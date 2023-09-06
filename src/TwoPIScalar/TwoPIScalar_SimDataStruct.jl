###############################################################################################################################
# Implementation of SimData for TwoPIScalar
#   - TwoPIScalarSimDataCPUfull{sdim} 
#   - TwoPIScalarSimDataGPUreduced{sdim} 
#
# Comments: 
#   * initialize SimData, expand it ("add a row"), populate the columns of the new row, usw 
#   * F/r implementation:
#       - set values: F/r[t,tp,k] where k referes to the lattice site to be set and t>=tp
#       - get values: F/r[t,tp,k]
#       - get lattice: F/r[t,tp]
#   * SigF/Sigr implementation
#       - set lattice: SigF/Sigr[t,tp]
#       - get lattice: SigF/Sigr[t,tp]
#   
# Note: For multiple component fields use array of F -> [F0,F1,F2,F4],
#       this is faster than SVector in FFTW 
#
###############################################################################################################################
using CUDA
using FFTW
export TwoPIScalarSimData
export expandSimData!
using SparseArrays

#
# helpers
#


# gpubatch -> real
const gpubatchlattice{sdimpone} = CuArray{Float64,sdimpone}
function getGPUMemory(Nx,sdim,NstepsinMemory)
    if sdim==1 return CuArray{Float64}(undef,Nx,NstepsinMemory) end
    if sdim==2 return CuArray{Float64}(undef,Nx,Nx,NstepsinMemory) end
    if sdim==3 return CuArray{Float64}(undef,Nx,Nx,Nx,NstepsinMemory) end
end

# gputmpbatch -> complex reduced
const gputmpbatchlattice{sdimpone} = CuArray{ComplexF64,sdimpone}
function getGPUtmpMemory(Nx,sdim,NstepsinMemory)
    if sdim==1 return CuArray{ComplexF64}(undef,Int(floor(Nx/2)+1),NstepsinMemory) end
    if sdim==2 return CuArray{ComplexF64}(undef,Int(floor(Nx/2)+1),Nx,NstepsinMemory) end
    if sdim==3 return CuArray{ComplexF64}(undef,Int(floor(Nx/2)+1),Nx,Nx,NstepsinMemory) end
end

const gpuinterfacelattice{sdimpone} = Array{Float64,sdimpone}
function creategpuinterface(Nx,sdim,NstepsinMemory)
    return Mem.alloc(Mem.Host, Nx^sdim*NstepsinMemory*sizeof(Float64))
end

#
# Define SimData 
#
abstract type TwoPIScalarSimData <: AbstractSimData end

abstract type TwoPIScalarSimDataCPU <: TwoPIScalarSimData end
abstract type TwoPIScalarSimDataGPU <: TwoPIScalarSimData end

##################################################################################
# for CPUfull computation
##################################################################################
struct TwoPIScalarSimDataCPUfull{sdim} <: TwoPIScalarSimDataCPU
    # quantities to evolve
    F::SymMatrix{lattice{sdim}}
    r::SymMatrix{lattice{sdim}}
    # self Energy calculation
    SigF::CyclVector{lattice{sdim}}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)
    Sigr::CyclVector{lattice{sdim}}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)
    IFx::CyclVector{rfftlattice{sdim}} #SelfEnergy computation
    Irx::CyclVector{rfftlattice{sdim}} #SelfEnergy computation
    ### tmp storage for computations
    hartreemass2::Vector{Float64}
    kL2values::lattice{sdim}
    omega2values::lattice{sdim}
    # helpers of SimData structures
    NstepsinMemory::Integer
    Nx::Integer
    sdim::Integer
    indices::Array{Int64} # [tmin, tmax]

    function TwoPIScalarSimDataCPUfull(NstepsinMemory, Nx, sdim)
        return new{ sdim }(SymMatrix(NstepsinMemory, Nx, sdim), SymMatrix(NstepsinMemory, Nx, sdim), CyclVector(NstepsinMemory, Nx, sdim), CyclVector(NstepsinMemory, Nx, sdim),
        CyclVectorx(NstepsinMemory, Nx, sdim), #IFx - not used in lambda expansion
        CyclVectorx(NstepsinMemory, Nx, sdim), #IFr - not used in lambda expansion
        [0.], createlattice(Nx, sdim), createlattice(Nx, sdim),
        NstepsinMemory, Nx, sdim, [0,0])
    end
end

export MemoryUseofTwoPIScalarSimDataCPUfull
function MemoryUseofTwoPIScalarSimDataCPUfull(NstepsinMemory, Nx, sdim)
    Mem = 0
    # F and rho
    Mem += NstepsinMemory^2 * Nx^sdim * 8 * 10^(-9)
    # SigF and Sigr
    Mem += 2 * NstepsinMemory * Nx^sdim * 8 * 10^(-9)
    # IFx and Irx (roughly half the size, but complex)
    Mem += 2 * NstepsinMemory * (Nx^sdim)* 8 * 10^(-9)
    return Mem
end

##################################################################################
# for CPUcont computation)
##################################################################################
struct TwoPIScalarSimDataCPUcont <: TwoPIScalarSimDataCPU
    # quantities to evolve
    F::SymMatrix{reducedlattice}
    r::SymMatrix{reducedlattice}
    Fr::SymMatrix{reducedlattice}
    F2kr2::SymMatrix{reducedlattice}
    # self Energy calculation
    SigF::CyclVector{reducedlattice}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)
    Sigr::CyclVector{reducedlattice}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)
    #SigFfull::CyclVector{lattice{sdim}}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)
    #Sigrfull::CyclVector{lattice{sdim}}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)
    IFx::CyclVector{reducedlattice} #SelfEnergy computation
    Irx::CyclVector{reducedlattice} #SelfEnergy computation
    IFk::CyclVector{reducedlattice}
    Irk::CyclVector{reducedlattice}
    Fx::CyclVector{reducedlattice} #SelfEnergy computation
    rx::CyclVector{reducedlattice} #SelfEnergy computation
    ### tmp storage for computations
    hartreemass2::Vector{Float64}
    kL2values::reducedlattice
    omega2values::reducedlattice
    # helpers of SimData structures
    NstepsinMemory::Integer
#    Nx::Integer
#    sdim::Integer
    indices::Array{Int64} # [tmin, tmax]

    function TwoPIScalarSimDataCPUcont(NstepsinMemory, nmomenta)
        return new(SymMatrix(NstepsinMemory, nmomenta), SymMatrix(NstepsinMemory, nmomenta), SymMatrix(NstepsinMemory, nmomenta), SymMatrix(NstepsinMemory, nmomenta), CyclVector(NstepsinMemory, nmomenta), CyclVector(NstepsinMemory, nmomenta),
        CyclVector(NstepsinMemory, nmomenta), #IFx - not used in lambda expansion
        CyclVector(NstepsinMemory, nmomenta), #IFr - not used in lambda expansion
        CyclVector(NstepsinMemory, nmomenta), CyclVector(NstepsinMemory, nmomenta), #IFk, Irk
        CyclVector(NstepsinMemory, nmomenta), 
        CyclVector(NstepsinMemory, nmomenta), 
        [0.],
        createreducedlattice( nmomenta), createreducedlattice( nmomenta),
        NstepsinMemory, [0,0])
    end
end

##################################################################################
# for CPUred computation)
##################################################################################
struct TwoPIScalarSimDataCPUred2{sdim} <: TwoPIScalarSimDataCPU
    # quantities to evolve
    F::SymMatrix{reducedlattice}
    r::SymMatrix{reducedlattice}
    Fr::SymMatrix{reducedlattice}
    F2kr2::SymMatrix{reducedlattice}
    # self Energy calculation
    SigF::CyclVector{reducedlattice}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)
    Sigr::CyclVector{reducedlattice}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)
    #SigFfull::CyclVector{lattice{sdim}}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)
    #Sigrfull::CyclVector{lattice{sdim}}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)
    IFx::CyclVector{rfftlattice{sdim}} #SelfEnergy computation
    Irx::CyclVector{rfftlattice{sdim}} #SelfEnergy computation
    IFk::CyclVector{reducedlattice}
    Irk::CyclVector{reducedlattice}
    Fx::CyclVector{rfftlattice{sdim}} #SelfEnergy computation
    rx::CyclVector{rfftlattice{sdim}} #SelfEnergy computation
    #Frxloc::CyclVector{rfftlattice{sdim}} #SelfEnergy computation
    #F2kr2xloc::CyclVector{rfftlattice{sdim}} #SelfEnergy computation
    ### tmp storage for computations
    hartreemass2::Vector{Float64}
    kL2values::reducedlattice
    omega2values::reducedlattice
    # helpers of SimData structures
    NstepsinMemory::Integer
    Nx::Integer
    sdim::Integer
    indices::Array{Int64} # [tmin, tmax]

    function TwoPIScalarSimDataCPUred2(NstepsinMemory, Nx, sdim, nmomenta)
        return new{ sdim }(SymMatrix(NstepsinMemory, nmomenta), SymMatrix(NstepsinMemory, nmomenta), SymMatrix(NstepsinMemory, nmomenta), SymMatrix(NstepsinMemory, nmomenta), CyclVector(NstepsinMemory, nmomenta), CyclVector(NstepsinMemory, nmomenta),
        CyclVectorx(NstepsinMemory, Nx, sdim), #IFx - not used in lambda expansion
        CyclVectorx(NstepsinMemory, Nx, sdim), #IFr - not used in lambda expansion
        CyclVector(NstepsinMemory, nmomenta), CyclVector(NstepsinMemory, nmomenta), #IFk, Irk
        CyclVectorx(NstepsinMemory, Nx, sdim), 
        CyclVectorx(NstepsinMemory, Nx, sdim), 
        [0.],
        createreducedlattice( nmomenta), createreducedlattice( nmomenta),
        NstepsinMemory, Nx, sdim, [0,0])
    end
end

export MemoryUseofTwoPIScalarSimDataCPUred2
function MemoryUseofTwoPIScalarSimDataCPUred2(NstepsinMemory, Nx, sdim)
    Mem = 0
    nmom = length(getfftwhelper(Nx,sdim))
    # F and rho, Fr, F2kr2
    Mem += 2 * NstepsinMemory^2 * nmom * 8 * 10^(-9)
    # SigF and Sigr, Irk and IFk
    Mem += 4 * NstepsinMemory * nmom * 8 * 10^(-9)
    # Fx, rx, IFx, Irx
    Mem += 4 * NstepsinMemory * (Nx^sdim) * 8 * 10^(-9)
    return Mem
end

struct TwoPIScalarSimDataCPUred{sdim} <: TwoPIScalarSimDataCPU
    # quantities to evolve
    F::SymMatrix{reducedlattice}
    r::SymMatrix{reducedlattice}
    # self Energy calculation
    SigF::CyclVector{reducedlattice}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)
    Sigr::CyclVector{reducedlattice}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)
    SigFfull::CyclVector{lattice{sdim}}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)
    Sigrfull::CyclVector{lattice{sdim}}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)

    IFx::CyclVector{rfftlattice{sdim}} #SelfEnergy computation
    Irx::CyclVector{rfftlattice{sdim}} #SelfEnergy computation
    IFk::CyclVector{reducedlattice}
    Irk::CyclVector{reducedlattice}

    Fxloc::CyclVector{rfftlattice{sdim}} #SelfEnergy computation
    rxloc::CyclVector{rfftlattice{sdim}} #SelfEnergy computation
    Frxloc::CyclVector{rfftlattice{sdim}} #SelfEnergy computation
    F2kr2xloc::CyclVector{rfftlattice{sdim}} #SelfEnergy computation
    ### tmp storage for computations
    hartreemass2::Vector{Float64}
    kL2values::reducedlattice
    omega2values::reducedlattice
    # helpers of SimData structures
    NstepsinMemory::Integer
    Nx::Integer
    sdim::Integer
    indices::Array{Int64} # [tmin, tmax]

    function TwoPIScalarSimDataCPUred(NstepsinMemory, Nx, sdim, nmomenta)
        return new{ sdim }(SymMatrix(NstepsinMemory, nmomenta), SymMatrix(NstepsinMemory, nmomenta), 
        CyclVector(NstepsinMemory, nmomenta), CyclVector(NstepsinMemory, nmomenta), CyclVector(NstepsinMemory, Nx, sdim), CyclVector(NstepsinMemory, Nx, sdim),
        CyclVectorx(NstepsinMemory, Nx, sdim), #IFx - not used in lambda expansion
        CyclVectorx(NstepsinMemory, Nx, sdim), #IFr - not used in lambda expansion
        CyclVector(NstepsinMemory, nmomenta), CyclVector(NstepsinMemory, nmomenta), #IFk, Irk
        CyclVectorx(NstepsinMemory, Nx, sdim), 
        CyclVectorx(NstepsinMemory, Nx, sdim), 
        CyclVectorx(NstepsinMemory, Nx, sdim), 
        CyclVectorx(NstepsinMemory, Nx, sdim), 
        [0.],
        createreducedlattice( nmomenta), createreducedlattice( nmomenta),
        NstepsinMemory, Nx, sdim, [0,0])
    end
end

export MemoryUseofTwoPIScalarSimDataCPUred
function MemoryUseofTwoPIScalarSimDataCPUred(NstepsinMemory, Nx, sdim)
    Mem = 0
    nmom = length(getfftwhelper(Nx,sdim))
    # F and rho
    Mem += NstepsinMemory^2 * nmom * 8 * 10^(-9)
    # reduced
    Mem += 4 * NstepsinMemory * nmom * 8 * 10^(-9)
    # full and rfft
    Mem += 8 * NstepsinMemory * Nx^sdim * 8 * 10^(-9)
    return Mem
end

function expandSimData!(x::TwoPIScalarSimDataCPUfull)
    # change indices list of SimData
    if (x.indices[2] + 1) > x.NstepsinMemory
        x.indices[1] += 1
        x.indices[2] += 1
    else
        x.indices[1] = 1
        x.indices[2] += 1
    end
    # expand individual objects
    expand!(x.F)
    expand!(x.r)
    expand!(x.SigF)
    expand!(x.Sigr)
    expand!(x.IFx)
    expand!(x.Irx)
end

function expandSimData!(x::TwoPIScalarSimDataCPUred)
    # change indices list of SimData
    if (x.indices[2] + 1) > x.NstepsinMemory
        x.indices[1] += 1
        x.indices[2] += 1
    else
        x.indices[1] = 1
        x.indices[2] += 1
    end
    # expand individual objects
    expand!(x.F)
    expand!(x.r)
    expand!(x.SigF)
    expand!(x.Sigr)
    expand!(x.SigFfull)
    expand!(x.Sigrfull)
    expand!(x.IFx)
    expand!(x.Irx)
    expand!(x.IFk)
    expand!(x.Irk)
end

function expandSimData!(x::TwoPIScalarSimDataCPUred2)
    # change indices list of SimData
    if (x.indices[2] + 1) > x.NstepsinMemory
        x.indices[1] += 1
        x.indices[2] += 1
    else
        x.indices[1] = 1
        x.indices[2] += 1
    end
    # expand individual objects
    expand!(x.F)
    expand!(x.r)
    expand!(x.Fr)
    expand!(x.F2kr2)
    expand!(x.SigF)
    expand!(x.Sigr)
    expand!(x.IFx)
    expand!(x.Irx)
    expand!(x.IFk)
    expand!(x.Irk)
    expand!(x.Fx)
    expand!(x.rx)
end

function expandSimData!(x::TwoPIScalarSimDataCPUcont)
    # change indices list of SimData
    if (x.indices[2] + 1) > x.NstepsinMemory
        x.indices[1] += 1
        x.indices[2] += 1
    else
        x.indices[1] = 1
        x.indices[2] += 1
    end
    # expand individual objects
    expand!(x.F)
    expand!(x.r)
    expand!(x.Fr)
    expand!(x.F2kr2)
    expand!(x.SigF)
    expand!(x.Sigr)
    expand!(x.IFx)
    expand!(x.Irx)
    expand!(x.IFk)
    expand!(x.Irk)
    expand!(x.Fx)
    expand!(x.rx)
end

export AbstractTmpData
export TwoPIScalarTmpData
export TwoPIScalarTmpDataCPU
abstract type AbstractTmpData end
abstract type TwoPIScalarTmpData <: AbstractTmpData end
abstract type TwoPIScalarTmpDataCPU <: TwoPIScalarTmpData end

export TwoPIScalarTmpDataCPUfull
struct TwoPIScalarTmpDataCPUfull{sdim} <: TwoPIScalarTmpDataCPU
    # FFT plans
    ftplan::FFTW.rFFTWPlan
    iftplan::FFTW.rFFTWPlan
    # ranges for threads
    threadranges::Vector{UnitRange{Int64}}
    threadranges2::Vector{UnitRange{Int64}}
    threadranges3::Vector{UnitRange{Int64}}
    # tmp storage for computations
    RHS::Vector{lattice{sdim}} #RHS
    RHS2::Vector{lattice{sdim}} #RHS
    RHS3::Vector{lattice{sdim}} #RHS
    rx::Vector{rfftlattice{sdim}} #SelfEnergy computation
    Fx::Vector{rfftlattice{sdim}} #SelfEnergy computation
    Fx2::Vector{rfftlattice{sdim}} #SelfEnergy computation
    # helpers of SimData structures
    NstepsinMemory::Integer
    Nx::Integer
    sdim::Integer
    nchunks::Integer

    function TwoPIScalarTmpDataCPUfull(simdata::TwoPIScalarSimDataCPU, nchunks)
        return new{ simdata.sdim }(
        plan_rfft(createlattice(simdata.Nx, simdata.sdim); flags=FFTW.EXHAUSTIVE), # just one representative 
        plan_brfft(createrfftlattice(simdata.Nx, simdata.sdim) + im*createrfftlattice(simdata.Nx, simdata.sdim), simdata.Nx; flags=FFTW.EXHAUSTIVE),
        [0:0 for i in 1:nchunks], [0:0 for i in 1:nchunks], [0:0 for i in 1:nchunks],
        [createlattice(simdata.Nx, simdata.sdim)     for i in 1:nchunks],
        [createlattice(simdata.Nx, simdata.sdim)     for i in 1:nchunks],
        [createlattice(simdata.Nx, simdata.sdim)     for i in 1:nchunks],
        [createrfftlattice(simdata.Nx, simdata.sdim) for i in 1:nchunks], 
        [createrfftlattice(simdata.Nx, simdata.sdim) for i in 1:nchunks], 
        [createrfftlattice(simdata.Nx, simdata.sdim) for i in 1:nchunks], 
        simdata.NstepsinMemory, simdata.Nx, simdata.sdim, nchunks)
    end
end

export MemoryUseofTwoPIScalarTmpDataCPUfull
function MemoryUseofTwoPIScalarTmpDataCPUfull(NstepsinMemory, Nx, sdim, nchunks)
    Mem = 0
    # RHS
    Mem += nchunks * Nx^sdim * 8 * 10^(-9)
    # Fx, rx
    Mem += 2 * nchunks * Nx^sdim * 0.5 * 8 * 10^(-9)
    return Mem
end

export TwoPIScalarTmpDataCPUred
struct TwoPIScalarTmpDataCPUred{sdim} <: TwoPIScalarTmpDataCPU
    # FFT plans
    ftplan::FFTW.rFFTWPlan
    iftplan::FFTW.rFFTWPlan
    # ranges for threads
    threadranges::Vector{UnitRange{Int64}}
    threadranges2::Vector{UnitRange{Int64}}
    threadranges3::Vector{UnitRange{Int64}}
    # tmp storage for computations
    RHS::Vector{reducedlattice} #RHS
    tmpfull::Vector{lattice{sdim}} # for Irs/IFs calc, for CPUfull i could use RHS; here that does not work
    tmpx::Vector{rfftlattice{sdim}} #SelfEnergy computation
################################################################################################################
    RHS2::Vector{reducedlattice} #RHS
    #RHS3::Vector{reducedlattice} #RHS
    tmpfull2::Vector{lattice{sdim}} # for Irs/IFs calc, for CPUfull i could use RHS; here that does not work
    tmpfull3::Vector{lattice{sdim}} # for Irs/IFs calc, for CPUfull i could use RHS; here that does not work
    r::Vector{lattice{sdim}} #for reconstruction of lattice from reduced one
    F::Vector{lattice{sdim}} #for reconstruction of lattice from reduced one
    rx::Vector{rfftlattice{sdim}} #SelfEnergy computation
    Fx::Vector{rfftlattice{sdim}} #SelfEnergy computation
    Fx2::Vector{rfftlattice{sdim}} #SelfEnergy computation
    # helpers for copying to lattice / reduced lattice
    tofulllattice::SparseMatrixCSC{Float64, Int64}
    toredlattice::SparseMatrixCSC{Float64, Int64}
    tmpvec::Vector{Float64}
################################################################################################################
    # helpers of SimData structures
    NstepsinMemory::Integer
    Nx::Integer
    sdim::Integer
    nchunks::Integer

    function TwoPIScalarTmpDataCPUred(simdata::TwoPIScalarSimDataCPU, fftwhelper, nchunks)
################################################################################################################
        toredlattice = spzeros(length(fftwhelper), Int(simdata.Nx^simdata.sdim))
        for i in 1:length(fftwhelper)
            for j in 1:Int(simdata.Nx^simdata.sdim)
                toredlattice[i,fftwhelper[i].ind[1]] = 1
            end
        end
        tofulllattice = spzeros(Int(simdata.Nx^simdata.sdim), length(fftwhelper))
        for i in 1:length(fftwhelper)
            for j in fftwhelper[i].ind
                tofulllattice[j,i] = 1
            end
        end
        tmpvec = zeros(Int(simdata.Nx^simdata.sdim))
################################################################################################################
        return new{ simdata.sdim }(
        plan_rfft(createlattice(simdata.Nx, simdata.sdim); flags=FFTW.EXHAUSTIVE), # just one representative 
        plan_brfft(createrfftlattice(simdata.Nx, simdata.sdim) + im*createrfftlattice(simdata.Nx, simdata.sdim), simdata.Nx; flags=FFTW.EXHAUSTIVE),
        [0:0 for i in 1:nchunks],[0:0 for i in 1:nchunks],[0:0 for i in 1:nchunks],
        [createreducedlattice( length(simdata.kL2values) )  for i in 1:nchunks], # RHS
        [createlattice(simdata.Nx, simdata.sdim)            for i in 1:nchunks], # tmpfull
        [createrfftlattice(simdata.Nx, simdata.sdim) for i in 1:nchunks], # tmpx
################################################################################################################
        [createreducedlattice( length(simdata.kL2values) )  for i in 1:nchunks],
        [createlattice(simdata.Nx, simdata.sdim)            for i in 1:nchunks],
        [createlattice(simdata.Nx, simdata.sdim)            for i in 1:nchunks],
        [createlattice(simdata.Nx, simdata.sdim)            for i in 1:nchunks],
        [createlattice(simdata.Nx, simdata.sdim)            for i in 1:nchunks],
        #[createclattice(simdata.Nx, simdata.sdim)           for i in 1:nchunks], 
        #[createclattice(simdata.Nx, simdata.sdim)           for i in 1:nchunks], 
        #[createclattice(simdata.Nx, simdata.sdim)           for i in 1:nchunks], 
        [createrfftlattice(simdata.Nx, simdata.sdim) for i in 1:nchunks], 
        [createrfftlattice(simdata.Nx, simdata.sdim) for i in 1:nchunks], 
        [createrfftlattice(simdata.Nx, simdata.sdim) for i in 1:nchunks], 
        tofulllattice, toredlattice, tmpvec,
################################################################################################################
        simdata.NstepsinMemory, simdata.Nx, simdata.sdim, nchunks)
    end
end

export MemoryUseofTwoPIScalarTmpDataCPUred
function MemoryUseofTwoPIScalarTmpDataCPUred(NstepsinMemory, Nx, sdim, nchunks)
    Mem = 0
    nmom = length(getfftwhelper(Nx,sdim))
    # RHS
    Mem += nchunks * nmom * 8 * 10^(-9)
    # tmpfull and tmpx
    Mem += 2 * nchunks * Nx^sdim * 8 * 10^(-9)
    return Mem
end

# CPUcont 
#function calcFx(F, nrmomenta)
#    isofactor = [i for i in 1:nrmomenta]
#    normfactor = (1/(4*nrmomenta^2)) * [1/(i-0.5) for i in 1:nrmomenta]
#    result = FFTW.r2r(isofactor .* F,FFTW.RODFT01)
#    return normfactor .* result
#end
#
#function calcF(Fx, nrmomenta)
#    isofactor = [(i-0.5) for i in 1:nrmomenta]
#    normfactor = 2*nrmomenta * [1/i for i in 1:nrmomenta]
#    result = FFTW.r2r(isofactor .* Fx,FFTW.RODFT10)
#    return normfactor .* result
#end

export TwoPIScalarTmpDataCPUcont
struct TwoPIScalarTmpDataCPUcont <: TwoPIScalarTmpDataCPU
    # FFT plans
    ftplan::FFTW.r2rFFTWPlan  # from p to x: FFTW.plan_r2r(r, FFTW.RODFT01),
    iftplan::FFTW.r2rFFTWPlan # from x to p: FFTW.plan_r2r(r, FFTW.RODFT10),
    # ranges for threads
    threadranges::Vector{UnitRange{Int64}}
    threadranges2::Vector{UnitRange{Int64}}
    threadranges3::Vector{UnitRange{Int64}}
    # tmp storage for computations
    RHS::Vector{reducedlattice} #for reconstruction of lattice from reduced one

#    RHS::Vector{reducedlattice} #RHS
#    tmpfull::Vector{lattice{sdim}} # for Irs/IFs calc, for CPUfull i could use RHS; here that does not work
#    tmpx::Vector{rfftlattice{sdim}} #SelfEnergy computation
#
#
#
#    RHS2::Vector{reducedlattice} #RHS
#    #RHS3::Vector{reducedlattice} #RHS
#    tmpfull2::Vector{lattice{sdim}} # for Irs/IFs calc, for CPUfull i could use RHS; here that does not work
#    tmpfull3::Vector{lattice{sdim}} # for Irs/IFs calc, for CPUfull i could use RHS; here that does not work
#    r::Vector{lattice{sdim}} #for reconstruction of lattice from reduced one
#    F::Vector{lattice{sdim}} #for reconstruction of lattice from reduced one
#    rx::Vector{rfftlattice{sdim}} #SelfEnergy computation
#    Fx::Vector{rfftlattice{sdim}} #SelfEnergy computation
#    Fx2::Vector{rfftlattice{sdim}} #SelfEnergy computation
    # helpers for cont fft
    isofactortox::Vector
    normfactortox::Vector
    isofactortop::Vector
    normfactortop::Vector
    # helpers of SimData structures
    NstepsinMemory::Integer
#    Nx::Integer
#    sdim::Integer
    nchunks::Integer

    function TwoPIScalarTmpDataCPUcont(simdata::TwoPIScalarSimDataCPU, fftwhelper, nchunks)
        r=rand(length(fftwhelper))

        return new(
        FFTW.plan_r2r(r, FFTW.RODFT01),
        FFTW.plan_r2r(r, FFTW.RODFT10),
        [0:0 for i in 1:nchunks],[0:0 for i in 1:nchunks],[0:0 for i in 1:nchunks],
        [createreducedlattice( length(simdata.kL2values) )  for i in 1:nchunks], # RHS

#        [createreducedlattice( length(simdata.kL2values) )  for i in 1:nchunks], # RHS
#        [createlattice(simdata.Nx, simdata.sdim)            for i in 1:nchunks], # tmpfull
#        [createrfftlattice(simdata.Nx, simdata.sdim) for i in 1:nchunks], # tmpx
#        [createreducedlattice( length(simdata.kL2values) )  for i in 1:nchunks],
#        [createlattice(simdata.Nx, simdata.sdim)            for i in 1:nchunks],
#        [createlattice(simdata.Nx, simdata.sdim)            for i in 1:nchunks],
#        [createlattice(simdata.Nx, simdata.sdim)            for i in 1:nchunks],
#        [createlattice(simdata.Nx, simdata.sdim)            for i in 1:nchunks],
#        #[createclattice(simdata.Nx, simdata.sdim)           for i in 1:nchunks], 
#        #[createclattice(simdata.Nx, simdata.sdim)           for i in 1:nchunks], 
#        #[createclattice(simdata.Nx, simdata.sdim)           for i in 1:nchunks], 
#        [createrfftlattice(simdata.Nx, simdata.sdim) for i in 1:nchunks], 
#        [createrfftlattice(simdata.Nx, simdata.sdim) for i in 1:nchunks], 
#        [createrfftlattice(simdata.Nx, simdata.sdim) for i in 1:nchunks], 
        # helpers for cont fft
        [i for i in 1:length(fftwhelper)], #isofactortox
        (1/(4*length(fftwhelper)^2)) * [1/(i-0.5) for i in 1:length(fftwhelper)], #normfactortox
        [(i-0.5) for i in 1:length(fftwhelper)], #isofactortop 
        2*length(fftwhelper) * [1/i for i in 1:length(fftwhelper)], #normfactortop 
        # general helpers
        simdata.NstepsinMemory, nchunks)
    end
end

export MemoryUseofTwoPIScalarTmpDataCPUred
function MemoryUseofTwoPIScalarTmpDataCPUred(NstepsinMemory, Nx, sdim, nchunks)
    Mem = 0
    nmom = length(getfftwhelper(Nx,sdim))
    # RHS
    Mem += nchunks * nmom * 8 * 10^(-9)
    # r, F, tmpfull
    Mem += 3 * nchunks * Nx^sdim * 8 * 10^(-9)
    # Fx, rx
    Mem += 2 * nchunks * Nx^sdim * 0.5 * 8 * 10^(-9)
    return Mem
end

##################################################################################
# for GPUreduced computation
##################################################################################
struct TwoPIScalarSimDataGPUreduced{sdim,sdimpone} <: TwoPIScalarSimDataGPU
    # quantities to evolve
    F::SymMatrix{reducedlattice}
    r::SymMatrix{reducedlattice}
    # self Energy calculation
    SigF::CyclVector{reducedlattice}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)
    Sigr::CyclVector{reducedlattice}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)
    F_GPUinterface::lattice{sdimpone} #continous memoryblock for transfer to/from GPU - one dimension higher than sdim
    r_GPUinterface::lattice{sdimpone} #continous memoryblock for transfer to/from GPU - one dimension higher than sdim
    F_GPUinterface_ptr::Ptr{Float64}
    r_GPUinterface_ptr::Ptr{Float64}
    F_GPUinterface_buf::CUDA.Mem.HostBuffer
    r_GPUinterface_buf::CUDA.Mem.HostBuffer
    # FFT plans
    ftplan::CUDA.CUFFT.rCuFFTPlan
    iftplan::AbstractFFTs.ScaledPlan
    # GPU Memory
    Fbatch_gpu::gpubatchlattice{sdimpone} #real -> Nx^sdim * NstepsinMemory 
    rbatch_gpu::gpubatchlattice{sdimpone} #real -> Nx^sdim * NstepsinMemory 
    Fxbatch_gpu::gputmpbatchlattice{sdimpone} #complex -> 2 *(Nx/2 + 1) * Nx^(sdim-1) * NstepsinMemory
    rxbatch_gpu::gputmpbatchlattice{sdimpone} #complex -> 2 *(Nx/2 + 1) * Nx^(sdim-1) * NstepsinMemory
    Fbatch_gpu_ptr::CuPtr{Float64}
    rbatch_gpu_ptr::CuPtr{Float64}
    Fxbatch_gpu_ptr::CuPtr{Float64}
    rxbatch_gpu_ptr::CuPtr{Float64}
    # tmp storage for computations
    tmplattice::Vector{reducedlattice}
    tmpscalar::Vector{Float64}
    kL2values::Vector{Float64}
    omega2values::Vector{Float64}
    # helpers of SimData structures
    NstepsinMemory::Integer
    Nx::Integer
    sdim::Integer
    indices::Array{Int64} # [tmin, tmax]

    function TwoPIScalarSimDataGPUreduced(NstepsinMemory, Nx, sdim, nmomenta)
        F_GPUinterface_buf = creategpuinterface(Nx, sdim, NstepsinMemory)
        r_GPUinterface_buf = creategpuinterface(Nx, sdim, NstepsinMemory)
        F_GPUinterface_ptr = convert(Ptr{Float64}, F_GPUinterface_buf)
        r_GPUinterface_ptr = convert(Ptr{Float64}, r_GPUinterface_buf)
        #get first GPU memory since it is needed for fft plan
        Fbatch_gpu =  getGPUMemory(Nx,sdim,NstepsinMemory) #real -> Nx^sdim * NstepsinMemory 
        rbatch_gpu =  getGPUMemory(Nx,sdim,NstepsinMemory) #real -> Nx^sdim * NstepsinMemory 
        Fxbatch_gpu = getGPUtmpMemory(Nx,sdim,NstepsinMemory) #complex -> 2 *(Nx/2 + 1) * Nx^(sdim-1) * NstepsinMemory
        rxbatch_gpu = getGPUtmpMemory(Nx,sdim,NstepsinMemory) #complex -> 2 *(Nx/2 + 1) * Nx^(sdim-1) * NstepsinMemory
    return new{ sdim, sdim+1 }(SymMatrix(NstepsinMemory, nmomenta), SymMatrix(NstepsinMemory, nmomenta), CyclVector(NstepsinMemory, nmomenta), CyclVector(NstepsinMemory, nmomenta),
        unsafe_wrap(Array, F_GPUinterface_ptr, tuple( map(i->Nx,1:sdim)...,NstepsinMemory) ), unsafe_wrap(Array, r_GPUinterface_ptr, tuple( map(i->Nx,1:sdim)...,NstepsinMemory) ),
        F_GPUinterface_ptr, r_GPUinterface_ptr,
        F_GPUinterface_buf, r_GPUinterface_buf,
        plan_rfft(rbatch_gpu,       map(i->i,1:sdim)), # just one representative 
        plan_irfft(rxbatch_gpu, Nx, map(i->i,1:sdim)), # operates on complex one
        Fbatch_gpu,          rbatch_gpu,            Fxbatch_gpu,            rxbatch_gpu,          
        pointer(Fbatch_gpu), pointer(rbatch_gpu),   pointer(Fxbatch_gpu),   pointer(rxbatch_gpu), 
        [createreducedlattice(nmomenta) for i in 1:Threads.nthreads()], [0.], 
        zeros(nmomenta), zeros(nmomenta),
        NstepsinMemory, Nx, sdim, [0,0])
    end
end

function expandSimData!(x::TwoPIScalarSimDataGPUreduced)
    # change indices list of SimData
    if (x.indices[2] + 1) > x.NstepsinMemory
        x.indices[1] += 1
        x.indices[2] += 1
    else
        x.indices[1] = 1
        x.indices[2] += 1
    end
    # expand individual objects
    expand!(x.F)
    expand!(x.r)
    expand!(x.SigF)
    expand!(x.Sigr)
end

##########################################################################################################################
function getTwoPISimData(problem::QFTdynamicsProblem)
    if typeof(problem.num) == QFTdynamics.TwoPIScalarCPUfull
        # same simdata for lambda and 1/N expansion!
        simdata = TwoPIScalarSimDataCPUfull(problem.simsetup.NstepsinMemory, problem.disc.Nx, problem.disc.sdim)
        for element in problem.disc.fftwhelper
            for idx in element.ind
                simdata.kL2values[idx] = element.lev2
                simdata.omega2values[idx] = element.lev2 + problem.disc.Mass^2
                simdata.omega2values[idx] *= problem.disc.dt^2
            end
        end
    end
    if typeof(problem.num) == QFTdynamics.TwoPIScalarCPUred
        # same simdata for lambda and 1/N expansion!
        simdata = TwoPIScalarSimDataCPUred(problem.simsetup.NstepsinMemory, problem.disc.Nx, problem.disc.sdim, length(problem.disc.fftwhelper))
        simdata.kL2values .= [element.lev2 for element in problem.disc.fftwhelper]
        simdata.omega2values .= simdata.kL2values .+ problem.disc.Mass^2   
        simdata.omega2values .*= problem.disc.dt^2
    end
    if typeof(problem.num) == QFTdynamics.TwoPIScalarCPUred2
        # same simdata for lambda and 1/N expansion!
        simdata = TwoPIScalarSimDataCPUred2(problem.simsetup.NstepsinMemory, problem.disc.Nx, problem.disc.sdim, length(problem.disc.fftwhelper))
        simdata.kL2values .= [element.lev2 for element in problem.disc.fftwhelper]
        simdata.omega2values .= simdata.kL2values .+ problem.disc.Mass^2   
        simdata.omega2values .*= problem.disc.dt^2
    end
    if typeof(problem.num) == QFTdynamics.TwoPIScalarCPUcont
        # same simdata for lambda and 1/N expansion!
        simdata = TwoPIScalarSimDataCPUcont(problem.simsetup.NstepsinMemory, length(problem.disc.fftwhelper))
        simdata.kL2values .= [element.lev2 for element in problem.disc.fftwhelper]
        simdata.omega2values .= simdata.kL2values .+ problem.disc.Mass^2   
        simdata.omega2values .*= problem.disc.dt^2
    end
    #if typeof(problem.num) == QFTdynamics.TwoPIScalarGPUreduced
    #    simdata = TwoPIScalarSimDataGPUreduced(problem.simsetup.NstepsinMemory, problem.disc.Nx, problem.disc.sdim, length(problem.disc.fftwhelper))
    #    simdata.kL2values .= [element.lev2 for element in problem.disc.fftwhelper]
    #    simdata.omega2values .= simdata.kL2values .+ problem.disc.Mass^2   
    #end
    return simdata
end

function getTwoPIScalarTmpData(simdata::TwoPIScalarSimDataCPUfull, problem::QFTdynamicsProblem)
    return TwoPIScalarTmpDataCPUfull(simdata, problem.num.nchunks)
end
function getTwoPIScalarTmpData(simdata::TwoPIScalarSimDataCPUred, problem::QFTdynamicsProblem)
    return TwoPIScalarTmpDataCPUred(simdata, problem.disc.fftwhelper, problem.num.nchunks)
end
function getTwoPIScalarTmpData(simdata::TwoPIScalarSimDataCPUred2, problem::QFTdynamicsProblem)
    return TwoPIScalarTmpDataCPUred(simdata, problem.disc.fftwhelper, problem.num.nchunks)
end
function getTwoPIScalarTmpData(simdata::TwoPIScalarSimDataCPUcont, problem::QFTdynamicsProblem)
    return TwoPIScalarTmpDataCPUcont(simdata, problem.disc.fftwhelper, problem.num.nchunks)
end
##########################################################################################################################
export calc_Fr_F2kr2!
function calc_Fr_F2kr2!(tmone::Int,tp::Int,simdata::TwoPIScalarSimDataCPUfull, tmpdata::TwoPIScalarTmpDataCPUfull, disc, ichunk::Int) end
function calc_Fr_F2kr2!(tmone::Int,tp::Int,simdata::TwoPIScalarSimDataCPUred, tmpdata::TwoPIScalarTmpDataCPUred, disc, ichunk::Int) end

function calc_Fr_F2kr2!(tmone::Int,tp::Int,simdata::TwoPIScalarSimDataCPUred2, tmpdata::TwoPIScalarTmpDataCPUred, disc, ichunk::Int)
    # calc Fx  
    copytolattice!( tmpdata.tmpfull[ichunk], simdata.F[tmone,tp], disc.fftwhelper)
    mul!(simdata.Fx[tp], tmpdata.ftplan, tmpdata.tmpfull[ichunk])
    simdata.Fx[tp] .*= disc.ivol
    simdata.Fx[tp] .= real.(simdata.Fx[tp])

    # calc rx  
    copytolattice!( tmpdata.tmpfull[ichunk], simdata.r[tmone,tp], disc.fftwhelper)
    mul!(simdata.rx[tp], tmpdata.ftplan, tmpdata.tmpfull[ichunk])
    simdata.rx[tp] .*= disc.ivol
    simdata.rx[tp] .= real.(simdata.rx[tp])

    # calc Fx * rx
    tmpdata.tmpx[ichunk] .=  simdata.Fx[tp] .* simdata.rx[tp] # no thesign(t,tp,simdata.NstepsinMemory) here, use it when i fetch Fr[t,tp]?
    # calc (Fx * rx)k
    mul!(tmpdata.tmpfull[ichunk], tmpdata.iftplan, tmpdata.tmpx[ichunk])
    copytoreducedlattice!( simdata.Fr[tmone,tp], tmpdata.tmpfull[ichunk], disc.fftwhelper)
    #simdata.Fr[tmone,tp] .= real.(simdata.Fr[tmone,tp])

    # calc Fx^2 - 0.25 * rx^2
    tmpdata.tmpx[ichunk] .=  simdata.Fx[tp].^2 .- 0.25 * simdata.rx[tp].^2
    # calc (Fx^2 - 0.25 * rx^2)k
    mul!(tmpdata.tmpfull[ichunk], tmpdata.iftplan, tmpdata.tmpx[ichunk])
    copytoreducedlattice!( simdata.F2kr2[tmone,tp], tmpdata.tmpfull[ichunk], disc.fftwhelper)
    #simdata.F2kr2[tmone,tp] .= real.(simdata.F2kr2[tmone,tp])
end

# CPUcont 
#function calc_Fr_F2kr2!(tmone::Int,tp::Int,simdata::TwoPIScalarSimDataCPUcont, tmpdata::TwoPIScalarTmpDataCPUcont, disc, ichunk::Int)
#    # calc Fx  
#    tmpdata.RHS[ichunk] .= tmpdata.isofactortox .* simdata.F[tmone,tp]
#    mul!(simdata.Fx[tp], tmpdata.ftplan, tmpdata.RHS[ichunk])
#    simdata.Fx[tp] .*= tmpdata.normfactortox
#
#    # calc rx  
#    tmpdata.RHS[ichunk] .= tmpdata.isofactortox .* simdata.r[tmone,tp]
#    mul!(simdata.rx[tp], tmpdata.ftplan, tmpdata.RHS[ichunk])
#    simdata.rx[tp] .*= tmpdata.normfactortox
#
#    # calc Fx * rx
#    tmpdata.RHS[ichunk] .=  simdata.Fx[tp] .* simdata.rx[tp] # no thesign(t,tp,simdata.NstepsinMemory) here, use it when i fetch Fr[t,tp]?
#    # calc (Fx * rx)k
#    tmpdata.RHS[ichunk] .*= tmpdata.isofactortop
#    mul!(simdata.Fr[tmone,tp], tmpdata.iftplan, tmpdata.RHS[ichunk])
#    simdata.Fr[tmone,tp] .*= tmpdata.normfactortop
#
#    # calc Fx^2 - 0.25 * rx^2
#    tmpdata.RHS[ichunk] .=  simdata.Fx[tp].^2 .- 0.25 * simdata.rx[tp].^2
#    # calc (Fx^2 - 0.25 * rx^2)k
#    tmpdata.RHS[ichunk] .*= tmpdata.isofactortop
#    mul!(simdata.F2kr2[tmone,tp], tmpdata.iftplan, tmpdata.RHS[ichunk])
#    simdata.F2kr2[tmone,tp] .*= tmpdata.normfactortop
#end

# copy to full/reducedlattice using LinAlg
export copytofulllattice!
function copytofulllattice!(alattice::Array{Float64}, areducedlattice::reducedlattice, tmpdata::TwoPIScalarTmpDataCPUred)
    #tmpdata.tmpvec .= tmpdata.tofulllattice * areducedlattice
    #mul!(tmpdata.tmpvec, tmpdata.tofulllattice, areducedlattice) #allocates 0bytes, 15 ns
    #alattice .= reshape(tmpdata.tmpvec, tmpdata.Nx, tmpdata.Nx) #allocates 96bytes, 25ns
    # sligthly better option: 
    vecalattice = @view alattice[:]
    mul!(vecalattice, tmpdata.tofulllattice, areducedlattice)
end

export copytoreducedlattice!
function copytoreducedlattice!(areducedlattice::reducedlattice, alattice, tmpdata::TwoPIScalarTmpDataCPUred)
    mul!(areducedlattice, tmpdata.toredlattice, @view alattice[:])
end
