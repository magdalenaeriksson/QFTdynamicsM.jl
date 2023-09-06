###############################################################################################################################
# Implementation of SimData for 2PI -> F(t,t'), r(t,t'), SigF(t'), Sigr(t') for N^d
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
using FFTW
using StaticArrays
export TwoPIGaugeScalarSimData
export expandSimData!

function rangesplitter(idxmin::Int,n::Int, k::Int)
    xz = Int.(ceil.(range(0, n, length = k+1)))
    return [(xz[i]+idxmin):(xz[i+1]+(idxmin-1)) for i in 1:k]
end
#
# Define SimData 
#
abstract type TwoPIGaugeScalarSimData <: AbstractSimData end
abstract type TwoPIGaugeScalarSimDataCPU <: TwoPIGaugeScalarSimData end

export TwoPIGaugeScalarSimDataCPUfull
struct TwoPIGaugeScalarSimDataCPUfull{sdim} <: TwoPIGaugeScalarSimDataCPU
    FS::SymMatrix{lattice{sdim}}
    FT::SymMatrix{lattice{sdim}}
    FL::SymMatrix{lattice{sdim}}
    rS::SymMatrix{lattice{sdim}}
    rT::SymMatrix{lattice{sdim}}
    rL::SymMatrix{lattice{sdim}}
    # self-energies
    SigFS::CyclVector{lattice{sdim}}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)
    SigFT::CyclVector{lattice{sdim}}
    SigFL::CyclVector{lattice{sdim}}
    SigrS::CyclVector{lattice{sdim}}
    SigrT::CyclVector{lattice{sdim}}
    SigrL::CyclVector{lattice{sdim}}
    # gauge sunset self-energies -- for energy computation
    SigFTgaugesunset::CyclVector{lattice{sdim}}
    SigFLgaugesunset::CyclVector{lattice{sdim}}
    SigrTgaugesunset::CyclVector{lattice{sdim}}
    SigrLgaugesunset::CyclVector{lattice{sdim}}
    # masses
    scalarmass2::Vector{Float64}
    transvmass2::lattice{sdim}
    longitmass2::lattice{sdim}
    #EfromSelf::CyclVector{Float64}#{sdim}
    # various momenta
    k2::lattice{sdim}
    k4::lattice{sdim}
    k6::lattice{sdim}
    k2inv::lattice{sdim}
    # omega values
    omega2Svalues::lattice{sdim}
    omega2Tvalues::lattice{sdim}
    omega2Lvalues::lattice{sdim}
    # lattice dimensions and numbers
    NstepsinMemory::Integer
    Nx::Integer
    sdim::Integer
    indices::Array{Int64} # [tmin, tmax]
    # FFT plans
    fftplan::FFTW.FFTWPlan
    #rfftplan::FFTW.rFFTWPlan#AbstractFFTs.Plan
    bfftplan::AbstractFFTs.Plan
    #brfftplan::FFTW.rFFTWPlan
    #irfftplan::AbstractFFTs.ScaledPlan
end

function TwoPIGaugeScalarSimDataCPUfull(NstepsinMemory, Nx, sdim)
    return TwoPIGaugeScalarSimDataCPUfull{ sdim }(
    #propagators
    SymMatrix(NstepsinMemory, Nx, sdim), SymMatrix(NstepsinMemory, Nx, sdim), SymMatrix(NstepsinMemory, Nx, sdim),
    SymMatrix(NstepsinMemory, Nx, sdim), SymMatrix(NstepsinMemory, Nx, sdim), SymMatrix(NstepsinMemory, Nx, sdim),
    # self-energies
    CyclVector(NstepsinMemory, Nx, sdim), CyclVector(NstepsinMemory, Nx, sdim), CyclVector(NstepsinMemory, Nx, sdim),
    CyclVector(NstepsinMemory, Nx, sdim), CyclVector(NstepsinMemory, Nx, sdim), CyclVector(NstepsinMemory, Nx, sdim),
    # gauge sunset self-energies
    CyclVector(NstepsinMemory, Nx, sdim), CyclVector(NstepsinMemory, Nx, sdim), CyclVector(NstepsinMemory, Nx, sdim), CyclVector(NstepsinMemory, Nx, sdim),
    # masses
    [0.], 0. * createlattice(Nx, sdim), 0. * createlattice(Nx, sdim),
    # various momenta
    createlattice(Nx, sdim), createlattice(Nx, sdim), createlattice(Nx, sdim), createlattice(Nx, sdim),
    # omega values
    createlattice(Nx, sdim), createlattice(Nx, sdim), createlattice(Nx, sdim),
    # lattice numbers
    NstepsinMemory, Nx, sdim, [0,0],
    # FFT plans
    plan_fft(createlattice(Nx, sdim)),
    #plan_rfft(createlattice(Nx, sdim)),
    plan_bfft(createlattice(Nx, sdim))
    #plan_brfft(createrfftlattice(Nx, sdim)+ im*createrfftlattice(Nx, sdim), Nx ), #NB may need to change this to plan_brfft( createrfftlattice(Nx, sdim) + im*createrfftlattice(Nx, sdim) , Nx ) 
                                                 # if Gerhard made createrfftlattice real-valued
    #plan_irfft(createrfftlattice(Nx, sdim), Nx)
    )
end

function expandSimData!(x::TwoPIGaugeScalarSimDataCPU)
    # change indices list of SimData
    if (x.indices[2] + 1) > x.NstepsinMemory
        x.indices[1] += 1
        x.indices[2] += 1
    else
        x.indices[1] = 1
        x.indices[2] += 1
    end
    # expand individual objects
    expand!(x.FS)
    expand!(x.FT)
    expand!(x.FL)
    expand!(x.rS)
    expand!(x.rT)
    expand!(x.rL)
    
    expand!(x.SigFS)
    expand!(x.SigFT)
    expand!(x.SigFL)
    expand!(x.SigrS)
    expand!(x.SigrT)
    expand!(x.SigrL)

    expand!(x.SigFTgaugesunset)
    expand!(x.SigFLgaugesunset)
    expand!(x.SigrTgaugesunset)
    expand!(x.SigrLgaugesunset)

    #expand!(x.EfromSelf)
end


export AbstractTmpData
export TwoPIGaugeScalarTmpData
abstract type AbstractTmpData end
abstract type TwoPIGaugeScalarTmpData <: AbstractTmpData end
abstract type TwoPIGaugeScalarTmpData2 <: AbstractTmpData end

export TwoPIGaugeScalarTmpDataCPUfullSimple
struct TwoPIGaugeScalarTmpDataCPUfullSimple{sdim} <: TwoPIGaugeScalarTmpData 
    # ranges for threads
    threadranges::Vector{UnitRange{Int64}}
    # tmp storage for computations
    scalarRHS::Vector{lattice{sdim}}
    transvRHS::Vector{lattice{sdim}}
    longitRHS::Vector{lattice{sdim}}
    nchunks::Integer
    # temp lattice vectors
    tmp::Array
    res::Array
    #
    G1T::Array
    G1L::Array
    G2T::Array
    G2L::Array
    function TwoPIGaugeScalarTmpDataCPUfullSimple(simdata::TwoPIGaugeScalarSimDataCPU, nchunks)
        @unpack Nx, sdim = simdata

        return new{sdim}(
            [0:0 for i in 1:nchunks],
            [createlattice(Nx, sdim)     for i in 1:nchunks], 
            [createlattice(Nx, sdim)     for i in 1:nchunks], 
            [createlattice(Nx, sdim)     for i in 1:nchunks], 
            nchunks,
            0*createlattice(Nx, sdim), 0*createlattice(Nx, sdim), 
            0*createlattice(Nx, sdim), 0*createlattice(Nx, sdim),
            0*createlattice(Nx, sdim), 0*createlattice(Nx, sdim))
    end
end

export TwoPIGaugeScalarTmpDataCPUfull
struct TwoPIGaugeScalarTmpDataCPUfull{sdim} <: TwoPIGaugeScalarTmpData 
    # ranges for threads
    threadranges::Vector{UnitRange{Int64}}
    # tmp storage for computations
    scalarRHS::Vector{lattice{sdim}}
    transvRHS::Vector{lattice{sdim}}
    longitRHS::Vector{lattice{sdim}}
    #tmplattice::Vector{lattice{sdim}} # vector containing 2 elements (lattices)
    nchunks::Integer
    # Gauss constraint matrices
    rhok::Vector{clattice{sdim}}
    chik::Vector{clattice{sdim}}
    rhox::Vector{clattice{sdim}}
    chix::Vector{clattice{sdim}}
    #
    # Initial real field scalar modes in x-space
    phix::Vector
    piix::Vector
    # Initial scalar field propagator helpers
    Pix::Array
    Phix::Array
    Pik::Array
    Phik::Array
    PiPik::Array
    PiPhik::Array
    PhiPhik::Array
    Ea::Array
    # Pauli matrix vector
    pauli::Vector
    # nvalues 
    nvalues::Vector
    # Initial gauge field propagator helpers
    Aaix::Vector{Vector{clattice}}#Vector{clattice{sdim}}
    Aaik::Vector{Vector{clattice}}#Vector{clattice{sdim}}
    Dk::Matrix
    DLk::clattice{sdim}
    DTk::clattice{sdim}
    # longitudinal projection operator
    P_L::Matrix
    # temp lattice vectors
    tmp::Array
    res::Array
    #
    G1T::Array
    G1L::Array
    G2T::Array
    G2L::Array
    #tmplatvec::Vector
    # CS averaged initial propagators
    AvrgPhiPhik::Vector
    AvrgPiPhik::Vector
    AvrgPiPik::Vector
    AvrgDTk::Vector
    AvrgDLk::Vector
    function TwoPIGaugeScalarTmpDataCPUfull(simdata::TwoPIGaugeScalarSimDataCPU, nchunks)
        @unpack Nx, sdim = simdata
        SU2 = SMatrix{2,2,ComplexF64}
        colouridx = SVector{sdim,Float64}   # a number for each value of spacetime, direction and adjoint index
        EaObj = SVector{sdim,colouridx}     # each directional index carries three values of adjoint indices
        #
        rhok = 0*[createclattice(Nx, sdim),createclattice(Nx, sdim),createclattice(Nx, sdim),createclattice(Nx, sdim)]
        chik = 0*[createclattice(Nx, sdim),createclattice(Nx, sdim),createclattice(Nx, sdim)]
        rhox = 0*[createlattice(Nx, sdim),createlattice(Nx, sdim),createlattice(Nx, sdim),createlattice(Nx, sdim)]
        chix = 0*[createlattice(Nx, sdim),createlattice(Nx, sdim),createlattice(Nx, sdim)]
        #
        phix = 0*[createlattice(Nx, sdim), createlattice(Nx, sdim), createlattice(Nx, sdim), createlattice(Nx, sdim)]
        piix = 0*[createlattice(Nx, sdim), createlattice(Nx, sdim), createlattice(Nx, sdim), createlattice(Nx, sdim)]
        #
        Pix   = Array{SU2}(undef,Nx,Nx,Nx)
        Phix  = Array{SU2}(undef,Nx,Nx,Nx)
        Pik   = Array{SU2}(undef,Nx,Nx,Nx)
        Phik  = Array{SU2}(undef,Nx,Nx,Nx)
        PiPik = createlattice(Nx, sdim)
        PiPhik = createlattice(Nx, sdim)
        PhiPhik = createlattice(Nx, sdim)
        Ea   = Array{EaObj}(undef,Nx,Nx,Nx)
        #
        pauli1 =  [0.0 1.0; 1.0 0.0]
        pauli2 =  [0.0 -1.0*im; 1.0*im 0.0]
        pauli3 =  [1.0 0.0; 0.0 -1.0]
        pauli =  [pauli1, pauli2, pauli3]
        #
        nvalues = zeros(Int64,Nx)
	    for i in 0:(Nx-1)
	        if (i>Nx/2)
	            nvalues[i+1] = -(Nx-i) # julia indexing
	        else
	            nvalues[i+1] = i  # julia indexing
	        end
	    end
        #
        Aaix   = Vector{Vector{clattice}}(undef,3)
        Aaik   = Vector{Vector{clattice}}(undef,3)
        Dk     = Matrix{clattice}(undef,3,3)
        DLk    = 0*createclattice(Nx, sdim)
        DTk    = 0*createclattice(Nx, sdim)
        P_L    = Matrix{clattice}(undef,3,3)
        for a in 1:3
            Aaix[a] = 0 * [ createclattice(Nx, sdim), createclattice(Nx, sdim), createclattice(Nx, sdim) ]
            Aaik[a] = 0 * [ createclattice(Nx, sdim), createclattice(Nx, sdim), createclattice(Nx, sdim) ]
        end
        for i in 1:3
            for j in 1:3
                Dk[i,j]  = 0 * createclattice(Nx, sdim) 
                P_L[i,j] = 0 * createclattice(Nx, sdim)
            end
        end
        #
        # fill P_L
        for i in 1:3
            for j in 1:3
                for nx in 1:Nx
                    for ny in 1:Nx
                        for nz in 1:Nx
                            p_phys = [ 2*sin(pi*nvalues[nx]/Nx), 2*sin(pi*nvalues[ny]/Nx), 2*sin(pi*nvalues[nz]/Nx)]
                            p2 = sum( p_phys .^2 )
                            
                            if p2 == 0
                                P_L[i,j][nx,ny,nz] = 1
                            else
                                P_L[i,j][nx,ny,nz] =  p_phys[i] * p_phys[j] / p2
                            end
                        end
                    end
                end
            end
         end
        tmp = 0*createlattice(Nx, sdim)
        res = 0*createlattice(Nx, sdim)
        G1T = 0*createlattice(Nx, sdim)
        G1L = 0*createlattice(Nx, sdim)
        G2T = 0*createlattice(Nx, sdim)
        G2L = 0*createlattice(Nx, sdim)
        #tmplatvec = 0*im*[createlattice(Nx, sdim), createlattice(Nx, sdim), createlattice(Nx, sdim)]
        #
        AvrgPhiPhik = 0 * [PhiPhik]
        AvrgPiPhik = 0 * [PhiPhik]
        AvrgPiPik = 0 * [PhiPhik]
        AvrgDTk = 0 * [PhiPhik]
        AvrgDLk = 0 * [PhiPhik]
        return new{sdim}(
            [0:0 for i in 1:nchunks],
            [createlattice(Nx, sdim)     for i in 1:nchunks], 
            [createlattice(Nx, sdim)     for i in 1:nchunks], 
            [createlattice(Nx, sdim)     for i in 1:nchunks], nchunks,
            #[0. * createlattice(Nx, sdim), 0. * createlattice(Nx, sdim)], nchunks,
            rhok, chik, rhox, chix, phix, piix, Pix, Phix, Pik, Phik, PiPik, PiPhik, PhiPhik, Ea, pauli, nvalues,
            Aaix, Aaik, Dk, DLk, DTk, P_L, tmp, res, G1T, G1L, G2T, G2L,
            AvrgPhiPhik, AvrgPiPhik, AvrgPiPik, AvrgDTk, AvrgDLk )
    end
end

##########################################################################################################################
export getTwoPIGaugeSimData
function getTwoPIGaugeSimData(problem::QFTdynamicsProblem)
    if typeof(problem.num) == QFTdynamics.TwoPIGaugeScalarCPUfull
        simdata = TwoPIGaugeScalarSimDataCPUfull(problem.simsetup.NstepsinMemory, problem.disc.Nx, problem.disc.sdim)
        for element in problem.disc.fftwhelper
            for idx in element.ind
                simdata.k2[idx] = element.lev2
                simdata.k4[idx] = element.lev2^2
                simdata.k6[idx] = element.lev2^3
                simdata.k2inv[idx] = 1/element.lev2
                #simdata.omega2Svalues[idx] = element.lev2 + problem.disc.Mass^2
            end
        end
        simdata.k2inv[1] = 0.
    end
    return simdata
end
##########################################################################################################################
