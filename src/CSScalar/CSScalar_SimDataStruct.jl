using FFTW

#
# helpers to copy between SimData and TmpData
#
function copyON_toSimdata!(model::CSScalarModel,field::Array{Float64, 1},field_cpu::Array{ComplexF64, 1}) field[2:end-1] = real(field_cpu) end
function copyON_toSimdata!(model::CSScalarModel,field::Array{Float64, 2},field_cpu::Array{ComplexF64, 2}) field[2:end-1,2:end-1] = real(field_cpu) end
function copyON_toSimdata!(model::CSScalarModel,field::Array{Float64, 3},field_cpu::Array{ComplexF64, 3}) field[2:end-1,2:end-1,2:end-1] = real(field_cpu) end
export copy_toSimdata!
function copy_toSimdata!(model::CSScalarModel,field,field_cpu)
    for i in 1:model.ONgroup copyON_toSimdata!(model,field[i],field_cpu[i]) end
end

function copyON_toTmpdata!(model::CSScalarModel,field_cpu::Array{ComplexF64, 1},field::Array{Float64, 1}) field_cpu .= field[2:end-1] end
function copyON_toTmpdata!(model::CSScalarModel,field_cpu::Array{ComplexF64, 2},field::Array{Float64, 2}) field_cpu .= field[2:end-1,2:end-1] end
function copyON_toTmpdata!(model::CSScalarModel,field_cpu::Array{ComplexF64, 3},field::Array{Float64, 3}) field_cpu .= field[2:end-1,2:end-1,2:end-1] end
export copy_toTmpdata!
function copy_toTmpdata!(model::CSScalarModel,field_cpu,field)
    for i in 1:model.ONgroup copyON_toTmpdata!(model,field_cpu[i],field[i]) end
end

function myupdate_ONhalo!(model::CSScalarModel,field::Array{Float64, 1})
    field[1] = field[end-1]; field[end] = field[2] #x dir
end
function myupdate_ONhalo!(model::CSScalarModel,field::Array{Float64, 2})
    field[1,:] .= field[end-1,:]; field[end,:] .= field[2,:] #x dir
    field[:,1] .= field[:,end-1]; field[:,end] .= field[:,2] #y dir
end
function myupdate_ONhalo!(model::CSScalarModel,field::Array{Float64, 3})
        field[1,:,:] .= field[end-1,:,:]; field[end,:,:] .= field[2,:,:] #x dir
        field[:,1,:] .= field[:,end-1,:]; field[:,end,:] .= field[:,2,:] #y dir
        field[:,:,1] .= field[:,:,end-1]; field[:,:,end] .= field[:,:,2] #z dir
end
export myupdate_halo!
function myupdate_halo!(model::CSScalarModel,field)
    for i in 1:model.ONgroup myupdate_ONhalo!(model,field[i]) end
end

#
# Define SimData 
#
abstract type CSScalarSimData <: AbstractSimData end

export CSScalarSimDataCPU
struct CSScalarSimDataCPU <: CSScalarSimData
    # quantities to evolve
    phi::Array
    pii::Array
    # helpers of SimData structures
    Nx::Integer
    sdim::Integer
    ONgroup::Integer

    function CSScalarSimDataCPU(Nx, sdim, ONgroup)
        if sdim==3
            phi  	= Array{typeof(zeros(Nx+2,Nx+2,Nx+2)), 1}(undef, ONgroup)
            pii  	= Array{typeof(zeros(Nx+2,Nx+2,Nx+2)), 1}(undef, ONgroup)
            for i in 1:ONgroup
                phi[i] = zeros(Nx+2,Nx+2,Nx+2)
                pii[i] = zeros(Nx+2,Nx+2,Nx+2)
            end
        end
        if sdim==2
            phi  	= Array{typeof(zeros(Nx+2,Nx+2)), 1}(undef, ONgroup)
            pii  	= Array{typeof(zeros(Nx+2,Nx+2)), 1}(undef, ONgroup)
            for i in 1:ONgroup
                phi[i] = zeros(Nx+2,Nx+2)
                pii[i] = zeros(Nx+2,Nx+2)
            end
        end
        if sdim==1
            phi  	= Array{typeof(zeros(Nx+2)), 1}(undef, ONgroup)
            pii  	= Array{typeof(zeros(Nx+2)), 1}(undef, ONgroup)
            for i in 1:ONgroup
                phi[i] = zeros(Nx+2)
                pii[i] = zeros(Nx+2)
            end
        end
        return new(phi,pii,Nx,sdim,ONgroup)
    end
end

export MemoryUseofCSScalarSimDataCPU
function MemoryUseofCSScalarSimDataCPU(Nx, sdim, ONgroup)
    Mem = 0
    # phi and pi 
    Mem += 2 * (Nx+2)^sdim * ONgroup * 8 * 10^(-9)
    return Mem
end

abstract type CSScalarTmpData <: AbstractTmpData end

export CSScalarTmpDataCPU
struct CSScalarTmpDataCPU <: CSScalarTmpData
    #
    phi_cpu::Array
    pii_cpu::Array
    ## FFT plans
    ftplan::FFTW.cFFTWPlan
    iftplan::AbstractFFTs.ScaledPlan
    # ranges for threads
    update::Array
    # helpers of SimData structures
    Nx::Integer
    sdim::Integer
    ONgroup::Integer

    function CSScalarTmpDataCPU(simdata::CSScalarSimDataCPU)
        ##################
        if simdata.sdim==3
            phi_cpu  	= Array{typeof(zeros(ComplexF64,simdata.Nx,simdata.Nx,simdata.Nx)), 1}(undef, simdata.ONgroup)
            pii_cpu  	= Array{typeof(zeros(ComplexF64,simdata.Nx,simdata.Nx,simdata.Nx)), 1}(undef, simdata.ONgroup)
            for i in 1:simdata.ONgroup
                phi_cpu[i] = zeros(ComplexF64,simdata.Nx,simdata.Nx,simdata.Nx)
                pii_cpu[i] = zeros(ComplexF64,simdata.Nx,simdata.Nx,simdata.Nx)
            end
        end
        if simdata.sdim==2
            phi_cpu  	= Array{typeof(zeros(ComplexF64,simdata.Nx,simdata.Nx)), 1}(undef, simdata.ONgroup)
            pii_cpu  	= Array{typeof(zeros(ComplexF64,simdata.Nx,simdata.Nx)), 1}(undef, simdata.ONgroup)
            for i in 1:simdata.ONgroup
                phi_cpu[i] = zeros(ComplexF64,simdata.Nx,simdata.Nx)
                pii_cpu[i] = zeros(ComplexF64,simdata.Nx,simdata.Nx)
            end
        end
        if simdata.sdim==1
            phi_cpu  	= Array{typeof(zeros(ComplexF64,simdata.Nx)), 1}(undef, simdata.ONgroup)
            pii_cpu  	= Array{typeof(zeros(ComplexF64,simdata.Nx)), 1}(undef, simdata.ONgroup)
            for i in 1:simdata.ONgroup
                phi_cpu[i] = zeros(ComplexF64,simdata.Nx)
                pii_cpu[i] = zeros(ComplexF64,simdata.Nx)
            end
        end
        ##################
        if simdata.sdim == 3
            if simdata.ONgroup>1
                update	= zeros(simdata.Nx+2,simdata.Nx+2,simdata.Nx+2)
            else
                update	= zeros(1)
            end
        end
        if simdata.sdim == 2
            if simdata.ONgroup>1
                update	= zeros(simdata.Nx+2,simdata.Nx+2)
            else
                update	= zeros(1)
            end
        end
        if simdata.sdim == 1
            if simdata.ONgroup>1
                update	= zeros(simdata.Nx+2)
            else
                update	= zeros(1)
            end
        end
        return new(phi_cpu,pii_cpu, plan_fft!(phi_cpu[1]; flags=FFTW.EXHAUSTIVE), plan_ifft!(phi_cpu[1]; flags=FFTW.EXHAUSTIVE),update,simdata.Nx,simdata.sdim,simdata.ONgroup)
    end
end

export MemoryUseofCSScalarTmpDataCPU
function MemoryUseofTwoPIScalarTmpDataCPU(Nx, sdim, ONgroup)
    Mem = 0
    Mem += 4 * (Nx)^sdim * ONgroup * 8 * 10^(-9)
    if ONgroup > 1
        Mem += Nx^sdim * 8 * 10^(-9)
    end
    return Mem
end