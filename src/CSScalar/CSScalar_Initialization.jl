using Parameters
using FFTW
using Roots
using ForwardDiff
using Random

# Initialization
function initphi!(initmass::Float64, tmpdata::CSScalarTmpDataCPU, model::CSScalarPhi4, disc::CSScalarDiscretization, n::Vector, rng::Random.MersenneTwister)
    norm = sqrt(disc.Nx^disc.sdim)
    for j in 1:model.ONgroup
        tmpdata.phi_cpu[j] .= randn(rng, Float64, size(tmpdata.phi_cpu[j])) # Specified Float64 to get No(0,1), imaginary part is 0 (from initmemory)
        # inplace FFTW  
        #fft!(tmpdata.phi_cpu[j]) # creates asym one
        tmpdata.ftplan * tmpdata.phi_cpu[j] # creates asym one
        tmpdata.phi_cpu[j] ./= norm # now we have sym
        # modify FFTWs
        @inbounds for i in 1:length(disc.fftwhelper)
            omega = sqrt(disc.fftwhelper[i].lev2 + initmass^2)
            for deg in 1:disc.fftwhelper[i].deg
                tmpdata.phi_cpu[j][disc.fftwhelper[i].ind[deg]] *= sqrt( (n[i]+0.5) /omega)
            end
        end
        # inplace inverse FFTW  
        #ifft!(tmpdata.phi_cpu[j])
        tmpdata.iftplan * tmpdata.phi_cpu[j] # creates asym one
        tmpdata.phi_cpu[j] .*= norm # now we have sym
        # get rid of eventual imaginary parts
        tmpdata.phi_cpu[j] .= real(tmpdata.phi_cpu[j])
    end
    return
end

function initpii!(initmass::Float64, tmpdata::CSScalarTmpDataCPU, model::CSScalarPhi4, disc::CSScalarDiscretization, n::Vector, rng::Random.MersenneTwister)
    norm = sqrt(disc.Nx^disc.sdim)
    for j in 1:model.ONgroup
        tmpdata.pii_cpu[j] .= randn(rng, Float64, size(tmpdata.pii_cpu[j])) # Specified Float64 to get No(0,1), imaginary part is 0 (from initmemory)
        # inplace FFTW  
        #fft!(tmpdata.pii_cpu[j]) # creates asym one
        tmpdata.ftplan * tmpdata.pii_cpu[j] # creates asym one
        tmpdata.pii_cpu[j] ./= norm # now we have sym
        # modify FFTWs
        @inbounds for i in 1:length(disc.fftwhelper)
            omega = sqrt(disc.fftwhelper[i].lev2 + initmass^2)
            for deg in 1:disc.fftwhelper[i].deg
                tmpdata.pii_cpu[j][disc.fftwhelper[i].ind[deg]] *= sqrt( (n[i] + 0.5 )*omega)
            end
        end
        # inplace inverse FFTW  
        #ifft!(tmpdata.pii_cpu[j])
        tmpdata.iftplan * tmpdata.pii_cpu[j] # creates asym one
        tmpdata.pii_cpu[j] .*= norm # now we have sym
        # get rid of eventual imaginary parts
        tmpdata.pii_cpu[j] .= real(tmpdata.pii_cpu[j])
    end
    return
end

# CSPhi4Tachyonic -> only initialize unstable modes
function initphi!(initmass::Float64, tmpdata::CSScalarTmpDataCPU, model::CSPhi4TachyonicUNS, disc::CSScalarDiscretization, n::Vector, rng::Random.MersenneTwister)
    norm = sqrt(disc.Nx^disc.sdim)
    for j in 1:model.ONgroup
        tmpdata.phi_cpu[j] .= randn(rng, Float64, size(tmpdata.phi_cpu[j])) # Specified Float64 to get No(0,1), imaginary part is 0 (from initmemory)
        # inplace FFTW  
        #fft!(tmpdata.phi_cpu[j]) # creates asym one
        tmpdata.ftplan * tmpdata.phi_cpu[j] # creates asym one
        tmpdata.phi_cpu[j] ./= norm # now we have sym
        # modify FFTWs
        @inbounds for i in 1:length(disc.fftwhelper)
            if disc.fftwhelper[i].lev < disc.Mass
                omega = sqrt(disc.fftwhelper[i].lev2 + initmass^2)
                for deg in 1:disc.fftwhelper[i].deg
                    tmpdata.phi_cpu[j][disc.fftwhelper[i].ind[deg]] *= sqrt( (n[i]+0.5) /omega)
                end
            else
                for deg in 1:disc.fftwhelper[i].deg
                    tmpdata.phi_cpu[j][disc.fftwhelper[i].ind[deg]] *= 0
                end
            end
        end
        # inplace inverse FFTW  
        #ifft!(tmpdata.phi_cpu[j])
        tmpdata.iftplan * tmpdata.phi_cpu[j] # creates asym one
        tmpdata.phi_cpu[j] .*= norm # now we have sym
        # get rid of eventual imaginary parts
        tmpdata.phi_cpu[j] .= real(tmpdata.phi_cpu[j])
    end
    return
end

# CSPhi4Tachyonic -> only initialize unstable modes
function initpii!(initmass::Float64, tmpdata::CSScalarTmpDataCPU, model::CSPhi4TachyonicUNS, disc::CSScalarDiscretization, n::Vector, rng::Random.MersenneTwister)
    norm = sqrt(disc.Nx^disc.sdim)
    for j in 1:model.ONgroup
        tmpdata.pii_cpu[j] .= randn(rng, Float64, size(tmpdata.pii_cpu[j])) # Specified Float64 to get No(0,1), imaginary part is 0 (from initmemory)
        # inplace FFTW  
        #fft!(tmpdata.pii_cpu[j]) # creates asym one
        tmpdata.ftplan * tmpdata.pii_cpu[j] # creates asym one
        tmpdata.pii_cpu[j] ./= norm # now we have sym
        # modify FFTWs
        @inbounds for i in 1:length(disc.fftwhelper)
            if disc.fftwhelper[i].lev < disc.Mass
                omega = sqrt(disc.fftwhelper[i].lev2 + initmass^2)
                for deg in 1:disc.fftwhelper[i].deg
                    tmpdata.pii_cpu[j][disc.fftwhelper[i].ind[deg]] *= sqrt( (n[i] + 0.5 )*omega)
                end
            else
                for deg in 1:disc.fftwhelper[i].deg
                    tmpdata.pii_cpu[j][disc.fftwhelper[i].ind[deg]] *= 0
                end
            end
        end
        # inplace inverse FFTW  
        #ifft!(tmpdata.pii_cpu[j])
        tmpdata.iftplan * tmpdata.pii_cpu[j] # creates asym one
        tmpdata.pii_cpu[j] .*= norm # now we have sym
        # get rid of eventual imaginary parts
        tmpdata.pii_cpu[j] .= real(tmpdata.pii_cpu[j])
    end
    return
end

# CSPhi4TachyonicALL -> initialize all modes
function initphi!(initmass::Float64, tmpdata::CSScalarTmpDataCPU, model::CSPhi4TachyonicALL, disc::CSScalarDiscretization, n::Vector, rng::Random.MersenneTwister)
    norm = sqrt(disc.Nx^disc.sdim)
    for j in 1:model.ONgroup
        tmpdata.phi_cpu[j] .= randn(rng, Float64, size(tmpdata.phi_cpu[j])) # Specified Float64 to get No(0,1), imaginary part is 0 (from initmemory)
        # inplace FFTW  
        #fft!(tmpdata.phi_cpu[j]) # creates asym one
        tmpdata.ftplan * tmpdata.phi_cpu[j] # creates asym one
        tmpdata.phi_cpu[j] ./= norm # now we have sym
        # modify FFTWs
        @inbounds for i in 1:length(disc.fftwhelper)
            omega = sqrt(disc.fftwhelper[i].lev2 + initmass^2)
            for deg in 1:disc.fftwhelper[i].deg
                tmpdata.phi_cpu[j][disc.fftwhelper[i].ind[deg]] *= sqrt( (n[i]+0.5) /omega)
            end
        end
        # inplace inverse FFTW  
        #ifft!(tmpdata.phi_cpu[j])
        tmpdata.iftplan * tmpdata.phi_cpu[j] # creates asym one
        tmpdata.phi_cpu[j] .*= norm # now we have sym
        # get rid of eventual imaginary parts
        tmpdata.phi_cpu[j] .= real(tmpdata.phi_cpu[j])
    end
    return
end

# CSPhi4TachyonicALL -> initialize all modes
function initpii!(initmass::Float64, tmpdata::CSScalarTmpDataCPU, model::CSPhi4TachyonicALL, disc::CSScalarDiscretization, n::Vector, rng::Random.MersenneTwister)
    norm = sqrt(disc.Nx^disc.sdim)
    for j in 1:model.ONgroup
        tmpdata.pii_cpu[j] .= randn(rng, Float64, size(tmpdata.pii_cpu[j])) # Specified Float64 to get No(0,1), imaginary part is 0 (from initmemory)
        # inplace FFTW  
        #fft!(tmpdata.pii_cpu[j]) # creates asym one
        tmpdata.ftplan * tmpdata.pii_cpu[j] # creates asym one
        tmpdata.pii_cpu[j] ./= norm # now we have sym
        # modify FFTWs
        @inbounds for i in 1:length(disc.fftwhelper)
            omega = sqrt(disc.fftwhelper[i].lev2 + initmass^2)
            for deg in 1:disc.fftwhelper[i].deg
                tmpdata.pii_cpu[j][disc.fftwhelper[i].ind[deg]] *= sqrt( (n[i] + 0.5 )*omega)
            end
        end
        # inplace inverse FFTW  
        #ifft!(tmpdata.pii_cpu[j])
        tmpdata.iftplan * tmpdata.pii_cpu[j] # creates asym one
        tmpdata.pii_cpu[j] .*= norm # now we have sym
        # get rid of eventual imaginary parts
        tmpdata.pii_cpu[j] .= real(tmpdata.pii_cpu[j])
    end
    return
end


export initSimDataseriell!
function initSimDataseriell!(initmass::Float64, simdata::Vector{CSScalarSimDataCPU}, tmpdata::Vector{CSScalarTmpDataCPU}, model::CSScalarPhi4, num::CSScalarCPU, disc::CSScalarDiscretization, init::CSScalarInitialization, rng::Random.MersenneTwister)
    println("Initialzising with mass: ",  initmass); flush(stdout)
    # get n
    n = getparticlenr(init, disc)
    for i in 1:num.Runs 
        # init fields in tmpdata _cpu (sequentially!)
        initphi!(initmass, tmpdata[1], model, disc, n, rng)
        initpii!(initmass, tmpdata[1], model, disc, n, rng)
        # copy to fields to be evolved in simdata
        copy_toSimdata!(model, simdata[i].phi, tmpdata[1].phi_cpu)
        copy_toSimdata!(model, simdata[i].pii, tmpdata[1].pii_cpu)
    end
end

export initSimDataparallel!
function initSimDataparallel!(initmass::Float64, simdata::Vector{CSScalarSimDataCPU}, tmpdata::Vector{CSScalarTmpDataCPU}, model::CSScalarPhi4, num::CSScalarCPU, disc::CSScalarDiscretization, init::CSScalarInitialization, rng::Random.MersenneTwister)
    println("Initialzising with mass: ",  initmass); flush(stdout)
    # get n
    n = getparticlenr(init, disc)

    @Threads.threads for ichunk in 1:num.threads
        for i in num.threadranges[ichunk]
            # init fields in tmpdata _cpu (sequentially!)
            initphi!(initmass, tmpdata[ichunk], model, disc, n, rng)
            initpii!(initmass, tmpdata[ichunk], model, disc, n, rng)
            # copy to fields to be evolved in simdata
            copy_toSimdata!(model, simdata[i].phi, tmpdata[ichunk].phi_cpu)
            copy_toSimdata!(model, simdata[i].pii, tmpdata[ichunk].pii_cpu)
        end
    end
end

export getparticlenr
function getparticlenr(init::CSScalarParticle, disc::CSScalarDiscretization)
    return init.n .* ones(length(disc.fftwhelper))
end

function getparticlenr(init::CSScalarThermal, disc::CSScalarDiscretization)
    return [ 1/(exp(sqrt(disc.fftwhelper[i].lev2 + disc.Mass^2)/(init.T*disc.Mass))-1) for i in 1:length(disc.fftwhelper)]
end

function getparticlenr(init::CSScalarTopHatT1, disc::CSScalarDiscretization)
    n = zeros(length(disc.fftwhelper))
    k2min = 0.68 * disc.sdim * disc.Mass^2
    idxmin = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2min))[2]
    k2max = 2.04 * disc.sdim * disc.Mass^2
    idxmax = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2max))[2]
    eta = 2
    for i in idxmin:idxmax
        n[i] = eta 
    end
    return n
end

function getparticlenr(init::CSScalarTopHatT2, disc::CSScalarDiscretization)
    n = zeros(length(disc.fftwhelper))
    k2min = 0 * disc.sdim * disc.Mass^2
    idxmin = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2min))[2]
    k2max = 1.9 * disc.sdim * disc.Mass^2
    idxmax = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2max))[2]
    eta = 1.85
    for i in idxmin:idxmax
        n[i] = eta 
    end
    return n
end

function getparticlenr(init::CSScalarTopHatT3, disc::CSScalarDiscretization)
    n = zeros(length(disc.fftwhelper))
    k2min = 2.04 * disc.sdim * disc.Mass^2
    idxmin = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2min))[2]
    k2max = 2.72 * disc.sdim * disc.Mass^2
    idxmax = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2max))[2]
    eta = 1.6
    for i in idxmin:idxmax
        n[i] = eta 
    end
    return n
end

export initialize!
function initialize!( thesolution::QFTdynamicsSolutionCSScalar, tmpdata::Vector{CSScalarTmpDataCPU})
    @unpack problem, simdata, measurearray, measurearrayofruns =thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem

    if num.seed == 0 
        rng = MersenneTwister( )
    else
        rng = MersenneTwister( num.seed )
    end
    println("used seed: ", rng.seed)

    initSimDataseriell!(disc.Mass, simdata, tmpdata, model, num, disc, init, rng)
    #initSimDataparallel(disc.Mass, simdata, tmpdata, model, num, disc, init, rng)

    simsetup.lastEvolStep = 2
    measure!(thesolution,tmpdata,2)
end