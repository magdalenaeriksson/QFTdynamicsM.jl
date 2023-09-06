using Parameters
using FFTW
using Roots
using ForwardDiff

# Initialization
export initSimData!
function initSimData!(initmass::Float64, simdata::TwoPIScalarSimDataCPUfull, tmpdata::TwoPIScalarTmpDataCPUfull, model::TwoPIScalarPhi4, pexp::TwoPIScalarPertExpansion, disc::TwoPIScalarDiscretization, init::TwoPIScalarInitialization)
    println("Initialzising with mass: ",  initmass); flush(stdout)
    # get n
    n = getparticlenr(init, disc)
    # init simdata with n
    expandSimData!(simdata) # expand simdataata twice since we populate up to index 2,2
    expandSimData!(simdata) # expand simdataata twice since we populate up to index 2,2
    for i in 1:length(disc.fftwhelper)
        omega = sqrt(disc.fftwhelper[i].lev2 + initmass^2)
	    for deg in 1:disc.fftwhelper[i].deg
            idx = disc.fftwhelper[i].ind[deg]
            simdata.F[1,1, idx] = (1/omega) * (n[i] + 0.5)
            simdata.F[2,1, idx] = (1/omega) * (n[i] + 0.5) 
            simdata.F[2,2, idx] = (  omega  * (n[i] + 0.5) * disc.dt^2 + (1/omega) * (n[i] + 0.5))
            simdata.r[1,1, idx] = 0 
            simdata.r[2,1, idx] = disc.dt * thesign(2,1, simdata.NstepsinMemory)
            simdata.r[2,2, idx] = 0
        end
    end
    simdata.hartreemass2[1] = getHartreeMass2(model, pexp, disc, simdata, 2, 2)
    println("HartreeMass2 at t=2: ",  simdata.hartreemass2[1]); flush(stdout)
end

function initSimData!(initmass::Float64, simdata::TwoPIScalarSimDataCPUred, tmpdata::TwoPIScalarTmpDataCPUred, model::TwoPIScalarPhi4, pexp::TwoPIScalarPertExpansion, disc::TwoPIScalarDiscretization, init::TwoPIScalarInitialization)
    println("Initialzising with mass: ",  initmass); flush(stdout)
    # get n
    n = getparticlenr(init, disc)
    # init simdata with n
    expandSimData!(simdata) # expand simdataata twice since we populate up to index 2,2
    expandSimData!(simdata) # expand simdataata twice since we populate up to index 2,2
    omega = sqrt.(simdata.kL2values .+ initmass^2) 
    simdata.F[1,1] .= [ 1/o * (n + 0.5) for (o,n) in zip(omega,n)]
    simdata.F[2,1] .= [ 1/o * (n + 0.5) for (o,n) in zip(omega,n)]
    simdata.F[2,2] .= [ o * (n + 0.5) * disc.dt^2 + 1/o * (n + 0.5) for (o,n) in zip(omega,n)]
    simdata.r[1,1] .= 0
    simdata.r[2,1] .= disc.dt * thesign(2,1, simdata.NstepsinMemory)
    simdata.r[2,2] .= 0
    simdata.hartreemass2[1] = getHartreeMass2(model, pexp, disc, simdata, 2, 2)
    println("HartreeMass2 at t=2: ",  simdata.hartreemass2[1]); flush(stdout)
end

function initSimData!(initmass::Float64, simdata::TwoPIScalarSimDataCPUred2, tmpdata::TwoPIScalarTmpDataCPUred, model::TwoPIScalarPhi4, pexp::TwoPIScalarPertExpansion, disc::TwoPIScalarDiscretization, init::TwoPIScalarInitialization)
    println("Initialzising with mass: ",  initmass); flush(stdout)
    # get n
    n = getparticlenr(init, disc)
    # init simdata with n
    expandSimData!(simdata) # expand simdataata twice since we populate up to index 2,2
    expandSimData!(simdata) # expand simdataata twice since we populate up to index 2,2
    omega = sqrt.(simdata.kL2values .+ initmass^2) 
    simdata.F[1,1] .= [ 1/o * (n + 0.5) for (o,n) in zip(omega,n)]
    simdata.F[2,1] .= [ 1/o * (n + 0.5) for (o,n) in zip(omega,n)]
    simdata.F[2,2] .= [ o * (n + 0.5) * disc.dt^2 + 1/o * (n + 0.5) for (o,n) in zip(omega,n)]
    simdata.r[1,1] .= 0
    simdata.r[2,1] .= disc.dt * thesign(2,1, simdata.NstepsinMemory)
    simdata.r[2,2] .= 0

    ichunk = 1
    calc_Fr_F2kr2!(1,1,simdata, tmpdata, disc, ichunk)
    calc_Fr_F2kr2!(2,1,simdata, tmpdata, disc, ichunk)
    calc_Fr_F2kr2!(2,2,simdata, tmpdata, disc, ichunk)
    simdata.hartreemass2[1] = getHartreeMass2(model, pexp, disc, simdata, 2, 2)
    println("HartreeMass2 at t=2: ",  simdata.hartreemass2[1]); flush(stdout)
end

function initSimData!(initmass::Float64, simdata::TwoPIScalarSimDataCPUcont, tmpdata::TwoPIScalarTmpDataCPUcont, model::TwoPIScalarPhi4, pexp::TwoPIScalarPertExpansion, disc::TwoPIScalarDiscretization, init::TwoPIScalarInitialization)
    println("Initialzising with mass: ",  initmass); flush(stdout)
    # get n
    n = getparticlenr(init, disc)
    # init simdata with n
    expandSimData!(simdata) # expand simdataata twice since we populate up to index 2,2
    expandSimData!(simdata) # expand simdataata twice since we populate up to index 2,2
    omega = sqrt.(simdata.kL2values .+ initmass^2) 
    simdata.F[1,1] .= [ 1/o * (n + 0.5) for (o,n) in zip(omega,n)]
    simdata.F[2,1] .= [ 1/o * (n + 0.5) for (o,n) in zip(omega,n)]
    simdata.F[2,2] .= [ o * (n + 0.5) * disc.dt^2 + 1/o * (n + 0.5) for (o,n) in zip(omega,n)]
    simdata.r[1,1] .= 0
    simdata.r[2,1] .= disc.dt * thesign(2,1, simdata.NstepsinMemory)
    simdata.r[2,2] .= 0

    ichunk = 1
    calc_Fr_F2kr2!(1,1,simdata, tmpdata, disc, ichunk)
    calc_Fr_F2kr2!(2,1,simdata, tmpdata, disc, ichunk)
    calc_Fr_F2kr2!(2,2,simdata, tmpdata, disc, ichunk)
    simdata.hartreemass2[1] = getHartreeMass2(model, pexp, disc, simdata, 2, 2)
    println("HartreeMass2 at t=2: ",  simdata.hartreemass2[1]); flush(stdout)
end

function initSimData!(initmass::Float64, simdata::TwoPIScalarSimDataCPUfull, tmpdata::TwoPIScalarTmpDataCPUfull, model::TwoPIScalarPhi4, pexp::TwoPIScalarPertExpansion, disc::TwoPIScalarDiscretization, init::TwoPIScalarGauss)
    expandSimData!(simdata) # expand simdataata twice since we populate up to index 2,2
    expandSimData!(simdata) # expand simdataata twice since we populate up to index 2,2
    for i in 1:length(disc.fftwhelper)
	    for deg in 1:disc.fftwhelper[i].deg
            idx = disc.fftwhelper[i].ind[deg]
            simdata.F[1,1, idx] = init.xi^2
            simdata.F[2,1, idx] = init.xi * init.eta * disc.dt + simdata.F[1,1, idx]
            simdata.F[2,2, idx] = (init.eta^2 + init.sig^2/(4*init.xi^2) ) * disc.dt^2 - simdata.F[1,1, idx] + 2 * simdata.F[2,1, idx]
            simdata.r[1,1, idx] = 0 
            simdata.r[2,1, idx] = disc.dt * thesign(2,1, simdata.NstepsinMemory)
            simdata.r[2,2, idx] = 0
        end
    end
    println("Conserved Energy: ",  0.5 * init.eta^2 + init.sig^2/(8*init.xi^2) + 0.5* model.Mass2 * init.xi^2 + (1/24) * model.Lambda * (model.ONgroup+2)/model.ONgroup * init.xi^4 ); flush(stdout)
    simdata.hartreemass2[1] = getHartreeMass2(model, pexp, disc, simdata, 2, 2)
    println("HartreeMass2 at t=2: ",  simdata.hartreemass2[1]); flush(stdout)
end

#function initSimData!(simdata::TwoPIScalarSimDataGPUreduced, model::Phi4, pexp::TwoPIScalarPertExpansion, disc::TwoPIScalarDiscretization, init::TwoPIScalarInitialization)
#    # get n
#    n = getparticlenr(init, disc)
#    # init simdata with n
#    expandSimData!(simdata) # expand simdataata twice since we populate up to index 2,2
#    expandSimData!(simdata) # expand simdataata twice since we populate up to index 2,2
#    for i in 1:length(disc.fftwhelper)
#        omega = sqrt(disc.fftwhelper[i].lev2 + disc.Mass^2)
#        simdata.F[1,1,i] = (1/omega) * (n[i] + 0.5)
#        simdata.F[2,1,i] = (1/omega) * (n[i] + 0.5) 
#        simdata.F[2,2,i] = (  omega  * (n[i] + 0.5) * disc.dt^2 + (1/omega) * (n[i] + 0.5))
#        simdata.r[1,1,i] = 0
#        simdata.r[2,1,i] = disc.dt * thesign(2,1, simdata.NstepsinMemory)
#        simdata.r[2,2,i] = 0
#    end
#    simdata.hartreemass2[1] = getHartreeMass2(model, pexp, disc, simdata, 2, 2)
#    println("HartreeMass2 at t=2: ",  simdata.hartreemass2[1]); flush(stdout)
#end

export getparticlenr
function getparticlenr(init::TwoPIScalarParticle, disc::TwoPIScalarDiscretization)
    return init.n .* ones(length(disc.fftwhelper))
end

function getparticlenr(init::TwoPIScalarThermal, disc::TwoPIScalarDiscretization)
    return [ 1/(exp(sqrt(disc.fftwhelper[i].lev2 + disc.Mass^2)/(init.T*disc.Mass))-1) for i in 1:length(disc.fftwhelper)]
end

function getparticlenr(init::TwoPIScalarTopHatT1, disc::TwoPIScalarDiscretizationLattice)
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

function getparticlenr(init::TwoPIScalarTopHatT1, disc::TwoPIScalarDiscretizationCont)
    n = zeros(length(disc.fftwhelper))
    k2min = 0.68 * 3 * disc.Mass^2
    idxmin = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2min))[2]
    k2max = 2.04 * 3 * disc.Mass^2
    idxmax = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2max))[2]
    eta = 2
    for i in idxmin:idxmax
        n[i] = eta 
    end
    return n
end

function getparticlenr(init::TwoPIScalarTopHatT2, disc::TwoPIScalarDiscretization)
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

function getparticlenr(init::TwoPIScalarTopHatT3, disc::TwoPIScalarDiscretization)
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

function getparticlenr(init::TwoPIScalarBox, disc::TwoPIScalarDiscretizationCont)
    n = zeros(length(disc.fftwhelper))
    k2max = 1.0
    idxmax = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2max))[2]
    eta = 1000
    for i in 1:idxmax
        n[i] = eta 
    end
    return n
end

# Hartreemass
export getHartreeMass2
function getHartreeMass2(model::TwoPIScalarPhi4, pexp::TwoPIScalarPertExpansion, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUfull, idxone::Int64, idxtwo::Int64)
    #Hartreemass of F[idxone,idxtwo]; in sumulation this is F[t-1,t-1]
    HMass2 = 0
    for i in 1:length(disc.fftwhelper)
	    for deg in 1:disc.fftwhelper[i].deg
            HMass2 += simdata.F[idxone, idxtwo, disc.fftwhelper[i].ind[deg]] 
        end
    end
    HMass2 *= Hartreeprefactor(model, pexp) * disc.ivol
    HMass2 += model.Mass2
    return HMass2
end

function getHartreeMass2(model::TwoPIScalarPhi4, pexp::TwoPIScalarPertExpansion, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUred, idxone::Int64, idxtwo::Int64)
    #Hartreemass of F[idxone,idxtwo]; in sumulation this is F[t-1,t-1]
    HMass2 = 0
    for i in 1:length(disc.fftwhelper)
        HMass2 += simdata.F[idxone, idxtwo, i] * disc.fftwhelper[i].deg 
    end
    HMass2 *= Hartreeprefactor(model, pexp) * disc.ivol
    HMass2 += model.Mass2
    return HMass2
end

function getHartreeMass2(model::TwoPIScalarPhi4, pexp::TwoPIScalarPertExpansion, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUred2, idxone::Int64, idxtwo::Int64)
    #Hartreemass of F[idxone,idxtwo]; in sumulation this is F[t-1,t-1]
    HMass2 = Hartreeprefactor(model, pexp) * simdata.Fx[idxone, idxtwo][1] # index 1 -> x=0 contribution
    HMass2 += model.Mass2
    return HMass2
end

function getHartreeMass2(model::TwoPIScalarPhi4, pexp::TwoPIScalarPertExpansion, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUcont, idxone::Int64, idxtwo::Int64)
    #Hartreemass of F[idxone,idxtwo]; in sumulation this is F[t-1,t-1]
    HMass2 = Hartreeprefactor(model, pexp) * simdata.Fx[idxone, idxtwo][1] # index 1 -> x=0 contribution
    HMass2 += model.Mass2
    return HMass2
end

function Hartreeprefactor(model::TwoPIScalarPhi4, pexp::TwoPIScalarLambdaEx)    return model.Lambda * (model.ONgroup+2)/(6*model.ONgroup) end
function Hartreeprefactor(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseLO)  return model.Lambda * (1/6)   end
function Hartreeprefactor(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO) return model.Lambda * (model.ONgroup+2)/(6*model.ONgroup) end

export getGapMass
function getGapMass(model::TwoPIScalarPhi4, pexp::TwoPIScalarPertExpansion, disc::TwoPIScalarDiscretization, init::TwoPIScalarInitialization)
    n = getparticlenr(init, disc) # note: this contains potentially the gapmass
    # define Gapequation
    gapequation(x) = begin
        result = 0 
        for i in 1:length(disc.fftwhelper)
            result += (1/sqrt(disc.fftwhelper[i].lev2 + x^2)) * (n[i] + 0.5) * disc.fftwhelper[i].deg
        end
        return result*disc.ivol*Hartreeprefactor(model, pexp) + model.Mass2 - x^2
    end
    # define derivative of Gapequation
    D(gapequation) = x -> ForwardDiff.derivative(gapequation,float(x))
    # find root of Gap equation
    gapMass = find_zero((gapequation,D(gapequation)), 0.5, Roots.Newton()) 
    println("GapMass2 (used for initialisation at LO fixed point): ",  gapMass^2); flush(stdout)
    return gapMass
end

export getThermalGapMass
function getThermalGapMass(model::TwoPIScalarPhi4, pexp::TwoPIScalarPertExpansion, disc::TwoPIScalarDiscretization, init::TwoPIScalarInitialization, mu)
    # mu is dim full!!
    println("getThermalGapMass with selfconsistent n(MH), CAUTION: dimension of mu")
    # define Gapequation
    gapequation(x) = begin
        result = 0 
        n = [ 1/(exp( (sqrt(disc.fftwhelper[i].lev2 + x^2) - mu)/(init.T*disc.Mass) ) - 1 ) for i in 1:length(disc.fftwhelper)]
        for i in 1:length(disc.fftwhelper)
            result += (1/sqrt(disc.fftwhelper[i].lev2 + x^2)) * (n[i]) * disc.fftwhelper[i].deg
        end
        return result*disc.ivol*Hartreeprefactor(model, pexp) + disc.Mass^2 - x^2
    end
    # define derivative of Gapequation
    D(gapequation) = x -> ForwardDiff.derivative(gapequation,float(x))
    # find root of Gap equation
    gapMass = find_zero((gapequation,D(gapequation)), disc.Mass, Roots.Newton()) 
    return gapMass/disc.Mass
end

export initialize!
function initialize!( thesolution::QFTdynamicsSolutionTwoPIScalar, tmpdata::TwoPIScalarTmpDataCPU)
    @unpack problem, simdata, measurearray =thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem
    
    if isnewsimulation(simsetup) == true
        # mass to initialize with
        if typeof(model) == QFTdynamics.Phi4Tachyonic
            initSimData!(disc.Mass, simdata, tmpdata, model, pexp, disc, init)
        elseif typeof(init) == QFTdynamics.TwoPIScalarGauss
            initSimData!(disc.Mass, simdata, tmpdata, model, pexp, disc, init)
        elseif typeof(init) == QFTdynamics.TwoPIScalarQuench
            initSimData!(disc.Mass, simdata, tmpdata, model, pexp, disc, init)
            #println("Special: Quench Initialzising with mass: ",  sqrt(2)*disc.Mass); flush(stdout)
            #println("Special: Quench Initialzising with mass2: ",  2*disc.Mass^2); flush(stdout)
            gapMass = getGapMass(model, pexp, disc, init)
            println("Quench: M(init)^2=, ", 2*disc.Mass^2, ", then drop to M^2=",disc.Mass^2, " , LO dynamics should bring it to ", gapMass^2)
        elseif typeof(init) == QFTdynamics.TwoPIScalarTopHatT1
            initSimData!(disc.Mass, simdata, tmpdata, model, pexp, disc, init)
        elseif typeof(init) == QFTdynamics.TwoPIScalarThermal
            initSimData!(disc.Mass, simdata, tmpdata, model, pexp, disc, init)
        elseif typeof(init) == QFTdynamics.TwoPIScalarTsunami
            println("reno-mass: ", disc.Mass); flush(stdout)
            println("bare-mass2: ", model.Mass2); flush(stdout)
            initmass = getGapMass(model, pexp, disc, init)
            println("init-mass: ", initmass); flush(stdout)
            initSimData!(initmass, simdata, tmpdata, model, pexp, disc, init)
        elseif typeof(init) == QFTdynamics.TwoPIScalarBox
            initSimData!(disc.Mass, simdata, tmpdata, model, pexp, disc, init)
        else
            initmass = getGapMass(model, pexp, disc, init)
            initSimData!(initmass, simdata, tmpdata, model, pexp, disc, init)
        end
    end

    simsetup.lastEvolStep = 2
    measure!(thesolution,2)
end