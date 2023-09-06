using Parameters
using FFTW

export initialize!

# Initialization
function initSimData!(simdata::GenericModelSimData, model::GenericModelA, pexp::anexpansion, disc::DiscretizationGenericModel, init::aninitialization)
    # init SimD with n
    expandSimData!(simdata) # expand SimData twice since we populate up to index 2,2
    expandSimData!(simdata) # expand SimData twice since we populate up to index 2,2

    for idx in 1:(disc.Nx^disc.sdim)
        simdata.quantitytoevolve[1][idx] = 0
        simdata.quantitytoevolve[2][idx] = 0.1
    end
    
end

function initialize!( thesolution::QFTdynamicsSolutionGenericModel)
    @unpack problem, simdata, measurearray = thesolution
    @unpack model, pexp, disc, init, reno, simsetup = problem
    if isnewsimulation(simsetup) == true
        initSimData!(simdata, model, pexp, disc, init)
    end

    simsetup.lastEvolStep = 2
    measure!(thesolution,2)
end