using Parameters
using FFTW

export funa!
"""
    funa!(a,b)
    Thats a doc string
"""
function funa!(a,b)
    return a+b+1
end

"""
    funa!(a::Float64,b)
    Thats a doc string for Float64
"""
function funa!(a::Float64,b)
    return a+b+1
end



export evolve!
"""
    evolve!(thesolution::QFTdynamicsSolutionGenericModel,t)
    evolve function
"""
function evolve!(thesolution::QFTdynamicsSolutionGenericModel, t)
    @unpack problem, simdata, measurearray = thesolution
    @unpack model, pexp, disc, init, reno, simsetup = problem

    println("Working on:" * string(t)); flush(stdout)
    expandSimData!(simdata) 

    # solve x_i''(t) = -omega^2 x_i(t)
    for idx in 1:(disc.Nx^disc.sdim)
        simdata.quantitytoevolve[t][idx] = (2 - model.omega^2*disc.dt^2) * simdata.quantitytoevolve[t-1][idx] - simdata.quantitytoevolve[t-2][idx]
    end

    #update last EvolStep
    simsetup.lastEvolStep = t 
end