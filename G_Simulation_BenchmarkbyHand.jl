###################################################################################################################################
using QFTdynamics
using Serialization
using BenchmarkTools
using Parameters

###################################################################################################################################
# Set up
###################################################################################################################################
# get parameters as dict from parameter file and specify modelfile (a type)
#include("src/parseParameterfile.jl") # julia Simulation.jl --Parameterfile <<myparameterfile located in QFTdynamics>>
#include("Parfile_TwoPIScalar.jl") 
#include("src/GenericModel/GenericModel_Parfile.jl") 
include("src/TwoPIScalar/TwoPIScalar_Parfile.jl") 
parameters["Pexp"] = "NinverseNLO"
parameters["Num"] = "CPUred2"

parameters["Nsteps"] = 400
parameters["Nx"] = 32
parameters["sdim"] = 3
parameters["NstepsinMemory"] = 400
parameters["Nchunks"] = 1

# create Problem
theproblem = QFTdynamicsProblem(modelfile, parameters)

thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)

initialize!(thesolution, tmpdata)

#elapsedtime = @elapsed for t in (thesolution.problem.simsetup.lastEvolStep+1):(thesolution.problem.simsetup.Nsteps+2)
#    #@show thesolution.simdata.indices, ichunk::Int64
#    println("Working on:[" * string(t) * "," * string(t) * "] "); flush(stdout)
#
#    evolve!(thesolution, tmpdata, t)
#    if (t-2)%thesolution.problem.simsetup.NstepsbetweenNmeas == 0
#        measure!(thesolution,t)
#    end
#
#    #thesolution.problem.simsetup.lastEvolStep = t # update last EvolStep -> NOw insede Evolve
#end
for t in (thesolution.problem.simsetup.lastEvolStep+1):(thesolution.simdata.NstepsinMemory)
    expandSimData!(thesolution.simdata) # now indices range from 1:NstepsinMemory
    nrofindices = (thesolution.simdata.indices[2]-1) - thesolution.simdata.indices[1] + 1
    tmpdata.threadranges .= splitter(thesolution.simdata.indices[1], nrofindices, thesolution.problem.num.nchunks)
    fakecalcSelfEnergies!(thesolution.problem.model, thesolution.problem.pexp, thesolution.problem.disc, thesolution.simdata, tmpdata, 2)
end

# benchmarking

# evolve function
#@benchmark evolve!($thesolution, $tmpdata, $(thesolution.problem.simsetup.Nsteps+2)) seconds=100
@benchmark evolve!($thesolution, $tmpdata, $(thesolution.problem.simsetup.Nsteps+2))

## individual ones
@unpack problem, simdata, measurearray = thesolution
@unpack model, pexp, disc, init, reno, simsetup = problem

t=thesolution.problem.simsetup.Nsteps
tp=thesolution.problem.simsetup.Nsteps-5
#tp=150

# all functions in evolve!
#@benchmark calcSelfEnergies!($model, $pexp, $disc, $simdata, $tmpdata, t-1)
#@benchmark getHartreeMass2(  $model, $pexp, $disc, $simdata, t-1, t-1)
#@benchmark setomega2values!( $model, $simdata, $disc)
#@benchmark getRHS_r!($model, $pexp, $simdata, $disc, $tmpdata.RHS[1], t, tp)
#@benchmark evolve_r!($model, $simdata, $disc, $tmpdata.RHS[1], t, tp)
#@benchmark getRHS_F!($model, $pexp, $simdata, $disc, $tmpdata.RHS[1], t, tp)
#@benchmark evolve_F!($model, $simdata, $disc, $tmpdata.RHS[1], t, tp)
#@benchmark calc_Fr_F2kr2!(t,tp, $simdata, $tmpdata, $disc, 1)

## evolve_F
#@benchmark getRHS_r!($model, $pexp, $simdata, $disc, $tmpdata.RHS[1], t, tp)
#@benchmark optimgetRHS_r!($model, $pexp, $simdata, $disc, $tmpdata.RHS[1], t, tp)
#
#function optimgetRHS_r!(model::QFTdynamics.TwoPIScalarPhi4, pexp::QFTdynamics.TwoPIScalarNLO, simdata::QFTdynamics.TwoPIScalarSimDataCPUred, disc::QFTdynamics.TwoPIScalarDiscretization, Mem::QFTdynamics.reducedlattice, t::Int64, tp::Int64)
#    ## TO DO: REDUCE SIG
#    # calculates the RHS for solving for r(t,tp). Gives R_r(t-1,tp)
#    # Brute force implementation - no symmetries taken into account
#    # Note: r(t,t) is set by hand -> R_r(t-1,t) does not need to be calculated (tp=t > t-1)
#    Mem .= 0
#    for z in tp+1:t-2
#        Mem .-= simdata.Sigr[z] .* simdata.r[z,tp] .* thesign(z, tp, simdata.NstepsinMemory) 
#    end
#end
