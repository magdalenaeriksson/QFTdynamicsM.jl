###################################################################################################################################
using QFTdynamics
using Serialization

###################################################################################################################################
# Set up
###################################################################################################################################
# get parameters as dict from parameter file and specify modelfile (a type)
#include("src/parseParameterfile.jl") # julia Simulation.jl --Parameterfile <<myparameterfile located in QFTdynamics>>
#include("Parfile_TwoPIScalar.jl") 
#include("src/GenericModel/GenericModel_Parfile.jl") 
include("src/TwoPIScalar/TwoPIScalar_Parfile.jl") 

# create Problem
theproblem = QFTdynamicsProblem(modelfile, parameters)

thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)

initialize!(thesolution, tmpdata)

elapsedtime = @elapsed for t in (thesolution.problem.simsetup.lastEvolStep+1):(thesolution.problem.simsetup.Nsteps+2)
    #@show thesolution.simdata.indices, ichunk::Int64
    println("Working on:[" * string(t) * "," * string(t) * "] "); flush(stdout)

    evolve!(thesolution, tmpdata, t)
    if (t-2)%thesolution.problem.simsetup.NstepsbetweenNmeas == 0
        measure!(thesolution,t)
    end

    #thesolution.problem.simsetup.lastEvolStep = t # update last EvolStep -> NOw insede Evolve
end

# benchmarking
#using BenchmarkTools
#benchmarkfunction(thesolution)

## store solution - Serialization broke because of pointers in simdata
#println("Storing the solution in: " * thesolution.problem.simsetup.datapath * "/SimSolution.jld"); flush(stdout)
#serialize( thesolution.problem.simsetup.datapath * "/SimSolution.jld", thesolution)
#
### plotting
##include("src/plot2pi.jl")
#plots = Dict()
##plotdata(plots, thesolution, "c")  #create canvas
##plotdata(plots, thesolution, "h", 1, "myplot") # fill canvas
#plotsimdataF(plots, thesolution, "c")  #create canvas
#plotsimdataF(plots, thesolution, "h", 1, "myplot") # fill canvas
##addLOplot(plots, thesolution)
#if isinteractive() == true
#    displayplots(plots)
#else
#    saveplots(plots, thesolution.problem.simsetup.plotpath)
#end

#saveplots(plots, thesolution.problem.simsetup.plotpath)
#displayplots(plots)
#
#plots2 = Dict()
#plotdata_anders(plots2, thesolution, "c")  #create canvas
#plotdata_anders(plots2, thesolution, "h", 1, "T1") # fill canvas
##
#if isinteractive() == true
#    displayplots(plots2)
#else
#    saveplots(plots2, thesolution.problem.simsetup.plotpath)
#end
#
#displayplots(plots2)