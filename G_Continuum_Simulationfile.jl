###################################################################################################################################
using QFTdynamics
using Serialization
using ArgParse

###################################################################################################################################
# Set up
###################################################################################################################################
# get parameters as dict from parameter file and specify modelfile (a type)
modelfile = TwoPIScalarfile()

s = ArgParseSettings()
@add_arg_table s begin
    "--Mass"
        arg_type = Float64
        default = 0.7
    "--Nmom"
        arg_type = Int
        default = 200
    "--Nx"
        arg_type = Int
        default = 16
    "--sdim"
        arg_type = Int
        default = 3
    "--dt"
        arg_type = Float64
        default = 0.1
    "--Nsteps"
        arg_type = Int
        default = 100
    "--NstepsinMemory"# 3 - min for LO
        arg_type = Int
        default = 0
    "--Nmeas"
        arg_type = Int
        default = 0
    "--Mod" 
        arg_type = String
        default = "phi4"
    "--Lambda"
        arg_type = Float64
        default = 5.
    "--ONgroup"
        arg_type = Int
        default = 4
    #initialisation
    "--init" 
        arg_type = String
        default = "TopHatT1"
    "--xi"
        arg_type = Float64
        default = 1.
    "--eta"
        arg_type = Float64
        default = 1.
    "--sig"
        arg_type = Float64
        default = 1.
    "--T"
        arg_type = Float64
        default = 1.
    "--n"
        arg_type = Float64
        default = 1.
    #Reno
    "--Reno"
        arg_type = String
        default = "RNone"
    #PertExp
    "--Pexp"
        arg_type = String
        default = "NinverseNLO"
    # Numerics
    "--Num"
        arg_type = String
        default = "CPUred" #"CPUcont"
    "--Nchunks"
        arg_type = Int
        default = 0
    #extra
    "--tag"
        arg_type = String
        default = "test"
    "--ov"
        action = :store_true
        default = true
end
parameters = parse_args(ARGS, s)
parameters["ov"] = true
#parameters["Pexp"] = "NinverseLO"
#parameters["Pexp"] = "LambdaLO"
parameters["Pexp"] = "LambdaNLO"
#parameters["Num"] = "CPUred2"
parameters["Num"] = "CPUred"

# create Problem
theproblem = QFTdynamicsProblem(modelfile, parameters)

thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)

initialize!(thesolution, tmpdata)
#tmpdata = TwoPIScalarTmpDataCPUfull(thesolution.simdata, thesolution.problem.num.nchunks)

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

# store solution - Serialization broke because of pointers in simdata
#println("Storing the solution in: " * thesolution.problem.simsetup.datapath * "/SimSolution.jld"); flush(stdout)
#serialize( thesolution.problem.simsetup.datapath * "/SimSolution.jld", thesolution)
println("Storing the measurearray in: " * thesolution.problem.simsetup.datapath * "/Measurearray.jld"); flush(stdout)
serialize( thesolution.problem.simsetup.datapath * "/Measurearray.jld", thesolution.measurearray)
println("Storing the problem in: " * thesolution.problem.simsetup.datapath * "/Problem.jld"); flush(stdout)
serialize( thesolution.problem.simsetup.datapath * "/Problem.jld", thesolution.problem)

## plotting
#include("src/plot2pi.jl")
plots = Dict()
plotdata(plots, thesolution, "c")  #create canvas
plotdata(plots, thesolution, "h", 1, "myplot") # fill canvas
#addLOplot(plots, thesolution)
if isinteractive() == true
    displayplots(plots)
else
    saveplots(plots, thesolution.problem.simsetup.plotpath)
end

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

#reducedlattice = thesolution.simdata.F[49,20]
#fulllattice = tmpdata.F[1]
#
#copytolattice!(     fulllattice, reducedlattice, thesolution.problem.disc.fftwhelper)
#@show fulllattice
#copytofulllattice!( fulllattice, reducedlattice, tmpdata)
#@show fulllattice
#
#copytoreducedlattice!( reducedlattice, fulllattice, thesolution.problem.disc.fftwhelper)
#copytoreducedlattice!( reducedlattice, fulllattice, tmpdata)