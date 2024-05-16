###################################################################################################################################
using QFTdynamics
using Serialization
using ArgParse

###################################################################################################################################
# Set up
###################################################################################################################################
modelfile = TwoPIGaugeScalarfile()
s = ArgParseSettings()
@add_arg_table s begin
    "--Mass"
        arg_type = Float64
        default = 1.0
    "--Nx"
        arg_type = Int
        default = 32
    "--sdim"
        arg_type = Int
        default = 3
    "--dt"
        arg_type = Float64
        default = 0.05
    "--Nsteps"
        arg_type = Int
        default = 100
    "--NstepsinMemory"
        arg_type = Int
        default = 0
    "--Nmeas"
        arg_type = Int
        default = 0
    "--Mod" 
        arg_type = String
        default = "SUNgaugeScalar"
    "--Lambda"
        arg_type = Float64
        default = 0.125
    "--g"
        arg_type = Float64
        default = 1.0
    "--N"
        arg_type = Int
        default = 2
    #initialisation
    "--init" 
        arg_type = String
        default = "Pnr"
    "--T"
        arg_type = Float64
        default = 5.
    "--n"
        arg_type = Float64
        default = 1.
    #Reno
    "--Reno"
        arg_type = String
        default = "RNone"
    "--Pexp" 
        arg_type = String
        default = "LOloop"
    # Numerics
    "--Num"
        arg_type = String
        default = "CPUfull"
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
@show parameters

# create Problem
theproblem = QFTdynamicsProblem(modelfile, parameters)

thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem) 

initialize!(thesolution, tmpdata)

elapsedtime = @elapsed for t in (thesolution.problem.simsetup.lastEvolStep+1):(thesolution.problem.simsetup.Nsteps+2)
   #@show thesolution.simdata.indices
   println("Working on:[" * string(t) * "," * string(t) * "] "); flush(stdout)

   evolve!(thesolution, tmpdata, t)
   if (t-2)%thesolution.problem.simsetup.NstepsbetweenNmeas == 0
       measure!(thesolution, tmpdata, t)
   end
   #@show thesolution.simdata.F[t,t]
   #thesolution.problem.simsetup.lastEvolStep = t # update last EvolStep -> NOw insede Evolve
end

# benchmarking
#using BenchmarkTools
#benchmarkfunction(thesolution)

# store solution - Serialization broke because of pointers in simdata
#println("Storing the solution in: " * thesolution.problem.simsetup.datapath * "/SimSolution.jld"); flush(stdout)
#serialize( thesolution.problem.simsetup.datapath * "/SimSolution.jld", thesolution)

## plotting
#include("src/plot2pi.jl")
plots = Dict()
plotdata(plots, thesolution, "c")  #create canvas
plotdata(plots, thesolution, "h", 1, "myplot") # fill canvas
#addLOplot(plots, thesolution)
if isinteractive() == true
    displayplots(plots)
    #saveplots(plots, thesolution.problem.simsetup.plotpath)
else
    saveplots(plots, thesolution.problem.simsetup.plotpath)
end

#saveplots(plots, thesolution.problem.simsetup.plotpath)
#displayplots(plots)

#plots2 = Dict()
#plotdata_anders(plots2, thesolution, "c")  #create canvas
#plotdata_anders(plots2, thesolution, "h", 1, "T1") # fill canvas
#
#if REPL == true
#    displayplots(plots2)
#else
#    saveplots(plots2, thesolution.problem.simsetup.plotpath)
#end

#displayplots(plots2)