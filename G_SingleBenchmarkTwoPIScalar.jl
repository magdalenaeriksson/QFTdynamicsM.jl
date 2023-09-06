using QFTdynamics
using BenchmarkTools
using Parameters
using Plots
using CSV
using DataFrames
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
    "--Nx"
        arg_type = Int
        default = 16
    "--sdim"
        arg_type = Int
        default = 2
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
        default = 1.
    "--ONgroup"
        arg_type = Int
        default = 1
    #initialisation
    "--init" 
        arg_type = String
        default = "Thermal"
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
        default = "RMass"
    #PertExp
    "--Pexp"
        arg_type = String
        default = "LambdaLO"
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
end
parameters = parse_args(ARGS, s)

# load benchmark
include("src/TwoPIScalar/TwoPIScalar_Bechmark.jl")

Memlimit = 245
# rule of thumb: benchmark time should be 7 times a single execution
timeforbenchmark = 0. #seconds

approxmemgb = approximateMemoryusage( parameters, Threads.nthreads())
println("Approximate memory usage / GB :",approxmemgb)
if approxmemgb < Memlimit
    (time, memgb) = TwoPIScalar_bechmarkevolution(modelfile, parameters, timeforbenchmark)
    println("Actual memory usage / GB :", memgb)
    println("Time per time step / sec:", time)

    # save to DataFrame
    df = DataFrame( Nx = [ parameters["Nx"] ], NstepsinMemory=[ parameters["NstepsinMemory"] ], threads=[Threads.nthreads()],  time = [time], Mem=[memgb],Num=[ parameters["Num"] ], Pexp = [ parameters["Pexp"] ] )

    #
    # Save the df
    #
    lines = readlines("etc/paths.txt")
    plotpath = chop(lines[4],head=5,tail=0)
    dfname = "2PIScalarBenchmark_" * parameters["Pexp"] * "_nThreads" * string(Threads.nthreads()) * "_" * parameters["tag"] * ".csv"
    CSV.write(joinpath(plotpath * "/Benchmark", dfname), df, append=true)
    println("Dataframe saved as " * joinpath(plotpath * "/Benchmark", dfname) ); flush(stdout)
    chmod(plotpath * "/Benchmark", 0o777, recursive=true)
end
