using QFTdynamics
using BenchmarkTools
using Parameters
using Plots
using CSV
using DataFrames
using ArgParse

# Argparser - which expansion to benchmark
s = ArgParseSettings()
@add_arg_table s begin
    "--Pexp"
        arg_type = String
        default = "LambdaNLO"
end
parameters = parse_args(ARGS, s)

# load benchmark
include("src/TwoPIScalar/TwoPIScalar_Bechmark.jl")
# define memory limit
Memlimit = 200
# tag for files to write
tag = "200GB"

# lets go! Which expansion to benchmark "LambdaNLO" or "NinverseNLO"
PertExpansion = parameters["Pexp"]

# rule of thumb: benchmark time should be 7 times a single execution
#timeforbenchmark = 10. #seconds
timeforbenchmark = 0. # --> run the simple benchmark 

modelfile = TwoPIScalarfile()
# Dictionary
parameters = Dict(  "Mass"  =>  0.7, "sdim"  =>  3, "dt"    =>  0.1, "Nsteps"=>  1000, # Discretisation
                    "Nmeas" =>  0, # either 0 (all timesteps are measured) or >=3 (because the plot script wants to draw 3 graphs)
                    "Mod"  =>  "phi4", "Lambda"=>  6, "ONgroup"=> 1, #model
                    "init"  =>  "TopHatT1", "T"     =>  5, "n"     =>  2, #initialisation
                    "Reno"  =>  "RNone", "Pexp"  =>  PertExpansion, #Reno end Expansion type
                    "tag"   =>  "benchmark", "ov"    =>  true, #extr#
                    "Num"   =>  "CPUfull",
                    "Nchunks"   =>  0,
                    "NstepsinMemory"  =>  3,
                    "Nx"    =>  2,)

# Numrange
Numrange = ["CPUfull", "CPUred"]
# Nxrange
#Nxrange = [16,24,]
Nxrange = [16,32,42,54,64,112,128,144,162,180]
# NstepsinMemoryrange
NstepsinMemoryrange = [400,600]
#NstepsinMemoryrange = [50,100]

df = DataFrame()
for Num in Numrange 
    for Nx in Nxrange
        for NstepsinMemory in NstepsinMemoryrange
            println(" *** Working on Nx=", Nx, ", NstepsinMemory=", NstepsinMemory, ", Num=" * Num)
            # update parameters
            parameters["Num"]=Num
            parameters["Nx"]=Nx
            parameters["NstepsinMemory"]=NstepsinMemory
            # set up the problem and "fill up the kernel" (just changing the index range)
            approxmemgb = approximateMemoryusage( parameters, Threads.nthreads())
            @show approxmemgb 
            if approxmemgb < Memlimit
                (time, memgb) = TwoPIScalar_bechmarkevolution(modelfile, parameters, timeforbenchmark)
                GC.gc() # run GC collector 
                tmpdf = DataFrame( Nx = [Nx], NstepsinMemory=[NstepsinMemory], threads=[Threads.nthreads()],
                               time = [time], Mem=[memgb],Num=[Num])
                append!(df, tmpdf)
            end
        end    
    end
end
#
# Save the df
#
lines = readlines("etc/paths.txt")
plotpath = chop(lines[4],head=5,tail=0)
dfname = "2PIScalarBenchmark_" * PertExpansion * "_nThreads" * string(Threads.nthreads()) * "_" * tag * ".csv"
CSV.write(joinpath(plotpath * "/Benchmark", dfname), df)
println("Dataframe saved as " * joinpath(plotpath * "/Benchmark", dfname) ); flush(stdout)

#
# Create Plots
#
# computation time 
cplot=plot(xlabel="Nx", ylabel="time / step in seconds", legend=:topleft, title="2PI phi4 in 3d, " * PertExpansion * ", nthreads=" * string(Threads.nthreads()) )
dfCPUfull = filter(:Num => ==("CPUfull"), df)
dfCPUred = filter(:Num => ==("CPUred"), df)

label = ["CPUfull", "CPUred"]
color = [:blue,:green]

nNsteps = length(NstepsinMemoryrange)
ls = [:solid,:dash,:dot]

for (j, dfs) in enumerate([dfCPUfull,dfCPUred])
    for (i, NstepsinMemory) in enumerate(NstepsinMemoryrange)
        tmpdf = filter(:NstepsinMemory => ==(NstepsinMemory), dfs)
        plot!(cplot, tmpdf.Nx, tmpdf.time, label= label[j] * " Kernel=" * string(NstepsinMemory), line=(ls[i],2), color=color[j], marker=:x)
    end
end
plotname = "2PIScalarSimpleBenchmark_" * PertExpansion * "_nThreads" * string(Threads.nthreads()) * "_" * tag * ".png"
savefig( cplot, joinpath(plotpath * "/Benchmark", plotname) )
println("Saved as " * joinpath(plotpath, plotname) ); flush(stdout)
chmod(plotpath * "/Benchmark", 0o777, recursive=true)
