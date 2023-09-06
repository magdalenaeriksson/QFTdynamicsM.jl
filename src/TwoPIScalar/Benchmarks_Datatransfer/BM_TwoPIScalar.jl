using QFTdynamics
using BenchmarkTools
using Parameters
using Plots
using CSV
using DataFrames

# rule of thumb: benchmark time should be 7 times a single execution
timeforbenchmark = 10 #seconds

modelfile = TwoPIScalarfile()
# Dictionary
parameters = Dict(  "Mass"  =>  0.7, "sdim"  =>  3, "dt"    =>  0.1, "Nsteps"=>  100, # Discretisation
                    "Nmeas" =>  0, # either 0 (all timesteps are measured) or >=3 (because the plot script wants to draw 3 graphs)
                    "Mod"  =>  "phi4", "Lambda"=>  6, "ONgroup"=> 1, #model
                    "init"  =>  "TopHatT1", "T"     =>  5, "n"     =>  2, #initialisation
                    "Reno"  =>  "RNone", "Pexp"  =>  "LambdaNLO", #Reno end Expansion type
                    "tag"   =>  "benchmark", "ov"    =>  true, #extr#
                    "Nchunks"   =>  0,
################################# what matters ################################
                    "Num"   =>  "GPUreduced",
                    "NstepsinMemory"  =>  3,
                    "Nx"    =>  2,)

### loop through
# Numrange
Numrange = ["CPUfull", "CPUred"]
# Nxrange
Nxrange = [16,24,]
#Nmin = 16
#Nmax = 100
#NrofNvalues = 5
#function optimalNvalues(Nmin, Nmax)
#    N = zeros(11*7*5*4)
#    i = 1
#    for a=0:10, b=0:6, c=0:4, d=0:3
#        N[i] = 2^a*3^b*5^c*7^d ## FFTW works best if N can be factorized like that
#        i+=1
#    end 
#    sort!(N)
#    idxNmin = findmin(abs.(N.-Nmin))[2]
#    idxNmax = findmin(abs.(N.-Nmax))[2]
#    return N[idxNmin:idxNmax]
#end
#optimalN = optimalNvalues(Nmin, Nmax)
#indicesofNvalues = 1:Int(round(length(optimalN)/NrofNvalues)):length(optimalN)
#Nxrange = Int.(optimalN[indicesofNvalues])

#NstepsinMemory
#NstepsinMemoryrange = [400,800,1200]
NstepsinMemoryrange = [100,200]

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
            theproblem = QFTdynamicsProblem(modelfile, parameters)
            thesolution = QFTdynamicsSolution(modelfile, theproblem)

            #PrintMemoryofSimData!(simdata)
            for i in 1:thesolution.simdata.NstepsinMemory
                expandSimData!(thesolution.simdata) # now indices range from 1:NstepsinMemory
            end
            #@show simdata.indices

            #@benchmark evolve!($thesolution,t) samples= seconds= evals=
            #@benchmark evolve!($thesolution,$thesolution.simdata.indices[2]) seconds=10
            time = @belapsed evolve!($thesolution,$thesolution.simdata.indices[2]) seconds=timeforbenchmark

            tmpdf = DataFrame( Nx = [Nx], NstepsinMemory=[NstepsinMemory], threads=[Threads.nthreads()],
                               time =   [time], Num=[Num])
            append!(df, tmpdf)
        end    
    end
end

# plot
theplot=plot(xlabel="Nx", ylabel="time / step in seconds", title="2PI phi4/lambda in 3d, nthreads=" * string(Threads.nthreads()) )
dfCPU = filter(:Num => ==("CPUfull"), df)
dfGPU = filter(:Num => ==("GPUreduced"), df)

label = ["CPU", "GPU"]
color = [:blue,:green]

nNsteps = length(NstepsinMemoryrange)
ls = [:solid,:dash,:dot]

for (j, dfs) in enumerate([dfCPU,dfGPU])
    for (i, NstepsinMemory) in enumerate(NstepsinMemoryrange)
        tmpdf = filter(:NstepsinMemory => ==(NstepsinMemory), dfs)
        plot!(theplot, tmpdf.Nx, tmpdf.time, label= label[j] * " Kernel=" * string(NstepsinMemory), line=(ls[i],2), color=color[j], marker=:x)
    end
end
#display(theplot)
lines = readlines("etc/paths.txt")
plotpath = chop(lines[4],head=5,tail=0)
plotname = "Benchmark/2PIScalar_Benchmark_nThreads" * string(Threads.nthreads()) * ".png"
savefig( theplot, joinpath(plotpath, plotname) )
println("Saved as " * joinpath(plotpath, plotname) ); flush(stdout)