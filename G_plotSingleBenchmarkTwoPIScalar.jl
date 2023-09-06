# useage
# julia --project=. G_plotSingleBenchmarkTwoPIScalar.jl --file 2PIScalarBenchmark_LambdaNLO_nThreads28_benchmark.csv
using Plots
using CSV
using DataFrames
using ArgParse
using StatsBase

s = ArgParseSettings()
@add_arg_table s begin
    "--file"
        arg_type = String
end
parameters = parse_args(ARGS, s)

filename = parameters["file"]
filename = "2PIScalarBenchmark_LambdaNLO_nThreads28_benchmark.csv"
println("Working on ", filename); flush(stdout)

lines = readlines("etc/paths.txt")
plotpath = chop(lines[4],head=5,tail=0)

#
# read df
#
println("Reading df..."); flush(stdout)
df = CSV.read(joinpath(plotpath * "/Benchmark", filename), DataFrame; header=["Nx","NstepsinMemory", "threads", "time", "Mem", "Num", "Pexp",])

#
# create average Dataframe
#

#
# plot
#
cplot=plot(xlabel="Nx", ylabel="time / step in seconds", legend=:topleft, title="2PI phi4 in 3d")
dfCPUfull = filter(:Num => ==("CPUfull"), df)
dfCPUred = filter(:Num => ==("CPUred"), df)

label = ["CPUfull", "CPUred"]
color = [:blue,:green]

NstepsinMemoryrange = unique(df.NstepsinMemory) 
println("NstepsinMemory values:", length( NstepsinMemoryrange )); flush(stdout)

ls = [:solid,:dash,:dot]

for (j, dfs) in enumerate([dfCPUfull,dfCPUred])
    for (i, NstepsinMemory) in enumerate(NstepsinMemoryrange)
        tmpdf = filter(:NstepsinMemory => ==(NstepsinMemory), dfs)
        #plot!(cplot, tmpdf.Nx, tmpdf.time, label= label[j] * " Kernel=" * string(NstepsinMemory), line=(ls[i],2), color=color[j], marker=:x)
        # calx mean and std for given Nx
        gdf = groupby(tmpdf, Not(r"time") )
        gdf2 = combine(gdf, :time => mean, :time => sem)
        sort!(gdf2,:Nx)
        #plot!(cplot, gdf2.Nx, gdf2.time_mean, yerr=gdf2.time_sem, label= label[j] * " Kernel=" * string(NstepsinMemory), line=(ls[i],2), color=color[j], marker=:x)
        plot!(cplot, gdf2.Nx, gdf2.time_mean, yerr=gdf2.time_sem, label= label[j] * " Kernel=" * string(NstepsinMemory),)
    end
end
plotname = filename[1:end-4] * ".pdf"
savefig( cplot, joinpath(plotpath * "/Benchmark", plotname) )
println("Saved as " * joinpath(plotpath, plotname) ); flush(stdout)
chmod(plotpath * "/Benchmark", 0o777, recursive=true)
