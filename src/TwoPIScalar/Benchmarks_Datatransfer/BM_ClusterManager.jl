cd(@__DIR__)
using Pkg
Pkg.activate(".") 
using Distributed
using ClusterManagers

addprocs(SlurmManager(parse(Int, ENV["SLURM_NTASKS"])-1); exeflags="--project")
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere @show Base.load_path()
@everywhere using BenchmarkTools

# Hello!
for i in procs()
    host, pid, id = fetch(@spawnat i (gethostname(), getpid(), myid()))
    println("Hello from  ", id , "/", pid, " located at ", host)
end

## copying data
## master (id=1) asks processor "theid" to create a 2x2 matrix
#theid=5
#r = remotecall(rand, theid, 2, 2) # this returns a Future
## tell the master to wait until "theid" completet its task
#wait(r)
## thell the id to add 1 to the random matrix 
#s = @spawnat theid 1 .+ fetch(r)
## now i want to have the result
#fetch(s)
#@show fetch(s)


# 
#Nx = 128
#Memt = 100
Nx = 200
Memt = 100
println("datasize: ", Base.format_bytes( Nx^3*Memt*8) )

println("****** @belapsed *******")
theid1 = 5
println("1. starting remote call:")
r1 = remotecall(rand, theid1, Nx, Nx, Nx, Memt) # this returns a Future
println("waiting for it:")
wait(r1)
println("fetching data:")
time = @belapsed fetch(r1)
println("id: ", theid1, " / rate: ", Base.format_bytes( Nx^3*Memt*8 / time) * "/s")

theid2 = 35
println("2. starting remote call:")
r2 = remotecall(rand, theid2, Nx, Nx, Nx, Memt) # this returns a Future
println("waiting for it:")
wait(r2)
println("fetching data:")
time = @belapsed fetch(r2)
println("id: ", theid2, " / rate: ", Base.format_bytes( Nx^3*Memt*8 / time) * "/s")

println("****** @elapsed *******")
theid3 = 5
println("3. starting remote call:")
r3 = remotecall(rand, theid3, Nx, Nx, Nx, Memt) # this returns a Future
println("waiting for it:")
wait(r3)
println("fetching data:")
time = @elapsed fetch(r3)
println("id: ", theid3, " / rate: ", Base.format_bytes( Nx^3*Memt*8 / time) * "/s")

theid4 = 35
println("4. starting remote call:")
r4 = remotecall(rand, theid4, Nx, Nx, Nx, Memt) # this returns a Future
println("waiting for it:")
wait(r4)
println("fetching data:")
time = @elapsed fetch(r4)
println("id: ", theid4, " / rate: ", Base.format_bytes( Nx^3*Memt*8 / time) * "/s")



######################### SOME EXTRA STUFF #####################################

# 1)EXTRA: something i found online...could be useful at some point
#addprocs(SlurmManager(512), N=16, topology=:master_worker, exeflags="--project=.")
#pmap(x -> println(x), 1:10)
#
# 2)EXTRA: something i found online...could be useful at some point
#cd(@__DIR__)
#using Pkg
#Pkg.activate(".")
#using Distributed, ClusterManagers
#subs = Dict("x"=>"*", "(" => "", ")" => "");
#np = sum(eval(Meta.parse(replace(ENV["SLURM_JOB_CPUS_PER_NODE"], r"x|\(|\)" => s -> subs[s]))))
#addprocs(SlurmManager(np); exeflags="--project")
