using Plots
using LaTeXStrings

Nmin=30
Nmax=160

GPUMemt(N,d,GPUMem) = GPUMem/(2*(N^d+2*(N/2+1)*N^(d-1))*8*10^(-9)) # F,rho and tmp storage
CPUMemt(N,d,CPUMem) = sqrt( (48*CPUMem)/(N^d*8*10^(-9)) ) # reduced lattice for F, rho
a(N,d) = N^d*8*10^(-9)
b(N,d) = (2*(N^d+2*(N/2+1)*N^(d-1))*8*10^(-9))
onlyGPUMemt(N,d,Mem) = (-b(N,d)+sqrt(b(N,d)^2+4*Mem*a(N,d)) )/(2*a(N,d)) 

function optimalNvalues(Nmin, Nmax)
    N = zeros(11*7*5*4)
    i = 1
    for a=0:10, b=0:6, c=0:4, d=0:3
        N[i] = 2^a*3^b*5^c*7^d ## FFTW works best if N can be factorized like that
        i+=1
    end 
    sort!(N)
    idxNmin = findmin(abs.(N.-Nmin))[2]
    idxNmax = findmin(abs.(N.-Nmax))[2]
    return N[idxNmin:idxNmax]
end

N = optimalNvalues(Nmin, Nmax)
Nticks = [32,64,128,256,512]
kticks=[500,1000,5000,10000,50000,100000]


p = plot(xlabel = "Nx", ylabel = "K", title="2PI Scalarfield in 3 dimensions - Memory constraint", xaxis=:log, yaxis=:log, xticks=(Nticks,Nticks), yticks=(kticks,kticks), ylims=[300,100000])

# GPU
plot!(p,N,map(x->GPUMemt(x,3,  40),N), label="1xGPU (40GB)", color="green",  line=(:solid, 3) )
plot!(p,N,map(x->GPUMemt(x,3,2*40),N), label="2xGPU (a 40GB)", color="green",  line=(:dash,  3))
#plot!(p,N,map(x->GPUMemt(x,3,3*40),N), label="3xGPU")
plot!(p,N,map(x->GPUMemt(x,3,4*40),N), label="4xGPU (a 40 GB)", color="green",  line=(:dot,   3)   )
# CPU
plot!(p,N,map(x->CPUMemt(x,3,1*256),N), label="1xCPU (256 GB)", color="blue",  line=(:solid, 3) )
plot!(p,N,map(x->CPUMemt(x,3,2*256),N), label="2xCPU (a 256 GB)", color="blue",  line=(:dash,  3))
#plot!(p,N,map(x->CPUMemt(x,3,3*256),N), label="3xCPU")")
plot!(p,N,map(x->CPUMemt(x,3,4*256),N), label="4xCPU (a 256 GB)", color="blue",  line=(:dot,   3)   )
# entire simulation on GPU
plot!(p,N,map(x->onlyGPUMemt(x,3,40),N), label="only 1xGPU (40 GB)", color="black",  line=(:solid,   3)   )
plot!(p,N,map(x->onlyGPUMemt(x,3,80),N), label="only 1xGPU (80 GB)", color="black",  line=(:dot,   3)   )
#show
#display(p)
lines = readlines("../../etc/paths.txt")
plotpath = chop(lines[4],head=5,tail=0)
plotname = "Benchmark/2PIScalar_Memoryconstraints.png"
savefig( p, joinpath(plotpath, plotname) )
println("Saved as " * joinpath(plotpath, plotname) ); flush(stdout)
#savefig(p,"2PIScalar_Memory.png")