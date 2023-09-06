using CUDA, BenchmarkTools

#GPUMemorylimit = 5.5 #in GiB
GPUMemorylimit = 35 #in GiB
mytype = Float64
Nx=128

Memtlist = []
timelist = []

datasizelattice = Nx^3*8*10^(-9)
#datasize(Memt) = datasizelattice*Memt
nrMemtinGPUlimit = GPUMemorylimit/datasizelattice
println("set GPU limit: ", GPUMemorylimit," Maximum nr of Memt ", nrMemtinGPUlimit, " for a 3dim lattice with Nx: ", Nx ); flush(stdout)

#Memtlist = range(1,Int(round(nrMemtinGPUlimit)),step=75)
#Memtlist = range(1,Int(round(nrMemtinGPUlimit)),step=75)
Memtlist = Int.(round.(Array(LinRange(1,Int(round(nrMemtinGPUlimit)), 6))))
@show Memtlist
timelist = zero(Memtlist)

for (i,Memt) in enumerate(Memtlist)
    # first version
    global data = zeros(mytype, Nx,Nx,Nx,Memt)
    #println("before alloc. "); flush(stdout)
    #@show CUDA.memory_status()
    gpu_data = CuArray(data);
    #println("after alloc. "); flush(stdout)
    #@show CUDA.memory_status()
    time = @belapsed copyto!($gpu_data, $data)
    CUDA.unsafe_free!(gpu_data) # free the momory
    CUDA.reclaim() # reclaims the "unused" memory
    #println("after free. "); flush(stdout)
    #@show CUDA.memory_status()
    println("1. Memt:", Memt, " / Size:", Base.format_bytes(sizeof(data)), " / Speed:", Base.format_bytes(sizeof(data) / time) * "/s"); flush(stdout)
    # second version
    #println("before alloc. "); flush(stdout)
    #@show CUDA.memory_status()
    gpu = Mem.alloc(Mem.Device, sizeof(data)) #gpu magic
    gpu_ptr = convert(CuPtr{mytype}, gpu)  #gpu magic
    #println("after alloc. "); flush(stdout)
    #@show CUDA.memory_status()
    cpu = Mem.alloc(Mem.Host, length(data)*sizeof(mytype)) #cpu magic
    cpu_ptr = convert(Ptr{mytype}, cpu) # cpu magic
    time = @belapsed unsafe_copyto!($gpu_ptr, $cpu_ptr, length(data))
    println("2. Memt:", Memt, " / Size:", Base.format_bytes(sizeof(data)), " / Speed:", Base.format_bytes(sizeof(data) / time) * "/s"); flush(stdout)
    Mem.free(cpu)
    Mem.free(gpu)
    CUDA.reclaim() # reclaims the "unused" memory
    #println("after free. "); flush(stdout)
    #@show CUDA.memory_status()
end
