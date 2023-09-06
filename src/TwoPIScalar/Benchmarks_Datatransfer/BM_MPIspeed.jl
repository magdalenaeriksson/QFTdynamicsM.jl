# preparation: use cluster MPI
# julia --project -e 'ENV["JULIA_MPI_BINARY"]="system"; using Pkg; Pkg.build("MPI"; verbose=true)'
# create mpiexecjl and add it to path
# julia --project -e 'using MPI; MPI.install_mpiexecjl()'
# export PATH=$PATH:/home/stud/gerhard/.julia/bin
cd(@__DIR__)
using Pkg
Pkg.activate(".") 
using MPI
MPI.Init()

# print ids
comm = MPI.COMM_WORLD
println("Hello world, I am $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm)) operating on $(gethostname())")
MPI.Barrier(comm)

# send data
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

#N = 4048
N = 125000000 #max number otherwise problem....3d would be better

if rank==1 println("Sending between nodes: ") end
dst = 35 
src = 5 

time = @elapsed begin 
    if rank==src
        send_mesg = Array{Float64}(undef, N)
        fill!(send_mesg, 666)
        sreq = MPI.Isend(send_mesg, dst, rank+(N*8), comm)
        stats = MPI.Waitall!([sreq])
        #print("$rank: Sending   $rank -> $dst = $send_mesg\n")
    end
    
    #size of data
    if rank==dst
        recv_mesg = Array{Float64}(undef, N)
        rreq = MPI.Irecv!(recv_mesg, src,  src+(N*8), comm)
        stats = MPI.Waitall!([rreq])
        #print("$rank: Received $src -> $rank = $recv_mesg\n")
    end
end
if rank==1
    println("Elapsed time: ", time)
    println("datasize: ", Base.format_bytes( N*8) )
    println("Rate: ", Base.format_bytes( N*8 / time) * "/s")
end

MPI.Barrier(comm)

if rank==1 println("Sending inbetween node: ") end
dst = 6 
src = 5 

time = @elapsed begin 
    if rank==src
        send_mesg = Array{Float64}(undef, N)
        fill!(send_mesg, 666)
        sreq = MPI.Isend(send_mesg, dst, rank+(N*8), comm)
        stats = MPI.Waitall!([sreq])
        #print("$rank: Sending   $rank -> $dst = $send_mesg\n")
    end
    
    #size of data
    if rank==dst
        recv_mesg = Array{Float64}(undef, N)
        rreq = MPI.Irecv!(recv_mesg, src,  src+(N*8), comm)
        stats = MPI.Waitall!([rreq])
        #print("$rank: Received $src -> $rank = $recv_mesg\n")
    end
end
if rank==1
    println("Elapsed time: ", time)
    println("datasize: ", Base.format_bytes( N*8) )
    println("Rate: ", Base.format_bytes( N*8 / time) * "/s")
end

MPI.Barrier(comm)

#MPI.Finalize() automatically called by julia
#
