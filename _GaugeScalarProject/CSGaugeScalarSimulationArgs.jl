###################################################################################################################################
using QFTdynamics
using Serialization
using ArgParse
using Plots # for additional plot
using LaTeXStrings  # for additional plot
###################################################################################################################################
# Set up
###################################################################################################################################
modelfile = CSGaugeScalarFile()

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
        default = 300
    "--Nmeas"
        arg_type = Int
        default = 0
    "--Mod" 
        arg_type = String
        default = "CS_SUNgaugeScalar"
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
        default = "CS"
    # Numerics
    "--Num"
        arg_type = String
        default = "CPU"
    "--Runs"
        arg_type = Int
        default = 8
    "--seed"
        arg_type = Int
        default = 1
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

## create problem
theproblem = QFTdynamicsProblem(modelfile, parameters)
## create solution
thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)

println("...Initializing..."), flush(stdout)
initialize!(thesolution, tmpdata)

elapsedtime = @elapsed for t in (thesolution.problem.simsetup.lastEvolStep+1):(thesolution.problem.simsetup.Nsteps+2)
    println("Working on timestep: " * string(t)); flush(stdout)

    evolve!(thesolution, tmpdata, t)
    #evolveSerial!(thesolution, tmpdata, t)

    if (t-2)%thesolution.problem.simsetup.NstepsbetweenNmeas == 0
        measure!(thesolution, tmpdata, t)
    end
end

println("Simulation time in h: ", elapsedtime/3600); flush(stdout)
# store solution - Serialization broke because of pointers in simdata
#println("Storing the solution in: " * thesolution.problem.simsetup.datapath * "/SimSolution.jld"); flush(stdout)
#serialize( thesolution.problem.simsetup.datapath * "/SimSolution.jld", thesolution)

## plotting
plots = Dict()
plotdata(plots, thesolution, "c")  #create canvas
plotdata(plots, thesolution, "h", 1, "myplot") # fill canvas
# plot scalar compontents
plotScalarcomponentdata(plots, thesolution, 1, "myplot") # fill canvas

# additional plot (phi without zero mode)
#
# Phi2 snapshots
#
timeidx = [1, Int64(floor( length(thesolution.measurearray)/2) ), length(thesolution.measurearray)]
kvals = [el.lev for el in thesolution.problem.disc.fftwhelper]
plots["Phi2_nozero.png"] = plot(xlabel = L"k_L", ylabel =L"<|\Phi(t;\mathbf{p}=0)|^2>/V", xlim=(kvals[1],kvals[end]))
for (i, tidx) in enumerate(timeidx)
    plot!( plots["Phi2_nozero.png"],kvals[2:end], thesolution.measurearray[tidx].Phi2k[2:end], label = "tm=" * string(round(thesolution.measurearray[tidx].time,digits=1)))
end
#
# Phi2 timeevolution
#
nkLmomenta = 3
kLmax = thesolution.problem.disc.fftwhelper[end].lev
kLidx = zeros(Int64, nkLmomenta)
for i in 0:(length(kLidx)-1)
    kLidx[i+1] = Int64(findmin(abs.([element.lev for element in thesolution.problem.disc.fftwhelper] .- ( (i/(nkLmomenta-1))*kLmax ) ))[2])
end
kLidx[1] = 2 # the fist entry should be the lowest IR mode, not the 0 mode
plots["Phi2_nozero_t.png"] = plot(xlabel = L"tm",ylabel =L"<|\Phi(t;\mathbf{p}=0)|^2>/V")
for (i, kidx) in enumerate(kLidx)
    plot!(plots["Phi2_nozero_t.png"],[element.time for element in thesolution.measurearray], [element.Phi2k[kidx] for element in thesolution.measurearray], label="k=" * string(round(thesolution.problem.disc.fftwhelper[kidx].lev, digits=1)))
end
 
if isinteractive() == true
    displayplots(plots)
    saveplots(plots, thesolution.problem.simsetup.plotpath)
else
    saveplots(plots, thesolution.problem.simsetup.plotpath)
end
