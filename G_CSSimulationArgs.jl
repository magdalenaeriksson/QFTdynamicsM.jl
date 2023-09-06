###################################################################################################################################
#julia -t 10 --project=. G_CSSimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.01 --Nsteps 50 --Nmeas 0 --Mod CSphi4 --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --Reno RNone --Pexp NinverseNLO --Num CPU --Runs 4 --tag tere --ov
###################################################################################################################################
using QFTdynamics
using Serialization
using ArgParse

###################################################################################################################################
# Set up
###################################################################################################################################
# get parameters as dict from parameter file and specify modelfile (a type)
modelfile = CSScalarfile()

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
        default = 3
    "--dt"
        arg_type = Float64
        default = 0.05
    "--Nsteps"
        arg_type = Int
        default = 200
    "--Nmeas"
        arg_type = Int
        default = 0
    "--Mod" 
        arg_type = String
        default = "CSphi4"
    "--Lambda"
        arg_type = Float64
        default = 5.
    "--ONgroup"
        arg_type = Int
        default = 4
    #initialisation
    "--init" 
        arg_type = String
        default = "Pnr"
    "--T"
        arg_type = Float64
        default = 1.
    "--n"
        arg_type = Float64
        default = 0.
    #Reno
    "--Reno"
        arg_type = String
        default = "RNone"
    # Numerics
    "--Num"
        arg_type = String
        default = "CPU"
    "--Runs"
        arg_type = Int
        default = 2
    "--seed"
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

if typeof(modelfile) == CSScalarfile include("src/CSScalar/CSScalar_Evolution.jl") end

# create Problem
theproblem = QFTdynamicsProblem(modelfile, parameters)

thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)

initialize!(thesolution, tmpdata)

kickoff!(thesolution, tmpdata)

elapsedtime = @elapsed for t in (thesolution.problem.simsetup.lastEvolStep+1):(thesolution.problem.simsetup.Nsteps+2)
    if t%100 == 0 println("Working on:" * string(t)); flush(stdout) end

    evolve!(thesolution, tmpdata, t)
    if (t-2)%thesolution.problem.simsetup.NstepsbetweenNmeas == 0
        measure!(thesolution,tmpdata,t)
    end
end
@show elapsedtime

# store solution - Serialization broke because of pointers in simdata
#println("Storing the solution in: " * thesolution.problem.simsetup.datapath * "/SimSolution.jld"); flush(stdout)
#serialize( thesolution.problem.simsetup.datapath * "/SimSolution.jld", thesolution)

# plotting
plots = Dict()
plotdata(plots, thesolution, "c")  #create canvas
plotdata(plots, thesolution, "h", 1, "data") # fill canvas

if isinteractive() == true
    displayplots(plots)
else
    saveplots(plots, thesolution.problem.simsetup.plotpath)
end