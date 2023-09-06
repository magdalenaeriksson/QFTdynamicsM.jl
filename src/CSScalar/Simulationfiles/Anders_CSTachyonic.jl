###################################################################################################################################
using QFTdynamics
using Serialization

###################################################################################################################################
# Set up
###################################################################################################################################
modelfile = CSScalarfile()

parameters = Dict(
    # Discretisation
    "Mass"  =>  0.7, "Nx"    =>  32, "sdim"  =>  3, "dt"    =>  0.05, "Nsteps"=>  300, "Nmeas" =>  0,
    #model
    "Mod"  =>  "CSphi4tachyonicUNS", "Lambda"=>  1, "ONgroup"=> 4,
    #initialisation
    "init"  =>  "Pnr", "xi"     =>  0.4, "sig"     =>  1, "eta"     =>  5, "T"     =>  5, "n"     =>  0, 
    #Reno #PertExp
    "Reno"  =>  "RNone",
    # Numerics
    "Num"   =>  "CPU", "Runs" => 5, "seed" => 0,
    #extra
    "tag"   =>  "AndersTachyonicpaperCS", "ov"    =>  true)
# create Problem

if typeof(modelfile) == CSScalarfile include("../CSScalar_Evolution.jl") end

# create Problem
theproblem = QFTdynamicsProblem(modelfile, parameters)

thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)

initialize!(thesolution, tmpdata)

kickoff!(thesolution, tmpdata)
elapsedtime = @elapsed for t in (thesolution.problem.simsetup.lastEvolStep+1):(thesolution.problem.simsetup.Nsteps+2)
    println("Working on:" * string(t)); flush(stdout)
    evolve!(thesolution, tmpdata, t)
    if (t-2)%thesolution.problem.simsetup.NstepsbetweenNmeas == 0
        measure!(thesolution,tmpdata,t)
    end
end
@show elapsedtime

## store solution - Serialization broke because of pointers in simdata
#println("Storing the solution in: " * thesolution.problem.simsetup.datapath * "/SimSolution.jld"); flush(stdout)
#serialize( thesolution.problem.simsetup.datapath * "/SimSolution.jld", thesolution)

## plotting
plots = Dict()
plotdata(plots, thesolution, "c")  #create canvas
plotdata(plots, thesolution, "h", 1, "myplot") # fill canvas
if isinteractive() == true
    displayplots(plots)
else
    saveplots(plots, thesolution.problem.simsetup.plotpath)
end