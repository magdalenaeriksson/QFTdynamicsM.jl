###################################################################################################################################
using QFTdynamics
using Serialization

###################################################################################################################################
# Set up
###################################################################################################################################
modelfile = TwoPIScalarfile()

parameters = Dict(
    # Discretisation
    "Mass"  =>  0.7, "Nx"    =>  16, "sdim"  =>  3, "dt"    =>  0.1, "Nsteps"=>  1000, "NstepsinMemory"  =>  400, "Nmeas" =>  0, 
    #model
    "Mod"  =>  "phi4", "Lambda"=>  6, "ONgroup"=> 1,
    #initialisation
    "init"  =>  "TopHatT1", "xi"     =>  0.4, "sig"     =>  1, "eta"     =>  5, "T"     =>  5, "n"     =>  2, 
    #Reno #PertExp
    "Reno"  =>  "RMass", "Pexp"  =>  "LambdaLO",
    # Numerics
    "Num"   =>  "CPUfull", "Nchunks"   =>  0, # if 0 -> set to @Threads.nthreads
    #extra
    "tag"   =>  "Tophatcor", "ov"    =>  true)
# create Problem

#############################################################################
# What to vary
#############################################################################
argument = "Pexp"
values = ["LambdaNLO","LambdaLO"]#, "NinverseNLO"]
#values = ["LambdaNLO"]#, "NinverseNLO"]

for (i,value) in enumerate(values)
    parameters[argument]=value

    global theproblem = QFTdynamicsProblem(modelfile, parameters)
    global thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)

    initialize!(thesolution)

    for t in (thesolution.problem.simsetup.lastEvolStep+1):(thesolution.problem.simsetup.Nsteps+2)
        println("Working on:[" * string(t) * "," * string(t) * "] "); flush(stdout)
        evolve!(thesolution, tmpdata, t)
        if (t-2)%thesolution.problem.simsetup.NstepsbetweenNmeas == 0
            measure!(thesolution,t)
        end
    end
    #println("Storing the solution in: " * thesolution.problem.simsetup.datapath * "/SimSolution.jld"); flush(stdout)
    #serialize( thesolution.problem.simsetup.datapath * "/SimSolution.jld", thesolution)

    ## plotting
    if i==1
        global plots = Dict()
        plotdata(plots, thesolution, "c")  #create canvas
        plotdata_anders(plots, thesolution, "c")
    end

    plotdata(plots, thesolution, "h", i, string(value)) # fill canvas
    plotdata_anders(plots, thesolution, "h", i, string(value)) # fill canvas

    # free memory
    #global plotpath = thesolution.problem.simsetup.plotpath
    #thesolution = 0
end

if isinteractive() == true
    displayplots(plots)
else
    #saveplots(plots, plotpath)
    saveplots(plots, thesolution.problem.simsetup.plotpath)
end
