###################################################################################################################################
using QFTdynamics
using Serialization

###################################################################################################################################
# Set up
###################################################################################################################################
modelfile = TwoPIScalarfile()

parameters = Dict(
    # Discretisation
    "Mass"  =>  0.5, "Nx"    =>  16, "sdim"  =>  2, "dt"    =>  0.1, "Nsteps"=>  100, "NstepsinMemory"  =>  0, "Nmeas" =>  0, 
    #model
    "Mod"  =>  "phi4", "Lambda"=>  10, "ONgroup"=> 4,
    #initialisation
    "init"  =>  "TopHatT1", "T"     =>  5, "n"     =>  2,
    #Reno #PertExp
    "Reno"  =>  "RNone", "Pexp"  =>  "LambdaLO",
    # Numerics
    "Num"   =>  "CPUfull", "Nchunks"   =>  0, # if 0 -> set to @Threads.nthreads
    #extra
    "tag"   =>  "testing", "ov"    =>  true)
# create Problem

#############################################################################
# What to vary
#############################################################################
argument = "Pexp"
values = ["NinverseLO","NinverseNLO"]#, "NinverseNLO"]

#argument = "Nx"
#values = [32,128]#, "NinverseNLO"]

for (i,value) in enumerate(values)
    parameters[argument]=value

    global theproblem = QFTdynamicsProblem(modelfile, parameters)
    global thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)

    initialize!(thesolution, tmpdata)

    for t in (thesolution.problem.simsetup.lastEvolStep+1):(thesolution.problem.simsetup.Nsteps+2)
        println("Working on:[" * string(t) * "," * string(t) * "] "); flush(stdout)
        evolve!(thesolution, tmpdata, t)
        if (t-2)%thesolution.problem.simsetup.NstepsbetweenNmeas == 0
            measure!(thesolution,t)
        end
    end
    #println("Storing the solution in: " * thesolution.problem.simsetup.datapath * "/SimSolution.jld"); flush(stdout)
    #serialize( thesolution.problem.simsetup.datapath * "/SimSolution.jld", thesolution.measurearray)
    #corr = (el.corr(0) for el in thesolution.measurearray)


    ## plotting
    if i==1
        global plots = Dict()
        plotdata(plots, thesolution, "c")  #create canvas
    end

    plotdata(plots, thesolution, "h", i, string(value)) # fill canvas
    if isinteractive() == true
        displayplots(plots)
    else
        saveplots(plots, thesolution.problem.simsetup.plotpath)
    end
    saveplots(plots, thesolution.problem.simsetup.plotpath)
end
