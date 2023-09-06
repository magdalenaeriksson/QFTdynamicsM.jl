###################################################################################################################################
using QFTdynamics
using Serialization

###################################################################################################################################
# Set up
###################################################################################################################################
modelfile = TwoPIScalarfile()

parameters = Dict(
    # Discretisation
    "Mass"  =>  1.0, "Nx"    =>  1, "sdim"  =>  1, "dt"    =>  0.01, "Nsteps"=>  600, "NstepsinMemory"  =>  0, "Nmeas" =>  0, 
    #model
    "Mod"  =>  "phi4", "Lambda"=>  1, "ONgroup"=> 8,
    #initialisation
    "init"  =>  "Gauss", "xi"     =>  0.4, "sig"     =>  1, "eta"     =>  5, "T"     =>  5, "n"     =>  2, 
    #Reno #PertExp
    "Reno"  =>  "RNone", "Pexp"  =>  "LambdaLO",
    # Numerics
    "Num"   =>  "CPUfull", "Nchunks"   =>  0, # if 0 -> set to @Threads.nthreads
    #extra
    "tag"   =>  "AndersNNLOpaperALL", "ov"    =>  true)
# create Problem

#############################################################################
# What to vary
#############################################################################
argument = "Pexp"
values = ["LambdaLO","NinverseLO", "NinverseNLO"]

for (i,value) in enumerate(values)
    parameters[argument]=value

    global theproblem = QFTdynamicsProblem(modelfile, parameters)
    global thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)

    initialize!(thesolution, tmpdata)

    for t in (thesolution.problem.simsetup.lastEvolStep+1):(thesolution.problem.simsetup.Nsteps+2)
        #println("Working on:[" * string(t) * "," * string(t) * "] "); flush(stdout)
        print("*"); flush(stdout)
        evolve!(thesolution, tmpdata, t)
        if (t-2)%thesolution.problem.simsetup.NstepsbetweenNmeas == 0
            measure!(thesolution,t)
        end
    end
    println("Storing the solution in: " * thesolution.problem.simsetup.datapath * "/SimSolution.jld"); flush(stdout)
    serialize( thesolution.problem.simsetup.datapath * "/SimSolution.jld", thesolution)

    ## plotting
    if i==1
        global plots = Dict()
        plotdata(plots, thesolution, "c")  #create canvas
    end

    plotdata(plots, thesolution, "h", i, string(value)) # fill canvas

end

if isinteractive() == true
    displayplots(plots)
else
    saveplots(plots, thesolution.problem.simsetup.plotpath)
end
