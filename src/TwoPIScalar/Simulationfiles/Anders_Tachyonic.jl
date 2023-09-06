###################################################################################################################################
using QFTdynamics
using Serialization

###################################################################################################################################
# Set up
###################################################################################################################################
modelfile = TwoPIScalarfile()

parameters = Dict(
    # Discretisation
    "Mass"  =>  0.7, "Nx"    =>  32, "sdim"  =>  3, "dt"    =>  0.1, "Nsteps"=>  150, "NstepsinMemory"  =>  0, "Nmeas" =>  0, 
    #model
    "Mod"  =>  "phi4tachyonic", "Lambda"=>  1, "ONgroup"=> 4,
    #initialisation
    "init"  =>  "Pnr", "xi"     =>  0.4, "sig"     =>  1, "eta"     =>  5, "T"     =>  5, "n"     =>  0, 
    #Reno #PertExp
    "Reno"  =>  "RMass", "Pexp"  =>  "NinverseNLO",
    # Numerics
    "Num"   =>  "CPUred2", "Nchunks"   =>  0, # if 0 -> set to @Threads.nthreads
    #extra
    "tag"   =>  "AndersTachyonicpaperLO", "ov"    =>  true)
# create Problem

#############################################################################
# What to vary
#############################################################################
argument = "Num"
#values = ["LambdaLO","NinverseLO"]#, "NinverseNLO"]
#values = ["LambdaLO","NinverseLO"]
#values = ["LambdaNLO"]
values = ["CPUred","CPUred2"]

for (i,value) in enumerate(values)
    println("Starting loop..."); flush(stdout)
    parameters[argument]=value

    global theproblem = QFTdynamicsProblem(modelfile, parameters)
    global thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)

    println("Initialization ..."); flush(stdout)
    initialize!(thesolution, tmpdata)

    elapsedtime = @elapsed for t in (thesolution.problem.simsetup.lastEvolStep+1):(thesolution.problem.simsetup.Nsteps+2)
        #println("Working on:[" * string(t) * "," * string(t) * "] "); flush(stdout)
        #if (t-2)%10 == 0 print("* ",t ,"*"); flush(stdout) end
        evolve!(thesolution, tmpdata, t)
        if (t-2)%thesolution.problem.simsetup.NstepsbetweenNmeas == 0
            measure!(thesolution,t)
        end
        if (t-2)%1 == 0
            val = 0
            for j in 1:length(thesolution.problem.disc.fftwhelper)
                val += thesolution.measurearray[t-1].F[j] * thesolution.problem.disc.fftwhelper[j].deg
            end
            val *= thesolution.problem.disc.ivol/thesolution.problem.disc.Mass^2
            println("step:", t , ", time:",thesolution.measurearray[t-1].time ,", F[1]:",thesolution.measurearray[t-1].F[1], ", F[x=0]:", val); flush(stdout)
        end
    end

    println("Simulation time ", value ," [sec] ", elapsedtime)
    #println("Storing the solution in: " * thesolution.problem.simsetup.datapath * "/SimSolution.jld"); flush(stdout)
    #serialize( thesolution.problem.simsetup.datapath * "/SimSolution.jld", thesolution)

    ## plotting
    if i==1
        global plots = Dict()
        plotdata(plots, thesolution, "c")  #create canvas
    end

    plotdata(plots, thesolution, "h", i, string(value)) # fill canvas

end

display(plots["statpropat0.png"])
#if isinteractive() == true
#    displayplots(plots)
#else
#    saveplots(plots, thesolution.problem.simsetup.plotpath)
#end
