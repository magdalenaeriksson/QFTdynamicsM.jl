###################################################################################################################################
using QFTdynamics
using Serialization
using Plots

###################################################################################################################################
# Set up
###################################################################################################################################
modelfile = TwoPIScalarfile()

parameters = Dict(
    # Discretisation
    "Mass"  =>  0.5, "Nx"    => 16, "sdim"  =>  2, "dt"    =>  0.1, "Nsteps"=>  200, "NstepsinMemory"  =>  0, "Nmeas" =>  0, 
    #model
    "Mod"  =>  "phi4", "Lambda"=>  5, "ONgroup"=> 8,
    #initialisation
    "init"  =>  "Quench", "xi"     =>  0.4, "sig"     =>  1, "eta"     =>  5, "T"     =>  0, "n"     =>  2, 
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
argument = "Lambda"
values = [0.5,5,20]
#argument = "Pexp"
#values = ["LambdaLO","NinverseNLO"]

for (i,value) in enumerate(values)
    parameters[argument]=value

    global theproblem = QFTdynamicsProblem(modelfile, parameters)
    global thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)

    initialize!(thesolution)

    for t in (thesolution.problem.simsetup.lastEvolStep+1):(thesolution.problem.simsetup.Nsteps+2)
        #println("Working on:[" * string(t) * "," * string(t) * "] "); flush(stdout)
        if (t-2)%10 == 0 print("*") end
        evolve!(thesolution, tmpdata, t)
        if (t-2)%thesolution.problem.simsetup.NstepsbetweenNmeas == 0
            measure!(thesolution,t)
        end
    end
    #println("Storing the solution in: " * thesolution.problem.simsetup.datapath * "/SimSolution.jld"); flush(stdout)
    #serialize( thesolution.problem.simsetup.datapath * "/SimSolution.jld", thesolution)

    local avgM2 = sum([thesolution.measurearray[i].M2 for i in 1:length(thesolution.measurearray)])/length(thesolution.measurearray)
    println("average Hartree mass:", avgM2)

    ## plotting
    if i==1
        global plots = Dict()
        plotdata(plots, thesolution, "c")  #create canvas
    end

    plotdata(plots, thesolution, "h", i, string(value)) # fill canvas
    println(" *********************************************************** ")

end

if isinteractive() == true
    display(plots["Hartreemass.png"])
    display(plots["stat0prop.png"])
    display(plots["Sigmar.png"])
    display(plots["SigmaF.png"])
    display(plots["cnumber.png"])
else
    saveplots(plots, thesolution.problem.simsetup.plotpath)
end
