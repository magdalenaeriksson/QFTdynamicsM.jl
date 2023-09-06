###################################################################################################################################
using QFTdynamics
using Serialization
###################################################################################################################################
# Set up
###################################################################################################################################
# get parameters as dict from parameter file and specify modelfile (a type)
include("src/CSGaugeScalar/CSGaugeScalar_Parfile.jl") # Define modelfile and parameters
#thesamples = deserialize("/Users/magdalenaeriksson/code/2PIcode/data/MCSampledIC_Nx32_sdim2_Mass100_n0_Samples796_B1_ith5_test/Samples.jld")
#@show thesamples[1].phi1x


## create problem
theproblem = QFTdynamicsProblem(modelfile, parameters)
## create solution
thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)

initialize!(thesolution, tmpdata)
#initializeGaugeProps!(thesolution, tmpdata)
#evolve!(thesolution, tmpdata, 3)
#measure!(thesolution, tmpdata, 3)
#evolve!(thesolution, tmpdata, 4)
#measure!(thesolution, tmpdata, 4)

elapsedtime = @elapsed for t in (thesolution.problem.simsetup.lastEvolStep+1):(thesolution.problem.simsetup.Nsteps+2)
    println("Working on timestep: " * string(t)); flush(stdout)

    evolve!(thesolution, tmpdata, t)
    #evolveSerial!(thesolution, tmpdata, t)

    if (t-2)%thesolution.problem.simsetup.NstepsbetweenNmeas == 0
        measure!(thesolution, tmpdata, t)
    end
end
@show elapsedtime

## plotting
plots = Dict()
plotdata(plots, thesolution, "c")  #create canvas
plotdata(plots, thesolution, "h", 1, "myplot") # fill canvas
if isinteractive() == true
    displayplots(plots)
    saveplots(plots, thesolution.problem.simsetup.plotpath)
else
    saveplots(plots, thesolution.problem.simsetup.plotpath)
end
