using QFTdynamics
using Serialization
using Parameters
using Plots
using LaTeXStrings

################################################################################################################
# Load stuff
################################################################################################################
#path1 = "bhome/blabla/CSsimulation/SimSolution.jld"
#path2 = "bhome/blabla/2PIsimulation/SimSolution.jld"

# path 1: CS code
# path 2: 2PI code 
path1 = "/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarCS|Lambda=0.125|g=0.667|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=0.05|m=1|seed=1|Runs=4|test/SimSolution.jld"
path2 = "/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarLOloop|Lambda=0.125|g=667|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=50|m=100|NstepsinMemory=3|CSnbrOfRuns=4/SimSolution.jld"

#path1 = "/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarCS|Lambda=0.125|g=0.667|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=0.05|m=1|seed=1|Runs=10|test/SimSolution.jld"
#path2 = "/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarLOloop|Lambda=0.125|g=667|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=50|m=100|NstepsinMemory=3|CSnbrOfRuns=10/SimSolution.jld"

#zero coupling: 3 Runs
#path1 = "/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarCS|Lambda=0|g=0.002|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=0.05|m=1|seed=1|Runs=3|test/SimSolution.jld"
#path2 = "/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarLOloop|Lambda=0|g=2|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=50|m=100|NstepsinMemory=3|CSnbrOfRuns=3/SimSolution.jld"

# zero coupling: 6 runs
#path1 = "/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarCS|Lambda=0|g=0.002|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=0.05|m=1|seed=1|Runs=6|test/SimSolution.jld"
#path2 = "/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarLOloop|Lambda=0|g=2|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=50|m=100|NstepsinMemory=3|CSnbrOfRuns=6/SimSolution.jld"

# zero coupling: 12 runs
#path1 = "/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarCS|Lambda=0|g=0.002|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=0.05|m=1|seed=1|Runs=12|test/SimSolution.jld"
#path2 = "/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarLOloop|Lambda=0|g=2|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=50|m=100|NstepsinMemory=3|CSnbrOfRuns=12/SimSolution.jld"

# zero coupling with new initial conds: 3 runs
#path1 = "/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarCS|Lambda=0|g=0.002|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=0.05|m=1|seed=1|Runs=3|newInitalCondsTest/SimSolution.jld"
#path2 = "/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarLOloop|Lambda=0|g=2|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=50|m=100|NstepsinMemory=3|CSnbrOfRuns=3NewInits/SimSolution.jld"

# with new initial conds: 130 runs, standard coupling
#path1 = "/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarCS|Lambda=0.125|g=0.667|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=0.05|m=1|seed=1|Runs=130|newInitalCondsTest/SimSolution.jld"
#path2 = "/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarLOloop|Lambda=0.125|g=667|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=50|m=100|NstepsinMemory=3|CSnbrOfRuns=130_NewInits/SimSolution.jld"

# zero coupling with new initial conds: 8 runs
#path1 = "/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarCS|Lambda=0|g=0.002|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=0.05|m=1|seed=1|Runs=8|newInitalCondsTest/SimSolution.jld"
#path2 = "/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarLOloop|Lambda=0|g=2|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=50|m=100|NstepsinMemory=3|CSnbrOfRuns=8_NewInits/SimSolution.jld"

# with new initial conds: 8 runs, lambda=0.125, g=0.3
#path1="/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarCS|Lambda=0.125|g=0.3|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=0.05|m=1|seed=1|Runs=8|newInitalCondsTest/SimSolution.jld"
#path2="/Users/magdalenaeriksson/code/2PIcode/data/SU(2)gaugeScalarLOloop|Lambda=0.125|g=300|RNone|n=1|Nx=32|sdim=3|Nsteps=200|dt=50|m=100|NstepsinMemory=3|CSnbrOfRuns=8_NewInits/SimSolution.jld"
thesolution1 = deserialize(  path1  )
#@unpack problem1, simdata1, measurearray1 = thesolution1

thesolution2 = deserialize( path2  )

#thesolution3 = deserialize( path3  )
#@unpack problem2, simdata2, measurearray2 = thesolution2

################################################################################################################
# Extract whatever you need
################################################################################################################
time1 = [el.time for el in thesolution1.measurearray]
Prop1 = [el.Phi2k[1] for el in thesolution1.measurearray]

time2 = [el.time for el in thesolution2.measurearray]
Prop2 = [el.FS[1] for el in thesolution2.measurearray]

#time3 = [el.time for el in thesolution3.measurearray]
#Prop3 = [el.Phi2k[123] for el in thesolution3.measurearray]

myplot = plot(xlabel=L"tm", ylabel=L"G(t,t';k=0)")

plot!(myplot, time1, Prop1, label="CS")
plot!(myplot, time2, Prop2, label="2PI")
#plot!(myplot, time3, Prop3, label="CS 10")

#savefig(myplot, <path/filename.png or .pdg>)