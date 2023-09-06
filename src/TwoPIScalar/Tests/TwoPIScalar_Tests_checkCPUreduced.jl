using QFTdynamics
using Serialization
using Parameters
using FFTWhelper

include("../TwoPIScalar_Testhelpers.jl")

function TestTwoPIScalar_CPUred(Nx, sdim, Nchunks, Pexpstring, evolveto)
    #
    # set up parameters
    #
    parameters = Dict(
        # Discretisation
        "Mass"  =>  1,
        "Nx"    =>  Nx, 
        "sdim"  =>  sdim, 
        "dt"    =>  0.01, #0.05
        "Nsteps"=>  evolveto,
        "NstepsinMemory"  =>  0, # 3 - min for LO
        "Nmeas" =>  0, # either 0 (all timesteps are measured) or >=3 (because the plot script wants to draw 3 graphs)
        #model
        "Mod"  =>  "phi4",
        "Lambda"=>  10,
        "ONgroup"=> 8,
        #initialisation
        "init"  =>  "Thermal",
        "xi"    => 0.4,
        "eta"   => 5.0,
        "sig"   => 1.0,
        "T"     =>  5,
        "n"     =>  2,
        #Reno
        "Reno"  =>  "RNone",
        #"Reno"  =>  "RMass",
        #PertExp
        "Pexp"  =>  Pexpstring,
        # Numerics
        "Num"   =>  "CPUfull",
        "Nchunks"   =>  Nchunks, # if 0 -> set to @Threads.nthreads
        #extra
        "tag"   =>  "testing", #dont use "_" in the tag
        "ov"    =>  true)

    parameters["Num"] = "CPUfull"
    parameters["NstepsinMemory"] = 0
    parameters["Nmeas"] = 0
    (simdata, measurearray) = TwoPIScalar_testevolution(TwoPIScalarfile(), parameters, evolveto)

    parameters["Num"] = "CPUred"
    parameters["NstepsinMemory"] = 0
    parameters["Nmeas"] = 0
    (simdata2, measurearray2) = TwoPIScalar_testevolution(TwoPIScalarfile(), parameters, evolveto)

    #
    # comparison on level of F 
    #
    # copy full to reduced opject
    F = zero(simdata2.F[evolveto,evolveto]) # create emty reduced lattice
    copytoreducedlattice!( F, simdata.F[evolveto,evolveto], getfftwhelper(parameters["Nx"],parameters["sdim"] ) )

    # compare F (result of CPUfull calc, but copied to reduced), and simdata.F2[evolveto,evolveto] (result of CPUred calc)
    boolarray = isapprox.(F,simdata2.F[evolveto,evolveto])
    if sum( boolarray ) != length(F)
        println("Error in comparison of F from CPUfull vs CPUred")
        println("correct vs total number of elements: ", sum(boolarray), "/", length(F))
        println("Note: For evolveto > 40 i get different results - they are minor but still puzzling?!")
        @show F[1]
        @show simdata2.F[evolveto,evolveto][1]
        return false
    end

    #
    # comparison on level of measurearray F
    #
    # -> not necessary
    return true
end

#TestTwoPIScalar_CPUred(4, 2, 2, "LambdaNLO", 2)
TestTwoPIScalar_CPUred(64,  2, 1, "NinverseNLO", 20)
#TestTwoPIScalar_CPUred(12,  2, 1, "NinverseNLO", 2)

#@test TestTwoPIScalar_CPUred(512, 1, 4, "NinverseNLO", 50)
#@test TestTwoPIScalar_CPUred(64,  2, 1, "NinverseNLO", 50)
#@test TestTwoPIScalar_CPUred(16,  3, 5, "NinverseNLO", 50)
#@test TestTwoPIScalar_CPUred(512, 1, 4, "LambdaNLO", 50)
#@test TestTwoPIScalar_CPUred(64,  2, 1, "LambdaNLO", 50)
#@test TestTwoPIScalar_CPUred(16,  3, 5, "LambdaNLO", 50)