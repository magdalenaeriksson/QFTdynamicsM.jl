using QFTdynamics
using Test

include("../src/TwoPIScalar/Tests/TwoPIScalar_Tests_NinverseNLO_checkfastimplementation.jl")
#include("../src/TwoPIScalar/Tests/TwoPIScalar_Tests_NinverseNLO_checkfastimplementation_old.jl") #TestTwoPIScalar_NinverseNLO_full works in 1,2,3 sdim

#include("../src/TwoPIScalar/Tests/TwoPIScalar_Tests_LambdaNLO_checkfastimplementation.jl")
#include("../src/TwoPIScalar/Tests/TwoPIScalar_Tests_checkCPUreduced.jl")

@testset "QFTdynamics.jl" begin
    # TwoPIScalar
    # fast vs simple implementation, Args: Nx, sdim, nchunks
    #@test TestTwoPIScalar_NinverseNLO_full(64,1,3)
    #@test TestTwoPIScalar_NinverseNLO_full(64,2,3)
    #@test TestTwoPIScalar_NinverseNLO_full(32,3,3)
    #@test   TestTwoPIScalar_LambdaNLO_full(1024,1,1)
    #@test   TestTwoPIScalar_LambdaNLO_full(512,2,5)
    #@test   TestTwoPIScalar_LambdaNLO_full(128,3,6)
    # CPUred vs CPUfull, Arg: Nx, sdim, nchunks, Pexpstring, evolveto
    #@test TestTwoPIScalar_CPUred(512, 1, 4, "NinverseNLO", 50)
    #@test TestTwoPIScalar_CPUred(64,  2, 1, "NinverseNLO", 50)
    #@test TestTwoPIScalar_CPUred(16,  3, 5, "NinverseNLO", 50)
    #@test TestTwoPIScalar_CPUred(512, 1, 4, "LambdaNLO", 50)
    #@test TestTwoPIScalar_CPUred(64,  2, 1, "LambdaNLO", 50)
    #@test TestTwoPIScalar_CPUred(16,  3, 5, "LambdaNLO", 50)
end
