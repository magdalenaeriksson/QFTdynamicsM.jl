################################################################################################################################
# defines test function for TwoPIScalar models
################################################################################################################################
using QFTdynamics
using BenchmarkTools

function TwoPIScalar_exampletest()
  evolveto = 10
  # testparameterfile
  parameters = Dict(
      # Discretisation
      "Mass"  =>  1.0, "Nx"    =>  6, "sdim"  =>  2, "dt"    =>  0.05, "Nsteps"=>  evolveto, "NstepsinMemory"  =>  0, "Nmeas" =>  0, 
      #model
      "Mod"  =>  "phi4", "Lambda"=>  1, "ONgroup"=> 8,
      #initialisation
      "init"  =>  "Pnr", "xi"     =>  0.4, "sig"     =>  1, "eta"     =>  5, "T"     =>  5, "n"     =>  2, 
      #Reno #PertExp
      "Reno"  =>  "RNone", "Pexp"  =>  "LambdaNLO",
      # Numerics
      "Num"   =>  "CPUfull", "Nchunks"   =>  0, # if 0 -> set to @Threads.nthreads
      #extra
      "tag"   =>  "AndersNNLOpaper0dim", "ov"    =>  true)

  return TwoPIScalar_testevolution(TwoPIScalarfile(), parameters, evolveto)
end

function TwoPIScalar_testevolution(modelfile::TwoPIScalarfile, parameters::Dict, evolveto::Int64)
    parameters["tag"] = "testevolution" 

    #
    # create Problem
    #
    theproblem = QFTdynamicsProblem(modelfile, parameters)
    thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)
    initialize!(thesolution, tmpdata)

    #
    # evolve
    #
    for t in (thesolution.problem.simsetup.lastEvolStep+1):evolveto
      print("*"); flush(stdout)
      evolve!(thesolution, tmpdata, t)
      if (t-2)%thesolution.problem.simsetup.NstepsbetweenNmeas == 0
          measure!(thesolution,t)
      end
    end
    return (thesolution.simdata, thesolution.measurearray)
end

function TwoPIScalar_testevolution_returnall(modelfile::TwoPIScalarfile, parameters::Dict, evolveto::Int64)
    parameters["tag"] = "testevolution" 

    #
    # create Problem
    #
    theproblem = QFTdynamicsProblem(modelfile, parameters)
    thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)
    initialize!(thesolution, tmpdata)

    #
    # evolve
    #
    for t in (thesolution.problem.simsetup.lastEvolStep+1):evolveto
      print("*"); flush(stdout)
      evolve!(thesolution, tmpdata, t)
      if (t-2)%thesolution.problem.simsetup.NstepsbetweenNmeas == 0
          measure!(thesolution,t)
      end
    end
    return (thesolution, tmpdata)
end
