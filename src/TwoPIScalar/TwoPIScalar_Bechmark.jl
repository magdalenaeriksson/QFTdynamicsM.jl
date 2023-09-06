################################################################################################################################
# defines benchmark function for TwoPIScalar models
################################################################################################################################
using QFTdynamics
using BenchmarkTools

function approximateMemoryusage( parameters::Dict, nthreads)
  if parameters["NstepsinMemory"] == 0
    NstepsinMemory = parameters["Nsteps"] + 2
  else
    NstepsinMemory = parameters["NstepsinMemory"] 
  end
  GB = 0
  if parameters["Num"] == "CPUfull"
    GB += MemoryUseofTwoPIScalarSimDataCPUfull( NstepsinMemory, parameters["Nx"], parameters["sdim"])
    GB += MemoryUseofTwoPIScalarTmpDataCPUfull( NstepsinMemory, parameters["Nx"], parameters["sdim"], nthreads)
  end
  if parameters["Num"] == "CPUred"
    GB += MemoryUseofTwoPIScalarSimDataCPUred( NstepsinMemory, parameters["Nx"], parameters["sdim"])
    GB += MemoryUseofTwoPIScalarTmpDataCPUred( NstepsinMemory, parameters["Nx"], parameters["sdim"], nthreads)
  end
  return GB
end

function TwoPIScalar_examplebenchmark()
  # testparameterfile
  parameters = Dict(
      # Discretisation
      "Mass"  =>  1.0, "Nx"    =>  32, "sdim"  =>  3, "dt"    =>  0.05, "Nsteps"=>  200, "NstepsinMemory"  =>  0, "Nmeas" =>  0, 
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

  approxmemgb = approximateMemoryusage( parameters, Threads.nthreads())
  println("approximate Memory/GB:", approxmemgb); flush(stdout)
  (time, memgb) = TwoPIScalar_bechmarkevolution(TwoPIScalarfile(), parameters, 0.)
  println("SimpleBenchmark: Time/s per evolution step: ", string(time), ", Memory/GB of simulation: ", string(memgb))
  (time, memgb) = TwoPIScalar_bechmarkevolution(TwoPIScalarfile(), parameters, 10.)
  println("ToolsBenchmark: Time/s per evolution step: ", string(time), ", Memory/GB of simulation: ", string(memgb))
end

function TwoPIScalar_bechmarkevolution(modelfile::TwoPIScalarfile, parameters::Dict, timeforbenchmark::Float64)
    # if timeforbenchmark ==0, then it simply takes 10 evolution steps and measures the time it takes
    # The issue is that for Memory intensive configurations Benchmark tools does not work since it think it passes it a copy!
    parameters["tag"] = "benchmark" 

    #
    # create Problem
    #
    theproblem = QFTdynamicsProblem(modelfile, parameters)
    thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)
    initialize!(thesolution,tmpdata)

    #
    # fill up memory - actually we only need to expandSimData - benchmarking does not depend on which values we ve stored in the objects
    #
    # Lambda expansion: expandSimData is enough
    if typeof(thesolution.problem.pexp) == QFTdynamics.TwoPIScalarLambdaNLO
      for t in (thesolution.problem.simsetup.lastEvolStep+1):(thesolution.simdata.NstepsinMemory)
        expandSimData!(thesolution.simdata) # now indices range from 1:NstepsinMemory
      end
    end

    # 1/N expansion: expandSimData is not enough - dunno why
    if typeof(thesolution.problem.pexp) == QFTdynamics.TwoPIScalarNinverseNLO
      for t in (thesolution.problem.simsetup.lastEvolStep+1):(thesolution.simdata.NstepsinMemory)
        expandSimData!(thesolution.simdata) # now indices range from 1:NstepsinMemory
        nrofindices = (thesolution.simdata.indices[2]-1) - thesolution.simdata.indices[1] + 1
        tmpdata.threadranges .= splitter(thesolution.simdata.indices[1], nrofindices, thesolution.problem.num.nchunks)
        fakecalcSelfEnergies!(thesolution.problem.model, thesolution.problem.pexp, thesolution.problem.disc, thesolution.simdata, tmpdata, 2)
      end
    end

    # 
    # measure memory usage
    # 
    membytes_simdata = Base.summarysize(thesolution.simdata)
    membytes_tmpdata = Base.summarysize(tmpdata)
    membytes_problem = Base.summarysize(theproblem)

    membytes = membytes_simdata + membytes_tmpdata + membytes_problem 
    memgb = membytes*10^(-9)

    #
    # benchmark one step
    #
    if timeforbenchmark == 0
      print("Running simple benchmark...");flush(stdout)
      # call evole once to get rid off compilation time
      evolve!(thesolution,tmpdata,thesolution.simdata.indices[2])
      nsteps = 5
      time = @elapsed for i in (thesolution.simdata.NstepsinMemory+1):(thesolution.simdata.NstepsinMemory+1+nsteps-1)
        print("*");flush(stdout)
        evolve!(thesolution,tmpdata,thesolution.simdata.indices[2])
      end
      time = time / nsteps
      println("time [s] per step: ", time);flush(stdout)
    else
      print("Running BenchmarkTools benchmark...");flush(stdout)
      time = @belapsed evolve!($thesolution,$tmpdata,$thesolution.simdata.indices[2]) seconds=timeforbenchmark
      # rule of thumb: benchmark time should be 7 times a single execution
      # run again if this is the case
      if time*7 > timeforbenchmark 
        print("Rerunning benchmark for higher accuracy...");flush(stdout)
        timeforbenchmark = 10*time
        time = @belapsed evolve!($thesolution,$tmpdata,$thesolution.simdata.indices[2]) seconds=timeforbenchmark
      end
      println("");flush(stdout)
    end

    # free memory?
    theproblem = 0
    thesolution = 0
    tmpdata = 0 
    return (time, memgb)
end

