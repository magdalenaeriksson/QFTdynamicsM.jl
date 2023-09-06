#"""
# this file defines abstract types common to all QFTdynamics Simulations
#"""
export AbstractDiscretization
export AbstractModel
export AbstractPertExpansion
export AbstractRenormalization
export AbstractInitialization
export AbstractNumerics

abstract type AbstractDiscretization end
abstract type AbstractModel end
abstract type AbstractPertExpansion end
abstract type AbstractRenormalization end
abstract type AbstractInitialization end
abstract type AbstractNumerics end

export QFTdynamicsProblem
mutable struct QFTdynamicsProblem
  model::AbstractModel
  pexp::AbstractPertExpansion
  disc::AbstractDiscretization
  init::AbstractInitialization
  reno::AbstractRenormalization
  num::AbstractNumerics
  simsetup::SimSetup
end

export Measurement
export QFTdynamicsSolution
abstract type Measurement end
abstract type QFTdynamicsSolution end

# function to save timing
export saveelapsedtime
function saveelapsedtime(thesolution::QFTdynamicsSolution, elapsedtime)
  @unpack problem, simdata, measurearray = thesolution
  @unpack model, pexp, disc, init, reno, num, simsetup = problem

  if split(simsetup.parameterstring,"_")[end] == "benchmark"
    filename = join(split(simsetup.parameterstring,"_")[1:end-1],"_") * ".txt"
    if Threads.nthreads()==1 # when benchmarking we go from 1 to some nr. Starting from 1 -> init new file
      rm( filename , force=true)
    end
    timingfile = open( filename, "a") 
    println("nrofThreads/Time in seconds: " * string(Threads.nthreads()) * " " * string(elapsedtime) * "\n") 
    write(timingfile, string(Threads.nthreads()) * " " * string(elapsedtime) * "\n") 
    close(timingfile)
  end 

end

using BenchmarkTools
export benchmarkfunction

function benchmarkfunction(thesolution::QFTdynamicsSolution)

  if split(thesolution.problem.simsetup.parameterstring,"_")[end] == "benchmark"

    lastevolvedstep = (thesolution.problem.simsetup.Nsteps+2)

    # eventually let each sample work on the same thesolution thing -> create a copy and use
    # evolve!($thesolution, lastevolvedstep) setup=(x = deepcopy($thesolution))
    evolvebenchmarkable = @benchmarkable evolve!(x, $lastevolvedstep) setup=(x = deepcopy($thesolution))
    #evolvebenchmarkable = @benchmarkable evolve!($thesolution, $lastevolvedstep) samples = thesolution.problem.simsetup.NstepsinMemory
    evolvebenchmark = run(evolvebenchmarkable)
    #measurebenchmarkable = @benchmarkable measure!($thesolution, $lastevolvedstep)  samples = thesolution.problem.simsetup.NstepsinMemory
    measurebenchmark = @benchmark measure!($thesolution, $lastevolvedstep)
    evolvetime = median(evolvebenchmark).time
    measuretime = median(measurebenchmark).time

    @unpack problem, simdata, measurearray = thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem
    filename = join(split(simsetup.parameterstring,"_")[1:end-1],"_") * ".txt"
    if Threads.nthreads()==1 # when benchmarking we go from 1 to some nr. Starting from 1 -> init new file
      rm( filename , force=true)
    end
    timingfile = open( filename, "a") 
    println("nrofThreads/Evolvetime/Measuretime: " * string(Threads.nthreads()) * " " * string(evolvetime) * " " * string(measuretime) * "\n") 
    write(timingfile, string(Threads.nthreads()) * " " * string(evolvetime) * " " * string(measuretime) * "\n") 
    close(timingfile)
  end 

end
