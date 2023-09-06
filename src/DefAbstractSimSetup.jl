#"""
#this file defines SimSetupuration used to control the Simulation steps and Measurement intervalls.
#"""
export SimSetup

export splitter
function splitter(idxmin::Int,n::Int, k::Int)
    xz = Int.(ceil.(range(0, n, length = k+1)))
    return [(xz[i]+idxmin):(xz[i+1]+(idxmin-1)) for i in 1:k]
end

function testsplitter(idxmin::Int, n::Int, nchunks::Int)
    a = splitter(idxmin, n, nchunks)
    res = 0
    for i in 1:length(a)
        res += sum(a[i])
    end
    if res != Int((n/2)*(2*idxmin+(n-1))) print("***idx splitter is wrong") end
end

function check_parameters!(parameters)
  # basic checks
  if parameters["NstepsinMemory"] == 0
    println("* Changing NstepsinMemory to (Nsteps+2)"); flush(stdout)
    parameters["NstepsinMemory"] = parameters["Nsteps"] + 2
  end
  if parameters["Nmeas"] == 0
    println("* Changing Nmeas to Nsteps"); flush(stdout)
    parameters["Nmeas"] = parameters["Nsteps"] 
  end
end

Base.@kwdef mutable struct SimSetup
  # Parameterstring
  parameterstring::String=""
  # Paths
  datapath::String=""
  postpath::String=""
  plotpath::String=""
  # Steps between Measurements
  Nsteps::Int64=0
  Nmeas::Int64=0
  NstepsbetweenNmeas::Int64=0
  #newsimulation::Bool=true # maybe merge with lastEvolStep
  lastEvolStep::Int64=0
  # Memory
  NstepsinMemory::Int64=0
end

#function Base.copy(simsetup::SimSetup)
#  return SimSetup(simsetup.parameterstring,
#                  simsetup.datapath,
#                  simsetup.postpath,
#                  simsetup.plotpath,
#                  simsetup.Nsteps,
#                  simsetup.Nmeas,
#                  simsetup.NstepsbetweenNmeas,
#                  simsetup.lastEvolStep,
#                  simsetup.NstepsinMemory)
#end

function isnewsimulation(simsetup::SimSetup) return simsetup.lastEvolStep == 0 end

function createfilestructandSimSetup(model, parameters)
  # create parameterstring for runpath
	parameterstring = getParameterstring(model, parameters)
  
  # get datapath
	lines = readlines("etc/paths.txt")
	datapath = chop(lines[2],head=5,tail=0)
	postpath = chop(lines[3],head=5,tail=0)
	plotpath = chop(lines[4],head=5,tail=0)
   
  # define dirs
	runpath = datapath * "/" * parameterstring
	runpostpath = postpath * "/" * parameterstring
	runplotpath = plotpath * "/" * parameterstring

  # create SimSetup
  simsetup = SimSetup()
  simsetup.parameterstring = parameterstring
  simsetup.datapath = runpath
  simsetup.postpath = runpostpath
  simsetup.plotpath = runplotpath
  simsetup.NstepsinMemory = parameters["NstepsinMemory"]

  if parameters["ov"] == true
    # setting everything up - Note this lines are also used for a similar use case, see below
	  println("Purging and creating files in... ", runpath); flush(stdout)
    rm(runpath, force=true, recursive=true)
    rm(runplotpath, force=true, recursive=true)
    rm(runplotpath * "/pdf", force=true, recursive=true)
    rm(runpostpath, force=true, recursive=true)
		mkpath(runpath)
		mkpath(runpostpath)
		mkpath(runplotpath)
		mkpath(runplotpath* "/pdf")
    
    cp("etc/index.php", joinpath( join(split(runplotpath,"/")[1:end-1],"/"),"index.php"), force=true )
	  cp("etc/index.php", joinpath(runplotpath,"index.php"), force=true )
	  cp("etc/index.php", joinpath(runplotpath * "/pdf" ,"index.php"), force=true )
    chmod( join(split(runplotpath,"/")[1:end-1],"/") , 0o777, recursive=true)

    simsetup.lastEvolStep=0
    simsetup.NstepsbetweenNmeas = floor( parameters["Nsteps"] / parameters["Nmeas"])
    simsetup.Nsteps = simsetup.NstepsbetweenNmeas * parameters["Nmeas"]
    simsetup.Nmeas = parameters["Nmeas"]
    println("Nsteps: ", simsetup.Nsteps, ", Nmeas: ", simsetup.Nmeas, ", NstepsbetweenNmeas: ", simsetup.NstepsbetweenNmeas )
  #else	    
    #if isdir(SimSetup.datapath)
    #  # datapath already exists -
    #  println("Datapath already exists...Loading SimData from: ", SimSetup.datapath)
    #  SimResult = deserialize( SimSetup.datapath * "/SimResult.jld")
    #  println("Checking if SimPar is compatible with current arguments...")
    #  if !(SimResult.SimPar == SimPar)
    #    # this should actually never be the case, since if the SimPar are not compatible then we have a different runpath!?
    #    println("Loaded SimSetup is not compatible with current arguments...aborting: ")
    #    return 0
    #  else
    #    println("Loaded SimSetup is compatible with current arguments...continuing the simulation..."); flush(stdout)
    #    # use saved SimSetup and SimPar from arguments
    #    SimSetup = SimResult.SimSetup
    #    # update stuff caused by a larger SimPar.Nsteps
    #    SimSetup.contofsim=true 
    #    SimSetup.Nmeas = floor(SimPar.Nsteps / SimSetup.NstepsbetwNmeas)
    #    SimSetup.Nsteps = SimSetup.Nmeas * SimSetup.NstepsbetwNmeas

    #    if !isfile(SimSetup.datapath * "/fields.jld")
    #      println("fields.jld does not exist -> CANNOT CONTINUE SIMULATION")
    #    end

    #  end
    #else
	   # println("Creating files in... ", runpath); flush(stdout)
    #  rm(runpath, force=true, recursive=true)
    #  rm(runplotpath, force=true, recursive=true)
    #  rm(runpostpath, force=true, recursive=true)
		  #mkpath(runpath)
		  #mkpath(runpostpath)
		  #mkpath(runplotpath)

    #  SimSetup.contofsim=false 
    #  SimSetup.lastEvolStep=0
    #  SimSetup.NstepsbetwNmeas = floor(SimPar.Nsteps / SimPar.Nmeas)
    #  SimSetup.Nsteps = SimSetup.NstepsbetwNmeas * SimPar.Nmeas
    #  SimSetup.Nmeas = SimPar.Nmeas
    #end
  end
   
  return simsetup
end