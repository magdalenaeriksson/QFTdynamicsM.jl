using FFTWhelper

export TwoPIScalarfile
struct TwoPIScalarfile <: AbstractModelfile end # defnine modelfile as type
include("TwoPIScalar_Problemdefinition.jl") #this is a "header", defining all the abstract types for Model, PertExpansion, Renormalization ...

include("TwoPIScalar_SimDataStruct.jl")
# actual implementation of types (Model, PertExpansion, Renormalization ...) defined in "header"
include("TwoPIScalar_Discretization.jl")
include("TwoPIScalar_ModelExpansionRenormalization.jl")

#
# Create Problem
#
function getParameterstring(modelfile::TwoPIScalarfile, parameters)
  parameterstring =  parameters["Mod"]
  parameterstring *= parameters["Pexp"]
  parameterstring *= "_L"   * string( Int64(floor(parameters["Lambda"]*100))	) 
  parameterstring *= "_ON"  * string( parameters["ONgroup"] 	)
  parameterstring *= "_"   * parameters["Reno"]
  if parameters["init"]=="Thermal"
    parameterstring *= "_T" * string( Int64(parameters["T"]*10) )
  elseif parameters["init"]=="Pnr" 
    parameterstring *= "_n" * string( Int64(parameters["n"]*10))
  else
    parameterstring *= "_" * parameters["init"]
  end
  # Discretization
  parameterstring *= "_"
  if parameters["Num"] != "CPUcont"
    parameterstring *=  "Nx" 		  * string( parameters["Nx"]      	) 
    parameterstring *= "_sdim" 	  * string( parameters["sdim"]     ) 
  else
    parameterstring *= "_Nmom" 	  * string( parameters["Nmom"]     ) 
  end
  parameterstring *= "_Nsteps"  * string( parameters["Nsteps"]		)  
  parameterstring *= "_dt"    	* string( Int64(parameters["dt"]*1000)	)  
  parameterstring *= "_M"     	* string( Int64(parameters["Mass"]*100)	)
  parameterstring *= "_NstepsinMemory"    * string( parameters["NstepsinMemory"]		  )
  parameterstring *= "_Num"     * parameters["Num"]		 
  if parameters["tag"]!=""
    parameterstring *= "_"* parameters["tag"]
  end
end

function QFTdynamicsProblem(modelfile::TwoPIScalarfile, parameters)
  # check parameters
  check_parameters!(parameters)

  ## simsetup
  simsetup = createfilestructandSimSetup(modelfile, parameters)

  # Discretization
  if parameters["Num"] != "CPUcont"
    disc = TwoPIScalarDiscretizationLattice(  parameters["Mass"],
                          parameters["Nx"],
                          parameters["sdim"],
                          parameters["Nsteps"],
                          parameters["dt"],
                          parameters["Nx"]^parameters["sdim"],
                          1/(parameters["Nx"]^parameters["sdim"]),
                          parameters["NstepsinMemory"],
                          parameters["Nmeas"],
                          getfftwhelper(parameters["Nx"], parameters["sdim"]) )
  else 
    disc = TwoPIScalarDiscretizationCont(  parameters["Mass"],
                          parameters["Nmom"],
                          parameters["Nsteps"],
                          parameters["dt"],
                          parameters["NstepsinMemory"],
                          parameters["Nmeas"],
                          getfftwhelper(parameters["Nmom"]) )
  end

  # Model
  if parameters["Mod"] == "phi4" model = Phi4( parameters["Lambda"], parameters["ONgroup"], parameters["Mass"]^2) end
  if parameters["Mod"] == "phi4tachyonic" model = Phi4Tachyonic( parameters["Lambda"], parameters["ONgroup"], parameters["Mass"]^2) end

  # Initialization
  if parameters["init"] == "Thermal"  init = TwoPIScalarThermal(parameters["T"]) end
  if parameters["init"] == "Pnr"      init = TwoPIScalarParticle(parameters["n"]) end
  if parameters["init"] == "TopHatT1" init = TwoPIScalarTopHatT1() end
  if parameters["init"] == "TopHatT2" init = TwoPIScalarTopHatT2() end
  if parameters["init"] == "TopHatT3" init = TwoPIScalarTopHatT3() end
  if parameters["init"] == "Quench"   init = TwoPIScalarQuench() end
  if parameters["init"] == "Tsunami"  init = TwoPIScalarTsunami() end
  if parameters["init"] == "Gauss"    init = TwoPIScalarGauss(parameters["xi"],parameters["eta"],parameters["sig"]) end
  if parameters["init"] == "Box"      init = TwoPIScalarBox() end

  # Perturbative Expansion
  if parameters["Pexp"] == "LambdaLO" pexp =  TwoPIScalarLambdaLO() end
  if parameters["Pexp"] == "LambdaNLO" pexp = TwoPIScalarLambdaNLOquantum() end
  if parameters["Pexp"] == "LambdaNLOclassical" pexp = TwoPIScalarLambdaNLOclassical() end
  if parameters["Pexp"] == "NinverseLO" pexp =  TwoPIScalarNinverseLO() end
  if parameters["Pexp"] == "NinverseNLO" pexp = TwoPIScalarNinverseNLO() end

  # Renormalization
  if parameters["Reno"] == "RNone" reno = TwoPIScalarRNone() end
  if parameters["Reno"] == "RMass" reno = TwoPIScalarRMass() end

  # Numerics
  if parameters["Nchunks"] == 0 parameters["Nchunks"]=Threads.nthreads() end

  if parameters["Num"] == "CPUfull"     num = TwoPIScalarCPUfull(parameters["Nchunks"]) end
  if parameters["Num"] == "CPUred"      num = TwoPIScalarCPUred(parameters["Nchunks"]) end
  if parameters["Num"] == "CPUred2"     num = TwoPIScalarCPUred2(parameters["Nchunks"]) end
  if parameters["Num"] == "CPUcont"     num = TwoPIScalarCPUcont(parameters["Nchunks"]) end
  #if parameters["Num"] == "GPUreduced"  num = TwoPIScalarGPUreduced(parameters["Nchunks"]) end

  # Renormalize
  renomodel = Renormalize(model, pexp, reno, disc)

  # initialize the Problem
  theproblem = QFTdynamicsProblem(
    renomodel,
    pexp,
    disc,
    init,
    reno,
    num,
    simsetup)

  return theproblem
end

export QFTdynamicsSolutionTwoPIScalar
mutable struct QFTdynamicsSolutionTwoPIScalar <: QFTdynamicsSolution
  problem::QFTdynamicsProblem
  simdata::TwoPIScalarSimData
  measurearray::Vector{Measurement}
end

function QFTdynamicsSolution(modelfile::TwoPIScalarfile, problem::QFTdynamicsProblem)
  # Get Memory
  simdata = getTwoPISimData(problem)
  PrintMemoryofSimData!(simdata)
  #if SimPar.binning != 0 rfftwhelper = rebin(fftwhelper, SimPar.binning) else rfftwhelper=fftwhelper end # bins includes 0 bin!
  measurearray = Vector{Measurement}(undef, problem.simsetup.Nmeas+1) # includes 0 measurement
    
  tmpdata = getTwoPIScalarTmpData(simdata, problem) # second argument solely for nr of chunks

  return QFTdynamicsSolutionTwoPIScalar(problem, simdata, measurearray), tmpdata
end

include("TwoPIScalar_Initialization.jl")
include("TwoPIScalar_Evolution.jl")
include("TwoPIScalar_Measurement.jl")
include("TwoPIScalar_Plotting.jl")
