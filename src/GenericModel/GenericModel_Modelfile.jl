export GenericModelFile
struct GenericModelFile <: AbstractModelfile end # defnine modelfile as type
include("GenericModel_Problemdefinition.jl") #this is a "header", defining all the abstract types for Model, PertExpansion, Renormalization ...

include("GenericModel_SimDataStruct.jl")
# actual GenericModelion of types (Model, PertExpansion, Renormalization ...) defined in "header"
include("GenericModel_ModelExpansionRenormalization.jl")
include("GenericModel_Discretization.jl")

#
# Create Problem
#
function getParameterstring(modelfile::GenericModelFile, parameters)
  # maybe add if statements depening on parameters, eg parameters["Model"]
  parameterstring =  "GenericModel"
  parameterstring *= "_Omega"   * string( Int64(parameters["Omega"]*100)	) 
  # Discretization
  parameterstring *= "/"
  parameterstring *=  "Nx" 		  * string( parameters["Nx"]      	) 
  parameterstring *= "_sdim" 	  * string( parameters["sdim"]     ) 
  parameterstring *= "_Nsteps"  * string( parameters["Nsteps"]		)  
  parameterstring *= "_dt"    	* string( Int64(parameters["dt"]*1000)	)  
  if parameters["tag"]!=""
    parameterstring *= "_"* parameters["tag"]
  end
end

function QFTdynamicsProblem(modelfile::GenericModelFile, parameters)
  # check parameters
  check_parameters!(parameters)

  ## SimConfig
  simsetup = createfilestructandSimSetup(modelfile, parameters)

  # Discretization
  disc = DiscretizationGenericModel(  parameters["Mass"],
                          parameters["Nx"],
                          parameters["sdim"],
                          parameters["Nsteps"],
                          parameters["dt"],
                          parameters["NstepsinMemory"],
                          parameters["Nmeas"],)
  # Model
  if parameters["Mod"] == "GenericModelA" 
      model = GenericModelA( parameters["Omega"])
  end

  # Initialization
  #if parameters["init"] == "Thermal"  init = aninitialisation() end
  init = aninitialization()

  # Perturbative Expansion
  #if parameters["Pexp"] == "LambdaLO" pexp =  anexpansion() end
  pexp =  anexpansion() 

  # Renormalization
  #if parameters["Reno"] == "RNone" reno = anrenormalization() end
  reno = arenormalization()

  # initialize the Problem
  theproblem = QFTdynamicsProblem(
    model,
    pexp,
    disc,
    init,
    reno,
    simsetup)

  return theproblem
end

export QFTdynamicsSolutionGenericModel
mutable struct QFTdynamicsSolutionGenericModel <: QFTdynamicsSolution
  problem::QFTdynamicsProblem
  simdata::AbstractSimData
  measurearray::Vector{Measurement}
end

function QFTdynamicsSolution(modelfile::GenericModelFile, problem::QFTdynamicsProblem)
  # Get Memory
  simdata = GenericModelSimData(problem.simsetup.NstepsinMemory, problem.disc.Nx, problem.disc.sdim)
  PrintMemoryofSimData!(simdata)
  #if SimPar.binning != 0 rfftwhelper = rebin(fftwhelper, SimPar.binning) else rfftwhelper=fftwhelper end # bins includes 0 bin!
  measurearray = Vector{Measurement}(undef, problem.simsetup.Nmeas+1) # includes 0 measurement
    
  return QFTdynamicsSolutionGenericModel(problem, simdata, measurearray)
end

include("GenericModel_Initialization.jl")
include("GenericModel_Evolution.jl")
include("GenericModel_Measurement.jl")
include("GenericModel_Plotting.jl")