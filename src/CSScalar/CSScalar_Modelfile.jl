using FFTWhelper

export CSScalarfile
struct CSScalarfile <: AbstractModelfile end # defnine modelfile as type

include("CSScalar_Problemdefinition.jl") #this is a "header", defining all the abstract types for Model, PertExpansion, Renormalization ...

include("CSScalar_SimDataStruct.jl")
include("CSScalar_Discretization.jl")
include("CSScalar_ModelExpansionRenormalization.jl")

# actual implementation of types (Model, PertExpansion, Renormalization ...) defined in "header"


#
# Create Problem
#
function getParameterstring(modelfile::CSScalarfile, parameters)
  parameterstring =  parameters["Mod"]
  parameterstring *= "CS"
  parameterstring *= "_L"   * string( Int64(floor(parameters["Lambda"]*100))	) 
  parameterstring *= "_ON"  * string( parameters["ONgroup"] 	)
  parameterstring *= "_"   * parameters["Reno"]
    if parameters["init"]=="Thermal" parameterstring *= "_T" * string( Int64(parameters["T"]*10) ) end
    if parameters["init"]=="Pnr" parameterstring *= "_n" * string( Int64(parameters["n"]*10)) end
  # Discretization
  parameterstring *= "_"
  parameterstring *=  "Nx" 		  * string( parameters["Nx"]      	) 
  parameterstring *= "_sdim" 	  * string( parameters["sdim"]     ) 
  parameterstring *= "_Nsteps"  * string( parameters["Nsteps"]		)  
  parameterstring *= "_dt"    	* string( Int64(parameters["dt"]*1000)	)  
  parameterstring *= "_M"     	* string( Int64(parameters["Mass"]*100)	)
  parameterstring *= "_Runs"    * string( Int64(parameters["Runs"])	)
  if parameters["tag"]!=""
    parameterstring *= "_"* parameters["tag"]
  end

end

function QFTdynamicsProblem(modelfile::CSScalarfile, parameters)
  # add NstepsinMemory otherwise problem with 2PI code
  parameters["NstepsinMemory"] = 0
  if parameters["Nmeas"] == 0
    parameters["Nmeas"] = parameters["Nsteps"]
  end
  @show parameters
  ## simsetup
  simsetup = createfilestructandSimSetup(modelfile, parameters)
  @show simsetup.Nsteps

  # Discretization
  helper = getfftwhelper(parameters["Nx"], parameters["sdim"])
  disc = CSScalarDiscretizationLattice(  parameters["Mass"],
                          parameters["Nx"],
                          parameters["sdim"],
                          parameters["Nsteps"],
                          parameters["dt"],
                          parameters["Nx"]^parameters["sdim"],
                          1/(parameters["Nx"]^parameters["sdim"]),
                          parameters["Nmeas"],
                          helper,
                          [lev.deg for lev in helper])

  # Model
  if parameters["Mod"] == "CSphi4" model = CSPhi4( parameters["Lambda"], parameters["ONgroup"], parameters["Mass"]^2) end
  if parameters["Mod"] == "CSphi4tachyonicUNS" model = CSPhi4TachyonicUNS( parameters["Lambda"], parameters["ONgroup"], -parameters["Mass"]^2) end
  if parameters["Mod"] == "CSphi4tachyonicALL" model = CSPhi4TachyonicALL( parameters["Lambda"], parameters["ONgroup"], -parameters["Mass"]^2) end

  # Initialization
  if parameters["init"] == "Thermal"  init = CSScalarThermal(parameters["T"]) end
  if parameters["init"] == "Pnr"      init = CSScalarParticle(parameters["n"]) end
  if parameters["init"] == "TopHatT1" init = CSScalarTopHatT1() end
  if parameters["init"] == "TopHatT2" init = CSScalarTopHatT2() end
  if parameters["init"] == "TopHatT3" init = CSScalarTopHatT3() end

  # Perturbative espansion
  pexp = CSScalarPertExpansion()

  # Renormalization
  if parameters["Reno"] == "RNone" reno = CSScalarRNone() end
  #if parameters["Reno"] == "RMass" reno = CSScalarRMass() end

  # Numerics
  if parameters["Num"] == "CPU"    num = CSScalarCPU(parameters["Runs"], parameters["seed"], Threads.nthreads(),  splitter( 1, parameters["Runs"], Threads.nthreads()) ) end

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

export QFTdynamicsSolutionCSScalar
mutable struct QFTdynamicsSolutionCSScalar <: QFTdynamicsSolution
  problem::QFTdynamicsProblem
  simdata::Vector{CSScalarSimDataCPU}
  measurearray::Vector{Measurement} # elements are timeslices
  measurearrayofruns::Vector{Vector{Measurement}} # elements are timeslices
end


function QFTdynamicsSolution(modelfile::CSScalarfile, problem::QFTdynamicsProblem)
  # Get Memory
  simdata = Vector{CSScalarSimDataCPU}(undef,problem.num.Runs)
  for i in 1:problem.num.Runs
    simdata[i] = CSScalarSimDataCPU(problem.disc.Nx, problem.disc.sdim, problem.model.ONgroup)
  end
  
  # get tmpdata  
  tmpdata = Vector{CSScalarTmpDataCPU}(undef, problem.num.threads)
  for i in 1:problem.num.threads
    tmpdata[i] = CSScalarTmpDataCPU(simdata[1])
  end

  # get measurearray
  measurearray = Vector{Measurement}(undef, problem.simsetup.Nmeas+1) # includes 0 measurement

  # get measurearrayofruns
  measurearrayofruns = Vector{Vector{Measurement}}(undef, problem.simsetup.Nmeas+1) # includes 0 measurement
  for i in 1:problem.simsetup.Nmeas+1
    measurearrayofruns[i] = Vector{Measurement}(undef, problem.num.Runs)
  end

  return QFTdynamicsSolutionCSScalar(problem, simdata, measurearray, measurearrayofruns), tmpdata
end


include("CSScalar_Initialization.jl")
#include("CSScalar_Evolution.jl") # not part of the package
include("CSScalar_Measurement.jl")
include("CSScalar_Plotting.jl")