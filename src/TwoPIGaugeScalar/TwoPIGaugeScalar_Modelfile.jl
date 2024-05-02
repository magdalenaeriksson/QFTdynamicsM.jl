using FFTWhelper
using DataFrames
using CSV

export TwoPIGaugeScalarfile
struct TwoPIGaugeScalarfile <: AbstractModelfile end # defnine modelfile as type
include("TwoPIGaugeScalar_Problemdefinition.jl") #this is a "header", defining all the abstract types for Model, PertExpansion, Renormalization ...
include("TwoPIGaugeScalar_SimDataStruct.jl")
# actual implementation of types (Model, PertExpansion, Renormalization ...) defined in "header"
include("TwoPIGaugeScalar_Discretization.jl")
include("TwoPIGaugeScalar_ModelExpansionRenormalization.jl")

#
# Create Problem
#
function getParameterstring(modelfile::TwoPIGaugeScalarfile, parameters)
    # Model
    if parameters["Mod"]=="SUNgaugeScalar" parameterstring = "SU" * string(Int64(parameters["N"])) * "gaugeScalar"  end
    if parameters["Mod"]=="U1gaugeScalar" parameterstring = "U1gaugeScalar" end
    parameterstring *= parameters["Pexp"]
    parameterstring *= "_Lambda"   * string( parameters["Lambda"]*100 ) 
    parameterstring *= "_g"        * string( Int64(parameters["g"]*1000))
    parameterstring *= "_"          * parameters["Reno"]
    if parameters["init"]=="Thermo" 
        parameterstring *= "_T" * string(Int64(parameters["T"])) 
    elseif parameters["init"]=="Pnr" 
        parameterstring *= "_n" * string(Int64(parameters["n"]))
    else
        parameterstring *= "_" * parameters["init"]
    end
    parameterstring *= "_Nx"       * string(parameters["Nx"]) 
    parameterstring *= "_sdim" 	* string(parameters["sdim"]) 
    parameterstring *= "_Nsteps"   * string(parameters["Nsteps"])  
    parameterstring *= "_dt"    	* string(Int64(parameters["dt"]*1000))  
    parameterstring *= "_m"     	* string(Int64(parameters["Mass"]*100))
    parameterstring *= "_NstepsinMemory"    * string(parameters["NstepsinMemory"])
    if parameters["tag"]!= "" parameterstring *= "_" * parameters["tag"] end
end

function QFTdynamicsProblem(modelfile::TwoPIGaugeScalarfile, parameters)
    # Check parameters
    check_parameters!(parameters)

    # Simsetup
    simsetup = createfilestructandSimSetup(modelfile, parameters)

    # Discretisation
    disc = TwoPIGaugeScalarDiscretizationLattice(
                        parameters["Mass"],
                        parameters["Nx"],
                        parameters["sdim"],
                        parameters["sdim"]+1,
                        parameters["Nsteps"],
                        parameters["dt"],
                        parameters["Nx"]^parameters["sdim"],
                        1/(parameters["Nx"]^parameters["sdim"]),
                        parameters["NstepsinMemory"],
                        parameters["Nmeas"],
                        getfftwhelper(parameters["Nx"], parameters["sdim"]))

    # Models
    if parameters["Mod"] == "SUNgaugeScalar" 
        model = SUNgaugeScalar( parameters["Lambda"], 
                                parameters["Mass"], 
                                parameters["g"], 
                                parameters["N"],
                                (parameters["N"]^2-1)/(2*parameters["N"])) #CF = (N^2-1)/(2N)
    end  
    
    if parameters["Mod"] == "U1gaugeScalar" 
        model = U1gaugeScalar( parameters["Lambda"], parameters["Mass"], parameters["g"])
    end

    # Initialization
    if parameters["init"] == "Thermal"  init = TwoPIGaugeScalarThermal(parameters["T"])  end
    if parameters["init"] == "Pnr"      init = TwoPIGaugeScalarParticle(parameters["n"]) end
    if parameters["init"] == "TopHatT1" init = TwoPIGaugeScalarTopHatT1() end
    if parameters["init"] == "TopHatT2" init = TwoPIGaugeScalarTopHatT2() end
    if parameters["init"] == "TopHatT3" init = TwoPIGaugeScalarTopHatT3() end

    # Perturbative expansion
    if parameters["Pexp"] == "LOloop"   pexp = TwoPIGaugeScalarLOloop()  end
    if parameters["Pexp"] == "NLOloop"  pexp = TwoPIGaugeScalarNLOloop() end

    # Renormalization
    if parameters["Reno"] == "RNone" reno = TwoPIGaugeScalarRNone() end
    if parameters["Reno"] == "RMass" reno = TwoPIGaugeScalarRMass() end

      # Numerics
    if parameters["Nchunks"] == 0 parameters["Nchunks"]=Threads.nthreads() end
    if parameters["Num"] == "CPUfull"     num = TwoPIGaugeScalarCPUfull(parameters["Nchunks"]) end

    # initialize the problem
    theproblem = QFTdynamicsProblem(
        model,
        pexp, 
        disc, 
        init, 
        reno, 
        num, 
        simsetup)

    return theproblem
end

export QFTdynamicsSolutionTwoPIGaugeScalar
mutable struct QFTdynamicsSolutionTwoPIGaugeScalar <: QFTdynamicsSolution
    problem::QFTdynamicsProblem
    simdata::TwoPIGaugeScalarSimData
    measurearray::Vector{Measurement}
end

function QFTdynamicsSolution(modelfile::TwoPIGaugeScalarfile, problem::QFTdynamicsProblem)
    # Get Memory
    simdata = getTwoPIGaugeSimData(problem)
    PrintMemoryofSimData!(simdata)
    #if SimPar.binning != 0 rfftwhelper = rebin(fftwhelper, SimPar.binning) else rfftwhelper=fftwhelper end # bins includes 0 bin!
    measurearray = Vector{Measurement}(undef, problem.simsetup.Nmeas+1) # includes 0 measurement
    tmpdata = TwoPIGaugeScalarTmpDataCPUfull(simdata, problem.num.nchunks)
    
    return QFTdynamicsSolutionTwoPIGaugeScalar(problem, simdata, measurearray), tmpdata
end

include("TwoPIGaugeScalar_Initialization.jl")
include("TwoPIGaugeScalar_Evolution.jl")
include("TwoPIGaugeScalar_Measurement.jl")
include("TwoPIGaugeScalar_Plotting.jl")