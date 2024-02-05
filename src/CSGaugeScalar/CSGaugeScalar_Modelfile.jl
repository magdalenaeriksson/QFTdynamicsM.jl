using FFTWhelper
using DataFrames
using CSV

export CSGaugeScalarFile
struct CSGaugeScalarFile <: AbstractModelfile end

include("CSGaugeScalar_Problemdefinition.jl")
include("CSGaugeScalar_SimDataStruct.jl")
include("CSGaugeScalar_Discretization.jl")
include("CSGaugeScalar_ModelExpansionRenormalization.jl")
#
# Create Problem
#
function getParameterstring(modelfile::CSGaugeScalarFile, parameters)
    # Model
    if parameters["Mod"]=="CS_SUNgaugeScalar" parameterstring = "SU" * string(Int64(parameters["N"])) * "gaugeScalar"  end
    if parameters["Mod"]=="CS_U1gaugeScalar" parameterstring = "U1gaugeScalar" end
    parameterstring *= parameters["Pexp"]
    parameterstring *= "_Lambda"   * string( Int64(parameters["Lambda"] *1000 )	) 
    parameterstring *= "_g"        * string( Int64(parameters["g"] *1000      )   )
    parameterstring *= "_"          * parameters["Reno"]
    if parameters["init"]=="Thermo" 
        parameterstring *= "_T" * string( Int64(parameters["T"]*10)) 
    elseif parameters["init"]=="Pnr" 
        parameterstring *= "_n" * string( Int64(parameters["n"]*10))
    elseif parameters["init"]=="MCThermo" 
        parameterstring *= "_MCT" * string(Int64(parameters["T"]*10)) 
    elseif parameters["init"]=="MCPnr" 
        parameterstring *= "_MCn" * string(Int64(parameters["n"]*10))
    else
        parameterstring *= "_" * parameters["init"]
    end
    parameterstring *= "_Nx"        * string(parameters["Nx"]) 
    parameterstring *= "_sdim" 	    * string(parameters["sdim"]) 
    parameterstring *= "_Nsteps"    * string(parameters["Nsteps"])  
    parameterstring *= "_dt"    	* string(Int64(parameters["dt"]   * 1000)) 
    parameterstring *= "_m"     	* string(Int64(parameters["Mass"] * 100) )
    #parameterstring *= "|NstepsinMemory="    * string(parameters["NstepsinMemory"])
    parameterstring *= "_seed"     * string(Int64(parameters["seed"]))
    parameterstring *= "_Runs"    * string( Int64(parameters["Runs"])	)
    if parameters["tag"]!= "" parameterstring *= "_" * parameters["tag"] end
end

function QFTdynamicsProblem(modelfile::CSGaugeScalarFile, parameters)
    # Check parameters
    # add NstepsinMemory otherwise problem with 2PI code
    parameters["NstepsinMemory"] = 0
    if parameters["Nmeas"] == 0
        parameters["Nmeas"] = parameters["Nsteps"]
    end

    # Simsetup
    simsetup = createfilestructandSimSetup(modelfile, parameters)

    # Discretisation
    fftwhelper = getfftwhelper(parameters["Nx"], parameters["sdim"])
    
    disc = CSGaugeScalarDiscretizationLattice(
                        parameters["Mass"],
                        parameters["Nx"],
                        parameters["sdim"],
                        parameters["Nsteps"],
                        parameters["dt"],
                        parameters["Nx"]^parameters["sdim"],
                        1/(parameters["Nx"]^parameters["sdim"]),
                        parameters["Nmeas"],
                        fftwhelper,
                        [lev.deg for lev in fftwhelper])

    # Models
    if parameters["Mod"] == "CS_SUNgaugeScalar" 
        model = CS_SUNgaugeScalar( parameters["Lambda"], 
                                parameters["Mass"], 
                                parameters["g"], 
                                parameters["N"])
                                #(parameters["N"]^2-1)/(2*parameters["N"])) #CF = (N^2-1)/(2N)
    elseif parameters["Mod"] == "CS_U1gaugeScalar" 
        model = CS_U1gaugeScalar( parameters["Lambda"], parameters["Mass"], parameters["g"] )
    end

    # Initialization
    if parameters["init"] == "Thermal"  init = CSGaugeScalarThermal(parameters["T"])  end
    if parameters["init"] == "Pnr"      init = CSGaugeScalarParticle(parameters["n"]) end
    if parameters["init"] == "MCThermal"  init = CSGaugeScalarThermalMC(parameters["T"])  end
    if parameters["init"] == "MCPnr"      init = CSGaugeScalarParticleMC(parameters["n"]) end
    if parameters["init"] == "TopHatT1" init = CSGaugeScalarTopHatT1() end
    if parameters["init"] == "TopHatT2" init = CSGaugeScalarTopHatT2() end
    if parameters["init"] == "TopHatT3" init = CSGaugeScalarTopHatT3() end

    # Perturbative expansion
    pexp = CSGaugeScalarPertExpansion()

    # Renormalization
    if parameters["Reno"] == "RNone" reno = CSGaugeScalarRNone() end

    # Numerics
    if parameters["Num"] == "CPU" 
        num = CSNumericsGaugeScalarCPU(parameters["Runs"], parameters["seed"], Threads.nthreads(),  splitter( 1, parameters["Runs"], Threads.nthreads()))
    end

    # Renormalised model
    renomodel = Renormalize(model, pexp, reno, disc)

    # initialize the problem
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

export QFTdynamicsSolutionCSGaugeScalar
#mutable struct QFTdynamicsSolutionCSGaugeScalar <: QFTdynamicsSolution
mutable struct QFTdynamicsSolutionCSGaugeScalar{SU2HiggsSimData <: CSGaugeScalarSimData} <: QFTdynamicsSolution    
    problem::QFTdynamicsProblem
    #simdata::CSGaugeScalarSimData
    simdata::Vector{SU2HiggsSimData}
    measurearray::Vector{Measurement} # elements are timeslices
    measurearrayofruns::Vector{Vector{Measurement}} # collection of measurearray's for each run in vector "measurearrayofruns"
end

function QFTdynamicsSolution(modelfile::CSGaugeScalarFile, problem::QFTdynamicsProblem)

        # set simdata
        #simdata = Vector{CSGaugeScalarSimData}(undef,problem.num.Runs)
        simdata = Vector{SU2HiggsSimData}(undef,problem.num.Runs)
        for i in 1:problem.num.Runs
            simdata[i] = SU2HiggsSimData(problem.disc.Nx, problem.disc.sdim)
        end
        # set tmpdata
        #tmpdata = Vector{CSGaugeScalarTmpData}(undef,problem.num.Runs)
        tmpdata = Vector{SU2HiggsTmpData}(undef,problem.num.Runs)
        for i in 1:problem.num.Runs
            tmpdata[i] = SU2HiggsTmpData(simdata[1])
            
            for element in problem.disc.fftwhelper
                for idx in element.ind
                    tmpdata[i].k2values[idx] = element.lev2
                end
            end
        end

    # get measurearray (containing total, averaged quantities)
    measurearray = Vector{Measurement}(undef, problem.simsetup.Nmeas+1) # includes 0 measurement

    # get measurearrayofruns
    measurearrayofruns = Vector{Vector{Measurement}}(undef, problem.simsetup.Nmeas+1) # includes 0 measurement
    for i in 1:problem.simsetup.Nmeas+1
        measurearrayofruns[i] = Vector{Measurement}(undef, problem.num.Runs)
    end
    
    return QFTdynamicsSolutionCSGaugeScalar(problem, simdata, measurearray, measurearrayofruns), tmpdata
end

include("CSGaugeScalar_Initialization.jl")
include("CSGaugeScalar_Evolution.jl")
include("CSGaugeScalar_Measurement.jl")
include("CSGaugeScalar_Plotting.jl")