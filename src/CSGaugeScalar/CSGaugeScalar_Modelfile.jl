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
    if parameters["Mod"]=="CS_SUNgaugeScalar" parameterstring = "SU(" * string(Int64(parameters["N"])) * ")gaugeScalar"  end
    if parameters["Mod"]=="CS_U1gaugeScalar" parameterstring = "U(1)gaugeScalar" end
    parameterstring *= parameters["Pexp"]
    parameterstring *= "|Lambda="   * string( parameters["Lambda"]	) 
    parameterstring *= "|g="        * string( parameters["g"])
    parameterstring *= "|"          * parameters["Reno"]
    if parameters["init"]=="Thermo" 
        parameterstring *= "|T=" * string(Int64(parameters["T"])) 
    elseif parameters["init"]=="Pnr" 
        parameterstring *= "|n=" * string(Int64(parameters["n"]))
    else
        parameterstring *= "|" * parameters["init"]
    end
    parameterstring *= "|Nx="       * string(parameters["Nx"]) 
    parameterstring *= "|sdim=" 	* string(parameters["sdim"]) 
    parameterstring *= "|Nsteps="   * string(parameters["Nsteps"])  
    parameterstring *= "|dt="    	* string(parameters["dt"])  
    parameterstring *= "|m="     	* string(parameters["Mass"])
    #parameterstring *= "|NstepsinMemory="    * string(parameters["NstepsinMemory"])
    parameterstring *= "|seed="     * string(Int64(parameters["seed"]))
    parameterstring *= "|Runs="    * string( Int64(parameters["Runs"])	)
    if parameters["tag"]!= "" parameterstring *= "|" * parameters["tag"] end
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
        # load samples from csv, we need one sample for each run
        #location = "/Users/magdalenaeriksson/code/2PIcode/data/MCSampledIC_Nx32_sdim3_Mass100_n0_Samples90_B50_ith5_test"
        location = "/Users/magdalenaeriksson/code/2PIcode/data/MCSampledIC_Nx32_sdim3_Mass100_n0_Samples916_B42_ith2_test"
        for i in 1:problem.num.Runs
            simdata[i] = SU2HiggsSimData(problem.disc.Nx, problem.disc.sdim)
            df = CSV.read(location * "/Sample_" * string(i) * ".csv", DataFrame)
            for idx in 1:problem.disc.vol
                # phi
                simdata[i].phix[1][idx] = df.phi1x[idx]
                simdata[i].phix[2][idx] = df.phi2x[idx]
                simdata[i].phix[3][idx] = df.phi3x[idx]
                simdata[i].phix[4][idx] = df.phi4x[idx]
                # pi
                simdata[i].piix[1][idx] = df.pii1x[idx]
                simdata[i].piix[2][idx] = df.pii2x[idx]
                simdata[i].piix[3][idx] = df.pii3x[idx]
                simdata[i].piix[4][idx] = df.pii4x[idx]
            end
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