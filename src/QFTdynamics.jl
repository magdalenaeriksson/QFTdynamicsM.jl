module QFTdynamics

using Parameters
using FFTWhelper
# common to all types of QFTdynamics simulations
include("DefSimDataElements.jl") # defines all the Elements (SymMartix, AsymMatrix, CyclVector..) to build model specific SimData
include("DefAbstractSimSetup.jl") # defines SimConfig and createfilestructandSimConfig, relies on model specific "getParameterstring" that is used for path creation
include("AbstractTypesSimulation.jl") # defines abstract types Discretisation, Initialisation etc.; and QFTdynamicsProblem which is a structure combining all of them
include("DefPlotting.jl")

# include models
export AbstractModelfile
"""
    AbstractModelfile 
    asdl;fjslkdf;.
"""
abstract type AbstractModelfile end
#include("GenericModel/GenericModel_Modelfile.jl") 
#include("TwoPIScalar/TwoPIScalar_Modelfile.jl") 
#include("TwoPIScalar/AnalysisTools/TwoPIScalar_AnalysisTools.jl")
include("TwoPIGaugeScalar/TwoPIGaugeScalar_Modelfile.jl") 
#include("CSScalar/CSScalar_Modelfile.jl") 
#
include("CSGaugeScalar/CSGaugeScalar_Modelfile.jl") 


export fun
"""
    fun(x)
    Thats a doc string
"""
function fun(x)
    return x+1
end

end