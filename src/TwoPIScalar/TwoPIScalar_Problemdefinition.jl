################################################################################################################################
# TwoPIProblem contains different aspects of a model
# - physical Model: Parameters: lambda, ON, Mass
# - Perturbative Expansion: lambda, 1/N, and Order
# - discretization: all lattice information, including scale (M at t=0), time and fftwhelper and measurement
# - initialization: 1/2, Thermal, Tophat 
# - simsetup: restart, paths, Simtime,
#  OPEN: OPEN: OPEN: OPEN: OPEN: OPEN: OPEN: OPEN: OPEN: OPEN: OPEN: OPEN: OPEN: OPEN: OPEN: OPEN: OPEN: OPEN: OPEN: OPEN: OPEN:
# - struct Metric:  changes evolution equation (dalambert!)
#     * Minkowsi -> as is
#     * de Sitter -> Parameter H, maybe also evolved with Einstein Equation
#     * FRLW -> a(t) evolution with Friedman Equation (rho and p?)
#     * Black hole-> spherical symetric spacetime: Einstein Equation
# - Energy Renormalization -> required for solution of Einstein Equation -> substract something from T00
# - lattice type -> full or reduced lattice 
# - parallisation
################################################################################################################################
using FFTWhelper

#
# define physicsal model: Phi4
#
export TwoPIScalarModel
abstract type TwoPIScalarModel <: AbstractModel end

#
# Perturbative Expansion - selects Selfenergy and HartreeMass2
#
export TwoPIScalarPertExpansion
export TwoPIScalarLambdaLO
export TwoPIScalarLambdaNLO
export TwoPIScalarLambdaNLOquantum
export TwoPIScalarLambdaNLOclassical
abstract type TwoPIScalarPertExpansion <: AbstractPertExpansion end
abstract type TwoPIScalarLambdaEx <: TwoPIScalarPertExpansion end
struct TwoPIScalarLambdaLO <:      TwoPIScalarLambdaEx end
abstract type TwoPIScalarLambdaNLO <:     TwoPIScalarLambdaEx end
struct TwoPIScalarLambdaNLOquantum <: TwoPIScalarLambdaNLO end
struct TwoPIScalarLambdaNLOclassical <: TwoPIScalarLambdaNLO end

export TwoPIScalarNinverseLO
export TwoPIScalarNinverseNLO
abstract type TwoPIScalarNinverseEx <: TwoPIScalarPertExpansion end
struct TwoPIScalarNinverseLO <: TwoPIScalarNinverseEx  end
struct TwoPIScalarNinverseNLO <: TwoPIScalarNinverseEx  end
export TwoPIScalarNLO
export TwoPIScalarLO
TwoPIScalarNLO = Union{TwoPIScalarLambdaNLO,TwoPIScalarNinverseNLO}
TwoPIScalarLO = Union{TwoPIScalarLambdaLO,TwoPIScalarNinverseLO}
#
# Renormalization, done at t=0
#
export TwoPIScalarRenormalization
export TwoPIScalarRNone
export TwoPIScalarRMass
abstract type TwoPIScalarRenormalization <: AbstractRenormalization end
struct TwoPIScalarRNone <:           TwoPIScalarRenormalization end
struct TwoPIScalarRMass <:           TwoPIScalarRenormalization end

#
# Initialization 
#
export TwoPIScalarInitialization
export TwoPIScalarParticle
export TwoPIScalarThermal
export TwoPIScalarTopHatT1
export TwoPIScalarTopHatT2
export TwoPIScalarTopHatT3
export TwoPIScalarQuench
abstract type TwoPIScalarInitialization <: AbstractInitialization end
struct TwoPIScalarParticle <: TwoPIScalarInitialization
  n::Float64
end
struct TwoPIScalarThermal <:  TwoPIScalarInitialization
  T::Float64
end
struct TwoPIScalarTopHatT1 <: TwoPIScalarInitialization end
struct TwoPIScalarTopHatT2 <: TwoPIScalarInitialization end
struct TwoPIScalarTopHatT3 <: TwoPIScalarInitialization end
struct TwoPIScalarQuench <: TwoPIScalarInitialization end
struct TwoPIScalarTsunami <: TwoPIScalarInitialization end
struct TwoPIScalarGauss <:  TwoPIScalarInitialization
  xi::Float64
  eta::Float64
  sig::Float64
end
struct TwoPIScalarBox <: TwoPIScalarInitialization end

#
# Numerics 
#
abstract type TwoPIScalarNumerics <: AbstractNumerics end
struct TwoPIScalarCPUfull <: TwoPIScalarNumerics
  nchunks::Int
end
struct TwoPIScalarCPUred <: TwoPIScalarNumerics
  nchunks::Int
end
struct TwoPIScalarCPUred2 <: TwoPIScalarNumerics
  nchunks::Int
end
struct TwoPIScalarCPUcont <: TwoPIScalarNumerics
  nchunks::Int
end
struct TwoPIScalarGPUreduced <: TwoPIScalarNumerics
  nchunks::Int
end