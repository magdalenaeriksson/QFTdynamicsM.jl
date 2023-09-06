using FFTWhelper

#
# define physicsal model: Phi4
#
export CSScalarModel
abstract type CSScalarModel <: AbstractModel end

#
# Perturbative Expansion - selects Selfenergy and HartreeMass2
#
export CSScalarPertExpansion
struct CSScalarPertExpansion <: AbstractPertExpansion end 

#
# Renormalization, done at t=0
#
export CSScalarRenormalization
export CSScalarRNone
export CSScalarRMass
abstract type CSScalarRenormalization <: AbstractRenormalization end
struct CSScalarRNone <:           CSScalarRenormalization end
struct CSScalarRMass <:           CSScalarRenormalization end

#
# Initialization 
#
export CSScalarInitialization
export CSScalarParticle
export CSScalarThermal
export CSScalarTopHatT1
export CSScalarTopHatT2
export CSScalarTopHatT3
abstract type CSScalarInitialization <: AbstractInitialization end
struct CSScalarParticle <: CSScalarInitialization
  n::Float64
end
struct CSScalarThermal <:  CSScalarInitialization
  T::Float64
end
struct CSScalarTopHatT1 <: CSScalarInitialization end
struct CSScalarTopHatT2 <: CSScalarInitialization end
struct CSScalarTopHatT3 <: CSScalarInitialization end

#
# Numerics 
#
export CSScalarCPU
export CSScalarNumerics 
abstract type CSScalarNumerics <: AbstractNumerics end
struct CSScalarCPU <: CSScalarNumerics
  Runs::Int
  seed::Int
  threads::Int
  threadranges::Vector{UnitRange{Int64}}
end