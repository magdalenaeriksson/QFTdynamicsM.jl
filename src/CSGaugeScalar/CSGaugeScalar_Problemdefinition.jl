using FFTWhelper

#
# define physicsal model: Phi4
#
export CSGaugeScalarModel
abstract type CSGaugeScalarModel <: AbstractModel end

#
# Perturbative Expansion - selects Selfenergy and HartreeMass2
#
export CSGaugeScalarPertExpansion
struct CSGaugeScalarPertExpansion <: AbstractPertExpansion end 

#
# Renormalization, done at t=0
#
export CSGaugeScalarRenormalization
export CSGaugeScalarRNone
export CSGaugeScalarRMass
abstract type CSGaugeScalarRenormalization <: AbstractRenormalization end
struct CSGaugeScalarRNone <:           CSGaugeScalarRenormalization end
struct CSGaugeScalarRMass <:           CSGaugeScalarRenormalization end

#
# Initialization 
#
export CSGaugeScalarInitialization
export CSGaugeScalarParticle
export CSGaugeScalarThermal
export CSGaugeScalarTopHatT1
export CSGaugeScalarTopHatT2
export CSGaugeScalarTopHatT3
abstract type CSGaugeScalarInitialization <: AbstractInitialization end
struct CSGaugeScalarParticle <: CSGaugeScalarInitialization
  n::Float64
end
struct CSGaugeScalarThermal <:  CSGaugeScalarInitialization
  T::Float64
end
struct CSGaugeScalarTopHatT1 <: CSGaugeScalarInitialization end
struct CSGaugeScalarTopHatT2 <: CSGaugeScalarInitialization end
struct CSGaugeScalarTopHatT3 <: CSGaugeScalarInitialization end

#
# Numerics 
#
export CSNumericsGaugeScalarCPU
export CSGaugeScalarNumerics 
abstract type CSGaugeScalarNumerics <: AbstractNumerics end
struct CSNumericsGaugeScalarCPU <: CSGaugeScalarNumerics
  Runs::Int     # number of individual runs N (total of [run#1, run#2, ..., run#N])
  seed::Int
  threads::Int  # number of threads available to do runs on 
  threadranges::Vector{UnitRange{Int64}} # range over run number/identity, ie [run#1, run#2, run#3] only
end