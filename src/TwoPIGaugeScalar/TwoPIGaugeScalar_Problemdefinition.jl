using FFTWhelper

#
# define umbrella physicsal model: contains SU(N) gauge-scalar and O(N) gauge-scalar theories
#
export TwoPIGaugeScalarModel
abstract type TwoPIGaugeScalarModel <: AbstractModel end


#
# Perturbative Expansion - selects Selfenergy and HartreeMass2
#
export TwoPIGaugeScalarPertExpansion
export TwoPIGaugeScalarLOloop
export TwoPIGaugeScalarNLOloop
abstract type TwoPIGaugeScalarPertExpansion <: AbstractPertExpansion end
abstract type TwoPIGaugeScalarLoopEx <: TwoPIGaugeScalarPertExpansion end
struct TwoPIGaugeScalarLOloop <:      TwoPIGaugeScalarLoopEx end
struct TwoPIGaugeScalarNLOloop <:     TwoPIGaugeScalarLoopEx end

#
# Renormalization, done at t=0
#
export TwoPIGaugeScalarRenormalization
export TwoPIGaugeScalarRNone
export TwoPIGaugeScalarRMass
abstract type TwoPIGaugeScalarRenormalization <: AbstractRenormalization end
struct TwoPIGaugeScalarRNone <:           TwoPIGaugeScalarRenormalization end
struct TwoPIGaugeScalarRMass <:           TwoPIGaugeScalarRenormalization end

#
# Initialization 
#
export TwoPIGaugeScalarInitialization
export TwoPIGaugeScalarParticle
export TwoPIGaugeScalarThermal
export TwoPIGaugeScalarTopHatT1
export TwoPIGaugeScalarTopHatT2
export TwoPIGaugeScalarTopHatT3
abstract type TwoPIGaugeScalarInitialization <: AbstractInitialization end
struct TwoPIGaugeScalarParticle <: TwoPIGaugeScalarInitialization
  n::Float64
end
struct TwoPIGaugeScalarThermal <:  TwoPIGaugeScalarInitialization
  T::Float64
end
struct TwoPIGaugeScalarTopHatT1 <: TwoPIGaugeScalarInitialization end
struct TwoPIGaugeScalarTopHatT2 <: TwoPIGaugeScalarInitialization end
struct TwoPIGaugeScalarTopHatT3 <: TwoPIGaugeScalarInitialization end

#
# Numerics 
#
abstract type TwoPIGaugeScalarNumerics <: AbstractNumerics end
struct TwoPIGaugeScalarCPUfull <: TwoPIGaugeScalarNumerics
  nchunks::Int
end
