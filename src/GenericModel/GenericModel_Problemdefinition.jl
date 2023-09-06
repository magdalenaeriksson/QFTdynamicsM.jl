#
# define physicsal model: Here I create the abstract type GenericModel, and in "GenericModel_ModelExpansionRenormalisation" I 
# create an sturcture GenericModelA, containing all parameters of this model (eg. Phi4) 
#
export AbstractGenericModel
abstract type AbstractGenericModel <: AbstractModel end

#
# Perturbative Expansion - selects Selfenergy and HartreeMass2
#
export anexpansion
struct anexpansion <: AbstractPertExpansion end

#
# Renormalization, done at t=0
#
export arenormalization
struct arenormalization <: AbstractRenormalization end

#
# Initialization 
#
export aninitialization
struct aninitialization <: AbstractInitialization end