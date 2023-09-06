
abstract type TwoPIScalarDiscretization <: AbstractDiscretization end

struct TwoPIScalarDiscretizationLattice <: TwoPIScalarDiscretization 
  Mass            ::Float64
  Nx	            ::Int64
  sdim	          ::Int64
  Nsteps          ::Int64
  dt	            ::Float64
  vol             ::Int64
  ivol            ::Float64
  NstepsinMemory  ::Int64
  # Measurement
  Nmeas           ::Int64
  # fftwhelper
  fftwhelper      ::Vector{FFTWhelper.fftwlevel}
end

struct TwoPIScalarDiscretizationCont <: TwoPIScalarDiscretization 
  Mass            ::Float64
  Nmom	          ::Int64
  Nsteps          ::Int64
  dt	            ::Float64
  NstepsinMemory  ::Int64
  # Measurement
  Nmeas           ::Int64
  # fftwhelper
  fftwhelper      ::Vector{FFTWhelper.fftwlevel}
end