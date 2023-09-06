
abstract type CSScalarDiscretization <: AbstractDiscretization end

export CSScalarDiscretizationLattice
struct CSScalarDiscretizationLattice <: CSScalarDiscretization 
  Mass            ::Float64
  Nx	            ::Int64
  sdim	          ::Int64
  Nsteps          ::Int64
  dt	            ::Float64
  vol             ::Int64
  ivol            ::Float64
  # Measurement
  Nmeas           ::Int64
  # fftwhelper
  fftwhelper      ::Vector{FFTWhelper.fftwlevel}
  deg             ::Vector{Int}
end