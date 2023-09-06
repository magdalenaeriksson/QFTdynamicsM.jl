abstract type CSGaugeScalarDiscretization <: AbstractDiscretization end

export CSGaugeScalarDiscretizationLattice
struct CSGaugeScalarDiscretizationLattice <: CSGaugeScalarDiscretization 
  Mass            ::Float64
  Nx	            ::Int64
  sdim	          ::Int64
  Nsteps          ::Int64
  dt	            ::Float64
  vol             ::Int64
  ivol            ::Float64
  Nmeas           ::Int64
  fftwhelper      ::Vector{FFTWhelper.fftwlevel}
  deg             ::Vector{Int}
end