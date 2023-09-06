
abstract type TwoPIGaugeScalarDiscretization <: AbstractDiscretization end

struct TwoPIGaugeScalarDiscretizationLattice <: TwoPIGaugeScalarDiscretization 
    Mass          ::Float64 # scalar mass
    Nx	          ::Int64
    sdim        	::Int64
    dim           ::Int64
    Nsteps        ::Int64
    dt	          ::Float64
    vol           ::Int64
    ivol          ::Float64
    NstepsinMemory::Int64
    Nmeas         ::Int64
    fftwhelper    ::Vector{FFTWhelper.fftwlevel}
  end