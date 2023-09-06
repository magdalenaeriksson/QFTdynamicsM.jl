struct DiscretizationGenericModel <: AbstractDiscretization 
    Mass    ::Float64
    Nx	    ::Int64
    sdim	  ::Int64
    Nsteps  ::Int64
    dt	    ::Float64
    NstepsinMemory::Int64
  # Measurement
    Nmeas   ::Int64
  # fftwhelper
    #fftwhelper::Vector{FFTWhelper.fftwlevel}
  end