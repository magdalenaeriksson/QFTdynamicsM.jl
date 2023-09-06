export GenericModelSimData
export expandSimData!

#
# Define SimData 
#
struct GenericModelSimData{sdim} <: AbstractSimData
    quantitytoevolve::CyclVector{lattice{sdim}}
    #other useful objects:
    #F::SymMatrix{sdim}
    #r::AsymMatrix{sdim}
    #SigF::CyclVector{lattice{sdim}}#{sdim} #Note: lenght is NstepsinMemory (actually required is only NstepsinMemory-1, but this way the modulo operation is wrong)
    #tmplattice::Vector{lattice{sdim}}
    #tmpscalar::Vector{Float64}

    NstepsinMemory::Integer
    Nx::Integer
    sdim::Integer
    indices::Array{Int64} # [tmin, tmax]

    function GenericModelSimData(NstepsinMemory, Nx, sdim)
        return new{ sdim }(CyclVector(NstepsinMemory, Nx, sdim), NstepsinMemory, Nx, sdim, [0,0])
    end
end

function expandSimData!(x::GenericModelSimData)
    # change indices list of SimData
    if (x.indices[2] + 1) > x.NstepsinMemory
        x.indices[1] += 1
        x.indices[2] += 1
    else
        x.indices[1] = 1
        x.indices[2] += 1
    end
    # expand individual objects
    expand!(x.quantitytoevolve)
end