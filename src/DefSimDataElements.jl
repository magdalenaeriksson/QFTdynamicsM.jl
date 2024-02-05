using Printf
using FFTWhelper

#"""
#This file provides atomic structures to build a custom SimData structure.
#"""
#export lattice
#export createlattice
#
#export CyclVector
#export SymMatrix

export SimData
""" 
    SimData

Partent type of all modelspecific SimData (need to be defined in corresponding modelfiles).
"""

#
# constructor for a lattice and reduced lattice
#
const lattice{sdim} = Array{Float64,sdim}
function createlattice(Nx,sdim)
    if sdim==1 return ones(Nx) end
    if sdim==2 return ones(Nx,Nx) end
    if sdim==3 return ones(Nx,Nx,Nx) end
end

export clattice
const clattice{sdim} = Array{ComplexF64,sdim}
export createclattice
function createclattice(Nx,sdim)
    if sdim==1 return ones(ComplexF64, Nx) end
    if sdim==2 return ones(ComplexF64, Nx,Nx) end
    if sdim==3 return ones(ComplexF64, Nx,Nx,Nx) end
end

#const rfftlattice{sdim} = Array{ComplexF64,sdim}
#function createrfftlattice(Nx,sdim)
#    if sdim==1 return rand(ComplexF64, Int(floor(Nx/2)+1)) end
#    if sdim==2 return rand(ComplexF64, Int(floor(Nx/2)+1),Nx) end
#    if sdim==3 return rand(ComplexF64, Int(floor(Nx/2)+1),Nx,Nx) end
#end

#const rfftlattice{sdim} = Array{Float64,sdim}
#function createrfftlattice(Nx,sdim)
#    if sdim==1 return rand(Float64, Int(floor(Nx/2)+1)) end
#    if sdim==2 return rand(Float64, Int(floor(Nx/2)+1),Nx) end
#    if sdim==3 return rand(Float64, Int(floor(Nx/2)+1),Nx,Nx) end
#end

const rfftlattice{sdim} = Array{ComplexF64,sdim}
function createrfftlattice(Nx,sdim)
    if sdim==1 return rand(ComplexF64, Int(floor(Nx/2)+1)) end
    if sdim==2 return rand(ComplexF64, Int(floor(Nx/2)+1),Nx) end
    if sdim==3 return rand(ComplexF64, Int(floor(Nx/2)+1),Nx,Nx) end
end



const reducedlattice = Vector{Float64}
function createreducedlattice(nmomenta)
    return ones(nmomenta)
end


#
# some functions to copy data bewteen lattice and reducedlattices 
#
function copytolattice!(alattice, areducedlattice::reducedlattice, fftwhelper::Vector{FFTWhelper.fftwlevel})
	for i in 1:length(fftwhelper)
	    for idx in fftwhelper[i].ind
            alattice[idx] = areducedlattice[i]
        end
    end
end

function copytoreducedlattice!(areducedlattice::reducedlattice, alattice, fftwhelper::Vector{FFTWhelper.fftwlevel})
	for i in 1:length(fftwhelper)
        areducedlattice[i] = alattice[fftwhelper[i].ind[1]]
    end
end

#
# create CyclVector
#
struct CyclVector{Sigtype}
    Sig::Vector{Sigtype}
    NstepsinMemory::Integer
    indices::Array{Int64} # [tmin, tmax]

#    function CyclVector(NstepsinMemory, Nx, sdim)
#        Sig = [ 0 .* createlattice(Nx, sdim) for i in 1:NstepsinMemory]
#        return new{eltype(Sig)}(Sig,NstepsinMemory,[0,0])
#    end
#    function CyclVector(NstepsinMemory, nmomenta)
#        Sig = [ 0 .* createreducedlattice(nmomenta) for i in 1:NstepsinMemory]
#        return new{eltype(Sig)}(Sig,NstepsinMemory,[0,0])
#    end
#    function CyclVector(NstepsinMemory)
#        return new{Float64}(zeros(NstepsinMemory),NstepsinMemory,[0,0])
#    end
end

function CyclVector(NstepsinMemory, Nx, sdim)
    Sig = [ 0 .* createlattice(Nx, sdim) for i in 1:NstepsinMemory]
    #return new{eltype(Sig)}(Sig,NstepsinMemory,[0,0])
    return CyclVector{eltype(Sig)}(Sig,NstepsinMemory,[0,0])
end
function CyclVector(NstepsinMemory, nmomenta)
    Sig = [ 0 .* createreducedlattice(nmomenta) for i in 1:NstepsinMemory]
    #return new{eltype(Sig)}(Sig,NstepsinMemory,[0,0])
    return CyclVector{eltype(Sig)}(Sig,NstepsinMemory,[0,0])
end
function CyclVector(NstepsinMemory)
    #return new{Float64}(zeros(NstepsinMemory),NstepsinMemory,[0,0])
    return CyclVector{Float64}(zeros(NstepsinMemory),NstepsinMemory,[0,0])
end
function CyclVectorx(NstepsinMemory, Nx, sdim)
    Sig = [ 0 .* createrfftlattice(Nx, sdim) for i in 1:NstepsinMemory]
    return CyclVector{eltype(Sig)}(Sig,NstepsinMemory,[0,0])
end
function CyclVectorc(NstepsinMemory, Nx, sdim)
    Sig = [ 0 .* createclattice(Nx, sdim) for i in 1:NstepsinMemory]
    return CyclVector{eltype(Sig)}(Sig,NstepsinMemory,[0,0])
end

function Base.getindex(x::CyclVector,i...)
    if !(x.indices[1] <= i[1] <= x.indices[2]) throw(BoundsError(x,"Indices not in range!") ) end
    return x.Sig[mod1(i[1],x.NstepsinMemory)]
end

function Base.setindex!(x::CyclVector,v,i)
    # Takes lattice from [i,j] and assigns the site k the value v
    # check range of indices
    if !(x.indices[1] <= i <= x.indices[2]) throw(BoundsError(x,"Indices not in range!") ) end
    x.Sig[mod1(i,x.NstepsinMemory)] = v
end

function expand!(x::CyclVector)
    if (x.indices[2] + 1) > x.NstepsinMemory
        x.indices[1] += 1
        x.indices[2] += 1
    else
        x.indices[1] = 1
        x.indices[2] += 1
    end
end

#
# Define abstract Matrix type and functions
#
abstract type TwoPIMatrix end

export expand!
function expand!(x::TwoPIMatrix)
    if (x.indices[2] + 1) > x.NstepsinMemory
        x.indices[1] += 1
        x.indices[2] += 1
    else
        x.indices[1] = 1
        x.indices[2] += 1
    end
end

# create SymMatrix type with constructor and geindex and setindex function
export SymMatrix
struct SymMatrix{Ftype} <: TwoPIMatrix
    F::Array{Array{Ftype}}
    NstepsinMemory::Integer
    indices::Array{Int64} # [tmin, tmax]
    # eventually add range object

    function SymMatrix(NstepsinMemory, Nx, sdim)
        F = [ [ createlattice(Nx, sdim) for i in 1:j] for j in 1:NstepsinMemory] # F[i][j] for i>= j is populated
        return new{eltype(eltype(F))}(F,NstepsinMemory, [0, 0])
    end
    function SymMatrix(NstepsinMemory, nmomenta)
        F = [ [ createreducedlattice(nmomenta) for i in 1:j] for j in 1:NstepsinMemory] # F[i][j] for i>= j is populated
        return new{eltype(eltype(F))}(F,NstepsinMemory, [0, 0])
    end
end

function Base.getindex(x::SymMatrix,i,j)
    # check range of indeices
    if !(x.indices[1] <= i <= x.indices[2]) throw(BoundsError(x,"Indices not in range!") ) end
    if !(x.indices[1] <= j <= x.indices[2]) throw(BoundsError(x,"Indices not in range!") ) end
    # get index 
    if mod1(i,x.NstepsinMemory) >= mod1(j,x.NstepsinMemory)
        return x.F[mod1(i,x.NstepsinMemory)][mod1(j,x.NstepsinMemory)]
    else
        return x.F[mod1(j,x.NstepsinMemory)][mod1(i,x.NstepsinMemory)]
    end
end

function Base.getindex(x::SymMatrix,i,j,k)
    # Gets lattice from [i,j] and returns value of lattice site k 
    # check range of indices
    if !(x.indices[1] <= i <= x.indices[2]) throw(BoundsError(x,"First index (" * string(i) * ") not in range!") ) end
    if !(x.indices[1] <= j <= x.indices[2]) throw(BoundsError(x,"Second index (" * string(j) * ") not in range!") ) end
    # return
    if mod1(i,x.NstepsinMemory) >= mod1(j,x.NstepsinMemory)
        return x.F[mod1(i,x.NstepsinMemory)][mod1(j,x.NstepsinMemory)][k]
    else
        return x.F[mod1(j,x.NstepsinMemory)][mod1(i,x.NstepsinMemory)][k]
    end
end

function Base.setindex!(x::SymMatrix,v,i,j,k)
    # Takes lattice from [i,j] and assigns the site k the value v
    # check range of indices
    if !(x.indices[1] <= i <= x.indices[2]) throw(BoundsError(x,"Indices not in range!") ) end
    if !(x.indices[1] <= j <= x.indices[2]) throw(BoundsError(x,"Indices not in range!") ) end
    # return
    if mod1(i,x.NstepsinMemory) >= mod1(j,x.NstepsinMemory)
        x.F[mod1(i,x.NstepsinMemory)][mod1(j,x.NstepsinMemory)][k] = v
    else
        x.F[mod1(j,x.NstepsinMemory)][mod1(i,x.NstepsinMemory)][k] = v
    end
end

export thesign
function thesign(i, j, NstepsinMemory)
    if mod1(i,NstepsinMemory) >= mod1(j,NstepsinMemory)
        return 1.
    else
        return -1.
    end
end

# create AsymMatrix type with constructor and geindex and setindex function
#struct AsymMatrix{rtype} <: TwoPIMatrix
#    r::Array{Array{rtype}}
#    NstepsinMemory::Integer
#    indices::Array{Int64} # [tmin, tmax]
#    # eventually add range object
#
#    function AsymMatrix(NstepsinMemory, Nx, sdim)
#        r = [ [ createlattice(Nx, sdim) for i in 1:j] for j in 1:NstepsinMemory] # r[i][j] for i>= j is populated
#        return new{eltype(eltype(r))}(r,NstepsinMemory, [0, 0])
#    end
#    function AsymMatrix(NstepsinMemory, nmomenta)
#        r = [ [ createreducedlattice(nmomenta) for i in 1:j] for j in 1:NstepsinMemory] # r[i][j] for i>= j is populated
#        return new{eltype(eltype(r))}(r,NstepsinMemory, [0, 0])
#    end
#end
#
#function Base.getindex(x::AsymMatrix,i,j)
#    # check range of indices
#    if !(x.indices[1] <= i <= x.indices[2]) throw(BoundsError(x,"Indices not in range!") ) end
#    if !(x.indices[1] <= j <= x.indices[2]) throw(BoundsError(x,"Indices not in range!") ) end
#    # return
#    if mod1(i,x.NstepsinMemory) >= mod1(j,x.NstepsinMemory)
#        return x.r[mod1(i,x.NstepsinMemory)][mod1(j,x.NstepsinMemory)]
#    else
#        return x.r[mod1(j,x.NstepsinMemory)][mod1(i,x.NstepsinMemory)]
#    end
#
#end

#function Base.getindex(x::AsymMatrix,i,j,k)
#    # Gets lattice from [i,j] and returns value of lattice site k 
#    # check range of indices
#    if !(x.indices[1] <= i <= x.indices[2]) throw(BoundsError(x,"Indices not in range!") ) end
#    if !(x.indices[1] <= j <= x.indices[2]) throw(BoundsError(x,"Indices not in range!") ) end
#    # return
#    if mod1(i,x.NstepsinMemory) >= mod1(j,x.NstepsinMemory)
#        return x.r[mod1(i,x.NstepsinMemory)][mod1(j,x.NstepsinMemory)][k]
#    else
#        return -x.r[mod1(j,x.NstepsinMemory)][mod1(i,x.NstepsinMemory)][k]
#    end
#end
#
#function Base.setindex!(x::AsymMatrix,v,i,j,k)
#    # Takes lattice from [i,j] and assigns the site k the value v
#    # check range of indices
#    if !(x.indices[1] <= i <= x.indices[2]) throw(BoundsError(x,"Indices not in range!") ) end
#    if !(x.indices[1] <= j <= x.indices[2]) throw(BoundsError(x,"Indices not in range!") ) end
#    # return
#    if mod1(i,x.NstepsinMemory) >= mod1(j,x.NstepsinMemory)
#        x.r[mod1(i,x.NstepsinMemory)][mod1(j,x.NstepsinMemory)][k] = v
#    else
#        x.r[mod1(j,x.NstepsinMemory)][mod1(i,x.NstepsinMemory)][k] = -v
#    end
#end

export AbstractSimData
abstract type AbstractSimData end

export PrintMemoryofSimData!
""" 
    PrintMemoryofSimData!

Bla Bla Bla.
"""

function PrintMemoryofSimData!(simdata::AbstractSimData)
    membytes = Base.summarysize(simdata)
    if round(membytes*10^(-9)) == 0
        println("Memory of simdata [MB]: ", round(membytes*10^(-6))); flush(stdout)
    else
        println("Memory of simdata [GB]: ", round( membytes*10^(-9),digits=1)); flush(stdout)
    end
end