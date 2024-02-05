using FFTW
using StaticArrays

#
# Define SimData 
#

#abstract type AbstractSimData already defined in general DefSimDataElements.jl 
abstract type AbstractTmpData end
abstract type CSGaugeScalarSimData <: AbstractSimData end
abstract type CSGaugeScalarTmpData <: AbstractTmpData end

export SU2HiggsSimData
struct SU2HiggsSimData <: CSGaugeScalarSimData
    # Initial real field scalar modes in x-space
    phix::Vector
    piix::Vector
    # Fields
    Pi::Array
    Phi::Array
    U::Array
    Ea::Array
    Ectd::Array
    Eaσa::Array
    E2i::Array
    E0i::Array
    E2::Vector
    E0::Vector
    # SimData helpers
    pauli::Vector
    Nx::Integer
    sdim::Integer
    vol::Integer
    #
    function SU2HiggsSimData(Nx::Int64, sdim::Int64)
        phix = 0*[createlattice(Nx, sdim),createlattice(Nx, sdim),createlattice(Nx, sdim),createlattice(Nx, sdim)]
        piix = 0*[createlattice(Nx, sdim),createlattice(Nx, sdim),createlattice(Nx, sdim),createlattice(Nx, sdim)]
        #
        SU2 = SMatrix{2,2,ComplexF64}
        SU2Object = SVector{sdim,SU2}
        colouridx = SVector{sdim,Float64}   # a number for each value of spacetime, direction and adjoint index
        EaObj = SVector{sdim,colouridx}     # each directional index carries three values of adjoint indices
        # sdim == 3
        Pi   = Array{SU2}(undef,Nx,Nx,Nx)
        Phi  = Array{SU2}(undef,Nx,Nx,Nx)
        U    = Array{SU2Object}(undef,Nx,Nx,Nx)
        Ea   = Array{EaObj}(undef,Nx,Nx,Nx)
        Ectd = Array{SU2Object}(undef,Nx,Nx,Nx)
        Eaσa = Array{SU2Object}(undef,Nx,Nx,Nx)
        E2i  = Array{Vector{Float64}}(undef,Nx,Nx,Nx)
        E0i  = Array{Vector{Float64}}(undef,Nx,Nx,Nx)
        #
        E2   = [ 0. ]
        E0   = [ 0. ]
        #
        pauli1 = SU2(0, 1, 1, 0)
        pauli2 = SU2(0, im, -im, 0)
        pauli3 = SU2(1, 0, 0, -1)
        pauli = [pauli1, pauli2, pauli3]
        #T^a = (1/2)*pauli^a

        return new(phix, piix, Pi, Phi, U, Ea, Ectd, Eaσa, E2i, E0i, E2, E0, pauli, Nx, sdim, Nx^sdim)
    end
end

 export SU2HiggsTmpData
 struct SU2HiggsTmpData{sdim} <: CSGaugeScalarTmpData
    # Gauss constraint matrices
    rhok::Vector{clattice{sdim}}
    chik::Vector{clattice{sdim}}
    rhox::Vector{clattice{sdim}}
    chix::Vector{clattice{sdim}}
    # Energies
    EelecLatt::Array
    EmagnLatt::Array
    EscalLattKin::Array
    EscalLattPot::Array
    ## Gauss constraint quantities
    Ga::Array
    GaEcont::Array
    GaRelNumerator::Array
    GaRel::Array
    Ga2Latt::clattice{sdim}
    GaRelLatt::clattice{sdim}
    ## Higgs field propagator helpers
    Phi_xcomp::Vector{clattice{sdim}}
    Phi_kcomp::Vector{clattice{sdim}}
    Phi_k::Array
    Phi2ktmp::Array
    Phi4ktmp::Array
    # Electric field propagator helpers
    Eaix::Vector{Vector{clattice{sdim}}}
    Eaik::Vector{Vector{clattice{sdim}}}
    E2k::Matrix
    E2L::clattice{sdim}
    E2T::clattice{sdim}
    # Gauge field propagator helpers
    Aaix::Vector{Vector{clattice{sdim}}}
    Aaik::Vector{Vector{clattice{sdim}}}
    Dk::Matrix
    DL::clattice{sdim}
    DT::clattice{sdim}
    ## FFT plans
    ftplan::AbstractFFTs.Plan#FFTW.cFFTWPlan
    iftplan::AbstractFFTs.Plan#AbstractFFTs.ScaledPlan
    ## Simdata helpers
    plaq::Vector{SMatrix{2,2,ComplexF64}}
    plaqsum::Vector{SMatrix{2,2,ComplexF64}}   # Matrix{ComplexF64} + return zeros(ComplexF64,2,2)
    k2values::Array
    nvalues::Vector
    Nx::Int64
    sdim::Int64
    # longitudinal/transverse projection operators
    P_L::Matrix
    P_T::Matrix
    function SU2HiggsTmpData(simdata::CSGaugeScalarSimData)
        @unpack Nx, sdim = simdata
        rhok = 0*[createclattice(Nx, sdim),createclattice(Nx, sdim),createclattice(Nx, sdim),createclattice(Nx, sdim)]
        chik = 0*[createclattice(Nx, sdim),createclattice(Nx, sdim),createclattice(Nx, sdim)]
        rhox = 0*[createlattice(Nx, sdim),createlattice(Nx, sdim),createlattice(Nx, sdim),createlattice(Nx, sdim)]
        chix = 0*[createlattice(Nx, sdim),createlattice(Nx, sdim),createlattice(Nx, sdim)]
        #
        EelecLatt    = 0*createclattice(Nx,sdim)
        EmagnLatt    = 0*createclattice(Nx,sdim)
        EscalLattKin = 0*createclattice(Nx,sdim)
        EscalLattPot = 0*createclattice(Nx,sdim)
        #
        Ga               = Array{Vector{Float64}}(undef,Nx,Nx,Nx)
        GaEcont          = Array{Vector{Float64}}(undef,Nx,Nx,Nx)
        GaRelNumerator   = Array{Vector{Float64}}(undef,Nx,Nx,Nx)
        GaRel            = Array{Vector{Float64}}(undef,Nx,Nx,Nx)
        Ga2Latt          = 0*createclattice(Nx, sdim)
        GaRelLatt        = 0*createclattice(Nx, sdim)
        #
        Phi_xcomp = 0*[createclattice(Nx, sdim),createclattice(Nx, sdim),createclattice(Nx, sdim),createclattice(Nx, sdim)]
        Phi_kcomp = 0*[createclattice(Nx, sdim),createclattice(Nx, sdim),createclattice(Nx, sdim),createclattice(Nx, sdim)]
        Phi_k     = Array{SMatrix{2,2,ComplexF64}}(undef,Nx,Nx,Nx)
        Phi2ktmp  = 0*createlattice(Nx, sdim)
        Phi4ktmp  = 0*createlattice(Nx, sdim)
        #
        Eaix   = Vector{Vector{clattice}}(undef,3)
        Eaik   = Vector{Vector{clattice}}(undef,3)
        E2k    = Matrix{clattice}(undef,3,3)
        E2L    = 0*createclattice(Nx, sdim)
        E2T    = 0*createclattice(Nx, sdim)
        #
        Aaix   = Vector{Vector{clattice}}(undef,3)
        Aaik   = Vector{Vector{clattice}}(undef,3)
        Dk     = Matrix{clattice}(undef,3,3)
        DL     = 0*createclattice(Nx, sdim)
        DT     = 0*createclattice(Nx, sdim)
        P_L    = Matrix{clattice}(undef,3,3)
        P_T    = Matrix{clattice}(undef,3,3)
        for a in 1:3
            Eaix[a] = 0 * [ createclattice(Nx, sdim), createclattice(Nx, sdim), createclattice(Nx, sdim) ]
            Eaik[a] = 0 * [ createclattice(Nx, sdim), createclattice(Nx, sdim), createclattice(Nx, sdim) ]
            Aaix[a] = 0 * [ createclattice(Nx, sdim), createclattice(Nx, sdim), createclattice(Nx, sdim) ]
            Aaik[a] = 0 * [ createclattice(Nx, sdim), createclattice(Nx, sdim), createclattice(Nx, sdim) ]
        end
        for i in 1:3
            for j in 1:3
                Dk[i,j]  = 0 * createclattice(Nx, sdim)
                E2k[i,j] = 0 * createclattice(Nx, sdim)
                P_L[i,j] = 0 * createclattice(Nx, sdim)
                P_T[i,j] = 0 * createclattice(Nx, sdim)
            end
        end
        #
        ftplan  = plan_fft(E2L) # NB: Not in-place FFTs
        iftplan = plan_ifft(E2L)
        #
        plaq     = [SMatrix{2,2}(0,0,0,0), SMatrix{2,2}(0,0,0,0)]
        plaqsum  = [SMatrix{2,2}(0,0,0,0)]
        k2values = createlattice(Nx, sdim)
        nvalues = zeros(Int64,Nx)
	    for i in 0:(Nx-1)
	        if (i>Nx/2)
	            nvalues[i+1] = -(Nx-i) # julia indexing
	        else
	            nvalues[i+1] = i  # julia indexing
	        end
	    end
        #
        # fill P_L & P_T
        for i in 1:3 # i,j directional indices
           for j in 1:3
               for nx in 1:Nx
                   for ny in 1:Nx
                       for nz in 1:Nx
                           #p_phys = [ 2*sin(pi*nvalues[nx]/Nx), 2*sin(pi*nvalues[ny]/Nx), 2*sin(pi*nvalues[nz]/Nx)]
                           #p_phys = -im * [ exp(2*pi*nvalues[nx]/Nx)-1, exp(2*pi*nvalues[ny]/Nx)-1, exp(2*pi*nvalues[nz]/Nx)-1 ]
                           p_phys = [ -im * exp(2*pi*im*(nx-1)/Nx)-1, -im * exp(2*pi*im*(ny-1)/Nx)-1, -im * exp(2*pi*im*(nz-1)/Nx)-1 ]
                           p2 = sum( abs2.(p_phys) .^2 )

                            if p2 == 0
                               P_L[i,j][nx,ny,nz] = 1
                            else
                               P_L[i,j][nx,ny,nz] =  p_phys[i] * conj(p_phys[j]) / p2
                            end
                        end
                    end
                end
                if i==j
                    P_T[i,j] .= (1 .- P_L[i,j])
                else
                    P_T[i,j] .= - P_L[i,j]
                end
           end
        end
        #@show sum(P_T[2,2] + P_L[2,2])
        #@show sum(P_T[2,1] + P_L[2,1])
        #@show sum(P_T[2,3] + P_L[2,3])
        return new{sdim}(   rhok, chik, rhox, chix,
                            EelecLatt, EmagnLatt, EscalLattKin, EscalLattPot,
                            Ga, GaEcont, GaRelNumerator, GaRel, Ga2Latt, GaRelLatt, 
                            Phi_xcomp, Phi_kcomp, Phi_k, Phi2ktmp, Phi4ktmp,
                            Eaix, Eaik, E2k, E2L, E2T, Aaix, Aaik, Dk, DL, DT,
                            ftplan, iftplan, plaq, plaqsum, k2values, nvalues, Nx, sdim, P_L, P_T )
    end
 end

