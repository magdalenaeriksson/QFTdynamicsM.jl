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
    Eaix::Vector{Vector{clattice{sdim}}} # Acttually no clattice required, it is a real field : Aaix[a][i][idx] = tr( im * pauli[a] * (simdata.U[idx][i] - adjoint(simdata.U[idx][i])) ) / (2*model.g)
    Eaik::Vector{Vector{clattice{sdim}}}
    E2k::Array{Array{ComplexF64}, 3} #Matrix
    E2L::Array{Float64, 3} 
    E2T::Array{Float64, 3}
    # Gauge field propagator helpers
    Aaix::Vector{Vector{clattice{sdim}}} # Acttually no clattice required, it is a real field : Aaix[a][i][idx] = tr( im * pauli[a] * (simdata.U[idx][i] - adjoint(simdata.U[idx][i])) ) / (2*model.g)
    Aaik::Vector{Vector{clattice{sdim}}}
    Dk::Array{Array{ComplexF64}, 3} #Matrix
    DL::Array{Float64, 3}
    DT::Array{Float64, 3}
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
    P_L::Array{Array{ComplexF64}, 3} #Matrix
    P_T::Array{Array{ComplexF64}, 3} #Matrix
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
        E2k = Array{Array{ComplexF64}}(undef,Nx,Nx,Nx)
        E2L = zeros(Nx,Nx,Nx)
        E2T = zeros(Nx,Nx,Nx)
        #
        Aaix   = Vector{Vector{clattice}}(undef,3)
        Aaik   = Vector{Vector{clattice}}(undef,3)
        Dk = Array{Array{ComplexF64}}(undef,Nx,Nx,Nx)
        DL = zeros(Nx,Nx,Nx)
        DT = zeros(Nx,Nx,Nx)

        for a in 1:3
            Eaix[a] = 0 * [ createclattice(Nx, sdim), createclattice(Nx, sdim), createclattice(Nx, sdim) ]
            Eaik[a] = 0 * [ createclattice(Nx, sdim), createclattice(Nx, sdim), createclattice(Nx, sdim) ]
            Aaix[a] = 0 * [ createclattice(Nx, sdim), createclattice(Nx, sdim), createclattice(Nx, sdim) ]
            Aaik[a] = 0 * [ createclattice(Nx, sdim), createclattice(Nx, sdim), createclattice(Nx, sdim) ]
        end
       #
        ftplan  = plan_fft(Ga2Latt)  # NB: Not in-place FFTs
        iftplan = plan_ifft(Ga2Latt)
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

        # projection operators - at each point in k space a 3x3 projection operator
        P_L = Array{Array{ComplexF64}}(undef,Nx,Nx,Nx) #sdim
        P_T = Array{Array{ComplexF64}}(undef,Nx,Nx,Nx) #sdim
        for i in 1:Nx # "in k space"
            for j in 1:Nx # "in k space"
                for k in 1:Nx # "in k space"
                    P_L[i,j,k] = zeros(ComplexF64,3,3) 
                    P_T[i,j,k] = zeros(ComplexF64,3,3)
                    E2k[i,j,k] = zeros(ComplexF64,3,3)
                    Dk[i,j,k]  = zeros(ComplexF64,3,3)
                end
            end
        end
        # momenta based on central derivative
        momentavalues = -im .* ( exp.(2*pi*im*nvalues/Nx) .- 1 ) #2*sin.(nvalues*(pi)/Nx) 
        # mometa at index 0 and Int(Nx/2) + 1 are 0 (yes, there are 2).Set them manually (otherwise the Nx/2+1 is like 10^-16)
        momentavalues[1] = 0
        #momentavalues[ Int(Nx/2) + 1 ] = 0

        for i in 1:Nx # "in k space" # x component
            for j in 1:Nx # "in k space" # y component
                for k in 1:Nx # "in k space" # y component
                    #p2 = momentavalues[i]^2 + momentavalues[j]^2 + momentavalues[k]^2
                    p2 = abs2(momentavalues[i]) + abs2(momentavalues[j]) + abs2(momentavalues[k])
                    if p2 == 0
                        P_T[i,j,k][1,1] = 1 #xx component
                        P_T[i,j,k][2,2] = 1 #yy component
                        P_T[i,j,k][3,3] = 1 #zz component
                                            
                        P_T[i,j,k][1,2] = 0 #xy component
                        P_T[i,j,k][1,3] = 0 #xz component
                        P_T[i,j,k][2,1] = 0 #yx component
                        P_T[i,j,k][2,3] = 0 #yz component
                        P_T[i,j,k][3,1] = 0 #zx component
                        P_T[i,j,k][3,2] = 0 #zy component
                
                        P_L[i,j,k][1,1] = 0 #xx component
                        P_L[i,j,k][2,2] = 0 #yy component
                        P_L[i,j,k][3,3] = 0 #zz component
                                            
                        P_L[i,j,k][1,2] = 0 #xy component
                        P_L[i,j,k][1,3] = 0 #xz component
                        P_L[i,j,k][2,1] = 0 #yx component
                        P_L[i,j,k][2,3] = 0 #yz component
                        P_L[i,j,k][3,1] = 0 #zx component
                        P_L[i,j,k][3,2] = 0 #zy component
                   else
                        # P_T
                        P_T[i,j,k][1,1] = 1 - momentavalues[i] * conj(momentavalues[i]) / p2 #xx component
                        P_T[i,j,k][2,2] = 1 - momentavalues[j] * conj(momentavalues[j]) / p2 #yy component
                        P_T[i,j,k][3,3] = 1 - momentavalues[k] * conj(momentavalues[k]) / p2 #zz component
                        P_T[i,j,k][1,2] =   - momentavalues[i] * conj(momentavalues[j]) / p2 #xy component
                        P_T[i,j,k][1,3] =   - momentavalues[i] * conj(momentavalues[k]) / p2 #xz component
                        P_T[i,j,k][2,1] =   - momentavalues[j] * conj(momentavalues[i]) / p2 #yx component
                        P_T[i,j,k][2,3] =   - momentavalues[j] * conj(momentavalues[k]) / p2 #yz component
                        P_T[i,j,k][3,1] =   - momentavalues[k] * conj(momentavalues[i]) / p2 #zx component
                        P_T[i,j,k][3,2] =   - momentavalues[k] * conj(momentavalues[j]) / p2 #zy component
                        # P_L
                        P_L[i,j,k][1,1] = momentavalues[i] * conj(momentavalues[i]) / p2 #xx component
                        P_L[i,j,k][2,2] = momentavalues[j] * conj(momentavalues[j]) / p2 #yy component
                        P_L[i,j,k][3,3] = momentavalues[k] * conj(momentavalues[k]) / p2 #zz component
                        P_L[i,j,k][1,2] = momentavalues[i] * conj(momentavalues[j]) / p2 #xy component
                        P_L[i,j,k][1,3] = momentavalues[i] * conj(momentavalues[k]) / p2 #xz component
                        P_L[i,j,k][2,1] = momentavalues[j] * conj(momentavalues[i]) / p2 #yx component
                        P_L[i,j,k][2,3] = momentavalues[j] * conj(momentavalues[k]) / p2 #yz component
                        P_L[i,j,k][3,1] = momentavalues[k] * conj(momentavalues[i]) / p2 #zx component
                        P_L[i,j,k][3,2] = momentavalues[k] * conj(momentavalues[j]) / p2 #zy component
                   end
                end
            end
        end

        return new{sdim}(   rhok, chik, rhox, chix,
                            EelecLatt, EmagnLatt, EscalLattKin, EscalLattPot,
                            Ga, GaEcont, GaRelNumerator, GaRel, Ga2Latt, GaRelLatt, 
                            Phi_xcomp, Phi_kcomp, Phi_k, Phi2ktmp, Phi4ktmp,
                            Eaix, Eaik, E2k, E2L, E2T, Aaix, Aaik, Dk, DL, DT,
                            ftplan, iftplan, plaq, plaqsum, k2values, nvalues, Nx, sdim, P_L, P_T )
    end
 end