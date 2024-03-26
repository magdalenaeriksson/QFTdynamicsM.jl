using StatsBase

## measurement struct for single run
export MeasurementCSGaugeScalar
Base.@kwdef mutable struct MeasurementCSGaugeScalar{N} <: Measurement
	time::Float64=0
    Etot::Float64=0
    Eelec::Float64=0
    Emagn::Float64=0
    EscalKin::Float64=0
    EscalPot::Float64=0
    GaussTot::Float64=0
    GaussRel::Float64=0
    MeanPhix2::Float64=0
    Phi2k::Vector{Float64}=zeros(Float64, N )
    Phi4k::Vector{Float64}=zeros(Float64, N )
    DTk::Vector{Float64}=zeros(Float64, N )
    DLk::Vector{Float64}=zeros(Float64, N )
    E2T::Vector{Float64}=zeros(Float64, N )
    E2L::Vector{Float64}=zeros(Float64, N )
    # Scalar components
    phi2k::Vector{Vector{Float64}}=[zeros(Float64, N ) for i in 1:4]
    phi4k::Vector{Vector{Float64}}=[zeros(Float64, N ) for i in 1:4]
    pi2k::Vector{Vector{Float64}}=[zeros(Float64, N ) for i in 1:4]
    pi4k::Vector{Vector{Float64}}=[zeros(Float64, N ) for i in 1:4]
end

# add function to add single-run measurement arrays
import Base: +
function +(A::MeasurementCSGaugeScalar,B::MeasurementCSGaugeScalar)
	return MeasurementCSGaugeScalar{length(A.Phi2k)}( A.time,
				                	A.Etot     	+B.Etot     ,
				                	A.Eelec   	+B.Eelec    ,
				                	A.Emagn   	+B.Emagn   	,
				                	A.EscalKin  +B.EscalKin ,
				                	A.EscalPot  +B.EscalPot ,
				                	A.GaussTot  +B.GaussTot ,
                                    A.GaussRel  +B.GaussRel ,
                                    A.MeanPhix2 +B.MeanPhix2,
                                    A.Phi2k    .+B.Phi2k    ,
                                    A.Phi4k    .+B.Phi4k    ,
                                    A.DTk      .+B.DTk      ,
                                    A.DLk      .+B.DLk      ,
                                    A.E2T      .+B.E2T      ,
                                    A.E2L      .+B.E2L      ,
                                    # Scalar components
                                    [ A.phi2k[1] .+ B.phi2k[1], A.phi2k[2] .+ B.phi2k[2] , A.phi2k[3] .+ B.phi2k[3] , A.phi2k[4] .+ B.phi2k[4] ],
                                    [ A.phi4k[1] .+ B.phi4k[1], A.phi4k[2] .+ B.phi4k[2] , A.phi4k[3] .+ B.phi4k[3] , A.phi4k[4] .+ B.phi4k[4] ],
                                    [ A.pi2k[1]  .+ B.pi2k[1],  A.pi2k[2]  .+ B.pi2k[2]  , A.pi2k[3]  .+ B.pi2k[3]  , A.pi2k[4]  .+ B.pi2k[4]  ],
                                    [ A.pi4k[1]  .+ B.pi4k[1],  A.pi4k[2]  .+ B.pi4k[2]  , A.pi4k[3]  .+ B.pi4k[3]  , A.pi4k[4]  .+ B.pi4k[4]  ],)
end

# total (final) measurement array containing multiple averaged runs
export TotalMeasurementCSGaugeScalar
Base.@kwdef mutable struct TotalMeasurementCSGaugeScalar{N} <: Measurement
    time::Float64=0
    #
    Etot::Float64=0
    Eelec::Float64=0
    Emagn::Float64=0
    EscalKin::Float64=0
    EscalPot::Float64=0
    #
    GaussTot::Float64=0
    GaussRel::Float64=0
    # 
    MeanPhix2::Float64=0
    Phi2k::Vector{Float64}=zeros(Float64, N )
    Phi4k::Vector{Float64}=zeros(Float64, N )
    DTk::Vector{Float64}=zeros(Float64, N ) 
    DLk::Vector{Float64}=zeros(Float64, N ) 
    E2T::Vector{Float64}=zeros(Float64, N ) 
    E2L::Vector{Float64}=zeros(Float64, N ) 
	# 
    phi2k::Vector{Vector{Float64}}=[zeros(Float64, N ) for i in 1:4]
    phi4k::Vector{Vector{Float64}}=[zeros(Float64, N ) for i in 1:4]
    phi2k_err::Vector{Vector{Float64}}=[zeros(Float64, N ) for i in 1:4]
    pi2k::Vector{Vector{Float64}}=[zeros(Float64, N ) for i in 1:4]
    pi4k::Vector{Vector{Float64}}=[zeros(Float64, N ) for i in 1:4]
    pi2k_err::Vector{Vector{Float64}}=[zeros(Float64, N ) for i in 1:4]
    n::Vector{Vector{Float64}}=[zeros(Float64, N ) for i in 1:4]
    n_err::Vector{Vector{Float64}}=[zeros(Float64, N ) for i in 1:4]
    omega::Vector{Vector{Float64}}=[zeros(Float64, N ) for i in 1:4]
    omega_err::Vector{Vector{Float64}}=[zeros(Float64, N ) for i in 1:4]
    #
	Phi2k_err::Vector{Float64}=zeros(Float64, N )
end


# function to constuct total measurement array by combining and averaging single-run measurement arrays
function combinedMeasurements(model::CSGaugeScalarModel, MeasArrOfRuns::Vector{Measurement}, num::CSNumericsGaugeScalarCPU, disc::CSGaugeScalarDiscretization)
	# sum all measurments of runs
	meassum = MeasurementCSGaugeScalar{length(disc.fftwhelper)}()
	for i in 1:num.Runs
		meassum += MeasArrOfRuns[i]
	end
    # average over number of runs
    totalmeas = TotalMeasurementCSGaugeScalar{length(disc.fftwhelper)}()
    totalmeas.time      = MeasArrOfRuns[1].time
	totalmeas.Etot      = meassum.Etot      / num.Runs
    totalmeas.Eelec     = meassum.Eelec     / num.Runs
    totalmeas.Emagn     = meassum.Emagn     / num.Runs
    totalmeas.EscalKin  = meassum.EscalKin  / num.Runs
    totalmeas.EscalPot  = meassum.EscalPot  / num.Runs
    totalmeas.GaussTot  = meassum.GaussTot  / num.Runs
    totalmeas.GaussRel  = meassum.GaussRel  / num.Runs
    totalmeas.MeanPhix2 = meassum.MeanPhix2 / num.Runs
    totalmeas.Phi2k    .= meassum.Phi2k     / num.Runs
    totalmeas.Phi4k    .= meassum.Phi4k     / num.Runs
    totalmeas.DTk      .= meassum.DTk       / num.Runs
    totalmeas.DLk      .= meassum.DLk       / num.Runs
    totalmeas.E2T      .= meassum.E2T       / num.Runs
    totalmeas.E2L      .= meassum.E2L       / num.Runs
    for i in 1:4
        totalmeas.phi2k[i] .= meassum.phi2k[i]  / num.Runs
        totalmeas.phi4k[i] .= meassum.phi4k[i]  / num.Runs
        totalmeas.pi2k[i]  .= meassum.pi2k[i]  / num.Runs
        totalmeas.pi4k[i]  .= meassum.pi4k[i]  / num.Runs
    end
    # new quantities
    for i in 1:4
        # phi2k_err: first without sqrt, then check if any comp is <0, then sqrt
	    totalmeas.phi2k_err[i]  .= (totalmeas.phi4k[i] .- totalmeas.phi2k[i].^2 ) ./ (num.Runs*disc.deg)
        if sum(totalmeas.phi2k_err[i] .< 0) > 0 # if any of the components is smaller than 0
            println("Error is complex?!")
            totalmeas.phi2k_err[i]  .*= 0 # otherwise plotting does not work
        else
	        totalmeas.phi2k_err[i]  .= sqrt.( totalmeas.phi2k_err[i] )
        end
        # pi2k_err: first without sqrt, then check if any comp is <0, then sqrt
	    totalmeas.pi2k_err[i]   .= ( totalmeas.pi4k[i]  .- totalmeas.pi2k[i].^2 ) ./ (num.Runs*disc.deg)
        if sum(totalmeas.pi2k_err[i] .< 0) > 0 # if any of the components is smaller than 0
            println("Error is complex?!")
            totalmeas.pi2k_err[i]  .*= 0 # otherwise plotting does not work
        else
	        totalmeas.pi2k_err[i]  .= sqrt.( totalmeas.pi2k_err[i] )
        end
        # n: first without sqrt, then check if any comp is <0, then sqrt
        totalmeas.n[i]          .= totalmeas.phi2k[i] .* totalmeas.pi2k[i] 
        if sum(totalmeas.n[i] .< 0) > 0 # if any of the components is smaller than 0
            println("Error is complex?!")
            totalmeas.n[i]  .*= 0 # otherwise plotting does not work
        else
	        totalmeas.n[i]  .= sqrt.( totalmeas.n[i] ) .- 0.5
        end
        # n_err: first without sqrt, then check if any comp is <0, then sqrt
	    totalmeas.n_err[i]      .= 0.25 * (totalmeas.pi2k[i] ./ totalmeas.phi2k[i]) .* totalmeas.phi2k_err[i].^2 + 0.25 * (totalmeas.phi2k[i] ./ totalmeas.phi2k[i]) .* totalmeas.pi2k_err[i].^2
        if sum(totalmeas.n_err[i] .< 0) > 0 # if any of the components is smaller than 0
            println("Error is complex?!")
            totalmeas.n_err[i]  .*= 0 # otherwise plotting does not work
        else
	        totalmeas.n_err[i]  .= sqrt.( totalmeas.n_err[i] )
        end
        # omega: first without sqrt, then check if any comp is <0, then sqrt
	    totalmeas.omega[i]      .= totalmeas.phi2k[i] ./ totalmeas.phi2k[i] 
        if sum(totalmeas.omega[i] .< 0) > 0 # if any of the components is smaller than 0
            println("Error is complex?!")
            totalmeas.omega[i]  .*= 0 # otherwise plotting does not work
        else
	        totalmeas.omega[i]  .= sqrt.( totalmeas.omega[i] )
        end
	    totalmeas.omega_err[i]  .= 0.25 *  (totalmeas.omega[i] ./ totalmeas.phi2k[i]) .* totalmeas.phi2k_err[i].^2 + 0.25 * ( totalmeas.omega[i] ./ totalmeas.pi2k[i])  .* totalmeas.pi2k_err[i].^2 
        if sum(totalmeas.omega_err[i] .< 0) > 0 # if any of the components is smaller than 0
            println("Error is complex?!")
            totalmeas.omega_err[i]  .*= 0 # otherwise plotting does not work
        else
	        totalmeas.omega_err[i]  .= sqrt.( totalmeas.omega_err[i] )
        end
    end

    totalmeas.Phi2k_err = sqrt.( totalmeas.Phi4k .- totalmeas.Phi2k.^2 ) ./ (num.Runs*disc.deg)

    return totalmeas
end

export createMeasurement
function createMeasurement(t::Int64, model::CSGaugeScalarModel, simdata::CSGaugeScalarSimData, tmpdata::CSGaugeScalarTmpData, disc::CSGaugeScalarDiscretization)
    @unpack pauli, sdim, vol, Nx = simdata
    @unpack Phi_xcomp, Phi_kcomp, Phi_k, Phi2ktmp, Phi4ktmp, Eaix, Eaik, E2k, E2L, E2T, Aaix, Aaik, Dk, DL, DT, ftplan, iftplan, nvalues, P_L, P_T = tmpdata
    #
    # init measurement object
    #
	meas = MeasurementCSGaugeScalar{length(disc.fftwhelper)}()
    #
    # scalar measurements
    #
	meas.time = t * disc.dt * disc.Mass
    meas.Eelec = getEelec(model, simdata, tmpdata, disc)
    meas.Emagn = getEmagn(model, simdata, tmpdata, disc)
    meas.EscalKin = getEscalKin(model, simdata, tmpdata, disc)
    meas.EscalPot = getEscalPot(model, simdata, tmpdata, disc)
    meas.Etot = meas.Eelec + meas.Emagn + meas.EscalKin + meas.EscalPot

    if t == 2
        meas.GaussTot = sum(tmpdata.Ga2Latt)/vol
        meas.GaussRel = sum(tmpdata.GaRelLatt)/vol
    else
        meas.GaussTot = getGauss(model, simdata, tmpdata, disc)[1]
        meas.GaussRel = getGauss(model, simdata, tmpdata, disc)[2]
    end
    # 1/vol sum_x Phi(x)^dagger Phi(x)
    for idx in 1:vol
        meas.MeanPhix2 += real( tr(adjoint(simdata.Phi[idx])*simdata.Phi[idx]) )
    end
    meas.MeanPhix2 /= vol
   
    #
    # scalar propagator measurements
    #
    # Measure Pi components (using Phi storage in tmpdata) - Actually not components but small pi!
    for idx in 1:vol
        # Pi
        Phi_xcomp[1][idx] = -real(2 * simdata.Pi[idx][2]) # or real(2 * simdata.Pi[idx][3])
        Phi_xcomp[2][idx] =  imag(2 * simdata.Pi[idx][2]) # or imag(2 * simdata.Pi[idx][3])
        Phi_xcomp[3][idx] =  real(2 * simdata.Pi[idx][1]) # or real(2 * simdata.Pi[idx][4])
        Phi_xcomp[4][idx] = -imag(2 * simdata.Pi[idx][1]) # or imag(2 * simdata.Pi[idx][4])
    end
    Phi_kcomp[1] .= ftplan * Phi_xcomp[1]
    Phi_kcomp[2] .= ftplan * Phi_xcomp[2]
    Phi_kcomp[3] .= ftplan * Phi_xcomp[3]
    Phi_kcomp[4] .= ftplan * Phi_xcomp[4]
    for i in 1:length(disc.fftwhelper)
        for j in 1:disc.fftwhelper[i].deg
            idx = disc.fftwhelper[i].ind[j]
            for c in 1:4
                meas.pi2k[c][i] += real(Phi_kcomp[c][idx])^2 + imag(Phi_kcomp[c][idx])^2
                meas.pi4k[c][i] += real(Phi_kcomp[c][idx])^4 + imag(Phi_kcomp[c][idx])^4 + 2*real(Phi_kcomp[c][idx])^2 * imag(Phi_kcomp[c][idx])^2
            end
        end
        for c in 1:4
            meas.pi2k[c][i]  /= disc.fftwhelper[i].deg
            meas.pi4k[c][i]  /= disc.fftwhelper[i].deg
        end
    end
    # Normalize
    for c in 1:4
        meas.pi2k[c] /= (disc.Nx^disc.sdim) # norm = 1/sqrt(N^3) & each field gets one
        meas.pi4k[c] /= (disc.Nx^disc.sdim)^2 # norm = 1/sqrt(N^3) & each field gets one
    end
    # Measure Phi components (using Phi storage in tmpdata) - Actually not components but small phi!
    for idx in 1:vol
        # Phi
        Phi_xcomp[1][idx] = -real(2 * simdata.Phi[idx][2]) # or real(2 * simdata.Phi[idx][3])
        Phi_xcomp[2][idx] = imag(2 * simdata.Phi[idx][2]) # or imag(2 * simdata.Phi[idx][3])
        Phi_xcomp[3][idx] = real(2 * simdata.Phi[idx][1]) # or real(2 * simdata.Phi[idx][4])
        Phi_xcomp[4][idx] = -imag(2 * simdata.Phi[idx][1]) # or imag(2 * simdata.Phi[idx][4])
    end
    Phi_kcomp[1] .= ftplan * Phi_xcomp[1]
    Phi_kcomp[2] .= ftplan * Phi_xcomp[2]
    Phi_kcomp[3] .= ftplan * Phi_xcomp[3]
    Phi_kcomp[4] .= ftplan * Phi_xcomp[4]
    for i in 1:length(disc.fftwhelper)
        for j in 1:disc.fftwhelper[i].deg
            idx = disc.fftwhelper[i].ind[j]
            for c in 1:4
                meas.phi2k[c][i] += real(Phi_kcomp[c][idx])^2 + imag(Phi_kcomp[c][idx])^2
                meas.phi4k[c][i] += real(Phi_kcomp[c][idx])^4 + imag(Phi_kcomp[c][idx])^4 + 2*real(Phi_kcomp[c][idx])^2 * imag(Phi_kcomp[c][idx])^2
            end
        end
        for c in 1:4
            meas.phi2k[c][i]  /= disc.fftwhelper[i].deg
            meas.phi4k[c][i]  /= disc.fftwhelper[i].deg
        end
    end
    for c in 1:4
        meas.phi2k[c] /= (disc.Nx^disc.sdim) # norm = 1/sqrt(N^3) & each field gets one
        meas.phi4k[c] /= (disc.Nx^disc.sdim)^2 # norm = 1/sqrt(N^3) & each field gets one
    end
    # Measure Phi
    for idx in 1:vol
        Phi_xcomp[1][idx] = simdata.Phi[idx][1]
        Phi_xcomp[2][idx] = simdata.Phi[idx][2]
        Phi_xcomp[3][idx] = simdata.Phi[idx][3]
        Phi_xcomp[4][idx] = simdata.Phi[idx][4]
    end
    Phi_kcomp[1] .= ftplan * Phi_xcomp[1]
    Phi_kcomp[2] .= ftplan * Phi_xcomp[2]
    Phi_kcomp[3] .= ftplan * Phi_xcomp[3]
    Phi_kcomp[4] .= ftplan * Phi_xcomp[4]
    for idx in 1:vol
        Phi_k[idx] = SMatrix{2,2}( Phi_kcomp[1][idx], Phi_kcomp[2][idx], Phi_kcomp[3][idx], Phi_kcomp[4][idx] )
        Phi2ktmp[idx] = real(tr( adjoint(Phi_k[idx]) * Phi_k[idx] )/2)
        Phi4ktmp[idx] = real(tr( adjoint(Phi_k[idx]) * Phi_k[idx] )/4)^2 # factor of 1/4?
    end
    Phi2ktmp /= vol
    Phi4ktmp /= vol

    #
    # gauge proagator measurements
    #
    # extract gauge field in coordinate space
    for a in 1:3
        for i in 1:3
            for idx in 1:vol
                Aaix[a][i][idx] = tr( im * pauli[a] * (simdata.U[idx][i] - adjoint(simdata.U[idx][i])) ) / (2*model.g) # real data stored in complex array
                Eaix[a][i][idx] = simdata.Ea[idx][i][a]  # real data stored in complex array
            end
        end
    end
  
    # construct physical gauge field in momentum space
    for a in 1:3
        for i in 1:3
            Aaik[a][i] .= (ftplan * Aaix[a][i]) # hermitian array b(k) = b(-k)*
            Eaik[a][i] .= (ftplan * Eaix[a][i]) # hermitian array b(k) = b(-k)*

            # add phase for physical gauge field located between lattice sites:
            #for nx in 1:Nx
            #    for ny in 1:Nx
            #        for nz in 1:Nx
            #            p = [ 2*pi*nvalues[nx]/Nx, 2*pi*nvalues[ny]/Nx, 2*pi*nvalues[nz]/Nx ]
            #            Aaik[a][i][nx,ny,nz] *= exp(-im * p[i]/2)
            #            Eaik[a][i][nx,ny,nz] *= exp(-im * p[i]/2)
            #        end
            #    end
            #end
        end
    end
    
    # construct physical propagators (already color averaged)
    tmplattice = zeros(ComplexF64,3,3)
    for i in 1:Nx
        for j in 1:Nx
            for k in 1:Nx
                # set up propagator
                for idir in 1:3
                    for jdir in 1:3
                        Dk[i,j,k][idir,jdir] = Aaik[1][idir][i,j,k] * conj(Aaik[1][jdir][i,j,k])  +  Aaik[2][idir][i,j,k] * conj(Aaik[2][jdir][i,j,k])  +  Aaik[3][idir][i,j,k] * conj(Aaik[3][jdir][i,j,k]) 
                        Dk[i,j,k][idir,jdir] /= 3*vol # hermitian for every k: b(x,y) = b(y,x)*
                        E2k[i,j,k][idir,jdir] = Eaik[1][idir][i,j,k] * conj(Eaik[1][jdir][i,j,k])  +  Eaik[2][idir][i,j,k] * conj(Eaik[2][jdir][i,j,k])  +  Eaik[3][idir][i,j,k] * conj(Eaik[3][jdir][i,j,k]) 
                        E2k[i,j,k][idir,jdir] /= 3*vol # hermitian for every k: b(x,y) = b(y,x)*
                    end
                end
                # if PL contains a zero
                if sum( iszero.(P_L[i,j,k]) ) > 0 
                    DL[i,j,k]  = 0
                    E2L[i,j,k] = 0
                    tmplattice .= P_T[i,j,k] *  Dk[i,j,k] # Not hermitian
                    DT[i,j,k]  = 1/3 * real(tmplattice[1,1] + tmplattice[2,2] + tmplattice[3,3] )
                    tmplattice .= P_T[i,j,k] *  E2k[i,j,k] # Not hermitian
                    E2T[i,j,k] = 1/3 * real(tmplattice[1,1] + tmplattice[2,2] + tmplattice[3,3] )
                else
                    # compute longitudinal and transverse propagators (P_L & P_T pre-calculated in tmpdata)
                    tmplattice .= P_L[i,j,k] *  Dk[i,j,k] # Not hermitian
                    DL[i,j,k] = mean( real.( tmplattice ./ P_L[i,j,k] ) ) # if perfectly isotropic then all entries are the same
                    tmplattice .= P_T[i,j,k] *  Dk[i,j,k] # Not hermitian
                    DT[i,j,k] = mean( real.( tmplattice ./ P_T[i,j,k] ) ) # if perfectly isotropic then all entries are the same
                    tmplattice .= P_L[i,j,k] *  E2k[i,j,k] # Not hermitian
                    E2L[i,j,k] = mean( real.( tmplattice ./ P_L[i,j,k] ) ) # if perfectly isotropic then all entries are the same
                    tmplattice .= P_T[i,j,k] *  E2k[i,j,k] # Not hermitian
                    E2T[i,j,k] = mean( real.( tmplattice ./ P_T[i,j,k] ) ) # if perfectly isotropic then all entries are the same
                end
            end
        end
    end

    # Assign meas values 
    for i in 1:length(disc.fftwhelper)
        for j in 1:disc.fftwhelper[i].deg
            idx = disc.fftwhelper[i].ind[j]
            meas.Phi2k[i] += Phi2ktmp[idx]
            meas.Phi4k[i] += Phi4ktmp[idx]
            meas.DTk[i]   += DT[idx] 
            meas.DLk[i]   += DL[idx]
            meas.E2T[i]   += E2T[idx]
            meas.E2L[i]   += E2L[idx] 
        end
        meas.Phi2k[i] /= disc.fftwhelper[i].deg
        meas.Phi4k[i] /= disc.fftwhelper[i].deg
        meas.DTk[i]   /= disc.fftwhelper[i].deg
        meas.DLk[i]   /= disc.fftwhelper[i].deg
        meas.E2T[i]   /= disc.fftwhelper[i].deg
        meas.E2L[i]   /= disc.fftwhelper[i].deg 
    end
    return meas
end
    

####################################################################################################
## Calculate energy components
##
function getEelec(model::CSSUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization)
    @unpack Ea, vol = simdata
    for idx in 1:vol
        tmpdata.EelecLatt[idx] = sum(Ea[idx][1] .^2) + sum(Ea[idx][2] .^2) + sum(Ea[idx][3] .^2)
    end
    #Eelec = real( 0.5*sum(EelecLatt)/vol )
    return real( 0.5*sum(tmpdata.EelecLatt)/vol )
end
#
function getEmagn(model::CSSUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization)
    @unpack U, Nx, vol = simdata
    for x in 1:Nx
        for y in 1:Nx
            for z in 1:Nx
                tmpdata.EmagnLatt[x,y,z] = 3 - 0.5 * tr(  U[x,y,z][1] * U[mod1(x+1,Nx),y,z][2] * adjoint(U[x,mod1(y+1,Nx),z][1]) * adjoint(U[x,y,z][2])
                                                        + U[x,y,z][1] * U[mod1(x+1,Nx),y,z][3] * adjoint(U[x,y,mod1(z+1,Nx)][1]) * adjoint(U[x,y,z][3])
                                                        + U[x,y,z][2] * U[x,mod1(y+1,Nx),z][3] * adjoint(U[x,y,mod1(z+1,Nx)][2]) * adjoint(U[x,y,z][3]) )
            end
        end
    end
    #@show (4/(model.g^2) * sum(real(EmagnLatt))/vol)
    #Emagn = real( 4/(model.g^2) * sum(EmagnLatt)/vol )
    #@show sum(EmagnLatt)
    return (4/(model.g^2) * sum(real(tmpdata.EmagnLatt))/vol)
end
#
#function getEscalTot(model::CS_SUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization)
#    @unpack Pi, Phi, U, EscalLattTot, Nx, vol = simdata
#    for x in 1:Nx
#        for y in 1:Nx
#            for z in 1:Nx
#                EscalLattTot[x,y,z] = tr(adjoint(Pi[x,y,z])*Pi[x,y,z] 
#                                    + (6 + model.Mass^2 + model.Lambda*tr(adjoint(Phi[x,y,z])*Phi[x,y,z])) * adjoint(Phi[x,y,z])*Phi[x,y,z]
#                                    - 2*(   adjoint(Phi[x,y,z]) * U[x,y,z][1] * Phi[mod1(x+1,Nx),y,z] 
#                                          + adjoint(Phi[x,y,z]) * U[x,y,z][2] * Phi[x,mod1(y+1,Nx),z]
#                                          + adjoint(Phi[x,y,z]) * U[x,y,z][3] * Phi[x,y,mod1(z+1,Nx)] ) )
#            end
#        end
#    end
#    #Escal = real( sum(EscalLatt)/vol )
#    return real( sum(EscalLattTot)/vol )
#end
function getEscalPot(model::CSSUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization)
    for idx in 1:simdata.vol
        tmpdata.EscalLattPot[idx] = ( model.Mass2 + model.Lambda*tr(adjoint(simdata.Phi[idx]) * simdata.Phi[idx]) ) * tr(adjoint(simdata.Phi[idx]) * simdata.Phi[idx])
    end
    return real( sum(tmpdata.EscalLattPot)/simdata.vol )
end
function getEscalKin(model::CSSUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization)
    @unpack Pi, Phi, U, Nx, vol = simdata
    for x in 1:Nx
        for y in 1:Nx
            for z in 1:Nx
                tmpdata.EscalLattKin[x,y,z] = tr( adjoint(Pi[x,y,z]) * Pi[x,y,z] + 6 * adjoint(Phi[x,y,z]) * Phi[x,y,z]
                                            - 2*(   adjoint(Phi[x,y,z]) * U[x,y,z][1] * Phi[mod1(x+1,Nx),y,z] 
                                                  + adjoint(Phi[x,y,z]) * U[x,y,z][2] * Phi[x,mod1(y+1,Nx),z]
                                                  + adjoint(Phi[x,y,z]) * U[x,y,z][3] * Phi[x,y,mod1(z+1,Nx)] ) )
            end
        end
    end
    return real( sum(tmpdata.EscalLattKin)/vol ) 
end

#################################################################################################### 
## Calculate Gauss constraint
##
function getGauss(model::CSSUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization)
    @unpack Pi, Phi, U, Ea, Ectd, E0, pauli, Nx, vol = simdata
    @unpack Ga, GaEcont, GaRelNumerator, GaRel, Ga2Latt, GaRelLatt, rhok, chik, rhox, chix = tmpdata
    
    # compute rho
    for idx in 1:disc.vol
        rhox[1][idx] = -imag(tr(adjoint(Pi[idx]) * pauli[1] * Phi[idx]))
        rhox[2][idx] = -imag(tr(adjoint(Pi[idx]) * pauli[2] * Phi[idx]))
        rhox[3][idx] = -imag(tr(adjoint(Pi[idx]) * pauli[3] * Phi[idx]))
    end
    rhox .*= -model.g
    
    # compute E contributions
    for x in 1:Nx
        for y in 1:Nx
            for z in 1:Nx
                GaEcont[x,y,z] = real(im*tr.( pauli .* tuple(sum(Ectd[x,y,z]) - ( adjoint(U[mod1(x-1,Nx),y,z][1]) * Ectd[mod1(x-1,Nx),y,z][1] * U[mod1(x-1,Nx),y,z][1]
                                                                        +  adjoint(U[x,mod1(y-1,Nx),z][2]) * Ectd[x,mod1(y-1,Nx),z][2] * U[x,mod1(y-1,Nx),z][2]
                                                                        +  adjoint(U[x,y,mod1(z-1,Nx)][3]) * Ectd[x,y,mod1(z-1,Nx)][3] * U[x,y,mod1(z-1,Nx)][3] ))))
                Ga[x,y,z] = [ GaEcont[x,y,z][1] + rhox[1][x,y,z], GaEcont[x,y,z][2] + rhox[2][x,y,z], GaEcont[x,y,z][3] + rhox[3][x,y,z] ]
                GaRelNumerator[x,y,z] = [ GaEcont[x,y,z][1] - rhox[1][x,y,z], GaEcont[x,y,z][2] - rhox[2][x,y,z], GaEcont[x,y,z][3] - rhox[3][x,y,z] ]
                GaRel[x,y,z] = Ga[x,y,z] ./ GaRelNumerator[x,y,z]
                Ga2Latt[x,y,z] = sum(Ga[x,y,z] .^2)
                GaRelLatt[x,y,z] = sum(GaRel[x,y,z] .^2)
            end
        end
    end
    return [ sum(Ga2Latt)/vol, sum(GaRelLatt)/vol ]
    # returns total + relative vol-averaged Gauss constraint
end

####################################################################################################
####################################################################################################
export measure!
function measure!(thesolution::QFTdynamicsSolutionCSGaugeScalar, tmpdata::Vector{SU2HiggsTmpData}, t)
    @unpack problem, simdata, measurearray, measurearrayofruns = thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem


    @Threads.threads for ichunk in 1:num.threads
        for i in num.threadranges[ichunk]
            measurearrayofruns[Int64((t-2)/simsetup.NstepsbetweenNmeas) + 1][i]= createMeasurement(t, model, simdata[i], tmpdata[ichunk], disc)
        end
    end

    # combine meas from all runs into final measarr, which goes into plotting
    measurearray[Int64((t-2)/simsetup.NstepsbetweenNmeas) + 1] = combinedMeasurements(model, measurearrayofruns[Int64((t-2)/simsetup.NstepsbetweenNmeas) + 1], num, disc)
end

export measureSingle!
function measureSingle!(thesolution::QFTdynamicsSolutionCSGaugeScalar, tmpdata::CSGaugeScalarTmpData, t)
    @unpack problem, simdata, measurearray, measurearrayofruns = thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem

    measurearray[Int64((t-2)/simsetup.NstepsbetweenNmeas) + 1] = createMeasurement(t, model, simdata, tmpdata, disc)
end