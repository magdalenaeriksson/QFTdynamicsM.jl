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
    E2k::Vector{Float64}=zeros(Float64, N )
    EpropT::Vector{Float64}=zeros(Float64, N )
    EpropL::Vector{Float64}=zeros(Float64, N )
    #Phi2x::Vector{Float64}=zeros(Float64, N )
    Phi2k::Vector{Float64}=zeros(Float64, N )
    Phi4k::Vector{Float64}=zeros(Float64, N )
    DLk::Vector{Float64}=zeros(Float64, N )
    DTk::Vector{Float64}=zeros(Float64, N )
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
                                    A.E2k       +B.E2k      ,
                                    A.EpropT    +B.EpropT   ,
                                    A.EpropL    +B.EpropL   ,
                                    A.Phi2k     +B.Phi2k    ,
                                    A.Phi4k     +B.Phi4k    ,
                                    A.DLk       +B.DLk      ,
                                    A.DTk       +B.DTk      )
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
	E2k::Vector{Float64}=zeros(Float64, N )
    EpropT::Vector{Float64}=zeros(Float64, N )
    EpropL::Vector{Float64}=zeros(Float64, N )
    #Phi2x::Vector{Float64}=zeros(Float64, N )
    Phi2k::Vector{Float64}=zeros(Float64, N )
    Phi4k::Vector{Float64}=zeros(Float64, N )
    DLk::Vector{Float64}=zeros(Float64, N )
    DTk::Vector{Float64}=zeros(Float64, N )
	# 
    n::Vector{Float64}=zeros(Float64, N )
    omega::Vector{Float64}=zeros(Float64, N )
    #
	Phi2k_err::Vector{Float64}=zeros(Float64, N )
	n_err::Vector{Float64}=zeros(Float64, N )
	omega_err::Vector{Float64}=zeros(Float64, N )
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
    totalmeas.E2k      .= meassum.E2k       / num.Runs
    totalmeas.EpropT   .= meassum.EpropT    / num.Runs
    totalmeas.EpropL   .= meassum.EpropL    / num.Runs
    totalmeas.Phi2k    .= meassum.Phi2k     / num.Runs
    totalmeas.Phi4k    .= meassum.Phi4k     / num.Runs
    totalmeas.DLk      .= meassum.DLk       / num.Runs
    totalmeas.DTk      .= meassum.DTk       / num.Runs
    # new quantities
    #totalmeas.n = 
    #totalmeas.omega = 
    totalmeas.Phi2k_err = sqrt.( totalmeas.Phi4k .- totalmeas.Phi2k.^2 ) ./ (num.Runs*disc.deg)
    #totalmeas.n_err = 
    #totalmeas.omega_err = 

    return totalmeas
end

export createMeasurement
function createMeasurement(t::Int64, model::CSGaugeScalarModel, simdata::CSGaugeScalarSimData, tmpdata::CSGaugeScalarTmpData, disc::CSGaugeScalarDiscretization)
    @unpack pauli, sdim, vol, Nx = simdata
    @unpack Phi_xcomp, Phi_kcomp, Phi_k, Phi2ktmp, Phi4ktmp, Eaix, Eaik, E2k, E2L, E2T, Aaix, Aaik, Dk, DL, DT, ftplan, iftplan, nvalues, P_L, P_T = tmpdata
    # init measurement object
	meas = MeasurementCSGaugeScalar{length(disc.fftwhelper)}()
    # scalar measurements
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
    #@show meas.GaussTot
    #@show meas.GaussRel

    ######################################################################################################
    ## intial field propagator measurements
    #phikTest1 = fft.(tmpdata.phix) ./ sqrt(4*vol)
    #phikTest2 = conj.(phikTest1)
    #Phi2ktmp = phikTest1[1] .* phikTest2[1] + phikTest1[2] .* phikTest2[2] + phikTest1[3] .* phikTest2[3] + phikTest1[4] .* phikTest2[4]
    ## Phi2 from initial conditions: normalisation check
    #@show sum(tmpdata.phix[1] .^2 + tmpdata.phix[2] .^2 + tmpdata.phix[3] .^2 + tmpdata.phix[4] .^2)/vol
    #@show sum(simdata.PhiPhi)/vol
    ######################################################################################################

    ## scalar propagator measurements
    for idx in 1:vol
        Phi_xcomp[1][idx] = simdata.Phi[idx][1]
        Phi_xcomp[2][idx] = simdata.Phi[idx][2]
        Phi_xcomp[3][idx] = simdata.Phi[idx][3]
        Phi_xcomp[4][idx] = simdata.Phi[idx][4]
    end
    #Phi_kcomp .= fft.(Phi_xcomp)
    Phi_kcomp .= tuple(ftplan) .* Phi_xcomp

    for idx in 1:vol
        Phi_k[idx] = SMatrix{2,2}( Phi_kcomp[1][idx], Phi_kcomp[2][idx], Phi_kcomp[3][idx], Phi_kcomp[4][idx] )
        Phi2ktmp[idx] = real(tr( adjoint(Phi_k[idx]) * Phi_k[idx] )/2)
        Phi4ktmp[idx] = real(tr( adjoint(Phi_k[idx]) * Phi_k[idx] )/4)^2
    end
    Phi2ktmp /= vol
    Phi4ktmp /= vol

    ## gauge proagator measurements
    # extract gauge field in coordinate space
    for a in 1:3
        for i in 1:3
            for idx in 1:vol
                #Aaix[a][i][idx] = tr( im * pauli[a] * log.(simdata.U[idx][i]) ) / model.g
                #Aaix[a][i][idx] = tr( im * pauli[a] * simdata.U[idx][i] ) / model.g
                Aaix[a][i][idx] = tr( im * pauli[a] * (simdata.U[idx][i] - adjoint(simdata.U[idx][i])) ) / (2*model.g)
                Eaix[a][i][idx] = simdata.Ea[idx][i][a]
            end
        end
    end
    #@show sum(Aaix[1][1])
  
    # construct physical gauge field in momentum space
    for a in 1:3
        for i in 1:3
            Aaik[a][i] .= (ftplan * Aaix[a][i])
            Eaik[a][i] .= (ftplan * Eaix[a][i])

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
    
    # construct physical propagators
    DL .*= 0
    DT .*= 0
    E2L .*=0
    E2T .*=0
    for i in 1:3
        for j in 1:3
            Dk[i,j]  .= Aaik[1][i] .* conj(Aaik[1][j])  +  Aaik[2][i] .* conj(Aaik[2][j])  +  Aaik[3][i] .* conj(Aaik[3][j]) 
            E2k[i,j] .= Eaik[1][i] .* conj(Eaik[1][j])  +  Eaik[2][i] .* conj(Eaik[2][j])  +  Eaik[3][i] .* conj(Eaik[3][j]) 
            Dk[i,j]  ./= 3*vol
            E2k[i,j] ./= 3*vol
            # compute longitudinal and transverse propagators (P_L & P_T pre-calculated in tmpdata)
            DL  .+= (P_L[i,j] .* Dk[i,j] )
            E2L .+= (P_L[i,j] .* E2k[i,j])
            DT  .+= (P_T[i,j] .* Dk[i,j] )
            E2T .+= (P_T[i,j] .* E2k[i,j])
        end
    end

    #@show sum(Phi2ktmp)/vol
    #@show sum(P_L[1,1] * Dk[1,1])
    #@show sum(DT)/vol
    #@show sum(real(E2L))
    #@show sum(real(E2T))

    # Assign meas values 
    for i in 1:length(disc.fftwhelper)
        for j in 1:disc.fftwhelper[i].deg
            idx = disc.fftwhelper[i].ind[j]
            meas.Phi2k[i] += Phi2ktmp[idx]
            meas.Phi4k[i] += Phi4ktmp[idx]
            meas.E2k[i] += (E2k[1,1][idx]+E2k[2,2][idx]+E2k[3,3][idx])/3 # for example, just a measarr to check stuff
            meas.DLk[i] += abs( DL[idx])
            meas.DTk[i] += abs( DT[idx])
            meas.EpropL[i] += abs( E2L[idx] )
            meas.EpropT[i] += abs( E2T[idx] )
        end
        meas.Phi2k[i] /= disc.fftwhelper[i].deg
        meas.Phi4k[i] /= disc.fftwhelper[i].deg
        meas.E2k[i] /= disc.fftwhelper[i].deg
        meas.DLk[i] /= disc.fftwhelper[i].deg
        meas.DTk[i] /= disc.fftwhelper[i].deg
        meas.EpropL[i] /= disc.fftwhelper[i].deg
        meas.EpropT[i] /= disc.fftwhelper[i].deg
    end
    return meas
end
    

####################################################################################################
## Calculate energy components
##
function getEelec(model::CS_SUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization)
    @unpack Ea, vol = simdata
    for idx in 1:vol
        tmpdata.EelecLatt[idx] = sum(Ea[idx][1] .^2) + sum(Ea[idx][2] .^2) + sum(Ea[idx][3] .^2)
    end
    #Eelec = real( 0.5*sum(EelecLatt)/vol )
    return real( 0.5*sum(tmpdata.EelecLatt)/vol )
end
#
function getEmagn(model::CS_SUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization)
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
function getEscalPot(model::CS_SUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization)
    for idx in 1:simdata.vol
        tmpdata.EscalLattPot[idx] = ( model.Mass^2 + model.Lambda*tr(adjoint(simdata.Phi[idx]) * simdata.Phi[idx]) ) * tr(adjoint(simdata.Phi[idx]) * simdata.Phi[idx])
    end
    return real( sum(tmpdata.EscalLattPot)/simdata.vol )
end
function getEscalKin(model::CS_SUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization)
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
function getGauss(model::CS_SUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization)
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