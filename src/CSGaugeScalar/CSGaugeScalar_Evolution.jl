using FFTW
using Parameters

## get plaquette and plaquette sum
function getplaq(simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, x::Int64, y::Int64, z::Int64, i::Int64, j::Int64) # U_ij(x,y,z) with i,j = 1,2,3 being directional indices
    @unpack U, Nx = simdata
    #tmpdata.plaq .*= 0
    #P = zeros(2,2)
    if !(1 <= abs(i) <= 3 && 1 <= abs(j) <= 3 && abs(i) != abs(j))
        error("getplaq: indices either not in range or equal")
    elseif (i >= 1 && j >= 1) # postitive-postitive plaquette indices
        if (i == 1 && j == 2)
            tmpdata.plaq[1] = U[x,y,z][i] * U[mod1(x+1,Nx),y,z][j] * adjoint(U[x,mod1(y+1,Nx),z][i]) * adjoint(U[x,y,z][j])
        elseif (i == 1 && j == 3)
            tmpdata.plaq[1] = U[x,y,z][i] * U[mod1(x+1,Nx),y,z][j] * adjoint(U[x,y,mod1(z+1,Nx)][i]) * adjoint(U[x,y,z][j])
        elseif (i == 2 && j == 1)
            tmpdata.plaq[1] = U[x,y,z][i] * U[x,mod1(y+1,Nx),z][j] * adjoint(U[mod1(x+1,Nx),y,z][i]) * adjoint(U[x,y,z][j])
        elseif (i == 2 && j == 3)
            tmpdata.plaq[1] = U[x,y,z][i] * U[x,mod1(y+1,Nx),z][j] * adjoint(U[x,y,mod1(z+1,Nx)][i]) * adjoint(U[x,y,z][j])
        elseif (i == 3 && j == 1)
            tmpdata.plaq[1] = U[x,y,z][i] * U[x,y,mod1(z+1,Nx)][j] * adjoint(U[mod1(x+1,Nx),y,z][i]) * adjoint(U[x,y,z][j])
        elseif (i == 3 && j == 2)
            tmpdata.plaq[1] = U[x,y,z][i] * U[x,y,mod1(z+1,Nx)][j] * adjoint(U[x,mod1(y+1,Nx),z][i]) * adjoint(U[x,y,z][j])
        end
    elseif (i >= 1 && j < 0) # postitive-negative plaquette indices
        if (i == 1 && j == -2)
            tmpdata.plaq[1] = U[x,y,z][1] * adjoint(U[mod1(x+1,Nx),mod1(y-1,Nx),z][2]) * adjoint(U[x,mod1(y-1,Nx),z][1]) * U[x,mod1(y-1,Nx),z][2]
        elseif (i == 1 && j == -3)
            tmpdata.plaq[1] = U[x,y,z][1] * adjoint(U[mod1(x+1,Nx),y,mod1(z-1,Nx)][3]) * adjoint(U[x,y,mod1(z-1,Nx)][1]) * U[x,y,mod1(z-1,Nx)][3]
        elseif (i == 2 && j == -1)
            tmpdata.plaq[1] = U[x,y,z][2] * adjoint(U[mod1(x-1,Nx),mod1(y+1,Nx),z][1]) * adjoint(U[mod1(x-1,Nx),y,z][2]) * U[mod1(x-1,Nx),y,z][1]
        elseif (i == 2 && j == -3)
            tmpdata.plaq[1] = U[x,y,z][2] * adjoint(U[x,mod1(y+1,Nx),mod1(z-1,Nx)][3]) * adjoint(U[x,y,mod1(z-1,Nx)][2]) * U[x,y,mod1(z-1,Nx)][3]
        elseif (i == 3 && j == -1)
            tmpdata.plaq[1] = U[x,y,z][3] * adjoint(U[mod1(x-1,Nx),y,mod1(z+1,Nx)][1]) * adjoint(U[mod1(x-1,Nx),y,z][3]) * U[mod1(x-1,Nx),y,z][1]
        elseif (i == 3 && j == -2)
            tmpdata.plaq[1] = U[x,y,z][3] * adjoint(U[x,mod1(y-1,Nx),mod1(z+1,Nx)][2]) * adjoint(U[x,mod1(y-1,Nx),z][3]) * U[x,mod1(y-1,Nx),z][2]
        end
    end
    return tmpdata.plaq[1]
end
###################################################################################################################################################
##
## Compute plaquettes and plaquette sums
##
function getP(simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, x::Int64, y::Int64, z::Int64, i::Int64, j::Int64) # U_ij(x,y,z) with i,j = 1,2,3 being directional indices
    @unpack U, Nx = simdata
    if (i == j)
        tmpdata.plaq[1] = [1 0; 0 1]
    elseif (i == 1 && j == 2)
        tmpdata.plaq[1] = U[x,y,z][i] * U[mod1(x+1,Nx),y,z][j] * adjoint(U[x,mod1(y+1,Nx),z][i]) * adjoint(U[x,y,z][j])
    elseif (i == 1 && j == 3)
        tmpdata.plaq[1] = U[x,y,z][i] * U[mod1(x+1,Nx),y,z][j] * adjoint(U[x,y,mod1(z+1,Nx)][i]) * adjoint(U[x,y,z][j])
    elseif (i == 2 && j == 1)
        tmpdata.plaq[1] = U[x,y,z][i] * U[x,mod1(y+1,Nx),z][j] * adjoint(U[mod1(x+1,Nx),y,z][i]) * adjoint(U[x,y,z][j])
    elseif (i == 2 && j == 3)
        tmpdata.plaq[1] = U[x,y,z][i] * U[x,mod1(y+1,Nx),z][j] * adjoint(U[x,y,mod1(z+1,Nx)][i]) * adjoint(U[x,y,z][j])
    elseif (i == 3 && j == 1)
        tmpdata.plaq[1] = U[x,y,z][i] * U[x,y,mod1(z+1,Nx)][j] * adjoint(U[mod1(x+1,Nx),y,z][i]) * adjoint(U[x,y,z][j])
    elseif (i == 3 && j == 2)
        tmpdata.plaq[1] = U[x,y,z][i] * U[x,y,mod1(z+1,Nx)][j] * adjoint(U[x,mod1(y+1,Nx),z][i]) * adjoint(U[x,y,z][j])
    end
    return tmpdata.plaq[1]
end

function getOtherP(simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, x::Int64, y::Int64, z::Int64, i::Int64, j::Int64)
    @unpack U, Nx = simdata
    if (i == j)
        tmpdata.plaq[2] = [1 0; 0 1]
    elseif (i == 1 && j == 2)
        tmpdata.plaq[2] = U[x,y,z][i] * adjoint(U[mod1(x+1,Nx),mod1(y-1,Nx),z][j]) * adjoint(U[x,mod1(y-1,Nx),z][i]) * U[x,mod1(y-1,Nx),z][j]
    elseif (i == 1 && j == 3)
        tmpdata.plaq[2] = U[x,y,z][i] * adjoint(U[mod1(x+1,Nx),y,mod1(z-1,Nx)][j]) * adjoint(U[x,y,mod1(z-1,Nx)][i]) * U[x,y,mod1(z-1,Nx)][j]
    elseif (i == 2 && j == 1)
        tmpdata.plaq[2] = U[x,y,z][i] * adjoint(U[mod1(x-1,Nx),mod1(y+1,Nx),z][j]) * adjoint(U[mod1(x-1,Nx),y,z][i]) * U[mod1(x-1,Nx),y,z][j]
    elseif (i == 2 && j == 3)
        tmpdata.plaq[2] = U[x,y,z][i] * adjoint(U[x,mod1(y+1,Nx),mod1(z-1,Nx)][j]) * adjoint(U[x,y,mod1(z-1,Nx)][i]) * U[x,y,mod1(z-1,Nx)][j]
    elseif (i == 3 && j == 1)
        tmpdata.plaq[2] = U[x,y,z][i] * adjoint(U[mod1(x-1,Nx),y,mod1(z+1,Nx)][j]) * adjoint(U[mod1(x-1,Nx),y,z][i]) * U[mod1(x-1,Nx),y,z][j]
    elseif (i == 3 && j == 2)
        tmpdata.plaq[2] = U[x,y,z][i] * adjoint(U[x,mod1(y-1,Nx),mod1(z+1,Nx)][j]) * adjoint(U[x,mod1(y-1,Nx),z][i]) * U[x,mod1(y-1,Nx),z][j]
    end
    return tmpdata.plaq[2]
end

function getPsums(simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, x::Int64, y::Int64, z::Int64, i::Int64)
    tmpdata.plaqsum[1] *= 0
    for j in 1:3
        tmpdata.plaqsum[1] += getP(simdata,tmpdata,x,y,z,i,j)
        tmpdata.plaqsum[1] += getOtherP(simdata,tmpdata,x,y,z,i,j)
    end
    return tmpdata.plaqsum[1] # returns a 2x2 matrix 
end
###################################################################################################################################################
function getplaqsum(simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, x::Int64, y::Int64, z::Int64, i::Int64) # sum plaquettes at lattice site (x,y,z) given a directional index i=1,2,3
    #tmpdata.plaqsum .*= 0
    if !(1 <= abs(i) <= 3)
        return println("(getplaqsum): direction index not in range") 
    elseif i == 1
        tmpdata.plaqsum[1] = getplaq(simdata, tmpdata, x, y, z, i, 2) + getplaq(simdata, tmpdata, x, y, z, i, -2) + getplaq(simdata, tmpdata, x, y, z, i, 3) + getplaq(simdata, tmpdata, x, y, z, i, -3) 
    elseif i == 2
        tmpdata.plaqsum[1] = getplaq(simdata, tmpdata, x, y, z, i, 1) + getplaq(simdata, tmpdata, x, y, z, i, -1) + getplaq(simdata, tmpdata, x, y, z, i, 3) + getplaq(simdata, tmpdata, x, y, z, i, -3) 
    elseif i == 3
        tmpdata.plaqsum[1] = getplaq(simdata, tmpdata, x, y, z, i, 1) + getplaq(simdata, tmpdata, x, y, z, i, -1) + getplaq(simdata, tmpdata, x, y, z, i, 2) + getplaq(simdata, tmpdata, x, y, z, i, -2) 
    end
    return tmpdata.plaqsum[1] # returns a 2x2 matrix 
end
###############################################################################################################################################
##
## Field update/evolution functions
##
function updateEpuregauge!(model::CS_SUNgaugeScalar,simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, dt::Float64)
    @unpack Pi, Phi, U, Ea, Ectd, E0, Ta, Nx = simdata
    Ea_old = deepcopy(simdata.Ea)
    C1 = dt/model.g
    for x in 1:Nx
        for y in 1:Nx
            for z in 1:Nx
                simdata.Ea[x,y,z] = SVector( Ea_old[x,y,z][1] .+ C1*imag(tr.( Ta .* tuple(getPsums(simdata,tmpdata,x,y,z,1)) )),
                                             Ea_old[x,y,z][2] .+ C1*imag(tr.( Ta .* tuple(getPsums(simdata,tmpdata,x,y,z,2)) )),
                                             Ea_old[x,y,z][3] .+ C1*imag(tr.( Ta .* tuple(getPsums(simdata,tmpdata,x,y,z,3)) )))
            end
        end
    end
    return
end

function updateEa!(model::CS_SUNgaugeScalar,simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, dt::Float64)
    @unpack Pi, Phi, U, Ea, Ectd, E0, pauli, Nx = simdata
    Ea_old = deepcopy(simdata.Ea)
    #C1 = dt/model.g # NOTE: -im*tr(..) == +imag( tr(..) ) (Finns notation), see below 
    #C2 = dt*model.g # NOTE: see below
    for x in 1:Nx
        for y in 1:Nx
            for z in 1:Nx
                simdata.Ea[x,y,z] = SVector( Ea_old[x,y,z][1] .+ dt/model.g*imag(tr.( pauli .* tuple(getPsums(simdata,tmpdata,x,y,z,1)) )) .+ dt*model.g*imag(tr.( pauli .* tuple(U[x,y,z][1] * Phi[mod1(x+1,Nx),y,z] * adjoint(Phi[x,y,z])) )),
                                             Ea_old[x,y,z][2] .+ dt/model.g*imag(tr.( pauli .* tuple(getPsums(simdata,tmpdata,x,y,z,2)) )) .+ dt*model.g*imag(tr.( pauli .* tuple(U[x,y,z][2] * Phi[x,mod1(y+1,Nx),z] * adjoint(Phi[x,y,z])) )),
                                             Ea_old[x,y,z][3] .+ dt/model.g*imag(tr.( pauli .* tuple(getPsums(simdata,tmpdata,x,y,z,3)) )) .+ dt*model.g*imag(tr.( pauli .* tuple(U[x,y,z][3] * Phi[x,y,mod1(z+1,Nx)] * adjoint(Phi[x,y,z])) )) )
            end
        end
    end
    # convention checks of first (gauge) term and second (phi) term :
    #@show -im*(tr.( pauli .* tuple(getPsums(simdata,tmpdata,2,2,2,1)) ))                                             # Finns
    #@show +imag(tr.( pauli .* tuple(getPsums(simdata,tmpdata,2,2,2,1)) ))                                            # me == Finns (OK)
    #@show +tr.( tuple(adjoint(Phi[mod1(1+1,Nx),1,1]) * adjoint(U[1,1,1][1])) .* pauli .* tuple(im * Phi[1,1,1]))      # Finns
    #@show +tr.( pauli .* tuple(im * Phi[1,1,1] * adjoint(Phi[mod1(1+1,Nx),1,1]) * adjoint(U[1,1,1][1])) )             # Finns version 2
    #@show +imag(tr.( pauli .* tuple(U[1,1,1][1] * Phi[mod1(1+1,Nx),1,1] * adjoint(Phi[1,1,1])) ))                     # me -imag(tr(Finns version/im)) == Finns (OK)
    
    return
end
function updateEi!(model::CS_SUNgaugeScalar,simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, dt::Float64)
    @unpack Ea, Ectd, Eaσa, E2i, E0i, E2, E0, pauli, vol, Nx = simdata
    #E2lat = 0*tmpdata.flatt[1] # sum_{i,a} E(x)ia^2
    #E0lat = 0*tmpdata.flatt[1]
    #E2i = 0*Ga
    #E0i = 0*Ga
    #for idx in 1:simdata.vol
    #    #E2lat[idx] = sum( Ea[idx][1] .^2 ) + sum( Ea[idx][2] .^2 ) + sum( Ea[idx][3] .^2 )
    #    #E0lat[idx] = 1 - (model.g*dt)^2*E2lat[idx] 
    #    E2ilat[idx] = [ sum( Ea[idx][1] .^2 ), sum( Ea[idx][2] .^2 ), sum( Ea[idx][3] .^2 ) ]
    #    E0ilat[idx] = sqrt.( 1 .- (model.g*dt)^2*E2ilat[idx]/4 ) 
    #end
    #E2[1] = sum(E2lat)/vol  # 1/V * sum_{x,i,a} E(x)ia^2
    #E0[1] = sqrt( 1 - (model.g*dt)^2*E2[1] )
    #@show E0ilat[23]
    #Eaσa = 0* Ectd#Array{Vector{Matrix}}(undef,Nx,Nx,Nx)
    for idx in 1:vol
        E2i[idx] = [ sum( Ea[idx][1] .^2 ), sum( Ea[idx][2] .^2 ), sum( Ea[idx][3] .^2 ) ]
        E0i[idx] = sqrt.( 1 .- (model.g*dt)^2*E2i[idx]/4 ) 
        Eaσa[idx] = SVector( sum( Ea[idx][1] .* pauli ), sum(Ea[idx][2] .* pauli), sum(Ea[idx][3] .* pauli) )
        Ectd[idx] = (1/(model.g*dt)) * SVector( E0i[idx][1]*[1 0; 0 1] - (im*dt*model.g/2) * Eaσa[idx][1], 
                                                E0i[idx][2]*[1 0; 0 1] - (im*dt*model.g/2) * Eaσa[idx][2], 
                                                E0i[idx][3]*[1 0; 0 1] - (im*dt*model.g/2) * Eaσa[idx][3] )
    end
    #@show simdata.E0i[2][2]
    #@show simdata.Eaσa[2][2]
    #@show Ectd[54][2]*model.g*dt
    #@show simdata.U[1]
    #@show real(im*tr.(pauli .* tuple(Ectd[54][2])))
    #@show Ea[54][2]
end

function updateU!(model::CS_SUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, dt::Float64)
    U_old = deepcopy(simdata.U)
    for idx in 1:simdata.vol
        #simdata.U[idx] = SVector(model.g*dt * simdata.Ectd[idx][1] * U_old[idx][1], 
        #                         model.g*dt * simdata.Ectd[idx][2] * U_old[idx][2], 
        #                         model.g*dt * simdata.Ectd[idx][3] * U_old[idx][3] )
        simdata.U[idx] = model.g*dt .* simdata.Ectd[idx] .* U_old[idx]                  
    end
    #@show simdata.U[1]
    #@show typeof(simdata.U[1])
    return
end

function updatePhi!(simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, dt::Float64)
    Phi_old = deepcopy(simdata.Phi)
    for idx in 1:simdata.vol
        simdata.Phi[idx] = Phi_old[idx] + dt * simdata.Pi[idx]
        #simdata.PhiPhi[idx] = real(tr( adjoint(simdata.Phi[idx]) * simdata.Phi[idx] ))   # PhiPhi is Phi propagator -- a lattice
    end
    return
end

function updatePi!(model::CS_SUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, dt::Float64)
    # NOTE: double check constants in scalar potential
    @unpack Pi, Phi, U, Nx = simdata
    oldPi = deepcopy(simdata.Pi)
    for x in 1:Nx
        for y in 1:Nx
            for z in 1:Nx
                Pi[x,y,z] = oldPi[x,y,z] + dt*( - ( 6 + model.Mass^2 + 2*model.Lambda*tr(adjoint(Phi[x,y,z])*Phi[x,y,z]) ) * Phi[x,y,z]
                                                + U[x,y,z][1] * Phi[mod1(x+1,Nx),y,z] + adjoint(U[mod1(x-1,Nx),y,z][1]) * Phi[mod1(x-1,Nx),y,z] 
                                                + U[x,y,z][2] * Phi[x,mod1(y+1,Nx),z] + adjoint(U[x,mod1(y-1,Nx),z][2]) * Phi[x,mod1(y-1,Nx),z]
                                                + U[x,y,z][3] * Phi[x,y,mod1(z+1,Nx)] + adjoint(U[x,y,mod1(z-1,Nx)][3]) * Phi[x,y,mod1(z-1,Nx)] )
                #Pi[x,y,z] = oldPi[x,y,z] + dt*( -2*Phi[x,y,z] + C3*( 1/2*model.Lambda*tr(adjoint(Phi[x,y,z])*Phi[x,y,z]) + v2 )*Phi[x,y,z] 
                #                                + adjoint(U[mod1(x-1,Nx),y,z][1]) * Phi[mod1(x-1,Nx),y,z] + U[x,y,z][1] * Phi[mod1(x+1,Nx),y,z]
                #                                + adjoint(U[x,mod1(y-1,Nx),z][2]) * Phi[x,mod1(y-1,Nx),z] + U[x,y,z][2] * Phi[x,mod1(y+1,Nx),z]
                #                                + adjoint(U[x,y,mod1(z-1,Nx)][3]) * Phi[x,y,mod1(z-1,Nx)] + U[x,y,z][3] * Phi[x,y,mod1(z+1,Nx)] )
            end
        end
    end
    return
end

#
## Compute Gauss constraint
function calcGausspuregauge!(model::CS_SUNgaugeScalar,simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, dt::Float64)
    @unpack Pi, Phi, U, Ea, Ectd, E0, Ta, Nx = simdata
    for idx in 1:simdata.vol
        Ectd[idx] = SVector( sum( simdata.Ea[idx][1] .* simdata.Ta ), sum(simdata.Ea[idx][2] .* simdata.Ta), sum(simdata.Ea[idx][3] .* simdata.Ta) )
    end
    G = im*0*tmpdata.phix[1]
    for x in 1:Nx
        for y in 1:Nx
            for z in 1:Nx
                G[x,y,z] = tr( sum(Ectd[x,y,z]) - (   adjoint(U[mod1(x-1,Nx),y,z][1]) * Ectd[mod1(x-1,Nx),y,z][1] * U[mod1(x-1,Nx),y,z][1]
                                                    + adjoint(U[x,mod1(y-1,Nx),z][2]) * Ectd[x,mod1(y-1,Nx),z][2] * U[x,mod1(y-1,Nx),z][2]
                                                    + adjoint(U[x,y,mod1(z-1,Nx)][3]) * Ectd[x,y,mod1(z-1,Nx)][3] * U[x,y,mod1(z-1,Nx)][3] ))
            end
        end
    end
    @show G[1]
    @show G[132]
end
function calcGaussEa(model::CS_SUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization)
    @unpack Pi, Phi, U, Ea, Ectd, E0, pauli, Ta, Nx, vol = simdata
    @unpack Ga, rhok, chik, rhox, chix = tmpdata
    #
    ## compute rho
    for idx in 1:disc.vol
        rhox[1][idx] = -imag(tr(adjoint(Pi[idx]) * pauli[1] * Phi[idx]))
        rhox[2][idx] = -imag(tr(adjoint(Pi[idx]) * pauli[2] * Phi[idx]))
        rhox[3][idx] = -imag(tr(adjoint(Pi[idx]) * pauli[3] * Phi[idx]))
    end
    rhox .*= -model.g
    #
    ## compute E contributions
    EpartGauss1 = 0*Ga
    EpartGauss2x = 0*Ga
    EpartGauss2y = 0*Ga
    EpartGauss2z = 0*Ga
    EpartGauss = 0*Ga
    Ga2latt = 0*tmpdata.flatt[1]
    Garel = 0*Ga
    Ga2rellatt = 0*tmpdata.flatt[1]
    for x in 1:Nx
        for y in 1:Nx
            for z in 1:Nx
                EpartGauss1[x,y,z] = Ea[x,y,z][1]+Ea[x,y,z][2]+Ea[x,y,z][3]
                EpartGauss2x[x,y,z] = real(tr.( pauli .* tuple(adjoint(U[mod1(x-1,Nx),y,z][1]) * pauli[1] * U[mod1(x-1,Nx),y,z][1]) ) * Ea[mod1(x-1,Nx),y,z][1][1]
                                         + tr.( pauli .* tuple(adjoint(U[mod1(x-1,Nx),y,z][1]) * pauli[2] * U[mod1(x-1,Nx),y,z][1]) ) * Ea[mod1(x-1,Nx),y,z][1][2]
                                         + tr.( pauli .* tuple(adjoint(U[mod1(x-1,Nx),y,z][1]) * pauli[3] * U[mod1(x-1,Nx),y,z][1]) ) * Ea[mod1(x-1,Nx),y,z][1][3] )
                #EpartGauss2x[x,y,z] = 
                EpartGauss2y[x,y,z] = real(tr.( pauli .* tuple(adjoint(U[x,mod1(y-1,Nx),z][2]) * pauli[1] * U[x,mod1(y-1,Nx),z][2]) ) * Ea[x,mod1(y-1,Nx),z][2][1]
                                         + tr.( pauli .* tuple(adjoint(U[x,mod1(y-1,Nx),z][2]) * pauli[2] * U[x,mod1(y-1,Nx),z][2]) ) * Ea[x,mod1(y-1,Nx),z][2][2]
                                         + tr.( pauli .* tuple(adjoint(U[x,mod1(y-1,Nx),z][2]) * pauli[3] * U[x,mod1(y-1,Nx),z][2]) ) * Ea[x,mod1(y-1,Nx),z][2][3] )
                EpartGauss2z[x,y,z] = real(tr.( pauli .* tuple(adjoint(U[x,y,mod1(z-1,Nx)][3]) * pauli[1] * U[x,y,mod1(z-1,Nx)][3]) ) * Ea[x,y,mod1(z-1,Nx)][3][1]
                                         + tr.( pauli .* tuple(adjoint(U[x,y,mod1(z-1,Nx)][3]) * pauli[2] * U[x,y,mod1(z-1,Nx)][3]) ) * Ea[x,y,mod1(z-1,Nx)][3][2]
                                         + tr.( pauli .* tuple(adjoint(U[x,y,mod1(z-1,Nx)][3]) * pauli[3] * U[x,y,mod1(z-1,Nx)][3]) ) * Ea[x,y,mod1(z-1,Nx)][3][3] )
                EpartGauss[x,y,z] = EpartGauss1[x,y,z] - 0.5*(EpartGauss2x[x,y,z]+EpartGauss2y[x,y,z]+EpartGauss2z[x,y,z])
                Ga[x,y,z] = [ EpartGauss[x,y,z][1] + rhox[1][x,y,z], EpartGauss[x,y,z][2] + rhox[2][x,y,z], EpartGauss[x,y,z][3] + rhox[3][x,y,z] ]
                Garel[x,y,z] = [ (EpartGauss[x,y,z][1] + rhox[1][x,y,z])/(abs(EpartGauss[x,y,z][1] - rhox[1][x,y,z])), 
                                 (EpartGauss[x,y,z][2] + rhox[2][x,y,z])/(abs(EpartGauss[x,y,z][2] - rhox[2][x,y,z])),
                                 (EpartGauss[x,y,z][3] + rhox[3][x,y,z])/(abs(EpartGauss[x,y,z][3] - rhox[3][x,y,z])) ]
                Ga2latt[x,y,z] = sum(Ga[x,y,z] .^2)
                Ga2rellatt[x,y,z] = sum(Garel[x,y,z])
            end
        end
    end
    Gavg[1] = sum(Ga2latt)/vol
    Grelavg = sum(Ga2rellatt)

    #@show Gavg[1]
    #@show Grelavg[1]
    #@show Ga[223]
    #@show Ga[20]
    @show Ga[1]
    @show Ga[3]
    @show Ga[802]
    
end

function calcGaussEi(model::CS_SUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization)
    @unpack Pi, Phi, U, Ea, Ectd, E0, pauli, Ta, Nx, vol = simdata
    @unpack Ga, rhok, chik, rhox, chix = tmpdata
    #
    ## compute rho
    for idx in 1:disc.vol
        rhox[1][idx] = -imag(tr(adjoint(Pi[idx]) * pauli[1] * Phi[idx]))
        rhox[2][idx] = -imag(tr(adjoint(Pi[idx]) * pauli[2] * Phi[idx]))
        rhox[3][idx] = -imag(tr(adjoint(Pi[idx]) * pauli[3] * Phi[idx]))
    end
    rhox .*= -model.g
    #
    ## compute E contributions
    Econt = 0*Ga
    GaRelNumerator = 0*Ga
    GaRel = 0*Ga
    Ga2Latt = 0*tmpdata.flatt[1]
    GaRelLat = 0*tmpdata.flatt[1]
    #G = im*0*tmpdata.phix[1]
    ## NOTE : Econt in terms of Ectd works AFTER first time step due to definition of Ectd. Econt in terms of Ea works at all time steps. After first time step the two expressions are equivalent!
    for x in 1:Nx
        for y in 1:Nx
            for z in 1:Nx
                Econt[x,y,z] = real(im*tr.( pauli .* tuple(sum(Ectd[x,y,z]) - ( adjoint(U[mod1(x-1,Nx),y,z][1]) * Ectd[mod1(x-1,Nx),y,z][1] * U[mod1(x-1,Nx),y,z][1]
                                                                        +  adjoint(U[x,mod1(y-1,Nx),z][2]) * Ectd[x,mod1(y-1,Nx),z][2] * U[x,mod1(y-1,Nx),z][2]
                                                                        +  adjoint(U[x,y,mod1(z-1,Nx)][3]) * Ectd[x,y,mod1(z-1,Nx)][3] * U[x,y,mod1(z-1,Nx)][3] ))))
                Ga[x,y,z] = SVector( Econt[x,y,z][1] + rhox[1][x,y,z], Econt[x,y,z][2] + rhox[2][x,y,z], Econt[x,y,z][3] + rhox[3][x,y,z] )
                GaRelNumerator[x,y,z] = SVector( Econt[x,y,z][1] - rhox[1][x,y,z], Econt[x,y,z][2] - rhox[2][x,y,z], Econt[x,y,z][3] - rhox[3][x,y,z] )
                GaRel[x,y,z] = Ga[x,y,z] ./ GaRelNumerator[x,y,z]
                Ga2Latt[x,y,z] = sum(Ga[x,y,z] .^2)
                GaRelLat[x,y,z] = sum(GaRel[x,y,z] .^2)
            end
        end
    end
    @show sum(Ga2Latt)/vol
    @show sum(GaRelLat)/vol
    #return [ sum(Ga2Latt)/vol, sum(GaRelLat)/vol ]
    # return total + relative vol-averaged Gauss constraint
end

###############################################################################################################################################
# evolve-function for system
export evolveSingle!
function evolveSingle!(thesolution::QFTdynamicsSolutionCSGaugeScalar, tmpdata::CSGaugeScalarTmpData, t::Int64)
    @unpack problem, simdata, measurearray, measurearrayofruns = thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem

    #@show disc.lev
    #
    ## First update t+dt/2 functions
    updatePi!(model, simdata, tmpdata, disc.dt)
    updateEa!(model, simdata, tmpdata, disc.dt)
    #
    ## Then update t+dt functions
    updatePhi!(simdata, tmpdata, disc.dt)
    updateEi!(model, simdata, tmpdata, disc.dt)
    updateU!(model, simdata, tmpdata, disc.dt)
    #
    ## Finally compute Gauss constraint with updated fields
    #calcGaussEa(model, simdata, tmpdata, disc)
    #calcGaussEi(model, simdata, tmpdata, disc)
    #calcPropagatrs(model, simdata, tmpdata, disc)
    #
    ## Pure gauge evolution
    #updateEpuregauge!(model, simdata, tmpdata, disc.dt)
    #updateEctd!(model, simdata, tmpdata, disc.dt)
    #updateUtest!(model, simdata, disc.dt)
    #calcGausspuregauge!(model, simdata, tmpdata, disc.dt)
    
    simsetup.lastEvolStep = t 
end

export evolve!
function evolve!(thesolution::QFTdynamicsSolutionCSGaugeScalar, tmpdata::Vector{SU2HiggsTmpData}, t::Int64)
    @unpack problem, simdata, measurearray, measurearrayofruns = thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem

    @Threads.threads for ichunk in 1:num.threads
        for i in num.threadranges[ichunk]
            ## First update t+dt/2 functions
            updatePi!(model, simdata[i], tmpdata[ichunk], disc.dt)
            updateEa!(model, simdata[i], tmpdata[ichunk], disc.dt)
            #
            ## Then update t+dt functions
            updatePhi!(simdata[i], tmpdata[ichunk], disc.dt)
            updateEi!(model, simdata[i], tmpdata[ichunk], disc.dt)
            updateU!(model, simdata[i], tmpdata[ichunk], disc.dt)
        end
    end

    simsetup.lastEvolStep = t 
end

export evolveSerial!
function evolveSerial!(thesolution::QFTdynamicsSolutionCSGaugeScalar, tmpdata::Vector{SU2HiggsTmpData}, t::Int64)
    @unpack problem, simdata, measurearray, measurearrayofruns = thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem

    for i in 1:num.Runs
        ## First update t+dt/2 functions
        updatePi!(model, simdata[i], tmpdata[1], disc.dt)
        updateEa!(model, simdata[i], tmpdata[1], disc.dt)
        #
        ## Then update t+dt functions
        updatePhi!(simdata[i], tmpdata[1], disc.dt)
        updateEi!(model, simdata[i], tmpdata[1], disc.dt)
        updateU!(model, simdata[i], tmpdata[1], disc.dt)
    end

    simsetup.lastEvolStep = t 
end