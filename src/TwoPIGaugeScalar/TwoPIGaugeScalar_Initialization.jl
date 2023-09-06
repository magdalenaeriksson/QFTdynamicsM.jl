using Parameters
using FFTW
using LinearAlgebra
using DataFrames
using CSV

############################################################################################################
# Initialisation
############################################################################################################

export CSinitTmpData!
function CSinitTmpData!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization)
    # t=0: determine scalar propagator according to Gauss constraint, compute initial E-field (gauge field A=0 here)
    @unpack fftplan = simdata
    @unpack phix, piix, Pix, Phix, Pik, Phik, PiPik, PiPhik, PhiPhik, Ea, pauli, nvalues, rhox, rhok, chik, chix, Ea, Aaix, Aaik, Dk, DLk, DTk, P_L = tmpdata
    @unpack vol, Nx = disc

    ##
    location = "/Users/magdalenaeriksson/code/2PIcode/data/MCSampledIC_Nx32_sdim3_Mass100_n0_Samples90_B50_ith5_test"
    CSnbrOfRuns = 10 # this must match CS code

    for sampleNbr in 1:CSnbrOfRuns
        ## load samples
        df = CSV.read(location * "/Sample_" * string(sampleNbr) * ".csv", DataFrame)
        for idx in 1:vol
            # phi
            phix[1][idx] = df.phi1x[idx]
            phix[2][idx] = df.phi2x[idx]
            phix[3][idx] = df.phi3x[idx]
            phix[4][idx] = df.phi4x[idx]
            #
            piix[1][idx] = df.pii1x[idx]
            piix[2][idx] = df.pii2x[idx]
            piix[3][idx] = df.pii3x[idx]
            piix[4][idx] = df.pii4x[idx]
        end

        ##
        ## construct Phix
        Phi1x =  phix[3] - im * phix[4]
        Phi2x =  phix[1] + im * phix[2]
        Phi3x = -phix[1] + im * phix[2]
        Phi4x =  phix[3] + im * phix[4]
        #
        ## construct Pix
        Pi1x =  piix[3] - im * piix[4]
        Pi2x =  piix[1] + im * piix[2]
        Pi3x = -piix[1] + im * piix[2]
        Pi4x =  piix[3] + im * piix[4]
        #
        ### construct Phik
        Phi1k = fftplan * Phi1x / sqrt(vol) #fft(Phi1x) / sqrt(vol) # divide by sqrt(vol) as we want volume-averaged correlators
        Phi2k = fftplan * Phi2x / sqrt(vol) #fft(Phi2x) / sqrt(vol)
        Phi3k = fftplan * Phi3x / sqrt(vol) #fft(Phi3x) / sqrt(vol)
        Phi4k = fftplan * Phi4x / sqrt(vol) #fft(Phi4x) / sqrt(vol)
        #
        ## construct Pik
        Pi1k = fftplan * Pi1x / sqrt(vol) #fft(Pi1x) / sqrt(vol) 
        Pi2k = fftplan * Pi2x / sqrt(vol) #fft(Pi2x) / sqrt(vol)
        Pi3k = fftplan * Pi3x / sqrt(vol) #fft(Pi3x) / sqrt(vol)
        Pi4k = fftplan * Pi4x / sqrt(vol) #fft(Pi4x) / sqrt(vol)
        #
        ## construct matrix Phi and Pi
        for idx in 1:disc.vol
            Phix[idx] = SMatrix{2,2}( Phi1x[idx], Phi3x[idx], Phi2x[idx], Phi4x[idx] ) /2
            Pix[idx]  = SMatrix{2,2}( Pi1x[idx], Pi3x[idx], Pi2x[idx], Pi4x[idx] ) /2
            Phik[idx] = SMatrix{2,2}( Phi1k[idx], Phi3k[idx], Phi2k[idx], Phi4k[idx] ) /2
            Pik[idx]  = SMatrix{2,2}( Pi1k[idx], Pi3k[idx], Pi2k[idx], Pi4k[idx] ) /2
            PiPik[idx] = real(tr( adjoint(Pik[idx]) * Pik[idx] )) /2 # -> dd F -> F[t=2,tp=2]
            PiPhik[idx] = real(tr( adjoint(Pik[idx]) * Phik[idx] )) /2 # -> d F -> F[t=1,tp=0] and F[t=0,tp=1]
            PhiPhik[idx] = real(tr( adjoint(Phik[idx]) * Phik[idx] )) /2 # -> F -> F[t=0,tp=0]
        end

        # compute Gauss contraint quantities
        ## calc rho 
        #pauli4 = [1.0*im 0.0; 0.0 1.0*im]
        for idx in 1:disc.vol
            rhox[1][idx] = -imag(tr(adjoint(Pix[idx]) * pauli[1] * Phix[idx]))
            rhox[2][idx] = -imag(tr(adjoint(Pix[idx]) * pauli[2] * Phix[idx]))
            rhox[3][idx] = -imag(tr(adjoint(Pix[idx]) * pauli[3] * Phix[idx]))
            rhox[4][idx] = -imag(tr(adjoint(Pix[idx]) * im * Phix[idx]))
        end
        #@show sum(rhox[1])
        #@show sum(rhox[2])
        #@show sum(rhox[3])
        #@show sum(rhox[4])
        rhox .*= -model.g
        
        for a in 1:3
            rhok[a] = fftplan * rhox[a] #fft(rhox[a])
            #@show rhok[a][1]               # this should be zero
            chik[a] = -rhok[a] ./ simdata.k2  # real coefficients here --> chix = ifft(chik) completely real
            chik[a][1] = 0                  # set first entry to zero manually to avoid 0/0 = NaN element 
            chix[a] = real( ifft(chik[a]) )
        end
        ##
        ## compute Ea (Ea: lattice < directional index vector < colour index vector)
        for x in 1:Nx
            for y in 1:Nx
                for z in 1:Nx
                    Ea[x,y,z] = SVector(SVector( chix[1][x,y,z]-chix[1][mod1(x+1,Nx),y,z], chix[2][x,y,z]-chix[2][mod1(x+1,Nx),y,z], chix[3][x,y,z]-chix[3][mod1(x+1,Nx),y,z] ),
                                        SVector( chix[1][x,y,z]-chix[1][x,mod1(y+1,Nx),z], chix[2][x,y,z]-chix[2][x,mod1(y+1,Nx),z], chix[3][x,y,z]-chix[3][x,mod1(y+1,Nx),z] ),
                                        SVector( chix[1][x,y,z]-chix[1][x,y,mod1(z+1,Nx)], chix[2][x,y,z]-chix[2][x,y,mod1(z+1,Nx)], chix[3][x,y,z]-chix[3][x,y,mod1(z+1,Nx)] ) )
                end
            end
        end
        ##
        ## compute gauge field at t = dt
        #summ = 0.
        for a in 1:3
            for i in 1:3
                for idx in 1:vol
                    Aaix[a][i][idx] = disc.dt * Ea[idx][i][a]
                end
                #summ += sum(Aaix[a][i])
            end
        end
        #@show summ/disc.dt
        # construct physical gauge field in momentum space
        for a in 1:3
            for i in 1:3
                Aaik[a][i] .= fftplan * Aaix[a][i]

                # add phase for physical gauge field located between lattice sites:
                for nx in 1:Nx
                    for ny in 1:Nx
                        for nz in 1:Nx
                            p = [ 2*pi*nvalues[nx]/Nx, 2*pi*nvalues[ny]/Nx, 2*pi*nvalues[nz]/Nx ]
                            Aaik[a][i][nx,ny,nz] *= exp(-im * p[i]/2)
                        end
                    end
                end
            end
        end

        # construct physical propagator
        for i in 1:3
            for j in 1:3
                Dk[i,j] .= ( Aaik[1][i] .* conj(Aaik[1][i])  +  Aaik[2][i] .* conj(Aaik[2][i])  +  Aaik[3][i] .* conj(Aaik[3][i]) )
                Dk[i,j] /= 3*vol

                # compute longitudinal propagator
                DLk  .= P_L[i,j] .* Dk[i,j]
            end
        end
        # compute transverse propagators
        DTk .= Dk[1,1] + Dk[2,2] + Dk[3,3] - DLk


        #######################################################
        # add to averaged propagators
        #######################################################
        tmpdata.AvrgPhiPhik[1] += PhiPhik
        tmpdata.AvrgPiPhik[1] += PiPhik
        tmpdata.AvrgPiPik[1] += PiPik
        tmpdata.AvrgDTk[1] += DTk
        tmpdata.AvrgDLk[1] += DLk
        #######################################################

    end


    # average propagators over number of runs
    tmpdata.AvrgPhiPhik[1] ./= CSnbrOfRuns
    tmpdata.AvrgPiPhik[1] ./= CSnbrOfRuns
    tmpdata.AvrgPiPik[1] ./= CSnbrOfRuns
    tmpdata.AvrgDTk[1] ./= CSnbrOfRuns
    tmpdata.AvrgDLk[1] ./= CSnbrOfRuns

end

export CSinitSimData!
function CSinitSimData!(simdata::TwoPIGaugeScalarSimData, model::SUNgaugeScalar, tmpdata::TwoPIGaugeScalarTmpData, disc::TwoPIGaugeScalarDiscretization)
    # init simdata with CS initial conds collected in tmpdata
    expandSimData!(simdata) # expand simdataata twice since we populate up to index 2,2
    expandSimData!(simdata) # expand simdataata twice since we populate up to index 2,2
    for i in 1:length(disc.fftwhelper)
	    for j in 1:disc.fftwhelper[i].deg
            idx = disc.fftwhelper[i].ind[j]
            simdata.FS[1,1][idx] = real(tmpdata.AvrgPhiPhik[1][idx])
            simdata.FS[2,1][idx] = disc.dt*real(tmpdata.AvrgPiPhik[1][idx]) + real(tmpdata.AvrgPhiPhik[1][idx])
            simdata.FS[2,2][idx] = disc.dt^2*real(tmpdata.AvrgPiPik[1][idx]) + 2*disc.dt*real(tmpdata.AvrgPiPhik[1][idx]) + real(tmpdata.AvrgPhiPhik[1][idx])
            simdata.rS[1,1][idx] *= 0
            simdata.rS[2,1][idx] *= disc.dt * thesign(2,1, simdata.NstepsinMemory)
            simdata.rS[2,2][idx] *= 0
            #
            simdata.FT[1,1][idx] = 0
            simdata.FT[2,1][idx] = 0
            simdata.FT[2,2][idx] = disc.dt^2 * real(tmpdata.AvrgDTk[1][idx])
            simdata.rT[1,1][idx] *= 0
            simdata.rT[2,1][idx] *= disc.dt * thesign(2,1,simdata.NstepsinMemory)
            simdata.rT[2,2][idx] *= 0
            #
            simdata.FL[1,1][idx] = 0
            simdata.FL[2,1][idx] = 0
            simdata.FL[2,2][idx] = disc.dt^2 * real(tmpdata.AvrgDLk[1][idx])
            simdata.rL[1,1][idx] *= 0
            simdata.rL[2,1][idx] *= disc.dt * thesign(2,1,simdata.NstepsinMemory)
            simdata.rL[2,2][idx] *= 0
        end
    end
    
    
    setmass2values!(model, disc, simdata, tmpdata, 2)

    # supply SelfEnergies -> purely for Energycalculation
    #@show simdata.SigFS[1][1]
    #calcSelfEnergies!(model, pexp, simdata, tmpdata, disc, 2, 1) #
    #@show simdata.SigFS[1][1]
end

export getparticlenr
function getparticlenr(init::TwoPIGaugeScalarParticle, disc::TwoPIGaugeScalarDiscretization)
    return init.n .* ones(length(disc.fftwhelper))
end

export initSimData!
function initSimData!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, model::SUNgaugeScalar, pexp::TwoPIGaugeScalarLoopEx, disc::TwoPIGaugeScalarDiscretization, init::TwoPIGaugeScalarInitialization)
    # get n TwoPIGaugeScalarLoopEx
    n = getparticlenr(init, disc)
    # init simdata with n
    expandSimData!(simdata) # expand simdataata twice since we populate up to index 2,2
    expandSimData!(simdata) # expand simdataata twice since we populate up to index 2,2
    for i in 1:length(disc.fftwhelper)
        omegaS = sqrt(disc.fftwhelper[i].lev2 + disc.Mass^2)
        omegaT = sqrt(disc.fftwhelper[i].lev2 + disc.Mass^2)
        omegaL = sqrt(disc.fftwhelper[i].lev2 + disc.Mass^2)
	    for j in 1:disc.fftwhelper[i].deg
            idx = disc.fftwhelper[i].ind[j]
            simdata.FS[1,1][idx] *= (1/omegaS) * (n[i] + 0.5)
            simdata.FS[2,1][idx] *= (1/omegaS) * (n[i] + 0.5) 
            simdata.FS[2,2][idx] *= (omegaS  * (n[i] + 0.5) * disc.dt^2 + (1/omegaS) * (n[i] + 0.5))
            simdata.rS[1,1][idx] *= 0
            simdata.rS[2,1][idx] *= disc.dt * thesign(2,1,simdata.NstepsinMemory)
            simdata.rS[2,2][idx] *= 0
            simdata.FT[1,1][idx] *= (1/omegaT) * (n[i] + 0.5)
            simdata.FT[2,1][idx] *= (1/omegaT) * (n[i] + 0.5) 
            simdata.FT[2,2][idx] *= (omegaT  * (n[i] + 0.5) * disc.dt^2 + (1/omegaT) * (n[i] + 0.5))
            simdata.rT[1,1][idx] *= 0
            simdata.rT[2,1][idx] *= disc.dt * thesign(2,1,simdata.NstepsinMemory)
            simdata.rT[2,2][idx] *= 0
            simdata.FL[1,1][idx] *= (1/omegaL) * (n[i] + 0.5)
            simdata.FL[2,1][idx] *= (1/omegaL) * (n[i] + 0.5) 
            simdata.FL[2,2][idx] *= (omegaL  * (n[i] + 0.5) * disc.dt^2 + (1/omegaL) * (n[i] + 0.5))
            simdata.rL[1,1][idx] *= 0
            simdata.rL[2,1][idx] *= disc.dt * thesign(2,1,simdata.NstepsinMemory)
            simdata.rL[2,2][idx] *= 0
        end
    end
    
    setmass2values!(model, disc, simdata, tmpdata, 2)

    # supply self-energies -> purely for energy calculation
    calcSelfEnergies!(model, pexp, simdata, tmpdata, disc, 2, 1) 
   
end

export calcAngprdterm!
function calcAngprdterm!(model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, t::Int64)
    #@unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack k2, k4, k6, k2inv, fftplan, bfftplan = simdata
    @unpack G1T = tmpdata
    
    # calculate angular product terms (introduce p,q,k for pedagogical labelling/motivation)
    #p2 = k2 # represents external momentum
    #q2 = k2 # loop lattice momentum
    #k2 = k2 # k==p-q --> separately handled loop momentum
    #k4 = k4
    #p2inv = k2inv
    #q2inv = k2inv
    #k2inv = k2inv

    G1T .= simdata.FT[t,t] .- simdata.FL[t,t]
    G1T[1] = 0.

    #res .=  p2    .* (fftplan * (     (bfftplan*(G1T .* q2inv)) .* (bfftplan * ones(disc.Nx,disc.Nx,disc.Nx))  ) )
    #res .+=           fftplan * ( 2*( (bfftplan*(G1T))          .* (bfftplan * ones(disc.Nx,disc.Nx,disc.Nx))  .- ((bfftplan*(G1T .* q2inv)) .* (bfftplan*(k2))) ) )
    #res .+= p2inv .* (fftplan * (     (bfftplan*(q2 .* G1T))    .* (bfftplan * ones(disc.Nx,disc.Nx,disc.Nx))  .+ ((bfftplan*(G1T .* q2inv)) .* (bfftplan*(k4))) .- 2*(bfftplan*(G1T)) .* (bfftplan*(k2))  ) )
    #res .*= 0.25 / disc.vol^2 # 1/V for two bfft products
    # below should be equal to above since eg: fftplan * ( (bfftplan*(simdata.FT[t,t])) .* (bfftplan * ones(disc.Nx,disc.Nx,disc.Nx))) == (bfftplan*(simdata.FT[t,t]))[1]*disc.vol

    #res .=  p2    .* ( (bfftplan*(G1T .* q2inv))[1] * disc.vol )
    #res .+= ( 2*(bfftplan*(G1T))[1]*disc.vol )  .-  2*(fftplan * ((bfftplan*(G1T .* q2inv)) .* (bfftplan*(k2)) ))
    #res .+= p2inv .* ( ((bfftplan*(q2 .* G1T))[1]*disc.vol)  .+ ( fftplan * (((bfftplan*(G1T .* q2inv)) .* (bfftplan*(k4))) .- 2*(bfftplan*(G1T)) .* (bfftplan*(k2))  ) ))
    #res .*= 0.25 / disc.vol^2 # 1/V for two bfft products

    
    #tmpdata.res .= 0
    #tmpdata.res .+=  k2 .* (   (fftplan*(G1T .* k2inv))[1] * disc.vol )
    #tmpdata.res .+=        ( 2*(fftplan*(G1T))[1]*disc.vol )  .-  2*(bfftplan * ((fftplan*(G1T .* k2inv)) .* (fftplan*(k2)) ))
    #tmpdata.res .+= k2inv .* ( ((fftplan*(k2 .* G1T))[1]*disc.vol)  .+ ( bfftplan * (((fftplan*(G1T .* k2inv)) .* (fftplan*(k4))) .- 2*(fftplan*(G1T)) .* (fftplan*(k2))  ) ))
    #tmpdata.res .*= 0.25 / disc.vol^2 # 1/V for bffts

    tmpdata.res .= 0
    tmpdata.res .+=  k2 .* (   (bfftplan*(G1T .* k2inv))[1] * disc.vol )
    tmpdata.res .+=        ( 2*(bfftplan*(G1T))[1]*disc.vol )  .-  2*(fftplan * ((bfftplan*(G1T .* k2inv)) .* (bfftplan*(k2)) ))
    tmpdata.res .+= k2inv .* ( ((bfftplan*(k2 .* G1T))[1]*disc.vol)  .+ ( fftplan * (((bfftplan*(G1T .* k2inv)) .* (bfftplan*(k4))) .- 2*(bfftplan*(G1T)) .* (bfftplan*(k2))  ) ))
    tmpdata.res .*= 0.25 / disc.vol^2 # 1/V for bffts

    # exception handling
    #res .+= 0  # q=0 (zero loop momentum)
    #res .+= 0  # k=0 (loop == external momentum)
    tmpdata.res[1] = 0  # p=0 (zero external momentum)
    
end

export setmass2values!
function setmass2values!(model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, t::Int64)

    ## Higgs tree-level mass
    simdata.scalarmass2[1] = model.Mass^2

    ## Higgs Hartree contribution
    simdata.scalarmass2[1] += 8 * model.Lambda * sum( simdata.FS[t,t] ) / disc.vol 
#
    ### Mixed Hartree contributions
    simdata.scalarmass2[1] += model.g^2 * sum( (disc.sdim-1)*simdata.FT[t,t] + simdata.FL[t,t] ) / disc.vol
    simdata.transvmass2    .= model.g^2 * sum( simdata.FS[t,t] ) / disc.vol 
    simdata.longitmass2    .= model.g^2 * sum( simdata.FS[t,t] ) / disc.vol 
    
    ## Gauge Hartree contributions
    #
    calcAngprdterm!(model, disc, simdata, tmpdata, t) # sets angprdterm == tmpdata.res
    simdata.transvmass2 .= model.N * model.g^2 * sum( (disc.sdim - 2 + (disc.sdim-1)^(-1)) .* simdata.FT[t,t] .+ (1 - (disc.sdim-1)^(-1)) .* simdata.FL[t,t] )/disc.vol
    simdata.transvmass2 .+= -(model.N * model.g^2 * (disc.sdim-1)^(-1)) .* tmpdata.res
    #
    simdata.longitmass2 .= model.N * model.g^2 * sum( (disc.sdim - 2) .* simdata.FT[t,t] .+ simdata.FL[t,t] )/disc.vol 
    simdata.longitmass2 .+= (model.N * model.g^2) .* tmpdata.res
   
end

# Initialise
############################################################################################################
export initialize!
function initialize!(thesolution::QFTdynamicsSolutionTwoPIGaugeScalar, tmpdata::TwoPIGaugeScalarTmpDataCPUfull)
    @unpack problem, simdata, measurearray = thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem

    if isnewsimulation(simsetup) == true
        println("Initialising tmpdata structs: averaging samples"); flush(stdout)
        CSinitTmpData!(simdata, tmpdata, model, disc)
        println("Initialising simdata structs: initialising propagators"); flush(stdout)
        CSinitSimData!(simdata, model, tmpdata, disc)
        #initSimData!(simdata, tmpdata, model, pexp, disc, init)
    end
    #@show tmpdata.tmplatvec[1]
    simsetup.lastEvolStep = 2
    measure!(thesolution, tmpdata, 2)
end