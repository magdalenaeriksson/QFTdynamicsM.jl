export getangprdterm
function getangprdterm(model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, time::Int64)
    @unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack res, G1T = tmpdata
    
    # calculate angular product terms (introduce p,q,k for pedagogical labelling/motivation)
    p2 = k2values # represents external momentum
    q2 = k2values # loop lattice momentum
    k2 = k2values # k==p-q --> separately handled loop momentum
    k4 = k4values
    p2inv = invk2values
    q2inv = invk2values
    k2inv = invk2values

    
    #vec = [createclattice(disc.Nx,disc.sdim), createclattice(disc.Nx,disc.sdim), createclattice(disc.Nx,disc.sdim)]
    #vec[1] .= (bfftplan*(fun .* q2inv)) * (disc.vol) 
    #vec[2] .= 2*( (bfftplan*(fun)) * (disc.vol) - (bfftplan*(fun .* q2inv)) .* (bfftplan*(k2)) )
    #vec[3] .= (bfftplan*(q2 .* fun)) * (disc.vol) + (bfftplan*(fun .* k2inv)) .* (bfftplan*(k4)) + (-2)*(bfftplan*(fun)) .* (bfftplan*(k2))
    #angprdterm = 0.25 * real( p2 .* (fftplan*vec[1]) + fftplan*vec[2] + p2inv .* (fftplan*vec[3]) ) / disc.vol^3 # 1/V for two bfft products, and one 1/V for loop summation

    #vec = 0*im*[p2,p2,p2]
    #fun = simdata.FT[time,time] - simdata.FL[time,time]
    #fun[1] = 0. # this is not true in general, but we want to construct tmpvec such that the zero loop momentum does not contribute (we set q2inv[1]=0) to the values of tmpvec
    #vec[1] .= p2    .* (fftplan * (     (bfftplan*(fun .* q2inv)) .* (bfftplan * ones(disc.Nx,disc.Nx,disc.Nx)) ))
    #vec[2] .=           fftplan * ( 2*( (bfftplan*(fun)) .* (bfftplan * ones(disc.Nx,disc.Nx,disc.Nx))  .- ((bfftplan*(fun .* q2inv)) .* (bfftplan*(k2))) ))
    #vec[3] .= p2inv .* (fftplan* (      (bfftplan*(q2 .* fun)) .* (bfftplan * ones(disc.Nx,disc.Nx,disc.Nx))  + ((bfftplan*(fun .* k2inv)) .* (bfftplan*(k4))) - 2*(bfftplan*(fun)) .* (bfftplan*(k2))))
    #angprdterm = 0.25 * real( vec[1] + vec[2] + vec[3])  / disc.vol^2 # 1/V for two bfft products, and one 1/V for loop summation
    #angprdterm[1] = 0. # vanishing external momentum
    #
    #return angprdterm
    #@show (fftplan * (     (bfftplan*(G1T .* q2inv)) .* (bfftplan * ones(disc.Nx,disc.Nx,disc.Nx))  ) )
    
    

    G1T .= simdata.FT[time,time] .- simdata.FL[time,time]
    G1T[1] = 0.

    res .=  p2    .* (fftplan * (     (bfftplan*(G1T .* q2inv)) .* (bfftplan * ones(disc.Nx,disc.Nx,disc.Nx))  ) )
    res .+=           fftplan * ( 2*( (bfftplan*(G1T))          .* (bfftplan * ones(disc.Nx,disc.Nx,disc.Nx))  .- ((bfftplan*(G1T .* q2inv)) .* (bfftplan*(k2))) ) )
    res .+= p2inv .* (fftplan * (     (bfftplan*(q2 .* G1T))    .* (bfftplan * ones(disc.Nx,disc.Nx,disc.Nx))  .+ ((bfftplan*(G1T .* q2inv)) .* (bfftplan*(k4))) .- 2*(bfftplan*(G1T)) .* (bfftplan*(k2))  ) )
    res *= 0.25 / disc.vol^2 # 1/V for two bfft products
    res[1] = 0

    return real(res)
end

export setmass2values!
function setmass2values!(model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, t::Int64)
    @unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, scalarmass2, transvmass2, longitmass2, fftplan, bfftplan = simdata
    #@unpack angprdterm, fun, tmplatvec = tmpdata

    

    #@show  (bfftplan * unitlattice)[1] 
    #@show bfftplan*(k2values)
    #@show disc.vol
    # calculate loops
    #scalarloop = sum( simdata.FS[t,t] ) / disc.vol 
    #gaugeloopT = sum( simdata.FT[t,t] ) / disc.vol 
    #gaugeloopL = sum( simdata.FL[t,t] ) / disc.vol 

    ## Higgs tree-level mass
    scalarmass2[1] = model.Mass^2

    ## Higgs Hartree contribution
    scalarmass2[1] += 8 * model.Lambda * sum( simdata.FS[t,t] ) / disc.vol 
#
    ### Mixed Hartree contributions
    #scalarmass2[1] += model.g^2 * ( (disc.sdim-1)*gaugeloopT + gaugeloopL )
    scalarmass2[1] += model.g^2 * sum( (disc.sdim-1)*simdata.FT[t,t] + simdata.FL[t,t] ) / disc.vol
    transvmass2    .= model.g^2 * sum( simdata.FS[t,t] ) / disc.vol 
    longitmass2    .= model.g^2 * sum( simdata.FS[t,t] ) / disc.vol 
    
    ## Gauge Hartree contributions
    angprdterm = getangprdterm(model, disc, simdata, tmpdata, t)
#
    #transvmass2 .+= model.N * model.g^2 * ( ((disc.sdim-2)*gaugeloopT + gaugeloopL + (gaugeloopT - gaugeloopL)/(disc.sdim-1)) .- (angprdterm/(disc.sdim-1)) )
    #longitmass2 .+= model.N * model.g^2 * ( ((disc.sdim-2)*gaugeloopT + gaugeloopL) .- angprdterm )

    transvmass2 .+= model.N * model.g^2 * sum((disc.sdim-2 + 1/(disc.sdim-1))*simdata.FT[t,t] + (1 - 1/(disc.sdim-1))simdata.FL[t,t] ) /disc.vol 
    transvmass2 .-= (model.N * model.g^2 / (disc.sdim-1))  * angprdterm

    longitmass2 .+= model.N * model.g^2 * sum((disc.sdim-2)*simdata.FT[t,t] + simdata.FL[t,t])/disc.vol 
    longitmass2 .-= model.N * model.g^2 * angprdterm
end


################################################################################################################################################################################
################################################################################################################################################################################
################################################################################################################################################################################

function getE_pot(t::Int64, model::SUNgaugeScalar, meas::MeasurementTwoPIGaugeScalar, pexp::TwoPIGaugeScalarLoopEx, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, disc::TwoPIGaugeScalarDiscretization)
    
    ### Trlog contributions
    E_pot =  0.5 * sum( (model.Mass^2 .+ simdata.k2values) .* simdata.FS[t-1,t-1] )# / disc.vol   # OK
    E_pot += 0.5 * (disc.sdim - 1) * sum( simdata.k2values .* simdata.FT[t-1,t-1] )# / disc.vol   # OK

    ### Hartree contributions
    # calculate loops
    scalarloop = sum( simdata.FS[t-1,t-1] ) / disc.vol 
    gaugeloopT = sum( simdata.FT[t-1,t-1] ) / disc.vol 
    gaugeloopL = sum( simdata.FL[t-1,t-1] ) / disc.vol 

    ## Higgs Hartree contribution 
    #E_pot += 2 * model.Lambda * scalarloop^2 
    E_pot += 2 * model.Lambda * scalarloop *  sum( simdata.FS[t-1,t-1] )

    ## Mixed Hartree contribution
    #E_pot += 0.5 * model.g^2 * ( (disc.sdim-1)*gaugeloopT + gaugeloopL) * scalarloop
    E_pot += 0.5 * model.g^2 * ( (disc.sdim-1)*gaugeloopT + gaugeloopL) * sum( simdata.FS[t-1,t-1] )
    
    ## Gauge Hartree contributions
    angprdterm = getangprdterm(model, disc, simdata, tmpdata, t-1)

    #transvmass2 .+= model.N * model.g^2 * ( ( sum( (disc.sdim-2 + 1/(disc.sdim-1))*simdata.FT[t,t] + (1 - 1/(disc.sdim-1))simdata.FL[t,t] ) /disc.vol ) .- angprdterm/(disc.sdim-1) ) .* simdata.FT[t-1,t-1]
    #transvmass2 .-= (model.N * model.g^2 / (disc.sdim-1))  * angprdterm

    #longitmass2 .+= model.N * model.g^2 * sum((disc.sdim-2)*simdata.FT[t,t] + simdata.FL[t,t])/disc.vol 
    #longitmass2 .-= model.N * model.g^2 * angprdterm

    transmassterm = model.N * model.g^2 * ( ((disc.sdim-2)*gaugeloopT + gaugeloopL + (gaugeloopT - gaugeloopL)/(disc.sdim-1)) .- (angprdterm/(disc.sdim-1)) )
    longimassterm = model.N * model.g^2 * ( ((disc.sdim-2)*gaugeloopT + gaugeloopL) .- angprdterm )

    E_pot += 0.25 * (disc.sdim-1) * sum( transmassterm .* simdata.FT[t-1,t-1] )#/disc.vol
    E_pot += 0.25 * sum( longimassterm .* simdata.FL[t-1,t-1] )#/disc.vol

    return E_pot
end


function getE_mem(t::Int64, model::SUNgaugeScalar, meas::MeasurementTwoPIGaugeScalar, pexp::TwoPIGaugeScalarLoopEx, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, disc::TwoPIGaugeScalarDiscretization)
    
    if t == 2 # first measurement after initialisation, before evolve!
        E_mem = 0
    else # t == 3: first step in trapez integration
        # Mixed sunset contribution
        #E_mem = 0.5 * 0.5 * disc.dt * sum( simdata.SigrS[simdata.indices[1]].*simdata.FS[simdata.indices[1],t-1] - simdata.SigFS[simdata.indices[1]].*simdata.rS[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) )
        # Gauge sunset contributions
        E_mem += 0.5 * disc.dt * (1/3) * sum( Sig_rT .* simdata.FT[simdata.indices[1],t-1] - Sig_FT .* simdata.rT[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) )
        #E_mem += 0.5 * disc.dt * (1/3) * sum( Sig_rL .* simdata.FL[simdata.indices[1],t-1] - Sig_FL .* simdata.rL[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) )
        #E_mem = 0.5 * disc.dt * (1/3) *(disc.sdim-1)* sum( simdata.SigrT[simdata.indices[1]] .* simdata.FT[simdata.indices[1],t-1] - simdata.SigFT[simdata.indices[1]] .* simdata.rT[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) )
        #E_mem += 0.5 * disc.dt * (1/3) * sum( simdata.SigrL[simdata.indices[1]] .* simdata.FL[simdata.indices[1],t-1] - simdata.SigFL[simdata.indices[1]] .* simdata.rL[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) )
        if t > 3 # two or more steps in trapez integration
            # Mixed sunset contribution
            E_mem +=  0.5 * disc.dt * sum( [ sum(simdata.SigrS[tp].*simdata.FS[tp,t-1] - simdata.SigFS[tp].*simdata.rS[tp,t-1]*thesign(tp, t-1, simdata.NstepsinMemory)) for tp in simdata.SigFS.indices[1]+1:(t-2) ] )
            # Gauge sunset contributions
            #E_mem +=  disc.dt * (1/3) *(disc.sdim-1)* sum( [ sum(simdata.SigrT[simdata.indices[1]] .* simdata.FT[simdata.indices[1],t-1] - simdata.SigFT[simdata.indices[1]] .* simdata.rT[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory)) for tp in simdata.SigFT.indices[1]+1:(t-2) ] )
            #E_mem +=  disc.dt * (1/3) * sum( [ sum(simdata.SigrL[simdata.indices[1]] .* simdata.FL[simdata.indices[1],t-1] - simdata.SigFL[simdata.indices[1]] .* simdata.rL[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory)) for tp in simdata.SigFL.indices[1]+1:(t-2) ] )
        end
    end

    return E_mem
end


################################################################################################################################################################################
################################################################################################################################################################################
################################################################################################################################################################################

function calcSelfEnergies2!(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarNLOloop, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    #calcrSMixedSunset!(model, simdata, tmpdata, disc, t, tp)
    #simdata.SigrS[tp] .= tmpdata.res

    ## statistical self-energy contributions  
    #@show tmpdata.tmplatvec[1]
    #simdata.SigFS[tp] = getFSMixedSunset(model,simdata,tmpdata,disc,t,tp)
    #simdata.SigFT[tp] = getFTMixedSunset(model,simdata,tmpdata,disc,t,tp) #+ getFTGaugeSunset(model,simdata,tmpdata,disc,t,tp)
    #simdata.SigFL[tp] = getFLMixedSunset(model,simdata,tmpdata,disc,t,tp) #+ getFLGaugeSunset(model,simdata,tmpdata,disc,t,tp)
#
    ## spectral self-energy contributions  
    #simdata.SigrS[tp] = getrSMixedSunset(model,simdata,tmpdata,disc,t,tp)
    #simdata.SigrT[tp] = getrTMixedSunset(model,simdata,tmpdata,disc,t,tp) #+ getrTGaugeSunset(model,simdata,tmpdata,disc,t,tp)
    #simdata.SigrL[tp] = getrLMixedSunset(model,simdata,tmpdata,disc,t,tp) #+ getrLGaugeSunset(model,simdata,tmpdata,disc,t,tp)

    #tmpdata.tmplatvec[1] =  cT_scalar(simdata.FT[t,tp], simdata.FS[t,tp], simdata, tmpdata, disc)  -  0.25 * cT_scalar(simdata.rT[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata, tmpdata, disc)
    #tmpdata.tmplatvec[1] += cL_scalar(simdata.FL[t,tp], simdata.FS[t,tp], simdata, tmpdata, disc)  -  0.25 * cL_scalar(simdata.rL[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata, tmpdata, disc)
    #tmpdata.tmplatvec[1] *= (-0.25) * model.g^2 # check coefficient 

    #simdata.SigFT[tp] .= 0#(-model.g^2 * ( NT("FF", simdata, tmpdata, disc, t, tp) - 0.25 * NT("rr", simdata, tmpdata, disc, t, tp) ))#getFTGaugeSunset(model,simdata,tmpdata,disc,t,tp)
    #simdata.SigrT[tp] .= 0#getrTGaugeSunset(model,simdata,tmpdata,disc,t,tp)
    #simdata.SigFL[tp] .= 0#getFLGaugeSunset(model,simdata,tmpdata,disc,t,tp)
    #simdata.SigrL[tp] .= 0#getrLGaugeSunset(model,simdata,tmpdata,disc,t,tp)
end

export cT_scalar
function cT_scalar(prop1::lattice, prop2::lattice, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    @unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    p2 = k2values   # pedagogical book-keeping
    q2 = k2values
    k2 = k2values
    p4 = k4values
    k4 = k4values
    q2inv = invk2values

    G1T .= prop1
    G2T .= prop2

    G1T[1] = 0.  # set q=0 contribution to zero (handled in exception)
    G2T[1] = 0.  # set k=0 contribution to zero (handled in exception)

    #tmpvec = [(-1)*fftplan*( (bfftplan*(q2inv .* G1)) .* (bfftplan*(G2)) ),
    #            2*fftplan*( (bfftplan*(q2inv .* G1)) .* (bfftplan*(k2 .* G2)) + (bfftplan*(G1)) .* (bfftplan*(G2)) ),
    #              fftplan*( 2*(bfftplan*(G1)) .* (bfftplan*(k2 .* G2))  -  (bfftplan*(q2inv .* G1)) .* (bfftplan*(k4 .* G2))  -  (bfftplan*(q2 .* G1)) .* (bfftplan*(G2)) )
    #        ]
    #tmp = p4 .* tmpvec[1] + p2 .* tmpvec[2] + tmpvec[3]
    tmpdata.res .*= 0
    tmpdata.res .+= -p4 .* (  fftplan*( (bfftplan*(q2inv .* G1T)) .* (bfftplan*(G2T)) ) )
    tmpdata.res .+=  p2 .* (2*fftplan*( (bfftplan*(q2inv .* G1T)) .* (bfftplan*(k2 .* G2T)) + (bfftplan*(G1T)) .* (bfftplan*(G2T)) ) )
    tmpdata.res .+=           fftplan*( 2*(bfftplan*(G1T)) .* (bfftplan*(k2 .* G2T))  -  (bfftplan*(q2inv .* G1T)) .* (bfftplan*(k4 .* G2T))  -  (bfftplan*(q2 .* G1T)) .* (bfftplan*(G2T)) )

    tmpdata.res .*= disc.ivol^2 # one volume factor for each bfft
    
    # exception handling: 
    tmpdata.res .+= 4*p4 .* G1T[1] .* G2T  # q=0
    tmpdata.res .+= 0                    # k=0
    tmpdata.res[1] = 0                   # p=0

    return real(tmpdata.res)
end

export cL_scalar
function cL_scalar(prop1::lattice, prop2::lattice, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    @unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    p2 = k2values   # pedagogical book-keeping
    q2 = k2values
    k2 = k2values
    p4 = k4values
    k4 = k4values
    q2inv = invk2values

    G1T .= prop1
    G2T .= prop2
    G1T[1] = 0.  # set q=0 contribution to zero (handled in exception)
    G2T[1] = 0.  # set k=0 contribution to zero (handled in exception)
    # tmpvec = [    fftplan*( (bfftplan*(q2inv .* G1)) .* (bfftplan*(G2)) ),
    #          (-2)*fftplan*( (bfftplan*(q2inv .* G1)) .* (bfftplan*(k2 .* G2)) ),
    #               fftplan*( (bfftplan*(q2inv .* G1)) .* (bfftplan*(k4 .* G2)) )
    #         ]
    # tmp = p4 .* tmpvec[1] + p2 .* tmpvec[2] + tmpvec[3]
    tmpdata.res .*= 0
    tmpdata.res .+= p4 .* (      fftplan*( (bfftplan*(q2inv .* G1T)) .* (bfftplan*(G2T)) ))
    tmpdata.res .+= p2 .* ( (-2)*fftplan*( (bfftplan*(q2inv .* G1T)) .* (bfftplan*(k2 .* G2T)) ) )
    tmpdata.res .+=             (fftplan*( (bfftplan*(q2inv .* G1T)) .* (bfftplan*(k4 .* G2T)) ))

    tmpdata.res .*= disc.ivol^2 # one volume factor for each bfft
    
    # exception handling:
    tmpdata.res .+= 0                    # q=0
    tmpdata.res .+= p2 .* G1T .* G2T[1]    # k=0
    tmpdata.res[1] = sum( q2 .* G1T .* (G2T) )/disc.vol   # p=0
    
    return real(tmpdata.res)
end

export CT
function CT(prop1::lattice, prop2::lattice, model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    @unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    p2 = k2values   # pedagogical book-keeping
    q2 = k2values
    k2 = k2values
    q4 = k4values
    k4 = k4values
    p2inv = invk2values

    G1T .= prop1
    G2T .= prop2
    G1T[1] = 0.  # set q=0 contribution to zero (handled in exception)
    G2T[1] = 0.  # set k=0 contribution to zero (handled in exception)

    #tmpvec = [(-1)*fftplan*( (bfftplan*(G1T)) .* (bfftplan*(G2T)) ),
    #             2*fftplan*( (bfftplan*(q2 .* G1T)) .* (bfftplan*(G2T))  .+  (bfftplan*(G1T)) .* (bfftplan*(k2 .* G2T))),
    #               fftplan*((-1)*(bfftplan*(q4 .* G1T)) .* (bfftplan*(G2T))  .+  2*(bfftplan*(q2 .* G1T)) .* (bfftplan*(k2 .* G2T))  .-  (bfftplan*(G1T)) .* (bfftplan*(k4 .* G2T)))
    #        ]
    #tmp = p2 .* tmpvec[1] + tmpvec[2] + p2inv .* tmpvec[3]

    res .*= 0
    res .+= p2 .* ((-1)*fftplan*( (bfftplan*(G1T)) .* (bfftplan*(G2T)) ))
    res .+=          (2*fftplan*( (bfftplan*(q2 .* G1T)) .* (bfftplan*(G2T))  .+  (bfftplan*(G1T)) .* (bfftplan*(k2 .* G2T))))
    res .+= p2inv .*   (fftplan*((-1)*(bfftplan*(q4 .* G1T)) .* (bfftplan*(G2T))  .+  2*(bfftplan*(q2 .* G1T)) .* (bfftplan*(k2 .* G2T))  .-  (bfftplan*(G1T)) .* (bfftplan*(k4 .* G2T))))

    res .*= disc.ivol^2 # one volume factor for each bfft
    
    # exception handling:
    res .+= 0                    # q=0
    res .+= 0                    # k=0
    res[1] = 4 * sum( q2 .* G1T .* G2T )/disc.vol   # p=0
    #@show tmp[1]

    #return real(res)
    res .*= (-0.5) * model.g^2
end

export CL
function CL(prop1::lattice, prop2::lattice, model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    @unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    p2 = k2values   # pedagogical book-keeping
    q2 = k2values
    k2 = k2values
    q4 = k4values
    k4 = k4values
    p2inv = invk2values

    G1T .= prop1
    G2T .= prop2
    G1T[1] = 0.  # set q=0 contribution to zero (handled in exception)
    G2T[1] = 0.  # set k=0 contribution to zero (handled in exception)
    # tmpvec = [
    #             fftplan*((bfftplan*(q4 .* G1)) .* (bfftplan*(G2))  .-  2*(bfftplan*(q2 .* G1)) .* (bfftplan*(k2 .* G2))  .+  (bfftplan*(G1)) .* (bfftplan*(k4 .* G2)))
    #         ]
    # tmp = p2inv .* tmpvec[1]

    res .*= 0
    res .+= p2inv .* ( fftplan * ( (bfftplan*(q4 .* G1T)) .* (bfftplan*(G2T))  .-  2*(bfftplan*(q2 .* G1T)) .* (bfftplan*(k2 .* G2T))  .+  (bfftplan*(G1T)) .* (bfftplan*(k4 .* G2T)) ) )
    
    res .*= disc.ivol^2 # one volume factor for each bfft
    
    # exception handling:
    res .+= p2 .* G1T[1] .* G2T    # q=0
    res .+= p2 .* G1T .* G2T[1]    # k=0
    res[1] = 0                   # p=0

    #return real(res)
    res .*= (-0.5) * model.g^2
end
export getFSMixedSunset
function getFSMixedSunset(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    #tmpdata.tmplatvec[1] *= 0
    tmpdata.tmplatvec[1] =  cT_scalar(simdata.FT[t,tp], simdata.FS[t,tp], simdata, tmpdata, disc)  -  0.25 * cT_scalar(simdata.rT[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata, tmpdata, disc)
    tmpdata.tmplatvec[1] += cL_scalar(simdata.FL[t,tp], simdata.FS[t,tp], simdata, tmpdata, disc)  -  0.25 * cL_scalar(simdata.rL[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata, tmpdata, disc)
    tmpdata.tmplatvec[1] *= (-0.25) * model.g^2 # check coefficient 

    return tmpdata.tmplatvec[1]
end
export getrSMixedSunset
function getrSMixedSunset(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    #tmpdata.tmplatvec[1] *= 0
    tmpdata.tmplatvec[1] =  cT_scalar(simdata.FT[t,tp], simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata, tmpdata, disc)  +  cT_scalar(simdata.rT[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.FS[t,tp], simdata, tmpdata, disc)
    tmpdata.tmplatvec[1] += cL_scalar(simdata.FL[t,tp], simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata, tmpdata, disc)  +  cT_scalar(simdata.rL[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.FS[t,tp], simdata, tmpdata, disc)
    tmpdata.tmplatvec[1] *= (-0.25) * model.g^2 # check coefficient 
    #@show sum(tmpdata.tmplatvec[1])

    return real(tmpdata.tmplatvec[1])
end

## Mixed sunset -> transverse 
#
export getFTMixedSunset
function getFTMixedSunset(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    #tmpdata.tmplatvec[1] *= 0
    tmpdata.tmplatvec[1] =  CT(simdata.FS[t,tp], simdata.FS[t,tp], simdata, tmpdata, disc)  -  0.25 * CT(simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata, tmpdata, disc)
    tmpdata.tmplatvec[1] *= (-0.5) * model.g^2 # check coefficient 

    return tmpdata.tmplatvec[1]
end
export getrTMixedSunset
function getrTMixedSunset(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    #tmpdata.tmplatvec[1] *= 0
    tmpdata.tmplatvec[1] =  CT(simdata.FS[t,tp], simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata, tmpdata, disc)  +  CT(simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.FS[t,tp], simdata, tmpdata, disc)
    tmpdata.tmplatvec[1] *= (-0.5) * model.g^2 # check coefficient 

    return tmpdata.tmplatvec[1]
end

## Mixed sunset -> longitudinal 
#
export getFLMixedSunset
function getFLMixedSunset(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    #tmpdata.tmplatvec[1] *= 0
    tmpdata.tmplatvec[1] =  CL(simdata.FS[t,tp], simdata.FS[t,tp], simdata, tmpdata, disc)  -  0.25 * CL(simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata, tmpdata, disc)
    tmpdata.tmplatvec[1] *= (-0.5) * model.g^2 # check coefficient 

    return tmpdata.tmplatvec[1]
end
export getrLMixedSunset
function getrLMixedSunset(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    #tmpdata.tmplatvec[1] *= 0
    tmpdata.tmplatvec[1] =  CL(simdata.FS[t,tp], simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata, tmpdata, disc)  +  CL(simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.FS[t,tp], simdata, tmpdata, disc)
    tmpdata.tmplatvec[1] *= (-0.5) * model.g^2 # check coefficient 

    return tmpdata.tmplatvec[1]
end
export getFTGaugeSunset
function getFTGaugeSunset(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    #lat =  NT("FF", simdata, tmpdata, disc, t, tp) #- 0.25 * NT("rr", simdata, tmpdata, disc, t, tp)
    #@unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    #G1T = deepcopy(simdata.FT[t,tp])
    #G1L = deepcopy(simdata.FL[t,tp])
    #G2T = deepcopy(simdata.FT[t,tp])
    #G2L = deepcopy(simdata.FL[t,tp])
    #
#
    #G1T[1] = 0.  # set q=0 contribution to zero (handled in exception)
    #G1L[1] = 0.
    #G2T[1] = 0.  # set k=0 contribution to zero (handled in exception)
    #G2L[1] = 0. 
#
#
    #p6 = k6values
    #q2inv = invk2values
    #k2inv = invk2values
#
   #
    #tmpdata.tmplatvec[1] .=  (bfftplan*(q2inv .* G1T)) .* (bfftplan*(k2inv .* G2T))
    #tmpdata.tmplatvec[1] .-= (bfftplan*(q2inv .* G1T)) .* (bfftplan*(k2inv .* G2L)) 
    #tmpdata.tmplatvec[1] .-= (bfftplan*(q2inv .* G1L)) .* (bfftplan*(k2inv .* G2T)) 
    #tmpdata.tmplatvec[1] .+= (bfftplan*(q2inv .* G1L)) .* (bfftplan*(k2inv .* G2L))
#
    #fft!(tmpdata.tmplatvec[1])
    #tmpdata.tmplatvec[1] .*= p6 
    #tmpdata.tmplatvec[1] .*= (-0.5)/disc.vol^2
    

    #return real(tmpdata.tmplatvec[1])#model.g^2 * ( NT("FF", simdata, tmpdata, disc, t, tp) - 0.25 * NT("rr", simdata, tmpdata, disc, t, tp) )
    return -model.g^2 * ( NT("FF", simdata, tmpdata, disc, t, tp) - 0.25 * NT("rr", simdata, tmpdata, disc, t, tp) )
end
export getrTGaugeSunset
function getrTGaugeSunset(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    lat = NT("Fr", simdata, tmpdata, disc, t, tp) + NT("rF", simdata, tmpdata, disc, t, tp)
    lat *= -model.g^2
    return lat#(-model.g^2) * ( NT("Fr", simdata, tmpdata, disc, t, tp) + NT("rF", simdata, tmpdata, disc, t, tp) )
end

## Gauge sunset -> Longitudinal
# 
export getFLGaugeSunset
function getFLGaugeSunset(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    lat = NL("FF", simdata, tmpdata, disc, t, tp) - 0.25 * NL("rr", simdata, tmpdata, disc, t, tp)
    lat *= -model.g^2
    return lat#(-model.g^2) * ( NL("FF", simdata, tmpdata, disc, t, tp) - 0.25 * NL("rr", simdata, tmpdata, disc, t, tp) )
end
export getrLGaugeSunset
function getrLGaugeSunset(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    lat = NL("Fr", simdata, tmpdata, disc, t, tp) + NL("rF", simdata, tmpdata, disc, t, tp)
    lat *= -model.g^2
    return lat#(-model.g^2) * ( NL("Fr", simdata, tmpdata, disc, t, tp) + NL("rF", simdata, tmpdata, disc, t, tp) )
end
## Gauge sunset contributions
#
export NT
function NT(string::String, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    #@unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    #@unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    p2 = simdata.k2values
    q2 = simdata.k2values
    
    if string == "FF"
        tmpdata.G1T .= simdata.FT[t,tp]
        tmpdata.G1L .= simdata.FL[t,tp]
        tmpdata.G2T .= simdata.FT[t,tp]
        tmpdata.G2L .= simdata.FL[t,tp]
    elseif string == "rF"
        tmpdata.G1T .= simdata.rT[t,tp]*thesign(t, tp, simdata.NstepsinMemory)
        tmpdata.G1L .= simdata.rL[t,tp]*thesign(t, tp, simdata.NstepsinMemory)
        tmpdata.G2T .= simdata.FT[t,tp]
        tmpdata.G2L .= simdata.FL[t,tp]
    elseif string == "Fr"
        tmpdata.G1T .= simdata.FT[t,tp]
        tmpdata.G1L .= simdata.FL[t,tp]
        tmpdata.G2T .= simdata.rT[t,tp]*thesign(t, tp, simdata.NstepsinMemory)
        tmpdata.G2L .= simdata.rL[t,tp]*thesign(t, tp, simdata.NstepsinMemory)
    elseif string == "rr"
        tmpdata.G1T .= simdata.rT[t,tp]*thesign(t, tp, simdata.NstepsinMemory)
        tmpdata.G1L .= simdata.rL[t,tp]*thesign(t, tp, simdata.NstepsinMemory)
        tmpdata.G2T .= simdata.rT[t,tp]*thesign(t, tp, simdata.NstepsinMemory)
        tmpdata.G2L .= simdata.rL[t,tp]*thesign(t, tp, simdata.NstepsinMemory)
    end
    tmpdata.G1T[1] = 0.  # set q=0 contribution to zero (handled in exception)
    tmpdata.G1L[1] = 0.
    tmpdata.G2T[1] = 0.  # set k=0 contribution to zero (handled in exception)
    tmpdata.G2L[1] = 0. 

    tmpdata.res .*= 0

    calcNTp6term!(simdata, tmpdata, disc)
    calcNTp4term!(simdata, tmpdata, disc)
    calcNTp2term!(simdata, tmpdata, disc)
    calcNTp0term!(simdata, tmpdata, disc)
    calcNTp2invterm!(simdata, tmpdata, disc)
    tmpdata.res .*= disc.ivol^2 # one volume factor for each bfft

    # exception handling:
    tmpdata.res .+=               (-2) * p2 .* ((disc.sdim-2)*tmpdata.G1T[1] - tmpdata.G1L[1]) .* tmpdata.G2L .+ (disc.sdim-1)*tmpdata.G1T[1] .* tmpdata.G2T                          # q=0
    tmpdata.res .+= (-2)*(disc.sdim-1) * p2 .* (tmpdata.G1T .+ tmpdata.G1L) .* tmpdata.G2T[1]                                                           # k=0
    tmpdata.res[1] = sum(                q2 .* ((disc.sdim-2)*tmpdata.G1L .* tmpdata.G2T .+ tmpdata.G1T .* ((disc.sdim-2)*tmpdata.G2L .+ 4*(disc.sdim-1)*tmpdata.G2T) ))/disc.vol    # p=0
    
    return real(tmpdata.res)
end

export calcNTp6term!
function calcNTp6term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    @unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    p6 = k6values
    q2inv = invk2values
    k2inv = invk2values

    tmpdata.tmp .*= 0 
    tmpdata.tmp .+= (bfftplan*(q2inv .* G1T)) .* (bfftplan*(k2inv .* G2T))
    tmpdata.tmp .-= (bfftplan*(q2inv .* G1T)) .* (bfftplan*(k2inv .* G2L))
    tmpdata.tmp .-= (bfftplan*(q2inv .* G1L)) .* (bfftplan*(k2inv .* G2T))
    tmpdata.tmp .+= (bfftplan*(q2inv .* G1L)) .* (bfftplan*(k2inv .* G2L))

    tmpdata.res .+= 0.5 * p6 .* ( fftplan * tmpdata.tmp )

    #@show sum(tmpdata.tmp)
    #tmpdata.res .+= 0.5 * p6 .* (fftplan * (  (bfftplan*(q2inv .* G1T)) .* (bfftplan*(k2inv .* G2T)) - (bfftplan*(q2inv .* G1T)) .* (bfftplan*(k2inv .* G2L))
    #                                        - (bfftplan*(q2inv .* G1L)) .* (bfftplan*(k2inv .* G2T)) + (bfftplan*(q2inv .* G1L)) .* (bfftplan*(k2inv .* G2L)) ))

    
end

export calcNTp4term!
function calcNTp4term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    @unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    p4 = k4values
    q2inv = invk2values
    k2inv = invk2values

    tmpdata.tmp .*= 0
    tmpdata.tmp .+=  0.5*(disc.sdim-3.5)*( (bfftplan*(q2inv .* G1T)).*(bfftplan*(G2T)) + (bfftplan*(G1T)).*(bfftplan*(k2inv .* G2T)) )
    tmpdata.tmp .+= (bfftplan*(q2inv .* G1T)) .* (bfftplan*(G2L))  +  0.5*(3.5-disc.sdim)*(bfftplan*(G1T)) .* (bfftplan*(k2inv .* G2L))
    tmpdata.tmp .+= 0.5*(3.5-disc.sdim)*(bfftplan*(q2inv .* G1L)) .* (bfftplan*(G2T))  +  (bfftplan*(G1L)) .* (bfftplan*(k2inv .* G2T))
    tmpdata.tmp .-= ( (bfftplan*(q2inv .* G1L)) .* (bfftplan*(G2L)) + (bfftplan*(G1L)) .* (bfftplan*(k2inv .* G2L)) )
    
    tmpdata.res .+=  p4 .* (fftplan * tmpdata.tmp)

    #@show sum(tmpdata.tmp)

    #tmpdata.res .+=  p4 .* (fftplan * ( 0.5*(disc.sdim-3.5)*( (bfftplan*(q2inv .* G1T)) .* (bfftplan*(G2T))  +                      (bfftplan*(G1T)) .* (bfftplan*(k2inv .* G2T)) )
    #                                    +                     (bfftplan*(q2inv .* G1T)) .* (bfftplan*(G2L))  +  0.5*(3.5-disc.sdim)*(bfftplan*(G1T)) .* (bfftplan*(k2inv .* G2L))
    #                                    + 0.5*(3.5-disc.sdim)*(bfftplan*(q2inv .* G1L)) .* (bfftplan*(G2T))  +                      (bfftplan*(G1L)) .* (bfftplan*(k2inv .* G2T))
    #                                    -                     (bfftplan*(q2inv .* G1L)) .* (bfftplan*(G2L))  -                      (bfftplan*(G1L)) .* (bfftplan*(k2inv .* G2L)) 
    #                                ))
#

    
    #tmpdata.res .+= (bfftplan*(q2inv .* G1T)) .* (bfftplan*(G2L))  +  0.5*(3.5-disc.sdim)*(bfftplan*(G1T)) .* (bfftplan*(k2inv .* G2L))
    #tmpdata.res .+= 0.5*(3.5-disc.sdim)*(bfftplan*(q2inv .* G1L)) .* (bfftplan*(G2T))  +  (bfftplan*(G1L)) .* (bfftplan*(k2inv .* G2T))
    #tmpdata.res .-= ( (bfftplan*(q2inv .* G1L)) .* (bfftplan*(G2L)) + (bfftplan*(G1L)) .* (bfftplan*(k2inv .* G2L)) )
    #
    #return (p4 .* (fftplan * tmpdata.tmplatvec[2]))

end

export calcNTp2term!
function calcNTp2term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    @unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    p2 = k2values
    q2 = k2values
    k2 = k2values
    q2inv = invk2values
    k2inv = invk2values
    
    # tmpdata.tmplatvec[2] *= 0
    # tmpdata.tmplatvec[2] .+=  -(bfftplan*(G1L)).*(bfftplan*(G2L)) 
    # tmpdata.tmplatvec[2] .+= (0.5*( (bfftplan*(q2inv .* G1L)).*(bfftplan*(k2 .* G2L)) + (bfftplan*(q2 .* G1L)).*(bfftplan*(k2inv .* G2L)) ) )
    # tmpdata.tmplatvec[2] .+= ((2-1.5*disc.sdim)*(bfftplan*(G1L)).*(bfftplan*(G2T)) + (disc.sdim-1.75)*(bfftplan*(q2inv .* G1L)).*(bfftplan*(k2 .* G2T)) - 0.5*(bfftplan*(q2 .* G1L)).*(bfftplan*(k2inv .* G2T)) )
    # tmpdata.tmplatvec[2] .+= ((2-1.5*disc.sdim)*(bfftplan*(G1T)).*(bfftplan*(G2L)) - 0.5*(bfftplan*(q2inv .* G1T)).*(bfftplan*(k2 .* G2L)) + (disc.sdim-1.75)*(bfftplan*(q2 .* G1T)).*(bfftplan*(k2inv .* G2L)) )
    # tmpdata.tmplatvec[2] .+= ((3.5-3*disc.sdim)*(bfftplan*(G1T)).*(bfftplan*(G2T)) + (1.75-disc.sdim)*( (bfftplan*(q2inv .* G1T)).*(bfftplan*(k2 .* G2T)) + (bfftplan*(q2 .* G1T)).*(bfftplan*(k2inv .* G2T)) ))

    # return (p2 .* (fftplan * tmpdata.tmplatvec[2]))

    tmpdata.tmp .*= 0
    tmpdata.tmp .+=  -(bfftplan*(G1L)).*(bfftplan*(G2L)) 
    tmpdata.tmp .+= (0.5*( (bfftplan*(q2inv .* G1L)).*(bfftplan*(k2 .* G2L)) + (bfftplan*(q2 .* G1L)).*(bfftplan*(k2inv .* G2L)) ) )
    tmpdata.tmp .+= ((2-1.5*disc.sdim)*(bfftplan*(G1L)).*(bfftplan*(G2T)) + (disc.sdim-1.75)*(bfftplan*(q2inv .* G1L)).*(bfftplan*(k2 .* G2T)) - 0.5*(bfftplan*(q2 .* G1L)).*(bfftplan*(k2inv .* G2T)) )
    tmpdata.tmp .+= ((2-1.5*disc.sdim)*(bfftplan*(G1T)).*(bfftplan*(G2L)) - 0.5*(bfftplan*(q2inv .* G1T)).*(bfftplan*(k2 .* G2L)) + (disc.sdim-1.75)*(bfftplan*(q2 .* G1T)).*(bfftplan*(k2inv .* G2L)) )
    tmpdata.tmp .+= ((3.5-3*disc.sdim)*(bfftplan*(G1T)).*(bfftplan*(G2T)) + (1.75-disc.sdim)*( (bfftplan*(q2inv .* G1T)).*(bfftplan*(k2 .* G2T)) + (bfftplan*(q2 .* G1T)).*(bfftplan*(k2inv .* G2T)) ))
    
    tmpdata.res .+= p2 .* (fftplan  * tmpdata.tmp)

end

export calcNTp0term!
function calcNTp0term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    @unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    q2 = k2values
    k2 = k2values
    q4 = k4values
    k4 = k4values
    q2inv = invk2values
    k2inv = invk2values
    
    tmpdata.tmp .*= 0
    tmpdata.tmp .+= ( 0.5*(3*disc.sdim-5) * (bfftplan*(G1L))  .* (bfftplan*(k2 .* G2T)) + 0.5*(disc.sdim-0.5) * (bfftplan*(q2inv .* G1L))  .* (bfftplan*(k4 .* G2T)) +                0.25 * (bfftplan*(q2 .* G1L)) .* (bfftplan*(G2T)))
    tmpdata.tmp .+= (                0.25 * (bfftplan*(G1T))  .* (bfftplan*(k2 .* G2L)) + 0.5*(3*disc.sdim-5) * (bfftplan*(G1T))           .* (bfftplan*(q2 .* G2L)) + 0.5*(0.5-disc.sdim) * (bfftplan*(q4 .* G1T)) .* (bfftplan*(k2inv .* G2L)) )
    tmpdata.tmp .+= (0.5*(3*disc.sdim-3.5)* ((bfftplan*(G1T)) .* (bfftplan*(k2 .* G2T)) +                       (bfftplan*(q2inv .* G1T))  .* (bfftplan*(k4 .* G2T)) +                       (bfftplan*(q2 .* G1T)) .* (bfftplan*(G2T))              + (bfftplan*(q4 .* G1T)) .* (bfftplan*(k2inv .* G2T)) ))

    tmpdata.res .+= (fftplan * tmpdata.tmp)
    #return (fftplan * tmpdata.tmplatvec[2])
end

export calcNTp2invterm!
function calcNTp2invterm!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    @unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    q2 = k2values
    k2 = k2values
    q4 = k4values
    k4 = k4values
    q6 = k6values
    k6 = k6values
    q2inv = invk2values
    k2inv = invk2values
    p2inv = invk2values

    tmpdata.tmp .*= 0
    tmpdata.tmp .+= (-0.5)*(bfftplan*(G1L)).*(bfftplan*(k4 .* G2T)) + 0.25*(bfftplan*(q2inv .* G1L)).*(bfftplan*(k6 .* G2T)) + 0.25*(bfftplan*(q2 .* G1L)).*(bfftplan*(k2 .* G2T))
    tmpdata.tmp .+= (-0.5)*(bfftplan*(q4 .* G1T)).*(bfftplan*(G2L)) + 0.25*(bfftplan*(q6 .* G1T)).*(bfftplan*(k2inv .* G2L)) + 0.25*(bfftplan*(q2 .* G1T)).*(bfftplan*(k2 .* G2L))
    tmpdata.tmp .+= (2-disc.sdim)*(bfftplan*(G1T)).*(bfftplan*(k4 .* G2T)) - 0.25*(bfftplan*(q2inv .* G1T)).*(bfftplan*(k6 .* G2T)) + (2*disc.sdim-3.5)*(bfftplan*(q2 .* G1T)).*(bfftplan*(k2 .* G2T)) + (2-disc.sdim)*(bfftplan*(q4 .* G1T)).*(bfftplan*(G2T)) - 0.25*(bfftplan*(q6 .* G1T)).*(bfftplan*(k2inv .* G2T))
    
    tmpdata.res .+= p2inv .* (fftplan * tmpdata.tmp) 
    #return ( p2inv .* (fftplan * tmpdata.tmplatvec[2]) )
end
##########################################################################################################################################################################

export NL
function NL(string::String, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    # pedagogical book-keeping
    p2 = simdata.k2values
    q2 = simdata.k2values
    
    if string == "FF"
        tmpdata.G1T .= simdata.FT[t,tp]
        tmpdata.G1L .= simdata.FL[t,tp]
        tmpdata.G2T .= simdata.FT[t,tp]
        tmpdata.G2L .= simdata.FL[t,tp]
    elseif string == "rF"
        tmpdata.G1T .= simdata.rT[t,tp]*thesign(t, tp, simdata.NstepsinMemory)
        tmpdata.G1L .= simdata.rL[t,tp]*thesign(t, tp, simdata.NstepsinMemory)
        tmpdata.G2T .= simdata.FT[t,tp]
        tmpdata.G2L .= simdata.FL[t,tp]
    elseif string == "Fr"
        tmpdata.G1T .= simdata.FT[t,tp]
        tmpdata.G1L .= simdata.FL[t,tp]
        tmpdata.G2T .= simdata.rT[t,tp]*thesign(t, tp, simdata.NstepsinMemory)
        tmpdata.G2L .= simdata.rL[t,tp]*thesign(t, tp, simdata.NstepsinMemory)
    elseif string == "rr"
        tmpdata.G1T .= simdata.rT[t,tp]*thesign(t, tp, simdata.NstepsinMemory)
        tmpdata.G1L .= simdata.rL[t,tp]*thesign(t, tp, simdata.NstepsinMemory)
        tmpdata.G2T .= simdata.rT[t,tp]*thesign(t, tp, simdata.NstepsinMemory)
        tmpdata.G2L .= simdata.rL[t,tp]*thesign(t, tp, simdata.NstepsinMemory)
    end
    tmpdata.G1T[1] = 0.  # set q=0 contribution to zero (handled in exception)
    tmpdata.G1L[1] = 0.
    tmpdata.G2T[1] = 0.  # set k=0 contribution to zero (handled in exception)
    tmpdata.G2L[1] = 0. 

    tmpdata.res .*= 0

    calcNLp2term!(simdata, tmpdata, disc)
    calcNLp0term!(simdata, tmpdata, disc)
    calcNLp2invterm!(simdata, tmpdata, disc)

    tmpdata.res .*= disc.ivol^2 # one volume factor for each bfft
    
    # exception handling:
    tmpdata.res .+=                 p2 .* ((disc.sdim-2)*tmpdata.G1T[1] + tmpdata.G1L[1]) .* tmpdata.G2T  # q=0
    tmpdata.res .+= (disc.sdim-1) * p2 .* tmpdata.G1T .* tmpdata.G2T[1]                                   # k=0
    tmpdata.res[1] = sum( q2 .* ((tmpdata.G1L .* tmpdata.G2T)  .+  (tmpdata.G1T .* tmpdata.G2L)) )/disc.vol    # p=0
    
    return real(tmpdata.res)
end

export calcNLp2term!
function calcNLp2term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    @unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    p2 = k2values
    q2 = k2values
    k2 = k2values
    q2inv = invk2values
    k2inv = invk2values

    tmpdata.tmp .*= 0
    tmpdata.tmp .+= (-0.25)*(bfftplan*(q2inv .* G1L)).*(bfftplan*(k2 .* G2T)) 
    tmpdata.tmp .+= (-0.25)*(bfftplan*(q2 .* G1T)).*(bfftplan*(k2inv .* G2L)) 
    tmpdata.tmp .+= (-0.5 )*(bfftplan*(G1T)).*(bfftplan*(G2T)) + 0.25*(bfftplan*(q2inv .* G1T)).*(bfftplan*(k2 .* G2T)) + 0.25*(bfftplan*(q2 .* G1T)).*(bfftplan*(k2inv .* G2T))

    tmpdata.res .+= p2 .* (fftplan * tmpdata.tmp)
    #return (p2 .* (fftplan * tmpdata.tmplatvec[3]))
end
export calcNLp0term!
function calcNLp0term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    @unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    q2 = k2values
    k2 = k2values
    q4 = k4values
    k4 = k4values
    q2inv = invk2values
    k2inv = invk2values

    tmpdata.tmp .*= 0
    tmpdata.tmp .+= 0.5 * ( (bfftplan*(G1L)).*(bfftplan*(k2 .* G2T)) .+ (bfftplan*(q2inv .* G1L)).*(bfftplan*(k4 .* G2T)) )
    tmpdata.tmp .+= 0.5 * ( (bfftplan*(q2 .* G1T)).*(bfftplan*(G2L)) .+ (bfftplan*(q4 .* G1T)).*(bfftplan*(k2inv .* G2L)) )
    tmpdata.tmp .+= 0.5 * ( (bfftplan*(G1T)).*(bfftplan*(k2 .* G2T)) .- (bfftplan*(q2inv .* G1T)).*(bfftplan*(k4 .* G2T)) .+ (bfftplan*(q2 .* G1T)).*(bfftplan*(G2T)) .- (bfftplan*(q4 .* G1T)).*(bfftplan*(k2inv .* G2T)) )
    
    tmpdata.res .+= (fftplan * tmpdata.tmp)
    #return (fftplan * tmpdata.tmplatvec[3])
end
export calcNLp2invterm!
function calcNLp2invterm!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    @unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    q2 = k2values
    k2 = k2values
    q4 = k4values
    k4 = k4values
    q6 = k6values
    k6 = k6values
    q2inv = invk2values
    k2inv = invk2values
    p2inv = invk2values

    tmpdata.tmp .*= 0
    tmpdata.tmp .+= 0.5*(bfftplan*(G1L)).*(bfftplan*(k4 .* G2T)) .- 0.25*(bfftplan*(q2inv .* G1L)).*(bfftplan*(k6 .* G2T)) .- 0.25*(bfftplan*(q2 .* G1L)).*(bfftplan*(k2 .* G2T))
    tmpdata.tmp .+= (-0.25)*(bfftplan*(q2 .* G1T)).*(bfftplan*(k2 .* G2L)) .+ 0.5*(bfftplan*(q4 .* G1T)).*(bfftplan*(G2L)) .- 0.25*(bfftplan*(q6 .* G1T)).*(bfftplan*(k2inv .* G2L))
    tmpdata.tmp .+= (disc.sdim-2)*( (bfftplan*(G1T)).*(bfftplan*(k4 .* G2T)) .+ (bfftplan*(q4 .* G1T)).*(bfftplan*(G2T)) ) .+ 0.25*( (bfftplan*(q2inv .* G1T)).*(bfftplan*(k6 .* G2T)) .+ (bfftplan*(q6 .* G1T)).*(bfftplan*(k2inv .* G2T))) .+ (2*disc.sdim-3.5)*(bfftplan*(q2 .* G1T)).*(bfftplan*(k2 .* G2T))
    
    tmpdata.res .+= p2inv .* (fftplan * tmpdata.tmp)
    #return (p2inv .* (fftplan * tmpdata.tmplatvec[3]))
end