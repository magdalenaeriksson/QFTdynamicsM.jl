# Mass terms
export calcAngprdterm!
function calcAngprdterm!(model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, t::Int64)
    #@unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack k2, k4, k6, k2inv, fftplan, bfftplan = simdata
    @unpack res, G1T = tmpdata
    
    # calculate angular product terms (introduce p,q,k for pedagogical labelling/motivation)
    #p2 = k2values # represents external momentum
    #q2 = k2values # loop lattice momentum
    #k2 = k2values # k==p-q --> separately handled loop momentum
    #k4 = k4values
    #p2inv = invk2values
    #q2inv = invk2values
    #k2inv = invk2values

    G1T .= simdata.FT[t,t] .- simdata.FL[t,t]
    G1T[1] = 0.
    # res .=  p2    .* (fftplan * (     (bfftplan*(G1T .* q2inv)) .* (bfftplan * ones(disc.Nx,disc.Nx,disc.Nx))  ) )
    # res .+=           fftplan * ( 2*( (bfftplan*(G1T))          .* (bfftplan * ones(disc.Nx,disc.Nx,disc.Nx))  .- ((bfftplan*(G1T .* q2inv)) .* (bfftplan*(k2))) ) )
    # res .+= p2inv .* (fftplan * (     (bfftplan*(q2 .* G1T))    .* (bfftplan * ones(disc.Nx,disc.Nx,disc.Nx))  .+ ((bfftplan*(G1T .* q2inv)) .* (bfftplan*(k4))) .- 2*(bfftplan*(G1T)) .* (bfftplan*(k2))  ) )
    # res .*= 0.25 / disc.vol^2 # 1/V for two bfft products


    # below should be equal to above since eg: fftplan * ( (bfftplan*(simdata.FT[t,t])) .* (bfftplan * ones(disc.Nx,disc.Nx,disc.Nx))) == (bfftplan*(simdata.FT[t,t]))[1]*disc.vol

    #res .=  p2    .* ( (bfftplan*(G1T .* q2inv))[1] * disc.vol )
    #res .+= ( 2*(bfftplan*(G1T))[1]*disc.vol )  .-  2*(fftplan * ((bfftplan*(G1T .* q2inv)) .* (bfftplan*(k2)) ))
    #res .+= p2inv .* ( ((bfftplan*(q2 .* G1T))[1]*disc.vol)  .+ ( fftplan * (((bfftplan*(G1T .* q2inv)) .* (bfftplan*(k4))) .- 2*(bfftplan*(G1T)) .* (bfftplan*(k2))  ) ))
    #res .*= 0.25 / disc.vol^2 # 1/V for two bfft products

    res .=  k2    .* ( (bfftplan*(G1T .* k2inv))[1] * disc.vol )
    res .+= ( 2*(bfftplan*(G1T))[1]*disc.vol )  .-  2*(fftplan * ((bfftplan*(G1T .* k2inv)) .* (bfftplan*(k2)) ))
    res .+= k2inv .* ( ((bfftplan*(k2 .* G1T))[1]*disc.vol)  .+ ( fftplan * (((bfftplan*(G1T .* k2inv)) .* (bfftplan*(k4))) .- 2*(bfftplan*(G1T)) .* (bfftplan*(k2))  ) ))
    res .*= 0.25 / disc.vol^2 # 1/V for two bfft products


    # exception handling
    #res .+= 0  # q=0 (zero loop momentum)
    #res .+= 0  # k=0 (loop == external momentum)
    res[1] = 0  # p=0 (zero external momentum)
    
end
export setmass2values!
function setmass2values!(model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, t::Int64)
    @unpack FS, FT, FL, rS, rT, rL, k2, k4, k6, k2inv, scalarmass2, transvmass2, longitmass2, fftplan, bfftplan = simdata
    #@unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, scalarmass2, transvmass2, longitmass2, fftplan, bfftplan = simdata
    #@unpack angprdterm, fun, tmplatvec = tmpdata

    ## Higgs tree-level mass
    scalarmass2[1] = model.Mass^2

    ## Higgs Hartree contribution
    scalarmass2[1] += 8 * model.Lambda * sum( simdata.FS[t,t] ) / disc.vol 
#
    ### Mixed Hartree contributions
    scalarmass2[1] += model.g^2 * sum( (disc.sdim-1)*simdata.FT[t,t] + simdata.FL[t,t] ) / disc.vol
    transvmass2    .= model.g^2 * sum( simdata.FS[t,t] ) / disc.vol 
    longitmass2    .= model.g^2 * sum( simdata.FS[t,t] ) / disc.vol 
    
    ## Gauge Hartree contributions
    calcAngprdterm!(model, disc, simdata, tmpdata, t) # sets angprdterm == tmpdata.res

    transvmass2 .+= model.N * model.g^2 * sum((disc.sdim-2 + 1/(disc.sdim-1))*simdata.FT[t,t] + (1 - 1/(disc.sdim-1))simdata.FL[t,t] ) /disc.vol 
    transvmass2 .-= (model.N * model.g^2 / (disc.sdim-1)) .* tmpdata.res

    longitmass2 .+= model.N * model.g^2 * sum((disc.sdim-2)*simdata.FT[t,t] + simdata.FL[t,t])/disc.vol 
    longitmass2 .-= model.N * model.g^2 .* tmpdata.res
end
########################################################################################################################################################################
export calc_cTHiggs!
function calc_cTHiggs!(prop1::lattice, prop2::lattice, model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    #@unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack FS, FT, FL, rS, rT, rL, k2, k4, k6, k2inv, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    #p2 = k2values   # pedagogical book-keeping
    #q2 = k2values
    #k2 = k2values
    #p4 = k4values
    #k4 = k4values
    #q2inv = invk2values

    G1T .= prop1
    G2T .= prop2

    G1T[1] = 0.  # set q=0 contribution to zero (handled in exception)
    G2T[1] = 0.  # set k=0 contribution to zero (handled in exception)

    #res .*= 0
    #res .+= -p4 .* (  fftplan*( (bfftplan*(q2inv .* G1T)) .* (bfftplan*(G2T)) ) )
    #res .+=  p2 .* (2*fftplan*( (bfftplan*(q2inv .* G1T)) .* (bfftplan*(k2 .* G2T)) + (bfftplan*(G1T)) .* (bfftplan*(G2T)) ) )
    #res .+=           fftplan*( 2*(bfftplan*(G1T)) .* (bfftplan*(k2 .* G2T))  -  (bfftplan*(q2inv .* G1T)) .* (bfftplan*(k4 .* G2T))  -  (bfftplan*(q2 .* G1T)) .* (bfftplan*(G2T)) )
    res .*= 0
    res .+= -k4 .* (  fftplan*( (bfftplan*(k2inv .* G1T)) .* (bfftplan*(G2T)) ) )
    res .+=  k2 .* (2*fftplan*( (bfftplan*(k2inv .* G1T)) .* (bfftplan*(k2 .* G2T)) + (bfftplan*(G1T)) .* (bfftplan*(G2T)) ) )
    res .+=           fftplan*( 2*(bfftplan*(G1T)) .* (bfftplan*(k2 .* G2T))  -  (bfftplan*(k2inv .* G1T)) .* (bfftplan*(k4 .* G2T))  -  (bfftplan*(k2 .* G1T)) .* (bfftplan*(G2T)) )

    res .*= disc.ivol^2 # one volume factor for each bfft
    
    # exception handling: 
    #res .+= 4*p4 .* G1T[1] .* G2T  # q=0
    #res .+= 0                      # k=0
    #res[1] = 0                     # p=0
    res .+= 4*k4 .* G1T[1] .* G2T  # q=0
    res .+= 0                      # k=0
    res[1] = 0                     # p=0

    res .*= (-0.25) * model.g^2

end
export calc_cLHiggs!
function calc_cLHiggs!(prop1::lattice, prop2::lattice, model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    #@unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack FS, FT, FL, rS, rT, rL, k2, k4, k6, k2inv, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    #p2 = k2values   # pedagogical book-keeping
    #q2 = k2values
    #k2 = k2values
    #p4 = k4values
    #k4 = k4values
    #q2inv = invk2values

    G1T .= prop1
    G2T .= prop2
    G1T[1] = 0.  # set q=0 contribution to zero (handled in exception)
    G2T[1] = 0.  # set k=0 contribution to zero (handled in exception)

    #res .*= 0
    #res .+= p4 .* (      fftplan*( (bfftplan*(q2inv .* G1T)) .* (bfftplan*(G2T)) ))
    #res .+= p2 .* ( (-2)*fftplan*( (bfftplan*(q2inv .* G1T)) .* (bfftplan*(k2 .* G2T)) ) )
    #res .+=             (fftplan*( (bfftplan*(q2inv .* G1T)) .* (bfftplan*(k4 .* G2T)) ))
    res .*= 0
    res .+= k4 .* (      fftplan*( (bfftplan*(k2inv .* G1T)) .* (bfftplan*(G2T)) ))
    res .+= k2 .* ( (-2)*fftplan*( (bfftplan*(k2inv .* G1T)) .* (bfftplan*(k2 .* G2T)) ) )
    res .+=             (fftplan*( (bfftplan*(k2inv .* G1T)) .* (bfftplan*(k4 .* G2T)) ))

    res .*= disc.ivol^2 # one volume factor for each bfft
    
    # exception handling:
    #res .+= 0                                     # q=0
    #res .+= p2 .* G1T .* G2T[1]                   # k=0
    #res[1] = sum( q2 .* G1T .* (G2T) )/disc.vol   # p=0
    res .+= 0                                     # q=0
    res .+= k2 .* G1T .* G2T[1]                   # k=0
    res[1] = sum( k2 .* G1T .* (G2T) )/disc.vol   # p=0

    res .*= (-0.25) * model.g^2
end
export calc_CTGauge!
function calc_CTGauge!(prop1::lattice, prop2::lattice, model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    #@unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack FS, FT, FL, rS, rT, rL, k2, k4, k6, k2inv, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    #p2 = k2values   # pedagogical book-keeping
    #q2 = k2values
    #k2 = k2values
    #q4 = k4values
    #k4 = k4values
    #p2inv = invk2values

    G1T .= prop1
    G2T .= prop2
    G1T[1] = 0.  # set q=0 contribution to zero (handled in exception)
    G2T[1] = 0.  # set k=0 contribution to zero (handled in exception)

    #res .*= 0
    #res .+= p2 .* ((-1)*fftplan*( (bfftplan*(G1T)) .* (bfftplan*(G2T)) ))
    #res .+= 2*fftplan*( (bfftplan*(q2 .* G1T)) .* (bfftplan*(G2T))  .+  (bfftplan*(G1T)) .* (bfftplan*(k2 .* G2T)))
    #res .+= p2inv .* (fftplan*((-1)*(bfftplan*(q4 .* G1T)) .* (bfftplan*(G2T))  .+  2*(bfftplan*(q2 .* G1T)) .* (bfftplan*(k2 .* G2T))  .-  (bfftplan*(G1T)) .* (bfftplan*(k4 .* G2T))))
    res .*= 0
    res .+= k2 .* ((-1)*fftplan*( (bfftplan*(G1T)) .* (bfftplan*(G2T)) ))
    res .+= 2*fftplan*( (bfftplan*(k2 .* G1T)) .* (bfftplan*(G2T))  .+  (bfftplan*(G1T)) .* (bfftplan*(k2 .* G2T)))
    res .+= k2inv .* (fftplan*((-1)*(bfftplan*(k4 .* G1T)) .* (bfftplan*(G2T))  .+  2*(bfftplan*(k2 .* G1T)) .* (bfftplan*(k2 .* G2T))  .-  (bfftplan*(G1T)) .* (bfftplan*(k4 .* G2T))))

    res .*= disc.ivol^2 # one volume factor for each bfft
    
    # exception handling:
    #res .+= 0                                       # q=0
    #res .+= 0                                       # k=0
    #res[1] = 4 * sum( q2 .* G1T .* G2T )/disc.vol   # p=0
    res .+= 0                                       # q=0
    res .+= 0                                       # k=0
    res[1] = 4 * sum( k2 .* G1T .* G2T )/disc.vol   # p=0

    res .*= (-0.5) * model.g^2
end
export calc_CLGauge!
function calc_CLGauge!(prop1::lattice, prop2::lattice, model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    #@unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack FS, FT, FL, rS, rT, rL, k2, k4, k6, k2inv, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    #p2 = k2values   # pedagogical book-keeping
    #q2 = k2values
    #k2 = k2values
    #q4 = k4values
    #k4 = k4values
    #p2inv = invk2values

    G1T .= prop1
    G2T .= prop2
    G1T[1] = 0.  # set q=0 contribution to zero (handled in exception)
    G2T[1] = 0.  # set k=0 contribution to zero (handled in exception)

    #res .*= 0
    #res .+= p2inv .* ( fftplan * ( (bfftplan*(q4 .* G1T)) .* (bfftplan*(G2T))  .-  2*(bfftplan*(q2 .* G1T)) .* (bfftplan*(k2 .* G2T))  .+  (bfftplan*(G1T)) .* (bfftplan*(k4 .* G2T)) ) )
    res .*= 0
    res .+= k2inv .* ( fftplan * ( (bfftplan*(k4 .* G1T)) .* (bfftplan*(G2T))  .-  2*(bfftplan*(k2 .* G1T)) .* (bfftplan*(k2 .* G2T))  .+  (bfftplan*(G1T)) .* (bfftplan*(k4 .* G2T)) ) )
    
    res .*= disc.ivol^2 # one volume factor for each bfft
    
    # exception handling:
    # res .+= p2 .* G1T[1] .* G2T    # q=0
    # res .+= p2 .* G1T .* G2T[1]    # k=0
    # res[1] = 0                     # p=0
    res .+= k2 .* G1T[1] .* G2T    # q=0
    res .+= k2 .* G1T .* G2T[1]    # k=0
    res[1] = 0                     # p=0


    res .*= (-0.5) * model.g^2
end
export calcNTp6term!
function calcNTp6term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    #@unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack k2, k4, k6, k2inv, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    #p6 = k6values
    #q2inv = invk2values
    #k2inv = invk2values

    # tmpdata.tmp .*= 0 
    # tmpdata.tmp .+= (bfftplan*(q2inv .* G1T)) .* (bfftplan*(k2inv .* G2T))
    # tmpdata.tmp .-= (bfftplan*(q2inv .* G1T)) .* (bfftplan*(k2inv .* G2L))
    # tmpdata.tmp .-= (bfftplan*(q2inv .* G1L)) .* (bfftplan*(k2inv .* G2T))
    # tmpdata.tmp .+= (bfftplan*(q2inv .* G1L)) .* (bfftplan*(k2inv .* G2L))

    tmpdata.tmp .*= 0 
    tmpdata.tmp .+= (bfftplan*(k2inv .* G1T)) .* (bfftplan*(k2inv .* G2T))
    tmpdata.tmp .-= (bfftplan*(k2inv .* G1T)) .* (bfftplan*(k2inv .* G2L))
    tmpdata.tmp .-= (bfftplan*(k2inv .* G1L)) .* (bfftplan*(k2inv .* G2T))
    tmpdata.tmp .+= (bfftplan*(k2inv .* G1L)) .* (bfftplan*(k2inv .* G2L))

    #tmpdata.res .+= 0.5 * p6 .* ( fftplan * tmpdata.tmp )
    tmpdata.res .+= 0.5 * k6 .* ( fftplan * tmpdata.tmp )
end
export calcNTp4term!
function calcNTp4term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    #@unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack k2, k4, k6, k2inv, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    #p4 = k4values
    #q2inv = invk2values
    #k2inv = invk2values

    # tmpdata.tmp .*= 0
    # tmpdata.tmp .+=  0.5*(disc.sdim-3.5)*( (bfftplan*(q2inv .* G1T)).*(bfftplan*(G2T)) + (bfftplan*(G1T)).*(bfftplan*(k2inv .* G2T)) )
    # tmpdata.tmp .+= (bfftplan*(q2inv .* G1T)) .* (bfftplan*(G2L))  +  0.5*(3.5-disc.sdim)*(bfftplan*(G1T)) .* (bfftplan*(k2inv .* G2L))
    # tmpdata.tmp .+= 0.5*(3.5-disc.sdim)*(bfftplan*(q2inv .* G1L)) .* (bfftplan*(G2T))  +  (bfftplan*(G1L)) .* (bfftplan*(k2inv .* G2T))
    # tmpdata.tmp .-= ( (bfftplan*(q2inv .* G1L)) .* (bfftplan*(G2L)) + (bfftplan*(G1L)) .* (bfftplan*(k2inv .* G2L)) )
    
    # tmpdata.res .+=  p4 .* (fftplan * tmpdata.tmp)

    tmpdata.tmp .*= 0
    tmpdata.tmp .+= 0.5*(disc.sdim-3.5)*( (bfftplan*(k2inv .* G1T)).*(bfftplan*(G2T)) + (bfftplan*(G1T)).*(bfftplan*(k2inv .* G2T)) )
    tmpdata.tmp .+= (bfftplan*(k2inv .* G1T)) .* (bfftplan*(G2L))  +  0.5*(3.5-disc.sdim)*(bfftplan*(G1T)) .* (bfftplan*(k2inv .* G2L))
    tmpdata.tmp .+= 0.5*(3.5-disc.sdim)*(bfftplan*(k2inv .* G1L)) .* (bfftplan*(G2T))  +  (bfftplan*(G1L)) .* (bfftplan*(k2inv .* G2T))
    tmpdata.tmp .-= ( (bfftplan*(k2inv .* G1L)) .* (bfftplan*(G2L)) + (bfftplan*(G1L)) .* (bfftplan*(k2inv .* G2L)) )
    
    tmpdata.res .+=  k4 .* (fftplan * tmpdata.tmp)
end
export calcNTp2term!
function calcNTp2term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    #@unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack k2, k4, k6, k2inv, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    #p2 = k2values
    #q2 = k2values
    #k2 = k2values
    #q2inv = invk2values
    #k2inv = invk2values

    # tmpdata.tmp .*= 0
    # tmpdata.tmp .+=  -(bfftplan*(G1L)).*(bfftplan*(G2L)) 
    # tmpdata.tmp .+= (0.5*( (bfftplan*(q2inv .* G1L)).*(bfftplan*(k2 .* G2L)) + (bfftplan*(q2 .* G1L)).*(bfftplan*(k2inv .* G2L)) ) )
    # tmpdata.tmp .+= ((2-1.5*disc.sdim)*(bfftplan*(G1L)).*(bfftplan*(G2T)) + (disc.sdim-1.75)*(bfftplan*(q2inv .* G1L)).*(bfftplan*(k2 .* G2T)) - 0.5*(bfftplan*(q2 .* G1L)).*(bfftplan*(k2inv .* G2T)) )
    # tmpdata.tmp .+= ((2-1.5*disc.sdim)*(bfftplan*(G1T)).*(bfftplan*(G2L)) - 0.5*(bfftplan*(q2inv .* G1T)).*(bfftplan*(k2 .* G2L)) + (disc.sdim-1.75)*(bfftplan*(q2 .* G1T)).*(bfftplan*(k2inv .* G2L)) )
    # tmpdata.tmp .+= ((3.5-3*disc.sdim)*(bfftplan*(G1T)).*(bfftplan*(G2T)) + (1.75-disc.sdim)*( (bfftplan*(q2inv .* G1T)).*(bfftplan*(k2 .* G2T)) + (bfftplan*(q2 .* G1T)).*(bfftplan*(k2inv .* G2T)) ))
    
    # tmpdata.res .+= p2 .* (fftplan  * tmpdata.tmp)

    tmpdata.tmp .*= 0
    tmpdata.tmp .+=  -(bfftplan*(G1L)).*(bfftplan*(G2L)) 
    tmpdata.tmp .+= (0.5*( (bfftplan*(k2inv .* G1L)).*(bfftplan*(k2 .* G2L)) + (bfftplan*(k2 .* G1L)).*(bfftplan*(k2inv .* G2L)) ) )
    tmpdata.tmp .+= ((2-1.5*disc.sdim)*(bfftplan*(G1L)).*(bfftplan*(G2T)) + (disc.sdim-1.75)*(bfftplan*(k2inv .* G1L)).*(bfftplan*(k2 .* G2T)) - 0.5*(bfftplan*(k2 .* G1L)).*(bfftplan*(k2inv .* G2T)) )
    tmpdata.tmp .+= ((2-1.5*disc.sdim)*(bfftplan*(G1T)).*(bfftplan*(G2L)) - 0.5*(bfftplan*(k2inv .* G1T)).*(bfftplan*(k2 .* G2L)) + (disc.sdim-1.75)*(bfftplan*(k2 .* G1T)).*(bfftplan*(k2inv .* G2L)) )
    tmpdata.tmp .+= ((3.5-3*disc.sdim)*(bfftplan*(G1T)).*(bfftplan*(G2T)) + (1.75-disc.sdim)*( (bfftplan*(k2inv .* G1T)).*(bfftplan*(k2 .* G2T)) + (bfftplan*(k2 .* G1T)).*(bfftplan*(k2inv .* G2T)) ))
    
    tmpdata.res .+= p2 .* (fftplan  * tmpdata.tmp)

end
export calcNTp0term!
function calcNTp0term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    #@unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack k2, k4, k6, k2inv, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    # q2 = k2values
    # k2 = k2values
    # q4 = k4values
    # k4 = k4values
    # q2inv = invk2values
    # k2inv = invk2values
    
    # tmpdata.tmp .*= 0
    # tmpdata.tmp .+= ( 0.5*(3*disc.sdim-5) * (bfftplan*(G1L))  .* (bfftplan*(k2 .* G2T)) + 0.5*(disc.sdim-0.5) * (bfftplan*(q2inv .* G1L))  .* (bfftplan*(k4 .* G2T)) +                0.25 * (bfftplan*(q2 .* G1L)) .* (bfftplan*(G2T)))
    # tmpdata.tmp .+= (                0.25 * (bfftplan*(G1T))  .* (bfftplan*(k2 .* G2L)) + 0.5*(3*disc.sdim-5) * (bfftplan*(G1T))           .* (bfftplan*(q2 .* G2L)) + 0.5*(0.5-disc.sdim) * (bfftplan*(q4 .* G1T)) .* (bfftplan*(k2inv .* G2L)) )
    # tmpdata.tmp .+= (0.5*(3*disc.sdim-3.5)* ((bfftplan*(G1T)) .* (bfftplan*(k2 .* G2T)) +                       (bfftplan*(q2inv .* G1T))  .* (bfftplan*(k4 .* G2T)) +                       (bfftplan*(q2 .* G1T)) .* (bfftplan*(G2T))              + (bfftplan*(q4 .* G1T)) .* (bfftplan*(k2inv .* G2T)) ))

    tmpdata.tmp .*= 0
    tmpdata.tmp .+= ( 0.5*(3*disc.sdim-5) * (bfftplan*(G1L))  .* (bfftplan*(k2 .* G2T)) + 0.5*(disc.sdim-0.5) * (bfftplan*(k2inv .* G1L))  .* (bfftplan*(k4 .* G2T)) +                0.25 * (bfftplan*(k2 .* G1L)) .* (bfftplan*(G2T)))
    tmpdata.tmp .+= (                0.25 * (bfftplan*(G1T))  .* (bfftplan*(k2 .* G2L)) + 0.5*(3*disc.sdim-5) * (bfftplan*(G1T))           .* (bfftplan*(k2 .* G2L)) + 0.5*(0.5-disc.sdim) * (bfftplan*(k4 .* G1T)) .* (bfftplan*(k2inv .* G2L)) )
    tmpdata.tmp .+= (0.5*(3*disc.sdim-3.5)* ((bfftplan*(G1T)) .* (bfftplan*(k2 .* G2T)) +                       (bfftplan*(k2inv .* G1T))  .* (bfftplan*(k4 .* G2T)) +                       (bfftplan*(k2 .* G1T)) .* (bfftplan*(G2T))              + (bfftplan*(k4 .* G1T)) .* (bfftplan*(k2inv .* G2T)) ))

    tmpdata.res .+= (fftplan * tmpdata.tmp)
end
export calcNTp2invterm!
function calcNTp2invterm!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    #@unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack k2, k4, k6, k2inv, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    # q2 = k2values
    # k2 = k2values
    # q4 = k4values
    # k4 = k4values
    # q6 = k6values
    # k6 = k6values
    # q2inv = invk2values
    # k2inv = invk2values
    # p2inv = invk2values

    # tmpdata.tmp .*= 0
    # tmpdata.tmp .+= (-0.5)*(bfftplan*(G1L)).*(bfftplan*(k4 .* G2T)) + 0.25*(bfftplan*(q2inv .* G1L)).*(bfftplan*(k6 .* G2T)) + 0.25*(bfftplan*(q2 .* G1L)).*(bfftplan*(k2 .* G2T))
    # tmpdata.tmp .+= (-0.5)*(bfftplan*(q4 .* G1T)).*(bfftplan*(G2L)) + 0.25*(bfftplan*(q6 .* G1T)).*(bfftplan*(k2inv .* G2L)) + 0.25*(bfftplan*(q2 .* G1T)).*(bfftplan*(k2 .* G2L))
    # tmpdata.tmp .+= (2-disc.sdim)*(bfftplan*(G1T)).*(bfftplan*(k4 .* G2T)) - 0.25*(bfftplan*(q2inv .* G1T)).*(bfftplan*(k6 .* G2T)) + (2*disc.sdim-3.5)*(bfftplan*(q2 .* G1T)).*(bfftplan*(k2 .* G2T)) + (2-disc.sdim)*(bfftplan*(q4 .* G1T)).*(bfftplan*(G2T)) - 0.25*(bfftplan*(q6 .* G1T)).*(bfftplan*(k2inv .* G2T))
    
    # tmpdata.res .+= p2inv .* (fftplan * tmpdata.tmp) 

    tmpdata.tmp .*= 0
    tmpdata.tmp .+= (-0.5)*(bfftplan*(G1L)).*(bfftplan*(k4 .* G2T)) + 0.25*(bfftplan*(k2inv .* G1L)).*(bfftplan*(k6 .* G2T)) + 0.25*(bfftplan*(k2 .* G1L)).*(bfftplan*(k2 .* G2T))
    tmpdata.tmp .+= (-0.5)*(bfftplan*(k4 .* G1T)).*(bfftplan*(G2L)) + 0.25*(bfftplan*(k6 .* G1T)).*(bfftplan*(k2inv .* G2L)) + 0.25*(bfftplan*(k2 .* G1T)).*(bfftplan*(k2 .* G2L))
    tmpdata.tmp .+= (2-disc.sdim)*(bfftplan*(G1T)).*(bfftplan*(k4 .* G2T)) - 0.25*(bfftplan*(k2inv .* G1T)).*(bfftplan*(k6 .* G2T)) + (2*disc.sdim-3.5)*(bfftplan*(k2 .* G1T)).*(bfftplan*(k2 .* G2T)) + (2-disc.sdim)*(bfftplan*(k4 .* G1T)).*(bfftplan*(G2T)) - 0.25*(bfftplan*(k6 .* G1T)).*(bfftplan*(k2inv .* G2T))
    
    tmpdata.res .+= k2inv .* (fftplan * tmpdata.tmp) 
end
##########################################################################################################################################################################
export calcNLp2term!
function calcNLp2term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    #@unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack k2, k4, k6, k2inv, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    #p2 = k2values
    #q2 = k2values
    #k2 = k2values
    #q2inv = invk2values
    #k2inv = invk2values

    # tmpdata.tmp .*= 0
    # tmpdata.tmp .+= (-0.25)*(bfftplan*(q2inv .* G1L)).*(bfftplan*(k2 .* G2T)) 
    # tmpdata.tmp .+= (-0.25)*(bfftplan*(q2 .* G1T)).*(bfftplan*(k2inv .* G2L)) 
    # tmpdata.tmp .+= (-0.5 )*(bfftplan*(G1T)).*(bfftplan*(G2T)) + 0.25*(bfftplan*(q2inv .* G1T)).*(bfftplan*(k2 .* G2T)) + 0.25*(bfftplan*(q2 .* G1T)).*(bfftplan*(k2inv .* G2T))

    # tmpdata.res .+= p2 .* (fftplan * tmpdata.tmp)

    tmpdata.tmp .*= 0
    tmpdata.tmp .+= (-0.25)*(bfftplan*(k2inv .* G1L)).*(bfftplan*(k2 .* G2T)) 
    tmpdata.tmp .+= (-0.25)*(bfftplan*(k2 .* G1T)).*(bfftplan*(k2inv .* G2L)) 
    tmpdata.tmp .+= (-0.5 )*(bfftplan*(G1T)).*(bfftplan*(G2T)) + 0.25*(bfftplan*(k2inv .* G1T)).*(bfftplan*(k2 .* G2T)) + 0.25*(bfftplan*(k2 .* G1T)).*(bfftplan*(k2inv .* G2T))

    tmpdata.res .+= k2 .* (fftplan * tmpdata.tmp)
end
export calcNLp0term!
function calcNLp0term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    #@unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack k2, k4, k6, k2inv, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    # q2 = k2values
    # k2 = k2values
    # q4 = k4values
    # k4 = k4values
    # q2inv = invk2values
    # k2inv = invk2values

    #tmpdata.tmp .*= 0
    #tmpdata.tmp .+= 0.5 * ( (bfftplan*(G1L)).*(bfftplan*(k2 .* G2T)) .+ (bfftplan*(q2inv .* G1L)).*(bfftplan*(k4 .* G2T)) )
    #tmpdata.tmp .+= 0.5 * ( (bfftplan*(q2 .* G1T)).*(bfftplan*(G2L)) .+ (bfftplan*(q4 .* G1T)).*(bfftplan*(k2inv .* G2L)) )
    #tmpdata.tmp .+= 0.5 * ( (bfftplan*(G1T)).*(bfftplan*(k2 .* G2T)) .- (bfftplan*(q2inv .* G1T)).*(bfftplan*(k4 .* G2T)) .+ (bfftplan*(q2 .* G1T)).*(bfftplan*(G2T)) .- (bfftplan*(q4 .* G1T)).*(bfftplan*(k2inv .* G2T)) )

    tmpdata.tmp .*= 0
    tmpdata.tmp .+= 0.5 * ( (bfftplan*(G1L)).*(bfftplan*(k2 .* G2T)) .+ (bfftplan*(k2inv .* G1L)).*(bfftplan*(k4 .* G2T)) )
    tmpdata.tmp .+= 0.5 * ( (bfftplan*(k2 .* G1T)).*(bfftplan*(G2L)) .+ (bfftplan*(k4 .* G1T)).*(bfftplan*(k2inv .* G2L)) )
    tmpdata.tmp .+= 0.5 * ( (bfftplan*(G1T)).*(bfftplan*(k2 .* G2T)) .- (bfftplan*(k2inv .* G1T)).*(bfftplan*(k4 .* G2T)) .+ (bfftplan*(k2 .* G1T)).*(bfftplan*(G2T)) .- (bfftplan*(k4 .* G1T)).*(bfftplan*(k2inv .* G2T)) )
    
    tmpdata.res .+= (fftplan * tmpdata.tmp)
end
export calcNLp2invterm!
function calcNLp2invterm!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
    #@unpack k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    @unpack k2, k4, k6, k2inv, fftplan, bfftplan = simdata
    @unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    # q2 = k2values
    # k2 = k2values
    # q4 = k4values
    # k4 = k4values
    # q6 = k6values
    # k6 = k6values
    # q2inv = invk2values
    # k2inv = invk2values
    # p2inv = invk2values

    # tmpdata.tmp .*= 0
    # tmpdata.tmp .+= 0.5*(bfftplan*(G1L)).*(bfftplan*(k4 .* G2T)) .- 0.25*(bfftplan*(q2inv .* G1L)).*(bfftplan*(k6 .* G2T)) .- 0.25*(bfftplan*(q2 .* G1L)).*(bfftplan*(k2 .* G2T))
    # tmpdata.tmp .+= (-0.25)*(bfftplan*(q2 .* G1T)).*(bfftplan*(k2 .* G2L)) .+ 0.5*(bfftplan*(q4 .* G1T)).*(bfftplan*(G2L)) .- 0.25*(bfftplan*(q6 .* G1T)).*(bfftplan*(k2inv .* G2L))
    # tmpdata.tmp .+= (disc.sdim-2)*( (bfftplan*(G1T)).*(bfftplan*(k4 .* G2T)) .+ (bfftplan*(q4 .* G1T)).*(bfftplan*(G2T)) ) .+ 0.25*( (bfftplan*(q2inv .* G1T)).*(bfftplan*(k6 .* G2T)) .+ (bfftplan*(q6 .* G1T)).*(bfftplan*(k2inv .* G2T))) .+ (2*disc.sdim-3.5)*(bfftplan*(q2 .* G1T)).*(bfftplan*(k2 .* G2T))
    
    # tmpdata.res .+= p2inv .* (fftplan * tmpdata.tmp)

    tmpdata.tmp .*= 0
    tmpdata.tmp .+= 0.5*(bfftplan*(G1L)).*(bfftplan*(k4 .* G2T)) .- 0.25*(bfftplan*(k2inv .* G1L)).*(bfftplan*(k6 .* G2T)) .- 0.25*(bfftplan*(k2 .* G1L)).*(bfftplan*(k2 .* G2T))
    tmpdata.tmp .+= (-0.25)*(bfftplan*(k2 .* G1T)).*(bfftplan*(k2 .* G2L)) .+ 0.5*(bfftplan*(k4 .* G1T)).*(bfftplan*(G2L)) .- 0.25*(bfftplan*(k6 .* G1T)).*(bfftplan*(k2inv .* G2L))
    tmpdata.tmp .+= (disc.sdim-2)*( (bfftplan*(G1T)).*(bfftplan*(k4 .* G2T)) .+ (bfftplan*(k4 .* G1T)).*(bfftplan*(G2T)) ) .+ 0.25*( (bfftplan*(k2inv .* G1T)).*(bfftplan*(k6 .* G2T)) .+ (bfftplan*(k6 .* G1T)).*(bfftplan*(k2inv .* G2T))) .+ (2*disc.sdim-3.5)*(bfftplan*(k2 .* G1T)).*(bfftplan*(k2 .* G2T))
    
    tmpdata.res .+= k2inv .* (fftplan * tmpdata.tmp)
end
export NT!
function NT!(string::String, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    #@unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    #@unpack tmp, res, G1T, G1L, G2T, G2L = tmpdata
    # pedagogical book-keeping
    #p2 = simdata.k2values
    #q2 = simdata.k2values
    
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
    #tmpdata.res .+=               (-2) * p2 .* ((disc.sdim-2)*tmpdata.G1T[1] - tmpdata.G1L[1]) .* tmpdata.G2L .+ (disc.sdim-1)*tmpdata.G1T[1] .* tmpdata.G2T                         # q=0
    #tmpdata.res .+= (-2)*(disc.sdim-1) * p2 .* (tmpdata.G1T .+ tmpdata.G1L) .* tmpdata.G2T[1]                                                                                        # k=0
    #tmpdata.res[1] = sum(                q2 .* ((disc.sdim-2)*tmpdata.G1L .* tmpdata.G2T .+ tmpdata.G1T .* ((disc.sdim-2)*tmpdata.G2L .+ 4*(disc.sdim-1)*tmpdata.G2T) ))/disc.vol    # p=0
    tmpdata.res .+=               (-2) * simdata.k2 .* ((disc.sdim-2)*tmpdata.G1T[1] - tmpdata.G1L[1]) .* tmpdata.G2L .+ (disc.sdim-1)*tmpdata.G1T[1] .* tmpdata.G2T                          # q=0
    tmpdata.res .+= (-2)*(disc.sdim-1) * simdata.k2 .* (tmpdata.G1T .+ tmpdata.G1L) .* tmpdata.G2T[1]                                                                                         # k=0
    tmpdata.res[1] = sum(                simdata.k2 .* ((disc.sdim-2)*tmpdata.G1L .* tmpdata.G2T .+ tmpdata.G1T .* ((disc.sdim-2)*tmpdata.G2L .+ 4*(disc.sdim-1)*tmpdata.G2T) ))/disc.vol     # p=0
    
    tmpdata.res .*= (-model.g^2)
end
export NL!
function NL!(string::String, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    # pedagogical book-keeping
    #p2 = simdata.k2values
    #q2 = simdata.k2values
    
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
    #tmpdata.res .+=                 p2 .* ((disc.sdim-2)*tmpdata.G1T[1] + tmpdata.G1L[1]) .* tmpdata.G2T  # q=0
    #tmpdata.res .+= (disc.sdim-1) * p2 .* tmpdata.G1T .* tmpdata.G2T[1]                                   # k=0
    #tmpdata.res[1] = sum( q2 .* ((tmpdata.G1L .* tmpdata.G2T)  .+  (tmpdata.G1T .* tmpdata.G2L)) )/disc.vol    # p=0
    tmpdata.res .+=                 simdata.k2 .* ((disc.sdim-2)*tmpdata.G1T[1] + tmpdata.G1L[1]) .* tmpdata.G2T  # q=0
    tmpdata.res .+= (disc.sdim-1) * simdata.k2 .* tmpdata.G1T .* tmpdata.G2T[1]                                   # k=0
    tmpdata.res[1] = sum( simdata.k2 .* ((tmpdata.G1L .* tmpdata.G2T)  .+  (tmpdata.G1T .* tmpdata.G2L)) )/disc.vol    # p=0
    
    tmpdata.res .*= (-model.g^2)
end
############################################################################################################################################################################################################
function getE_pot(t::Int64, model::SUNgaugeScalar, meas::MeasurementTwoPIGaugeScalar, pexp::TwoPIGaugeScalarLoopEx, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, disc::TwoPIGaugeScalarDiscretization)
    
    ### Trlog contributions
    #E_pot =  0.5 * sum( (model.Mass^2 .+ simdata.k2values) .* simdata.FS[t-1,t-1] )# / disc.vol   # OK
    #E_pot += 0.5 * (disc.sdim - 1) * sum( simdata.k2values .* simdata.FT[t-1,t-1] )# / disc.vol   # OK
    E_pot =  0.5 * sum( (model.Mass^2 .+ simdata.k2) .* simdata.FS[t-1,t-1] )# / disc.vol   # OK
    E_pot += 0.5 * (disc.sdim - 1) * sum( simdata.k2 .* simdata.FT[t-1,t-1] )# / disc.vol   # OK

    ### Hartree contributions
    E_pot += 2 * model.Lambda * sum( simdata.FS[t-1,t-1] )/disc.vol *  sum( simdata.FS[t-1,t-1] )

    ## Mixed Hartree contribution
    E_pot += 0.5 * model.g^2 * sum( (disc.sdim-1)*simdata.FT[t-1,t-1] + simdata.FL[t-1,t-1])/disc.vol * sum( simdata.FS[t-1,t-1] )
    
    ## Gauge Hartree contributions
    calcAngprdterm!(model, disc, simdata, tmpdata, t-1)
    transmassterm = model.N * model.g^2 * ( sum((disc.sdim-2 + 1/(disc.sdim-1))*simdata.FT[t-1,t-1] + (1 - 1/(disc.sdim-1))simdata.FL[t-1,t-1] )/disc.vol .- (tmpdata.res ./ (disc.sdim-1)) )
    longimassterm = model.N * model.g^2 * ( sum((disc.sdim-2)*simdata.FT[t-1,t-1] + simdata.FL[t-1,t-1])/disc.vol .- tmpdata.res )

    E_pot += 0.25 * (disc.sdim-1) * sum( transmassterm .* simdata.FT[t-1,t-1] )#/disc.vol
    E_pot += 0.25 * sum( longimassterm .* simdata.FL[t-1,t-1] )#/disc.vol

    return E_pot
end