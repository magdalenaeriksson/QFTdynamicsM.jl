using Parameters
using FFTW
##########################################################################################################################################################################
##########################################################################################################################################################################
## Mixed sunset contributions
#
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
    #res .= 0
    #res .+= (-1) .* k4 .* (fftplan*(  (bfftplan*(k2inv .* G1T)) .* (bfftplan*(G2T)) ))
    #res .+=    2 .* k2 .* (fftplan*( ((bfftplan*(k2inv .* G1T)) .* (bfftplan*(k2 .* G2T))) .+ ((bfftplan*(G1T)) .* (bfftplan*(G2T))) ))
    #res .+=          1 .* (fftplan*( 2 .* ((bfftplan*(G1T)) .* (bfftplan*(k2 .* G2T)))  .-  ((bfftplan*(k2inv .* G1T)) .* (bfftplan*(k4 .* G2T)))  .-  ((bfftplan*(k2 .* G1T)) .* (bfftplan*(G2T))) ))
#
    #res .*= disc.ivol^2 # one volume factor for each bfft

    res .= 0
    res .+= (-1) .* k4 .* (bfftplan*(      ( fftplan*(k2inv .* G1T)) .* (fftplan*(G2T)) ))
    res .+=    2 .* k2 .* (bfftplan*(      ((fftplan*(k2inv .* G1T)) .* (fftplan*(k2 .* G2T)))  .+  ((fftplan*(G1T))          .* (fftplan*(G2T))) ))
    res .+=          1 .* (bfftplan*( 2 .* ((fftplan*(G1T))          .* (fftplan*(k2 .* G2T)))  .-  ((fftplan*(k2inv .* G1T)) .* (fftplan*(k4 .* G2T)))  .-  ((fftplan*(k2 .* G1T)) .* (fftplan*(G2T))) ))

    res .*= disc.ivol^2 # one volume factor for each bfft
    
    # exception handling: 
    #res .+= 4*p4 .* prop1[1] .* prop2  # q=0
    #res .+= 0                          # k=0
    #res[1] = 0                         # p=0
    
    res .+= (4 .* k2 .* prop1[1] .* prop2)/disc.vol     # q=0
    #res .+= 0                                           # k=0
    res[1] = 0                                          # p=0

    res .*= -0.25 * model.g^2
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
    #res .= 0
    #res .+=         k4 .* (fftplan*( (bfftplan*(k2inv .* G1T)) .* (bfftplan*(G2T)) )       )
    #res .+= (-2) .* k2 .* (fftplan*( (bfftplan*(k2inv .* G1T)) .* (bfftplan*(k2 .* G2T)) ) )
    #res .+=          1 .* (fftplan*( (bfftplan*(k2inv .* G1T)) .* (bfftplan*(k4 .* G2T)) ) )

    res .= 0
    res .+=         k4 .* (bfftplan*( (fftplan*(k2inv .* G1T)) .* (fftplan*(G2T)) )       )
    res .+= (-2) .* k2 .* (bfftplan*( (fftplan*(k2inv .* G1T)) .* (fftplan*(k2 .* G2T)) ) )
    res .+=          1 .* (bfftplan*( (fftplan*(k2inv .* G1T)) .* (fftplan*(k4 .* G2T)) ) )
    
    res .*= disc.ivol^2 # one volume factor for each bfft
    
    # exception handling:

    #res .+= 0                                          # q=0
    #res .+= (p2 .* prop1 .* prop2[1])                  # k=0
    #res[1] = sum( q2 .* prop1 .* prop2 )/disc.vol      # p=0

    #@show sum(res)

    #res .+= 0                                                  # q=0
    res .+= (k2 .* prop1 .* prop2[1])/disc.vol                  # k=0
    res[1] = sum( k2 .* prop1 .* prop2 )/disc.vol               # p=0

    res .*= -0.25 * model.g^2
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
    res .= 0
    res .+= (-1) .* k2 .* (fftplan*( ((bfftplan*(G1T)) .* (bfftplan*(G2T)))   ))
    res .+= 2 .*          (fftplan*( ((bfftplan*(k2 .* G1T)) .* (bfftplan*(G2T)))  .+  ((bfftplan*(G1T)) .* (bfftplan*(k2 .* G2T)))   ))
    res .+= k2inv .*      (fftplan*(-((bfftplan*(k4 .* G1T)) .* (bfftplan*(G2T)))  .+  (2 .* (bfftplan*(k2 .* G1T)) .* (bfftplan*(k2 .* G2T)))  .-  ((bfftplan*(G1T)) .* (bfftplan*(k4 .* G2T)))   ))

    res .*= disc.ivol^2 # one volume factor for each bfft
    
    # exception handling:
    #res .+= 0                                          # q=0
    #res .+= 0                                          # k=0
    #res[1] = 4 * sum( q2 .* prop1 .* prop2 )/disc.vol   # p=0

    #res .+= 0                                          # q=0
    #res .+= 0                                          # k=0
    res[1] = 4 * sum( k2 .* prop1 .* prop2 )/disc.vol   # p=0

    res .*= -0.5 * model.g^2/(disc.sdim-1)
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
    res .= 0
    res .+= k2inv .* ( fftplan * ( ((bfftplan*(k4 .* G1T)) .* (bfftplan*(G2T)))  .-  (2*(bfftplan*(k2 .* G1T)) .* (bfftplan*(k2 .* G2T)))  .+  ((bfftplan*(G1T)) .* (bfftplan*(k4 .* G2T))) ) )
    
    res .*= disc.ivol^2 # one volume factor for each bfft
    
    # exception handling:
    #res .+= ( p2 .* prop1[1] .* prop2 )   # q=0
    #res .+= ( p2 .* prop1 .* prop2[1] )   # k=0
    #res[1] = 0                            # p=0

    res .+= ( k2 .* prop1[1] .* prop2 )/disc.vol   # q=0
    res .+= ( k2 .* prop1 .* prop2[1] )/disc.vol   # k=0
    #res[1] = 0                            # p=0


    res .*= -0.5 * model.g^2
end

##########################################################################################################################################################################
##########################################################################################################################################################################
## Gauge sunset contributions
#

export addNTp6term!
function addNTp6term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
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

    tmpdata.tmp .= 0 
    tmpdata.tmp .+= (bfftplan*(k2inv .* G1T)) .* (bfftplan*(k2inv .* G2T))
    tmpdata.tmp .-= (bfftplan*(k2inv .* G1T)) .* (bfftplan*(k2inv .* G2L))
    tmpdata.tmp .-= (bfftplan*(k2inv .* G1L)) .* (bfftplan*(k2inv .* G2T))
    tmpdata.tmp .+= (bfftplan*(k2inv .* G1L)) .* (bfftplan*(k2inv .* G2L))

    #tmpdata.res .+= 0.5 * p6 .* ( fftplan * tmpdata.tmp )
    tmpdata.res .+= 0.5 .* k6 .* ( fftplan * tmpdata.tmp )
#
    #tmp .= 0
    #tmp .+= k2 .* ((-1)*fftplan*( (bfftplan*(G1T)) .* (bfftplan*(G2T)) ))
    #tmp .+= 2*fftplan*( (bfftplan*(k2 .* G1T)) .* (bfftplan*(G2T))  .+  (bfftplan*(G1T)) .* (bfftplan*(k2 .* G2T)))
    #tmp .+= k2inv .* (fftplan*((-1)*(bfftplan*(k4 .* G1T)) .* (bfftplan*(G2T))  .+  2*(bfftplan*(k2 .* G1T)) .* (bfftplan*(k2 .* G2T))  .-  (bfftplan*(G1T)) .* (bfftplan*(k4 .* G2T))))
    
    #@show tmpdata.res[1]
    #@show sum(tmpdata.res)
    #@show tmpdata.res[1]
    #@show sum(tmpdata.res)
end

export addNTp4term!
function addNTp4term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
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

    tmpdata.tmp .= 0
    tmpdata.tmp .+= 0.5*(disc.sdim-3.5)*( (bfftplan*(k2inv .* G1T)).*(bfftplan*(G2T)) + (bfftplan*(G1T)).*(bfftplan*(k2inv .* G2T)) )
    tmpdata.tmp .+= (bfftplan*(k2inv .* G1T)) .* (bfftplan*(G2L))  +  0.5*(3.5-disc.sdim)*(bfftplan*(G1T)) .* (bfftplan*(k2inv .* G2L))
    tmpdata.tmp .+= 0.5*(3.5-disc.sdim)*(bfftplan*(k2inv .* G1L)) .* (bfftplan*(G2T))  +  (bfftplan*(G1L)) .* (bfftplan*(k2inv .* G2T))
    tmpdata.tmp .-= ( (bfftplan*(k2inv .* G1L)) .* (bfftplan*(G2L)) + (bfftplan*(G1L)) .* (bfftplan*(k2inv .* G2L)) )
    
    tmpdata.res .+=  k4 .* (fftplan * tmpdata.tmp)
end

export addNTp2term!
function addNTp2term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
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

    tmpdata.tmp .= 0
    tmpdata.tmp .+= (-1)*(bfftplan*(G1L)).*(bfftplan*(G2L)) 
    tmpdata.tmp .+= (0.5*( (bfftplan*(k2inv .* G1L)).*(bfftplan*(k2 .* G2L)) + (bfftplan*(k2 .* G1L)).*(bfftplan*(k2inv .* G2L)) ) )
    tmpdata.tmp .+= ((2-1.5*disc.sdim)*(bfftplan*(G1L)).*(bfftplan*(G2T)) + (disc.sdim-1.75)*(bfftplan*(k2inv .* G1L)).*(bfftplan*(k2 .* G2T)) - 0.5*(bfftplan*(k2 .* G1L)).*(bfftplan*(k2inv .* G2T)) )
    tmpdata.tmp .+= ((2-1.5*disc.sdim)*(bfftplan*(G1T)).*(bfftplan*(G2L)) - 0.5*(bfftplan*(k2inv .* G1T)).*(bfftplan*(k2 .* G2L)) + (disc.sdim-1.75)*(bfftplan*(k2 .* G1T)).*(bfftplan*(k2inv .* G2L)) )
    tmpdata.tmp .+= ((3.5-3*disc.sdim)*(bfftplan*(G1T)).*(bfftplan*(G2T)) + (1.75-disc.sdim)*( (bfftplan*(k2inv .* G1T)).*(bfftplan*(k2 .* G2T)) + (bfftplan*(k2 .* G1T)).*(bfftplan*(k2inv .* G2T)) ))
    
    tmpdata.res .+= k2 .* (fftplan  * tmpdata.tmp)

end

export addNTp0term!
function addNTp0term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
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

    tmpdata.tmp .= 0
    tmpdata.tmp .+= ( 0.5*(3*disc.sdim-5) * (bfftplan*(G1L))  .* (bfftplan*(k2 .* G2T)) + 0.5*(disc.sdim-0.5) * (bfftplan*(k2inv .* G1L))  .* (bfftplan*(k4 .* G2T)) +                0.25 * (bfftplan*(k2 .* G1L)) .* (bfftplan*(G2T)))
    tmpdata.tmp .+= (                0.25 * (bfftplan*(G1T))  .* (bfftplan*(k2 .* G2L)) + 0.5*(3*disc.sdim-5) * (bfftplan*(G1T))           .* (bfftplan*(k2 .* G2L)) + 0.5*(0.5-disc.sdim) * (bfftplan*(k4 .* G1T)) .* (bfftplan*(k2inv .* G2L)) )
    tmpdata.tmp .+= (0.5*(3*disc.sdim-3.5)* ((bfftplan*(G1T)) .* (bfftplan*(k2 .* G2T)) +                       (bfftplan*(k2inv .* G1T))  .* (bfftplan*(k4 .* G2T)) +                       (bfftplan*(k2 .* G1T)) .* (bfftplan*(G2T))              + (bfftplan*(k4 .* G1T)) .* (bfftplan*(k2inv .* G2T)) ))

    tmpdata.res .+= 1 .* (fftplan * tmpdata.tmp)
end

export addNTp2invterm!
function addNTp2invterm!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
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

    tmpdata.tmp .= 0
    tmpdata.tmp .+= (-0.5)*(bfftplan*(G1L)).*(bfftplan*(k4 .* G2T)) + 0.25*(bfftplan*(k2inv .* G1L)).*(bfftplan*(k6 .* G2T)) + 0.25*(bfftplan*(k2 .* G1L)).*(bfftplan*(k2 .* G2T))
    tmpdata.tmp .+= (-0.5)*(bfftplan*(k4 .* G1T)).*(bfftplan*(G2L)) + 0.25*(bfftplan*(k6 .* G1T)).*(bfftplan*(k2inv .* G2L)) + 0.25*(bfftplan*(k2 .* G1T)).*(bfftplan*(k2 .* G2L))
    tmpdata.tmp .+= (2-disc.sdim)*(bfftplan*(G1T)).*(bfftplan*(k4 .* G2T)) - 0.25*(bfftplan*(k2inv .* G1T)).*(bfftplan*(k6 .* G2T)) + (2*disc.sdim-3.5)*(bfftplan*(k2 .* G1T)).*(bfftplan*(k2 .* G2T)) + (2-disc.sdim)*(bfftplan*(k4 .* G1T)).*(bfftplan*(G2T)) - 0.25*(bfftplan*(k6 .* G1T)).*(bfftplan*(k2inv .* G2T))
    
    tmpdata.res .+= k2inv .* (fftplan * tmpdata.tmp) 
end

##########################################################################################################################################################################
export addNLp2term!
function addNLp2term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
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

    tmpdata.tmp .= 0
    tmpdata.tmp .+= (-0.25) .* ((bfftplan*(k2inv .* G1L)) .* (bfftplan*(k2 .* G2T)) )
    tmpdata.tmp .+= (-0.25) .* ((bfftplan*(k2 .* G1T))    .* (bfftplan*(k2inv .* G2L)) )
    tmpdata.tmp .+= (-0.5 ) .* ((bfftplan*(G1T)).*(bfftplan*(G2T)))  .+ 0.25 .* ((bfftplan*(k2inv .* G1T)) .* (bfftplan*(k2 .* G2T))) .+ 0.25 .* ((bfftplan*(k2 .* G1T)) .* (bfftplan*(k2inv .* G2T)))

    tmpdata.res .+= k2 .* (fftplan * tmpdata.tmp)
end
export addNLp0term!
function addNLp0term!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
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

    tmpdata.tmp .= 0
    tmpdata.tmp .+= 0.5 * ( (bfftplan*(G1L)).*(bfftplan*(k2 .* G2T)) .+ (bfftplan*(k2inv .* G1L)).*(bfftplan*(k4 .* G2T)) )
    tmpdata.tmp .+= 0.5 * ( (bfftplan*(k2 .* G1T)).*(bfftplan*(G2L)) .+ (bfftplan*(k4 .* G1T)).*(bfftplan*(k2inv .* G2L)) )
    tmpdata.tmp .+= 0.5 * ( (bfftplan*(G1T)).*(bfftplan*(k2 .* G2T)) .- (bfftplan*(k2inv .* G1T)).*(bfftplan*(k4 .* G2T)) .+ (bfftplan*(k2 .* G1T)).*(bfftplan*(G2T)) .- (bfftplan*(k4 .* G1T)).*(bfftplan*(k2inv .* G2T)) )
    
    tmpdata.res .+= 1 .* (fftplan * tmpdata.tmp)
end
export addNLp2invterm!
function addNLp2invterm!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization)
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

    tmpdata.tmp .= 0
    tmpdata.tmp .+= 0.5*(bfftplan*(G1L)).*(bfftplan*(k4 .* G2T)) .- 0.25*(bfftplan*(k2inv .* G1L)).*(bfftplan*(k6 .* G2T)) .- 0.25*(bfftplan*(k2 .* G1L)).*(bfftplan*(k2 .* G2T))
    tmpdata.tmp .+= (-0.25)*(bfftplan*(k2 .* G1T)).*(bfftplan*(k2 .* G2L)) .+ 0.5*(bfftplan*(k4 .* G1T)).*(bfftplan*(G2L)) .- 0.25*(bfftplan*(k6 .* G1T)).*(bfftplan*(k2inv .* G2L))
    tmpdata.tmp .+= (disc.sdim-2)*( (bfftplan*(G1T)).*(bfftplan*(k4 .* G2T)) .+ (bfftplan*(k4 .* G1T)).*(bfftplan*(G2T)) ) .+ 0.25*( (bfftplan*(k2inv .* G1T)).*(bfftplan*(k6 .* G2T)) .+ (bfftplan*(k6 .* G1T)).*(bfftplan*(k2inv .* G2T))) .+ (2*disc.sdim-3.5)*(bfftplan*(k2 .* G1T)).*(bfftplan*(k2 .* G2T))
    
    tmpdata.res .+= k2inv .* (fftplan * tmpdata.tmp)
end

##########################################################################################################################################################################
##########################################################################################################################################################################
export NT!
function NT!(string::String, model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
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
    tmpdata.G1T[1] = 0  # set q=0 contribution to zero (handled in exception)
    tmpdata.G1L[1] = 0
    tmpdata.G2T[1] = 0  # set k=0 contribution to zero (handled in exception)
    tmpdata.G2L[1] = 0 
    
    tmpdata.res .= 0

    addNTp6term!(simdata, tmpdata, disc)
    #addNTp4term!(simdata, tmpdata, disc)
    #addNTp2term!(simdata, tmpdata, disc)
    #addNTp0term!(simdata, tmpdata, disc)
    #addNTp2invterm!(simdata, tmpdata, disc)
    
    tmpdata.res .*= disc.ivol^2 # one volume factor for each bfft

    # exception handling:
    #tmpdata.res .+=               (-2) * p2 .* ((disc.sdim-2)*tmpdata.G1T[1] - tmpdata.G1L[1]) .* tmpdata.G2L .+ (disc.sdim-1)*tmpdata.G1T[1] .* tmpdata.G2T                         # q=0
    #tmpdata.res .+= (-2)*(disc.sdim-1) * p2 .* (tmpdata.G1T .+ tmpdata.G1L) .* tmpdata.G2T[1]                                                                                        # k=0
    #tmpdata.res[1] = sum(                q2 .* ((disc.sdim-2)*tmpdata.G1L .* tmpdata.G2T .+ tmpdata.G1T .* ((disc.sdim-2)*tmpdata.G2L .+ 4*(disc.sdim-1)*tmpdata.G2T) ))/disc.vol    # p=0
    tmpdata.res .+=               (-2) .* simdata.k2 .* ( ((disc.sdim-2)*tmpdata.G1T[1] + tmpdata.G1L[1]) .* tmpdata.G2L  .+  (disc.sdim-1) .* tmpdata.G1T[1] .* tmpdata.G2T )                        # q=0
    tmpdata.res .+= (-2)*(disc.sdim-1) .* simdata.k2 .* ( tmpdata.G1T .+ tmpdata.G1L ) .* tmpdata.G2T[1]                                                                                         # k=0
    tmpdata.res[1] = sum( simdata.k2 .* ((disc.sdim-2) .* tmpdata.G1L .* tmpdata.G2T  .+  tmpdata.G1T .* ((disc.sdim-2) .* tmpdata.G2L .+ 4*(disc.sdim-1) .* tmpdata.G2T) ))/disc.vol     # p=0
    
    tmpdata.res .*= model.g^2#(-model.g^2)

    #@show sum(tmpdata.res)
end
export NL!
function NL!(string::String, model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
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
    tmpdata.G1T[1] = 0  # set q=0 contribution to zero (handled in exception)
    tmpdata.G1L[1] = 0
    tmpdata.G2T[1] = 0  # set k=0 contribution to zero (handled in exception)
    tmpdata.G2L[1] = 0 

    tmpdata.res .= 0

    addNLp2term!(simdata, tmpdata, disc)
    #addNLp0term!(simdata, tmpdata, disc)
    #addNLp2invterm!(simdata, tmpdata, disc)

    tmpdata.res .*= disc.ivol^2 # one volume factor for each bfft
    
    # exception handling:
    #tmpdata.res .+=                 p2 .* ((disc.sdim-2)*tmpdata.G1T[1] + tmpdata.G1L[1]) .* tmpdata.G2T       # q=0
    #tmpdata.res .+= (disc.sdim-1) * p2 .* tmpdata.G1T .* tmpdata.G2T[1]                                        # k=0
    #tmpdata.res[1] = sum( q2 .* ((tmpdata.G1L .* tmpdata.G2T)  .+  (tmpdata.G1T .* tmpdata.G2L)) )/disc.vol    # p=0
    tmpdata.res .+=                 simdata.k2 .* ((disc.sdim-2)*tmpdata.G1T[1] + tmpdata.G1L[1]) .* tmpdata.G2T       # q=0
    tmpdata.res .+= (disc.sdim-1) * simdata.k2 .* tmpdata.G1T .* tmpdata.G2T[1]                                        # k=0
    tmpdata.res[1] = sum( simdata.k2 .* ((tmpdata.G1L .* tmpdata.G2T)  .+  (tmpdata.G1T .* tmpdata.G2L)) )/disc.vol    # p=0
    
    tmpdata.res .*= model.g^2#(-model.g^2)
    #@show sum(tmpdata.res) 
    #@show sum(simdata.SigFT[tp])
end
##########################################################################################################################################################################
export addMixedSunsetToSigFS!
function addMixedSunsetToSigFS!(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarNLOloop, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    # statistical scalar (mixed sunset)
    calc_cTHiggs!(simdata.FT[t,tp], simdata.FS[t,tp], model, simdata, tmpdata, disc)
    simdata.SigFS[tp] .+= tmpdata.res
    calc_cTHiggs!(simdata.rT[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), model, simdata, tmpdata, disc)
    simdata.SigFS[tp] .-= 0.25 .* tmpdata.res

    calc_cLHiggs!(simdata.FT[t,tp], simdata.FS[t,tp], model, simdata, tmpdata, disc)
    simdata.SigFS[tp] .+= tmpdata.res
    calc_cLHiggs!(simdata.rT[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), model, simdata, tmpdata, disc)
    simdata.SigFS[tp] .-= 0.25 .* tmpdata.res
end
export addMixedSunsetToSigrS!
function addMixedSunsetToSigrS!(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarNLOloop, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    # spectral scalar (mixed sunset)
    calc_cTHiggs!(simdata.FT[t,tp], simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), model, simdata, tmpdata, disc)
    simdata.SigrS[tp] .+= tmpdata.res
    calc_cTHiggs!(simdata.rT[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.FS[t,tp], model, simdata, tmpdata, disc)
    simdata.SigrS[tp] .+= tmpdata.res
    calc_cLHiggs!(simdata.FL[t,tp], simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), model, simdata, tmpdata, disc)
    simdata.SigrS[tp] .+= tmpdata.res
    calc_cLHiggs!(simdata.rL[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.FS[t,tp], model, simdata, tmpdata, disc)
    simdata.SigrS[tp] .+= tmpdata.res
end
export addMixedSunsetToSigFT!
function addMixedSunsetToSigFT!(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarNLOloop, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    # statistical transverse (mixed sunset)
    calc_CTGauge!(simdata.FS[t,tp], simdata.FS[t,tp], model, simdata, tmpdata, disc) 
    simdata.SigFT[tp] .+= tmpdata.res
    calc_CTGauge!(simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), model, simdata, tmpdata, disc)
    simdata.SigFT[tp] .-= 0.25 .* tmpdata.res
end
export addMixedSunsetToSigrT!
function addMixedSunsetToSigrT!(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarNLOloop, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    # spectral transverse (mixed sunset)
    calc_CTGauge!(simdata.FS[t,tp], simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), model, simdata, tmpdata, disc)  
    simdata.SigrT[tp] .+= tmpdata.res
    calc_CTGauge!(simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.FS[t,tp], model, simdata, tmpdata, disc)
    simdata.SigrT[tp] .+= tmpdata.res
end
export addMixedSunsetToSigFL!
function addMixedSunsetToSigFL!(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarNLOloop, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    # statistical longitudinal (mixed sunset)
    calc_CLGauge!(simdata.FS[t,tp], simdata.FS[t,tp], model, simdata, tmpdata, disc)
    simdata.SigFL[tp] .+= tmpdata.res
    calc_CLGauge!(simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), model, simdata, tmpdata, disc)
    simdata.SigFL[tp] .-= 0.25 .* tmpdata.res
end
export addMixedSunsetToSigrL!
function addMixedSunsetToSigrL!(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarNLOloop, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    # spectral longitudinal (mixed sunset)
    calc_CLGauge!(simdata.FS[t,tp], simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), model, simdata, tmpdata, disc)
    simdata.SigrL[tp] .+= tmpdata.res
    calc_CLGauge!(simdata.rS[t,tp]*thesign(t, tp, simdata.NstepsinMemory), simdata.FS[t,tp], model, simdata, tmpdata, disc)
    simdata.SigrL[tp] .+= tmpdata.res
end

#############################################

export addGaugeSunsetToSigFT!
function addGaugeSunsetToSigFT!(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarNLOloop, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    # statistical transverse (gauge sunset)
    NT!("FF", model, simdata, tmpdata, disc, t, tp)
    simdata.SigFT[tp] .+= tmpdata.res
    simdata.SigFTgaugesunset[tp] .+= tmpdata.res

    NT!("rr", model, simdata, tmpdata, disc, t, tp)
    simdata.SigFT[tp] .-= 0.25 .* tmpdata.res
    simdata.SigFTgaugesunset[tp] .-= 0.25 .* tmpdata.res
end

export addGaugeSunsetToSigrT!
function addGaugeSunsetToSigrT!(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarNLOloop, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    # spectral transverse (gauge sunset)
    NT!("Fr", model, simdata, tmpdata, disc, t, tp)
    simdata.SigrT[tp] .+= tmpdata.res
    simdata.SigrTgaugesunset[tp] .+= tmpdata.res

    NT!("rF", model, simdata, tmpdata, disc, t, tp)
    simdata.SigrT[tp] .+= tmpdata.res
    simdata.SigrTgaugesunset[tp] .+= tmpdata.res
end

export addGaugeSunsetToSigFL!
function addGaugeSunsetToSigFL!(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarNLOloop, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    # statistical longitudinal (gauge sunset)
    NL!("FF", model, simdata, tmpdata, disc, t, tp)
    simdata.SigFL[tp] .+= tmpdata.res
    simdata.SigFLgaugesunset[tp] .+= tmpdata.res

    NL!("rr", model, simdata, tmpdata, disc, t, tp)
    simdata.SigFL[tp] .-= 0.25 .* tmpdata.res
    simdata.SigFLgaugesunset[tp] .-= 0.25 .* tmpdata.res
end

export addGaugeSunsetToSigrL!
function addGaugeSunsetToSigrL!(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarNLOloop, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    # spectral longitudinal (gauge sunset)
    NL!("Fr", model, simdata, tmpdata, disc, t, tp)
    simdata.SigrL[tp] .+= tmpdata.res
    simdata.SigrLgaugesunset[tp] .+= tmpdata.res

    NL!("rF", model, simdata, tmpdata, disc, t, tp)
    simdata.SigrL[tp] .+= tmpdata.res
    simdata.SigrLgaugesunset[tp] .+= tmpdata.res
end


##########################################################################################################################################################################
##########################################################################################################################################################################
export calcSelfEnergies!
function calcSelfEnergies!(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarLOloop, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    simdata.SigFS[tp] .= 0
    simdata.SigFT[tp] .= 0
    simdata.SigFL[tp] .= 0
    simdata.SigrS[tp] .= 0
    simdata.SigrT[tp] .= 0
    simdata.SigrL[tp] .= 0
end

function calcSelfEnergies!(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarNLOloop, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    ## Set self-energies to zero
    simdata.SigFS[tp] .= 0
    simdata.SigrS[tp] .= 0
    simdata.SigFT[tp] .= 0
    simdata.SigrT[tp] .= 0
    simdata.SigFL[tp] .= 0
    simdata.SigrL[tp] .= 0

    simdata.SigFTgaugesunset[tp] .= 0
    simdata.SigrTgaugesunset[tp] .= 0
    simdata.SigFLgaugesunset[tp] .= 0
    simdata.SigrLgaugesunset[tp] .= 0

    ## Add self-energy contributions from each diagram
    tmpdata.res .= 0

    # Mixed sunset contributions
    addMixedSunsetToSigFS!(model, pexp, simdata, tmpdata, disc, t, tp)
    addMixedSunsetToSigrS!(model, pexp, simdata, tmpdata, disc, t, tp)
    addMixedSunsetToSigFT!(model, pexp, simdata, tmpdata, disc, t, tp)
    addMixedSunsetToSigrT!(model, pexp, simdata, tmpdata, disc, t, tp)
    addMixedSunsetToSigFL!(model, pexp, simdata, tmpdata, disc, t, tp)
    addMixedSunsetToSigrL!(model, pexp, simdata, tmpdata, disc, t, tp)

    # Pure gauge sunset contributions
    # these also compute SigFTgaugesunset, SigrTgaugesunset, SigFLgaugesunset and SigrLgaugesunset.
    #addGaugeSunsetToSigFT!(model, pexp, simdata, tmpdata, disc, t, tp)
    #addGaugeSunsetToSigrT!(model, pexp, simdata, tmpdata, disc, t, tp)
    #addGaugeSunsetToSigFL!(model, pexp, simdata, tmpdata, disc, t, tp)
    #addGaugeSunsetToSigrL!(model, pexp, simdata, tmpdata, disc, t, tp)
end

##########################################################################################################################################################################
##########################################################################################################################################################################

export calcRHS_statprops!
function calcRHS_statprops!(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    
    tpmin = simdata.indices[1]

    # first term
    tmpdata.scalarRHS[1] .= (-0.5) .* (simdata.SigrS[tpmin] .* simdata.FS[tpmin,tp])
    tmpdata.transvRHS[1] .= (-0.5) .* (simdata.SigrT[tpmin] .* simdata.FT[tpmin,tp])
    tmpdata.longitRHS[1] .= (-0.5) .* (simdata.SigrL[tpmin] .* simdata.FL[tpmin,tp])
    for i in (tpmin+1):t-2
        tmpdata.scalarRHS[1] .-= simdata.SigrS[i] .* simdata.FS[i,tp]
        tmpdata.transvRHS[1] .-= simdata.SigrT[i] .* simdata.FT[i,tp]
        tmpdata.longitRHS[1] .-= simdata.SigrL[i] .* simdata.FL[i,tp]
    end
    # second term
    tmpdata.scalarRHS[1] .+= 0.5 .* (simdata.SigFS[tpmin] .* simdata.rS[tpmin,tp] .* thesign(tpmin, tp, simdata.NstepsinMemory))
    tmpdata.transvRHS[1] .+= 0.5 .* (simdata.SigFT[tpmin] .* simdata.rT[tpmin,tp] .* thesign(tpmin, tp, simdata.NstepsinMemory))
    tmpdata.longitRHS[1] .+= 0.5 .* (simdata.SigFL[tpmin] .* simdata.rL[tpmin,tp] .* thesign(tpmin, tp, simdata.NstepsinMemory))
    for i in (tpmin+1):tp-1
        tmpdata.scalarRHS[1] .+= simdata.SigFS[i] .* simdata.rS[i,tp] .* thesign(i, tp, simdata.NstepsinMemory)
        tmpdata.transvRHS[1] .+= simdata.SigFT[i] .* simdata.rT[i,tp] .* thesign(i, tp, simdata.NstepsinMemory)
        tmpdata.longitRHS[1] .+= simdata.SigFL[i] .* simdata.rL[i,tp] .* thesign(i, tp, simdata.NstepsinMemory)
    end
    # if tp == tpmin
    #     tmpdata.scalarRHS[1] .+= 0
    #     tmpdata.transvRHS[1] .+= 0
    #     tmpdata.longitRHS[1] .+= 0
    # else
    #     tmpdata.scalarRHS[1] .+= 0.5 .* (simdata.SigFS[tpmin] .* simdata.rS[tpmin,tp] .* thesign(tpmin, tp, simdata.NstepsinMemory))
    #     tmpdata.transvRHS[1] .+= 0.5 .* (simdata.SigFT[tpmin] .* simdata.rT[tpmin,tp] .* thesign(tpmin, tp, simdata.NstepsinMemory))
    #     tmpdata.longitRHS[1] .+= 0.5 .* (simdata.SigFL[tpmin] .* simdata.rL[tpmin,tp] .* thesign(tpmin, tp, simdata.NstepsinMemory))
    #     for i in (tpmin+1):tp-1
    #         tmpdata.scalarRHS[1] .+= simdata.SigFS[i] .* simdata.rS[i,tp] .* thesign(i, tp, simdata.NstepsinMemory)
    #         tmpdata.transvRHS[1] .+= simdata.SigFT[i] .* simdata.rT[i,tp] .* thesign(i, tp, simdata.NstepsinMemory)
    #         tmpdata.longitRHS[1] .+= simdata.SigFL[i] .* simdata.rL[i,tp] .* thesign(i, tp, simdata.NstepsinMemory)
    #     end
    # end
    
end
##########################################################################################################################################################################
export calcRHS_specprops!
function calcRHS_specprops!(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    tmpdata.scalarRHS[1] .= 0
    tmpdata.transvRHS[1] .= 0
    tmpdata.longitRHS[1] .= 0
    #if ((t-1) - tp - 1) <= 0 # 0 points
    #    #println( "integration boundaries are empty" )
    #else
    #    for i in tp+1:t-2
    #        tmpdata.scalarRHS[1] .-= simdata.SigrS[i] .* simdata.rS[i,tp] * thesign(i, tp, simdata.NstepsinMemory)
    #        tmpdata.transvRHS[1] .-= simdata.SigrT[i] .* simdata.rT[i,tp] * thesign(i, tp, simdata.NstepsinMemory)
    #        tmpdata.longitRHS[1] .-= simdata.SigrL[i] .* simdata.rL[i,tp] * thesign(i, tp, simdata.NstepsinMemory)
    #    end
    #end
    for i in tp+1:t-2
        tmpdata.scalarRHS[1] .-= simdata.SigrS[i] .* simdata.rS[i,tp] .* thesign(i, tp, simdata.NstepsinMemory)
        tmpdata.transvRHS[1] .-= simdata.SigrT[i] .* simdata.rT[i,tp] .* thesign(i, tp, simdata.NstepsinMemory)
        tmpdata.longitRHS[1] .-= simdata.SigrL[i] .* simdata.rL[i,tp] .* thesign(i, tp, simdata.NstepsinMemory)
    end
    
end
##########################################################################################################################################################################
##########################################################################################################################################################################

export evolve_statprops!
function evolve_statprops!(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    # simdata.FS[t,tp] .= ( 2 .- disc.dt^2 .* simdata.omega2Svalues  ) .* simdata.FS[t-1,tp] .- simdata.FS[t-2,tp] .+ tmpdata.scalarRHS[1] .* disc.dt^3
    # simdata.FT[t,tp] .= ( 2 .- disc.dt^2 .* simdata.omega2Tvalues  ) .* simdata.FT[t-1,tp] .- simdata.FT[t-2,tp] .+ tmpdata.transvRHS[1] .* disc.dt^3
    # simdata.FL[t,tp] .= ( 2 .- disc.dt^2 .* simdata.omega2Lvalues  ) .* simdata.FL[t-1,tp] .- simdata.FL[t-2,tp] .+ tmpdata.longitRHS[1] .* disc.dt^3
    @. simdata.FS[t,tp] = ( 2 - disc.dt^2 * simdata.omega2Svalues ) * simdata.FS[t-1,tp]  -  simdata.FS[t-2,tp]  +  tmpdata.scalarRHS[1] * disc.dt^3
    @. simdata.FT[t,tp] = ( 2 - disc.dt^2 * simdata.omega2Tvalues ) * simdata.FT[t-1,tp]  -  simdata.FT[t-2,tp]  +  tmpdata.transvRHS[1] * disc.dt^3
    @. simdata.FL[t,tp] = ( 2 - disc.dt^2 * simdata.omega2Lvalues ) * simdata.FL[t-1,tp]  -  simdata.FL[t-2,tp]  +  tmpdata.longitRHS[1] * disc.dt^3
end
##########################################################################################################################################################################
export evolve_specprops!
function evolve_specprops!(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    if tp == t-1 # phi-pi commutator
        @. simdata.rS[t,tp] = disc.dt * thesign(t, tp, simdata.NstepsinMemory)
        @. simdata.rT[t,tp] = disc.dt * thesign(t, tp, simdata.NstepsinMemory)
        @. simdata.rL[t,tp] = disc.dt * thesign(t, tp, simdata.NstepsinMemory)
    elseif tp == t 
        @. simdata.rS[t,tp] = 0
        @. simdata.rT[t,tp] = 0
        @. simdata.rL[t,tp] = 0
    else
        #simdata.rS[t,tp] .= ( ( 2 .- disc.dt^2 .* simdata.omega2Svalues ) .* simdata.rS[t-1,tp]*thesign(t-1, tp, simdata.NstepsinMemory) 
        #                        - simdata.rS[t-2,tp]*thesign(t-2, tp, simdata.NstepsinMemory) + tmpdata.scalarRHS[1] * disc.dt^3 ) .* thesign(t, tp, simdata.NstepsinMemory) 
        ##
        #simdata.rT[t,tp] .= ( ( 2 .- disc.dt^2 .* simdata.omega2Tvalues ) .* simdata.rT[t-1,tp]*thesign(t-1, tp, simdata.NstepsinMemory) 
        #                        - simdata.rT[t-2,tp] .*thesign(t-2, tp, simdata.NstepsinMemory) + tmpdata.transvRHS[1] * disc.dt^3 ) .* thesign(t, tp, simdata.NstepsinMemory) 
        ##
        #simdata.rL[t,tp] .= ( ( 2 .- disc.dt^2 .* simdata.omega2Lvalues ) .* simdata.rL[t-1,tp]*thesign(t-1, tp, simdata.NstepsinMemory) 
        #                        - simdata.rL[t-2,tp].*thesign(t-2, tp, simdata.NstepsinMemory) + tmpdata.longitRHS[1] * disc.dt^3 ) .* thesign(t, tp, simdata.NstepsinMemory) 
        @. simdata.rS[t,tp] = ( 2 - disc.dt^2 * simdata.omega2Svalues ) * simdata.rS[t-1,tp] * thesign(t-1, tp, simdata.NstepsinMemory)  -  simdata.rS[t-2,tp] * thesign(t-2, tp, simdata.NstepsinMemory)  +  tmpdata.scalarRHS[1] * disc.dt^3 
        @. simdata.rT[t,tp] = ( 2 - disc.dt^2 * simdata.omega2Tvalues ) * simdata.rT[t-1,tp] * thesign(t-1, tp, simdata.NstepsinMemory)  -  simdata.rT[t-2,tp] * thesign(t-2, tp, simdata.NstepsinMemory)  +  tmpdata.transvRHS[1] * disc.dt^3
        @. simdata.rL[t,tp] = ( 2 - disc.dt^2 * simdata.omega2Lvalues ) * simdata.rL[t-1,tp] * thesign(t-1, tp, simdata.NstepsinMemory)  -  simdata.rL[t-2,tp] * thesign(t-2, tp, simdata.NstepsinMemory)  +  tmpdata.longitRHS[1] * disc.dt^3
        simdata.rS[t,tp] .*= thesign(t, tp, simdata.NstepsinMemory) 
        simdata.rT[t,tp] .*= thesign(t, tp, simdata.NstepsinMemory) 
        simdata.rL[t,tp] .*= thesign(t, tp, simdata.NstepsinMemory) 
    end
end
##########################################################################################################################################################################
##########################################################################################################################################################################
export setomega2values!
function setomega2values!(simdata::TwoPIGaugeScalarSimData)
    @. simdata.omega2Svalues = simdata.k2 + simdata.scalarmass2[1]
    @. simdata.omega2Tvalues = simdata.k2 + simdata.transvmass2
    @. simdata.omega2Lvalues = simdata.longitmass2
end
##########################################################################################################################################################################
##########################################################################################################################################################################
export evolve!
function evolve!(thesolution::QFTdynamicsSolutionTwoPIGaugeScalar, tmpdata::TwoPIGaugeScalarTmpData, t::Int64)
    @unpack problem, simdata, measurearray = thesolution
    @unpack model, pexp, disc, init, reno, simsetup = problem

    expandSimData!(simdata) 
    
    #Threads.@threads for tp in simdata.indices[1]:(simdata.indices[2]-1)
    for tp in simdata.indices[1]:(simdata.indices[2]-1)
        calcSelfEnergies!(model, pexp, simdata, tmpdata, disc, t-1, tp)
    end

    setmass2values!(model, disc, simdata, tmpdata, t-1)
    setomega2values!(simdata) 

    # Evolve off-diagonal (t =/= t') components
    for tp in simdata.indices[1]:simdata.indices[2]-1
        calcRHS_specprops!(model, simdata, tmpdata, disc, t, tp)
        evolve_specprops!(model, simdata, tmpdata, disc, t, tp)
        
        calcRHS_statprops!(model, simdata, tmpdata, disc, t, tp)
        evolve_statprops!(model, simdata, tmpdata, disc, t, tp)

    end
    
    # Evolve diagonal (t == t') components
    #calcRHS_specprops!(model, simdata, tmpdata, disc, t, t) # this is just zero anyway
    evolve_specprops!(model, simdata, tmpdata, disc, t, t)

    calcRHS_statprops!(model, simdata, tmpdata, disc, t, t)
    evolve_statprops!(model, simdata, tmpdata, disc, t, t)
    
    # Update lastEvolStep
    simsetup.lastEvolStep = t 
end