##########################################################################################################################################################################
export addTwoLatticeArrays
function addTwoLatticeArrays(vec1::Vector{Array{Float64, 3}}, vec2::Vector{Array{Float64, 3}})
    summedarray = Vector{Array{Float64, 3}}(undef,length(vec1))
    for i in 1:length(vec1)
        summedarray[i] = vec1[i] + vec2[i] # add each element (which is a lattice) element-wise of two vectors
    end
    return summedarray
end
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
export calcSelfEnergies!
function calcSelfEnergies!(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarLOloop, simdata::TwoPIGaugeScalarSimData, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64) end
function calcSelfEnergies!(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarNLOloop, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    # fix sign of rho in following expression
    simdata.rS[t,tp] .*= thesign(t, tp, simdata.NstepsinMemory)
    simdata.rT[t,tp] .*= thesign(t, tp, simdata.NstepsinMemory)
    simdata.rL[t,tp] .*= thesign(t, tp, simdata.NstepsinMemory)

    # get self-energy contributions from each diagram
    ## statistical self-energy contributions from first element in 2-element getter functions 
    simdata.SigFS[tp] .= getMixedSunsetContrToScalar(model,simdata,tmpdata,disc,t,tp)[1]
    simdata.SigFT[tp] .= getMixedSunsetContrToTransv(model,simdata,tmpdata,disc,t,tp)[1] + getGaugeSunsetContrToTransv(model,simdata,tmpdata,disc,t,tp)[1]
    simdata.SigFL[tp] .= getMixedSunsetContrToLongit(model,simdata,tmpdata,disc,t,tp)[1] + getGaugeSunsetContrToLongit(model,simdata,tmpdata,disc,t,tp)[1]

    ## spectral self-energy contributions from second element in 2-element getter functions 
    simdata.SigrS[tp] .= getMixedSunsetContrToScalar(model,simdata,tmpdata,disc,t,tp)[2]
    simdata.SigrT[tp] .= getMixedSunsetContrToTransv(model,simdata,tmpdata,disc,t,tp)[2] + getGaugeSunsetContrToTransv(model,simdata,tmpdata,disc,t,tp)[2]
    simdata.SigrL[tp] .= getMixedSunsetContrToLongit(model,simdata,tmpdata,disc,t,tp)[2] + getGaugeSunsetContrToLongit(model,simdata,tmpdata,disc,t,tp)[2]
end
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
export getMixedSunsetContrToScalar
function getMixedSunsetContrToScalar(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    @unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    p2 = k2values   # pedagogical book-keeping
    q2 = k2values
    k2 = k2values
    p4 = k4values
    q4 = k4values
    #k4 = k4values
    #p6 = k6values
    #q6 = k6values
    #k6 = k6values
    #p2inv = invk2values
    #q2inv = invk2values
    k2inv = invk2values
    function cTarr(prop1::lattice, prop2::lattice)
        G1 = deepcopy(prop1)
        G2 = deepcopy(prop2)
        G1[1] = 0.  # This way the contributions from zero modes (both q=0 and k=0) are manually put to zero (since q2inv[1]=k2inv[1]=q2[1]=...=k4[1]=...= 0).
        G2[1] = 0.
        tmp = [(-1/4) * fftplan*( (bfftplan*(k2inv .* G1)) .* (bfftplan*(G2)) ),
                (1/2) * fftplan*( (bfftplan*(k2inv .* G1)) .* (bfftplan*(q2 .* G2)) +   (bfftplan*(G1)) .* (bfftplan*(G2)) ),
               (-1/4) * fftplan*( (bfftplan*(k2inv .* G1)) .* (bfftplan*(q4 .* G2)) - 2*(bfftplan*(G1)) .* (bfftplan*(q2 .* G2)) + (bfftplan*(k2 .* G1)) .* (bfftplan*(G2)) )
                ]
        return real(tmp) * disc.ivol^3 # one volume factor for each bfft and one for the fft
    end
    function cLarr(prop1::lattice, prop2::lattice)
        G1 = deepcopy(prop1)
        G2 = deepcopy(prop2)
        G1[1] = 0.
        G2[1] = 0.
        tmp = [(1/4) * fftplan*(      (bfftplan*(k2inv .* G1)) .* (bfftplan*(G2)) ),
               (1/2) * fftplan*(    3*(bfftplan*(G1))          .* (bfftplan*(G2))       -       (bfftplan*(k2inv .* G1)) .* (bfftplan*(q2 .* G2)) ),
                       fftplan*((1/4)*(bfftplan*(k2inv .* G1)) .* (bfftplan*(q4 .* G2)) + (9/4)*(bfftplan*(k2 .* G1))    .* (bfftplan*(G2)) - (3/2)*(bfftplan*(G1)) .* (bfftplan*(q2 .* G2)) )
                ]
        return real(tmp) * disc.ivol^3 # one volume factor for each bfft and one for the fft
    end
    cTarr_F1 = addTwoLatticeArrays( cTarr( FT[t,tp], FS[t,tp] ), (-1/4)*cTarr( rT[t,tp], rS[t,tp] ) )
    cTarr_F2 = addTwoLatticeArrays( cLarr( FL[t,tp], FS[t,tp] ), (-1/4)*cLarr( rL[t,tp], rS[t,tp] ) )
    cTarr_F  = addTwoLatticeArrays(cTarr_F1, cTarr_F2) # see eq. (1.95)

    cTarr_r1 = addTwoLatticeArrays( cTarr( FT[t,tp], rS[t,tp] ), cTarr( rT[t,tp], FS[t,tp] ) )
    cTarr_r2 = addTwoLatticeArrays( cLarr( FL[t,tp], rS[t,tp] ), cLarr( rL[t,tp], FS[t,tp] ) )
    cTarr_r  = addTwoLatticeArrays(cTarr_r1, cTarr_r2) # see eq. (1.95)

    # vanishing loop momentum exception 
    zeroloopcontSigFS = ( 4*p2 .* (FL[t,tp]    * FS[t,tp][1]  - (1/4)*rL[t,tp]    * rS[t,tp][1])
                          + p2 .* (FT[t,tp][1] * FS[t,tp]     - (1/4)*rT[t,tp][1] * rS[t,tp]) )     # q=0 case : cT_qVanishing = 0,  cL_qVanishing = 4*p2,
                                                                                                    # k=0 case : cT_kVanishing = p2, cL_kVanishing = 0
    zeroloopcontSigrS = ( 4*p2 .* (FL[t,tp]    * rS[t,tp][1]  +       rL[t,tp]    * FS[t,tp][1])
                          + p2 .* (FT[t,tp][1] * rS[t,tp]     +       rT[t,tp][1] * FS[t,tp]) )        # q=0 case above, k=0 case below

    # addition to SigFS and SigrS respectively:
    tmpdata.tmplattice[1] .= ( p4 .* cTarr_F[1] + p2 .* cTarr_F[2] + cTarr_F[3] + zeroloopcontSigFS ) 
    tmpdata.tmplattice[2] .= ( p4 .* cTarr_r[1] + p2 .* cTarr_r[2] + cTarr_r[3] + zeroloopcontSigrS )

    # vanishing external momentum exception
    tmpdata.tmplattice[1][1] = sum( q2 .* ( FL[t,tp] .* FS[t,tp] - (1/4)*rL[t,tp] .* rS[t,tp] ) ) * disc.ivol # p=0 case : cT_pVanishing = 0 and cL_pVanishing = q2
    tmpdata.tmplattice[2][1] = sum( q2 .* ( FL[t,tp] .* rS[t,tp] +       rL[t,tp] .* FS[t,tp] ) ) * disc.ivol # p=0 case

    return tmpdata.tmplattice * (-2)*model.g^2*model.CF
end
##########################################################################################################################################################################
export getMixedSunsetContrToTransv
function getMixedSunsetContrToTransv(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    @unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    #d = disc.dim
    #CF = (model.N^2-1)/(2*model.N)
    p2 = k2values   # pedagogical book-keeping
    q2 = k2values
    k2 = k2values
    #p4 = k4values
    q4 = k4values
    k4 = k4values
    #p6 = k6values
    #q6 = k6values
    #k6 = k6values
    p2inv = invk2values
    #q2inv = invk2values
    #k2inv = invk2values
    function CTarr(prop1::lattice, prop2::lattice)
        G1 = deepcopy(prop1)
        G2 = deepcopy(prop2)
        G1[1] = 0. # This way the contributions from zero modes (both q=0 and k=0) are manually put to zero (since q2inv[1]=k2inv[1]=q2[1]=...=k4[1]=...= 0).
        G2[1] = 0.
        tmp = [ (-1)*fftplan*( (bfftplan*(G1))       .* (bfftplan*(G2)) ),
                   2*fftplan*( (bfftplan*(k2 .* G1)) .* (bfftplan*(G2)) +   (bfftplan*(G1))       .* (bfftplan*(q2 .* G2)) ),
                (-1)*fftplan*( (bfftplan*(k4 .* G1)) .* (bfftplan*(G2)) - 2*(bfftplan*(k2 .* G1)) .* (bfftplan*(q2 .* G2)) + (bfftplan*(G1)) .* (bfftplan*(q4 .* G2)) )
                ]
        return real(tmp) * disc.ivol^3
    end
    CTarr_F = addTwoLatticeArrays( CTarr( FS[t,tp], FS[t,tp] ), (-1/4)*CTarr( rS[t,tp], rS[t,tp] ) )
    CTarr_r = addTwoLatticeArrays( CTarr( FS[t,tp], rS[t,tp] ),        CTarr( rS[t,tp], FS[t,tp] ) )

    # vanishing loop momentum exception: CT_qVanishing = CT_kVanishing = 0 => zero loop momentum contribution vanishes

    # addition to SigFS and SigrS respectively:
    tmpdata.tmplattice[1] .= ( p2 .* CTarr_F[1] + CTarr_F[2] + p2inv .* CTarr_F[3] ) # MixedSunsetFT = p2 .* CTarr_F[1] + CTarr_F[2] + p2inv .* CTarr_F[3] 
    tmpdata.tmplattice[2] .= ( p2 .* CTarr_r[1] + CTarr_r[2] + p2inv .* CTarr_r[3] ) # MixedSunsetrT = p2 .* CTarr_r[1] + CTarr_r[2] + p2inv .* CTarr_r[3] 

    # vanishing external momentum exception (zero mode)
    tmpdata.tmplattice[1][1] = sum( 4*q2 .* (FS[t,tp] .* FS[t,tp] - (1/4)*rS[t,tp] .* rS[t,tp]) ) * disc.ivol # p=0 case : CT_p0 = 4*q2 (FT propagator)
    tmpdata.tmplattice[2][1] = sum( 4*q2 .* (FS[t,tp] .* rS[t,tp] +       rS[t,tp] .* FS[t,tp]) ) * disc.ivol # p=0 case : CT_p0 = 4*q2 (rT propagator)

    return tmpdata.tmplattice * (-1)*model.g^2*4/(2*(disc.dim-1))
end
##########################################################################################################################################################################
export getMixedSunsetContrToLongit
function getMixedSunsetContrToLongit(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    @unpack FS, FT, FL, rS, rT, rL, SigFS, SigFT, SigFL, SigrS, SigrT, SigrL, scalarmass2, transvmass2, 
    longitmass2, k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    #d = disc.dim
    #CF = (model.N^2-1)/(2*model.N)
    p2 = k2values   # pedagogical book-keeping
    q2 = k2values
    k2 = k2values
    #p4 = k4values
    q4 = k4values
    k4 = k4values
    #p6 = k6values
    #q6 = k6values
    #k6 = k6values
    p2inv = invk2values
    #q2inv = invk2values
    #k2inv = invk2values
    function CLarr(prop1::lattice, prop2::lattice)
        G1 = deepcopy(prop1)
        G2 = deepcopy(prop2)
        G1[1] = 0. # This way the contributions from zero modes (both q=0 and k=0) are manually put to zero (since q2inv[1]=k2inv[1]=q2[1]=...=k4[1]=...= 0).
        G2[1] = 0.
        tmp = [fftplan*( (bfftplan*(k4 .* G1)) .* (bfftplan*(G2)) - 2*(bfftplan*(k2 .* G1)) .* (bfftplan*(q2 .* G2)) + (bfftplan*(G1)) .* (bfftplan*(q4 .* G2)) )]
        return real(tmp) * disc.ivol^3
    end
    CLarr_F = addTwoLatticeArrays( CLarr( FS[t,tp], FS[t,tp] ), (-1/4)*CLarr( rS[t,tp], rS[t,tp] ) )
    CLarr_r = addTwoLatticeArrays( CLarr( FS[t,tp], rS[t,tp] ),        CLarr( rS[t,tp], FS[t,tp] ) )

    # vanishing loop momentum exception
    zeroloopcontMixedSunsetFL = ( p2 .* (FS[t,tp]    .* FS[t,tp][1] - (1/4)*rS[t,tp]    .* rS[t,tp][1])
                                + p2 .* (FS[t,tp][1] .* FS[t,tp]    - (1/4)*rS[t,tp][1] .* rS[t,tp]) )  # q=0 case : CL_q0 = p2, k=0 case : CL_k0 = p2
    zeroloopcontMixedSunsetrL = ( p2 .* (FS[t,tp]    .* rS[t,tp][1] +       rS[t,tp]    .* FS[t,tp][1])
                                + p2 .* (FS[t,tp][1] .* rS[t,tp]    +       rS[t,tp][1] .* FS[t,tp]) )

    tmpdata.tmplattice[1] .=  ( p2inv .* CLarr_F[1] + zeroloopcontMixedSunsetFL )# MixedSunsetFL
    tmpdata.tmplattice[2] .=  ( p2inv .* CLarr_r[1] + zeroloopcontMixedSunsetrL )# MixedSunsetrL

    # vanishing external momentum exception (zero mode)
    tmpdata.tmplattice[1][1] = sum( 4*q2 .* ( FS[t,tp] .* FS[t,tp] - (1/4)*rS[t,tp] .* rS[t,tp] ) ) * disc.ivol # p=0 case : CL_p0 = 4*q2 (FL propagator)
    tmpdata.tmplattice[2][1] = sum( 4*q2 .* ( FS[t,tp] .* rS[t,tp] +       rS[t,tp] .* FS[t,tp] ) ) * disc.ivol # p=0 case : CL_p0 = 4*q2 (rL propagator)

    return  tmpdata.tmplattice * (-1)*model.g^2/2 
end
##########################################################################################################################################################################
export getGaugeSunsetContrToTransv
function getGaugeSunsetContrToTransv(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    @unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    d = disc.dim
    #CF = (model.N^2-1)/(2*model.N)
    p2 = k2values   # pedagogical book-keeping
    q2 = k2values
    k2 = k2values
    p4 = k4values
    q4 = k4values
    k4 = k4values
    p6 = k6values
    q6 = k6values
    k6 = k6values
    p2inv = invk2values
    q2inv = invk2values
    k2inv = invk2values
    function MTTarr(prop1::lattice, prop2::lattice)
        # case of zero loop momentum has other coefficients (depend only on external momentum) whose corresponding contribution is added further down
        # => make copy of propagator to make manual exception for zero loop momentum mode
        G1 = deepcopy(prop1)
        G2 = deepcopy(prop2)
        G1[1] = 0. # This way the contributions from zero modes (both q=0 and k=0) are manually put to zero (since q2inv[1]=k2inv[1]=q2[1]=...=k4[1]=...= 0).
        G2[1] = 0. 
        tmp = [ (-1/4)*fftplan*(         (bfftplan*(k2inv .* G1)) .* (bfftplan*(q2inv .* G2))),
                -(d-2)*fftplan*(         (bfftplan*(k2inv .* G1)) .* (bfftplan*(G2)) + (bfftplan*(G1)) .* (bfftplan*(q2inv .* G2)) ),
                       fftplan*( (3*d-4)*(bfftplan*(G1))  .* (bfftplan*(G2)) - (1/2)*(7-4*d)*( (bfftplan*(k2 .* G1)) .* (bfftplan*(q2inv .* G2)) + (bfftplan*(k2inv .* G1)) .* (bfftplan*(q2 .* G2)) ) ),
                       fftplan*( -(d-2)*((bfftplan*(k4 .* G1)) .* (bfftplan*(q2inv .* G2)) + (bfftplan*(k2inv .* G1)) .* (bfftplan*(q4 .* G2))) + (3*d-4)*((bfftplan*(k2 .* G1)) .* (bfftplan*(G2)) + (bfftplan*(G1)) .* (bfftplan*(q2 .* G2)) ) ),
                  (-1)*fftplan*(  (1/4)*((bfftplan*(k6 .* G1)) .* (bfftplan*(q2inv .* G2)) + (bfftplan*(k2inv .* G1)) .* (bfftplan*(q6 .* G2)) ) + (1/2)*(7-4*d)*(bfftplan*(k2 .* G1)) .* (bfftplan*(q2 .* G2)) 
                                 +(d-2)*((bfftplan*(k4 .* G1)) .* (bfftplan*(G2)) + (bfftplan*(G1)) .* (bfftplan*(q4 .* G2))) )
                ]
        return real(tmp) * disc.ivol^3
    end
    function MLLarr(prop1::lattice, prop2::lattice)
        G1 = deepcopy(prop1)
        G2 = deepcopy(prop2)
        G1[1] = 0.
        G2[1] = 0.
        tmp = [(-1/4)*fftplan*(   (bfftplan*(k2inv .* G1)) .* (bfftplan*(q2inv .* G2)) ),
                (1/2)*fftplan*(   (bfftplan*(k2inv .* G1)) .* (bfftplan*(G2)) + (bfftplan*(G1))       .* (bfftplan*(q2inv .* G2)) ),
                (1/4)*fftplan*( 2*(bfftplan*(G1))          .* (bfftplan*(G2)) - (bfftplan*(k2 .* G1)) .* (bfftplan*(q2inv .* G2)) - (bfftplan*(k2inv .* G1)) .* (bfftplan*(q2 .* G2)) )
                ]
        return real(tmp) * disc.ivol^3
    end
    function MTLarr(prop1::lattice, prop2::lattice)
        G1 = deepcopy(prop1)
        G2 = deepcopy(prop2)
        G1[1] = 0.
        G2[1] = 0.
        tmp = [ (1/4)*fftplan*(         (bfftplan*(k2inv .* G1)) .* (bfftplan*(q2inv .* G2)) ),
                      fftplan*(   (d-2)*(bfftplan*(k2inv .* G1)) .* (bfftplan*(G2))          - (1/2)*(bfftplan*(G1))       .* (bfftplan*(q2inv .* G2)) ),
                (1/2)*fftplan*( (7-4*d)*(bfftplan*(k2inv .* G1)) .* (bfftplan*(q2 .* G2))    + (1/2)*(bfftplan*(k2 .* G1)) .* (bfftplan*(q2inv .* G2))  + (bfftplan*(G1))       .* (bfftplan*(G2)) ),
                      fftplan*(   (d-2)*(bfftplan*(k2inv .* G1)) .* (bfftplan*(q4 .* G2))    + (1/2)*((bfftplan*(G1))      .* (bfftplan*(q2 .* G2))     - (bfftplan*(k2 .* G1)) .* (bfftplan*(G2))) ),
                (1/4)*fftplan*(         (bfftplan*(k2inv .* G1)) .* (bfftplan*(q6 .* G2))    -     2*(bfftplan*(G1))       .* (bfftplan*(q4 .* G2))     + (bfftplan*(k2 .* G1)) .* (bfftplan*(q2 .* G2)) )
                ]
        return real(tmp) * disc.ivol^3
    end

    MTTarr_F = addTwoLatticeArrays( MTTarr( FT[t,tp], FT[t,tp] ), (-1/4)*MTTarr( rT[t,tp], rT[t,tp] ) )
    MTTarr_r = addTwoLatticeArrays( MTTarr( FT[t,tp], rT[t,tp] ),        MTTarr( rT[t,tp], FT[t,tp] ) )

    MLLarr_F = addTwoLatticeArrays( MLLarr( FL[t,tp], FL[t,tp] ), (-1/4)*MLLarr( rL[t,tp], rL[t,tp] ) )
    MLLarr_r = addTwoLatticeArrays( MLLarr( FL[t,tp], rL[t,tp] ),        MLLarr( rL[t,tp], FL[t,tp] ) )

    MTLarr_F = addTwoLatticeArrays( MTLarr( FT[t,tp], FL[t,tp] ), (-1/4)*MTLarr( rT[t,tp], rL[t,tp] ) )
    MTLarr_r = addTwoLatticeArrays( MTLarr( FT[t,tp], rL[t,tp] ),        MTLarr( rT[t,tp], FL[t,tp] ) )

    MLTarr_F = addTwoLatticeArrays( MTLarr( FL[t,tp], FT[t,tp] ), (-1/4)*MTLarr( rL[t,tp], rT[t,tp] ) ) # MLTarr symmetric in q and k
    MLTarr_r = addTwoLatticeArrays( MTLarr( FL[t,tp], rT[t,tp] ),        MTLarr( rL[t,tp], FT[t,tp] ) ) # MLTarr symmetric in q and k
    
    MTT_F = ( p6 .* MTTarr_F[1] + p4 .* MTTarr_F[2] + p2 .* MTTarr_F[3] + MTTarr_F[4] + p2inv .* MTTarr_F[5] )
    MTT_r = ( p6 .* MTTarr_r[1] + p4 .* MTTarr_r[2] + p2 .* MTTarr_r[3] + MTTarr_r[4] + p2inv .* MTTarr_r[5] )
    MLL_F = ( p6 .* MLLarr_F[1] + p4 .* MLLarr_F[2] + p2 .* MLLarr_F[3] )
    MLL_r = ( p6 .* MLLarr_r[1] + p4 .* MLLarr_r[2] + p2 .* MLLarr_r[3] )
    MTL_F = ( p6 .* MTLarr_F[1] + p4 .* MTLarr_F[2] + p2 .* MTLarr_F[3] + MTLarr_F[4] + p2inv .* MTLarr_F[5] )
    MTL_r = ( p6 .* MTLarr_r[1] + p4 .* MTLarr_r[2] + p2 .* MTLarr_r[3] + MTLarr_r[4] + p2inv .* MTLarr_r[5] )
    MLT_F = ( p6 .* MLTarr_F[1] + p4 .* MLTarr_F[2] + p2 .* MLTarr_F[3] + MLTarr_F[4] + p2inv .* MLTarr_F[5] )
    MLT_r = ( p6 .* MLTarr_r[1] + p4 .* MLTarr_r[2] + p2 .* MLTarr_r[3] + MLTarr_r[4] + p2inv .* MLTarr_r[5] )

    # vanishing loop momentum exception
    #MLT_q0 = 0.
    MTL_q0 = (d-2)*p2
    MLL_q0 = p2
    MTT_q0 = 4*(d-1)*p2
    ## k=0 contributions -- MTT_k0 == MTT_q0 ===>>> factor 2 and MLL_k0 == MLL_q0 ===>>> factor 2
    MLT_k0 = (d-2)*p2
    zeroloopcontPureSunsetFT = ( MTL_q0 .* (FT[t,tp]    .* FL[t,tp][1] - (1/4)*rT[t,tp]    .* rL[t,tp][1])
                             +   MLT_k0 .* (FL[t,tp][1] .* FT[t,tp]    - (1/4)*rL[t,tp][1] .* rT[t,tp])
                             + 2*MLL_q0 .* (FL[t,tp][1] .* FL[t,tp]    - (1/4)*rL[t,tp][1] .* rL[t,tp])
                             + 2*MTT_q0 .* (FT[t,tp][1] .* FT[t,tp]    - (1/4)*rT[t,tp][1] .* rT[t,tp]) )

    zeroloopcontPureSunsetrT = (MTL_q0 .* (FT[t,tp]    .* rL[t,tp][1] + rT[t,tp]    .* FL[t,tp][1])
                            +   MLT_k0 .* (FL[t,tp][1] .* rT[t,tp]    + rL[t,tp][1] .* FT[t,tp])
                            + 2*MLL_q0 .* (FL[t,tp][1] .* rL[t,tp]    + rL[t,tp][1] .* FL[t,tp])
                            + 2*MTT_q0 .* (FT[t,tp][1] .* rT[t,tp]    + rT[t,tp][1] .* FT[t,tp]) )

    tmpdata.tmplattice[1] .= ( MTT_F + MTL_F + MLT_F + MLL_F + zeroloopcontPureSunsetFT ) # == PureSunsetFT
    tmpdata.tmplattice[2] .= ( MTT_r + MTL_r + MLT_r + MLL_r + zeroloopcontPureSunsetrT ) # == PureSunsetrT
                                
    # vanishing external momentum exception (zero mode)
    MLT_p0 = (d-2)*q2 #  MTL_p0 = (d-2)*q2 AND stat expression symmetric ===>>> factor 2 (NB not true for spectral contribution)
    #MTL_p0 = (d-2)*q2
    #MLL_p0 = 0.
    MTT_p0 = 4*(d-1)*q2
    tmpdata.tmplattice[1][1] = sum( 2*MLT_p0 .* (FL[t,tp] .* FT[t,tp] - (1/4)*rL[t,tp] .* rT[t,tp]) 
                                    + MTT_p0 .* (FT[t,tp] .* FT[t,tp] - (1/4)*rT[t,tp] .* rT[t,tp]) ) * disc.ivol       # == p0contPureSunsetFT

    tmpdata.tmplattice[2][1] = sum(   MLT_p0 .* (FL[t,tp] .* rT[t,tp] +       rL[t,tp] .* FT[t,tp]) 
                                    + MLT_p0 .* (FL[t,tp] .* rL[t,tp] +       rT[t,tp] .* FL[t,tp])
                                    + MTT_p0 .* (FT[t,tp] .* rT[t,tp] +       rT[t,tp] .* FT[t,tp]) ) * disc.ivol       # == p0contPureSunsetrT 
    
    return tmpdata.tmplattice * (-1)*model.N*model.g^2/(2*(d-1))
end
##########################################################################################################################################################################
export getGaugeSunsetContrToLongit
function getGaugeSunsetContrToLongit(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    @unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, fftplan, bfftplan = simdata
    d = disc.dim
    #CF = (model.N^2-1)/(2*model.N)
    p2 = k2values   # pedagogical book-keeping
    q2 = k2values
    k2 = k2values
    p4 = k4values
    q4 = k4values
    k4 = k4values
    p6 = k6values
    q6 = k6values
    k6 = k6values
    p2inv = invk2values
    q2inv = invk2values
    k2inv = invk2values

    function NTTarr(prop1::lattice, prop2::lattice)
        G1 = deepcopy(prop1)
        G2 = deepcopy(prop2)
        G1[1] = 0.
        G2[1] = 0.
        tmp = [       fftplan*(  (1/4)*(bfftplan*(k2 .* G1)) .* (bfftplan*(q2inv .* G2)) + (1/4)*(bfftplan*(k2inv .* G1)) .* (bfftplan*(q2 .* G2))  - (1/2)*(bfftplan*(G1)) .* (bfftplan*(G2)) ),
                (1/2)*fftplan*(  (-1)*((bfftplan*(k4 .* G1)) .* (bfftplan*(q2inv .* G2)) +       (bfftplan*(k2inv .* G1)) .* (bfftplan*(q4 .* G2))) + (bfftplan*(k2 .* G1)) .* (bfftplan*(G2)) + (bfftplan*(G1)) .* (bfftplan*(q2 .* G2)) ),
                      fftplan*( (1/4)*((bfftplan*(k6 .* G1)) .* (bfftplan*(q2inv .* G2)) +       (bfftplan*(k2inv .* G1)) .* (bfftplan*(q6 .* G2))) 
                               +(d-2)*((bfftplan*(k4 .* G1)) .* (bfftplan*(G2))          +       (bfftplan*(G1))          .* (bfftplan*(q4 .* G2))) + (1/2)*(7-4*d)*(bfftplan*(k2 .* G1)) .* (bfftplan*(q2 .* G2)) )
                ]
        return real(tmp) * disc.ivol^3
    end
    function NTLarr(prop1::lattice, G2::lattice)
        G1 = deepcopy(prop1)
        G1[1] = 0.
        tmp = [ (-1/4)*fftplan*( (bfftplan*(k2inv .* G1)) .* (bfftplan*(q2 .* G2)) ),
                 (1/2)*fftplan*( (bfftplan*(k2inv .* G1)) .* (bfftplan*(q4 .* G2)) +   (bfftplan*(G1)) .* (bfftplan*(q2 .* G2)) ),
                (-1/4)*fftplan*( (bfftplan*(k2inv .* G1)) .* (bfftplan*(q6 .* G2)) - 2*(bfftplan*(G1)) .* (bfftplan*(q4 .* G2)) + (bfftplan*(k2 .* G1)) .* (bfftplan*(q2 .* G2)) )
                ]
        return real(tmp) * disc.ivol^3
    end
    NTTarr_F = addTwoLatticeArrays( NTTarr(FT[t,tp],FT[t,tp]), (-1/4)*NTTarr(rT[t,tp],rT[t,tp]) )
    NTTarr_r = addTwoLatticeArrays( NTTarr(FT[t,tp],rT[t,tp]),        NTTarr(rT[t,tp],FT[t,tp]) )

    NTLarr_F = addTwoLatticeArrays( NTLarr(FT[t,tp],FL[t,tp]), (-1/4)*NTLarr(rT[t,tp],rL[t,tp]) )
    NTLarr_r = addTwoLatticeArrays( NTLarr(FT[t,tp],rL[t,tp]),        NTLarr(rT[t,tp],FL[t,tp]) )

    NLTarr_F = addTwoLatticeArrays( NTLarr(FL[t,tp],FT[t,tp]), (-1/4)*NTLarr(rL[t,tp],rT[t,tp]) )
    NLTarr_r = addTwoLatticeArrays( NTLarr(FL[t,tp],rT[t,tp]),        NTLarr(rL[t,tp],FT[t,tp]) )
    #
    NTT_F = p2 .* NTTarr_F[1] + NTTarr_F[2] + p2inv .* NTTarr_F[3] 
    NTT_r = p2 .* NTTarr_r[1] + NTTarr_r[2] + p2inv .* NTTarr_r[3]
    #
    NTL_F = p2 .* NTLarr_F[1] + NTLarr_F[2] + p2inv .* NTLarr_F[3]
    NTL_r = p2 .* NTLarr_r[1] + NTLarr_r[2] + p2inv .* NTLarr_r[3]
    #
    NLT_F = p2 .* NLTarr_F[1] + NLTarr_F[2] + p2inv .* NLTarr_F[3]
    NLT_r = p2 .* NLTarr_r[1] + NLTarr_r[2] + p2inv .* NLTarr_r[3]
    
    # vanishing loop momentum exception 
    NLT_q0 = p2 # == NTL_k0
    #NTL_q0 = 0.
    #NLL_q0 = 0.
    NTT_q0 = (d-2)*p2 # == NTT_k0
    ##
    #NTL_k0 = p2
    zeroloopcontPureSunsetFL = ( NLT_q0 .* (FL[t,tp]    .* FT[t,tp][1] - (1/4)*rL[t,tp]    .* rT[t,tp][1])
                             +   NLT_q0 .* (FT[t,tp][1] .* FL[t,tp]    - (1/4)*rT[t,tp][1] .* rL[t,tp])
                             + 2*NTT_q0 .* (FT[t,tp][1] .* FT[t,tp]    - (1/4)*rT[t,tp][1] .* rT[t,tp]) )

    zeroloopcontPureSunsetrL = ( NLT_q0 .* (FL[t,tp]    .* rT[t,tp][1] + rL[t,tp]    .* FT[t,tp][1])
                               + NLT_q0 .* (FT[t,tp][1] .* rL[t,tp]    + rT[t,tp][1] .* FL[t,tp])
                               + NTT_q0 .* (FT[t,tp]    .* rT[t,tp][1] + rT[t,tp]    .* FT[t,tp][1])
                               + NTT_q0 .* (FT[t,tp][1] .* rT[t,tp]    + rT[t,tp][1] .* FT[t,tp]) )

    tmpdata.tmplattice[1] = ( NTT_F + NTL_F + NLT_F + zeroloopcontPureSunsetFL ) # == PureSunsetFL, NLL=0
    tmpdata.tmplattice[2] = ( NTT_r + NTL_r + NLT_r + zeroloopcontPureSunsetrL ) # == PureSunsetrL, NLL=0

    # vanishing external momentum exception (zero mode)
    NLT_p0 = q2 # == NTL_p0 ==> statistical contribution symmetric
    #NTL_p0 = q2
    tmpdata.tmplattice[1][1] = sum( 2*NLT_p0 .* (FL[t,tp] .* FT[t,tp] - (1/4)*rL[t,tp] .* rT[t,tp]) ) * disc.ivol # == p0contPureSunsetFL
    tmpdata.tmplattice[2][1] = sum( NLT_p0 .* (FL[t,tp] .* rT[t,tp] + rL[t,tp] .* FT[t,tp])
                                  + NLT_p0 .* (FT[t,tp] .* rL[t,tp] + rT[t,tp] .* FL[t,tp]) ) * disc.ivol         # == p0contPureSunsetrL

    return tmpdata.tmplattice * (-1)*model.N*model.g^2/2 
end
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
export calcRHS_statprops!
function calcRHS_statprops!(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    tpmin = simdata.indices[1]
    tmpdata.scalarRHS[1] .= 0.5 * (simdata.SigrS[tpmin] .* simdata.FS[tpmin,tp])
    tmpdata.transvRHS[1] .= 0.5 * (simdata.SigrT[tpmin] .* simdata.FT[tpmin,tp])
    tmpdata.longitRHS[1] .= 0.5 * (simdata.SigrL[tpmin] .* simdata.FL[tpmin,tp])
    for i in (tpmin+1):t-2
        tmpdata.scalarRHS[1] .+= (simdata.SigrS[i] .* simdata.FS[i,tp])
        tmpdata.transvRHS[1] .+= (simdata.SigrT[i] .* simdata.FT[i,tp])
        tmpdata.longitRHS[1] .+= (simdata.SigrL[i] .* simdata.FL[i,tp])
    end
    if tp == tpmin
        tmpdata.scalarRHS[1] .+= 0
        tmpdata.transvRHS[1] .+= 0
        tmpdata.longitRHS[1] .+= 0
    else
        tmpdata.scalarRHS[1] .-= 0.5 * (simdata.SigFS[tpmin] .* simdata.rS[tpmin,tp] * thesign(tpmin, tp, simdata.NstepsinMemory))
        tmpdata.transvRHS[1] .-= 0.5 * (simdata.SigFT[tpmin] .* simdata.rT[tpmin,tp] * thesign(tpmin, tp, simdata.NstepsinMemory))
        tmpdata.longitRHS[1] .-= 0.5 * (simdata.SigFL[tpmin] .* simdata.rL[tpmin,tp] * thesign(tpmin, tp, simdata.NstepsinMemory))
        for i in (tpmin+1):tp-1
            tmpdata.scalarRHS[1] .-= (simdata.SigFS[i] .* simdata.rS[i,tp] * thesign(i, tp, simdata.NstepsinMemory))
            tmpdata.transvRHS[1] .-= (simdata.SigFT[i] .* simdata.rT[i,tp] * thesign(i, tp, simdata.NstepsinMemory))
            tmpdata.longitRHS[1] .-= (simdata.SigFL[i] .* simdata.rL[i,tp] * thesign(i, tp, simdata.NstepsinMemory))
        end
    end
    tmpdata.scalarRHS[1] *= disc.ivol
    tmpdata.transvRHS[1] *= disc.ivol
    tmpdata.longitRHS[1] *= disc.ivol
end
##########################################################################################################################################################################
export calcRHS_specprops!
function calcRHS_specprops!(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    tmpdata.scalarRHS[1] .= 0.
    tmpdata.transvRHS[1] .= 0.
    tmpdata.longitRHS[1] .= 0.
    if ((t-1) - tp - 1) <= 0 # 0 points
        #println( "integration boundaries are empty" )
    else
        for i in tp+1:t-2
            tmpdata.scalarRHS[1] .+= (simdata.SigrS[i] .* simdata.rS[i,tp] * thesign(i, tp, simdata.NstepsinMemory))
            tmpdata.transvRHS[1] .+= (simdata.SigrT[i] .* simdata.rT[i,tp] * thesign(i, tp, simdata.NstepsinMemory))
            tmpdata.longitRHS[1] .+= (simdata.SigrL[i] .* simdata.rL[i,tp] * thesign(i, tp, simdata.NstepsinMemory))
        end
    end
    tmpdata.scalarRHS[1] *= disc.ivol
    tmpdata.transvRHS[1] *= disc.ivol
    tmpdata.longitRHS[1] *= disc.ivol
end
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
export evolve_statprops!
function evolve_statprops!(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    simdata.FS[t,tp] .= ( 2 .- disc.dt^2 * simdata.omega2Svalues  ) .* simdata.FS[t-1,tp] - simdata.FS[t-2,tp] #+ tmpdata.scalarRHS[1] * disc.dt^3
    simdata.FT[t,tp] .= ( 2 .- disc.dt^2 * simdata.omega2Tvalues  ) .* simdata.FT[t-1,tp] - simdata.FT[t-2,tp] #+ tmpdata.transvRHS[1] * disc.dt^3
    simdata.FL[t,tp] .= ( 2 .- disc.dt^2 * simdata.omega2Lvalues  ) .* simdata.FL[t-1,tp] - simdata.FL[t-2,tp] #+ tmpdata.longitRHS[1] * disc.dt^3
end
##########################################################################################################################################################################
export evolve_specprops!
function evolve_specprops!(model::SUNgaugeScalar, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    if tp == t-1 # phi-pi commutator
        simdata.rS[t,tp] .= disc.dt * thesign(t, tp, simdata.NstepsinMemory)
        simdata.rT[t,tp] .= disc.dt * thesign(t, tp, simdata.NstepsinMemory)
        simdata.rL[t,tp] .= disc.dt * thesign(t, tp, simdata.NstepsinMemory)
    elseif tp == t 
        simdata.rS[t,tp] .= 0.
        simdata.rT[t,tp] .= 0.
        simdata.rL[t,tp] .= 0.
    else
        simdata.rS[t,tp] .= ( ( 2 .- disc.dt^2 * simdata.omega2Svalues ) .* simdata.rS[t-1,tp]*thesign(t-1, tp, simdata.NstepsinMemory) 
                                - simdata.rS[t-2,tp]*thesign(t-2, tp, simdata.NstepsinMemory) )*thesign(t, tp, simdata.NstepsinMemory)#+ tmpdata.scalarRHS[1] * disc.dt^3 ) * thesign(t, tp, simdata.NstepsinMemory) 
        #
        simdata.rT[t,tp] .= ( ( 2 .- disc.dt^2 * simdata.omega2Tvalues ) .* simdata.rT[t-1,tp]*thesign(t-1, tp, simdata.NstepsinMemory) 
                                - simdata.rT[t-2,tp]*thesign(t-2, tp, simdata.NstepsinMemory) )*thesign(t, tp, simdata.NstepsinMemory)#+ tmpdata.transvRHS[1] * disc.dt^3 ) * thesign(t, tp, simdata.NstepsinMemory) 
        #
        simdata.rL[t,tp] .= ( ( 2 .- disc.dt^2 * simdata.omega2Lvalues ) .* simdata.rL[t-1,tp]*thesign(t-1, tp, simdata.NstepsinMemory) 
                                - simdata.rL[t-2,tp]*thesign(t-2, tp, simdata.NstepsinMemory) )* thesign(t, tp, simdata.NstepsinMemory) #+ tmpdata.longitRHS[1] * disc.dt^3 ) * thesign(t, tp, simdata.NstepsinMemory) 

    end
end
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
export evolve!
function evolve!(thesolution::QFTdynamicsSolutionTwoPIGaugeScalar, tmpdata::TwoPIGaugeScalarTmpData, t)
    @unpack problem, simdata, measurearray = thesolution
    @unpack model, pexp, disc, init, reno, simsetup = problem

    expandSimData!(simdata) 
    
    #Threads.@threads for tp in simdata.indices[1]:(simdata.indices[2]-1)
        #calcSelfEnergies!(model, pexp, simdata, tmpdata, disc, t-1, tp)
    #end

    setmass2values!(model, disc, simdata, tmpdata, t-1)
    setomega2values!(simdata) 

    # evolve off-diagonal (t =/= t') components
    for tp in simdata.indices[1]:simdata.indices[2]-1
        #calcRHS_specprops!(model, simdata, tmpdata, disc, t, tp)
        evolve_specprops!(model, simdata, tmpdata, disc, t, tp)
        
        #calcRHS_statprops!(model, simdata, tmpdata, disc, t, tp)
        evolve_statprops!(model, simdata, tmpdata, disc, t, tp)

    end
    #@show simdata.FS[2,1][1] 
    
    # evolve diagonal (t == t') components
    #calcRHS_specprops!(model, simdata, tmpdata, disc, t, t)
    evolve_specprops!(model, simdata, tmpdata, disc, t, t)
    #calcRHS_statprops!(model, simdata, tmpdata, disc, t, t)
    evolve_statprops!(model, simdata, tmpdata, disc, t, t)
    
    ##### Update lastEvolStep
    simsetup.lastEvolStep = t 
end