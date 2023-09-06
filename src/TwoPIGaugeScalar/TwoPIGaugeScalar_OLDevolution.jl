# calcSelfEnergies
export calcSelfEnergiesOLD!
function calcSelfEnergiesOLD!(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarLOloop, simdata::TwoPIGaugeScalarSimData, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64) end
function calcSelfEnergiesOLD!(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarNLOloop, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpDataCPUfull, disc::TwoPIGaugeScalarDiscretization, t::Int64, tp::Int64)
    # calc self-energies from known F,r (i.e. F[t-1,x], r[t-1,x])
    @unpack FS, FT, FL, rS, rT, rL, SigFS, SigFT, SigFL, SigrS, SigrT, SigrL, scalarmass2, transvmass2, longitmass2, k2values, k4values, k6values, invk2values, omega2Svalues, omega2Tvalues, omega2Lvalues, NstepsinMemory, Nx, sdim, indices = simdata
    
    #  I think that this is fine here -- not doing anything t,tp related
    #simdata.rS = rS * thesign(t, tp, NstepsinMemory)
    #rT .*= thesign(t, tp, NstepsinMemory)
    #rL .*= thesign(t, tp, NstepsinMemory)

    d = disc.sdim+1
    CF = (model.N^2-1)/(2*model.N)
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
    ### compute coefficients in terms of loop momentum that should be multiplied with external momentum p
    # statistical coefficients
    function MTTarr(prop1::lattice, prop2::lattice)
        # case of zero loop momentum has other coefficients (depend only on external momentum) whose corresponding contribution is added further down
        # => make copy of propagator to make manual exception for zero loop momentum mode
        G1 = deepcopy(prop1)
        G2 = deepcopy(prop2)
        G1[1] = 0. # This way the contributions from zero modes (both q=0 and k=0) are manually put to zero (since q2inv[1]=k2inv[1]=q2[1]=...=k4[1]=...= 0).
        G2[1] = 0. 
        tmp = [-(1/4)*fft(bfft(k2inv.*G1).*bfft(q2inv.*G2)),
                -(d-2)*fft(bfft(k2inv.*G1).*bfft(G2) + bfft(G1).*bfft(q2inv.*G2)),
                fft((3*d-4)*(bfft(G1).*bfft(G2)) - (1/2)*(7-4*d)*(bfft(k2.*G1).*bfft(q2inv.*G2) + bfft(k2inv.*G1).*bfft(q2.*G2))),
                fft(-(d-2)*(bfft(k4.*G1).*bfft(q2inv.*G2) + bfft(k2inv.*G1).*bfft(q4.*G2)) + (3*d-4)*(bfft(k2.*G1).*bfft(G2)+bfft(G1).*bfft(q2.*G2))),
                -fft((1/4)*(bfft(k6.*G1).*bfft(q2inv.*G2)+bfft(k2inv.*G1).*bfft(q6.*G2)) + (1/2)*(7-4*d)*bfft(k2.*G1).*bfft(q2.*G2) 
                + (d-2)*(bfft(k4.*G1).*bfft(G2)+bfft(G1).*bfft(q4.*G2)))
                ]
        return real(tmp)
    end
    function MLLarr(prop1::lattice, prop2::lattice)
        G1 = deepcopy(prop1)
        G2 = deepcopy(prop2)
        G1[1] = 0.
        G2[1] = 0.
        tmp = [-(1/4)*fft(bfft(k2inv.*G1).*bfft(q2inv.*G2)),
                (1/2)*fft(bfft(k2inv.*G1).*bfft(G2) + bfft(G1).*bfft(q2inv.*G2)),
                (1/4)*fft(2*bfft(G1).*bfft(G2) - bfft(k2.*G1).*bfft(q2inv.*G2) - bfft(k2inv.*G1).*bfft(q2.*G2))
                ]
        return real(tmp)
    end
    function MTLarr(prop1::lattice, prop2::lattice)
        G1 = deepcopy(prop1)
        G2 = deepcopy(prop2)
        G1[1] = 0.
        G2[1] = 0.
        tmp = [(1/4)*fft(bfft(k2inv.*G1).*bfft(q2inv.*G2)),
                fft((d-2)*bfft(k2inv.*G1).*bfft(G2) - (1/2)*bfft(G1).*bfft(q2inv.*G2)),
                (1/2)*fft((7-4*d)*bfft(k2inv.*G1).*bfft(q2.*G2) + (1/2)*bfft(k2.*G1).*bfft(q2inv.*G2) + bfft(G1).*bfft(G2)),
                fft((d-2)*bfft(k2inv.*G1).*bfft(q4.*G2) + (1/2)*(bfft(G1).*bfft(q2.*G2) - bfft(k2.*G1).*bfft(G2))),
                (1/4)*fft(bfft(k6.*G1).*bfft(k2inv.*G2) - 2*bfft(G1).*bfft(q4.*G2) + bfft(k2.*G1).*bfft(q2.*G2))
                ]
        return real(tmp)
    end
    # spectral coefficients
    function NTTarr(prop1::lattice, prop2::lattice)
        G1 = deepcopy(prop1)
        G2 = deepcopy(prop2)
        G1[1] = 0.
        G2[1] = 0.
        tmp = [fft((1/4)*bfft(k2.*G1).*bfft(q2inv.*G2) + (1/4)*bfft(k2inv.*G1).*bfft(q2.*G2)  - (1/2)*bfft(G1).*bfft(G2)),
                (1/2)*fft((-1)*(bfft(k4.*G1).*bfft(q2inv.*G2)+bfft(k2inv.*G1).*bfft(q4.*G2))  + bfft(k2.*G1).*bfft(G2) + bfft(G1).*bfft(q2.*G2)),
                fft((1/4)*(bfft(k6.*G1).*bfft(q2inv.*G2)+bfft(k2inv.*G1).*bfft(q6.*G2)) 
                    + (d-2)*(bfft(k4.*G1).*bfft(G2)+bfft(G1).*bfft(q4.*G2)) + (1/2)*(7-4*d)*bfft(k2.*G1).*bfft(q2.*G2))
                ]
        return real(tmp)
    end
    function NTLarr(G1::lattice, G2::lattice)
        # All coefficients zero for zero mode in products that are being Fourier transformed (no need to make copies and set propagators zero modes to vanish).
        tmp = [(-1/4)*fft(bfft(k2inv.*G1).*bfft(q2.*G2)),
                (1/2)*fft(bfft(k2inv.*G1).*bfft(q4.*G2) + bfft(G1).*bfft(q2.*G2)),
                (-1/4)*fft(bfft(k2inv.*G1).*bfft(q6.*G2) - 2*bfft(G1).*bfft(q4.*G2) + bfft(k2.*G1).*bfft(q2.*G2))
                ]
        return real(tmp)
    end
    # mixed sunset contributions to gauge prop
    function CTarr(prop1::lattice, prop2::lattice)
        G1 = deepcopy(prop1)
        G2 = deepcopy(prop2)
        G1[1] = 0.
        G2[1] = 0.
        tmp = [-fft(bfft(G1).*bfft(G2)),
                2*fft(bfft(k2.*G1).*bfft(G2) + bfft(G1).*bfft(q2.*G2)),
                -fft(bfft(k4.*G1).*bfft(G2) - 2*bfft(k2.*G1).*bfft(q2.*G2) + bfft(G1).*bfft(q4.*G2))
                ]
        return real(tmp)
    end
    function CLarr(G1::lattice, G2::lattice)
        tmp = [fft(bfft(k4.*G1).*bfft(G2) - 2*bfft(k2.*G1).*bfft(q2.*G2) + bfft(G1).*bfft(q4.*G2))]
        #@show (k4.*G1)[1]
        return real(tmp)
    end
    # Mixed sunset contributions to scalar prop
    function cTarr(prop1::lattice, prop2::lattice)
        G1 = deepcopy(prop1)
        G2 = deepcopy(prop2)
        G1[1] = 0.
        G2[1] = 0.
        tmp = [(-1/4)*fft(bfft(k2inv.*G1).*bfft(G2)),
                (1/2)*fft(bfft(k2inv.*G1).*bfft(q2.*G2) + bfft(G1).*bfft(G2)),
               (-1/4)*fft(bfft(k2inv.*G1).*bfft(q4.*G2) - 2*bfft(G1).*bfft(q2.*G2) + bfft(k2.*G1).*bfft(G2))
                ]
        return real(tmp)
    end
    function cLarr(prop1::lattice, prop2::lattice)
        G1 = deepcopy(prop1)
        G2 = deepcopy(prop2)
        G1[1] = 0.
        G2[1] = 0.
        tmp = [(1/4)*fft(bfft(k2inv.*G1).*bfft(G2)),
               (1/2)*fft(3*bfft(G1).*bfft(G2) - bfft(k2inv.*G1).*bfft(q2.*G2)),
               fft((1/4)*bfft(k2inv.*G1).*bfft(q4.*G2) + (9/4)*bfft(k2.*G1).*bfft(G2) - (3/2)*bfft(G1).*bfft(q2.*G2))
                ]
        return real(tmp)
    end
    function addTwoLatticeArrays(vec1::Vector{Array{Float64, 3}}, vec2::Vector{Array{Float64, 3}})
        summedarray = Vector{Array{Float64, 3}}(undef,length(vec1))
        for i in 1:length(vec1)
            summedarray[i] = vec1[i] + vec2[i] # add each element (which is a lattice) element-wise of two vectors
        end
        return summedarray
    end
    ## compute vectors to be mulitplied with external momenta
    MTTarr_F = addTwoLatticeArrays( MTTarr(FT[t,tp],FT[t,tp]), (-1/4).*MTTarr(rT[t,tp] * thesign(t, tp, NstepsinMemory),rT[t,tp]) * thesign(t, tp, NstepsinMemory))
    MTTarr_r = addTwoLatticeArrays( MTTarr(FT[t,tp],rT[t,tp]) * thesign(t, tp, NstepsinMemory), MTTarr(rT[t,tp] * thesign(t, tp, NstepsinMemory),FT[t,tp]))

    MLLarr_F = addTwoLatticeArrays( MLLarr(FL[t,tp],FL[t,tp]), (-1/4).*MLLarr(rL[t,tp] * thesign(t, tp, NstepsinMemory),rL[t,tp] * thesign(t, tp, NstepsinMemory)))
    MLLarr_r = addTwoLatticeArrays( MLLarr(FL[t,tp],rL[t,tp] * thesign(t, tp, NstepsinMemory)), MLLarr(rL[t,tp] * thesign(t, tp, NstepsinMemory),FL[t,tp]))

    MTLarr_F = addTwoLatticeArrays( MTLarr(FT[t,tp],FL[t,tp]), (-1/4).*MTLarr(rT[t,tp] * thesign(t, tp, NstepsinMemory),rL[t,tp] * thesign(t, tp, NstepsinMemory)))
    MTLarr_r = addTwoLatticeArrays( MTLarr(FT[t,tp],rL[t,tp] * thesign(t, tp, NstepsinMemory)), MTLarr(rT[t,tp] * thesign(t, tp, NstepsinMemory),FL[t,tp]))

    NTTarr_F = addTwoLatticeArrays( NTTarr(FT[t,tp],FT[t,tp]), (-1/4).*NTTarr(rT[t,tp] * thesign(t, tp, NstepsinMemory),rT[t,tp] * thesign(t, tp, NstepsinMemory)))
    NTTarr_r = addTwoLatticeArrays( NTTarr(FT[t,tp],rT[t,tp] * thesign(t, tp, NstepsinMemory)), NTTarr(rT[t,tp] * thesign(t, tp, NstepsinMemory),FT[t,tp]))

    NTLarr_F = addTwoLatticeArrays( NTLarr(FT[t,tp],FL[t,tp]), (-1/4).*NTLarr(rT[t,tp] * thesign(t, tp, NstepsinMemory),rL[t,tp] * thesign(t, tp, NstepsinMemory)))
    NTLarr_r = addTwoLatticeArrays( NTLarr(FT[t,tp],rL[t,tp] * thesign(t, tp, NstepsinMemory)), NTLarr(rT[t,tp] * thesign(t, tp, NstepsinMemory),FL[t,tp]))

    CTarr_F = addTwoLatticeArrays( CTarr(FS[t,tp],FS[t,tp]), (-1/4)*CTarr(rS[t,tp] * thesign(t, tp, NstepsinMemory),rS[t,tp] * thesign(t, tp, NstepsinMemory)) )
    CTarr_r = addTwoLatticeArrays( CTarr(FS[t,tp],rS[t,tp] * thesign(t, tp, NstepsinMemory)), CTarr(rS[t,tp],FS[t,tp]) )

    CLarr_F = addTwoLatticeArrays( CLarr(FS[t,tp],FS[t,tp]), (-1/4)*CLarr(rS[t,tp] * thesign(t, tp, NstepsinMemory),rS[t,tp]) * thesign(t, tp, NstepsinMemory) )
    CLarr_r = addTwoLatticeArrays( CLarr(FS[t,tp],rS[t,tp] * thesign(t, tp, NstepsinMemory)), CLarr(rS[t,tp] * thesign(t, tp, NstepsinMemory),FS[t,tp]) )

    cTarr_F1 = addTwoLatticeArrays( cTarr(FT[t,tp],FS[t,tp]), (-1/4)*cTarr(rT[t,tp] * thesign(t, tp, NstepsinMemory),rS[t,tp] * thesign(t, tp, NstepsinMemory)) )
    cTarr_F2 = addTwoLatticeArrays( cLarr(FL[t,tp],FS[t,tp]), (-1/4)*cLarr(rL[t,tp] * thesign(t, tp, NstepsinMemory),rS[t,tp] * thesign(t, tp, NstepsinMemory)) )
    cTarr_F = addTwoLatticeArrays(cTarr_F1, cTarr_F2) # see eq. (1.95)

    cTarr_r1 = addTwoLatticeArrays( cTarr(FT[t,tp],rS[t,tp] * thesign(t, tp, NstepsinMemory)), (-1/4)*cTarr(rT[t,tp] * thesign(t, tp, NstepsinMemory),FS[t,tp]) )
    cTarr_r2 = addTwoLatticeArrays( cLarr(FL[t,tp],rS[t,tp] * thesign(t, tp, NstepsinMemory)), (-1/4)*cLarr(rL[t,tp] * thesign(t, tp, NstepsinMemory),FS[t,tp]) )
    cTarr_r = addTwoLatticeArrays(cTarr_r1, cTarr_r2) # see eq. (1.95)

    ########################################################### compute ZERO loop momentum contribution to be added on top of other contributions (q=0, p=/=0)
    #MLT_q0 = 0.
    MTL_q0 = (d-2)*p2
    MLL_q0 = p2
    MTT_q0 = 4*(d-1)*p2
    ## NB: to add k=0 contributions -- MTT/MLL/NTT/NLL identical -> multiply these terms w factor 2
    MLT_k0 = (d-2)*p2
    #####
    NLT_q0 = p2
    #NTL_q0 = 0.
    #NLL_q0 = 0.
    NTT_q0 = (d-2)*p2
    ##
    NTL_k0 = p2
    #####
    ##CT_q0 = 0. k=0 coefficients identical -> add factor 2
    CL_q0 = p2
    #CL_k0 = p2
    #####
    ##cT_q0 = 0.
    cL_q0 = 4*p2
    cT_k0 = p2
    zeroloopcontPureSunsetFT = MTL_q0 .* (FT[t,tp][1] .* FL[t,tp] .- (1/4)*rT[t,tp][1] .* rL[t,tp]) 
                          .+   MLT_k0 .* (FL[t,tp][1] .* FT[t,tp] .- (1/4)*rL[t,tp][1] .* rT[t,tp])
                          .+ 2*MLL_q0 .* (FL[t,tp][1] .* FL[t,tp] .- (1/4)*rL[t,tp][1] .* rL[t,tp])
                          .+ 2*MTT_q0 .* (FT[t,tp][1] .* FT[t,tp] .- (1/4)*rT[t,tp][1] .* rT[t,tp])

    zeroloopcontPureSunsetrT = MTL_q0 .* (FT[t,tp][1] .* rL[t,tp] * thesign(t, tp, NstepsinMemory) .+ rT[t,tp][1] * thesign(t, tp, NstepsinMemory) .* FL[t,tp])
                          .+   MLT_k0 .* (FL[t,tp][1] .* rT[t,tp] * thesign(t, tp, NstepsinMemory) .+ rL[t,tp][1] * thesign(t, tp, NstepsinMemory) .* FT[t,tp])
                          .+ 2*MLL_q0 .* (FL[t,tp][1] .* rL[t,tp] * thesign(t, tp, NstepsinMemory) .+ rL[t,tp][1] * thesign(t, tp, NstepsinMemory) .* FL[t,tp])
                          .+ 2*MTT_q0 .* (FT[t,tp][1] .* rT[t,tp] * thesign(t, tp, NstepsinMemory) .+ rT[t,tp][1] * thesign(t, tp, NstepsinMemory) .* FT[t,tp])
                            
    zeroloopcontPureSunsetFL = NLT_q0 .* (FL[t,tp][1] .* FT[t,tp] .- (1/4)*rL[t,tp][1] .* rT[t,tp])
                          .+   NTL_k0 .* (FT[t,tp][1] .* FL[t,tp] .- (1/4)*rT[t,tp][1] .* rL[t,tp])
                          .+ 2*NTT_q0 .* (FT[t,tp][1] .* FT[t,tp] .- (1/4)*rT[t,tp][1] .* rT[t,tp])

    zeroloopcontPureSunsetrL = NLT_q0 .* (FL[t,tp][1] .* rT[t,tp] * thesign(t, tp, NstepsinMemory) .+ rL[t,tp][1] * thesign(t, tp, NstepsinMemory) .* FT[t,tp])
                          .+   NTL_k0 .* (FT[t,tp][1] .* rL[t,tp] * thesign(t, tp, NstepsinMemory) .+ rT[t,tp][1] * thesign(t, tp, NstepsinMemory) .* FL[t,tp])
                          .+ 2*NTT_q0 .* (FT[t,tp][1] .* rT[t,tp] * thesign(t, tp, NstepsinMemory) .+ rT[t,tp][1] * thesign(t, tp, NstepsinMemory) .* FT[t,tp])

    zeroloopcontMixedSunsetFL = 2*CL_q0 .* (FS[t,tp][1] .* FS[t,tp] .- (1/4)*rS[t,tp][1] .* rS[t,tp])
    zeroloopcontMixedSunsetrL = 2*CL_q0 .* (FS[t,tp][1] .* rS[t,tp] * thesign(t, tp, NstepsinMemory) .+ rS[t,tp][1] * thesign(t, tp, NstepsinMemory) .* FS[t,tp])

    zeroloopcontSigFS = cL_q0 .* (FL[t,tp] .* FS[t,tp][1] .- (1/4)*rL[t,tp] .* rS[t,tp][1]) 
                     .+ cT_k0 .* (FT[t,tp] .* FS[t,tp][1] .- (1/4)*rT[t,tp] .* rS[t,tp][1])

    zeroloopcontSigrS = cL_q0 .* (FL[t,tp] .* rS[t,tp][1] * thesign(t, tp, NstepsinMemory) + rL[t,tp] * thesign(t, tp, NstepsinMemory) .* FS[t,tp][1]) 
                     .+ cT_k0 .* (FT[t,tp] .* rS[t,tp][1] * thesign(t, tp, NstepsinMemory) + rT[t,tp] * thesign(t, tp, NstepsinMemory) .* FS[t,tp][1])
    #@show zeroloopcontPureSunsetrL[2]
    ###########################################################
    ########################################################### compute NON-ZERO loop momentum contribution (q=/=0, p=/=0)
    # pure sunset contribution to gauge props
    MTT_F = p6 .* MTTarr_F[1] + p4 .* MTTarr_F[2] + p2 .* MTTarr_F[3] + MTTarr_F[4] + p2inv .* MTTarr_F[5]
    MTT_r = p6 .* MTTarr_r[1] + p4 .* MTTarr_r[2] + p2 .* MTTarr_r[3] + MTTarr_r[4] + p2inv .* MTTarr_r[5]
    #
    MLL_F = p6 .* MLLarr_F[1] + p4 .* MLLarr_F[2] + p2 .* MLLarr_F[3]
    MLL_r = p6 .* MLLarr_r[1] + p4 .* MLLarr_r[2] + p2 .* MLLarr_r[3]
    #
    MTL_F = p6 .* MTLarr_F[1] + p4 .* MTLarr_F[2] + p2 .* MTLarr_F[3] + MTLarr_F[4] + p2inv .* MTLarr_F[5]
    MTL_r = p6 .* MTLarr_r[1] + p4 .* MTLarr_r[2] + p2 .* MTLarr_r[3] + MTLarr_r[4] + p2inv .* MTLarr_r[5]
    #
    NTT_F = p2 .* NTTarr_F[1] + NTTarr_F[2] + p2inv .* NTTarr_F[3]
    NTT_r = p2 .* NTTarr_r[1] + NTTarr_r[2] + p2inv .* NTTarr_r[3]
    #
    NTL_F = p2 .* NTLarr_F[1] + NTLarr_F[2] + p2inv .* NTLarr_F[3]
    NTL_r = p2 .* NTLarr_r[1] + NTLarr_r[2] + p2inv .* NTLarr_r[3]
    #
    MLT_F = MTL_F # check if this is true
    MLT_r = MTL_r
    NLT_F = NTL_F
    NLT_r = NTL_r

    PureSunsetFT = MTT_F + MTL_F + MLT_F + MLL_F + zeroloopcontPureSunsetFT
    PureSunsetrT = MTT_r + MTL_r + MLT_r + MLL_r + zeroloopcontPureSunsetrT
    PureSunsetFL = NTT_F + NTL_F + NLT_F + zeroloopcontPureSunsetFL # NLL=0
    PureSunsetrL = NTT_r + NTL_r + NLT_r + zeroloopcontPureSunsetrL      

    # mixed sunset contribution to gauge props

    MixedSunsetFT = p2 .* CTarr_F[1] + CTarr_F[2] + p2inv .* CTarr_F[3] 
    MixedSunsetrT = p2 .* CTarr_r[1] + CTarr_r[2] + p2inv .* CTarr_r[3] 

    MixedSunsetFL =  p2inv .* CLarr_F[1] + zeroloopcontMixedSunsetFL 
    MixedSunsetrL =  p2inv .* CLarr_r[1] + zeroloopcontMixedSunsetrL 
        
    PureSunsetCoeffT  = model.N*model.g^2/(2*(d-1)) / (disc.Nx^disc.sdim)
    PureSunsetCoeffL  = model.N*model.g^2/2         / (disc.Nx^disc.sdim)
    MixedSunsetCoeffT = model.g^2*4/(2*(d-1))       / (disc.Nx^disc.sdim)
    MixedSunsetCoeffL = model.g^2*4/2               / (disc.Nx^disc.sdim)
    MixedSunsetCoeffS = 2*model.g^2*CF              / (disc.Nx^disc.sdim)
    
    SigFT[tp] = MixedSunsetCoeffT*MixedSunsetFT + PureSunsetCoeffT*PureSunsetFT
    SigrT[tp] = MixedSunsetCoeffT*MixedSunsetrT + PureSunsetCoeffT*PureSunsetrT
    SigFL[tp] = MixedSunsetCoeffL*MixedSunsetFL + PureSunsetCoeffL*PureSunsetFL
    SigrL[tp] = MixedSunsetCoeffL*MixedSunsetrL + PureSunsetCoeffL*PureSunsetrL
    
    # mixed sunset (only) contribution to scalar prop
    SigFS[tp] = p4 .* cTarr_F[1] + p2 .* cTarr_F[2] + cTarr_F[3] + zeroloopcontSigFS
    SigrS[tp] = p4 .* cTarr_r[1] + p2 .* cTarr_r[2] + cTarr_r[3] + zeroloopcontSigrS

    ###########################################################
    ### compute self-energies for special case with zero external momentum (p=0 mode): idx==1
    MLT_p0 = (d-2)*q2
    MTL_p0 = (d-2)*q2
    #MLL_p0 = 0.
    MTT_p0 = 4*(d-1)*q2
    ####
    NLT_p0 = q2
    NTL_p0 = q2
    #NLL_p0 = 0.
    #NTT_p0 = 0.
    ####
    CT_p0 = 4*q2
    CL_p0 = 4*q2
    ####
    #cT_p0 = 0.
    cL_p0 = q2
    
    p0contPureSunsetFT = sum( MLT_p0 .* (FL[t,tp] .* FT[t,tp] - (1/4)*rL[t,tp] .* rT[t,tp]) # sum over all idx since arguments depend on loop momenta
                            + MTL_p0 .* (FT[t,tp] .* FL[t,tp] - (1/4)*rT[t,tp] .* rL[t,tp])
                            + MTT_p0 .* (FT[t,tp] .* FT[t,tp] - (1/4)*rT[t,tp] .* rT[t,tp]) )
    p0contPureSunsetrT = sum( MLT_p0 .* (FL[t,tp] .* rT[t,tp] +       rL[t,tp] .* FT[t,tp])
                            + MTL_p0 .* (FL[t,tp] .* rL[t,tp] +       rT[t,tp] .* FL[t,tp])
                            + MTT_p0 .* (FT[t,tp] .* rT[t,tp] +       rT[t,tp] .* FT[t,tp]) )
    p0contPureSunsetFL = sum( NLT_p0 .* (FL[t,tp] .* FT[t,tp] - (1/4)*rL[t,tp] .* rT[t,tp])
                            + NTL_p0 .* (FT[t,tp] .* FL[t,tp] - (1/4)*rT[t,tp] .* rL[t,tp]) )
    p0contPureSunsetrL = sum( NLT_p0 .* (FL[t,tp] .* rT[t,tp] +       rL[t,tp] .* FT[t,tp])
                            + NTL_p0 .* (FT[t,tp] .* rL[t,tp] +       rT[t,tp] .* FL[t,tp]) )
    #### mixed sunset contribution to gauge
    p0contMixedSunsetFT = sum( CT_p0 .* (FS[t,tp] .* FS[t,tp] - (1/4)*rS[t,tp] .* rS[t,tp]) )
    p0contMixedSunsetrT = sum( CT_p0 .* (FS[t,tp] .* rS[t,tp] +       rS[t,tp] .* FS[t,tp]) )
    p0contMixedSunsetFL = sum( CL_p0 .* (FS[t,tp] .* FS[t,tp] - (1/4)*rS[t,tp] .* rS[t,tp]) )
    p0contMixedSunsetrL = sum( CL_p0 .* (FS[t,tp] .* rS[t,tp] +       rS[t,tp] .* FS[t,tp]) )
    #### mixed sunset contribution to scalar
    SigFS[tp][1] = sum( cL_p0 .* (FL[t,tp] .* FS[t,tp] - (1/4)*rL[t,tp] .* rS[t,tp]) ) # + cT_p0(...)=0
    SigrS[tp][1] = sum( cL_p0 .* (FL[t,tp] .* rS[t,tp] +       rL[t,tp] .* FS[t,tp]) ) # + cT_p0(...)=0  

    SigFS[tp] *= MixedSunsetCoeffS
    SigrS[tp] *= MixedSunsetCoeffS
    @show SigFS[tp][1]
    @show SigrS[tp][1] 

    SigFT[tp][1] = MixedSunsetCoeffT*p0contMixedSunsetFT + PureSunsetCoeffT*p0contPureSunsetFT
    SigrT[tp][1] = MixedSunsetCoeffT*p0contMixedSunsetrT + PureSunsetCoeffT*p0contPureSunsetrT
    SigFL[tp][1] = MixedSunsetCoeffL*p0contMixedSunsetFL + PureSunsetCoeffL*p0contPureSunsetFL
    SigrL[tp][1] = MixedSunsetCoeffL*p0contMixedSunsetrL + PureSunsetCoeffL*p0contPureSunsetrL

    return nothing
end
export calcRHS_FS!
function calcRHS_FS!(model::TwoPIGaugeScalarModel, simdata::TwoPIGaugeScalarSimData, disc::TwoPIGaugeScalarDiscretization, RHS::lattice, t::Int64, tp::Int64)
    tpmin = simdata.indices[1]
    RHS .= 0.5 * (simdata.SigrS[tpmin] .* simdata.FS[tpmin,tp])
    for z0 in (tpmin+1):t-2
        RHS .+= (simdata.SigrS[z0] .* simdata.FS[z0,tp])
    end
    if tp == tpmin
        RHS .+= 0
    else
        RHS .-= 0.5 * (simdata.SigFS[tpmin] .* simdata.rS[tpmin,tp] * thesign(tpmin, tp, simdata.NstepsinMemory) )
        for z0 in (tpmin+1):tp-1
            RHS .-= (simdata.SigFS[z0] .* simdata.rS[z0,tp] * thesign(z0, tp, simdata.NstepsinMemory))
        end
    end
    RHS *= 1/(disc.Nx^disc.sdim)
    return nothing
end
export calcRHS_FT!
function calcRHS_FT!(model::TwoPIGaugeScalarModel, simdata::TwoPIGaugeScalarSimData, disc::TwoPIGaugeScalarDiscretization, RHS::lattice, t::Int64, tp::Int64)
    tpmin = simdata.indices[1]
    RHS .= 0.5 * (simdata.SigrT[tpmin] .* simdata.FT[tpmin,tp])
    for i in (tpmin+1):t-2
        RHS .+= (simdata.SigrT[i] .* simdata.FT[i,tp])
    end
    if tp == tpmin
        RHS .+= 0
    else
        RHS .-= 0.5 * (simdata.SigFT[tpmin] .* simdata.rT[tpmin,tp] * thesign(tpmin, tp, simdata.NstepsinMemory) )
        for i in (tpmin+1):tp-1
            RHS .-= (simdata.SigFT[i] .* simdata.rT[i,tp] * thesign(i, tp, simdata.NstepsinMemory))
        end
    end
    RHS *= 1/(disc.Nx^disc.sdim)
    return nothing
end
export calcRHS_FL!
function calcRHS_FL!(model::TwoPIGaugeScalarModel, simdata::TwoPIGaugeScalarSimData, disc::TwoPIGaugeScalarDiscretization, RHS::lattice, t::Int64, tp::Int64)
    tpmin = simdata.indices[1]
    RHS .= 0.5 * (simdata.SigrL[tpmin] .* simdata.FL[tpmin,tp])
    for i in (tpmin+1):t-2
        RHS .+= (simdata.SigrL[i] .* simdata.FL[i,tp])
    end
    if tp == tpmin
        RHS .+= 0
    else
        RHS .-= 0.5 * (simdata.SigFL[tpmin] .* simdata.rL[tpmin,tp] * thesign(tpmin, tp, simdata.NstepsinMemory))
        for i in (tpmin+1):tp-1
            RHS .-= (simdata.SigFL[i] .* simdata.rL[i,tp] * thesign(i, tp, simdata.NstepsinMemory))
        end
    end
    RHS *= 1/(disc.Nx^disc.sdim)
    return nothing
end


export calcRHS_rS!
function calcRHS_rS!(model::TwoPIGaugeScalarModel, simdata::TwoPIGaugeScalarSimData, disc::TwoPIGaugeScalarDiscretization, RHS::lattice, t::Int64, tp::Int64)
    RHS .= 0.
    if ((t-1) - tp - 1) <= 0 # 0 points
        #println( "integration boundaries are empty" )
    else
        for i in tp+1:t-2
            RHS .+= (simdata.SigrS[i] .* simdata.rS[i,tp] * thesign(i, tp, simdata.NstepsinMemory))
        end
    end
    RHS *= 1/(disc.Nx^disc.sdim)
    return nothing
end

export calcRHS_rT!
function calcRHS_rT!(model::TwoPIGaugeScalarModel, simdata::TwoPIGaugeScalarSimData, disc::TwoPIGaugeScalarDiscretization, RHS::lattice, t::Int64, tp::Int64)
    RHS .= 0.
    if ((t-1) - tp - 1) <= 0 # 0 points
        #println( "integration boundaries are empty" )
    else
        for i in tp+1:t-2
            RHS .+= (simdata.SigrT[i] .* simdata.rT[i,tp] * thesign(i, tp, simdata.NstepsinMemory))
        end
    end
    RHS *= 1/(disc.Nx^disc.sdim)
    return nothing
end
export calcRHS_rL!
function calcRHS_rL!(model::TwoPIGaugeScalarModel, simdata::TwoPIGaugeScalarSimData, disc::TwoPIGaugeScalarDiscretization, RHS::lattice, t::Int64, tp::Int64)
    RHS .= 0.
    if ((t-1) - tp - 1) <= 0 
        #
    else
        for i in tp+1:t-2
            RHS .+= (simdata.SigrL[i] .* simdata.rL[i,tp] * thesign(i, tp, simdata.NstepsinMemory))
        end
    end
    RHS *= 1/(disc.Nx^disc.sdim)
    return nothing
end

export evolve_FS!
function evolve_FS!(model::TwoPIGaugeScalarModel, simdata::TwoPIGaugeScalarSimData, disc::TwoPIGaugeScalarDiscretization, RHS::lattice, t::Int64, tp::Int64)
    simdata.FS[t,tp] .= ( 2 .- disc.dt^2 * simdata.omega2Svalues  ) .* simdata.FS[t-1,tp] - simdata.FS[t-2,tp] + RHS * disc.dt^3
end
export evolve_FT!
function evolve_FT!(model::TwoPIGaugeScalarModel, simdata::TwoPIGaugeScalarSimData, disc::TwoPIGaugeScalarDiscretization, RHS::lattice, t::Int64, tp::Int64)
    simdata.FT[t,tp] .= ( 2 .- disc.dt^2 * simdata.omega2Tvalues  ) .* simdata.FT[t-1,tp] - simdata.FT[t-2,tp] + RHS * disc.dt^3
end
export evolve_FL!
function evolve_FL!(model::TwoPIGaugeScalarModel, simdata::TwoPIGaugeScalarSimData, disc::TwoPIGaugeScalarDiscretization, RHS::lattice, t::Int64, tp::Int64)
    simdata.FL[t,tp] .= ( 2 .- disc.dt^2 * simdata.omega2Lvalues ) .* simdata.FL[t-1,tp] - simdata.FL[t-2,tp] + RHS * disc.dt^3
end


export evolve_rS!
function evolve_rS!(model::TwoPIGaugeScalarModel, simdata::TwoPIGaugeScalarSimData, disc::TwoPIGaugeScalarDiscretization, RHS::lattice, t::Int64, tp::Int64)
    if tp == t-1 # phi-pi commutator
        simdata.rS[t,tp] .= disc.dt * thesign(t, tp, simdata.NstepsinMemory)
    elseif tp == t 
        simdata.rS[t,tp] .= 0.
    else
        simdata.rS[t,tp] .= thesign(t, tp, simdata.NstepsinMemory) * ( ( 2 .- disc.dt^2 * simdata.omega2Svalues ) .* simdata.rS[t-1,tp]*thesign(t-1, tp, simdata.NstepsinMemory) 
                - simdata.rS[t-2,tp]*thesign(t-2, tp, simdata.NstepsinMemory) + RHS * disc.dt^3 )
    end
    return nothing
end
export evolve_rT!
function evolve_rT!(model::TwoPIGaugeScalarModel, simdata::TwoPIGaugeScalarSimData, disc::TwoPIGaugeScalarDiscretization, RHS::lattice, t::Int64, tp::Int64)
    if tp == t-1 # phi-pi commutator
        simdata.rT[t,tp] .= disc.dt
    elseif tp == t 
        simdata.rT[t,tp] .= 0.
    else
        simdata.rT[t,tp] .= thesign(t, tp, simdata.NstepsinMemory) * ( ( 2 .- disc.dt^2 * simdata.omega2Tvalues ) .* simdata.rT[t-1,tp]*thesign(t-1, tp, simdata.NstepsinMemory) 
            - simdata.rT[t-2,tp]*thesign(t-2, tp, simdata.NstepsinMemory) + RHS * disc.dt^3 )
    end
    return nothing
end
export evolve_rL!
function evolve_rL!(model::TwoPIGaugeScalarModel, simdata::TwoPIGaugeScalarSimData, disc::TwoPIGaugeScalarDiscretization, RHS::lattice, t::Int64, tp::Int64)
    if tp == t-1 # phi-pi commutator
        simdata.rL[t,tp] .= disc.dt
    elseif tp == t 
        simdata.rL[t,tp] .= 0.
    else
        simdata.rL[t,tp] .= thesign(t, tp, simdata.NstepsinMemory) * ( ( 2 .- disc.dt^2 * simdata.omega2Lvalues ) .* simdata.rL[t-1,tp]*thesign(t-1, tp, simdata.NstepsinMemory) 
                - simdata.rL[t-2,tp]*thesign(t-2, tp, simdata.NstepsinMemory) + RHS * disc.dt^3 )
    end
    return nothing
end
#export setomega2valuesOld!
#function setomega2valuesOld!(model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization, simdata::TwoPIGaugeScalarSimData, t::Int64)
#    simdata.omega2Svalues .= simdata.k2values .+ getScalarMass2(model,disc,simdata,t)
#    simdata.omega2Tvalues .= simdata.k2values + getTransvMass2(model,disc,simdata,t)
#    simdata.omega2Lvalues .= getLongitMass2(model,disc,simdata,t)
#end
export evolveOld!
function evolveOld!(thesolution::QFTdynamicsSolutionTwoPIGaugeScalar, tmpdata::TwoPIGaugeScalarTmpData, t)
    @unpack problem, simdata, measurearray = thesolution
    @unpack model, pexp, disc, init, reno, simsetup = problem

    expandSimData!(simdata) 
    
    Threads.@threads for tp in simdata.indices[1]:(simdata.indices[2]-1)
        #calcSelfEnergies!(model, pexp, simdata, tmpdata, disc, t-1, tp)
        #calcSelfEnergiesTest!(model, pexp, simdata, tmpdata, disc, t-1, tp)
    end

    setmass2values!(model, disc, simdata, t-1)

    setomega2values!(simdata)

    
    ##### Evolve off-diagonal components (non-equal time)
    #@Threads.threads for ichunk in 1:tmpdata.nchunks
        #for tp in tmpdata.threadranges[ichunk]
        for tp in simdata.indices[1]:simdata.indices[2]-1
            ichunk =1
            # scalar
            calcRHS_rS!(model, simdata, disc, tmpdata.scalarRHS[ichunk], t, tp)
            evolve_rS!(model, simdata, disc, tmpdata.scalarRHS[ichunk], t, tp)
            calcRHS_FS!(model, simdata, disc, tmpdata.scalarRHS[ichunk], t, tp)
            evolve_FS!(model, simdata, disc, tmpdata.scalarRHS[ichunk], t, tp)
            # transverse
            calcRHS_rT!(model, simdata, disc, tmpdata.scalarRHS[ichunk], t, tp)
            evolve_rT!(model, simdata, disc, tmpdata.scalarRHS[ichunk], t, tp)
            calcRHS_FT!(model, simdata, disc, tmpdata.scalarRHS[ichunk], t, tp)
            evolve_FT!(model, simdata, disc, tmpdata.scalarRHS[ichunk], t, tp)
            ## longitudinal
            calcRHS_rL!(model, simdata, disc, tmpdata.scalarRHS[ichunk], t, tp)
            evolve_rL!(model, simdata, disc, tmpdata.scalarRHS[ichunk], t, tp)
            calcRHS_FL!(model, simdata, disc, tmpdata.scalarRHS[ichunk], t, tp)
            evolve_FL!(model, simdata, disc, tmpdata.scalarRHS[ichunk], t, tp)
        end
    #end

    ##### Evolve diagonal components (equal time)
    # scalar
    calcRHS_rS!(model, simdata, disc, tmpdata.scalarRHS[1], t, t)
    evolve_rS!(model, simdata, disc, tmpdata.scalarRHS[1], t, t)
    calcRHS_FS!(model, simdata, disc, tmpdata.scalarRHS[1], t, t)
    evolve_FS!(model, simdata, disc, tmpdata.scalarRHS[1], t, t)
    # transverse
    calcRHS_rT!(model, simdata, disc, tmpdata.scalarRHS[1], t, t)
    evolve_rT!(model, simdata, disc, tmpdata.scalarRHS[1], t, t)
    calcRHS_FT!(model, simdata, disc, tmpdata.scalarRHS[1], t, t)
    evolve_FT!(model, simdata, disc, tmpdata.scalarRHS[1], t, t)
    ## longitudinal
    calcRHS_rL!(model, simdata, disc, tmpdata.scalarRHS[1], t, t)
    evolve_rL!(model, simdata, disc, tmpdata.scalarRHS[1], t, t)
    calcRHS_FL!(model, simdata, disc, tmpdata.scalarRHS[1], t, t)
    evolve_FL!(model, simdata, disc, tmpdata.scalarRHS[1], t, t)

    ##### Update lastEvolStep
    simsetup.lastEvolStep = t 
end