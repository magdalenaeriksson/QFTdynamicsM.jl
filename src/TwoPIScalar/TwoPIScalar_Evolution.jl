using Parameters
using FFTW
using LinearAlgebra
# calcSelfEnergies
export calcSelfEnergies!
export calcSelfEnergiessimple!
export calcSelfEnergyattp!
function calcSelfEnergies!(  model::TwoPIScalarPhi4, pexp::TwoPIScalarLambdaLO,   disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimData, tmpdata::TwoPIScalarTmpData, tmone::Int64) end
function calcSelfEnergyattp!(model::TwoPIScalarPhi4, pexp::TwoPIScalarLambdaLO,   disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimData, tmpdata::TwoPIScalarTmpData, tmone::Int64, tp::Int64, ichunk::Int64) end
function calcSelfEnergies!(  model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimData, tmpdata::TwoPIScalarTmpData, tmone::Int64) end
function calcSelfEnergyattp!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimData, tmpdata::TwoPIScalarTmpData, tmone::Int64, tp::Int64, ichunk::Int64) end

###########################################################################################################################
# Lambda expansion 
###########################################################################################################################
function calcSelfEnergies!(model::TwoPIScalarPhi4, pexp::TwoPIScalarLambdaNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPU, tmpdata::TwoPIScalarTmpDataCPU, tmone::Int64)
    @Threads.threads for ichunk in 1:tmpdata.nchunks
        for tp in tmpdata.threadranges[ichunk]
            calcSelfEnergyattp!(model, pexp, disc, simdata, tmpdata, tmone, tp, ichunk)
        end
    end
    return
end

function calcSelfEnergyattp!(model::TwoPIScalarPhi4, pexp::TwoPIScalarLambdaNLOquantum, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUfull, tmpdata::TwoPIScalarTmpDataCPUfull, tmone::Int64, tp::Int64, ichunk::Int64)
    tmpdata.rx[ichunk] = real( tmpdata.ftplan * (simdata.r[tmone,tp] * thesign(tmone, tp, simdata.NstepsinMemory)) )
    tmpdata.Fx[ichunk] = real( tmpdata.ftplan *  simdata.F[tmone,tp] )
    simdata.Sigr[tp] .= disc.ivol^3 * (-model.Lambda^2*(model.ONgroup+2)/(6 *model.ONgroup^2)) * real(tmpdata.iftplan * ( tmpdata.rx[ichunk] .* (  tmpdata.Fx[ichunk].^2 .- (1/12.) * tmpdata.rx[ichunk].^2)) ) # Normalisation: from exchange of fft/ifft and powers in selfenergy cal
    simdata.SigF[tp] .= disc.ivol^3 * (-model.Lambda^2*(model.ONgroup+2)/(18*model.ONgroup^2)) * real(tmpdata.iftplan * ( tmpdata.Fx[ichunk] .* (  tmpdata.Fx[ichunk].^2 .- 0.75    * tmpdata.rx[ichunk].^2)) ) # Normalisation: from exchange of fft/ifft and powers in selfenergy cal
    return
end

function calcSelfEnergyattp!(model::TwoPIScalarPhi4, pexp::TwoPIScalarLambdaNLO, disc::TwoPIScalarDiscretizationLattice, simdata::TwoPIScalarSimDataCPUred, tmpdata::TwoPIScalarTmpDataCPUred, tmone::Int64, tp::Int64, ichunk::Int64)
    copytolattice!( tmpdata.r[ichunk], simdata.r[tmone,tp] * thesign(tmone, tp, simdata.NstepsinMemory) , disc.fftwhelper)
    copytolattice!( tmpdata.F[ichunk], simdata.F[tmone,tp], disc.fftwhelper)
    #copytofulllattice!( tmpdata.r[ichunk], simdata.r[tmone,tp] * thesign(tmone, tp, simdata.NstepsinMemory) , tmpdata)
    #copytofulllattice!( tmpdata.F[ichunk], simdata.F[tmone,tp], tmpdata)
    tmpdata.rx[ichunk] = real( tmpdata.ftplan * tmpdata.r[ichunk])
    tmpdata.Fx[ichunk] = real( tmpdata.ftplan * tmpdata.F[ichunk])
    simdata.Sigrfull[tp] .= disc.ivol^3 * (-model.Lambda^2*(model.ONgroup+2)/(6 *model.ONgroup^2)) * real(tmpdata.iftplan * ( tmpdata.rx[ichunk] .* (  tmpdata.Fx[ichunk].^2 .- (1/12.) * tmpdata.rx[ichunk].^2)) ) # Normalisation: from exchange of fft/ifft and powers in selfenergy cal
    simdata.SigFfull[tp] .= disc.ivol^3 * (-model.Lambda^2*(model.ONgroup+2)/(18*model.ONgroup^2)) * real(tmpdata.iftplan * ( tmpdata.Fx[ichunk] .* (  tmpdata.Fx[ichunk].^2 .- 0.75    * tmpdata.rx[ichunk].^2)) ) # Normalisation: from exchange of fft/ifft and powers in selfenergy cal
    copytoreducedlattice!( simdata.Sigr[tp], simdata.Sigrfull[tp], disc.fftwhelper)
    copytoreducedlattice!( simdata.SigF[tp], simdata.SigFfull[tp], disc.fftwhelper)
    #copytoreducedlattice!( simdata.Sigr[tp], simdata.Sigrfull[tp], tmpdata)
    #copytoreducedlattice!( simdata.SigF[tp], simdata.SigFfull[tp], tmpdata)
    return
end

function calcSelfEnergyattp!(model::TwoPIScalarPhi4, pexp::TwoPIScalarLambdaNLO, disc::TwoPIScalarDiscretizationCont, simdata::TwoPIScalarSimDataCPUcont, tmpdata::TwoPIScalarTmpDataCPUcont, tmone::Int64, tp::Int64, ichunk::Int64)
    simdata.rx[tp] = real( tmpdata.ftplan * (tmpdata.isofactortox .* simdata.r[tmone,tp] * thesign(tmone, tp, simdata.NstepsinMemory)) )
    simdata.rx[tp] .*= tmpdata.normfactortox
    simdata.Fx[tp] = real( tmpdata.ftplan * (tmpdata.isofactortox .* simdata.F[tmone,tp] ) )
    simdata.Fx[tp] .*= tmpdata.normfactortox

    # use Irx, IFx as tmp storage for Sigrx,SigFx
    simdata.Irx[tp] .= simdata.rx[tp] .* (  simdata.Fx[tp].^2 .- (1/12.) * simdata.rx[tp].^2 ) 
    simdata.Irx[tp] .*= (-model.Lambda^2*(model.ONgroup+2)/(6 *model.ONgroup^2))
    simdata.IFx[tp] .= simdata.Fx[tp] .* (  simdata.Fx[tp].^2 .- 0.75   * simdata.rx[tp].^2 ) 
    simdata.IFx[tp] .*= (-model.Lambda^2*(model.ONgroup+2)/(18*model.ONgroup^2))

    simdata.Sigr[tp] = real( tmpdata.iftplan * (tmpdata.isofactortop .* simdata.Irx[tp]))
    simdata.Sigr[tp] .*= tmpdata.normfactortop
    simdata.SigF[tp] = real( tmpdata.iftplan * (tmpdata.isofactortop .* simdata.IFx[tp]))
    simdata.SigF[tp] .*= tmpdata.normfactortop
    return
end

#function calcSelfEnergyattp!(model::TwoPIScalarPhi4, pexp::TwoPIScalarLambdaNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUred, tmpdata::TwoPIScalarTmpDataCPUred, tmone::Int64, tp::Int64, ichunk::Int64)
#    copytolattice!( tmpdata.r[ichunk], simdata.r[tmone,tp], disc.fftwhelper)
#    copytolattice!( tmpdata.F[ichunk], simdata.F[tmone,tp], disc.fftwhelper)
#
#    mul!(tmpdata.Fx[ichunk], tmpdata.ftplan, tmpdata.F[ichunk] )
#    mul!(tmpdata.rx[ichunk], tmpdata.ftplan, tmpdata.r[ichunk] )
#
#    s = thesign(tmone, tp, simdata.NstepsinMemory)
#    @. tmpdata.Fx2[ichunk] = s * tmpdata.rx[ichunk] * (  tmpdata.Fx[ichunk]^2 - 1/12. * tmpdata.rx[ichunk]^2)
#    mul!(simdata.Sigr[tp], tmpdata.iftplan, tmpdata.Fx2[ichunk] )
#    factor1 = disc.ivol^3 * (-model.Lambda^2*(model.ONgroup+2)/(6 *model.ONgroup^2))
#    simdata.Sigr[tp] .*= factor1
#
#    @. tmpdata.Fx2[ichunk] = tmpdata.Fx[ichunk] * (  tmpdata.Fx[ichunk]^2 - 0.75   * tmpdata.rx[ichunk]^2 )
#    mul!(simdata.SigF[tp], tmpdata.iftplan, tmpdata.Fx2[ichunk] )
#    factor2 = disc.ivol^3 * (-model.Lambda^2*(model.ONgroup+2)/(18*model.ONgroup^2))
#    simdata.SigF[tp] .*= factor2
#
#    copytoreducedlattice!( simdata.Sigr[tp], simdata.Sigrfull[tp], disc.fftwhelper)
#    copytoreducedlattice!( simdata.SigF[tp], simdata.SigFfull[tp], disc.fftwhelper)
#    return
#end

function calcSelfEnergyattp!(model::TwoPIScalarPhi4, pexp::TwoPIScalarLambdaNLOquantum, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUred, tmpdata::TwoPIScalarTmpDataCPUred, tmone::Int64, tp::Int64, ichunk::Int64)
    copytolattice!( tmpdata.r[ichunk], simdata.r[tmone,tp] * thesign(tmone, tp, simdata.NstepsinMemory) , disc.fftwhelper)
    copytolattice!( tmpdata.F[ichunk], simdata.F[tmone,tp], disc.fftwhelper)
    #copytofulllattice!( tmpdata.r[ichunk], simdata.r[tmone,tp] * thesign(tmone, tp, simdata.NstepsinMemory) , tmpdata)
    #copytofulllattice!( tmpdata.F[ichunk], simdata.F[tmone,tp], tmpdata)
    tmpdata.rx[ichunk] = real( tmpdata.ftplan * tmpdata.r[ichunk])
    tmpdata.Fx[ichunk] = real( tmpdata.ftplan * tmpdata.F[ichunk])
    simdata.Sigrfull[tp] .= disc.ivol^3 * (-model.Lambda^2*(model.ONgroup+2)/(6 *model.ONgroup^2)) * real(tmpdata.iftplan * ( tmpdata.rx[ichunk] .* (  tmpdata.Fx[ichunk].^2 .- (1/12.) * tmpdata.rx[ichunk].^2)) ) # Normalisation: from exchange of fft/ifft and powers in selfenergy cal
    simdata.SigFfull[tp] .= disc.ivol^3 * (-model.Lambda^2*(model.ONgroup+2)/(18*model.ONgroup^2)) * real(tmpdata.iftplan * ( tmpdata.Fx[ichunk] .* (  tmpdata.Fx[ichunk].^2 .- 0.75    * tmpdata.rx[ichunk].^2)) ) # Normalisation: from exchange of fft/ifft and powers in selfenergy cal
    copytoreducedlattice!( simdata.Sigr[tp], simdata.Sigrfull[tp], disc.fftwhelper)
    copytoreducedlattice!( simdata.SigF[tp], simdata.SigFfull[tp], disc.fftwhelper)
    #copytoreducedlattice!( simdata.Sigr[tp], simdata.Sigrfull[tp], tmpdata)
    #copytoreducedlattice!( simdata.SigF[tp], simdata.SigFfull[tp], tmpdata)
    return
end

function calcSelfEnergyattp!(model::TwoPIScalarPhi4, pexp::TwoPIScalarLambdaNLOclassical, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUred, tmpdata::TwoPIScalarTmpDataCPUred, tmone::Int64, tp::Int64, ichunk::Int64)
    copytolattice!( tmpdata.r[ichunk], simdata.r[tmone,tp] * thesign(tmone, tp, simdata.NstepsinMemory) , disc.fftwhelper)
    copytolattice!( tmpdata.F[ichunk], simdata.F[tmone,tp], disc.fftwhelper)
    #copytofulllattice!( tmpdata.r[ichunk], simdata.r[tmone,tp] * thesign(tmone, tp, simdata.NstepsinMemory) , tmpdata)
    #copytofulllattice!( tmpdata.F[ichunk], simdata.F[tmone,tp], tmpdata)
    tmpdata.rx[ichunk] = real( tmpdata.ftplan * tmpdata.r[ichunk])
    tmpdata.Fx[ichunk] = real( tmpdata.ftplan * tmpdata.F[ichunk])
    simdata.Sigrfull[tp] .= disc.ivol^3 * (-model.Lambda^2*(model.ONgroup+2)/(6 *model.ONgroup^2)) * real(tmpdata.iftplan * ( tmpdata.rx[ichunk] .* (  tmpdata.Fx[ichunk].^2 )) ) # Normalisation: from exchange of fft/ifft and powers in selfenergy cal
    simdata.SigFfull[tp] .= disc.ivol^3 * (-model.Lambda^2*(model.ONgroup+2)/(18*model.ONgroup^2)) * real(tmpdata.iftplan * ( tmpdata.Fx[ichunk] .* (  tmpdata.Fx[ichunk].^2 )) ) # Normalisation: from exchange of fft/ifft and powers in selfenergy cal
    copytoreducedlattice!( simdata.Sigr[tp], simdata.Sigrfull[tp], disc.fftwhelper)
    copytoreducedlattice!( simdata.SigF[tp], simdata.SigFfull[tp], disc.fftwhelper)
    #copytoreducedlattice!( simdata.Sigr[tp], simdata.Sigrfull[tp], tmpdata)
    #copytoreducedlattice!( simdata.SigF[tp], simdata.SigFfull[tp], tmpdata)
    return
end

function calcSelfEnergiessimple!(model::TwoPIScalarPhi4, pexp::TwoPIScalarLambdaNLOquantum, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUfull, tmpdata::TwoPIScalarTmpDataCPUfull, tmone::Int64)
    for tp in simdata.indices[1]:(simdata.indices[2]-1) # [tmone,tp], the diagonal one is tmone,tmone
        calcSelfEnergyattpsimple!(model, pexp, disc, simdata, tmpdata, tmone, tp, 1)
    end
    return
end

function calcSelfEnergyattpsimple!(model::TwoPIScalarPhi4, pexp::TwoPIScalarLambdaNLOquantum, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUfull, tmpdata::TwoPIScalarTmpDataCPUfull, tmone::Int64, tp::Int64, ichunk::Int64)
    # brute force version, no prealloc
    rx = real( rfft(simdata.r[tmone,tp] * thesign(tmone, tp, simdata.NstepsinMemory)) )
    Fx = real( rfft(simdata.F[tmone,tp]) )
    Sigrx = rx .* (Fx.^2 .- (1/12.) * rx.^2) 
    SigFx = Fx .* (Fx.^2 .-    0.75 * rx.^2) 
    simdata.Sigr[tp] .= disc.ivol^3 * (-model.Lambda^2*(model.ONgroup+2)/(6 *model.ONgroup^2)) * real( brfft( Sigrx, disc.Nx) )
    simdata.SigF[tp] .= disc.ivol^3 * (-model.Lambda^2*(model.ONgroup+2)/(18*model.ONgroup^2)) * real( brfft( SigFx, disc.Nx) )
    return
end

###########################################################################################################################
# 1/N expansion 
###########################################################################################################################
function calcSelfEnergies!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPU, tmpdata::TwoPIScalarTmpDataCPU, tmone::Int64)
    calcIrs!(model, pexp, disc, simdata, tmpdata, tmone) # for all tp!
    calcIFs!(model, pexp, disc, simdata, tmpdata, tmone) # for all tp!
    @Threads.threads for ichunk in 1:tmpdata.nchunks
        for tp in tmpdata.threadranges[ichunk]
        #for tp in simdata.indices[1]:simdata.indices[2]-1
            calcSelfEnergyattp!(model, pexp, disc, simdata, tmpdata, tmone, tp, ichunk) # doesn t need tmpstorrange, therefore last index is not in use
        end
    end
    return
end

export fakecalcSelfEnergies!
function fakecalcSelfEnergies!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPU, tmpdata::TwoPIScalarTmpDataCPU, tmone::Int64)
    #calcIrs!(model, pexp, disc, simdata, tmpdata, tmone) # for all tp!
    #calcIFs!(model, pexp, disc, simdata, tmpdata, tmone) # for all tp!
    @Threads.threads for ichunk in 1:tmpdata.nchunks
        for tp in tmpdata.threadranges[ichunk]
        #for tp in simdata.indices[1]:simdata.indices[2]-1
            calcSelfEnergyattp!(model, pexp, disc, simdata, tmpdata, tmone, tp, ichunk) # doesn t need tmpstorrange, therefore last index is not in use
        end
    end
    return
end

export calcIrs!
function calcIrs!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUfull, tmpdata::TwoPIScalarTmpDataCPUfull, tmone::Int64)
    simdata.Irx[tmone] .= 0 #because of r[tmone,tmone]
    factor1 = (model.Lambda/3.) * disc.ivol^2
    factor2 = disc.dt * (model.Lambda/3.) * disc.ivol
    factor3 = disc.ivol^2
    for tp in (tmone-1):-1:simdata.indices[1] # counting down because of Ir dependence
        mul!(tmpdata.Fx[1], tmpdata.ftplan, simdata.F[tmone,tp] )
        mul!(tmpdata.rx[1], tmpdata.ftplan, simdata.r[tmone,tp] )
        s = thesign(tmone, tp, simdata.NstepsinMemory) 
        @. simdata.Irx[tp] = s * factor1 * real(tmpdata.Fx[1]) * real(tmpdata.rx[1])
        #simdata.Irx[tp] .= real( (model.Lambda/3.) * disc.ivol^2 * (tmpdata.ftplan * simdata.F[tmone,tp]) .* (tmpdata.ftplan * simdata.r[tmone,tp] * thesign(tmone, tp, simdata.NstepsinMemory)) )
        tmpdata.threadranges2 .= splitter(tp+1, (tmone-1) - (tp+1) + 1, tmpdata.nchunks) #smallest index, nr of indices
        @Threads.threads for ichunk in 1:tmpdata.nchunks
            tmpdata.RHS[ichunk] .= 0 # each thread writes to its tmpstorage
            for z in tmpdata.threadranges2[ichunk]
                mul!(tmpdata.Fx[ichunk], tmpdata.ftplan, simdata.F[z,tp] )
                mul!(tmpdata.rx[ichunk], tmpdata.ftplan, simdata.r[z,tp] )
                s = thesign(z, tp, simdata.NstepsinMemory)
                tmpdata.Fx2[ichunk] .=  real.(tmpdata.Fx[ichunk]) .* real.(tmpdata.rx[ichunk])
                mul!(tmpdata.RHS2[ichunk] , tmpdata.iftplan, tmpdata.Fx2[ichunk] )
                #that line causes problems in sdim 2 & 3:
                #mul!(tmpdata.RHS3[ichunk] , tmpdata.iftplan, simdata.Irx[z] )
                #- therefore I use the following two lines:
                tmpdata.Fx2[ichunk] .=  simdata.Irx[z] 
                mul!(tmpdata.RHS3[ichunk] , tmpdata.iftplan, tmpdata.Fx2[ichunk] )
                tmpdata.RHS[ichunk] .-= factor3 * s .* real.(tmpdata.RHS2[ichunk]) .* real.(tmpdata.RHS3[ichunk] )
                #tmpdata.RHS[ichunk] .-= real( (tmpdata.iftplan * simdata.Irx[z]) .* (tmpdata.iftplan * ( disc.ivol^2 * (tmpdata.ftplan * simdata.F[z,tp]) .* (tmpdata.ftplan * simdata.r[z,tp] * thesign(z, tp, simdata.NstepsinMemory))))  )
            end
        end
        for ichunk in 2:tmpdata.nchunks
            tmpdata.RHS[1] .+=  tmpdata.RHS[ichunk] 
        end
        mul!(tmpdata.Fx[1], tmpdata.ftplan, tmpdata.RHS[1] )
        @. simdata.Irx[tp] += factor2 * real(tmpdata.Fx[1])
        #simdata.Irx[tp] .+= disc.dt * (model.Lambda/3.) * disc.ivol * real( tmpdata.ftplan * tmpdata.RHS[1] ) 
    end
end

function calcIrs!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUred, tmpdata::TwoPIScalarTmpDataCPUred, tmone::Int64)
    simdata.Irx[tmone] .= 0 #because of r[tmone,tmone]
    factor1 = (model.Lambda/3.) * disc.ivol^2
    factor2 = disc.dt * (model.Lambda/3.) * disc.ivol
    factor3 = disc.ivol^2
    for tp in (tmone-1):-1:simdata.indices[1] # counting down because of Ir dependence
        copytolattice!( tmpdata.r[1], simdata.r[tmone,tp], disc.fftwhelper)
        copytolattice!( tmpdata.F[1], simdata.F[tmone,tp], disc.fftwhelper)
        #copytofulllattice!( tmpdata.r[1], simdata.r[tmone,tp], tmpdata)
        #copytofulllattice!( tmpdata.F[1], simdata.F[tmone,tp], tmpdata)
        mul!(tmpdata.Fx[1], tmpdata.ftplan, tmpdata.F[1] )
        mul!(tmpdata.rx[1], tmpdata.ftplan, tmpdata.r[1] )
        s = thesign(tmone, tp, simdata.NstepsinMemory) 
        @. simdata.Irx[tp] = s * factor1 * real(tmpdata.Fx[1]) * real(tmpdata.rx[1])
        tmpdata.threadranges2 .= splitter(tp+1, (tmone-1) - (tp+1) + 1, tmpdata.nchunks) #smallest index, nr of indices
        @Threads.threads for ichunk in 1:tmpdata.nchunks
            tmpdata.tmpfull[ichunk] .= 0
            for z in tmpdata.threadranges2[ichunk]
                copytolattice!( tmpdata.r[ichunk], simdata.r[z,tp], disc.fftwhelper)
                copytolattice!( tmpdata.F[ichunk], simdata.F[z,tp], disc.fftwhelper)
                #copytofulllattice!( tmpdata.r[ichunk], simdata.r[z,tp], tmpdata)
                #copytofulllattice!( tmpdata.F[ichunk], simdata.F[z,tp], tmpdata)
                mul!(tmpdata.Fx[ichunk], tmpdata.ftplan, tmpdata.F[ichunk] )
                mul!(tmpdata.rx[ichunk], tmpdata.ftplan, tmpdata.r[ichunk] )
                s = thesign(z, tp, simdata.NstepsinMemory)
                tmpdata.Fx[ichunk] .=  tmpdata.Fx[ichunk] .* tmpdata.rx[ichunk]
                mul!(tmpdata.tmpfull2[ichunk] , tmpdata.iftplan, tmpdata.Fx[ichunk] )
                #that line causes problems in sdim 2 & 3:
                #mul!(tmpdata.tmpfull3[ichunk] , tmpdata.iftplan, simdata.Irx[z] )
                #- therefore I use the following two lines:
                tmpdata.Fx2[ichunk] .=  simdata.Irx[z] 
                mul!(tmpdata.tmpfull3[ichunk], tmpdata.iftplan, tmpdata.Fx2[ichunk] )
                tmpdata.tmpfull[ichunk] .-= real.( factor3 * s .* tmpdata.tmpfull2[ichunk] .* tmpdata.tmpfull3[ichunk] )
                #tmpdata.tmpfull[ichunk] .-= real( (tmpdata.iftplan * simdata.Irx[z]) .* (tmpdata.iftplan * ( disc.ivol^2 * (tmpdata.ftplan * tmpdata.F[ichunk]) .* (tmpdata.ftplan * tmpdata.r[ichunk])))  )
            end
        end
        for ichunk in 2:tmpdata.nchunks
            tmpdata.tmpfull[1] .+=  tmpdata.tmpfull[ichunk] 
        end
        mul!(tmpdata.Fx[1], tmpdata.ftplan, tmpdata.tmpfull[1] )
        @. simdata.Irx[tp] += factor2 * tmpdata.Fx[1] 
    end
end

function calcIrs!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUred2, tmpdata::TwoPIScalarTmpDataCPUred, tmone::Int64)
    simdata.Irk[tmone] .= 0 #because of r[tmone,tmone]
    factor1 = (model.Lambda/3.)
    for tp in (tmone-1):-1:simdata.indices[1] # counting down because of Ir dependence
        # resum kernel part
        tmpdata.threadranges2 .= splitter(tp+1, (tmone-1) - (tp+1) + 1, tmpdata.nchunks) #smallest index, nr of indices
        @Threads.threads for ichunk in 1:tmpdata.nchunks
            tmpdata.RHS[ichunk] .= 0
            for z in tmpdata.threadranges2[ichunk]
                s = thesign(z, tp, simdata.NstepsinMemory) 
                tmpdata.RHS[ichunk] .-= real.( s.* simdata.Irk[z] .* simdata.Fr[z,tp] )
            end
        end
        # collect thread contribution
        simdata.Irk[tp] .= 0 
        for ichunk in 1:tmpdata.nchunks
            simdata.Irk[tp] .+=  tmpdata.RHS[ichunk] 
        end
        # multiply dt from time integral
        simdata.Irk[tp] .*= disc.dt
        # add local contribution 
        simdata.Irk[tp] .+= simdata.Fr[tmone,tp] * thesign(tmone, tp, simdata.NstepsinMemory) 
        # multiply factor
        simdata.Irk[tp] .*= factor1 
        # calc FT
#        copytolattice!(tmpdata.tmpfull[1], simdata.Irk[tp], disc.fftwhelper)
#        mul!(simdata.Irx[tp], tmpdata.ftplan, tmpdata.tmpfull[1] )
#        simdata.Irx[tp] .*= disc.ivol
    end
    # calc FT of Irk in parallel
    simdata.Irx[tmone] .= 0 #because of r[tmone,tmone]
    #for tp in simdata.indices[1]:(tmone-1)
    tmpdata.threadranges2 .= splitter(simdata.indices[1], (tmone-1) - simdata.indices[1] + 1, tmpdata.nchunks) #smallest index, nr of indices
    @Threads.threads for ichunk in 1:tmpdata.nchunks
        for tp in tmpdata.threadranges2[ichunk]
            copytolattice!(tmpdata.tmpfull[ichunk], simdata.Irk[tp], disc.fftwhelper)
            mul!(simdata.Irx[tp], tmpdata.ftplan, tmpdata.tmpfull[ichunk] )
            simdata.Irx[tp] .*= disc.ivol
        end
    end
end

function calcIrssimple!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUfull, tmpdata::TwoPIScalarTmpDataCPUfull, tmone::Int64)
    #cnt = zeros(tmone-simdata.indices[1]+1)
    simdata.Irx[tmone] .= 0 #because of r[tmone,tmone]
    for tp in (tmone-1):-1:simdata.indices[1] # counting down because of Ir dependence
        # add local contribution: Irx[tmone] will be 0 because of r[tmone,tmone]
        #simdata.Irx[tp] .= (model.Lambda/3.) * bfft(simdata.F[tmone,tp]) .* bfft(simdata.r[tmone,tp] * thesign(tmone, tp, simdata.NstepsinMemory))
        simdata.Irx[tp] .= real( (model.Lambda/3.) * disc.ivol^2 * rfft(simdata.F[tmone,tp]) .* rfft(simdata.r[tmone,tp] * thesign(tmone, tp, simdata.NstepsinMemory)) )
        #cnt[tp] = 1*tp #local -> one point
        # scan through the past
        for z in tp+1:tmone-1 # first contribution at tp=tmone-2
            #println("z:", z)
            #Bx = bfft(simdata.F[z,tp]) .* bfft(simdata.r[z,tp] * thesign(z, tp, simdata.NstepsinMemory))
            Bx = real( disc.ivol^2 * rfft(simdata.F[z,tp]) .* rfft(simdata.r[z,tp] * thesign(z, tp, simdata.NstepsinMemory)) )
            #simdata.Irx[tp] .-= disc.dt * (model.Lambda/3.) * bfft( fft(simdata.Irx[z]) .* fft(Bx) ) 
            simdata.Irx[tp] .-= real( disc.dt * (model.Lambda/3.) * disc.ivol * rfft( brfft(simdata.Irx[z],disc.Nx) .* brfft(Bx,disc.Nx) ) )
            #cnt[tp] = 5*z
        end
        #simdata.Irx[tp] .*= disc.ivol^2
    end
    #if tmone==37 @show cnt end
end

#function calcIrs!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUcont, tmpdata::TwoPIScalarTmpDataCPUcont, tmone::Int64)
#    # TODO
#    simdata.Irx[tmone] .= 0 #because of r[tmone,tmone]
#    simdata.Irk[tmone] .= 0 #because of r[tmone,tmone]
#    factor1 = (model.Lambda/3.)
#    for tp in (tmone-1):-1:simdata.indices[1] # counting down because of Ir dependence
#        # resum kernel part
#        tmpdata.threadranges2 .= splitter(tp+1, (tmone-1) - (tp+1) + 1, tmpdata.nchunks) #smallest index, nr of indices
#        @Threads.threads for ichunk in 1:tmpdata.nchunks
#            tmpdata.RHS[ichunk] .= 0
#            for z in tmpdata.threadranges2[ichunk]
#                s = thesign(z, tp, simdata.NstepsinMemory) 
#                tmpdata.RHS[ichunk] .-= real.( s.* simdata.Irk[z] .* simdata.Fr[z,tp] )
#            end
#        end
#        # collect thread contribution
#        simdata.Irk[tp] .= 0 
#        for ichunk in 1:tmpdata.nchunks
#            simdata.Irk[tp] .+=  tmpdata.RHS[ichunk] 
#        end
#        # multiply dt from time integral
#        simdata.Irk[tp] .*= disc.dt
#        # add local contribution 
#        simdata.Irk[tp] .+= simdata.Fr[tmone,tp] * thesign(tmone, tp, simdata.NstepsinMemory) 
#        # multiply factor
#        simdata.Irk[tp] .*= factor1 
#        # calc FT
#        #copytolattice!(tmpdata.tmpfull[1], simdata.Irk[tp], disc.fftwhelper)
#        #mul!(simdata.Irx[tp], tmpdata.ftplan, tmpdata.tmpfull[1] )
#        #simdata.Irx[tp] .*= disc.ivol
#        tmpdata.RHS[1] .= simdata.Irk[tp]
#        tmpdata.RHS[1] .*= tmpdata.isofactortox
#        mul!(simdata.Irx[tp], tmpdata.ftplan, tmpdata.RHS[1] )
#        simdata.Irx[tp] .*= tmpdata.normfactortox
#    end
#end

export calcIFs!
function calcIFs!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUfull, tmpdata::TwoPIScalarTmpDataCPUfull, tmone::Int64)
    factor1 = model.Lambda/6. * disc.ivol^2
    factor2 = model.Lambda/6. * disc.dt
    factor3 = model.Lambda/3. * disc.dt
    for tp in simdata.indices[1]:tmone
        mul!(tmpdata.Fx[1], tmpdata.ftplan, simdata.F[tmone,tp] )
        mul!(tmpdata.rx[1], tmpdata.ftplan, simdata.r[tmone,tp] )
        @. tmpdata.Fx[1] =  tmpdata.Fx[1]^2 - 0.25 * tmpdata.rx[1]^2
        @. simdata.IFx[tp] = real( factor1 * tmpdata.Fx[1])
        # first integral
        tmpdata.threadranges2 .= splitter(simdata.indices[1]+1, (tp-1) - (simdata.indices[1]+1) + 1, tmpdata.nchunks) #smallest index, nr of indices
        @Threads.threads for ichunk in 1:tmpdata.nchunks
            tmpdata.RHS[ichunk] .= 0
            for z in tmpdata.threadranges2[ichunk]
                mul!(tmpdata.Fx[ichunk], tmpdata.ftplan, simdata.F[z,tp] )
                mul!(tmpdata.rx[ichunk], tmpdata.ftplan, simdata.r[z,tp] )
                tmpdata.Fx[ichunk] .*= disc.ivol
                tmpdata.rx[ichunk] .*= disc.ivol
                # this line
                #mul!(tmpdata.RHS3[ichunk] , tmpdata.iftplan, simdata.Irx[z] )
                tmpdata.Fx2[ichunk] .= simdata.Irx[z]
                mul!(tmpdata.RHS3[ichunk] , tmpdata.iftplan, tmpdata.Fx2[ichunk])
                @. tmpdata.Fx2[ichunk] = tmpdata.Fx[ichunk]^2 - 0.25 * tmpdata.rx[ichunk]^2 
                mul!(tmpdata.RHS2[ichunk] , tmpdata.iftplan, tmpdata.Fx2[ichunk] )
                @. tmpdata.RHS[ichunk] -= real( factor2 * tmpdata.RHS3[ichunk] * tmpdata.RHS2[ichunk] ) #nr 1
                @. tmpdata.Fx2[ichunk] = tmpdata.Fx[ichunk] * tmpdata.rx[ichunk]
                s = thesign(z, tp, simdata.NstepsinMemory) 
                mul!(tmpdata.RHS3[ichunk] , tmpdata.iftplan, tmpdata.Fx2[ichunk] )
                # this line
                #mul!(tmpdata.RHS2[ichunk] , tmpdata.iftplan, simdata.IFx[z] ) 
                tmpdata.Fx2[ichunk] .= simdata.IFx[z]
                mul!(tmpdata.RHS2[ichunk] , tmpdata.iftplan, tmpdata.Fx2[ichunk]) 
                @. tmpdata.RHS[ichunk] += real( factor3 * s * tmpdata.RHS2[ichunk] * tmpdata.RHS3[ichunk] ) #nr 2
            end
        end
        # second integral
        if tp == simdata.indices[1] beginning = simdata.indices[1]+1 else beginning = tp end
        tmpdata.threadranges3 .= splitter(beginning, (tmone-1) - beginning + 1, tmpdata.nchunks) #smallest index, nr of indices
        @Threads.threads for ichunk in 1:tmpdata.nchunks
            for z in tmpdata.threadranges3[ichunk]
                mul!(tmpdata.Fx[ichunk], tmpdata.ftplan, simdata.F[z,tp] )
                mul!(tmpdata.rx[ichunk], tmpdata.ftplan, simdata.r[z,tp] )
                tmpdata.Fx[ichunk] .*= disc.ivol
                tmpdata.rx[ichunk] .*= disc.ivol
                @. tmpdata.Fx2[ichunk] = tmpdata.Fx[ichunk]^2 - 0.25 * tmpdata.rx[ichunk]^2 
                mul!(tmpdata.RHS3[ichunk] , tmpdata.iftplan, tmpdata.Fx2[ichunk] ) 
                # this line
                #mul!(tmpdata.RHS2[ichunk] , tmpdata.iftplan, simdata.Irx[z] ) 
                tmpdata.Fx2[ichunk] .= simdata.Irx[z]
                mul!(tmpdata.RHS2[ichunk] , tmpdata.iftplan, tmpdata.Fx2[ichunk] ) 
                @. tmpdata.RHS[ichunk] -= real( factor2 * tmpdata.RHS2[ichunk] * tmpdata.RHS3[ichunk] )
            end
        end
        for ichunk in 2:tmpdata.nchunks
            tmpdata.RHS[1] .+=  tmpdata.RHS[ichunk] 
        end
        # border part
        mul!(tmpdata.Fx[1], tmpdata.ftplan, simdata.F[simdata.indices[1],tp] )
        mul!(tmpdata.rx[1], tmpdata.ftplan, simdata.r[simdata.indices[1],tp] )
        @. tmpdata.Fx2[1] = tmpdata.Fx[1]^2 - 0.25 * tmpdata.rx[1]^2 
        mul!(tmpdata.RHS3[1] , tmpdata.iftplan, tmpdata.Fx2[1] ) 
        tmpdata.RHS3[1]  .*= disc.ivol^2
        # this line
        #mul!(tmpdata.RHS2[1] , tmpdata.iftplan, simdata.Irx[simdata.indices[1]] ) 
        tmpdata.Fx2[1] .= simdata.Irx[simdata.indices[1]]
        mul!(tmpdata.RHS2[1] , tmpdata.iftplan, tmpdata.Fx2[1] ) 
        @. tmpdata.RHS[1] -= real( 0.5*disc.dt*model.Lambda/6.  *  tmpdata.RHS2[1] * tmpdata.RHS3[1] )
        if tp>simdata.indices[1]
            mul!(tmpdata.Fx[1], tmpdata.ftplan, simdata.F[simdata.indices[1],tp] )
            mul!(tmpdata.rx[1], tmpdata.ftplan, simdata.r[simdata.indices[1],tp] )
            s = thesign(simdata.indices[1], tp, simdata.NstepsinMemory)
            @. tmpdata.Fx2[1] = s * tmpdata.Fx[1] * tmpdata.rx[1]
            mul!(tmpdata.RHS3[1], tmpdata.iftplan, tmpdata.Fx2[1] )
            # this line
            #mul!(tmpdata.RHS2[1] , tmpdata.iftplan, simdata.IFx[simdata.indices[1]] ) 
            tmpdata.Fx2[1] .= simdata.IFx[simdata.indices[1]]
            mul!(tmpdata.RHS2[1] , tmpdata.iftplan, tmpdata.Fx2[1] ) 
            @. tmpdata.RHS[1] += real( 0.5*disc.dt*model.Lambda/3. * disc.ivol^2 * tmpdata.RHS2[1] * tmpdata.RHS3[1] )
        end
        mul!(tmpdata.Fx[1], tmpdata.ftplan, tmpdata.RHS[1] )
        @. simdata.IFx[tp] += disc.ivol * tmpdata.Fx[1] 
    end
end

function calcIFs!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUred, tmpdata::TwoPIScalarTmpDataCPUred, tmone::Int64)
    factor1 = model.Lambda/6. * disc.ivol^2
    factor2 = model.Lambda/6. * disc.dt
    factor3 = model.Lambda/3. * disc.dt
    for tp in simdata.indices[1]:tmone
        copytolattice!( tmpdata.r[1], simdata.r[tmone,tp], disc.fftwhelper)
        copytolattice!( tmpdata.F[1], simdata.F[tmone,tp], disc.fftwhelper)
        #copytofulllattice!( tmpdata.r[1], simdata.r[tmone,tp], tmpdata)
        #copytofulllattice!( tmpdata.F[1], simdata.F[tmone,tp], tmpdata)
        mul!(tmpdata.Fx[1], tmpdata.ftplan, tmpdata.F[1] )
        mul!(tmpdata.rx[1], tmpdata.ftplan, tmpdata.r[1] )
        @. tmpdata.Fx[1] =  tmpdata.Fx[1]^2 - 0.25 * tmpdata.rx[1]^2
        @. simdata.IFx[tp] = real( factor1 * tmpdata.Fx[1])
        # first integral
        tmpdata.threadranges2 .= splitter(simdata.indices[1]+1, (tp-1) - (simdata.indices[1]+1) + 1, tmpdata.nchunks) #smallest index, nr of indices
        @Threads.threads for ichunk in 1:tmpdata.nchunks
            tmpdata.tmpfull[ichunk] .=0
            for z in tmpdata.threadranges2[ichunk]
                copytolattice!( tmpdata.r[ichunk], simdata.r[z,tp], disc.fftwhelper)
                copytolattice!( tmpdata.F[ichunk], simdata.F[z,tp], disc.fftwhelper)
                #copytofulllattice!( tmpdata.r[ichunk], simdata.r[z,tp], tmpdata)
                #copytofulllattice!( tmpdata.F[ichunk], simdata.F[z,tp], tmpdata)
                mul!(tmpdata.Fx[ichunk], tmpdata.ftplan, tmpdata.F[ichunk] )
                mul!(tmpdata.rx[ichunk], tmpdata.ftplan, tmpdata.r[ichunk] )
                tmpdata.Fx[ichunk] .*= disc.ivol
                tmpdata.rx[ichunk] .*= disc.ivol
                # this line
                #mul!(tmpdata.tmpfull3[ichunk] , tmpdata.iftplan, simdata.Irx[z] )
                tmpdata.Fx2[ichunk] .= simdata.Irx[z]
                mul!(tmpdata.tmpfull3[ichunk] , tmpdata.iftplan, tmpdata.Fx2[ichunk] )
                @. tmpdata.Fx2[ichunk] = tmpdata.Fx[ichunk]^2 - 0.25 * tmpdata.rx[ichunk]^2 
                mul!(tmpdata.tmpfull2[ichunk] , tmpdata.iftplan, tmpdata.Fx2[ichunk] )
                @. tmpdata.tmpfull[ichunk] -= real( factor2 * tmpdata.tmpfull3[ichunk] * tmpdata.tmpfull2[ichunk] ) #nr 1
                @. tmpdata.Fx2[ichunk] = tmpdata.Fx[ichunk] * tmpdata.rx[ichunk]
                s = thesign(z, tp, simdata.NstepsinMemory) 
                mul!(tmpdata.tmpfull3[ichunk] , tmpdata.iftplan, tmpdata.Fx2[ichunk] )
                # this line
                #mul!(tmpdata.tmpfull2[ichunk] , tmpdata.iftplan, simdata.IFx[z] ) 
                tmpdata.Fx2[ichunk] .= simdata.IFx[z]
                mul!(tmpdata.tmpfull2[ichunk] , tmpdata.iftplan, tmpdata.Fx2[ichunk] ) 
                @. tmpdata.tmpfull[ichunk] += real( factor3 * s * tmpdata.tmpfull2[ichunk] * tmpdata.tmpfull3[ichunk] ) #nr 2
                ## just for printing
                #copytoreducedlattice!( tmpdata.RHS2[ichunk], factor3 * s * tmpdata.tmpfull2[ichunk] .* tmpdata.tmpfull3[ichunk], disc.fftwhelper)
                #@show real.( tmpdata.RHS2[ichunk] ) #nr 2
            end
        end
        # second integral
        if tp == simdata.indices[1] beginning = simdata.indices[1]+1 else beginning = tp end
        tmpdata.threadranges3 .= splitter(beginning, (tmone-1) - beginning + 1, tmpdata.nchunks) #smallest index, nr of indices
        @Threads.threads for ichunk in 1:tmpdata.nchunks
            for z in tmpdata.threadranges3[ichunk]
                copytolattice!( tmpdata.r[ichunk], simdata.r[z,tp], disc.fftwhelper)
                copytolattice!( tmpdata.F[ichunk], simdata.F[z,tp], disc.fftwhelper)
                #copytofulllattice!( tmpdata.r[ichunk], simdata.r[z,tp], tmpdata)
                #copytofulllattice!( tmpdata.F[ichunk], simdata.F[z,tp], tmpdata)
                mul!(tmpdata.Fx[ichunk], tmpdata.ftplan, tmpdata.F[ichunk] )
                mul!(tmpdata.rx[ichunk], tmpdata.ftplan, tmpdata.r[ichunk] )
                tmpdata.Fx[ichunk] .*= disc.ivol
                tmpdata.rx[ichunk] .*= disc.ivol
                @. tmpdata.Fx2[ichunk] = tmpdata.Fx[ichunk]^2 - 0.25 * tmpdata.rx[ichunk]^2 
                mul!(tmpdata.tmpfull3[ichunk] , tmpdata.iftplan, tmpdata.Fx2[ichunk] ) 
                # this line
                #mul!(tmpdata.tmpfull2[ichunk] , tmpdata.iftplan, simdata.Irx[z] ) 
                tmpdata.Fx2[ichunk] .= simdata.Irx[z]
                mul!(tmpdata.tmpfull2[ichunk] , tmpdata.iftplan, tmpdata.Fx2[ichunk] ) 
                @. tmpdata.tmpfull[ichunk] -= real( factor2 * tmpdata.tmpfull2[ichunk] * tmpdata.tmpfull3[ichunk] )
            end
        end
        for ichunk in 2:tmpdata.nchunks
            tmpdata.tmpfull[1] .+=  tmpdata.tmpfull[ichunk] 
        end
        # border part
        copytolattice!( tmpdata.r[1], simdata.r[simdata.indices[1],tp], disc.fftwhelper)
        copytolattice!( tmpdata.F[1], simdata.F[simdata.indices[1],tp], disc.fftwhelper)
        #copytofulllattice!( tmpdata.r[1], simdata.r[simdata.indices[1],tp], tmpdata)
        #copytofulllattice!( tmpdata.F[1], simdata.F[simdata.indices[1],tp], tmpdata)
        mul!(tmpdata.Fx[1], tmpdata.ftplan, tmpdata.F[1] )
        mul!(tmpdata.rx[1], tmpdata.ftplan, tmpdata.r[1] )
        @. tmpdata.Fx2[1] = tmpdata.Fx[1]^2 - 0.25 * tmpdata.rx[1]^2 
        mul!(tmpdata.tmpfull3[1] , tmpdata.iftplan, tmpdata.Fx2[1] ) 
        tmpdata.tmpfull3[1]  .*= disc.ivol^2
        # this line
        #mul!(tmpdata.tmpfull2[1] , tmpdata.iftplan, simdata.Irx[simdata.indices[1]] ) 
        tmpdata.Fx2[1] .= simdata.Irx[simdata.indices[1]]
        mul!(tmpdata.tmpfull2[1] , tmpdata.iftplan, tmpdata.Fx2[1] ) 
        @. tmpdata.tmpfull[1] -= real( 0.5*disc.dt*model.Lambda/6.  *  tmpdata.tmpfull2[1] * tmpdata.tmpfull3[1] )
        if tp>simdata.indices[1]
            s = thesign(simdata.indices[1], tp, simdata.NstepsinMemory)
            @. tmpdata.Fx2[1] = s * tmpdata.Fx[1] * tmpdata.rx[1]
            mul!(tmpdata.tmpfull3[1], tmpdata.iftplan, tmpdata.Fx2[1] )
            # this line
            #mul!(tmpdata.tmpfull2[1] , tmpdata.iftplan, simdata.IFx[simdata.indices[1]] ) 
            tmpdata.Fx2[1] .= simdata.IFx[simdata.indices[1]]
            mul!(tmpdata.tmpfull2[1] , tmpdata.iftplan, tmpdata.Fx2[1] ) 
            @. tmpdata.tmpfull[1] += real( 0.5*disc.dt*model.Lambda/3. * disc.ivol^2 * tmpdata.tmpfull2[1] * tmpdata.tmpfull3[1] )
        end
        mul!(tmpdata.Fx[1], tmpdata.ftplan, tmpdata.tmpfull[1] )
        @. simdata.IFx[tp] += disc.ivol * tmpdata.Fx[1] 
    end
end

function calcIFs!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUred2, tmpdata::TwoPIScalarTmpDataCPUred, tmone::Int64)
    factor1 = model.Lambda/6.
    for tp in simdata.indices[1]:tmone
        # resum kernel
        # first integral
        tmpdata.threadranges2 .= splitter(simdata.indices[1]+1, (tp-1) - (simdata.indices[1]+1) + 1, tmpdata.nchunks) #smallest index, nr of indices
        @Threads.threads for ichunk in 1:tmpdata.nchunks
            tmpdata.RHS[ichunk] .=0
            for z in tmpdata.threadranges2[ichunk]
                # nr 1
                @. tmpdata.RHS[ichunk] += real( simdata.Irk[z] * simdata.F2kr2[z,tp]  )
                # nr 2
                s = 2 * thesign(z, tp, simdata.NstepsinMemory) 
                @. tmpdata.RHS[ichunk] -= real( s * simdata.IFk[z] * simdata.Fr[z,tp] )
            end
        end
        # second integral
        if tp == simdata.indices[1] beginning = simdata.indices[1]+1 else beginning = tp end
        tmpdata.threadranges3 .= splitter(beginning, (tmone-1) - beginning + 1, tmpdata.nchunks) #smallest index, nr of indices
        @Threads.threads for ichunk in 1:tmpdata.nchunks
            for z in tmpdata.threadranges3[ichunk]
                # only nr 1
                @. tmpdata.RHS[ichunk] += real( simdata.Irk[z] * simdata.F2kr2[z,tp] )
            end
        end
        # sum contributions from threads
        simdata.IFk[tp] .= 0
        for ichunk in 1:tmpdata.nchunks
            simdata.IFk[tp] .-=  tmpdata.RHS[ichunk] 
        end
        # border part of nr 1
        @. simdata.IFk[tp] -= real( 0.5 * simdata.Irk[simdata.indices[1]] * simdata.F2kr2[simdata.indices[1],tp])
        # border part of nr 2
        if tp>simdata.indices[1]
            # Note: Factor 2 from nr 2 and factor 0.5 of border integration cnacel out
            s = thesign(simdata.indices[1], tp, simdata.NstepsinMemory) 
            @. simdata.IFk[tp] += real( s * simdata.IFk[simdata.indices[1]] * simdata.Fr[simdata.indices[1], tp] )
        end
        # multiply dt factor
        simdata.IFk[tp] .*= disc.dt
        # add local contribution
        simdata.IFk[tp] .+= simdata.F2kr2[tmone,tp]
        # multiply expansion parameter
        simdata.IFk[tp] .*= factor1
        # make fft
#        copytolattice!( tmpdata.tmpfull[1], simdata.IFk[tp], disc.fftwhelper)
#        mul!(simdata.IFx[tp], tmpdata.ftplan, tmpdata.tmpfull[1] )
#        simdata.IFx[tp] .*= disc.ivol
    end
    #for tp in simdata.indices[1]:tmone
    tmpdata.threadranges2 .= splitter(simdata.indices[1], (tmone) - simdata.indices[1] + 1, tmpdata.nchunks) #smallest index, nr of indices
    @Threads.threads for ichunk in 1:tmpdata.nchunks
        for tp in tmpdata.threadranges2[ichunk]
            copytolattice!( tmpdata.tmpfull[ichunk], simdata.IFk[tp], disc.fftwhelper)
            mul!(simdata.IFx[tp], tmpdata.ftplan, tmpdata.tmpfull[ichunk] )
            simdata.IFx[tp] .*= disc.ivol
        end
    end
end

function calcIFssimple!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUfull, tmpdata::TwoPIScalarTmpDataCPUfull, tmone::Int64)
    #cnt = zeros(tmone-simdata.indices[1]+1)
    for tp in simdata.indices[1]:tmone
        # local contribution
        #simdata.IFx[tp] .= (model.Lambda/6.) * ( bfft(simdata.F[tmone,tp]).^2 .- 0.25 * bfft(simdata.r[tmone,tp]).^2 )
        simdata.IFx[tp] .= real( (model.Lambda/6.) * disc.ivol^2 * ( rfft(simdata.F[tmone,tp]).^2 .- 0.25 * rfft(simdata.r[tmone,tp]).^2 ) )
        #cnt[tp] = 0.001
        # scan through the past
        # add Ir part
        for z in simdata.indices[1]:tmone-1
            if z == simdata.indices[1] w = 0.5*disc.dt else w = 1.0*disc.dt end
            #simdata.IFx[tp] .-= w*(model.Lambda/6.) * bfft( fft(simdata.Irx[z]) .* fft(bfft(simdata.F[z,tp]).^2 .- 0.25 * bfft(simdata.r[z,tp]).^2) )
            simdata.IFx[tp] .-= w*(model.Lambda/6.) * real( disc.ivol * rfft( brfft(simdata.Irx[z],disc.Nx) .* brfft( disc.ivol^2*  ( rfft(simdata.F[z,tp]).^2 .- 0.25 * rfft(simdata.r[z,tp]).^2) ,disc.Nx) ) )
            #if z == simdata.indices[1] cnt[tp] += 500 else cnt[tp] += 1000 end
        end
        # add IF part
        for z in simdata.indices[1]:tp-1 # count up because of dependency
            if z == simdata.indices[1] w = 0.5*disc.dt else w = 1.0*disc.dt end
            #simdata.IFx[tp] .+= w*(model.Lambda/3.) * bfft( fft(simdata.IFx[z]) .* fft(bfft(simdata.F[z,tp]) .* bfft(simdata.r[z,tp] * thesign(z, tp, simdata.NstepsinMemory))) ) 
            simdata.IFx[tp] .+= w*(model.Lambda/3.) * real( disc.ivol * rfft( brfft(simdata.IFx[z],disc.Nx) .* brfft( disc.ivol^2*  ( rfft(simdata.F[z,tp]) .* rfft(simdata.r[z,tp] * thesign(z, tp, simdata.NstepsinMemory)) ) ,disc.Nx) )  )
            #if z == simdata.indices[1] cnt[tp] += 0.5 else cnt[tp] += 1.0 end
        end
        #simdata.IFx[tp] .*= disc.ivol^2
    end
    #if tmone==3 println("simple:", cnt) end
end

#function calcIFs!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUcont, tmpdata::TwoPIScalarTmpDataCPUcont, tmone::Int64)
#    # TODO
#    factor1 = model.Lambda/6.
#    for tp in simdata.indices[1]:tmone
#        # resum kernel
#        # first integral
#        tmpdata.threadranges2 .= splitter(simdata.indices[1]+1, (tp-1) - (simdata.indices[1]+1) + 1, tmpdata.nchunks) #smallest index, nr of indices
#        @Threads.threads for ichunk in 1:tmpdata.nchunks
#            tmpdata.RHS[ichunk] .=0
#            for z in tmpdata.threadranges2[ichunk]
#                # nr 1
#                @. tmpdata.RHS[ichunk] += real( simdata.Irk[z] * simdata.F2kr2[z,tp]  )
#                # nr 2
#                s = 2 * thesign(z, tp, simdata.NstepsinMemory) 
#                @. tmpdata.RHS[ichunk] -= real( s * simdata.IFk[z] * simdata.Fr[z,tp] )
#            end
#        end
#        # second integral
#        if tp == simdata.indices[1] beginning = simdata.indices[1]+1 else beginning = tp end
#        tmpdata.threadranges3 .= splitter(beginning, (tmone-1) - beginning + 1, tmpdata.nchunks) #smallest index, nr of indices
#        @Threads.threads for ichunk in 1:tmpdata.nchunks
#            for z in tmpdata.threadranges3[ichunk]
#                # only nr 1
#                @. tmpdata.RHS[ichunk] += real( simdata.Irk[z] * simdata.F2kr2[z,tp] )
#            end
#        end
#        # sum contributions from threads
#        simdata.IFk[tp] .= 0
#        for ichunk in 1:tmpdata.nchunks
#            simdata.IFk[tp] .-=  tmpdata.RHS[ichunk] 
#        end
#        # border part of nr 1
#        @. simdata.IFk[tp] -= real( 0.5 * simdata.Irk[simdata.indices[1]] * simdata.F2kr2[simdata.indices[1],tp])
#        # border part of nr 2
#        if tp>simdata.indices[1]
#            # Note: Factor 2 from nr 2 and factor 0.5 of border integration cnacel out
#            s = thesign(simdata.indices[1], tp, simdata.NstepsinMemory) 
#            @. simdata.IFk[tp] += real( s * simdata.IFk[simdata.indices[1]] * simdata.Fr[simdata.indices[1], tp] )
#        end
#        # multiply dt factor
#        simdata.IFk[tp] .*= disc.dt
#        # add local contribution
#        simdata.IFk[tp] .+= simdata.F2kr2[tmone,tp]
#        # multiply expansion parameter
#        simdata.IFk[tp] .*= factor1
#        # make fft
#        #copytolattice!( tmpdata.tmpfull[1], simdata.IFk[tp], disc.fftwhelper)
#        #mul!(simdata.IFx[tp], tmpdata.ftplan, tmpdata.tmpfull[1] )
#        #simdata.IFx[tp] .*= disc.ivol
#
#        tmpdata.RHS[1] .= simdata.IFk[tp]
#        tmpdata.RHS[1] .*= tmpdata.isofactortox
#        mul!(simdata.IFx[tp], tmpdata.ftplan, tmpdata.RHS[1] )
#        simdata.IFx[tp] .*= tmpdata.normfactortox
#    end
#end

function calcSelfEnergyattp!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUfull, tmpdata::TwoPIScalarTmpDataCPUfull, tmone::Int64, tp::Int64, ichunk::Int64)
    mul!(tmpdata.Fx[ichunk], tmpdata.ftplan, simdata.F[tmone,tp] )
    mul!(tmpdata.rx[ichunk], tmpdata.ftplan, simdata.r[tmone,tp] )
    s = thesign(tmone, tp, simdata.NstepsinMemory)
    tmpdata.rx[ichunk] .*= s

    @. tmpdata.Fx2[ichunk] = tmpdata.Fx[ichunk] * simdata.Irx[tp] +  tmpdata.rx[ichunk] * simdata.IFx[tp] 
    mul!(simdata.Sigr[tp], tmpdata.iftplan, tmpdata.Fx2[ichunk] )
    simdata.Sigr[tp] .*= -model.Lambda/(3. * model.ONgroup) * disc.ivol

    @. tmpdata.Fx2[ichunk] = tmpdata.Fx[ichunk] .* simdata.IFx[tp] .- 0.25 * tmpdata.rx[ichunk] .* simdata.Irx[tp]
    mul!(simdata.SigF[tp], tmpdata.iftplan, tmpdata.Fx2[ichunk] )
    simdata.SigF[tp] .*= -model.Lambda/(3. * model.ONgroup) * disc.ivol
    return
end


function calcSelfEnergyattp!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUred, tmpdata::TwoPIScalarTmpDataCPUred, tmone::Int64, tp::Int64, ichunk::Int64)
    copytolattice!( tmpdata.r[ichunk], simdata.r[tmone,tp], disc.fftwhelper)
    copytolattice!( tmpdata.F[ichunk], simdata.F[tmone,tp], disc.fftwhelper)

    mul!(tmpdata.Fx[ichunk], tmpdata.ftplan, tmpdata.F[ichunk] )
    mul!(tmpdata.rx[ichunk], tmpdata.ftplan, tmpdata.r[ichunk] )
    s = thesign(tmone, tp, simdata.NstepsinMemory)
    tmpdata.rx[ichunk] .*= s

    @. tmpdata.Fx2[ichunk] = tmpdata.Fx[ichunk] * simdata.Irx[tp] +  tmpdata.rx[ichunk] * simdata.IFx[tp] 
    mul!(tmpdata.tmpfull[ichunk], tmpdata.iftplan, tmpdata.Fx2[ichunk] )
    tmpdata.tmpfull[ichunk] .*= -model.Lambda/(3. * model.ONgroup) * disc.ivol
    copytoreducedlattice!( simdata.Sigr[tp], tmpdata.tmpfull[ichunk], disc.fftwhelper)

    @. tmpdata.Fx2[ichunk] = tmpdata.Fx[ichunk] * simdata.IFx[tp] - 0.25 * tmpdata.rx[ichunk] * simdata.Irx[tp]
    mul!(tmpdata.tmpfull[ichunk], tmpdata.iftplan, tmpdata.Fx2[ichunk] )
    tmpdata.tmpfull[ichunk] .*= -model.Lambda/(3. * model.ONgroup) * disc.ivol
    copytoreducedlattice!( simdata.SigF[tp], tmpdata.tmpfull[ichunk], disc.fftwhelper)
    return
end

function calcSelfEnergyattp!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUred2, tmpdata::TwoPIScalarTmpDataCPUred, tmone::Int64, tp::Int64, ichunk::Int64)
    factor1 = -model.Lambda/(3. * model.ONgroup)
    s = thesign(tmone, tp, simdata.NstepsinMemory)

    @. tmpdata.tmpx[ichunk] = simdata.Fx[tp] * simdata.Irx[tp] +  s*simdata.rx[tp] * simdata.IFx[tp] 
    mul!(tmpdata.tmpfull[ichunk], tmpdata.iftplan, tmpdata.tmpx[ichunk] )
    tmpdata.tmpfull[ichunk] .*=  factor1 
    copytoreducedlattice!( simdata.Sigr[tp], tmpdata.tmpfull[ichunk], disc.fftwhelper)

    @. tmpdata.tmpx[ichunk] = simdata.Fx[tp] * simdata.IFx[tp] - 0.25 * s* simdata.rx[tp] * simdata.Irx[tp]
    mul!(tmpdata.tmpfull[ichunk], tmpdata.iftplan, tmpdata.tmpx[ichunk] )
    tmpdata.tmpfull[ichunk] .*=  factor1 
    copytoreducedlattice!( simdata.SigF[tp], tmpdata.tmpfull[ichunk], disc.fftwhelper)
    return
end

function calcSelfEnergyattpsimple!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUfull, tmpdata::TwoPIScalarTmpDataCPUfull, tmone::Int64, tp::Int64, ichunk::Int64)
    simdata.Sigr[tp] .= -(model.Lambda/(3. * model.ONgroup)) * real(fft( (bfft(simdata.F[tmone,tp])) .* simdata.Irx[tp] .+        bfft(simdata.r[tmone,tp] * thesign(tmone, tp, simdata.NstepsinMemory)) .* simdata.IFx[tp] ))
    simdata.Sigr[tp] .*= disc.ivol
    simdata.SigF[tp] .= -(model.Lambda/(3. * model.ONgroup)) * real(fft( (bfft(simdata.F[tmone,tp])) .* simdata.IFx[tp] .- 0.25 * bfft(simdata.r[tmone,tp] * thesign(tmone, tp, simdata.NstepsinMemory)) .* simdata.Irx[tp] ))
    simdata.SigF[tp] .*= disc.ivol
end

#function calcSelfEnergyattp!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNinverseNLO, disc::TwoPIScalarDiscretization, simdata::TwoPIScalarSimDataCPUcont, tmpdata::TwoPIScalarTmpDataCPUcont, tmone::Int64, tp::Int64, ichunk::Int64)
#    factor1 = -model.Lambda/(3. * model.ONgroup)
#    s = thesign(tmone, tp, simdata.NstepsinMemory)
#
#    @. tmpdata.RHS[ichunk] = simdata.Fx[tp] * simdata.Irx[tp] +  s*simdata.rx[tp] * simdata.IFx[tp] 
#    tmpdata.RHS[ichunk] .*= tmpdata.isofactortop
#    mul!(simdata.Sigr[tp], tmpdata.iftplan, tmpdata.RHS[ichunk])
#    simdata.Sigr[tp] .*= factor1 * tmpdata.normfactortop
#
#    @. tmpdata.RHS[ichunk] = simdata.Fx[tp] * simdata.IFx[tp] - 0.25 * s* simdata.rx[tp] * simdata.Irx[tp]
#    tmpdata.RHS[ichunk] .*= tmpdata.isofactortop
#    mul!(simdata.Sigr[tp], tmpdata.iftplan, tmpdata.RHS[ichunk])
#    simdata.SigF[tp] .*= factor1 * tmpdata.normfactortop
#    return
#end
###############################################################################################################################################
# RHS functions
###############################################################################################################################################
export getRHS_F!
function getRHS_F!(model::TwoPIScalarPhi4, pexp::TwoPIScalarLO, simdata::TwoPIScalarSimData, disc::TwoPIScalarDiscretization, Mem::lattice, t::Int64, tp::Int64)
    Mem .= 0
end
function getRHS_F!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNLO, simdata::TwoPIScalarSimDataCPUfull, disc::TwoPIScalarDiscretization, Mem::lattice, t::Int64, tp::Int64)
    # calculates the RHS for solving for F(t,tp). Gives R_F(t-1,tp)
    # Brute force implementation - no symmetries taken into account
    # Note: For F(t,t) we need R_F(t-1,t). No special treatment required
    tpmin = simdata.indices[1]
    # first part
    Mem  .= simdata.Sigr[tpmin] .* simdata.F[tpmin,tp]
    Mem  .*= -0.5
    for z in tpmin+1:t-2
        Mem  .-= simdata.Sigr[z] .* simdata.F[z,tp]
    end
    # second part
    Mem  .+= 0.5 .* simdata.SigF[tpmin] .* simdata.r[tpmin,tp] .* thesign(tpmin, tp, simdata.NstepsinMemory) 
    for z in tpmin+1:tp-1
        Mem  .+= simdata.SigF[z] .* simdata.r[z,tp] .* thesign(z, tp, simdata.NstepsinMemory) 
    end
end

function getRHS_F!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNLO, simdata::TwoPIScalarSimDataCPUred, disc::TwoPIScalarDiscretization, Mem::reducedlattice, t::Int64, tp::Int64)
    # calculates the RHS for solving for F(t,tp). Gives R_F(t-1,tp)
    # Brute force implementation - no symmetries taken into account
    # Note: For F(t,t) we need R_F(t-1,t). No special treatment required
    ## TO DO: REDUCE SIG
    tpmin = simdata.indices[1]
    # first part
    Mem  .= simdata.Sigr[tpmin] .* simdata.F[tpmin,tp]
    Mem .*= -0.5
    for z in tpmin+1:t-2
        Mem  .-= simdata.Sigr[z] .* simdata.F[z,tp]
    end
    # second part
    Mem  .+= 0.5 .* simdata.SigF[tpmin] .* simdata.r[tpmin,tp] .* thesign(tpmin, tp, simdata.NstepsinMemory)
    for z in tpmin+1:tp-1
        Mem  .+= simdata.SigF[z] .* simdata.r[z,tp] .* thesign(z, tp, simdata.NstepsinMemory)
    end
end

function getRHS_F!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNLO, simdata::TwoPIScalarSimDataCPUred2, disc::TwoPIScalarDiscretization, Mem::reducedlattice, t::Int64, tp::Int64)
    # calculates the RHS for solving for F(t,tp). Gives R_F(t-1,tp)
    # Brute force implementation - no symmetries taken into account
    # Note: For F(t,t) we need R_F(t-1,t). No special treatment required
    ## TO DO: REDUCE SIG
    tpmin = simdata.indices[1]
    # first part
    Mem  .= simdata.Sigr[tpmin] .* simdata.F[tpmin,tp]
    Mem .*= -0.5
    for z in tpmin+1:t-2
        Mem  .-= simdata.Sigr[z] .* simdata.F[z,tp]
    end
    # second part
    Mem  .+= 0.5 .* simdata.SigF[tpmin] .* simdata.r[tpmin,tp] .* thesign(tpmin, tp, simdata.NstepsinMemory)
    for z in tpmin+1:tp-1
        Mem  .+= simdata.SigF[z] .* simdata.r[z,tp] .* thesign(z, tp, simdata.NstepsinMemory)
    end
end

function getRHS_F!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNLO, simdata::TwoPIScalarSimDataCPUcont, disc::TwoPIScalarDiscretization, Mem::reducedlattice, t::Int64, tp::Int64)
    # calculates the RHS for solving for F(t,tp). Gives R_F(t-1,tp)
    # Brute force implementation - no symmetries taken into account
    # Note: For F(t,t) we need R_F(t-1,t). No special treatment required
    ## TO DO: REDUCE SIG
    tpmin = simdata.indices[1]
    # first part
    Mem  .= simdata.Sigr[tpmin] .* simdata.F[tpmin,tp]
    Mem .*= -0.5
    for z in tpmin+1:t-2
        Mem  .-= simdata.Sigr[z] .* simdata.F[z,tp]
    end
    # second part
    Mem  .+= 0.5 .* simdata.SigF[tpmin] .* simdata.r[tpmin,tp] .* thesign(tpmin, tp, simdata.NstepsinMemory)
    for z in tpmin+1:tp-1
        Mem  .+= simdata.SigF[z] .* simdata.r[z,tp] .* thesign(z, tp, simdata.NstepsinMemory)
    end
end

export getRHS_r!
function getRHS_r!(model::TwoPIScalarPhi4, pexp::TwoPIScalarLO, simdata::TwoPIScalarSimData, disc::TwoPIScalarDiscretization, Mem::lattice, t::Int64, tp::Int64)
    Mem .= 0
end

function getRHS_r!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNLO, simdata::TwoPIScalarSimDataCPUfull, disc::TwoPIScalarDiscretization, Mem::lattice, t::Int64, tp::Int64)
    # calculates the RHS for solving for r(t,tp). Gives R_r(t-1,tp)
    # Brute force implementation - no symmetries taken into account
    # Note: r(t,t) is set by hand -> R_r(t-1,t) does not need to be calculated (tp=t > t-1)
    Mem .= 0
    for z in tp+1:t-2
        Mem .-= simdata.Sigr[z] .* simdata.r[z,tp] .* thesign(z, tp, simdata.NstepsinMemory)
    end
end

function getRHS_r!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNLO, simdata::TwoPIScalarSimDataCPUred, disc::TwoPIScalarDiscretization, Mem::reducedlattice, t::Int64, tp::Int64)
    ## TO DO: REDUCE SIG
    # calculates the RHS for solving for r(t,tp). Gives R_r(t-1,tp)
    # Brute force implementation - no symmetries taken into account
    # Note: r(t,t) is set by hand -> R_r(t-1,t) does not need to be calculated (tp=t > t-1)
    Mem .= 0
    for z in tp+1:t-2
        Mem .-= simdata.Sigr[z] .* simdata.r[z,tp] .* thesign(z, tp, simdata.NstepsinMemory) 
    end
end

function getRHS_r!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNLO, simdata::TwoPIScalarSimDataCPUred2, disc::TwoPIScalarDiscretization, Mem::reducedlattice, t::Int64, tp::Int64)
    ## TO DO: REDUCE SIG
    # calculates the RHS for solving for r(t,tp). Gives R_r(t-1,tp)
    # Brute force implementation - no symmetries taken into account
    # Note: r(t,t) is set by hand -> R_r(t-1,t) does not need to be calculated (tp=t > t-1)
    Mem .= 0
    for z in tp+1:t-2
        Mem .-= simdata.Sigr[z] .* simdata.r[z,tp] .* thesign(z, tp, simdata.NstepsinMemory) 
    end
end

function getRHS_r!(model::TwoPIScalarPhi4, pexp::TwoPIScalarNLO, simdata::TwoPIScalarSimDataCPUcont, disc::TwoPIScalarDiscretization, Mem::reducedlattice, t::Int64, tp::Int64)
    ## TO DO: REDUCE SIG
    # calculates the RHS for solving for r(t,tp). Gives R_r(t-1,tp)
    # Brute force implementation - no symmetries taken into account
    # Note: r(t,t) is set by hand -> R_r(t-1,t) does not need to be calculated (tp=t > t-1)
    Mem .= 0
    for z in tp+1:t-2
        Mem .-= simdata.Sigr[z] .* simdata.r[z,tp] .* thesign(z, tp, simdata.NstepsinMemory) 
    end
end

export evolve_r!
function evolve_r!(model::TwoPIScalarPhi4, simdata::TwoPIScalarSimData, disc::TwoPIScalarDiscretization, RHS, t::Int64, tp::Int64)
    if tp==t-1 # phi pi comutator
        #for i in 1:disc.Nx^disc.sdim simdata.r[t,tp,i] = thesign(t, tp, simdata.NstepsinMemory) * disc.dt end
        @. simdata.r[t,tp] = thesign(t, tp, simdata.NstepsinMemory)*disc.dt
    elseif tp==t # diagonal
        #for i in 1:disc.Nx^disc.sdim simdata.r[t,tp,i] = 0 end
        @. simdata.r[t,tp] = 0
    else
        #for i in 1:disc.Nx^disc.sdim simdata.r[t,tp,i] = thesign(t, tp, simdata.NstepsinMemory) * (( 2 - disc.dt^2 * simdata.omega2values[i] ) * simdata.r[t-1,tp,i] * thesign(t-1, tp, simdata.NstepsinMemory)- simdata.r[t-2,tp,i] * thesign(t-2, tp, simdata.NstepsinMemory) + RHS[i] * disc.dt^3 ) end
        #@. simdata.r[t,tp] = thesign(t, tp, simdata.NstepsinMemory) * ( ( 2 - disc.dt^2 * simdata.omega2values ) * simdata.r[t-1,tp] * thesign(t-1, tp, simdata.NstepsinMemory) - simdata.r[t-2,tp] * thesign(t-2, tp, simdata.NstepsinMemory) + RHS * disc.dt^3)
        # optimized
        RHS .*= disc.dt^3
        @. simdata.r[t,tp] =  ( 2 - simdata.omega2values ) * simdata.r[t-1,tp] * thesign(t-1, tp, simdata.NstepsinMemory) - simdata.r[t-2,tp] * thesign(t-2, tp, simdata.NstepsinMemory) + RHS
        simdata.r[t,tp] .*= thesign(t, tp, simdata.NstepsinMemory)
    end
end

export evolve_F!
function evolve_F!(model::TwoPIScalarPhi4,simdata::TwoPIScalarSimData,disc::TwoPIScalarDiscretization, RHS, t::Int64, tp::Int64)
    # to avoid setting values via r/F[t,tp,i] get pointers to lattice and then set the sites with bc
    #for i in 1:disc.Nx^disc.sdim simdata.F[t,tp,i] = ( 2 - disc.dt^2 * simdata.omega2values[i] ) * simdata.F[t-1,tp,i] - simdata.F[t-2,tp,i] + RHS[i] * disc.dt^3 end
    #@. simdata.F[t,tp] = ( 2 - disc.dt^2 * simdata.omega2values ) * simdata.F[t-1,tp] - simdata.F[t-2,tp] + RHS * disc.dt^3
    #optimized
    RHS .*= disc.dt^3
    @. simdata.F[t,tp] = ( 2 - simdata.omega2values ) * simdata.F[t-1,tp] - simdata.F[t-2,tp] + RHS
end

export setomega2values!
function setomega2values!(model::TwoPIScalarPhi4,simdata::TwoPIScalarSimData,disc::TwoPIScalarDiscretization)
    @. simdata.omega2values = simdata.kL2values + simdata.hartreemass2[1]
    simdata.omega2values .*= disc.dt^2
end

export evolve!
function evolve!( thesolution::QFTdynamicsSolutionTwoPIScalar, tmpdata::TwoPIScalarTmpData, t)
    @unpack problem, simdata, measurearray = thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem

    expandSimData!(simdata) 

    # calculate the range for individual threads
    nrofindices = (simdata.indices[2]-1) - simdata.indices[1] + 1
    tmpdata.threadranges .= splitter(simdata.indices[1], nrofindices, num.nchunks)

    # precalc stuff at tmone! - only for offdiagonal elements
    @Threads.threads for ichunk in 1:tmpdata.nchunks
        for tp in tmpdata.threadranges[ichunk]
            calc_Fr_F2kr2!(t-1, tp, simdata, tmpdata, disc, ichunk)
        end
    end

    # calcSelfEnergies at tmone
    calcSelfEnergies!(model, pexp, disc, simdata, tmpdata, t-1)

    # calc Mass
    simdata.hartreemass2[1] = getHartreeMass2(model, pexp, disc, simdata, t-1, t-1)
    # set omega2 values
    setomega2values!(model, simdata, disc)

    #print("Calculating offdiagonal F/g[", t, ",tp], where tp=")
    # Evolve off diagonal ones
    @Threads.threads for ichunk in 1:tmpdata.nchunks
        for tp in tmpdata.threadranges[ichunk]
            #print(tp,", ")
        #for tp in simdata.indices[1]:simdata.indices[2]-1
        #    ichunk =1
            # r
            getRHS_r!(model, pexp, simdata, disc, tmpdata.RHS[ichunk], t, tp)
            evolve_r!(model, simdata, disc, tmpdata.RHS[ichunk], t, tp)
            # F
            getRHS_F!(model, pexp, simdata, disc, tmpdata.RHS[ichunk], t, tp)
            evolve_F!(model, simdata, disc, tmpdata.RHS[ichunk], t, tp)

            # if CPUred2 calc Fr and F2kr2 
            #calc_Fr_F2kr2!(t,tp,simdata, tmpdata, disc, ichunk)

        end
    end

    # Evolve diagonal ones
    # r
    evolve_r!(model, simdata, disc, tmpdata.RHS[1], t, t)
    # F
    getRHS_F!(model, pexp, simdata, disc, tmpdata.RHS[1], t, t)
    evolve_F!(model, simdata, disc, tmpdata.RHS[1], t, t)
    
    # if CPUred2 calc Fr and F2kr2 
    #calc_Fr_F2kr2!(t,t,simdata, tmpdata, disc, 1)

    #update last EvolStep
    simsetup.lastEvolStep = t 
end
