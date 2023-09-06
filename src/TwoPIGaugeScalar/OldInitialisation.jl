# old mass initialisation file stuff
############################################################################################################
# Particle numbers
############################################################################################################
# export getparticlenr
# function getparticlenr(init::TwoPIGaugeScalarParticle, disc::TwoPIGaugeScalarDiscretization)
#     return init.n .* ones(length(disc.fftwhelper))
# end

# function getparticlenr(init::TwoPIGaugeScalarThermal, disc::TwoPIGaugeScalarDiscretization)
#     return [ 1/(exp(sqrt(disc.fftwhelper[i].lev2 + disc.Mass^2)/init.T)-1) for i in 1:length(disc.fftwhelper)]
# end

# function getparticlenr(init::TwoPIGaugeScalarTopHatT1, disc::TwoPIGaugeScalarDiscretization)
#     n = zeros(length(disc.fftwhelper))
#     k2min = 0.68 * disc.sdim
#     idxmin = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2min))[2]
#     k2max = 2.04 * disc.sdim
#     idxmax = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2max))[2]
#     eta = 2
#     for i in idxmin:idxmax
#         n[i] = eta 
#     end
#     return n
# end

# function getparticlenr(init::TwoPIGaugeScalarTopHatT2, disc::TwoPIGaugeScalarDiscretization)
#     n = zeros(length(disc.fftwhelper))
#     k2min = 0 * disc.sdim
#     idxmin = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2min))[2]
#     k2max = 1.9 * disc.sdim
#     idxmax = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2max))[2]
#     eta = 1.85
#     for i in idxmin:idxmax
#         n[i] = eta 
#     end
#     return n
# end

# function getparticlenr(init::TwoPIGaugeScalarTopHatT3, disc::TwoPIGaugeScalarDiscretization)
#     n = zeros(length(disc.fftwhelper))
#     k2min = 2.04 * disc.sdim
#     idxmin = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2min))[2]
#     k2max = 2.72 * disc.sdim
#     idxmax = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2max))[2]
#     eta = 1.6
#     for i in idxmin:idxmax
#         n[i] = eta 
#     end
#     return n
# end

############################################################################################################
# Masses
############################################################################################################
# export setmass2values!
# function setmass2values!(model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, t::Int64)
#     @unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, scalarmass2, transvmass2, longitmass2, fftplan, bfftplan = simdata
#     @unpack angprdterm, fun, tmplatvec = tmpdata

#     scalarloop = sum( simdata.FS[t,t] ) / sqrt(disc.vol)
#     gaugeloopT = sum( simdata.FT[t,t] ) / sqrt(disc.vol)
#     gaugeloopL = sum( simdata.FL[t,t] ) / sqrt(disc.vol)

#     # scalar mass
#     scalrcontr = 8 * model.Lambda * scalarloop                          # scalar contribution to scalar mass
#     gaugecontr = 4 * model.CF * model.g^2 * ( (disc.sdim-1)*gaugeloopT + gaugeloopL ) # gauge contribution to scalar mass

#     # gauge field masses
#     # ang product terms (introduce p,q,k for pedagogical labelling/motivation)
#     p2 = k2values # represents external momentum
#     q2 = k2values # loop lattice momentum
#     k2 = k2values # k==p-q --> separately handled loop momentum
#     k4 = k4values
#     p2inv = invk2values
#     q2inv = invk2values
#     k2inv = invk2values
#     fun .= FT[t,t] - FL[t,t]
#     @show fun 
#     fun[1] = 0. # this is not true in general, but we want to construct tmpvec such that the zero loop momentum does not contribute (we set q2inv[1]=0) to the values of tmpvec
#     tmplatvec[1] =   bfftplan*(fun .* q2inv) / sqrt(disc.vol)
#     tmplatvec[2] = 2*bfftplan*(fun)/sqrt(disc.vol) + (-2)*bfftplan*(fun .* q2inv)/sqrt(disc.vol) .* (bfftplan*(k2))/sqrt(disc.vol)
#     tmplatvec[3] =   bfftplan*(q2 .* fun)/sqrt(disc.vol) + (bfftplan*(fun .* k2inv))/sqrt(disc.vol) .* (bfftplan*(k4))/sqrt(disc.vol) + (-2)*(bfftplan*(fun))/sqrt(disc.vol) .* (bfftplan*(k2))/sqrt(disc.vol)
#     #tmpvec = real([ fftplan*(   (bfftplan*(fun .* q2inv))/ sqrt(disc.vol) ), 
#     #                fftplan*( 2*(bfftplan*(fun))/ sqrt(disc.vol) + (-2)*(bfftplan*(fun .* q2inv))/ sqrt(disc.vol) .* (bfftplan*(k2)) / sqrt(disc.vol)), 
#     #                fftplan*(   (bfftplan*(q2 .* fun))/ sqrt(disc.vol) + (bfftplan*(fun .* k2inv))/ sqrt(disc.vol) .* (bfftplan*(k4))/ sqrt(disc.vol) + (-2)*(bfftplan*(fun))/ sqrt(disc.vol) .* (bfftplan*(k2))/ sqrt(disc.vol) )
#     #                ]) /sqrt(vol)
#     # use of fft plan: fft(A) -->> simdata.fftplan * A
#     angprdterm .= (1/4) * real( p2 .* (fftplan*tmplatvec[1]) + fftplan*tmplatvec[2] + p2inv .* (fftplan*tmplatvec[3]) ) / sqrt(disc.vol)
#     #@show tmplatvec[1] 
#     #angprdterm = (1/4) * ( p2 .* tmpvec[1] + tmpvec[2] + p2inv .* tmpvec[3] )
#     angprdterm[1] = 0.
#     scalarmass2[1] = model.Mass^2 + scalrcontr  + gaugecontr    # total scalar mass
#     transvmass2 .= 2/model.N * scalarloop .+  model.N * ( (disc.sdim-2)*gaugeloopT + gaugeloopL ) .+ model.N/(disc.dim-1) .* (gaugeloopT .- gaugeloopL .- angprdterm)
#     longitmass2 .= 2/model.N * scalarloop .+  model.N * ( (disc.sdim-2)*gaugeloopT + gaugeloopL ) .+ model.N * angprdterm
#     transvmass2 .*= model.g^2
#     longitmass2 .*= model.g^2
#     #@show scalarmass2[1]
# end
# export setmass2valuesOld!
# function setmass2valuesOld!(model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization, simdata::TwoPIGaugeScalarSimData, t::Int64)
#     @unpack FS, FT, FL, rS, rT, rL, k2values, k4values, k6values, invk2values, scalarmass2, transvmass2, longitmass2, fftplan, bfftplan = simdata
#     scalarloop = sum( simdata.FS[t,t] ) * disc.ivol #/ (disc.Nx^disc.sdim)
#     gaugeloopT = sum( simdata.FT[t,t] ) * disc.ivol #/ (disc.Nx^disc.sdim)
#     gaugeloopL = sum( simdata.FL[t,t] ) * disc.ivol #/ (disc.Nx^disc.sdim)
#     # scalar mass
#     scalrcontr = 8 * model.Lambda * scalarloop                          # scalar contribution to scalar mass
#     gaugecontr = 4 * model.CF * model.g^2 * ( (disc.dim-1)*gaugeloopT + gaugeloopL ) # gauge contribution to scalar mass
#     #simdata.scalarmass2[1] = model.Mass^2 + scalrcontr  + gaugecontr    # total scalar mass
#     # gauge field masses
#     # ang product terms (introduce p,q,k for pedagogical labelling/motivation)
#     p2 = k2values # represents external momentum
#     q2 = k2values # loop lattice momentum
#     k2 = k2values # k==p-q --> separately handled loop momentum
#     k4 = k4values
#     p2inv = invk2values
#     q2inv = invk2values
#     k2inv = invk2values
#     fun = FT[t,t] - FL[t,t]
#     fun[1] = 0. # this is not true in general, but we want to construct tmpvec such that the zero loop momentum does not contribute (we set q2inv[1]=0) to the values of tmpvec
#     tmpvec = real([ fftplan*(   (bfftplan*(fun .* q2inv)) * disc.vol ), 
#                     fftplan*( 2*(bfftplan*(fun))          * disc.vol + (-2)*(bfftplan*(fun .* q2inv)) .* (bfftplan*(k2)) ), 
#                     fftplan*(   (bfftplan*(q2 .* fun))    * disc.vol +      (bfftplan*(fun .* k2inv)) .* (bfftplan*(k4)) + (-2)*(bfftplan*(fun)) .* (bfftplan*(k2)) )
#                     ]) * disc.ivol^3 #/ (disc.Nx^disc.sdim)^3
#     # use of fft plan: fft(A) -->> simdata.fftplan * A
    
#     angprdterm = (1/4) * ( p2 .* tmpvec[1] + tmpvec[2] + p2inv .* tmpvec[3] )
#     angprdterm[1] = 0.
#     scalarmass2[1] = model.Mass^2 + scalrcontr  + gaugecontr    # total scalar mass
#     transvmass2 .= 2/model.N * scalarloop .+  model.N * ( (disc.dim-2)*gaugeloopT + gaugeloopL ) .+ model.N/(disc.dim-1) .* (gaugeloopT .- gaugeloopL .- angprdterm)
#     longitmass2 .= 2/model.N * scalarloop .+  model.N * ( (disc.dim-2)*gaugeloopT + gaugeloopL ) .+ model.N * angprdterm
#     transvmass2 .*= model.g^2
#     longitmass2 .*= model.g^2
#     @show scalarmass2[1]
# end

##########################################################################################################
#   SEPARATE SETTERS EXPERIMENT
# export setSCALARmass2!
# function setSCALARmass2!(model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization, simdata::TwoPIGaugeScalarSimData, t::Int64)
#     #d = disc.sdim+1 
#     CF = (model.N^2-1)/(2*model.N) #3/4 # for SU(2) CF = (N^2-1)/(2N)
#     scalarloop = sum( simdata.FS[t,t] ) / (disc.Nx^disc.sdim)
#     gaugeloopT = sum( simdata.FT[t,t] ) / (disc.Nx^disc.sdim)
#     gaugeloopL = sum( simdata.FL[t,t] ) / (disc.Nx^disc.sdim)
#     # scalar mass
#     scalrcontr = 8 * model.Lambda * scalarloop                          # scalar contribution to scalar mass
#     gaugecontr = 4 * CF * model.g^2 * ( (disc.sdim-1)*gaugeloopT + gaugeloopL ) # gauge contribution to scalar mass
    
#     simdata.scalarmass2[1] = model.Mass^2 + scalrcontr  + gaugecontr    # total scalar mass
# end

# export setLONGITmass2!
# function setLONGITmass2!(model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization, simdata::TwoPIGaugeScalarSimData, t::Int64)# gauge field masses
#     #d = disc.sdim+1 
#     scalarloop = sum( simdata.FS[t,t] ) / (disc.Nx^disc.sdim)
#     gaugeloopT = sum( simdata.FT[t,t] ) / (disc.Nx^disc.sdim)
#     gaugeloopL = sum( simdata.FL[t,t] ) / (disc.Nx^disc.sdim)
#     # ang product terms (introduce p,q,k for pedagogical labelling/motivation)
#     p2 = simdata.k2values # represents external momentum
#     q2 = simdata.k2values # loop lattice momentum
#     k2 = simdata.k2values # k==p-q --> separately handled loop momentum
#     k4 = simdata.k4values
#     p2inv = simdata.invk2values
#     q2inv = simdata.invk2values
#     k2inv = simdata.invk2values
#     fun = simdata.FT[t,t] - simdata.FL[t,t]
#     fun[1] = 0. # this is not true in general, but we want to construct tmpvec such that the zero loop momentum does not contribute (we set q2inv[1]=0) to the values of tmpvec
#     tmpvec = real([fft(bfft(fun .* q2inv)*(disc.Nx^disc.sdim)), 
#                     fft(2*bfft(fun)*(disc.Nx^disc.sdim) + (-2)*bfft(fun .* q2inv) .* bfft(k2)), 
#                     fft(bfft(q2.*fun)*(disc.Nx^disc.sdim) + bfft(fun .* k2inv) .* bfft(k4) + (-2)*bfft(fun) .* bfft(k2))])/ (disc.Nx^disc.sdim)^2

#     tmp =  [fft(ifft(fun .* q2inv)),
#             2*fft(ifft(fun)) + (-2)*fft(ifft(fun .* q2inv) .* ifft(k2)),
#             fft(ifft(q2.*fun)) + fft(ifft(fun .* q2inv) .* ifft(k4)) + (-2)*fft(ifft(fun).*ifft(k2))]

#     angprdterm = (1/4) * real(p2 .* tmpvec[1] + tmpvec[2] + p2inv .* tmpvec[3])
#     angprdterm[1] = 0.
#     simdata.longitmass2 .= 2/model.N * scalarloop .+  model.N * ( (d-2)*gaugeloopT + gaugeloopL ) .+ model.N * angprdterm
#     simdata.longitmass2 .*= model.g^2
# end







############################################################################################################
# Masses (old version)
############################################################################################################
# export getScalarMass2
# function getScalarMass2(model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization, simdata::TwoPIGaugeScalarSimData, t::Int64)
#     Scalrcont = 0.
#     Gaugecont = 0.
#     d = disc.sdim+1 
#     CF = (model.N^2-1)/(2*model.N) #3/4 # for SU(2) CF = (N^2-1)/(2N)
#     Scalrcont = sum( simdata.FS[t,t] )
#     Gaugecont = sum(  (d-1)*simdata.FT[t,t] + simdata.FL[t,t] )
#     M2 = model.Mass^2 + ( 4*CF*model.g^2*Gaugecont + 8*model.Lambda*Scalrcont )/(disc.Nx^disc.sdim)
#     return M2
# end

# export setTRANSVmass2!
# function setTRANSVmass2!(model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization, simdata::TwoPIGaugeScalarSimData, t::Int64)# gauge field masses
#     d = disc.sdim+1 
#     scalarloop = sum( simdata.FS[t,t] ) / (disc.Nx^disc.sdim)
#     gaugeloopT = sum( simdata.FT[t,t] ) / (disc.Nx^disc.sdim)
#     gaugeloopL = sum( simdata.FL[t,t] ) / (disc.Nx^disc.sdim)
#     Scalrcont = 2/model.N * scalarloop
#     Gaugecont1 =  model.N * ( (d-2)*gaugeloopT + gaugeloopL )
#     # ang product terms (introduce p,q,k for pedagogical labelling/motivation)
#     p2 = simdata.k2values # represents external momentum
#     q2 = simdata.k2values # loop lattice momentum
#     k2 = simdata.k2values # k==p-q --> separately handled loop momentum
#     k4 = simdata.k4values
#     p2inv = simdata.invk2values
#     q2inv = simdata.invk2values
#     k2inv = simdata.invk2values
#     fun = simdata.FT[t,t] - simdata.FL[t,t]
#     fun[1] = 0. # this is not true in general, but we want to construct tmpvec such that the zero loop momentum does not contribute (we set q2inv[1]=0) to the values of tmpvec
#     tmpvec = real([fft(bfft(fun .* q2inv)*(disc.Nx^disc.sdim)), 
#                     fft(2*bfft(fun)*(disc.Nx^disc.sdim) + (-2)*bfft(fun .* q2inv) .* bfft(k2)), 
#                     fft(bfft(q2.*fun)*(disc.Nx^disc.sdim) + bfft(fun .* k2inv) .* bfft(k4) + (-2)*bfft(fun) .* bfft(k2))])/ (disc.Nx^disc.sdim)^2

#     tmp =   real([fft(ifft(fun .* q2inv)),
#             2*fft(ifft(fun)) + (-2)*fft(ifft(fun .* q2inv) .* ifft(k2)),
#             fft(ifft(q2.*fun)) + fft(ifft(fun .* q2inv) .* ifft(k4)) + (-2)*fft(ifft(fun).*ifft(k2))])

#     angprdterm = (1/4) * ( p2 .* tmpvec[1] .+ tmpvec[2] .+ p2inv .* tmpvec[3])
#     angprdterm[1] = 0.
#     #M2 = 0. * simdata.k2values 
#     #M2 .+= (-1)*model.N/(d-1) * angprdterm / (disc.Nx^disc.sdim)
#     simdata.transvmass2 .= (-1)*model.N/(d-1) * angprdterm / (disc.Nx^disc.sdim)                # momentum dependent part
#     simdata.transvmass2 .+= Scalrcont + Gaugecont1 + model.N/(d-1) * (gaugeloopT - gaugeloopL ) # scalar part
#     simdata.transvmass2 .*= model.g^2
#     #simdata.transvmass2 .= M2
# end



# export getTransvMass2NEW
# function getTransvMass2NEW(model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization, simdata::TwoPIGaugeScalarSimData, t::Int64)
#     d = disc.sdim+1 
#     scalarloop = sum( simdata.FS[t,t] ) / (disc.Nx^disc.sdim)
#     gaugeloopT = sum( simdata.FT[t,t] ) / (disc.Nx^disc.sdim)
#     gaugeloopL = sum( simdata.FL[t,t] ) / (disc.Nx^disc.sdim)
#     Scalrcont = 2/model.N * scalarloop
#     Gaugecont1 =  model.N * ( (d-2)*gaugeloopT + gaugeloopL )
#     # ang product terms (introduce p,q,k for pedagogical labelling/motivation)
#     p2 = simdata.k2values # represents external momentum
#     q2 = simdata.k2values # loop lattice momentum
#     k2 = simdata.k2values # k==p-q --> separately handled loop momentum
#     k4 = simdata.k4values
#     p2inv = simdata.invk2values
#     q2inv = simdata.invk2values
#     k2inv = simdata.invk2values
#     fun = simdata.FT[t,t] - simdata.FL[t,t]
#     fun[1] = 0. # this is not true in general, but we want to construct tmpvec such that the zero loop momentum does not contribute (we set q2inv[1]=0) to the values of tmpvec
#     tmpvec = real([fft(bfft(fun .* q2inv)*(disc.Nx^disc.sdim)), 
#                     fft(2*bfft(fun)*(disc.Nx^disc.sdim) + (-2)*bfft(fun .* q2inv) .* bfft(k2)), 
#                     fft(bfft(q2.*fun)*(disc.Nx^disc.sdim) + bfft(fun .* k2inv) .* bfft(k4) + (-2)*bfft(fun) .* bfft(k2))])/ (disc.Nx^disc.sdim)^2

#     tmp =   real([fft(ifft(fun .* q2inv)),
#             2*fft(ifft(fun)) + (-2)*fft(ifft(fun .* q2inv) .* ifft(k2)),
#             fft(ifft(q2.*fun)) + fft(ifft(fun .* q2inv) .* ifft(k4)) + (-2)*fft(ifft(fun).*ifft(k2))])

#     angprdterm = (1/4) * ( p2 .* tmp[1] .+ tmp[2] .+ p2inv .* tmp[3])
#     angprdterm[1] = 0.
#     M2 = 0. * simdata.k2values
#     #M2 = - model.N/(d-1) * angprdterm
#     M2 .+= Scalrcont + Gaugecont1 + model.N/(d-1) * (gaugeloopT - gaugeloopL ) 
#     M2 .+= (-1)*model.N/(d-1) * angprdterm / (disc.Nx^disc.sdim)
#     M2 .*= model.g^2
#     return M2
# end
# export getTransvMass2
# function getTransvMass2(model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization, simdata::TwoPIGaugeScalarSimData, t::Int64)
#     # pure loop terms #
#     d = disc.sdim+1
#     scalarloop = 0.
#     gaugeloopT = 0.
#     gaugeloopL = 0.
#     Gaugecont1 = 0.
#     for i in 1:length(disc.fftwhelper) # integrate over loop == sum over all lattice points (all modes with momentum)
#         for j in 1:disc.fftwhelper[i].deg
#             idx = disc.fftwhelper[i].ind[j]
#             scalarloop += simdata.FS[t,t][idx] # here we add all contributions (all values of FS) from all lattice sites to Scalarcont -> integral over q
#             gaugeloopT += simdata.FT[t,t][idx] # integrate over q
#             gaugeloopL += simdata.FL[t,t][idx] 
#         end
#     end
#     Scalrcont = 2/model.N * scalarloop
#     Gaugecont1 =  model.N * ( (d-2)*gaugeloopT + gaugeloopL )

#     # ang product term #
#     q2 = 0. * simdata.k2values
#     funoverq2 = 0. * simdata.k2values
#     fun = simdata.FT[t, t] - simdata.FL[t, t]
#     fun[1] = 0.
#     for i in 2:length(disc.fftwhelper)          #NB: start from first non-zero momentum!
#         for j in 1:disc.fftwhelper[i].deg
#             idx = disc.fftwhelper[i].ind[j]
#             q2[idx] = disc.fftwhelper[i].lev2
#             funoverq2[idx] = fun[idx] / disc.fftwhelper[i].lev2
#         end
#     end
#     q4 = q2.^2
#     k2 = 1. * q2
#     k4 = 1. * q4

#     # divide vector components into powers of p
#     tmp =  [fft(ifft(funoverq2)),
#             2*fft(ifft(fun)) + (-2)*fft(ifft(funoverq2) .* ifft(k2)),
#             fft(ifft(q2.*fun)) + fft(ifft(funoverq2) .* ifft(k4)) + (-2)*fft(ifft(fun).*ifft(k2))]
#     #@show tmp
#     ######
#     # make M2 a lattice instead
#     angprdterm = 0. * simdata.k2values
#     M2 = 0. * simdata.k2values
#     for i in 1:length(disc.fftwhelper)
#         p2 = disc.fftwhelper[i].lev2
#         for j in 1:disc.fftwhelper[i].deg
#             idx = disc.fftwhelper[i].ind[j]
#             if p2 == 0.
#                 angprdterm[idx] = 0.
#             else
#                 angprdterm[idx] = (1/4) * real(p2 * tmp[1][idx] + tmp[2][idx] + (1/p2) * tmp[3][idx])
#             end
#             M2[idx] = Scalrcont + Gaugecont1 + model.N/(d-1) * (gaugeloopT - gaugeloopL - angprdterm[idx])
#         end
#     end
#     return M2 .* model.g^2/(disc.Nx^disc.sdim)
# end

# export getLongitMass2
# function getLongitMass2(model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization, simdata::TwoPIGaugeScalarSimData, t::Int64)
#     # pure loop terms #
#     d = disc.sdim+1
#     scalarloop = 0
#     gaugeloopT = 0
#     gaugeloopL = 0
#     for i in 1:length(disc.fftwhelper) # integrate over loop == sum over all lattice points (all modes with momentum)
#         for j in 1:disc.fftwhelper[i].deg
#             idx = disc.fftwhelper[i].ind[j]
#             scalarloop += simdata.FS[t,t][idx] # here we add all contributions (all values of FS) from all lattice sites to Scalarcont -> integral over q
#             gaugeloopT += simdata.FT[t,t][idx] # integrate over q
#             gaugeloopL += simdata.FL[t,t][idx] 
#         end
#     end
#     Scalrcont = 2/model.N * scalarloop
#     Gaugecont1 =  model.N * ( (d-2)*gaugeloopT + gaugeloopL )

#     # ang product term #
#     q2 = 0. * simdata.k2values # loop momenta
#     funoverq2 = 0. * simdata.k2values
#     fun = simdata.FT[t, t] - simdata.FL[t, t]
#     for i in 1:length(disc.fftwhelper)
#         for j in 1:disc.fftwhelper[i].deg
#             idx = disc.fftwhelper[i].ind[j]
#             q2[idx] = disc.fftwhelper[i].lev2
#             if q2[idx] == 0.
#                 funoverq2[idx] = 0.
#             else
#                 funoverq2[idx] = fun[idx] / q2[idx]
#             end
#         end
#     end
#     q4 = q2.^2
#     k2 = 1. * q2
#     k4 = 1. * q4

#     # divide vector components into powers of p
#     tmp =  [fft(ifft(funoverq2)),
#             2*fft(ifft(fun)) + (-2)*fft(ifft(funoverq2) .* ifft(k2)),
#             fft(ifft(q2.*fun)) + fft(ifft(funoverq2) .* ifft(k4)) + (-2)*fft(ifft(fun).*ifft(k2))]
    
#     angprdterm = 0. * simdata.k2values
#     M2 = 0. * simdata.k2values
#     for i in 1:length(disc.fftwhelper)
#         p2 = disc.fftwhelper[i].lev2
#         for j in 1:disc.fftwhelper[i].deg
#             idx = disc.fftwhelper[i].ind[j]
#             if p2 == 0.
#                 angprdterm[idx] = 0.
#             else
#                 angprdterm[idx] = (1/4) * real(p2 * tmp[1][idx] + tmp[2][idx] + (1/p2) * tmp[3][idx])
#             end
#             M2[idx] =  Scalrcont + Gaugecont1 + model.N * angprdterm[i]
#         end
#     end
#     return M2 .* model.g^2/(disc.Nx^disc.sdim)
# end
############################################################################################################
export setCSinitTmpData1!
function setCSinitTmpData1!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, model::SUNgaugeScalar, disc::TwoPIGaugeScalarDiscretization) 
    # t=0: determine scalar propagator according to Gauss constraint, compute initial E-field (gauge field A=0 here)
    @unpack fftplan = simdata
    @unpack phix, piix, Pix, Phix, Pik, Phik, PiPik, PiPhik, PhiPhik, Ea, pauli, rhox, rhok, chik, chix = tmpdata
    @unpack vol, Nx = disc
    ### construct Phix
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
        PiPik[idx] = real(tr( adjoint(Pik[idx]) * Pik[idx] ))  # -> dd F -> F[t=2,tp=2]
        PiPhik[idx] = real(tr( adjoint(Pik[idx]) * Phik[idx] )) # -> d F -> F[t=1,tp=0] and F[t=0,tp=1]
        PhiPhik[idx] = real(tr( adjoint(Phik[idx]) * Phik[idx] )) # -> F -> F[t=0,tp=0]
    end
    #
    ## calc rho 
    #pauli4 = [1.0*im 0.0; 0.0 1.0*im]
    for idx in 1:disc.vol
        rhox[1][idx] = -imag(tr(adjoint(Pix[idx]) * pauli[1] * Phix[idx]))
        rhox[2][idx] = -imag(tr(adjoint(Pix[idx]) * pauli[2] * Phix[idx]))
        rhox[3][idx] = -imag(tr(adjoint(Pix[idx]) * pauli[3] * Phix[idx]))
        rhox[4][idx] = -imag(tr(adjoint(Pix[idx]) * im * Phix[idx]))
    end
    rhox .*= -model.g
    @show sum(rhox[1])
    @show sum(rhox[2])
    @show sum(rhox[3])
    @show sum(rhox[4])
    for a in 1:3
        rhok[a] = fftplan * rhox[a] #fft(rhox[a])
        #@show rhok[a][1]               # this should be zero
        chik[a] = -rhok[a] ./ simdata.k2values  # real coefficients here --> chix = ifft(chik) completely real
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

end

export setCSinitTmpData2!
function setCSinitTmpData2!(simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, disc::TwoPIGaugeScalarDiscretization)
    # t=dt (first time step): set A-field according to E-field in previous time step, E(t=0) = ( A(t=dt) - A(t=0) )/dt == A(t=dt)/dt 
    @unpack fftplan = simdata
    @unpack nvalues, Ea, Aaix, Aaik, Dk, DLk, DTk, P_L = tmpdata
    @unpack vol, Nx = disc 
    for a in 1:3
        for i in 1:3
            for idx in 1:vol
                Aaix[a][i][idx] = disc.dt * Ea[idx][i][a]
            end
        end
    end
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

end

## old modelfile construct
# function QFTdynamicsSolution(modelfile::TwoPIGaugeScalarfile, problem::QFTdynamicsProblem)
#     # Get Memory
#     simdata = getTwoPIGaugeSimData(problem)
#     PrintMemoryofSimData!(simdata)
#     #if SimPar.binning != 0 rfftwhelper = rebin(fftwhelper, SimPar.binning) else rfftwhelper=fftwhelper end # bins includes 0 bin!
#     measurearray = Vector{Measurement}(undef, problem.simsetup.Nmeas+1) # includes 0 measurement
#     tmpdata = TwoPIGaugeScalarTmpDataCPUfull(simdata, problem.num.nchunks)

#     # load CS initial conditions as an average over number of runs
#     location = "/Users/magdalenaeriksson/code/2PIcode/data/MCSampledIC_Nx32_sdim3_Mass100_n0_Samples90_B50_ith5_test"
#     CSnbrOfRuns = 50 # this must match CS code
#     for i in 1:CSnbrOfRuns
#         df = CSV.read(location * "/Sample_" * string(i) * ".csv", DataFrame)
#         for idx in 1:problem.disc.vol
#             # phi
#             tmpdata.phix[1][idx] += df.phi1x[idx]
#             tmpdata.phix[2][idx] += df.phi2x[idx]
#             tmpdata.phix[3][idx] += df.phi3x[idx]
#             tmpdata.phix[4][idx] += df.phi4x[idx]
#             #
#             tmpdata.piix[1][idx] += df.pii1x[idx]
#             tmpdata.piix[2][idx] += df.pii2x[idx]
#             tmpdata.piix[3][idx] += df.pii3x[idx]
#             tmpdata.piix[4][idx] += df.pii4x[idx]
#         end
#     end
#     # phi
#     tmpdata.phix[1] ./= CSnbrOfRuns
#     tmpdata.phix[2] ./= CSnbrOfRuns
#     tmpdata.phix[3] ./= CSnbrOfRuns
#     tmpdata.phix[4] ./= CSnbrOfRuns
#     #
#     tmpdata.piix[1] ./= CSnbrOfRuns
#     tmpdata.piix[2] ./= CSnbrOfRuns
#     tmpdata.piix[3] ./= CSnbrOfRuns
#     tmpdata.piix[4] ./= CSnbrOfRuns

    
#     return QFTdynamicsSolutionTwoPIGaugeScalar(problem, simdata, measurearray), tmpdata
# end