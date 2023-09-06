export MeasurementTwoPIGaugeScalar
Base.@kwdef mutable struct MeasurementTwoPIGaugeScalar{N,NstepsinMemory} <: Measurement
    # N = length(disc.fftwhelper) (nbr of distinct momenta)
	time::Float64=0
    E::Float64=0
	E_kinHiggs::Float64=0
    E_kinTrans::Float64=0
    E_kinLong::Float64=0
    E_TrLogHiggs::Float64=0
    E_TrLogGauge::Float64=0
    E_HiggsHartree::Float64=0
    E_MixedHartree::Float64=0
    E_GaugeHartree::Float64=0
    E_MixedSunset::Float64=0
	E_GaugeSunset::Float64=0
	MS2::Float64=0
    MT2::Vector=zeros(Float64, N) # gauge masses defined to be momentum dependent
    ML2::Vector=zeros(Float64, N)
	# from k space
	#phi::Float64=0
	FS::Vector=zeros(Float64, N)
    FT::Vector=zeros(Float64, N)
    FL::Vector=zeros(Float64, N)
	#ddFS::Vector=zeros(Float64, N )
    #ddFT::Vector=zeros(Float64, N )
    #ddFL::Vector=zeros(Float64, N )
	rS::Vector=zeros(Float64, N )
    rT::Vector=zeros(Float64, N )
    rL::Vector=zeros(Float64, N )
	#nS::Vector=zeros(Float64, N )
    #nT::Vector=zeros(Float64, N )
    #nL::Vector=zeros(Float64, N )
	#omegaS::Vector=zeros(Float64, N )
    #omegaT::Vector=zeros(Float64, N )
    #omegaL::Vector=zeros(Float64, N )
	#cS::Vector=zeros(Float64, N )
    #cT::Vector=zeros(Float64, N )
    #cL::Vector=zeros(Float64, N )
    # off-diagonal entries
    FSonet::Vector=zeros(Float64, NstepsinMemory )
    FTonet::Vector=zeros(Float64, NstepsinMemory )
    FLonet::Vector=zeros(Float64, NstepsinMemory )
    rSonet::Vector=zeros(Float64, NstepsinMemory )
    rTonet::Vector=zeros(Float64, NstepsinMemory )
    rLonet::Vector=zeros(Float64, NstepsinMemory )
    ## Sigma (self-energies)
	SigFS::Vector=zeros(Float64, NstepsinMemory )
    SigFT::Vector=zeros(Float64, NstepsinMemory )
    SigFL::Vector=zeros(Float64, NstepsinMemory )
	SigrS::Vector=zeros(Float64, NstepsinMemory )
    SigrT::Vector=zeros(Float64, NstepsinMemory )
    SigrL::Vector=zeros(Float64, NstepsinMemory )
end

export createMeasurement
function createMeasurement(t::Int64, model::SUNgaugeScalar, pexp::TwoPIGaugeScalarLoopEx, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, disc::TwoPIGaugeScalarDiscretization)
	# init measurement object
	meas = MeasurementTwoPIGaugeScalar{length(disc.fftwhelper),disc.NstepsinMemory}()
    # scalar measurements
	meas.time = (t-1) * disc.dt * disc.Mass
	meas.MS2 = simdata.scalarmass2[1]

    # vector measurements
	for i in 1:length(disc.fftwhelper)
	    idx = disc.fftwhelper[i].ind[1] # take degeneracy #1 
        meas.MT2[i] = simdata.transvmass2[idx]
        meas.ML2[i] = simdata.longitmass2[idx]
        # meas.ddFS[i] = ( simdata.FS[t,t][idx] - 2*simdata.FS[t-1,t][idx] + simdata.FS[t-1,t-1][idx] )/disc.dt^2
        # meas.ddFT[i] = ( simdata.FT[t,t][idx] - 2*simdata.FT[t-1,t][idx] + simdata.FT[t-1,t-1][idx] )/disc.dt^2
        # meas.ddFL[i] = ( simdata.FL[t,t][idx] - 2*simdata.FL[t-1,t][idx] + simdata.FL[t-1,t-1][idx] )/disc.dt^2
        meas.FS[i] = simdata.FS[t-1,t-1][idx] 
        meas.FT[i] = simdata.FT[t-1,t-1][idx]
        meas.FL[i] = simdata.FL[t-1,t-1][idx] 
        meas.rS[i] = simdata.rS[t-1,t-1][idx] * thesign(t-1, t-1, simdata.NstepsinMemory)
        meas.rT[i] = simdata.rT[t-1,t-1][idx] * thesign(t-1, t-1, simdata.NstepsinMemory)
        meas.rL[i] = simdata.rL[t-1,t-1][idx] * thesign(t-1, t-1, simdata.NstepsinMemory)
        #meas.omegaS[i] = sqrt( meas.ddFS[i] / simdata.FS[t-1,t-1][idx] )
        #meas.omegaT[i] = sqrt( meas.ddFT[i] / simdata.FT[t-1,t-1][idx] )
        #meas.omegaL[i] = sqrt( meas.ddFL[i] / simdata.FL[t-1,t-1][idx] )
        #meas.cS[i] = sqrt( 1. - 0.25*disc.dt^2*meas.omegaS[i]^2)
        #meas.cT[i] = sqrt( 1. - 0.25*disc.dt^2*meas.omegaT[i]^2)
        #meas.cL[i] = sqrt( 1. - 0.25*disc.dt^2*meas.omegaL[i]^2)
        #meas.nS[i] = meas.cS[i] * sqrt( meas.ddFS[i] * simdata.FS[t-1,t-1][idx] )
        #meas.nT[i] = meas.cS[i] * sqrt( meas.ddFT[i] * simdata.FT[t-1,t-1][idx] )
        #meas.nL[i] = meas.cS[i] * sqrt( meas.ddFL[i] * simdata.FL[t-1,t-1][idx] )
    end
    
    # Energy - depends on vector measurement
    meas.E_kinHiggs = getE_kinHiggs(t, model, meas, pexp, simdata, disc)
    meas.E_kinTrans = getE_kinTrans(t, model, meas, pexp, simdata, disc)
    meas.E_kinLong = getE_kinLong(t, model, meas, pexp, simdata, disc)
    meas.E_TrLogHiggs = getE_TrLogHiggs(t, model, meas, pexp, simdata, tmpdata, disc)
    meas.E_TrLogGauge = getE_TrLogGauge(t, model, meas, pexp, simdata, tmpdata, disc)
    meas.E_HiggsHartree = getE_HiggsHartree(t, model, meas, pexp, simdata, tmpdata, disc)
    meas.E_MixedHartree = getE_MixedHartree(t, model, meas, pexp, simdata, tmpdata, disc)
    meas.E_GaugeHartree = getE_GaugeHartree(t, model, meas, pexp, simdata, tmpdata, disc)
    #meas.E_MixedSunset = getE_MixedSunset(t, model, meas, pexp, simdata, tmpdata, disc)
    #meas.E_GaugeSunset = getE_GaugeSunset(t, model, meas, pexp, simdata, tmpdata, disc)

    meas.E = meas.E_kinHiggs + meas.E_kinTrans + meas.E_kinLong + meas.E_TrLogHiggs + meas.E_TrLogGauge + meas.E_HiggsHartree + meas.E_GaugeHartree #+ meas.E_MixedSunset + meas.E_GaugeSunset
    
    # Sig in Memory past
    whichindex = 2
    meas.SigFS = [ simdata.SigFS[tp][whichindex] for tp in (simdata.SigFS.indices[2]-1):-1:simdata.SigFS.indices[1] ]
    meas.SigFT = [ simdata.SigFT[tp][whichindex] for tp in (simdata.SigFT.indices[2]-1):-1:simdata.SigFT.indices[1] ] 
    meas.SigFL = [ simdata.SigFL[tp][whichindex] for tp in (simdata.SigFL.indices[2]-1):-1:simdata.SigFL.indices[1] ] 
    #@show meas.SigFS 
    meas.SigrS = [ simdata.SigrS[tp][whichindex] for tp in (simdata.SigrS.indices[2]-1):-1:simdata.SigrS.indices[1] ]
    meas.SigrT = [ simdata.SigrT[tp][whichindex] for tp in (simdata.SigrT.indices[2]-1):-1:simdata.SigrT.indices[1] ]
    meas.SigrL = [ simdata.SigrL[tp][whichindex] for tp in (simdata.SigrL.indices[2]-1):-1:simdata.SigrL.indices[1] ] 
    ## off-diagonal
    #meas.FSonet = simdata.FS[simdata.indices[1],simdata.indices[2]][whichindex] # wrong! 
    #meas.rSonet = simdata.rS[simdata.indices[1],simdata.indices[2]][whichindex] # wrong!


    # offdiagonal Gerhards --------------------------------------
    # correlation with oldest one
    #meas.Fonet = [ simdata.F[simdata.indices[1],tp,whichindex] for tp in simdata.indices[1]:(simdata.indices[2]-1)]
    #meas.ronet = [ simdata.r[simdata.indices[1],tp,whichindex]*thesign(simdata.indices[1], tp, simdata.NstepsinMemory) for tp in simdata.indices[1]:(simdata.indices[2]-1)]
    # Gerhards --------------------------------------

    #whichindex = 1

    meas.FSonet = [ simdata.FS[simdata.indices[1],tp][whichindex] for tp in simdata.indices[1]:(simdata.indices[2]-1)]
    meas.FTonet = [ simdata.FT[simdata.indices[1],tp][whichindex] for tp in simdata.indices[1]:(simdata.indices[2]-1)]
    meas.FLonet = [ simdata.FL[simdata.indices[1],tp][whichindex] for tp in simdata.indices[1]:(simdata.indices[2]-1)]
    #@show meas.FSonet 

    meas.rSonet = [ simdata.rS[simdata.indices[1],tp][whichindex]*thesign(simdata.indices[1], tp, simdata.NstepsinMemory) for tp in simdata.indices[1]:(simdata.indices[2]-1)]
    meas.rTonet = [ simdata.rT[simdata.indices[1],tp][whichindex]*thesign(simdata.indices[1], tp, simdata.NstepsinMemory) for tp in simdata.indices[1]:(simdata.indices[2]-1)]
    meas.rLonet = [ simdata.rL[simdata.indices[1],tp][whichindex]*thesign(simdata.indices[1], tp, simdata.NstepsinMemory) for tp in simdata.indices[1]:(simdata.indices[2]-1)]

    # meas.FTonet = [ simdata.FT[t-1,tp][whichindex] for tp in (simdata.indices[2]-1):-1:simdata.indices[1]]
    # meas.FLonet = [ simdata.FL[t-1,tp][whichindex] for tp in (simdata.indices[2]-1):-1:simdata.indices[1]]

    # meas.rSonet = [ simdata.rS[t-1,tp][whichindex] for tp in (simdata.indices[2]-1):-1:simdata.indices[1]]
    # meas.rTonet = [ simdata.rT[t-1,tp][whichindex] for tp in (simdata.indices[2]-1):-1:simdata.indices[1]]
    # meas.rLonet = [ simdata.rL[t-1,tp][whichindex] for tp in (simdata.indices[2]-1):-1:simdata.indices[1]]
    # meas.rSonet = [ simdata.rS[t-1,tp][whichindex]*thesign(t-1, tp, simdata.NstepsinMemory) for tp in (simdata.indices[2]-1):-1:simdata.indices[1]]
    # meas.rTonet = [ simdata.rT[t-1,tp][whichindex]*thesign(t-1, tp, simdata.NstepsinMemory) for tp in (simdata.indices[2]-1):-1:simdata.indices[1]]
    # meas.rLonet = [ simdata.rL[t-1,tp][whichindex]*thesign(t-1, tp, simdata.NstepsinMemory) for tp in (simdata.indices[2]-1):-1:simdata.indices[1]]

    return meas
end

function getE_kinHiggs(t::Int64, model::SUNgaugeScalar, meas::MeasurementTwoPIGaugeScalar, pexp::TwoPIGaugeScalarLoopEx, simdata::TwoPIGaugeScalarSimData, disc::TwoPIGaugeScalarDiscretization)
    ### Kinetic contributions
    E_kinHiggs =  0.5 * sum( simdata.FS[t,t] .- 2*simdata.FS[t-1,t] .+ simdata.FS[t-1,t-1] )/disc.dt^2
    return E_kinHiggs
end

function getE_kinTrans(t::Int64, model::SUNgaugeScalar, meas::MeasurementTwoPIGaugeScalar, pexp::TwoPIGaugeScalarLoopEx, simdata::TwoPIGaugeScalarSimData, disc::TwoPIGaugeScalarDiscretization)
    ### Kinetic contributions
    E_kinTrans = 0.5 * (disc.sdim-1) * sum( simdata.FT[t,t] .- 2*simdata.FT[t-1,t] .+ simdata.FT[t-1,t-1] )/disc.dt^2
    #E_kinTrans = 0.5 * sum( simdata.FT[t,t] .- 2*simdata.FT[t-1,t] .+ simdata.FT[t-1,t-1] )/disc.dt^2
    return E_kinTrans
end

function getE_kinLong(t::Int64, model::SUNgaugeScalar, meas::MeasurementTwoPIGaugeScalar, pexp::TwoPIGaugeScalarLoopEx, simdata::TwoPIGaugeScalarSimData, disc::TwoPIGaugeScalarDiscretization)
    ### Kinetic contributions
    E_kinLong = 0.5 * sum( simdata.FL[t,t] .- 2*simdata.FL[t-1,t] .+ simdata.FL[t-1,t-1] )/disc.dt^2
    return E_kinLong
end

function getE_TrLogHiggs(t::Int64, model::SUNgaugeScalar, meas::MeasurementTwoPIGaugeScalar, pexp::TwoPIGaugeScalarLoopEx, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, disc::TwoPIGaugeScalarDiscretization)
    ### Trlog Higgs contribution
    E_TrLogHiggs =  0.5 * sum( (model.Mass^2 .+ simdata.k2) .* simdata.FS[t-1,t-1] )
    return E_TrLogHiggs
end

function getE_TrLogGauge(t::Int64, model::SUNgaugeScalar, meas::MeasurementTwoPIGaugeScalar, pexp::TwoPIGaugeScalarLoopEx, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, disc::TwoPIGaugeScalarDiscretization)
    ### Trlog transverse contribution
    E_TrLogGauge = 0.5 * (disc.sdim - 1) * sum( simdata.k2 .* simdata.FT[t-1,t-1] )
    #E_TrLogGauge = 0.5 * sum( simdata.k2 .* simdata.FT[t-1,t-1] )
    return E_TrLogGauge
end

function getE_HiggsHartree(t::Int64, model::SUNgaugeScalar, meas::MeasurementTwoPIGaugeScalar, pexp::TwoPIGaugeScalarLoopEx, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, disc::TwoPIGaugeScalarDiscretization)
    ### Higgs Hartree contribution
    E_HiggsHartee = 2 * model.Lambda * sum( simdata.FS[t-1,t-1] )/disc.vol *  sum( simdata.FS[t-1,t-1] )
    return E_HiggsHartee
end

function getE_MixedHartree(t::Int64, model::SUNgaugeScalar, meas::MeasurementTwoPIGaugeScalar, pexp::TwoPIGaugeScalarLoopEx, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, disc::TwoPIGaugeScalarDiscretization)
    ### Higgs Hartree contribution
    E_MixedHartee = 0.5 * model.g^2 * sum( (disc.sdim-1)*simdata.FT[t-1,t-1] + simdata.FL[t-1,t-1])/disc.vol * sum( simdata.FS[t-1,t-1] )
    return E_MixedHartee
end

function getE_GaugeHartree(t::Int64, model::SUNgaugeScalar, meas::MeasurementTwoPIGaugeScalar, pexp::TwoPIGaugeScalarLoopEx, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, disc::TwoPIGaugeScalarDiscretization)
    ### Gauge Hartree contribution
    calcAngprdterm!(model, disc, simdata, tmpdata, t-1)
    
    transmassterm1 = model.N * model.g^2 * sum((disc.sdim-2 + (disc.sdim-1)^(-1)) .* simdata.FT[t-1,t-1] .+ (1 - (disc.sdim-1)^(-1)) .* simdata.FL[t-1,t-1])/disc.vol
    longimassterm1 = model.N * model.g^2 * sum((disc.sdim-2) .* simdata.FT[t-1,t-1] .+ simdata.FL[t-1,t-1])/disc.vol
    
    E_GaugeHartree =  0.25 * (disc.sdim-1) * sum(transmassterm1 * simdata.FT[t-1,t-1] )
    E_GaugeHartree += 0.25 *                 sum(longimassterm1 * simdata.FL[t-1,t-1] )

    transmassterm2 = (-model.N * model.g^2 * (disc.sdim-1)^(-1)) .* tmpdata.res
    longimassterm2 = (model.N * model.g^2) .* tmpdata.res
    E_GaugeHartree += 0.25 * (disc.sdim-1) * sum( transmassterm2 .* simdata.FT[t-1,t-1] )
    E_GaugeHartree += 0.25 *                 sum( longimassterm2 .* simdata.FL[t-1,t-1] )
    
    return E_GaugeHartree
end

function getE_MixedSunset(t::Int64, model::SUNgaugeScalar, meas::MeasurementTwoPIGaugeScalar, pexp::TwoPIGaugeScalarLoopEx, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, disc::TwoPIGaugeScalarDiscretization)
    ### Mixed sunset contribution
    if t == 2 # first measurement after initialisation, before evolve!
        E_MixedSunset = 0
    else # t == 3: first step in trapez integration
        #E_MixedSunset = 0.5 * 0.5 * disc.dt * sum( simdata.SigrS[simdata.indices[1]] .* simdata.FS[simdata.indices[1],t-1] - simdata.SigFS[simdata.indices[1]] .* simdata.rS[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) )
        E_MixedSunset = 0.5 * 0.25 * disc.dt *(disc.sdim-1)* sum( simdata.SigrT[simdata.indices[1]] .* simdata.FT[simdata.indices[1],t-1] - simdata.SigFT[simdata.indices[1]] .* simdata.rT[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) )
        E_MixedSunset += 0.5 * 0.25 * disc.dt * sum( simdata.SigrL[simdata.indices[1]] .* simdata.FL[simdata.indices[1],t-1] - simdata.SigFL[simdata.indices[1]] .* simdata.rL[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) )
        if t > 3 # two or more steps in trapez integration
            #E_MixedSunset +=  0.5 * disc.dt * sum( [ sum(simdata.SigrS[tp] .* simdata.FS[tp,t-1] - simdata.SigFS[tp] .* simdata.rS[tp,t-1]*thesign(tp, t-1, simdata.NstepsinMemory)) for tp in simdata.SigFS.indices[1]+1:(t-2) ] )
            E_MixedSunset += 0.25 * disc.dt * (disc.sdim-1)*sum([ sum( simdata.SigrT[simdata.indices[1]] .* simdata.FT[simdata.indices[1],t-1] - simdata.SigFT[simdata.indices[1]] .* simdata.rT[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) ) for tp in simdata.SigFT.indices[1]+1:(t-2) ] )
            E_MixedSunset += 0.25 * disc.dt *               sum([ sum( simdata.SigrL[simdata.indices[1]] .* simdata.FL[simdata.indices[1],t-1] - simdata.SigFL[simdata.indices[1]] .* simdata.rL[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) ) for tp in simdata.SigFL.indices[1]+1:(t-2) ] )
        end
    end
    return E_MixedSunset
end
function getE_GaugeSunset(t::Int64, model::SUNgaugeScalar, meas::MeasurementTwoPIGaugeScalar, pexp::TwoPIGaugeScalarLoopEx, simdata::TwoPIGaugeScalarSimData, tmpdata::TwoPIGaugeScalarTmpData, disc::TwoPIGaugeScalarDiscretization)
    ### Gauge sunset contribution
    if t == 2 # first measurement after initialisation, before evolve!
        E_GaugeSunset = 0
    else # t == 3: first step in trapez integration
        E_GaugeSunset = 0.5 * disc.dt * (1/3) *(disc.sdim-1)* sum( simdata.SigrTgaugesunset[simdata.indices[1]] .* simdata.FT[simdata.indices[1],t-1] - simdata.SigFTgaugesunset[simdata.indices[1]] .* simdata.rT[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) )
        E_GaugeSunset = 0.5 * disc.dt * (1/3) * sum( simdata.SigrLgaugesunset[simdata.indices[1]] .* simdata.FL[simdata.indices[1],t-1] - simdata.SigFLgaugesunset[simdata.indices[1]] .* simdata.rL[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) )
        if t > 3 # two or more steps in trapez integration
            E_GaugeSunset += disc.dt * (1/3) *(disc.sdim-1)* sum( [ sum(simdata.SigrTgaugesunset[simdata.indices[1]] .* simdata.FT[simdata.indices[1],t-1] - simdata.SigFTgaugesunset[simdata.indices[1]] .* simdata.rT[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory)) for tp in simdata.SigFT.indices[1]+1:(t-2) ] )
            E_GaugeSunset += disc.dt * (1/3) * sum( [ sum(simdata.SigrLgaugesunset[simdata.indices[1]] .* simdata.FL[simdata.indices[1],t-1] - simdata.SigFLgaugesunset[simdata.indices[1]] .* simdata.rL[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory)) for tp in simdata.SigFL.indices[1]+1:(t-2) ] )
        end
    end
    return E_GaugeSunset
end


export measure!
function measure!(thesolution::QFTdynamicsSolutionTwoPIGaugeScalar, tmpdata::TwoPIGaugeScalarTmpData, t)
    @unpack problem, simdata, measurearray = thesolution
    @unpack model, pexp, disc, init, reno, simsetup = problem
    measurearray[Int64((t-2)/simsetup.NstepsbetweenNmeas) + 1] = createMeasurement(t, model, pexp, simdata, tmpdata, disc)
end
