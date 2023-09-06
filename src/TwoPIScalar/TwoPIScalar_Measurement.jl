export MeasurementTwoPIScalar
Base.@kwdef mutable struct MeasurementTwoPIScalar{N,NstepsinMemory} <: Measurement
	time::Float64=0
	E::Float64=0
	E_pot::Float64=0
	E_kin::Float64=0
	E_mem::Float64=0
	M2::Float64=0
	# from k space
	phi::Float64=0
	F::Vector=zeros(Float64, N )
	ddF::Vector=zeros(Float64, N )
	r::Vector=zeros(Float64, N )
	n::Vector=zeros(Float64, N )
	omega::Vector=zeros(Float64, N )
	c::Vector=zeros(Float64, N )
    # offdiagonal entries
    Fonet::Vector=zeros(Float64, NstepsinMemory )
    ronet::Vector=zeros(Float64, NstepsinMemory )
    # Sigma
	SigF::Vector=zeros(Float64, NstepsinMemory )
	Sigr::Vector=zeros(Float64, NstepsinMemory )
end

function getmeasureidx(lev,simdata::TwoPIScalarSimDataCPUfull, disc::TwoPIScalarDiscretization) return disc.fftwhelper[lev].ind[1] end # take deg nr 1
function getmeasureidx(lev,simdata::TwoPIScalarSimDataCPUred, disc::TwoPIScalarDiscretization) return lev end
function getmeasureidx(lev,simdata::TwoPIScalarSimDataCPUred2, disc::TwoPIScalarDiscretization) return lev end
function getmeasureidx(lev,simdata::TwoPIScalarSimDataCPUcont, disc::TwoPIScalarDiscretization) return lev end

export createMeasurement
function createMeasurement(t::Int64, model::TwoPIScalarPhi4, simdata::TwoPIScalarSimDataCPU, disc::TwoPIScalarDiscretization, pexp::TwoPIScalarPertExpansion)
	# init measurement object
	meas = MeasurementTwoPIScalar{length(disc.fftwhelper),disc.NstepsinMemory}()
    # scalar Measurment
	meas.time = (t-1) * disc.dt * disc.Mass
	meas.M2 = simdata.hartreemass2[1]
    # vector Measurement
	for i in 1:length(disc.fftwhelper)
        idx = getmeasureidx(i,simdata, disc)
        meas.ddF[i] = (simdata.F[t,t,idx] - 2*simdata.F[t-1,t,idx] + simdata.F[t-1,t-1,idx] )/disc.dt^2
        meas.F[i] = simdata.F[t-1,t-1,idx] 
        meas.r[i] = simdata.r[t-1,t-1,idx] * thesign(t-1, t-1, simdata.NstepsinMemory)
        if (meas.ddF[i] / simdata.F[t-1,t-1,idx]) > 0
            meas.omega[i] = sqrt( meas.ddF[i] / simdata.F[t-1,t-1,idx])
        else
            meas.omega[i] = 0
        end
        if 0.25*disc.dt^2*meas.omega[i]^2 < 1
            meas.c[i] = sqrt( 1. - 0.25*disc.dt^2*meas.omega[i]^2)
        else
            meas.c[i] = 1
            print("E")
        end
        if meas.ddF[i] * simdata.F[t-1,t-1,idx] > 0
            meas.n[i] = (meas.c[i] * sqrt( meas.ddF[i] * simdata.F[t-1,t-1,idx])) - 0.5
        else 
            meas.n[i] =  0
        end
    end
    # Energy - depends on vector Measurement
    meas.E_kin = getE_kin(t, model, meas, simdata, disc, pexp)
    meas.E_pot = getE_pot(t, model, meas, simdata, disc, pexp)
    meas.E_mem = getE_mem(t, model, meas, simdata, disc, pexp)
    meas.E = meas.E_kin + meas.E_pot + meas.E_mem
    #hartreecont= 0
    #meas.E_kin = 0
    #meas.E_pot = 0
    #meas.E_mem = 0
	#for i in 1:length(disc.fftwhelper)
    #    idx = disc.fftwhelper[i].ind[1] # take deg nr 1
    #    hartreecont += simdata.F[t-1,t-1,idx] * disc.fftwhelper[i].deg 
    #end
	#for i in 1:length(disc.fftwhelper)
    #    idx = disc.fftwhelper[i].ind[1] # take deg nr 1
    #    meas.E_kin += 0.5 * meas.ddF[i] * disc.fftwhelper[i].deg # E_pi - count every degeneracy
    #    meas.E_pot += 0.5 * (model.Mass2 + disc.fftwhelper[i].lev2 + 0.5 * disc.ivol * Hartreeprefactor(model, pexp) * hartreecont ) * simdata.F[t-1,t-1,idx] * disc.fftwhelper[i].deg # E_phi - count every degeneracy
    #end
    ## trapezintegration: boundary
    #if t==2 # first measurment after initialization, before evolution
    #    meas.E_mem = 0
    #else
    #    meas.E_mem = 0.5 *0.25 * disc.dt * sum( simdata.Sigr[simdata.indices[1]].*simdata.F[simdata.indices[1],t-1] .- simdata.SigF[simdata.indices[1]].*simdata.r[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) )
    #    if t>3 # first case where time integration has more that two points is t=4 
    #    # trapezintegration: rest
    #        meas.E_mem += 0.25 * disc.dt * sum( [ sum(simdata.Sigr[tp].*simdata.F[tp,t-1] .- simdata.SigF[tp].*simdata.r[tp,t-1]*thesign(tp, t-1, simdata.NstepsinMemory)) for tp in simdata.SigF.indices[1]+1:(t-2) ] )
    #    end
    #end
    #meas.E = meas.E_kin + meas.E_pot + meas.E_mem
    # Sig in Memory past
    whichindex = 1
    meas.SigF = [ simdata.SigF[tp][whichindex] for tp in (simdata.SigF.indices[2]-1):-1:simdata.SigF.indices[1]] 
    meas.Sigr = [ simdata.Sigr[tp][whichindex] for tp in (simdata.Sigr.indices[2]-1):-1:simdata.Sigr.indices[1]] 
    # offdiagonal
    # correlation with oldest one
    meas.Fonet = [ simdata.F[simdata.indices[1],tp,whichindex] for tp in simdata.indices[1]:(simdata.indices[2]-1)]
    meas.ronet = [ simdata.r[simdata.indices[1],tp,whichindex]*thesign(simdata.indices[1], tp, simdata.NstepsinMemory) for tp in simdata.indices[1]:(simdata.indices[2]-1)]

    return meas
end

function getE_pot(t::Int64, model::TwoPIScalarPhi4, meas::MeasurementTwoPIScalar, simdata::TwoPIScalarSimDataCPU, disc::TwoPIScalarDiscretizationCont, pexp::TwoPIScalarPertExpansion)
    E_pot = 0
	for i in 1:length(disc.fftwhelper)
        #idx = disc.fftwhelper[i].ind[1] # take deg nr 1
        idx = getmeasureidx(i,simdata, disc)
        E_pot += 0.5 * (model.Mass2 + disc.fftwhelper[i].lev2 + 0.5 * Hartreeprefactor(model, pexp) * simdata.Fx[t-1, t-1][1] ) * simdata.F[t-1,t-1,idx] * (disc.fftwhelper[i].lev2 * 4*pi^2 / disc.Nmom) #last bit is FD!
    end
    return E_pot
end
function getE_kin(t::Int64, model::TwoPIScalarPhi4, meas::MeasurementTwoPIScalar, simdata::TwoPIScalarSimDataCPU, disc::TwoPIScalarDiscretizationCont, pexp::TwoPIScalarPertExpansion)
    E_kin = 0
	for i in 1:length(disc.fftwhelper)
        #idx = disc.fftwhelper[i].ind[1] # take deg nr 1
        idx = getmeasureidx(i,simdata, disc)
        E_kin += 0.5 * meas.ddF[i] *  (disc.fftwhelper[i].lev2 * 4*pi^2 / disc.Nmom)
    end
    return E_kin
end
function getE_mem(t::Int64, model::TwoPIScalarPhi4, meas::MeasurementTwoPIScalar, simdata::TwoPIScalarSimDataCPU, disc::TwoPIScalarDiscretizationCont, pexp::TwoPIScalarPertExpansion)
    E_mem = 0
    momenta2 = [a.lev2 for a in disc.fftwhelper]
    # trapezintegration: boundary
    if t==2 # first measurment after initialization, before evolution
        E_mem = 0
    else
        E_mem = 0.5 *0.25 * disc.dt * sum( (simdata.Sigr[simdata.indices[1]].*simdata.F[simdata.indices[1],t-1] .- simdata.SigF[simdata.indices[1]].*simdata.r[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory)) .* (momenta2 * 4*pi^2 / disc.Nmom))
        if t>3 # first case where time integration has more that two points is t=4 
        # trapezintegration: rest
            E_mem += 0.25 * disc.dt * sum( [ sum( (simdata.Sigr[tp].*simdata.F[tp,t-1] .- simdata.SigF[tp].*simdata.r[tp,t-1]*thesign(tp, t-1, simdata.NstepsinMemory) ) .*(momenta2 * 4*pi^2 / disc.Nmom) ) for tp in simdata.SigF.indices[1]+1:(t-2) ] )
        end
    end
    return E_mem
end

function getE_pot(t::Int64, model::TwoPIScalarPhi4, meas::MeasurementTwoPIScalar, simdata::TwoPIScalarSimDataCPU, disc::TwoPIScalarDiscretizationLattice, pexp::TwoPIScalarPertExpansion)
    E_pot = 0
    hartreecont= 0
	for i in 1:length(disc.fftwhelper)
        #idx = disc.fftwhelper[i].ind[1] # take deg nr 1
        idx = getmeasureidx(i,simdata, disc)
        hartreecont += simdata.F[t-1,t-1,idx] * disc.fftwhelper[i].deg 
    end
	for i in 1:length(disc.fftwhelper)
        #idx = disc.fftwhelper[i].ind[1] # take deg nr 1
        idx = getmeasureidx(i,simdata, disc)
        E_pot += 0.5 * (model.Mass2 + disc.fftwhelper[i].lev2 + 0.5 * disc.ivol * Hartreeprefactor(model, pexp) * hartreecont ) * simdata.F[t-1,t-1,idx] * disc.fftwhelper[i].deg # E_phi - count every degeneracy
    end
    return E_pot
end

function getE_kin(t::Int64, model::TwoPIScalarPhi4, meas::MeasurementTwoPIScalar, simdata::TwoPIScalarSimDataCPU, disc::TwoPIScalarDiscretizationLattice, pexp::TwoPIScalarPertExpansion)
    E_kin = 0
	for i in 1:length(disc.fftwhelper)
        #idx = disc.fftwhelper[i].ind[1] # take deg nr 1
        idx = getmeasureidx(i,simdata, disc)
        E_kin += 0.5 * meas.ddF[i] * disc.fftwhelper[i].deg # E_pi - count every degeneracy
    end
    return E_kin
end

function getE_mem(t::Int64, model::TwoPIScalarPhi4, meas::MeasurementTwoPIScalar, simdata::TwoPIScalarSimDataCPUfull, disc::TwoPIScalarDiscretizationLattice, pexp::TwoPIScalarPertExpansion)
    E_mem = 0
    # trapezintegration: boundary
    if t==2 # first measurment after initialization, before evolution
        E_mem = 0
    else
        E_mem = 0.5 *0.25 * disc.dt * sum( simdata.Sigr[simdata.indices[1]].*simdata.F[simdata.indices[1],t-1] .- simdata.SigF[simdata.indices[1]].*simdata.r[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) )
        if t>3 # first case where time integration has more that two points is t=4 
        # trapezintegration: rest
            E_mem += 0.25 * disc.dt * sum( [ sum(simdata.Sigr[tp].*simdata.F[tp,t-1] .- simdata.SigF[tp].*simdata.r[tp,t-1]*thesign(tp, t-1, simdata.NstepsinMemory)) for tp in simdata.SigF.indices[1]+1:(t-2) ] )
        end
    end
    return E_mem
end

function getE_mem(t::Int64, model::TwoPIScalarPhi4, meas::MeasurementTwoPIScalar, simdata::TwoPIScalarSimDataCPUred, disc::TwoPIScalarDiscretizationLattice, pexp::TwoPIScalarPertExpansion)
    E_mem = 0
    # trapezintegration: boundary
    if t==2 # first measurment after initialization, before evolution
        E_mem = 0
    else
	    for i in 1:length(disc.fftwhelper)
            E_mem += 0.5 *0.25 * disc.dt * disc.fftwhelper[i].deg * ( simdata.Sigr[simdata.indices[1]][i]*simdata.F[simdata.indices[1],t-1][i] - simdata.SigF[simdata.indices[1]][i]*simdata.r[simdata.indices[1],t-1][i]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) )
        end
        #E_mem = 0.5 *0.25 * disc.dt * sum( simdata.Sigr[simdata.indices[1]].*simdata.F[simdata.indices[1],t-1] .- simdata.SigF[simdata.indices[1]].*simdata.r[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) )
        if t>3 # first case where time integration has more that two points is t=4 
        # trapezintegration: rest
	        for i in 1:length(disc.fftwhelper)
                for tp in simdata.SigF.indices[1]+1:(t-2)
                    E_mem += 0.25 * disc.dt * disc.fftwhelper[i].deg * (simdata.Sigr[tp][i] * simdata.F[tp,t-1][i] - simdata.SigF[tp][i] * simdata.r[tp,t-1][i]*thesign(tp, t-1, simdata.NstepsinMemory))
                end
            end
            #E_mem += 0.25 * disc.dt * sum( [ sum(simdata.Sigr[tp].*simdata.F[tp,t-1] .- simdata.SigF[tp].*simdata.r[tp,t-1]*thesign(tp, t-1, simdata.NstepsinMemory)) for tp in simdata.SigF.indices[1]+1:(t-2) ] )
        end
    end
    return E_mem
end

function getE_mem(t::Int64, model::TwoPIScalarPhi4, meas::MeasurementTwoPIScalar, simdata::TwoPIScalarSimDataCPUred2, disc::TwoPIScalarDiscretizationLattice, pexp::TwoPIScalarPertExpansion)
    E_mem = 0
    # trapezintegration: boundary
    if t==2 # first measurment after initialization, before evolution
        E_mem = 0
    else
	    for i in 1:length(disc.fftwhelper)
            E_mem += 0.5 *0.25 * disc.dt * disc.fftwhelper[i].deg * ( simdata.Sigr[simdata.indices[1]][i]*simdata.F[simdata.indices[1],t-1][i] - simdata.SigF[simdata.indices[1]][i]*simdata.r[simdata.indices[1],t-1][i]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) )
        end
        #E_mem = 0.5 *0.25 * disc.dt * sum( simdata.Sigr[simdata.indices[1]].*simdata.F[simdata.indices[1],t-1] .- simdata.SigF[simdata.indices[1]].*simdata.r[simdata.indices[1],t-1]*thesign(simdata.indices[1], t-1, simdata.NstepsinMemory) )
        if t>3 # first case where time integration has more that two points is t=4 
        # trapezintegration: rest
	        for i in 1:length(disc.fftwhelper)
                for tp in simdata.SigF.indices[1]+1:(t-2)
                    E_mem += 0.25 * disc.dt * disc.fftwhelper[i].deg * (simdata.Sigr[tp][i] * simdata.F[tp,t-1][i] - simdata.SigF[tp][i] * simdata.r[tp,t-1][i]*thesign(tp, t-1, simdata.NstepsinMemory))
                end
            end
            #E_mem += 0.25 * disc.dt * sum( [ sum(simdata.Sigr[tp].*simdata.F[tp,t-1] .- simdata.SigF[tp].*simdata.r[tp,t-1]*thesign(tp, t-1, simdata.NstepsinMemory)) for tp in simdata.SigF.indices[1]+1:(t-2) ] )
        end
    end
    return E_mem
end

export measure!
function measure!( thesolution::QFTdynamicsSolutionTwoPIScalar, t)
    @unpack problem, simdata, measurearray =thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem

    measurearray[Int64((t-2)/simsetup.NstepsbetweenNmeas) + 1] = createMeasurement(t, model, simdata, disc, pexp)
end
