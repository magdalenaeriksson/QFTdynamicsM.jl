export MeasurementCSScalar
Base.@kwdef mutable struct MeasurementCSScalar{N} <: Measurement
	time::Float64=0
    # from x space
	E::Float64=0
	E_phi::Float64=0
	E_pii::Float64=0
	# from k space
	phi2::Vector{Float64}=zeros(Float64, N )
	phi2_err::Vector{Float64}=zeros(Float64, N )
	phi4::Vector{Float64}=zeros(Float64, N )
	pii2::Vector{Float64}=zeros(Float64, N )
	pii2_err::Vector{Float64}=zeros(Float64, N )
	pii4::Vector{Float64}=zeros(Float64, N )
	n::Vector{Float64}=zeros(Float64, N )
	n_err::Vector{Float64}=zeros(Float64, N )
	omega::Vector{Float64}=zeros(Float64, N )
	omega_err::Vector{Float64}=zeros(Float64, N )
end

export MeasurementCSScalarofRun
Base.@kwdef mutable struct MeasurementCSScalarofRun{N} <: Measurement
	time::Float64=0
    # from x space
	E::Float64=0
	E_phi::Float64=0
	E_pii::Float64=0
	# from k space
	phi2::Vector{Float64}=zeros(Float64, N )
	phi4::Vector{Float64}=zeros(Float64, N )
	pii2::Vector{Float64}=zeros(Float64, N )
	pii4::Vector{Float64}=zeros(Float64, N )
end

import Base: +
function +(A::MeasurementCSScalarofRun,B::MeasurementCSScalarofRun)
	return MeasurementCSScalarofRun{length(A.phi2)}(	A.time	,
					A.E     	+B.E     	,
					A.E_phi   	+B.E_phi    ,
					A.E_pii   	+B.E_pii   	,
					A.phi2  	+B.phi2  	,
					A.phi4  	+B.phi4  	,
					A.pii2  	+B.pii2  	,
					A.pii4  	+B.pii4  	, )
end

function combinemeasurements(model::CSScalarPhi4, MeasofRuns::Vector{Measurement}, num::CSScalarCPU, disc::CSScalarDiscretization)
	# sum all measurments of runs
	tmpmeas = MeasurementCSScalarofRun{length(disc.fftwhelper)}()
	for i in 1:num.Runs
		tmpmeas += MeasofRuns[i]
	end
	# divide by number of runs
	totalmeas = MeasurementCSScalar{length(disc.fftwhelper)}()
	totalmeas.time = MeasofRuns[1].time
	totalmeas.E = tmpmeas.E / num.Runs
	totalmeas.E_phi = tmpmeas.E_phi / num.Runs
	totalmeas.E_pii = tmpmeas.E_pii / num.Runs
	totalmeas.phi2 .= tmpmeas.phi2 / num.Runs
	totalmeas.phi4 .= tmpmeas.phi4 / num.Runs
	totalmeas.pii2 .= tmpmeas.pii2 / num.Runs
	totalmeas.pii4 .= tmpmeas.pii4 / num.Runs

	# create new things
	totalmeas.phi2_err .= sqrt.( totalmeas.phi4 .- totalmeas.phi2.^2 ) ./ (num.Runs*disc.deg)
	totalmeas.pii2_err .= sqrt.( totalmeas.pii4 .- totalmeas.pii2.^2 ) ./ (num.Runs*disc.deg)
	totalmeas.n 	   .= sqrt.( totalmeas.phi2 .* totalmeas.pii2 ) .- 0.5
	totalmeas.n_err    .= sqrt.( 0.25 * (totalmeas.pii2 ./ totalmeas.phi2) .* totalmeas.phi2_err.^2 + 0.25 * (totalmeas.phi2 ./ totalmeas.phi2) .* totalmeas.pii2_err.^2)
	totalmeas.omega	   .= sqrt.( totalmeas.pii2 ./ totalmeas.phi2 )
	totalmeas.omega_err.= sqrt.( 0.25 *  (totalmeas.omega ./ totalmeas.phi2) .* totalmeas.phi2_err.^2 + 0.25 * ( totalmeas.omega ./ totalmeas.pii2)  .* totalmeas.pii2_err.^2) 

	return totalmeas
end

export createMeasurementofRun
function createMeasurementofRun(t::Int64, model::CSScalarPhi4, tmpdata::CSScalarTmpDataCPU, disc::CSScalarDiscretization)
	# init measurement object
	meas = MeasurementCSScalarofRun{length(disc.fftwhelper)}()
    ## scalar Measurment
	meas.time = t * disc.dt * disc.Mass
    # from x space - Energy
   	meas.E_phi = getEphi(model, tmpdata)
  	meas.E_pii = getEpii(model, tmpdata)
   	meas.E = meas.E_phi + meas.E_pii
	# from k space
	for j in 1:model.ONgroup
		tmpmeas = MeasurementCSScalarofRun{length(disc.fftwhelper)}()
		# fft space stuff
		# inplace FFTW	
		#fft!(tmpdata.phi_cpu[j]) # creates asym one
        tmpdata.ftplan * tmpdata.phi_cpu[j] # creates asym one
		#fft!(tmpdata.pii_cpu[j]) # creates asym one
        tmpdata.ftplan * tmpdata.pii_cpu[j] # creates asym one
		norm = sqrt(disc.Nx^disc.sdim)
		tmpdata.phi_cpu[j] ./= norm # now we have sym
		tmpdata.pii_cpu[j] ./= norm # now we have sym
		for i in 1:length(disc.fftwhelper)
    		for deg in 1:disc.fftwhelper[i].deg
				index = disc.fftwhelper[i].ind[deg]
				# phi2
				tmpmeas.phi2[i]		+= real(tmpdata.phi_cpu[j][index])^2 + imag(tmpdata.phi_cpu[j][index])^2
				tmpmeas.phi4[i]		+= real(tmpdata.phi_cpu[j][index])^4 + imag(tmpdata.phi_cpu[j][index])^4 + 2*real(tmpdata.phi_cpu[j][index])^2 * imag(tmpdata.phi_cpu[j][index])^2
				# pii2
				tmpmeas.pii2[i]		+= real(tmpdata.pii_cpu[j][index])^2 + imag(tmpdata.pii_cpu[j][index])^2
				tmpmeas.pii4[i]		+= real(tmpdata.pii_cpu[j][index])^4 + imag(tmpdata.pii_cpu[j][index])^4 + 2*real(tmpdata.pii_cpu[j][index])^2 * imag(tmpdata.pii_cpu[j][index])^2
			end
			#phi2
			tmpmeas.phi2[i]		/= disc.fftwhelper[i].deg
			tmpmeas.phi4[i]		/= disc.fftwhelper[i].deg
			#pii2
			tmpmeas.pii2[i]		/= disc.fftwhelper[i].deg
			tmpmeas.pii4[i]		/= disc.fftwhelper[i].deg
		end
		# add each ON component up	
		meas += tmpmeas 
	end
	# divide by ON to get average 
	meas.phi2 ./= model.ONgroup
	meas.phi4 ./= model.ONgroup
	meas.pii2 ./= model.ONgroup
	meas.pii4 ./= model.ONgroup

   return meas
end

#
# Energy
#
function getEphi(model::CSScalarPhi4, tmpdata::CSScalarTmpDataCPU)
	Ephi = 0
	for j in 1:model.ONgroup
		# mass term
		Ephi += (0.5*model.Mass2*sum(real(tmpdata.phi_cpu[j]).^2))
	end
	# kinetic term
	Ephi += getEkin(model, tmpdata.phi_cpu)
	# interactionterm
	Ephi += model.Lambda/(24 * model.ONgroup) * sum(sum([real(tmpdata.phi_cpu[j]).^2 for j in 1:model.ONgroup]).^2)
	return Ephi
end

function getEkin(model::CSScalarPhi4, phi::Vector{Array{ComplexF64, 1}})
	Ekin = 0
	for j in 1:model.ONgroup
		Ekin += 0.5 * sum(real((circshift(phi[j], 1) - phi[j])).*real((circshift(phi[j], 1) - phi[j])))
	end
	return Ekin
end

function getEkin(model::CSScalarPhi4, phi::Vector{Array{ComplexF64, 2}})
	Ekin = 0
	for j in 1:model.ONgroup
		Ekin += 0.5 * sum(real((circshift(phi[j], (1,0)) - phi[j])).*real((circshift(phi[j], (1,0)) - phi[j])))
		Ekin += 0.5 * sum(real((circshift(phi[j], (0,1)) - phi[j])).*real((circshift(phi[j], (0,1)) - phi[j])))
	end
	return Ekin
end

function getEkin(model::CSScalarPhi4, phi::Vector{Array{ComplexF64, 3}})
	Ekin = 0
	for j in 1:model.ONgroup
		Ekin += 0.5 * sum(real((circshift(phi[j], (1,0,0)) - phi[j])).*real((circshift(phi[j], (1,0,0)) - phi[j])))
		Ekin += 0.5 * sum(real((circshift(phi[j], (0,1,0)) - phi[j])).*real((circshift(phi[j], (0,1,0)) - phi[j])))
		Ekin += 0.5 * sum(real((circshift(phi[j], (0,0,1)) - phi[j])).*real((circshift(phi[j], (0,0,1)) - phi[j])))
	end
	return Ekin
end

function getEpii(model::CSScalarPhi4, tmpdata::CSScalarTmpDataCPU)
	Epii = 0
	for j in 1:model.ONgroup
		Epii += 0.5* sum(real(tmpdata.pii_cpu[j]).*real(tmpdata.pii_cpu[j]))
	end
	return Epii
end

export measure!
function measure!( thesolution::QFTdynamicsSolutionCSScalar, tmpdata::Vector{CSScalarTmpDataCPU}, t)
    @unpack problem, simdata, measurearray, measurearrayofruns = thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem

	#MeasofRuns = Vector{MeasurementCSScalarofRun}(undef, num.Runs)
    @Threads.threads for ichunk in 1:num.threads
        for i in num.threadranges[ichunk]
			copy_toTmpdata!(model, tmpdata[ichunk].pii_cpu, simdata[i].pii)
			copy_toTmpdata!(model, tmpdata[ichunk].phi_cpu, simdata[i].phi)
			#MeasofRuns[i] = createMeasurementofRun(t, model, tmpdata[ichunk], disc)
			measurearrayofruns[Int64((t-2)/simsetup.NstepsbetweenNmeas) + 1][i]= createMeasurementofRun(t, model, tmpdata[ichunk], disc)
		end
	end
	# combine here the Measurments of the individual runs
    #measurearray[Int64((t-2)/simsetup.NstepsbetweenNmeas) + 1] = combinemeasurements(model, MeasofRuns, num, disc)
    measurearray[Int64((t-2)/simsetup.NstepsbetweenNmeas) + 1] = combinemeasurements(model, measurearrayofruns[Int64((t-2)/simsetup.NstepsbetweenNmeas) + 1], num, disc)
end