export MeasurementGenericModel
export measure!

Base.@kwdef mutable struct MeasurementGenericModel{N,Memt} <: Measurement
	time::Float64=0
	quantity::Float64=0
end

function createMeasurement(t::Int64, model::GenericModelA, pexp::anexpansion, simdata::GenericModelSimData, disc::DiscretizationGenericModel)
	# init measurement object
	#meas = MeasurementGenericModel{length(disc.fftwhelper),disc.Memt}()
	meas = MeasurementGenericModel{1,disc.NstepsinMemory}()

    # scalar Measurment
	meas.time = t * disc.dt
	meas.quantity = sum(simdata.quantitytoevolve[t]) 

    return meas
end

function measure!( thesolution::QFTdynamicsSolutionGenericModel, t)
    @unpack problem, simdata, measurearray = thesolution
    @unpack model, pexp, disc, init, reno, simsetup = problem

    measurearray[Int64((t-2)/simsetup.NstepsbetweenNmeas) + 1] = createMeasurement(t, model, pexp, simdata, disc)
end
