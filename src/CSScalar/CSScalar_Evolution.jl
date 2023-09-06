using Parameters
using FFTW
using ParallelStencil
@static if parameters["sdim"] == 1
    using ParallelStencil.FiniteDifferences1D
    @init_parallel_stencil(Threads, Float64, 1)
end
@static if parameters["sdim"] == 2
    using ParallelStencil.FiniteDifferences3D
    @init_parallel_stencil(Threads, Float64, 2)
end
@static if parameters["sdim"] == 3
    using ParallelStencil.FiniteDifferences3D
    @init_parallel_stencil(Threads, Float64, 3)
end

#
# stencil stuff
#
@static if parameters["sdim"] == 1
    #N=1
    @parallel function evolve_pii!(phi::Array{Float64, 1},pii::Array{Float64, 1}, update, dt, Mass2, Lambda) @inn(pii) = @inn(pii) + dt * (@d2(phi) - (Mass2 + Lambda/6. * @inn(phi).*@inn(phi)) * @inn(phi))
    return end
    #N=2
    @parallel function evolve_pii!(phi_1::Array{Float64, 1},phi_2::Array{Float64, 1},pii_1::Array{Float64, 1},pii_2::Array{Float64, 1}, update, dt, Mass2, Lambda)
	@inn(update) = @inn(phi_1)*@inn(phi_1) + @inn(phi_2)*@inn(phi_2)
	@inn(pii_1) = @inn(pii_1) + dt * (@d2(phi_1) - Mass2 * @inn(phi_1) - Lambda/12. * @inn(update) .* @inn(phi_1) )
	@inn(pii_2) = @inn(pii_2) + dt * (@d2(phi_2) - Mass2 * @inn(phi_2) - Lambda/12. * @inn(update) .* @inn(phi_2) )
    return end
    #N=4
    @parallel function evolve_pii!(phi_1::Array{Float64, 1},phi_2::Array{Float64, 1},phi_3::Array{Float64, 1},phi_4::Array{Float64, 1},pii_1::Array{Float64, 1},pii_2::Array{Float64, 1},pii_3::Array{Float64, 1},pii_4::Array{Float64, 1}, update, dt, Mass2, Lambda)
	@inn(update) = @inn(phi_1).*@inn(phi_1) + @inn(phi_2).*@inn(phi_2) + @inn(phi_3).*@inn(phi_3) + @inn(phi_4).*@inn(phi_4) 
	@inn(pii_1) = @inn(pii_1) + dt * (@d2(phi_1) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_1))
	@inn(pii_2) = @inn(pii_2) + dt * (@d2(phi_2) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_2))
	@inn(pii_3) = @inn(pii_3) + dt * (@d2(phi_3) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_3))
	@inn(pii_4) = @inn(pii_4) + dt * (@d2(phi_4) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_4))
    return end
    #N=8
    @parallel function evolve_pii!(phi_1::Array{Float64, 1},phi_2::Array{Float64, 1},phi_3::Array{Float64, 1},phi_4::Array{Float64, 1},phi_5::Array{Float64, 1},phi_6::Array{Float64, 1},phi_7::Array{Float64, 1},phi_8::Array{Float64, 1},pii_1::Array{Float64, 1},pii_2::Array{Float64, 1},pii_3::Array{Float64, 1},pii_4::Array{Float64, 1},pii_5::Array{Float64, 1},pii_6::Array{Float64, 1},pii_7::Array{Float64, 1},pii_8::Array{Float64, 1}, update, dt, Mass2, Lambda)
	@inn(update) = @inn(phi_1).*@inn(phi_1) + @inn(phi_2).*@inn(phi_2) + @inn(phi_3).*@inn(phi_3) + @inn(phi_4).*@inn(phi_4) + @inn(phi_5).*@inn(phi_5) + @inn(phi_6).*@inn(phi_6) + @inn(phi_7).*@inn(phi_7) + @inn(phi_8).*@inn(phi_8) 
	@inn(pii_1) = @inn(pii_1) + dt * (@d2(phi_1) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_1))
	@inn(pii_2) = @inn(pii_2) + dt * (@d2(phi_2) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_2))
	@inn(pii_3) = @inn(pii_3) + dt * (@d2(phi_3) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_3))
	@inn(pii_4) = @inn(pii_4) + dt * (@d2(phi_4) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_4))
	@inn(pii_5) = @inn(pii_5) + dt * (@d2(phi_5) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_5))
	@inn(pii_6) = @inn(pii_6) + dt * (@d2(phi_6) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_6))
	@inn(pii_7) = @inn(pii_7) + dt * (@d2(phi_7) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_7))
	@inn(pii_8) = @inn(pii_8) + dt * (@d2(phi_8) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_8))
    return end
end
@static if parameters["sdim"] == 2
    #N=1
    @parallel function evolve_pii!(phi::Array{Float64, 2},pii::Array{Float64, 2}, update, dt, Mass2, Lambda) @inn(pii) = @inn(pii) + dt * (@d2_xi(phi) + @d2_yi(phi) - (Mass2 + Lambda/6. * @inn(phi).*@inn(phi)) * @inn(phi))
    return end
    #N=2
    @parallel function evolve_pii!(phi_1::Array{Float64, 2},phi_2::Array{Float64, 2},pii_1::Array{Float64, 2},pii_2::Array{Float64, 2}, update, dt, Mass2, Lambda)
	@inn(update) = @inn(phi_1)*@inn(phi_1) + @inn(phi_2)*@inn(phi_2)
	@inn(pii_1) = @inn(pii_1) + dt * (@d2_xi(phi_1) + @d2_yi(phi_1) - Mass2 * @inn(phi_1) - Lambda/12. * @inn(update) .* @inn(phi_1) )
	@inn(pii_2) = @inn(pii_2) + dt * (@d2_xi(phi_2) + @d2_yi(phi_2) - Mass2 * @inn(phi_2) - Lambda/12. * @inn(update) .* @inn(phi_2) )
    return end
    #N=4
    @parallel function evolve_pii!(phi_1::Array{Float64, 2},phi_2::Array{Float64, 2},phi_3::Array{Float64, 2},phi_4::Array{Float64, 2},pii_1::Array{Float64, 2},pii_2::Array{Float64, 2},pii_3::Array{Float64, 2},pii_4::Array{Float64, 2}, update, dt, Mass2, Lambda)
	@inn(update) = @inn(phi_1).*@inn(phi_1) + @inn(phi_2).*@inn(phi_2) + @inn(phi_3).*@inn(phi_3) + @inn(phi_4).*@inn(phi_4) 
	@inn(pii_1) = @inn(pii_1) + dt * (@d2_xi(phi_1) + @d2_yi(phi_1) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_1))
	@inn(pii_2) = @inn(pii_2) + dt * (@d2_xi(phi_2) + @d2_yi(phi_2) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_2))
	@inn(pii_3) = @inn(pii_3) + dt * (@d2_xi(phi_3) + @d2_yi(phi_3) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_3))
	@inn(pii_4) = @inn(pii_4) + dt * (@d2_xi(phi_4) + @d2_yi(phi_4) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_4))
    return end
    #N=8
    @parallel function evolve_pii!(phi_1::Array{Float64, 2},phi_2::Array{Float64, 2},phi_3::Array{Float64, 2},phi_4::Array{Float64, 2},phi_5::Array{Float64, 2},phi_6::Array{Float64, 2},phi_7::Array{Float64, 2},phi_8::Array{Float64, 2},pii_1::Array{Float64, 2},pii_2::Array{Float64, 2},pii_3::Array{Float64, 2},pii_4::Array{Float64, 2},pii_5::Array{Float64, 2},pii_6::Array{Float64, 2},pii_7::Array{Float64, 2},pii_8::Array{Float64, 2}, update, dt, Mass2, Lambda)
	@inn(update) = @inn(phi_1).*@inn(phi_1) + @inn(phi_2).*@inn(phi_2) + @inn(phi_3).*@inn(phi_3) + @inn(phi_4).*@inn(phi_4) + @inn(phi_5).*@inn(phi_5) + @inn(phi_6).*@inn(phi_6) + @inn(phi_7).*@inn(phi_7) + @inn(phi_8).*@inn(phi_8) 
	@inn(pii_1) = @inn(pii_1) + dt * (@d2_xi(phi_1) + @d2_yi(phi_1) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_1))
	@inn(pii_2) = @inn(pii_2) + dt * (@d2_xi(phi_2) + @d2_yi(phi_2) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_2))
	@inn(pii_3) = @inn(pii_3) + dt * (@d2_xi(phi_3) + @d2_yi(phi_3) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_3))
	@inn(pii_4) = @inn(pii_4) + dt * (@d2_xi(phi_4) + @d2_yi(phi_4) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_4))
	@inn(pii_5) = @inn(pii_5) + dt * (@d2_xi(phi_5) + @d2_yi(phi_5) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_5))
	@inn(pii_6) = @inn(pii_6) + dt * (@d2_xi(phi_6) + @d2_yi(phi_6) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_6))
	@inn(pii_7) = @inn(pii_7) + dt * (@d2_xi(phi_7) + @d2_yi(phi_7) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_7))
	@inn(pii_8) = @inn(pii_8) + dt * (@d2_xi(phi_8) + @d2_yi(phi_8) - (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_8))
    return end
end
@static if parameters["sdim"] == 3
    #N=1
    @parallel function evolve_pii!(phi::Array{Float64, 3},pii::Array{Float64, 3}, update, dt, Mass2, Lambda) @inn(pii) = @inn(pii) + dt * (@d2_xi(phi) + @d2_yi(phi) + @d2_zi(phi) - (Mass2 + Lambda/6. * @inn(phi).*@inn(phi)) * @inn(phi))
    return end
    #N=2
    @parallel function evolve_pii!(phi_1::Array{Float64, 3},phi_2::Array{Float64, 3},pii_1::Array{Float64, 3},pii_2::Array{Float64, 3}, update, dt, Mass2, Lambda)
	@inn(update) = @inn(phi_1)*@inn(phi_1) + @inn(phi_2)*@inn(phi_2)
	@inn(pii_1) = @inn(pii_1) + dt * (@d2_xi(phi_1) + @d2_yi(phi_1) + @d2_zi(phi_1) - Mass2 * @inn(phi_1) - Lambda/12. * @inn(update) .* @inn(phi_1) )
	@inn(pii_2) = @inn(pii_2) + dt * (@d2_xi(phi_2) + @d2_yi(phi_2) + @d2_zi(phi_2) - Mass2 * @inn(phi_2) - Lambda/12. * @inn(update) .* @inn(phi_2) )
    return end
    #N=4
    @parallel function evolve_pii!(phi_1::Array{Float64, 3},phi_2::Array{Float64, 3},phi_3::Array{Float64, 3},phi_4::Array{Float64, 3},pii_1::Array{Float64, 3},pii_2::Array{Float64, 3},pii_3::Array{Float64, 3},pii_4::Array{Float64, 3}, update, dt, Mass2, Lambda)
	@inn(update) = @inn(phi_1).*@inn(phi_1) + @inn(phi_2).*@inn(phi_2) + @inn(phi_3).*@inn(phi_3) + @inn(phi_4).*@inn(phi_4) 
	@inn(pii_1) = @inn(pii_1) + dt * (@d2_xi(phi_1) + @d2_yi(phi_1) + @d2_zi(phi_1)- (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_1))
	@inn(pii_2) = @inn(pii_2) + dt * (@d2_xi(phi_2) + @d2_yi(phi_2) + @d2_zi(phi_2)- (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_2))
	@inn(pii_3) = @inn(pii_3) + dt * (@d2_xi(phi_3) + @d2_yi(phi_3) + @d2_zi(phi_3)- (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_3))
	@inn(pii_4) = @inn(pii_4) + dt * (@d2_xi(phi_4) + @d2_yi(phi_4) + @d2_zi(phi_4)- (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_4))
    return end
    #N=8
    @parallel function evolve_pii!(phi_1::Array{Float64, 3},phi_2::Array{Float64, 3},phi_3::Array{Float64, 3},phi_4::Array{Float64, 3},phi_5::Array{Float64, 3},phi_6::Array{Float64, 3},phi_7::Array{Float64, 3},phi_8::Array{Float64, 3},pii_1::Array{Float64, 3},pii_2::Array{Float64, 3},pii_3::Array{Float64, 3},pii_4::Array{Float64, 3},pii_5::Array{Float64, 3},pii_6::Array{Float64, 3},pii_7::Array{Float64, 3},pii_8::Array{Float64, 3}, update, dt, Mass2, Lambda)
	@inn(update) = @inn(phi_1).*@inn(phi_1) + @inn(phi_2).*@inn(phi_2) + @inn(phi_3).*@inn(phi_3) + @inn(phi_4).*@inn(phi_4) + @inn(phi_5).*@inn(phi_5) + @inn(phi_6).*@inn(phi_6) + @inn(phi_7).*@inn(phi_7) + @inn(phi_8).*@inn(phi_8) 
	@inn(pii_1) = @inn(pii_1) + dt * (@d2_xi(phi_1) + @d2_yi(phi_1) + @d2_zi(phi_1)- (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_1))
	@inn(pii_2) = @inn(pii_2) + dt * (@d2_xi(phi_2) + @d2_yi(phi_2) + @d2_zi(phi_2)- (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_2))
	@inn(pii_3) = @inn(pii_3) + dt * (@d2_xi(phi_3) + @d2_yi(phi_3) + @d2_zi(phi_3)- (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_3))
	@inn(pii_4) = @inn(pii_4) + dt * (@d2_xi(phi_4) + @d2_yi(phi_4) + @d2_zi(phi_4)- (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_4))
	@inn(pii_5) = @inn(pii_5) + dt * (@d2_xi(phi_5) + @d2_yi(phi_5) + @d2_zi(phi_5)- (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_5))
	@inn(pii_6) = @inn(pii_6) + dt * (@d2_xi(phi_6) + @d2_yi(phi_6) + @d2_zi(phi_6)- (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_6))
	@inn(pii_7) = @inn(pii_7) + dt * (@d2_xi(phi_7) + @d2_yi(phi_7) + @d2_zi(phi_7)- (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_7))
	@inn(pii_8) = @inn(pii_8) + dt * (@d2_xi(phi_8) + @d2_yi(phi_8) + @d2_zi(phi_8)- (Mass2 + Lambda/24. * @inn(update)) * @inn(phi_8))
    return end
end

#
# Time evolution functions, seriell
#
function seriell_laplace!(phi::Array{Float64, 1}) return circshift(phi, 1) .+ circshift(phi, -1)  .- 2 .* phi end
function seriell_laplace!(phi::Array{Float64, 2}) return circshift(phi, (1,0)) .+ circshift(phi, (-1,0)) .+ circshift(phi, (0,1)) .+ circshift(phi, (0,-1))  .- 4 .* phi end
function seriell_laplace!(phi::Array{Float64, 3}) return circshift(phi, (1,0,0)) .+ circshift(phi, (-1,0,0)) .+ circshift(phi, (0,1,0)) .+ circshift(phi, (0,-1,0)) .+ circshift(phi, (0,0,1)) .+ circshift(phi, (0,0,-1))  .- 6 .* phi end

function seriell_evolve_phi!(phi,pii, dt, ONgroup)
    for i in 1:ONgroup phi[i] .= phi[i] .+ dt .* pii[i] end
	return
end

function seriell_evolve_pii!(phi,pii, update, dt, Mass2, Lambda, ONgroup)
	update .= 0
    for i in 1:ONgroup update .+= phi[i].*phi[i] end
    for i in 1:ONgroup pii[i] .= pii[i] + dt * ( seriell_laplace!(phi[i]) .- Mass2 .* phi[i] .- Lambda/(6. * ONgroup) .* update .* phi[i] ) end
	return
end

#
# Time evolution functions, parallel
#
# N=1
@parallel function evolve_phi!(phi,pii,dt) @inn(phi) = @inn(phi) + dt * @inn(pii)
    return end

# for N=2
@parallel function evolve_phi!(phi_1,phi_2,pii_1,pii_2,dt)
	@inn(phi_1) = @inn(phi_1) + dt * @inn(pii_1)
	@inn(phi_2) = @inn(phi_2) + dt * @inn(pii_2)
    return end

# for N=4
@parallel function evolve_phi!(phi_1,phi_2,phi_3,phi_4,pii_1,pii_2,pii_3,pii_4,dt)
	@inn(phi_1) = @inn(phi_1) + dt * @inn(pii_1)
	@inn(phi_2) = @inn(phi_2) + dt * @inn(pii_2)
	@inn(phi_3) = @inn(phi_3) + dt * @inn(pii_3)
	@inn(phi_4) = @inn(phi_4) + dt * @inn(pii_4)
    return
end

# for N=8
@parallel function evolve_phi!(phi_1,phi_2,phi_3,phi_4,phi_5,phi_6,phi_7,phi_8,pii_1,pii_2,pii_3,pii_4,pii_5,pii_6,pii_7,pii_8,dt)
	@inn(phi_1) = @inn(phi_1) + dt * @inn(pii_1)
	@inn(phi_2) = @inn(phi_2) + dt * @inn(pii_2)
	@inn(phi_3) = @inn(phi_3) + dt * @inn(pii_3)
	@inn(phi_4) = @inn(phi_4) + dt * @inn(pii_4)
	@inn(phi_5) = @inn(phi_5) + dt * @inn(pii_5)
	@inn(phi_6) = @inn(phi_6) + dt * @inn(pii_6)
	@inn(phi_7) = @inn(phi_7) + dt * @inn(pii_7)
	@inn(phi_8) = @inn(phi_8) + dt * @inn(pii_8)
    return
end

#
# kickoff
# 
export kickoff!
function kickoff!( thesolution::QFTdynamicsSolutionCSScalar, tmpdata::Vector{CSScalarTmpDataCPU})
	println("kick off..."); flush(stdout)
    @unpack problem, simdata, measurearray, measurearrayofruns =thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem

	for i in 1:num.Runs
		myupdate_halo!(model,simdata[i].phi)
		@parallel evolve_pii!(simdata[i].phi...,simdata[i].pii..., tmpdata[1].update, 0.5*disc.dt, model.Mass2, model.Lambda)
	end
end

export kickoffsimple!
function kickoffsimple!( thesolution::QFTdynamicsSolutionCSScalar, tmpdata::Vector{CSScalarTmpDataCPU})
	println("kick off..."); flush(stdout)
    @unpack problem, simdata, measurearray, measurearrayofruns =thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem

	for i in 1:num.Runs
		myupdate_halo!(model,simdata[i].phi)
        seriell_evolve_pii!(simdata[i].phi, simdata[i].pii, tmpdata[1].update, 0.5*disc.dt, model.Mass2, model.Lambda, model.ONgroup)
	end
end

#
# evolve
# 
export evolve!
function evolve!( thesolution::QFTdynamicsSolutionCSScalar, tmpdata::Vector{CSScalarTmpDataCPU},t::Int)
    @unpack problem, simdata, measurearray, measurearrayofruns =thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem

	for i in 1:num.Runs
		@parallel evolve_phi!(simdata[i].phi..., simdata[i].pii..., disc.dt)
		myupdate_halo!(model, simdata[i].phi)
		@parallel evolve_pii!(simdata[i].phi..., simdata[i].pii..., tmpdata[1].update, disc.dt, model.Mass2, model.Lambda)
	end

    simsetup.lastEvolStep = t 
end

export evolvesimple!
function evolvesimple!( thesolution::QFTdynamicsSolutionCSScalar, tmpdata::Vector{CSScalarTmpDataCPU},t::Int)
    @unpack problem, simdata, measurearray, measurearrayofruns =thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem

	for i in 1:num.Runs
        seriell_evolve_phi!(simdata[i].phi, simdata[i].pii, disc.dt, model.ONgroup)
		myupdate_halo!(model, simdata[i].phi)
        seriell_evolve_pii!(simdata[i].phi, simdata[i].pii, tmpdata[1].update, disc.dt, model.Mass2, model.Lambda, model.ONgroup)
	end

    simsetup.lastEvolStep = t 
end