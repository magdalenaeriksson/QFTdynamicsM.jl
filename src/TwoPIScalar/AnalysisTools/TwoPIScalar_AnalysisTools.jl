using LsqFit
using Plots
using LaTeXStrings

##############################################################################################
# DampingRateFit
##############################################################################################
export DampingRateFit
Base.@kwdef mutable struct DampingRateFit
    # store the problem, @kwdef forces me to put a default thing - thats not nice but ok
    problem::QFTdynamicsProblem =  QFTdynamicsProblem( Phi4( 1, 1, 1), TwoPIScalarLambdaLO(), TwoPIScalarDiscretizationLattice(  1, 1, 1, 1, 0.02, 1, 1, 1, 1, getfftwhelper( 1, 1)), TwoPIScalarTopHatT1(), TwoPIScalarRNone(), TwoPIScalarCPUfull(1), SimSetup() )
    # r
    rMeff::Float64      = 0
    rMeff_std::Float64  = 0
    rgamma::Float64     = 0
    rgamma_std::Float64 = 0
    ronet_fit::Vector{Float64} = zeros(1)
    # F
    FMeff::Float64      = 0
    FMeff_std::Float64  = 0
    FAmpl::Float64      = 0
    FAmpl_std::Float64  = 0
    Fgamma::Float64     = 0 
    Fgamma_std::Float64 = 0
    Fonet_fit::Vector{Float64} = zeros(1)
end

export fitdampingrate
function fitdampingrate(thesolution::QFTdynamicsSolutionTwoPIScalar, fitplotpath::String)
    print("Performing fits on " * thesolution.problem.simsetup.parameterstring * "..." )
    dampingratefit = DampingRateFit()
    dampingratefit.problem = thesolution.problem 
    plots = Dict()
    linealpha = 1
    linewidth = 1.5
    #linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    # marker
    markersize = 3
    line = (:line, :solid, linealpha, linewidth); marker = (:x, markersize) 

    #
    # rho fit function
    #
    #rhofunction(x, p) = -1/p[2] * exp.(-p[1] * x) .* sin.(x*p[2])
    rhofunction(x, p) = -p[3] * exp.(-p[1] * x) .* sin.(x*p[2])
    function j_rhofunction(x,p)
        J = Matrix{Float64}(undef,length(x),length(p))
        J[:,1] = (1/p[2]) * x .* exp.(-p[1] * x) .* sin.(x*p[2])   #df/dp[1]
        J[:,2] = (1/p[2]^2) * exp.(-p[1] * x) .* (  sin.(x*p[2]) .- p[2]* x.* cos.(x*p[2]) ) #df/dp[2]
        J
    end
    #p0 = [0.02,1.2]
    #lb = [0, 0.1]
    #ub = [0.1, 5.]
    p0 = [0.02,1.2,0.8]
    lb = [0, 0.1,0.25]
    ub = [0.1, 5.,10]
    tdata =  [(j-1) * thesolution.problem.disc.dt * thesolution.problem.disc.Mass for j in 1:(thesolution.problem.disc.NstepsinMemory-1) ] 
    ydata = thesolution.measurearray[end].ronet
    #fit = curve_fit(rhofunction, j_rhofunction,tdata, ydata, p0)
    #fit = curve_fit(rhofunction, j_rhofunction,tdata, ydata, p0, lower=lb, upper=ub, maxIter=5000)
    fit = curve_fit(rhofunction, tdata, ydata, p0, lower=lb, upper=ub, maxIter=10000)
    cov = estimate_covar(fit)
    dampingratefit.rgamma       = fit.param[1]
    dampingratefit.rgamma_std   = sqrt(cov[1,1])
    dampingratefit.rMeff        = fit.param[2]
    dampingratefit.rMeff_std    = sqrt(cov[2,2])
    # evaluate fit 
    dampingratefit.ronet_fit = rhofunction(tdata ,fit.param)
    # Plot ronet and   
    #string_x = @sprintf("%5.3f",x)
    meanydata = 1/length(ydata) * sum(ydata)
    R2value = 1 - sum(fit.resid.^2) / sum( (ydata .- meanydata).^2 )
    #R2string = L"R^2 = %$( round(R2value, digits=2 ) )"
    R2string = L"1/M_{eff} = %$( round(fit.param[3], digits=4 ) )"
    Meffstring = L"M_{eff} = %$( round(dampingratefit.rMeff, digits=2 ) ) \pm %$( round(dampingratefit.rMeff_std, digits=2) )"
    gammastring = L"\gamma_{0} = %$( round(dampingratefit.rgamma, digits=4)) \pm %$( round(dampingratefit.rgamma_std, digits=4) )"
    plots[ thesolution.problem.simsetup.parameterstring * "_ronet.png"] = plot(title = Meffstring * ", " * gammastring * ", " * R2string, xlabel = L"tm", ylabel = L"\rho_0(t,t_{min})",)
    plot!(plots[ thesolution.problem.simsetup.parameterstring * "_ronet.png"], tdata, ydata, line=(line...,"red"), marker=(marker...,"red"), label="Simulation")
    plot!(plots[ thesolution.problem.simsetup.parameterstring * "_ronet.png"], tdata, dampingratefit.ronet_fit, line=(line...,"blue"), marker=(marker...,"blue"), label="fit")

    savefig( plots[ thesolution.problem.simsetup.parameterstring * "_ronet.png"] , joinpath(fitplotpath, thesolution.problem.simsetup.parameterstring * "_ronet.png") ) 
    println("Saved as " *  joinpath(fitplotpath, thesolution.problem.simsetup.parameterstring * "_ronet.png") ) ; flush(stdout)
 
    #
    # F fit function
    #
    Ffunction(x, p) = p[3] * exp.(-p[1] * x) .* cos.(x*p[2])
    function j_Ffunction(x,p)
        J = Matrix{Float64}(undef,length(x),length(p))
        J[:,1] = -x * p[3] .* exp.(-p[1] * x) .* cos.(x*p[2])  #df/dp[1]
        J[:,2] = -x * p[3] .* exp.(-p[1] * x) .* sin.(x*p[2])   #df/dp[2]
        J[:,3] = exp.(-p[1] * x) .* cos.(x*p[2])   #df/dp[2]
        J
    end
    p0 = [0.1,1,5]
    ydata = thesolution.measurearray[end].Fonet
    fit = curve_fit(Ffunction, j_Ffunction,tdata, ydata, p0)
    cov = estimate_covar(fit)
    dampingratefit.Fgamma       = fit.param[1]
    dampingratefit.Fgamma_std   = sqrt(cov[1,1])
    dampingratefit.FMeff        = fit.param[2]
    dampingratefit.FMeff_std    = sqrt(cov[2,2])
    dampingratefit.FAmpl        = fit.param[3]
    dampingratefit.FAmpl_std    = sqrt(cov[3,3])
    # evaluate fit 
    dampingratefit.Fonet_fit = Ffunction(tdata ,fit.param)
    # Plot ronet and   
    meanydata = 1/length(ydata) * sum(ydata)
    R2value = 1 - sum(fit.resid.^2) / sum( (ydata .- meanydata).^2 )
    R2string = L"R^2 = %$( round(R2value, digits=2 ) )"
    Meffstring = L"M_{eff} = %$( round(dampingratefit.FMeff, digits=2 ) ) \pm %$( round(dampingratefit.FMeff_std, digits=2) )"
    gammastring = L"\gamma_{0} = %$( round(dampingratefit.Fgamma, digits=4)) \pm %$( round(dampingratefit.Fgamma_std, digits=4) )"
    plots[ thesolution.problem.simsetup.parameterstring * "_Fonet.png"] = plot(title = Meffstring * ", " * gammastring * ", " * R2string, xlabel = L"tm", ylabel = L"F_0(t,t_{min})",)
    plot!(plots[ thesolution.problem.simsetup.parameterstring * "_Fonet.png"], tdata, ydata, line=(line...,"red"), marker=(marker...,"red"), label="Simulation")
    plot!(plots[ thesolution.problem.simsetup.parameterstring * "_Fonet.png"], tdata, dampingratefit.Fonet_fit, line=(line...,"blue"), marker=(marker...,"blue"), label="fit")

    savefig( plots[ thesolution.problem.simsetup.parameterstring * "_Fonet.png"] , joinpath(fitplotpath, thesolution.problem.simsetup.parameterstring * "_Fonet.png") ) 
    println("Saved as " *  joinpath(fitplotpath, thesolution.problem.simsetup.parameterstring * "_Fonet.png") ) ; flush(stdout)
 
    chmod(fitplotpath, 0o777, recursive=true)
    return dampingratefit
end

##############################################################################################
# DampingRateFit
##############################################################################################
export ChemPotFit
Base.@kwdef mutable struct ChemPotFit
    # store the problem, @kwdef forces me to put a default thing - thats not nice but ok
    problem::QFTdynamicsProblem =  QFTdynamicsProblem( Phi4( 1, 1, 1), TwoPIScalarLambdaLO(), TwoPIScalarDiscretizationLattice(  1, 1, 1, 1, 0.02, 1, 1, 1, 1, getfftwhelper( 1, 1)), TwoPIScalarTopHatT1(), TwoPIScalarRNone(), TwoPIScalarCPUfull(1), SimSetup() )
    # r
    Teff::Float64      = 0
    Teff_std::Float64  = 0
    mueff::Float64     = 0
    mueff_std::Float64 = 0
    BE_fit::Vector{Float64} = zeros(1)
end

export fitchempot
function fitchempot(thesolution::QFTdynamicsSolutionTwoPIScalar, fitplotpath::String)
    print("Performing fits on " * thesolution.problem.simsetup.parameterstring * "..." )
    chempotfit = ChemPotFit()
    chempotfit.problem = thesolution.problem 
    plots = Dict()
    linealpha = 1
    linewidth = 1.5
    #linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    # marker
    markersize = 3
    line = (:line, :solid, linealpha, linewidth); marker = (:x, markersize) 

    #
    # fit function
    #
    #rhofunction(x, p) = -1/p[2] * exp.(-p[1] * x) .* sin.(x*p[2])
    BEfunction(x, p) = p[1] .* (x .- p[2])
    p0 = [1.,1.]
    ##
    tmparray = (1 ./ thesolution.measurearray[end].n) .+ 1
    idxlog = tmparray .> 0.001 # for 1/n + 1 values for which the log function works
    tdata =  thesolution.measurearray[end].omega[idxlog]
    ydata = log.(1 ./ thesolution.measurearray[end].n[idxlog] .+ 1 )
    fit = curve_fit(BEfunction, tdata, ydata, p0, maxIter=500)
    cov = estimate_covar(fit)

    chempotfit.Teff       = 1/fit.param[1]
    chempotfit.Teff_std   = sqrt(cov[1,1])
    chempotfit.mueff      = fit.param[2]
    chempotfit.mueff_std  = sqrt(cov[2,2])
    # evaluate fit 
    chempotfit.BE_fit = BEfunction(tdata ,fit.param)
    # Plot ronet and   
    #string_x = @sprintf("%5.3f",x)
    meanydata = 1/length(ydata) * sum(ydata)
    R2value = 1 - sum(fit.resid.^2) / sum( (ydata .- meanydata).^2 )
    R2string = L"R^2 = %$( round(R2value, digits=2 ) )"
    Teffstring = L"T_{eff} = %$( round(chempotfit.Teff, digits=2 ) ) "
    mueffstring = L"\mu_{eff} = %$( round(chempotfit.mueff, digits=4)) \pm %$( round(chempotfit.mueff_std, digits=4) )"
    plots[ thesolution.problem.simsetup.parameterstring * "_BE.png"] = plot(title = Teffstring * ", " * mueffstring * ", " * R2string, xlabel = L"\omega_k", ylabel = L"Log[1+1/n_k(t)]",)
    plot!(plots[ thesolution.problem.simsetup.parameterstring * "_BE.png"], tdata, ydata, line=(line...,"red"), marker=(marker...,"red"), label="Simulation")
    plot!(plots[ thesolution.problem.simsetup.parameterstring * "_BE.png"], tdata, chempotfit.BE_fit, line=(line...,"blue"), marker=(marker...,"blue"), label="fit")

    savefig( plots[ thesolution.problem.simsetup.parameterstring * "_BE.png"] , joinpath(fitplotpath, thesolution.problem.simsetup.parameterstring * "_BE.png") ) 
    println("Saved as " *  joinpath(fitplotpath, thesolution.problem.simsetup.parameterstring * "_BE.png") ) ; flush(stdout)
 
    return chempotfit
end
