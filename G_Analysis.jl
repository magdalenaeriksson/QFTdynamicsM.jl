###################################################################################################################################
using QFTdynamics
using Serialization
using Plots
using LaTeXStrings
using PolyLog
using QuadGK
using Roots

###################################################################################################################################
# Set up
###################################################################################################################################
# dir name for fit in plotpath
fitname = "3p1EquilibrationTincor"
# simpaths to be analysed
simpaths = ["phi4LambdaNLO_L600_ON1_RMass_T27_Nx16_sdim3_Nsteps1400_dt100_M70_NstepsinMemory400_EquilibrationTincor",
            "phi4LambdaNLO_L600_ON1_RMass_T22_Nx16_sdim3_Nsteps1400_dt100_M70_NstepsinMemory400_EquilibrationTincor",
            "phi4LambdaNLO_L600_ON1_RMass_T18_Nx16_sdim3_Nsteps1400_dt100_M70_NstepsinMemory400_EquilibrationTincor",
            "phi4LambdaNLO_L600_ON1_RMass_T14_Nx16_sdim3_Nsteps1400_dt100_M70_NstepsinMemory800_EquilibrationTincor",
            "phi4LambdaNLO_L600_ON1_RMass_T11_Nx16_sdim3_Nsteps1400_dt100_M70_NstepsinMemory800_EquilibrationTincor",
            "phi4LambdaNLO_L600_ON1_RMass_T7_Nx16_sdim3_Nsteps1400_dt100_M70_NstepsinMemory800_EquilibrationTincor",
            "phi4LambdaNLO_L600_ON1_RMass_T3_Nx16_sdim3_Nsteps1400_dt100_M70_NstepsinMemory800_EquilibrationTincor",]
###################################################################################################################################

# read paths
lines = readlines("etc/paths.txt")
datapath = chop(lines[2],head=5,tail=0)
postpath = chop(lines[3],head=5,tail=0)
plotpath = chop(lines[4],head=5,tail=0)

# create pathtostore
fitpostpath = postpath * "/" * fitname
fitplotpath = plotpath * "/" * fitname

rm(fitpostpath, force=true, recursive=true)
rm(fitplotpath, force=true, recursive=true)

mkpath(fitpostpath)
mkpath(fitplotpath)
mkpath(fitplotpath * "/pdf")
   
cp("etc/index.php", joinpath( join(split(fitplotpath,"/")[1:end-1],"/"),"index.php"), force=true )
cp("etc/index.php", joinpath(fitplotpath,"index.php"), force=true )
cp("etc/index.php", joinpath(fitplotpath * "/pdf" ,"index.php"), force=true )
chmod( join(split(fitplotpath,"/")[1:end-1],"/") , 0o777, recursive=true)

###################################################################################################################################
# perform fits
###################################################################################################################################
fitresults = Vector{DampingRateFit}(undef, length(simpaths))
for (i,simpath) in enumerate(simpaths)
    #load solution
    thesolution = deserialize( datapath * "/" * simpath * "/SimSolution.jld"  )
    # make dampingrate fit
    fitresults[i] = fitdampingrate(thesolution, fitplotpath)
    # make chempot fit
    fitchempot(thesolution, fitplotpath)
end
###################################################################################################################################

###################################################################################################################################
# make the plot
###################################################################################################################################
plots = Dict()
plots["Meff.png"] = plot(xlabel = L"T/m", ylabel = L"M_{eff}/m")
plots["gamma.png"] = plot(xlabel = L"T/m", ylabel = L"\gamma/m")

linealpha = 1
linewidth = 1.5
#linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
# marker
markersize = 3
line = (:line, :solid, linealpha, linewidth); marker = (:x, markersize) 

#####################################################################################################
# Perturbative expressions from Continuum
#####################################################################################################
### dampingrate
function dampingrate(MH, lambda,T)
    println("calc dampingrate cont for T=", T)
    res =  ((lambda^2*T^2)/(128*pi^3*MH) * reli2(exp(-MH/T)))
    @show res
    return res
end

### Thermal Hartreemass
function ThermalHartreemass2_1loopcontinuum(MH, lambda,T,mu) 
    # lambda -> dimless
    # T -> 1/m
    # mu -> 1/m
    # MH2 -> 1/m^2
    # returns MH2/m^2 (not self consistent of course)
    integral, err = quadgk(x -> sqrt(x^2-(MH/T)^2)/(exp(x-mu)-1), MH/T, Inf, rtol=1e-8)
    return 1 + (lambda*T^2)/(4*pi^2) * integral
end

function getThermalHartreemass2_1loopcontinuum_selfconsistent(lambda,T,mu) 
    f(x) = x^2 - ThermalHartreemass2_1loopcontinuum(x, lambda,T,mu) 
    return find_zero(f, (0.1,5), Bisection())
end

###################################################################################################################################
# make and save plots
###################################################################################################################################
plot!(plots["gamma.png"], [ res.problem.init.T for res in fitresults], [res.rgamma for res in fitresults], yerr= [res.rgamma_std for res in fitresults], line=(line...,"red"), marker=(marker...,"red"), label=L"\rho_0(t,t^{\prime})")
plot!(plots["gamma.png"], [ res.problem.init.T for res in fitresults], [res.Fgamma for res in fitresults], yerr= [res.Fgamma_std for res in fitresults], line=(line...,"blue"), marker=(marker...,"blue"), label=L"F_0(t,t^{\prime})")
plot!(plots["gamma.png"], [ res.problem.init.T for res in fitresults], [ dampingrate(getThermalHartreemass2_1loopcontinuum_selfconsistent(res.problem.model.Lambda,res.problem.init.T, 0.),res.problem.model.Lambda, res.problem.init.T) for res in fitresults], line=(line...,"black"), marker=(marker...,"black"), label=L"Perturbative (Continuum)")

plot!(plots["Meff.png"],  [ res.problem.init.T for res in fitresults], [res.rMeff for res in fitresults], yerr= [res.rMeff_std for res in fitresults],line=(line...,"red"), marker=(marker...,"red"), label=L"\rho_0(t,t^{\prime})")
plot!(plots["Meff.png"],  [ res.problem.init.T for res in fitresults], [res.FMeff for res in fitresults], yerr= [res.FMeff_std for res in fitresults],line=(line...,"blue"), marker=(marker...,"blue"), label=L"F_0(t,t^{\prime})")
plot!(plots["Meff.png"],  [ res.problem.init.T for res in fitresults], [ getThermalGapMass(res.problem.model, res.problem.pexp, res.problem.disc, res.problem.init, 0.)  for res in fitresults], yerr= [res.FMeff_std for res in fitresults],line=(line...,"black"), marker=(marker...,"black"), label=L"Perturbative(Lattice, mu=0)")
plot!(plots["Meff.png"],  [ res.problem.init.T for res in fitresults], [ getThermalHartreemass2_1loopcontinuum_selfconsistent(res.problem.model.Lambda,res.problem.init.T, 0.)  for res in fitresults],line=(line...,"green"), marker=(marker...,"green"), label=L"Perturbative(Continuum, mu=0)")

savefig( plots[ "gamma.png"] , joinpath(fitplotpath, "gamma.png") ) 
savefig( plots[ "gamma.png"] , joinpath(fitplotpath, "pdf/gamma.pdf") ) 
println("Saved as " *  joinpath(fitplotpath, "gamma.png") ) ; flush(stdout)
savefig( plots[ "Meff.png"] , joinpath(fitplotpath, "Meff.png") ) 
savefig( plots[ "Meff.png"] , joinpath(fitplotpath, "pdf/Meff.pdf") ) 
println("Saved as " *  joinpath(fitplotpath, "Meff.png") ) ; flush(stdout)
chmod(fitplotpath, 0o777, recursive=true)
