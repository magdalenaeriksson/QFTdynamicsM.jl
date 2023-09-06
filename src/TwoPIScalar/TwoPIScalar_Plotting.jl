using Plots
using LaTeXStrings
using ColorSchemes

export plotdata
function plotdata(plots::Dict, sol::QFTdynamicsSolutionTwoPIScalar, mode::String, style::Int64=1, label::String="test")
    # linestyle
    # line
    linealpha = 1
    linewidth = 1.5
    #linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    # marker
    markersize = 3
    if style == 1 line = (:line, :solid, linealpha, linewidth); marker = (:x, markersize) end
    if style == 2 line = (:line, :dash, linealpha, linewidth); marker = (:o, markersize)  end
    if style == 3 line = (:line, :dot, linealpha, linewidth); marker = (:star5, markersize) end
    if style == 4 line = (:line, :dashdot, linealpha, linewidth); marker = (:circle, markersize) end
    # define colors
    color = ColorSchemes.Austria

    # unpack
    @unpack problem, simdata, measurearray = sol
    @unpack model, pexp, disc, init, reno, simsetup = problem

    # time snapshots to plot
    if length(measurearray) > 2 timeidx = [1, Int64(floor( length(measurearray)/2) ), length(measurearray)] end
    if length(measurearray) > 4 timeidx = [1, Int64(floor( length(measurearray)/4) ), Int64(floor( length(measurearray)/2) ), Int64(floor( (3*length(measurearray))/4) ), length(measurearray)] end

    # momenta to plot
    nkLmomenta = 5
    kLmax = disc.fftwhelper[end].lev
    kLidx = zeros(Int64, nkLmomenta)
    for i in 0:(length(kLidx)-1)
        kLidx[i+1] = Int64(findmin(abs.([element.lev for element in disc.fftwhelper] .- ( (i/(nkLmomenta-1))*kLmax ) ))[2])
    end

### TIME EVOLUTION PLOTS
    # plot EM(t)
    if mode == "c" plots["Energy.png"] = plot(xlabel = L"tm",ylabel = L"E")
    else
        plot!(plots["Energy.png"],[x.time for x in measurearray], [x.E for x in measurearray], line=(line...,color[1]), marker=(marker...,color[1]), label=label)
    end

    if mode == "c" plots["Energydetail.png"] = plot(xlabel = L"tm",ylabel = L"E")
    else
        plot!(plots["Energydetail.png"],[x.time for x in measurearray], [x.E_kin for x in measurearray], line=(line...,color[1]), marker=(marker...,"red"), label=label * L"E_{kin}")
        plot!(plots["Energydetail.png"],[x.time for x in measurearray], [x.E_pot for x in measurearray], line=(line...,color[1]), marker=(marker...,"blue"), label=label * L"E_{pot}")
        plot!(plots["Energydetail.png"],[x.time for x in measurearray], [x.E_mem for x in measurearray], line=(line...,color[1]), marker=(marker...,"black"), label=label * L"E_{mem}")
    end
    # plot M(t)
    if mode == "c" plots["Hartreemass.png"] = plot(xlabel = L"tm",ylabel = L"M^2(t)",)
    else
        plot!(plots["Hartreemass.png"],[element.time for element in measurearray], [element.M2 for element in measurearray], line=(line...,color[1]), marker=(marker...,color[1]), label=label)
    end
    # plot F_0(t,t)
    if mode == "c" plots["stat0prop.png"] = plot(xlabel = L"tm",ylabel = L"F_0(t,t)",)
    else
        plot!(plots["stat0prop.png"],[element.time for element in measurearray], [element.F[1] for element in measurearray], line=(line...,color[1]), marker=(marker...,color[1]), label=label)
    end
    # plot n_k(t)
    if mode == "c" plots["particlenumber_evol.png"] = plot(xlabel = L"tm",ylabel = L"n_k",)
    else
        for (i, kidx) in enumerate(kLidx)
        plot!(plots["particlenumber_evol.png"],[element.time for element in measurearray], [element.n[kidx] for element in measurearray], line=(line...,color[i]), marker=(marker...,color[i]), label=label * ": k=" * string(round(disc.fftwhelper[kidx].lev, digits=1)) )
        end
    end
    # plot F_x=0(t,t)
    if mode == "c" plots["statpropat0.png"] = plot(xlabel = L"tm",ylabel = L"F(t,t,x=0)/m^2",)
    else
        if typeof(problem.disc) == QFTdynamics.TwoPIScalarDiscretizationLattice
            tmparray = zeros(length(measurearray))
            for i in 1:length(measurearray)
                for j in 1:length(disc.fftwhelper)
                    tmparray[i] += measurearray[i].F[j] * disc.fftwhelper[j].deg
                end
            end
            tmparray .-= tmparray[1]
            tmparray .*= disc.ivol/disc.Mass^2
            plot!(plots["statpropat0.png"],[element.time for element in measurearray], tmparray, line=(line...,color[1]), marker=(marker...,color[1]), label=label)
        end
        if typeof(problem.disc) == QFTdynamics.TwoPIScalarDiscretizationCont
            # TODO
        end
    end
### SNAPSHOT PLOTS
    # plot c(k)
    if mode == "c" plots["cnumber.png"] = plot(xlabel = L"k_L",ylabel = L"c_k",)
    else
        for (i, tidx) in enumerate(timeidx)
        plot!(plots["cnumber.png"],[x.lev for x in disc.fftwhelper], measurearray[tidx].c, line=(line...,color[i]), marker=(marker...,color[i]), label=label * ": tm=" * string(round(measurearray[tidx].time,digits=1)))
        end
    end
    # plot n(k)
    if mode == "c" plots["particlenumber.png"] = plot(xlabel = L"k_L",ylabel = L"n_k",)
    else
        for (i, tidx) in enumerate(timeidx)
        plot!(plots["particlenumber.png"],[x.lev for x in disc.fftwhelper], measurearray[tidx].n, line=(line...,color[i]), marker=(marker...,color[i]), label=label * ": tm=" * string(round(measurearray[tidx].time,digits=1)))
        end
    end
    # plot r(k)
    #if mode == "c" plots["specprop.png"] = plot(xlabel = L"k_L",ylabel = L"r_k(t,t)",)
    #else
    #    for (i, tidx) in enumerate(timeidx)
    #    plot!(plots["specprop.png"],[x.lev for x in disc.fftwhelper], measurearray[tidx].r, line=(line...,color[i]), marker=(marker...,color[i]), label=label * ": tm=" * string(round(measurearray[tidx].time,digits=1)))
    #    end
    #end
    # plot F(k)
    if mode == "c" plots["statprop.png"] = plot(xlabel = L"k_L",ylabel = L"F_k(t,t)",)
    else
        for (i, tidx) in enumerate(timeidx)
        plot!(plots["statprop.png"],[x.lev for x in disc.fftwhelper], measurearray[tidx].F, line=(line...,color[i]), marker=(marker...,color[i]), label=label * ": tm=" * string(round(measurearray[tidx].time,digits=1)))
        end
    end
    # dispersion relation
    if mode == "c" plots["dispersionrelation.png"] = plot( xlabel = L"k_L", ylabel = L"\epsilon_k")
    else
        for (i, tidx) in enumerate(timeidx)
        plot!(plots["dispersionrelation.png"],[x.lev for x in disc.fftwhelper], measurearray[tidx].omega, line=(line...,color[i]), marker=(marker...,color[i]), label=label * ": tm=" * string(round(measurearray[tidx].time,digits=1)))
        end
    end
    # plot Log (1 + 1/np)
    if mode == "c" plots["BEdistribution.png"] = plot(xlabel = L"\epsilon_k",ylabel = L"Log[1+1/n_k]",)
    else
        for (i, tidx) in enumerate(timeidx)
        tmparray = (1 ./ measurearray[tidx].n) .+ 1
        idxlog = tmparray .> 0.001 # for 1/n + 1 values for which the log function works
        plot!(plots["BEdistribution.png"],measurearray[tidx].omega[idxlog], log.(tmparray[idxlog]), line=(line...,color[i]), marker=(marker...,color[i]), label=label * ": tm=" * string(round(measurearray[tidx].time,digits=1)))
        end
    end
    # plot np(omega)
    if mode == "c" plots["particledistribution.png"] = plot(xlabel = L"\epsilon_k",ylabel = L"n_k",)
    else
        for (i, tidx) in enumerate(timeidx)
        plot!(plots["particledistribution.png"],measurearray[tidx].omega, measurearray[tidx].n, line=(line...,color[i]), marker=(marker...,color[i]), label=label * ": tm=" * string(round(measurearray[tidx].time,digits=1)))
        end
    end
### MEMORY PLOT
    # Self Energy
    if mode == "c" plots["SigmaF.png"] = plot(xlabel = L"(t-t^\prime)m" ,ylabel = L"\Sigma^F_0",)
    else
    plot!(plots["SigmaF.png"], [i * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].SigF, line=(line...,color[1]), marker=(marker...,color[1]), label= label)
    end

    if mode == "c" plots["Sigmar.png"] = plot(xlabel = L"(t-t^\prime)m" ,ylabel = L"\Sigma^{\rho}_0",)
    else
    plot!(plots["Sigmar.png"], [i * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].Sigr, line=(line...,color[1]), marker=(marker...,color[1]), label= label)
    end
    # plot F_0(t,1)
    if mode == "c" plots["offdiagF.png"] = plot(xlabel = L"tm", ylabel = L"F_0(t,t_{min})",)
    else
    #plot!(plots["offdiagF.png"], [(i-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].Fonet, line=(line...,color[1]), marker=(marker...,color[1]), label=label)
    plot!(plots["offdiagF.png"], [(j-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].Fonet, line=(line...,color[1]), marker=(marker...,color[1]), label=label)
    end

    # plot r_0(t,1)
    if mode == "c" plots["offdiagr.png"] = plot(xlabel = L"tm", ylabel = L"\rho_0(t,t_{min})",)
    else
    #plot!(plots["offdiagr.png"], [(i-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].ronet, line=(line...,color[1]), marker=(marker...,color[1]), label=label)
    plot!(plots["offdiagr.png"], [(j-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].ronet, line=(line...,color[1]), marker=(marker...,color[1]), label=label)
    end
end

export plotsimdataF
function plotsimdataF(plots::Dict, sol::QFTdynamicsSolutionTwoPIScalar, mode::String, style::Int64=1, label::String="test")
    println("plotsimdataF: Works only for Memt=0 (=NstepsinMemory+2)!")
    # linestyle
    # line
    linealpha = 1
    linewidth = 1.5
    #linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    # marker
    markersize = 3
    if style == 1 line = (:line, :solid, linealpha, linewidth); marker = (:x, markersize) end
    if style == 2 line = (:line, :dash, linealpha, linewidth); marker = (:o, markersize)  end
    if style == 3 line = (:line, :dot, linealpha, linewidth); marker = (:star5, markersize) end
    if style == 4 line = (:line, :dashdot, linealpha, linewidth); marker = (:circle, markersize) end
    # define colors
    color = ColorSchemes.Austria

    # unpack
    @unpack problem, simdata, measurearray = sol
    @unpack model, pexp, disc, init, reno, simsetup = problem

    # time snapshots to plot
    if length(measurearray) > 2 timeidx = [1, Int64(floor( length(measurearray)/2) ), length(measurearray)] end
    if length(measurearray) > 4 timeidx = [1, Int64(floor( length(measurearray)/4) ), Int64(floor( length(measurearray)/2) ), Int64(floor( (3*length(measurearray))/4) ), length(measurearray)] end

    for time in timeidx
        if mode == "c" plots["simdataF_timeidx" * string(time) * ".png"] = plot(xlabel = L"k_L",ylabel = L"F_k(t,t)",)
        else
            # collect things to plot
            plot!(plots["simdataF_timeidx" * string(time) * ".png"],[x.lev for x in disc.fftwhelper], measurearray[time].F, line=(line...,color[1]), marker=(marker...,color[1]), label="Measurarray: tm=" * string(round(measurearray[time].time,digits=1)))
            plot!(plots["simdataF_timeidx" * string(time) * ".png"],[x.lev for x in disc.fftwhelper], [simdata.F[time,time,element.ind[end]] for element in disc.fftwhelper], line=(line...,color[1]), marker=(marker...,color[1]),                         label= "timeidx" * string(time) * ", last deg")
            plot!(plots["simdataF_timeidx" * string(time) * ".png"],[x.lev for x in disc.fftwhelper], [simdata.F[time,time,element.ind[round(Int, end/2, RoundUp)]] for element in disc.fftwhelper], line=(line...,color[3]), marker=(marker...,color[3]),  label= "timeidx" * string(time) * ", middle deg")
            plot!(plots["simdataF_timeidx" * string(time) * ".png"],[x.lev for x in disc.fftwhelper], [simdata.F[time,time,element.ind[1]] for element in disc.fftwhelper], line=(line...,color[2]), marker=(marker...,color[2]),                           label= "timeidx" * string(time) * ", first deg")
        end
    end
end

export addLOplot
function addLOplot(plots::Dict, sol::QFTdynamicsSolutionTwoPIScalar)
    label = "LO M(t)->M"
    # linestyle
    # line
    linealpha = 1
    linewidth = 1.5
    #linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    style = 1
    if style == 1 line = (:line, :solid, linealpha, linewidth) end
    if style == 2 line = (:line, :dash, linealpha, linewidth) end
    if style == 3 line = (:line, :dot, linealpha, linewidth) end
    if style == 4 line = (:line, :dashdot, linealpha, linewidth) end
    # marker
    markersize = 3
    marker = (:x, markersize)
    # define colors
    color = ColorSchemes.Austria

    # unpack
    @unpack problem, simdata, measurearray = sol
    @unpack model, pexp, disc, init, reno, simsetup = problem

    # particle number
    #function F(idx, t,tp, )
    #    n = getparticlenr(init, disc)
    #    return ((n[idx] + 0.5)/sqrt(simdata.kL2values[idx]+disc.Mass^2)) * cos( sqrt(simdata.kL2values[idx]+disc.Mass^2) * (t-tp) )
    #end
    #function r(idx, t,tp)
    #    return (1/sqrt(simdata.kL2values[idx]+disc.Mass^2)) * sin( sqrt(simdata.kL2values[idx]+disc.Mass^2) * (t-tp) )
    #end
    # gauss number
    function F(idx, t,tp)
        #return init.xi^2 * cos( sqrt(simdata.kL2values[idx]+disc.Mass^2) * t ) * cos( sqrt(simdata.kL2values[idx]+disc.Mass^2) * tp ) + 1/disc.Mass^2 * (init.eta^2 + init.sig^2/(4*init.xi^2) ) * sin( sqrt(simdata.kL2values[idx]+disc.Mass^2) * t ) * sin( sqrt(simdata.kL2values[idx]+disc.Mass^2) * tp ) + 1/disc.Mass * (init.eta *init.xi) * sin( sqrt(simdata.kL2values[idx]+disc.Mass^2) * (t+tp) )
        return init.xi^2*cos(disc.Mass*(t-disc.dt))*cos(disc.Mass*(tp-disc.dt))+1/disc.Mass^2*(init.eta^2 + init.sig^2/(4*init.xi^2))*sin(disc.Mass*(t-disc.dt))*sin(disc.Mass*(tp-disc.dt)) + 1/disc.Mass*(init.eta*init.xi)*sin(disc.Mass*(t+tp-disc.dt*2))
    end

    function r(idx, t,tp)
        #return (1/sqrt(simdata.kL2values[idx]+disc.Mass^2)) * sin( sqrt(simdata.kL2values[idx]+disc.Mass^2) * (t-tp) )
        return (1/disc.Mass) * sin(disc.Mass*(t-tp))
    end

    # plot F_0(t,tp)
    plot!(plots["stat0prop.png"],[element.time for element in measurearray], [F(1,element.time,element.time) for element in measurearray], line=(line...,"black"), marker=(marker...,"black"), label=label)
    # plot F_0(t,tp)
    # x axis: [0 ,dt , 2*dt,,,, NstepsinMemory], y axis: [F(last, last) , F(last-1, last) ]
    plot!(plots["offdiagF.png"], [(i-1)*disc.dt*disc.Mass for i in 1:(simdata.NstepsinMemory-1)], [F(1,(simdata.indices[2]-1)*disc.dt*disc.Mass, (j)*disc.dt*disc.Mass) for j in (simdata.indices[2]-1):-1:simdata.indices[1] ], line=(line...,"black"), marker=(marker...,"black"), label=label)
    # plot r_0(t,tp)
    # x axis: [0 ,dt , 2*dt,,,, NstepsinMemory], y axis: [F(last, last) , F(last-1, last) ]
    plot!(plots["offdiagr.png"], [(i-1)*disc.dt*disc.Mass for i in 1:(simdata.NstepsinMemory-1)], [r(1,(simdata.indices[2]-1)*disc.dt*disc.Mass, (j)*disc.dt*disc.Mass) for j in (simdata.indices[2]-1):-1:simdata.indices[1] ], line=(line...,"black"), marker=(marker...,"black"), label=label)
end

export plotdata_anders
function plotdata_anders(plots::Dict, sol::QFTdynamicsSolutionTwoPIScalar, mode::String, style::Int64=1, label::String="test")

    # linestyle
    # line
    linealpha = 1
    linewidth = 1.5
    #linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    if style == 1 line = (:line, :solid, linealpha, linewidth) end
    if style == 2 line = (:line, :dash, linealpha, linewidth) end
    if style == 3 line = (:line, :dot, linealpha, linewidth) end
    if style == 4 line = (:line, :dashdot, linealpha, linewidth) end
    # marker
    markersize = 3
    marker = (:x, markersize)
    # define colors
    color = ColorSchemes.Austria

    # unpack
    @unpack problem, simdata, measurearray = sol
    @unpack model, pexp, disc, init, reno, simsetup = problem

    stepstoplot = [1, Int(floor(length(measurearray)/2)), length(measurearray)] 
    mtvaluestoplot = [7,14,21,28,42,56,84,280,476]
    stepstoplot = zero(mtvaluestoplot)
    for i in 1:length(stepstoplot)
        stepstoplot[i] = Int(findmin(abs.([element.time for element in measurearray].-mtvaluestoplot[i]))[2])
    end
    # 1. nk vs omegak (NOTE: not sure which omegak to use!)
    if mode == "c" plots["Anders_particlenumber.png"] = scatter(layout=length(stepstoplot), xlim=[0,4], ylim=[0,5])
    else
        for (i, step) in enumerate(stepstoplot)
            #plot!(plots["Anders_particlenumber.png"][i], measurearray[step].omega, measurearray[step].n, label= "tm=" * string( round(measurearray[step].time,digits=1)), marker=:x)
            scatter!(plots["Anders_particlenumber.png"][i], [element.lev for element in disc.fftwhelper], measurearray[step].n, label= "tm=" * string( round(measurearray[step].time,digits=1)))
        end
    end
    # 2. omegak vs k2 (NOTE: not sure which omegak to use!)
    if mode == "c" plots["Anders_dispersionrelation.png"] = plot(layout=length(stepstoplot), xlim=[0,30], ylim=[0,30])
    else
        for (i, step) in enumerate(stepstoplot)
            scatter!(plots["Anders_dispersionrelation.png"][i], [element.lev2 for element in disc.fftwhelper], measurearray[step].omega.^2, label= "tm=" * string(round(measurearray[step].time,digits=1)))
        end
    end
    # 3. nk vs mt
    if mode == "c" plots["Anders_particlenumberevolution.png"] = plot(xlabel = L"mt",ylabel = L"n_k", legend=false)
    else
        #modestoplot = [1,Int(floor(length(disc.fftwhelper)/2)), length(disc.fftwhelper)]
        # 5 modes between kl = 2.04 and kl2 = 6.12
        k2min = 0.68 * disc.sdim
        idxmin = Int(findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2min))[2])
        k2max = 2.04 * disc.sdim
        idxmax = Int(findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2max))[2])
        modestoplot = [idxmin, Int64(floor(idxmin+(idxmax-idxmin)/2)) , idxmax]
        for mode in modestoplot
            plot!(plots["Anders_particlenumberevolution.png"],[element.time for element in measurearray],[element.n[mode] for element in measurearray],)
        end
    end
    # 4. ln(1+1/nk) at mt =1000 
    if mode == "c" plots["Anders_BEdistribution.png"] = scatter(xlabel = L"\omega_k",ylabel = L"Log[1+1/n_k(t)]")
    else
        tmparray = (1 ./ measurearray[end].n) .+ 1
        idxlog = tmparray .> 0.001 # for 1/n + 1 values for which the log function works
        scatter!(plots["Anders_BEdistribution.png"], measurearray[end].omega[idxlog], log.(tmparray[idxlog]), label= "tm=" * string(round(measurearray[end].time,digits=1)))
    end

end
