using Plots
using LaTeXStrings
using ColorSchemes

export plotdata
function plotdata(plots::Dict, thesolution::QFTdynamicsSolutionCSScalar, mode::String, style::Int64=1, label::String="test")
    # linestyle
    # line
    linealpha = 1
    linewidth = 1.5
    #linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    # marker
    markersize = 3
    #if style == 1 line = (:line, :solid, linealpha, linewidth); marker = (:x, markersize)           end
    #if style == 2 line = (:line, :dash, linealpha, linewidth); marker = (:o, markersize)            end
    #if style == 3 line = (:line, :dot, linealpha, linewidth); marker = (:star5, markersize)         end
    #if style == 4 line = (:line, :dashdot, linealpha, linewidth); marker = (:circle, markersize)    end
    if style == 1 line = :solid   ; marker = :x        end
    if style == 2 line = :dash    ; marker = :o        end
    if style == 3 line = :dot     ; marker = :star5    end
    if style == 4 line = :dashdot ; marker = :circle   end
    # define colors
    color = ColorSchemes.Austria

    # unpack
    @unpack problem, simdata, measurearray, measurearrayofruns = thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem

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
## degeneracy plot
    if mode == "c" plots["DegeneracySingleRun.png"] = plot(xlabel = L"k_L",ylabel = L"deg") end
    if mode == "h"
        #if hasproperty(disc,:deg)
        if haskey(plots, "DegeneracySingleRun.png")
            plot!(plots["DegeneracySingleRun.png"],[el.lev for el in disc.fftwhelper], disc.deg ,ls=line, lc=color[1], markerstrokecolor=:auto, label=label)
            annotate!(plots["DegeneracySingleRun.png"], 0, maximum(disc.deg), text( L"N_{Runs} = " * string(num.Runs), :black, :left, 10))
        end
    end

### TIME EVOLUTION PLOTS
    # plot EM(t)
    if mode == "c" plots["Energy.png"] = plot(xlabel = L"tm",ylabel = L"E")
    else
        plot!(plots["Energy.png"],[x.time for x in measurearray], [x.E for x in measurearray],ls=line, lc=color[1], markerstrokecolor=:auto, label=label)
    end
    if mode == "c" plots["deltaEnergy.png"] = plot(xlabel = L"tm",ylabel = L"\delta E/E (\%)")
    else
        meanEnergy = sum([x.E for x in measurearray])/ length(measurearray)
        plot!(plots["deltaEnergy.png"],[x.time for x in measurearray], [100*((x.E-meanEnergy)/meanEnergy) for x in measurearray],ls=line, lc=color[1], markerstrokecolor=:auto, label=label)
    end
    if mode == "c" plots["Energydetail.png"] = plot(xlabel = L"tm",ylabel = L"E")
    else
        plot!(plots["Energydetail.png"],[x.time for x in measurearray], [x.E_phi for x in measurearray], line=(line,"red"), marker=(marker,"red"), label=label * L"E_{\pi}")
        plot!(plots["Energydetail.png"],[x.time for x in measurearray], [x.E_pii for x in measurearray], line=(line,"blue"), marker=(marker,"blue"), label=label * L"E_{\phi}")
    end
    # plot F_0(t,t)
    if mode == "c" plots["stat0prop.png"] = plot(xlabel = L"tm",ylabel = L"F_0(t,t)",)
    else
        plot!(plots["stat0prop.png"],[element.time for element in measurearray], [element.phi2[1] for element in measurearray], yerr=[element.phi2_err[1] for element in measurearray], ls=line, lc=color[1], markerstrokecolor=:auto, label=label)
    end
    # plot n_k(t)
    if mode == "c" plots["particlenumber_evol.png"] = plot(xlabel = L"tm",ylabel = L"n_k",)
    else
        for (i, kidx) in enumerate(kLidx)
        plot!(plots["particlenumber_evol.png"],[element.time for element in measurearray], [element.n[kidx] for element in measurearray], yerr=[element.n_err[kidx] for element in measurearray], ls=line, lc=color[i], markerstrokecolor=:auto, label=label * ": k=" * string(round(disc.fftwhelper[kidx].lev, digits=1)) )
        end
    end
    # plot F_x=0(t,t)
    if mode == "c" plots["statpropat0.png"] = plot(xlabel = L"tm",ylabel = L"F(t,t,x=0)/m^2",)
    else
        tmparray = zeros(length(measurearray))
        tmparray_err = zeros(length(measurearray))
        for i in 1:length(measurearray)
            for j in 1:length(disc.fftwhelper)
                tmparray[i] += measurearray[i].phi2[j] * disc.fftwhelper[j].deg
                tmparray_err[i] += measurearray[i].phi2_err[j]^2 * disc.fftwhelper[j].deg
            end
        end
        tmparray .-= tmparray[1]
        tmparray .*= disc.ivol/disc.Mass^2
        tmparray_err .= sqrt.(tmparray_err)
        tmparray_err .*= disc.ivol/disc.Mass^2 
        plot!(plots["statpropat0.png"],[element.time for element in measurearray], tmparray, yerr=tmparray_err, line=(line,color[1]), marker=(marker,color[1]), label=label)
    end
## SNAPSHOT PLOTS
    # plot n(k)
    if mode == "c" plots["particlenumber.png"] = plot(xlabel = L"k_L",ylabel = L"n_k",)
    else
        for (i, tidx) in enumerate(timeidx)
        plot!(   plots["particlenumber.png"],[x.lev for x in disc.fftwhelper], measurearray[tidx].n, yerr=measurearray[tidx].n_err, ls=line, lc=color[i], markerstrokecolor=:auto, label=label * ": tm=" * string(round(measurearray[tidx].time,digits=1)))
        end
    end
    # plot F(k)
    if mode == "c" plots["statprop.png"] = plot(xlabel = L"k_L",ylabel = L"F_k(t,t)",)
    else
        for (i, tidx) in enumerate(timeidx)
        plot!(plots["statprop.png"],[x.lev for x in disc.fftwhelper], measurearray[tidx].phi2, yerr=measurearray[tidx].phi2_err,ls=line, lc=color[i], markerstrokecolor=:auto, label=label * ": tm=" * string(round(measurearray[tidx].time,digits=1)))
        end
    end
    # dispersion relation
    if mode == "c" plots["dispersionrelation.png"] = plot( xlabel = L"k_L", ylabel = L"\epsilon_k")
    else
        for (i, tidx) in enumerate(timeidx)
        plot!(plots["dispersionrelation.png"],[x.lev for x in disc.fftwhelper], measurearray[tidx].omega, yerr=measurearray[tidx].omega_err,ls=line, lc=color[i], markerstrokecolor=:auto, label=label * ": tm=" * string(round(measurearray[tidx].time,digits=1)))
        end
    end
    # plot Log (1 + 1/np)
    if mode == "c" plots["BEdistribution.png"] = plot(xlabel = L"\epsilon_k",ylabel = L"Log[1+1/n_k]",)
    else
        for (i, tidx) in enumerate(timeidx)
        tmparray = (1 ./ measurearray[tidx].n) .+ 1
        idxlog = tmparray .> 0.001 # for 1/n + 1 values for which the log function works
        plot!(plots["BEdistribution.png"],measurearray[tidx].omega[idxlog], log.(tmparray[idxlog]), line=(line,color[i]), marker=(marker,color[i]), label=label * ": tm=" * string(round(measurearray[tidx].time,digits=1)))
        end
    end
    # plot np(omega)
    if mode == "c" plots["particledistribution.png"] = plot(xlabel = L"\epsilon_k",ylabel = L"n_k",)
    else
        for (i, tidx) in enumerate(timeidx)
        plot!(plots["particledistribution.png"],measurearray[tidx].omega, measurearray[tidx].n, line=(line,color[i]), marker=(marker,color[i]), label=label * ": tm=" * string(round(measurearray[tidx].time,digits=1)))
        end
    end
end