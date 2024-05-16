using Plots
using LaTeXStrings
using ColorSchemes

export plotdata
function plotdata(plots::Dict, thesolution::QFTdynamicsSolutionCSGaugeScalar, mode::String, style::Int64=1, label::String="test")
    # linestyle
    # line
    linealpha = 1
    linewidth = 1.5
    linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    # marker
    markersize = 2
    marker = (:o, markersize)
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
    timeidx = [1, Int64(floor(length(measurearray)/3) ), Int64(floor(length(measurearray)*2/3) ), length(measurearray)]
    #if length(measurearray) > 2 timeidx = [1, Int64(floor( length(measurearray)/2) ), length(measurearray)] end
    #if length(measurearray) > 4 timeidx = [1, Int64(floor( length(measurearray)/4) ), Int64(floor( length(measurearray)/2) ), Int64(floor( (3*length(measurearray))/4) ), length(measurearray)] end

    # momenta to plot
    nkLmomenta = 3
    kLmax = disc.fftwhelper[end].lev
    kLidx = zeros(Int64, nkLmomenta)
    for i in 0:(length(kLidx)-1)
        kLidx[i+1] = Int64(findmin(abs.([element.lev for element in disc.fftwhelper] .- ( (i/(nkLmomenta-1))*kLmax ) ))[2])
    end

    # levels
    kvals = [el.lev for el in disc.fftwhelper]

    ####################################################################################################
    # TIME EVOLUTION PLOTS
    ####################################################################################################
    # Gauss constraint plots
    #if mode == "c" plots["GaussTot.png"] = plot(xlabel = L"tm",ylabel ="Gauss constraint")
    #else
    #    plot!(plots["GaussTot.png"],[element.time for element in measurearray], [element.GaussTot for element in measurearray], ls=line, lc=color[1], markerstrokecolor=:auto, label=L"\sum_x (G^a_x)^2/V")
    #end
    #if mode == "c" plots["GaussRel.png"] = plot(xlabel = L"tm",ylabel ="Gauss constraint")
    #else
    #    plot!(plots["GaussRel.png"],[element.time for element in measurearray], [element.GaussRel for element in measurearray], ls=line, lc=color[2], markerstrokecolor=:auto, label=L"\sum_x(G^a_{x,\text{rel}})^2/V")
    #end
    ## Energy Plots
    #if mode == "c" plots["EnergyFractions.png"] = plot(xlabel = L"tm",ylabel ="Energy fractions")
    # else
    #    plot!(plots["EnergyFractions.png"],[element.time for element in measurearray], [element.Eelec/element.Etot for element in measurearray], ls=line, lc=color[1],  markerstrokecolor=:auto, label="elec")
    #    plot!(plots["EnergyFractions.png"],[element.time for element in measurearray], [element.Emagn/element.Etot for element in measurearray], ls=line, lc=color[2], markerstrokecolor=:auto, label="magn")
    #    plot!(plots["EnergyFractions.png"],[element.time for element in measurearray], [element.EscalKin/element.Etot for element in measurearray], ls=line, lc=color[3], markerstrokecolor=:auto, label="scal kin")
    #    plot!(plots["EnergyFractions.png"],[element.time for element in measurearray], [element.EscalPot/element.Etot for element in measurearray], ls=line, lc=color[3], markerstrokecolor=:auto, label="scal pot")
    # end
    # if mode == "c" plots["EnergyTot.png"] = plot(xlabel = L"tm",ylabel ="Total energy")
    # else
    #    plot!(plots["EnergyTot.png"],[element.time for element in measurearray], [element.Etot for element in measurearray], ls=line,marker=(:x,2),markercolor=color[1], lc=color[1], markerstrokecolor=:auto, label="")
    # end
    # if mode == "c" plots["EnergyScalar.png"] = plot(xlabel = L"tm",ylabel ="Energy fractions")
    # else
    #    plot!(plots["EnergyScalar.png"],[element.time for element in measurearray], [element.EscalPot for element in measurearray], ls=line, lc=color[1], markerstrokecolor=:auto, label="scal pot")
    #    plot!(plots["EnergyScalar.png"],[element.time for element in measurearray], [element.EscalKin for element in measurearray], ls=line, lc=color[2], markerstrokecolor=:auto, label="scal kin")
    # end
    # if mode == "c" plots["EnergyGauge.png"] = plot(xlabel = L"tm",ylabel ="Energy components")
    # else
    #    plot!(plots["EnergyGauge.png"],[element.time for element in measurearray], [element.Eelec for element in measurearray], ls=line, lc=color[1], markerstrokecolor=:auto, label="elec")
    #    plot!(plots["EnergyGauge.png"],[element.time for element in measurearray], [element.Emagn for element in measurearray], ls=line, lc=color[2], markerstrokecolor=:auto, label="magn")
    # end
    # Scalar field
    if mode == "c" plots["MeanPhi2_t.png"] = plot(xlabel = L"tm",ylabel =L"1/V \sum_x |\Phi(t;x)|^2")
    else
        plot!(plots["MeanPhi2_t.png"],[element.time for element in measurearray], [element.MeanPhix2 for element in measurearray],)
    end
    if mode == "c" plots["Phi2_t.png"] = plot(xlabel = L"tm",ylabel =L"<|\Phi(t;\mathbf{k})|^2>/V")
    else
        for (i, kidx) in enumerate(kLidx)
        plot!(plots["Phi2_t.png"],[element.time for element in measurearray], [element.Phi2k[kidx] for element in measurearray], ls=line, lc=color[i], label="k=" * string(round(disc.fftwhelper[kidx].lev, digits=1)))
        end
    end
    # Gauge field
    if mode == "c" plots["AGauge_t.png"] = plot(xlabel = L"tm",ylabel =L"<A(t,k)A(t,-k)>/V")
    else
        for (i, kidx) in enumerate(kLidx)
        plot!(plots["AGauge_t.png"],[element.time for element in measurearray], [element.DLk[kidx] for element in measurearray], ls=:solid, lc=color[i], label="DL: k=" * string(round(disc.fftwhelper[kidx].lev, digits=1)))
        plot!(plots["AGauge_t.png"],[element.time for element in measurearray], [element.DTk[kidx] for element in measurearray], ls=:dot  , lc=color[i], label="DT: k=" * string(round(disc.fftwhelper[kidx].lev, digits=1)))
        end
    end
    if mode == "c" plots["EGauge_t.png"] = plot(xlabel = L"tm",ylabel =L"<E(t,k)E(t,-k)>/V")
    else
        for (i, kidx) in enumerate(kLidx)
        plot!(plots["EGauge_t.png"],[element.time for element in measurearray], [element.E2L[kidx] for element in measurearray], ls=:solid, lc=color[i], label="E2L: k=" * string(round(disc.fftwhelper[kidx].lev, digits=1)))
        plot!(plots["EGauge_t.png"],[element.time for element in measurearray], [element.E2T[kidx] for element in measurearray], ls=:dot  , lc=color[i], label="E2T: k=" * string(round(disc.fftwhelper[kidx].lev, digits=1)))
        end
    end
    ####################################################################################################
    ## SNAPSHOT PLOTS
    ####################################################################################################
    # Scalar field
    if mode == "c" plots["Phi2.png"] = plot(xlabel = L"k_L", ylabel =L"<|\Phi(t;\mathbf{p}=0)|^2>/V")
    else
        for (i, tidx) in enumerate(timeidx)
        plot!( plots["Phi2.png"],kvals, measurearray[tidx].Phi2k, ls=line, marker=(:x,2), markercolor=color[i], lc=color[i], markerstrokecolor=:auto, label="tm=" * string(round(measurearray[tidx].time,digits=1)))
        end
    end
    # Gauge field
    if mode == "c" plots["ALGauge.png"] = plot(xlabel = L"k_L", ylabel =L"P  <A(t,k)A(t,-k)>/V") 
    else
        for (i, tidx) in enumerate(timeidx)
        scatter!( plots["ALGauge.png"],kvals, measurearray[tidx].DLk, marker=(:x,2),      markercolor=color[i], markerstrokecolor=color[i], label="DL: tm=" * string(round(measurearray[tidx].time,digits=1)))
        end
    end
    if mode == "c" plots["ATGauge.png"] = plot(xlabel = L"k_L", ylabel =L"P  <A(t,k)A(t,-k)>/V") 
    else
        for (i, tidx) in enumerate(timeidx)
        scatter!( plots["ATGauge.png"],kvals, measurearray[tidx].DTk, marker=(:circle,2), markercolor=color[i], markerstrokecolor=color[i], label="DT: tm=" * string(round(measurearray[tidx].time,digits=1)))
        end
    end
    if mode == "c" plots["ELGauge.png"] = plot(xlabel = L"k_L", ylabel =L"P  <E(t,k)E(t,-k)>/V") 
    else
        for (i, tidx) in enumerate(timeidx)
        scatter!( plots["ELGauge.png"],kvals, measurearray[tidx].E2L, marker=(:x,2),      markercolor=color[i], markerstrokecolor=color[i], label="E2L: tm=" * string(round(measurearray[tidx].time,digits=1)))
        end
    end
    if mode == "c" plots["ETGauge.png"] = plot(xlabel = L"k_L", ylabel =L"P  <E(t,k)E(t,-k)>/V") 
    else
        for (i, tidx) in enumerate(timeidx)
        scatter!( plots["ETGauge.png"],kvals, measurearray[tidx].E2T, marker=(:circle,2), markercolor=color[i], markerstrokecolor=color[i], label="E2T: tm=" * string(round(measurearray[tidx].time,digits=1)))
        end
    end
    ####################################################################################################
    ####################################################################################################
    ####################################################################################################
    #if mode == "c" plots["Phi2t.png"] = plot(xlabel = L"tm",ylabel =L"<|\Phi(t;\mathbf{p}=0)|^2>/V")#yerr=[element.Phi2k_err[1] for element in measurearray]
    #else
    #    plot!(plots["Phi2t.png"],[element.time for element in measurearray], [element.Phi2k[1] for element in measurearray], ls=line, marker=(:x,2),markercolor=color[2], lc=color[2], markerstrokecolor=:auto, label="mode 1")
    #    plot!(plots["Phi2t.png"],[element.time for element in measurearray], [element.Phi2k[123] for element in measurearray],  ls=line, lc=color[3], markerstrokecolor=:auto, label="mode 123")
    #    plot!(plots["Phi2t.png"],[element.time for element in measurearray], [element.Phi2k[500] for element in measurearray], ls=line, lc=color[4], markerstrokecolor=:auto, label="mode 500")
    #    plot!(plots["Phi2t.png"],[element.time for element in measurearray], [element.E2k[1] for element in measurearray], ls=line, marker=(:x,1),markercolor=color[1], lc=color[1], markerstrokecolor=:auto, label="Emode 1")
    #    plot!(plots["Phi2t.png"],[element.time for element in measurearray], [element.E2k[123] for element in measurearray],  ls=line, lc=color[5], markerstrokecolor=:auto, label="Emode 123")
    #    plot!(plots["Phi2t.png"],[element.time for element in measurearray], [element.E2k[500] for element in measurearray], ls=line, lc=color[6], markerstrokecolor=:auto, label="Emode 500")
    #end
    #if mode == "c" plots["E11.png"] = plot(xlabel = L"tm",ylabel =L"<|E_{11}(t;\mathbf{p})|^2>/V")#yerr=[element.Phi2k_err[1] for element in measurearray]
    #else
    #    plot!(plots["E11.png"],[element.time for element in measurearray], [element.EpropL[1] for element in measurearray], ls=line, marker=(:x,1),markercolor=color[1], lc=color[1], markerstrokecolor=:auto, label="Emode 1")
    #    plot!(plots["E11.png"],[element.time for element in measurearray], [element.EpropL[123] for element in measurearray],  ls=line, lc=color[5], markerstrokecolor=:auto, label="Emode 123")
    #    plot!(plots["E11.png"],[element.time for element in measurearray], [element.EpropL[500] for element in measurearray], ls=line, lc=color[6], markerstrokecolor=:auto, label="Emode 500")
    #end
    #if mode == "c" plots["Es.png"] = plot(xlabel = L"tm",ylabel =L"<|E_{11}(t;\mathbf{p})|^2>/V")#yerr=[element.Phi2k_err[1] for element in measurearray]
    #else
    #    plot!(plots["Es.png"],[element.time for element in measurearray], [element.E2k[1] for element in measurearray], ls=line,  lc=color[1], markerstrokecolor=:auto, label="E11[0]")
    #    plot!(plots["Es.png"],[element.time for element in measurearray], [element.DLk[1] for element in measurearray],  ls=line, lc=color[2], markerstrokecolor=:auto, label="E22[0]")
    #    plot!(plots["Es.png"],[element.time for element in measurearray], [element.DTk[1] for element in measurearray], ls=line, lc=color[3], markerstrokecolor=:auto,  label="E33[0]")
    #end
    # if mode == "c" plots["Phi2t2.png"] = plot(xlabel = L"tm",ylabel =L"<|\Phi(t;\mathbf{p}=0)|^2>/V")
    # else
    #     #plot!(plots["Phi2t2.png"],[element.time for element in measurearray], [element.Phi2k[1] for element in measurearray], ls=line, marker=(:x,2),markercolor=color[2], lc=color[2], markerstrokecolor=:auto, label="mode 1")
    #     plot!(plots["Phi2t2.png"],[element.time for element in measurearray], [element.Phi2k[123] for element in measurearray], yerr=[element.Phi2k_err[123] for element in measurearray], ls=line, lc=color[3], markerstrokecolor=:auto, label="mode 2")
    #     plot!(plots["Phi2t2.png"],[element.time for element in measurearray], [element.Phi2k[500] for element in measurearray], yerr=[element.Phi2k_err[500] for element in measurearray], ls=line, lc=color[4], markerstrokecolor=:auto, label="mode 3")
    # end
    # if mode == "c" plots["EGauge_t.png"] = plot(xlabel = L"tm",ylabel =L"<E(t,k)E(t,-k)>/V")
    # else
    #     #plot!(plots["Gauge_t.png"],[element.time for element in measurearray], [element.DLk[1] for element in measurearray], ls=line, lc=color[1], markerstrokecolor=:auto, label=L"D^2")
    #     plot!(plots["EGauge_t.png"],[element.time for element in measurearray], [element.EpropL[1] for element in measurearray], ls=line, lc=color[1], markerstrokecolor=:auto, label="E2L(k=0)")
    #     plot!(plots["EGauge_t.png"],[element.time for element in measurearray], [element.EpropT[1] for element in measurearray], ls=line, lc=color[2], markerstrokecolor=:auto, label="E2T(k=0)")
    # end
    # if mode == "c" plots["EGauge_t2.png"] = plot(xlabel = L"tm",ylabel =L"<E(t,k)E(t,-k)>/V")
    # else
    #     #plot!(plots["Gauge_t.png"],[element.time for element in measurearray], [element.DLk[1] for element in measurearray], ls=line, lc=color[1], markerstrokecolor=:auto, label=L"D^2")
    #     plot!(plots["EGauge_t2.png"],[element.time for element in measurearray], [element.EpropL[123] for element in measurearray], ls=line, lc=color[1], markerstrokecolor=:auto, label="E2L(n=123)")
    #     plot!(plots["EGauge_t2.png"],[element.time for element in measurearray], [element.EpropT[123] for element in measurearray], ls=line, lc=color[2], markerstrokecolor=:auto, label="E2T(n=123)")
    #     plot!(plots["EGauge_t2.png"],[element.time for element in measurearray], [element.EpropL[500] for element in measurearray], ls=line, lc=color[3], markerstrokecolor=:auto, label="E2L(n=500)")
    #     plot!(plots["EGauge_t2.png"],[element.time for element in measurearray], [element.EpropT[500] for element in measurearray], ls=line, lc=color[4], markerstrokecolor=:auto, label="E2T(n=500)")
    # end
    #if mode == "c" plots["AGauge_t1.png"] = plot(xlabel = L"tm",ylabel =L"Dij(t,k)=<Ai(t,0)Aj(t,0)>/V")
    #else
    #    plot!(plots["AGauge_t1.png"],[element.time for element in measurearray], [element.DLk[1] for element in measurearray], ls=line, lc=color[1], markerstrokecolor=:auto, label="DL(k=0)")
    #    plot!(plots["AGauge_t1.png"],[element.time for element in measurearray], [element.DTk[1] for element in measurearray], ls=line, lc=color[2], markerstrokecolor=:auto, label="DT(k=0)")
    #end
    #if mode == "c" plots["AGauge_t.png"] = plot(xlabel = L"tm",ylabel =L"Dij(t,k)=<Ai(t,k)Aj(t,-k)>/V")
    #else
    #    #plot!(plots["Gauge_t.png"],[element.time for element in measurearray], [element.DLk[1] for element in measurearray], ls=line, lc=color[1], markerstrokecolor=:auto, label=L"D^2")
    #    plot!(plots["AGauge_t.png"],[element.time for element in measurearray], [element.DLk[123] for element in measurearray], ls=line, lc=color[1], markerstrokecolor=:auto, label="DL(n=123)")
    #    plot!(plots["AGauge_t.png"],[element.time for element in measurearray], [element.DTk[123] for element in measurearray], ls=line, lc=color[2], markerstrokecolor=:auto, label="DT(n=123)")
    #    plot!(plots["AGauge_t.png"],[element.time for element in measurearray], [element.DLk[500] for element in measurearray], ls=line, lc=color[3], markerstrokecolor=:auto, label="DL(n=500)")
    #    plot!(plots["AGauge_t.png"],[element.time for element in measurearray], [element.DTk[500] for element in measurearray], ls=line, lc=color[4], markerstrokecolor=:auto, label="DT(n=500)")
    #end
    #####################################################################################################
    ## momentum plots
    ##
    # if mode == "c" plots["Phi2k.png"] = plot(xlabel = L"\mathbf{|p}|",ylabel =L"<|\Phi(t;\mathbf{|p}|)|^2>/V")#yerr=measurearray[1].Phi2k_err,
    # else
    #     #scatter!(plots["Phi2k.png"],kvals, measurearray[1].Phi2k, mc=:blue, ms=1, ma=0.3, markerstrokecolor=:auto, label=L"\Phi^2(k,t=0)")
    #     plot!(plots["Phi2k.png"],kvals, measurearray[1].Phi2k, yerr=measurearray[1].Phi2k_err, ls=line, lc=color[2], markerstrokecolor=:auto, label=L"t=t_0")
    #     plot!(plots["Phi2k.png"],kvals, measurearray[end].Phi2k, yerr=measurearray[end].Phi2k_err, ls=line, lc=color[3], markerstrokecolor=:auto, label=L"t=t_e")
    #     plot!(plots["Phi2k.png"],kvals, [0.5 / sqrt(k^2 + disc.Mass^2) for k in kvals], ls=line, lc=color[1], markerstrokecolor=:auto, label=L"\frac{1/2}{\omega_\mathbf{p}}")
    # end
    # if mode == "c" plots["EGauge_k.png"] = plot(xlabel = L"k",ylabel ="E-field props; <Ei(t,k)Ej(t,-k)>/V")
    # else
    #     #scatter!(plots["Phi2k.png"],kvals, measurearray[1].Phi2k, mc=:blue, ms=1, ma=0.3, markerstrokecolor=:auto, label=L"\Phi^2(k,t=0)")
    #     plot!(plots["EGauge_k.png"],kvals, measurearray[1].EpropL, ls=line, lc=color[1], markerstrokecolor=:auto, label="E2L(k,t=t_0)")
    #     plot!(plots["EGauge_k.png"],kvals, measurearray[1].EpropT, ls=line, lc=color[2], markerstrokecolor=:auto, label="E2Tsum(k,t=t_0)")
    #     plot!(plots["EGauge_k.png"],kvals, measurearray[end].EpropL, ls=line, lc=color[4], markerstrokecolor=:auto, label="E2L(k,t=t_e)")
    #     plot!(plots["EGauge_k.png"],kvals, measurearray[end].EpropT, ls=line, lc=color[3], markerstrokecolor=:auto, label="E2Tsum(k,t=t_e)")
    #     #plot!(plots["Gauge_k.png"],kvals, [0.5 / sqrt(k^2 + disc.Mass^2) for k in kvals], ls=line, lc=color[1], markerstrokecolor=:auto, label=L"\frac{1/2}{\omega_k}")
    # end
    # if mode == "c" plots["AGauge_k.png"] = plot(xlabel = L"k",ylabel ="Dij(t,k)=<Ai(t,k)Aj(t,-k)>/V")
    # else
    #     #scatter!(plots["Phi2k.png"],kvals, measurearray[1].Phi2k, mc=:blue, ms=1, ma=0.3, markerstrokecolor=:auto, label=L"\Phi^2(k,t=0)")
    #     plot!(plots["AGauge_k.png"],kvals, measurearray[1].DLk, ls=line, lc=color[1], markerstrokecolor=:auto, label="DL(k,t=t_0)")
    #     plot!(plots["AGauge_k.png"],kvals, measurearray[1].DTk, ls=line, lc=color[2], markerstrokecolor=:auto, label="DT(k,t=t_0)")
    #     plot!(plots["AGauge_k.png"],kvals, measurearray[end].DLk, ls=line, lc=color[4], markerstrokecolor=:auto, label="DL(k,t=t_e)")
    #     plot!(plots["AGauge_k.png"],kvals, measurearray[end].DTk, ls=line, lc=color[3], markerstrokecolor=:auto, label="DT(k,t=t_e)")
    #     #plot!(plots["Gauge_k.png"],kvals, [0.5 / sqrt(k^2 + disc.Mass^2) for k in kvals], ls=line, lc=color[1], markerstrokecolor=:auto, label=L"\frac{1/2}{\omega_k}")
    # end
   ####################################################################################################
end

export plotScalarcomponentdata
function plotScalarcomponentdata(plots::Dict, thesolution::QFTdynamicsSolutionCSGaugeScalar, style::Int64=1, label::String="test")
    # linestyle
    # line
    linealpha = 1
    linewidth = 1.5
    linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    # marker
    markersize = 2
    marker = (:o, markersize)
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
    timeidx = [1, Int64(floor( length(measurearray)/2) ), length(measurearray)]
    #if length(measurearray) > 2 timeidx = [1, Int64(floor( length(measurearray)/2) ), length(measurearray)] end
    #if length(measurearray) > 4 timeidx = [1, Int64(floor( length(measurearray)/4) ), Int64(floor( length(measurearray)/2) ), Int64(floor( (3*length(measurearray))/4) ), length(measurearray)] end

    # momenta to plot
    nkLmomenta = 3
    kLmax = disc.fftwhelper[end].lev
    kLidx = zeros(Int64, nkLmomenta)
    for i in 0:(length(kLidx)-1)
        kLidx[i+1] = Int64(findmin(abs.([element.lev for element in disc.fftwhelper] .- ( (i/(nkLmomenta-1))*kLmax ) ))[2])
    end

    # levels
    kvals = [el.lev for el in disc.fftwhelper]
    tvals = [element.time for element in measurearray]

    ####################################################################################################

    ####################################################################################################
    #plot
    l = @layout [a b ; c d]
    #theplot = plot(phiplot, phi2plot, piiplot, pii2plot, layout=l)

   ####################################################################################################
    ## INITIAL TIME PLOTS
    ####################################################################################################
    # plot F_0(t,t)
    plotvector = [plot(size=(1200,800), xlabel = L"k_L",ylabel = L"F_k(t,t)", ylim=(0,5)) for i in 1:4]
    for c in 1:4 # iterte through components
        plot!(plotvector[c], kvals, measurearray[1].phi2k[c], yerr=measurearray[1].phi2k_err[c], label="comp=" * string(c) )
    end
    plots["phi2comp_initial.png"] = plot(plotvector..., layout=l)
    # plot \partial_t \partial_t F_0(t,t)
    plotvector = [plot(size=(1200,800), xlabel = L"k_L",ylabel = L"\partial_t \partial_t F_k(t,t)", ylim=(0,5)) for i in 1:4]
    for c in 1:4 # iterte through components
        plot!(plotvector[c], kvals, measurearray[1].pi2k[c], yerr=measurearray[1].pi2k_err[c], label="comp=" * string(c) )
    end
    plots["pi2comp_initial.png"] = plot(plotvector..., layout=l)
    # plot n_k(t)
    plotvector = [plot(size=(1200,800), xlabel = L"k_L",ylabel = L"n_k", ylim=(0,5)) for i in 1:4]
    for c in 1:4 # iterte through components
        plot!(plotvector[c], kvals, measurearray[1].n[c], yerr=measurearray[1].n_err[c], label="comp=" * string(c) )
    end
    plots["particlenumber_initial.png"] = plot(plotvector..., layout=l)
   ####################################################################################################
    # TIME EVOLUTION PLOTS
    ####################################################################################################
    # plot F_0(t,t)
    plotvector = [plot(size=(1200,800), xlabel = L"tm",ylabel = L"F_k(t,t)",) for i in 1:4]
    for c in 1:4 # iterte through components
        for (i, kidx) in enumerate(kLidx)
            #plot!(plotvector[c], tvals, [element.phi2k[c][i] for element in measurearray], yerr=[element.phi2k_err[c][i] for element in measurearray], label="k=" * string(round(disc.fftwhelper[kidx].lev, digits=1)) )
            plot!(plotvector[c], tvals, [element.phi2k[c][kidx] for element in measurearray], label="k=" * string(round(disc.fftwhelper[kidx].lev, digits=1)) )
        end
    end
    plots["phi2comp_t.png"] = plot(plotvector..., layout=l)
    # plot n_k(t)
    plotvector = [plot(size=(1200,800), xlabel = L"tm",ylabel = L"n_k",) for i in 1:4]
    for c in 1:4 # iterte through components
        for (i, kidx) in enumerate(kLidx)
            #plot!(plotvector[c], tvals, [element.n[c][i] for element in measurearray], yerr=[element.n_err[c][i] for element in measurearray], label="k=" * string(round(disc.fftwhelper[kidx].lev, digits=1)) )
            plot!(plotvector[c], tvals, [element.n[c][kidx] for element in measurearray], label="k=" * string(round(disc.fftwhelper[kidx].lev, digits=1)) )
        end
    end
    plots["particlenumber_t.png"] = plot(plotvector..., layout=l)
    # plot F_0(t,t)
    plotvector = [plot(size=(1200,800), xlabel = L"tm",ylabel = L"F_k(t,t)",) for i in 1:4]
    for c in 1:4 # iterte through components
        for (i, kidx) in enumerate(kLidx)
            #plot!(plotvector[c], tvals, [element.phi2k[c][i] for element in measurearray], yerr=[element.phi2k_err[c][i] for element in measurearray], label="k=" * string(round(disc.fftwhelper[kidx].lev, digits=1)) )
            if i==1
                plot!(plotvector[c], tvals, [element.phi2k[c][2] for element in measurearray], label="k=" * string(round(disc.fftwhelper[2].lev, digits=1)) )
            else
                plot!(plotvector[c], tvals, [element.phi2k[c][kidx] for element in measurearray], label="k=" * string(round(disc.fftwhelper[kidx].lev, digits=1)) )
            end
        end
    end
    plots["phi2comp_nozero_t.png"] = plot(plotvector..., layout=l)
    # plot n_k(t)
    plotvector = [plot(size=(1200,800), xlabel = L"tm",ylabel = L"n_k",) for i in 1:4]
    for c in 1:4 # iterte through components
        for (i, kidx) in enumerate(kLidx)
            #plot!(plotvector[c], tvals, [element.n[c][i] for element in measurearray], yerr=[element.n_err[c][i] for element in measurearray], label="k=" * string(round(disc.fftwhelper[kidx].lev, digits=1)) )
            if i==1
                plot!(plotvector[c], tvals, [element.n[c][2] for element in measurearray], label="k=" * string(round(disc.fftwhelper[2].lev, digits=1)) )
            else
                plot!(plotvector[c], tvals, [element.n[c][kidx] for element in measurearray], label="k=" * string(round(disc.fftwhelper[kidx].lev, digits=1)) )
            end
        end
    end
    plots["particlenumber_nozero_t.png"] = plot(plotvector..., layout=l)
    ####################################################################################################
    ## SNAPSHOT PLOTS
    ####################################################################################################
    # plot F_0(t,t)
    plotvector = [plot(size=(1200,800), xlabel = L"k_L",ylabel = L"F_k(t,t)",) for i in 1:4]
    for c in 1:4 # iterte through components
        for (i, tidx) in enumerate(timeidx)
            #plot!(plotvector[c], kvals, measurearray[tidx].phi2k[c], yerr=measurearray[tidx].phi2k_err[c], label="tm=" * string(round(measurearray[tidx].time,digits=1)) )
            plot!(plotvector[c], kvals, measurearray[tidx].phi2k[c], label="tm=" * string(round(measurearray[tidx].time,digits=1)) )
        end
    end
    plots["phi2comp.png"] = plot(plotvector..., layout=l)
    # plot n_k(t)
    plotvector = [plot(size=(1200,800), xlabel = L"k_L",ylabel = L"n_k",) for i in 1:4]
    for c in 1:4 # iterte through components
        for (i, tidx) in enumerate(timeidx)
            #plot!(plotvector[c], kvals, measurearray[tidx].n[c], yerr=measurearray[tidx].n_err[c], label="tm=" * string(round(measurearray[tidx].time,digits=1)) )
            plot!(plotvector[c], kvals, measurearray[tidx].n[c], label="tm=" * string(round(measurearray[tidx].time,digits=1)) )
        end
    end
    plots["particlenumber.png"] = plot(plotvector..., layout=l)
    # dispersion relation
    plotvector = [plot(size=(1200,800), xlabel = L"k_L",ylabel = L"\omega_k",) for i in 1:4]
    for c in 1:4 # iterte through components
        for (i, tidx) in enumerate(timeidx)
            #plot!(plotvector[c], kvals, measurearray[tidx].omega[c], yerr=measurearray[tidx].omega_err[c], label="tm=" * string(round(measurearray[tidx].time,digits=1)) )
            plot!(plotvector[c], kvals, measurearray[tidx].omega[c], label="tm=" * string(round(measurearray[tidx].time,digits=1)) )
        end
    end
    plots["dispersionrelation.png"] = plot(plotvector..., layout=l)
    # plot Log (1 + 1/np)
    plotvector = [plot(size=(1200,800), xlabel = L"\omega_k",ylabel = L"ln[1+1/n_k]",) for i in 1:4]
    for c in 1:4 # iterte through components
        for (i, tidx) in enumerate(timeidx)
            tmparray = (1 ./ measurearray[tidx].n[c]) .+ 1
            idxlog = tmparray .> 0.001 # for 1/n + 1 values for which the log function works
            #plot!(plotvector[c], kvals, measurearray[tidx].n[c][idxlog], yerr=measurearray[tidx].n_err[c][idxlog], label="tm=" * string(round(measurearray[tidx].time,digits=1)) )
            idxlog[1] = false # dont print the 0 mode 
            #plot!(plotvector[c], kvals[idxlog], measurearray[tidx].n[c][idxlog], yerr=measurearray[tidx].n_err[c][idxlog], label="tm=" * string(round(measurearray[tidx].time,digits=1)) )
            plot!(plotvector[c], kvals[idxlog], measurearray[tidx].n[c][idxlog], label="tm=" * string(round(measurearray[tidx].time,digits=1)) )
        end
    end
    plots["BEdistribution.png"] = plot(plotvector..., layout=l)
    # plot F_0(t,t)
    plotvector = [plot(size=(1200,800), xlabel = L"k_L",ylabel = L"F_k(t,t)",) for i in 1:4]
    for c in 1:4 # iterte through components
        for (i, tidx) in enumerate(timeidx)
            #plot!(plotvector[c], kvals, measurearray[tidx].phi2k[c], yerr=measurearray[tidx].phi2k_err[c], label="tm=" * string(round(measurearray[tidx].time,digits=1)) )
            plot!(plotvector[c], kvals[2:end], measurearray[tidx].phi2k[c][2:end], label="tm=" * string(round(measurearray[tidx].time,digits=1)) )
        end
    end
    plots["phi2comp_nozero.png"] = plot(plotvector..., layout=l)
    # plot n_k(t)
    plotvector = [plot(size=(1200,800), xlabel = L"k_L",ylabel = L"n_k",) for i in 1:4]
    for c in 1:4 # iterte through components
        for (i, tidx) in enumerate(timeidx)
            #plot!(plotvector[c], kvals, measurearray[tidx].n[c], yerr=measurearray[tidx].n_err[c], label="tm=" * string(round(measurearray[tidx].time,digits=1)) )
            plot!(plotvector[c], kvals[2:end], measurearray[tidx].n[c][2:end], label="tm=" * string(round(measurearray[tidx].time,digits=1)) )
        end
    end
    plots["particlenumber_nozero.png"] = plot(plotvector..., layout=l)
end