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

    # levels
    kvals = [el.lev for el in disc.fftwhelper]
    ####################################################################################################
    ## time plots
    ## 
    if mode == "c" plots["Phi2t.png"] = plot(xlabel = L"tm",ylabel =" ")#yerr=[element.Phi2k_err[1] for element in measurearray]
    else
        plot!(plots["Phi2t.png"],[element.time for element in measurearray], [element.Phi2k[1] for element in measurearray],    ls=line, lc=color[1], markerstrokecolor=:auto, label=L"\Phi^2[1]")
        plot!(plots["Phi2t.png"],[element.time for element in measurearray], [element.Phi2k[123] for element in measurearray],  ls=:dash, lc=color[2], markerstrokecolor=:auto, label=L"\Phi^2[123]")
        plot!(plots["Phi2t.png"],[element.time for element in measurearray], [element.Phi2k[500] for element in measurearray],  ls=:dashdot, lc=color[3], markerstrokecolor=:auto, label=L"\Phi^2[500]")
        plot!(plots["Phi2t.png"],[element.time for element in measurearray], [element.E2k[1] for element in measurearray],      ls=line, lc=color[4], markerstrokecolor=:auto, label=L"E^2[1]")
        plot!(plots["Phi2t.png"],[element.time for element in measurearray], [element.E2k[123] for element in measurearray],    ls=:dash, lc=color[5], markerstrokecolor=:auto, label=L"E^2[123]")
        plot!(plots["Phi2t.png"],[element.time for element in measurearray], [element.E2k[500] for element in measurearray],    ls=:dashdot, lc=color[6], markerstrokecolor=:auto, label=L"E^2[500]")
    end
    if mode == "c" plots["DLvsDT.png"] = plot(xlabel = L"tm",ylabel =" ")#yerr=[element.Phi2k_err[1] for element in measurearray]
    else
        plot!(plots["DLvsDT.png"],[element.time for element in measurearray], [element.DLk[1] for element in measurearray],     ls=line, lc=color[1], markerstrokecolor=:auto,    label="DL[1]")
        plot!(plots["DLvsDT.png"],[element.time for element in measurearray], [element.DTk[1] for element in measurearray],     ls=line, lc=color[2], markerstrokecolor=:auto,    label="DT[1]")
        plot!(plots["DLvsDT.png"],[element.time for element in measurearray], [element.DLk[123] for element in measurearray],   ls=:dash, lc=color[1], markerstrokecolor=:auto,   label="DL[123]")
        plot!(plots["DLvsDT.png"],[element.time for element in measurearray], [element.DTk[123] for element in measurearray],   ls=:dash, lc=color[2], markerstrokecolor=:auto,   label="DT[123]")
        plot!(plots["DLvsDT.png"],[element.time for element in measurearray], [element.DLk[500] for element in measurearray],   ls=:dashdot, lc=color[1], markerstrokecolor=:auto,   label="DL[500]")
        plot!(plots["DLvsDT.png"],[element.time for element in measurearray], [element.DTk[500] for element in measurearray],   ls=:dashdot, lc=color[2], markerstrokecolor=:auto,   label="DT[500]")
    end
    if mode == "c" plots["TvsL.png"] = plot(xlabel = L"tm",ylabel =" ")#yerr=[element.Phi2k_err[1] for element in measurearray]
    else
        plot!(plots["TvsL.png"],[element.time for element in measurearray], [element.EpropL[1] for element in measurearray],    ls=line, lc=color[1], markerstrokecolor=:auto,   label="EL[1]")
        plot!(plots["TvsL.png"],[element.time for element in measurearray], [element.EpropT[1] for element in measurearray],    ls=line, lc=color[2], markerstrokecolor=:auto,   label="ET[1]")
        plot!(plots["TvsL.png"],[element.time for element in measurearray], [element.DLk[1] for element in measurearray],       ls=line, lc=color[3], markerstrokecolor=:auto,   label="DL[1]")
        plot!(plots["TvsL.png"],[element.time for element in measurearray], [element.DTk[1] for element in measurearray],       ls=line, lc=color[4], markerstrokecolor=:auto,   label="DT[1]")
        plot!(plots["TvsL.png"],[element.time for element in measurearray], [element.EpropL[123] for element in measurearray],  ls=:dash, lc=color[1], markerstrokecolor=:auto,  label="EL[123]")
        plot!(plots["TvsL.png"],[element.time for element in measurearray], [element.EpropT[123] for element in measurearray],  ls=:dash, lc=color[2], markerstrokecolor=:auto,  label="ET[123]")
        plot!(plots["TvsL.png"],[element.time for element in measurearray], [element.DLk[123] for element in measurearray],     ls=:dash, lc=color[3], markerstrokecolor=:auto,  label="DL[123]")
        plot!(plots["TvsL.png"],[element.time for element in measurearray], [element.DTk[123] for element in measurearray],     ls=:dash, lc=color[4], markerstrokecolor=:auto,  label="DT[123]")
    end
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
    #####################################################################################################
    ### Gauss constraint plots
    ###
    #if mode == "c" plots["GaussTot.png"] = plot(xlabel = L"tm",ylabel ="Gauss constraint")
    #else
    #    plot!(plots["GaussTot.png"],[element.time for element in measurearray], [element.GaussTot for element in measurearray], ls=line, lc=color[1], markerstrokecolor=:auto, label=L"\sum_x (G^a_x)^2/V")
    #end
    #if mode == "c" plots["GaussRel.png"] = plot(xlabel = L"tm",ylabel ="Gauss constraint")
    #else
    #    plot!(plots["GaussRel.png"],[element.time for element in measurearray], [element.GaussRel for element in measurearray], ls=line, lc=color[2], markerstrokecolor=:auto, label=L"\sum_x(G^a_{x,\text{rel}})^2/V")
    #end
    #####################################################################################################
    ### Energy Plots
    ###
    # if mode == "c" plots["EnergyFractions.png"] = plot(xlabel = L"tm",ylabel ="Energy fractions")
    # else
    #    plot!(plots["EnergyFractions.png"],[element.time for element in measurearray], [element.Eelec/element.Etot for element in measurearray], ls=line, lc=color[1],  markerstrokecolor=:auto, label="elec")
    #    plot!(plots["EnergyFractions.png"],[element.time for element in measurearray], [element.Emagn/element.Etot for element in measurearray], ls=line, lc=color[2], markerstrokecolor=:auto, label="magn")
    #    #plot!(plots["EnergyFractions.png"],[element.time for element in measurearray], [element.EscalKin/element.Etot for element in measurearray], ls=line, lc=color[3], markerstrokecolor=:auto, label="scal kin")
    #    plot!(plots["EnergyFractions.png"],[element.time for element in measurearray], [element.EscalPot/element.Etot for element in measurearray], ls=line, lc=color[3], markerstrokecolor=:auto, label="scal pot")
    # end
    # if mode == "c" plots["EnergyTot.png"] = plot(xlabel = L"tm",ylabel ="Total energy")
    # else
    #    plot!(plots["EnergyTot.png"],[element.time for element in measurearray], [element.Etot for element in measurearray], ls=line,marker=(:x,2),markercolor=color[1], lc=color[1], markerstrokecolor=:auto, label="")
    # end
    # if mode == "c" plots["EnergyScalar.png"] = plot(xlabel = L"tm",ylabel ="Energy fractions")
    # else
    #    plot!(plots["EnergyScalar.png"],[element.time for element in measurearray], [element.Etot/element.Etot for element in measurearray], ls=line, lc=color[1], markerstrokecolor=:auto, label="total")
    #    plot!(plots["EnergyScalar.png"],[element.time for element in measurearray], [(element.EscalKin + element.EscalPot)/element.Etot for element in measurearray], ls=line, lc=color[2], markerstrokecolor=:auto, label="scal tot")
    #    plot!(plots["EnergyScalar.png"],[element.time for element in measurearray], [element.EscalKin/element.Etot for element in measurearray], ls=line, lc=color[3], markerstrokecolor=:auto, label="scal kin")
    #    #plot!(plots["EnergyScalar.png"],[element.time for element in measurearray], [10*element.EscalPot for element in measurearray], ls=line, lc=color[3], markerstrokecolor=:auto, label="10xpot")
    # end
    # if mode == "c" plots["EnergyGauge.png"] = plot(xlabel = L"tm",ylabel ="Energy components")
    # else
    #    plot!(plots["EnergyGauge.png"],[element.time for element in measurearray], [element.Eelec for element in measurearray], ls=line, lc=color[1], markerstrokecolor=:auto, label="elec")
    #    plot!(plots["EnergyGauge.png"],[element.time for element in measurearray], [element.Emagn for element in measurearray], ls=line, lc=color[2], markerstrokecolor=:auto, label="magn")
    #    plot!(plots["EnergyGauge.png"],[element.time for element in measurearray], [element.EscalPot for element in measurearray], ls=line, lc=color[3], markerstrokecolor=:auto, label="scal pot")
    # end
    #####################################################################################################
end