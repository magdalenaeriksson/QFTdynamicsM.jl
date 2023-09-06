using Plots
using LaTeXStrings
using ColorSchemes

export plotdata
function plotdata(plots::Dict, sol::QFTdynamicsSolutionTwoPIGaugeScalar, mode::String, style::Int64=1, label::String="test")
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

    # plot EM(t)
    if mode == "c" plots["Energy.png"] = plot(xlabel = L"tm",ylabel = L"E_\mathrm{tot}")
    else
        plot!(plots["Energy.png"],[x.time for x in measurearray], disc.ivol*[x.E for x in measurearray], line=(line...,color[1]), marker=(marker...,color[1]), label= L"E_\mathrm{tot}")
    end
    if mode == "c" plots["EnergyComps.png"] = plot(xlabel = L"tm",ylabel = L"E")
    else
        plot!(plots["EnergyComps.png"], [x.time for x in measurearray], [x.E_kinTrans for x in measurearray], line=(line...,color[1]), marker=marker=(marker...,color[1]), label= L"E_\mathrm{kinTrans}")
        plot!(plots["EnergyComps.png"], [x.time for x in measurearray], [x.E_kinLong for x in measurearray], line=(line...,color[2]), marker=marker=(marker...,color[2]), label= L"E_\mathrm{kinLong}")
        plot!(plots["EnergyComps.png"], [x.time for x in measurearray], [x.E_kinHiggs for x in measurearray], line=(line...,color[3]), marker=marker=(marker...,color[3]), label= L"E_\mathrm{kinHiggs}")
        plot!(plots["EnergyComps.png"], [x.time for x in measurearray], [x.E_TrLogHiggs for x in measurearray], line=(line...,color[4]), marker=marker=(marker...,color[4]), label= L"E_\mathrm{TrLogHiggs}")
        plot!(plots["EnergyComps.png"], [x.time for x in measurearray], [x.E_TrLogGauge for x in measurearray], line=(line...,color[5]), marker=marker=(marker...,color[5]), label= L"E_\mathrm{TrLogGauge}")
        plot!(plots["EnergyComps.png"], [x.time for x in measurearray], [x.E_HiggsHartree for x in measurearray], line=(line...,color[6]), marker=marker=(marker...,color[6]), label= L"E_\mathrm{HiggsHartee}")
        plot!(plots["EnergyComps.png"], [x.time for x in measurearray], [x.E_MixedHartree for x in measurearray], line=(line...,color[7]), marker=marker=(marker...,color[7]), label= L"E_\mathrm{MixedHartree}")
        plot!(plots["EnergyComps.png"], [x.time for x in measurearray], [x.E_GaugeHartree for x in measurearray], line=(line...,color[6]), marker=marker=(marker...,color[7]), label= L"E_\mathrm{GaugeHartree}")
        #plot!(plots["EnergyComps.png"], [x.time for x in measurearray], [x.E_MixedSunset for x in measurearray], line=(line...,color[6]), marker=marker=(marker...,color[6]), label= L"E_\mathrm{MixedSunset}")
    end
    #if mode == "c" plots["EnergyComps2.png"] = plot(xlabel = L"tm",ylabel = L"E")
    #else
    #    plot!(plots["EnergyComps2.png"], [x.time for x in measurearray], [x.E_MixedSunset for x in measurearray], line=(line...,color[6]), marker=marker=(marker...,color[6]), label= L"E_\mathrm{MixedSunset}")
    #end
    
    # Equal-time statistical propagators (zero mode)
    if mode == "c" plots["stat0props.png"] = plot(xlabel = L"tm",ylabel = "")
    else
        plot!(plots["stat0props.png"],[element.time for element in measurearray], [element.FS[1] for element in measurearray], line=(line...,color[1]), marker=(marker...,color[1]), label=L"G_F(t,t;\mathbf{p}=0)")
        #plot!(plots["stat0props.png"],[element.time for element in measurearray], [element.FT[1] for element in measurearray], line=(line...,color[2]), marker=(marker...,color[2]), label=L"D^\perp_F(t,t;\mathbf{p}=0)")
        #plot!(plots["stat0props.png"],[element.time for element in measurearray], [element.FL[1] for element in measurearray], line=(line...,color[3]), marker=(marker...,color[3]), label=L"D^\parallel_F(t,t;\mathbf{p}=0)")
    end

    if mode == "c" plots["Gaugestat0props.png"] = plot(xlabel = L"tm",ylabel = "")
    else
        #plot!(plots["Gaugestat0props.png"],[element.time for element in measurearray], [element.FS[1] for element in measurearray], line=(line...,color[1]), marker=(marker...,color[1]), label=L"G_F(t,t;\mathbf{p}=0)")
        plot!(plots["Gaugestat0props.png"],[element.time for element in measurearray], [element.FT[1] for element in measurearray], line=(line...,color[2]), marker=(marker...,color[2]), label=L"D^\perp_F(t,t;\mathbf{p}=0)")
        plot!(plots["Gaugestat0props.png"],[element.time for element in measurearray], [element.FL[1] for element in measurearray], line=(line...,color[3]), marker=(marker...,color[3]), label=L"D^\parallel_F(t,t;\mathbf{p}=0)")
    end

    ### Unequal-time statistical propagators
    #if mode == "c" plots["offdiagF.png"] = plot(xlabel = L"(t-t^\prime)m", ylabel = L"F(t,t^\prime;\mathbf{p}=0)")
    #else
    #  plot!(plots["offdiagF.png"], [(i-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].FSonet, line=(line...,color[1]), marker=(marker...,color[1]), label=L"G_F(t,t';\mathbf{p}=0)")
    #  plot!(plots["offdiagF.png"], [(j-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].FTonet, line=(line...,color[2]), marker=(marker...,color[2]), label=L"D^\perp_F(t,t';\mathbf{p}=0)")
    #  plot!(plots["offdiagF.png"], [(j-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].FLonet, line=(line...,color[3]), marker=(marker...,color[3]), label=L"D^\parallel_F(t,t';\mathbf{p}=0)")
    #end
#
    ### Unequal-time spectral propagators
    #if mode == "c" plots["offdiagr.png"] = plot(xlabel = L"(t-t^\prime)m", ylabel = L"\rho(t,t^\prime;\mathbf{p}=0)")
    #else
    #  plot!(plots["offdiagr.png"], [(i-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].rSonet, line=(line...,color[1]), marker=(marker...,color[1]), label=L"G_\rho(t,t';\mathbf{p}=0)")
    #  plot!(plots["offdiagr.png"], [(i-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].rTonet, line=(line...,color[2]), marker=(marker...,color[2]), label=L"D^\perp_\rho(t,t';\mathbf{p}=0)")
    #  plot!(plots["offdiagr.png"], [(i-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].rLonet, line=(line...,color[3]), marker=(marker...,color[3]), label=L"D^\parallel_\rho(t,t';\mathbf{p}=0)")
    #end

    ## Masses
    if mode == "c" plots["Hartreemasses.png"] = plot(xlabel = L"tm",ylabel = L"M^2(t)")
    else
        plot!(plots["Hartreemasses.png"],[element.time for element in measurearray], [element.MS2 for element in measurearray],    line=(line...,color[1]), marker=(marker...,color[1]), label=L"M_H^2")
        #plot!(plots["Hartreemasses.png"],[element.time for element in measurearray], [element.MT2[1] for element in measurearray], line=(line...,color[2]), marker=(marker...,color[2]), label=L"M_\perp^2(\mathbf{p}=0)")
        #plot!(plots["Hartreemasses.png"],[element.time for element in measurearray], [element.ML2[1] for element in measurearray], line=(line...,color[3]), marker=(marker...,color[3]), label=L"M_\parallel^2(\mathbf{p}=0)")
    end
    if mode == "c" plots["GaugeHartreemasses.png"] = plot(xlabel = L"tm",ylabel = L"M^2(t)")
    else
        plot!(plots["GaugeHartreemasses.png"],[element.time for element in measurearray], [element.MT2[1] for element in measurearray], line=(line...,color[2]), marker=(marker...,color[2]), label=L"M_\perp^2(\mathbf{p}=0)")
        plot!(plots["GaugeHartreemasses.png"],[element.time for element in measurearray], [element.ML2[1] for element in measurearray], line=(line...,color[3]), marker=(marker...,color[3]), label=L"M_\parallel^2(\mathbf{p}=0)")
    end

    ### Unequal-time statistical self-energies
    #if mode == "c" plots["StatSelfEnergy.png"] = plot(xlabel = L"(t-t^\prime)m" ,ylabel = " ")
    #else
    #   plot!(plots["StatSelfEnergy.png"], [i * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].SigFS, line=(line...,color[1]), marker=(marker...,color[1]), label=L"\Sigma_F(t,t';\mathbf{p}=0)")
    #   plot!(plots["StatSelfEnergy.png"], [i * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].SigFT, line=(line...,color[2]), marker=(marker...,color[2]), label=L"\Pi^\perp_F(t,t';\mathbf{p}=0)")
    #   plot!(plots["StatSelfEnergy.png"], [i * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].SigFL, line=(line...,color[3]), marker=(marker...,color[3]), label=L"\Pi^\parallel_F(t,t';\mathbf{p}=0)")
    #end  
#
    ### Unequal-time spectral self-energies
    #if mode == "c" plots["SpecSelfEnergy.png"] = plot(xlabel = L"(t-t^\prime)m" ,ylabel = " ")
    #else
    #   plot!(plots["SpecSelfEnergy.png"], [i * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].SigrS, line=(line...,color[1]), marker=(marker...,color[1]), label=L"\Sigma_\rho(t,t';\mathbf{p}=0)")
    #   plot!(plots["SpecSelfEnergy.png"], [i * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].SigrT, line=(line...,color[2]), marker=(marker...,color[2]), label=L"\Pi^\perp_\rho(t,t';\mathbf{p}=0)")
    #   plot!(plots["SpecSelfEnergy.png"], [i * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].SigrL, line=(line...,color[3]), marker=(marker...,color[3]), label=L"\Pi^\parallel_\rho(t,t';\mathbf{p}=0)")
    #end








    # if mode == "c" plots["offdiagScalar.png"] = plot(xlabel = L"(t-t^\prime)m", ylabel = "")
    # else
    #    plot!(plots["offdiagScalar.png"], [(j-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].FSonet, line=(line...,color[2]), marker=(marker...,color[2]), label=L"G_F(t,t';\mathbf{p}=0)")
    #    plot!(plots["offdiagScalar.png"], [(j-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].rSonet, line=(line...,color[1]), marker=(marker...,color[1]), label=L"G_\rho(t,t';\mathbf{p}=0)")
    #    #plot!(plots["offdiagF.png"], [(i-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].FTonet, line=(line...,color[2]), marker=(marker...,color[2]), label="trans")
    #    #plot!(plots["offdiagF.png"], [(i-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].FLonet, line=(line...,color[3]), marker=(marker...,color[3]), label="long")
    # end
    # if mode == "c" plots["offdiagGauge.png"] = plot(xlabel = L"(t-t^\prime)m", ylabel = "")
    # else
    #    plot!(plots["offdiagGauge.png"], [(j-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].FTonet, line=(line...,color[1]), marker=(marker...,color[1]), label=L"D^\perp_F(t,t';\mathbf{p}=0)")
    #    plot!(plots["offdiagGauge.png"], [(j-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].rTonet, line=(line...,color[2]), marker=(marker...,color[2]), label=L"D^\perp_\rho(t,t';\mathbf{p}=0)")
    #    plot!(plots["offdiagGauge.png"], [(j-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].FLonet, line=(line...,color[3]), marker=(marker...,color[3]), label=L"D^\parallel_F(t,t';\mathbf{p}=0)")
    #    plot!(plots["offdiagGauge.png"], [(j-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ], measurearray[end].rLonet, line=(line...,color[4]), marker=(marker...,color[4]), label=L"D^\parallel_\rho(t,t';\mathbf{p}=0)")
    # end
    
    # plot c(k)
    #if mode == "c" plots["cnumber.png"] = plot(xlabel = L"k_L",ylabel = L"c_k",)
    #else
    #    #for i in 1:simsetup.Nmeas+1
    #    for (i, time) in enumerate( [measurearray[1].time, measurearray[Int64(floor(end/2))].time, measurearray[end].time])
    #    plot!(plots["cnumber.png"],[x.lev for x in disc.fftwhelper], measurearray[i].c, line=(line...,color[i]), marker=(marker...,color[i]), label=label * ": tm=" * string(time))
    #    end
    #end
    # plot n(k)
    #if mode == "c" plots["particlenumber.png"] = plot(xlabel = L"k_L",ylabel = L"n_k")
    #else
    #    #for i in 1:simsetup.Nmeas+1
    #    for (i, time) in enumerate( [measurearray[1].time, measurearray[Int64(floor(end/2))].time, measurearray[end].time])
    #        plot!(plots["particlenumber.png"],[element.lev for element in disc.fftwhelper], measurearray[i].nS, line=(line...,color[i]), marker=(marker...,color[1]), label="scalar: tm=" * string(round(time,digits=2)))
    #        plot!(plots["particlenumber.png"],[element.lev for element in disc.fftwhelper], measurearray[i].nT, line=(line...,color[i]), marker=(marker...,color[2]), label="trans: tm=" * string(round(time,digits=2)))
    #        plot!(plots["particlenumber.png"],[element.lev for element in disc.fftwhelper], measurearray[i].nL, line=(line...,color[i]), marker=(marker...,color[3]), label="long: tm=" * string(round(time,digits=2)))
    #    end
    #end
    # plot r(k)
    #if mode == "c" plots["specprop.png"] = plot(xlabel = L"k_L",ylabel = L"r_k(t,t)",)
    #else
    #    #for i in 1:simsetup.Nmeas+1
    #    for (i, time) in enumerate( [measurearray[1].time, measurearray[Int64(floor(end/2))].time, measurearray[end].time])
    #       plot!(plots["specprop.png"],[element.lev for element in disc.fftwhelper], measurearray[i].rS, line=(line...,color[i]), marker=(marker...,color[i]), label="scalar: tm=" * string(time))
    #    end
    #end
    # plot F(k)
    # if mode == "c" plots["statprop.png"] = plot(xlabel = L"k_L",ylabel = L"F(t,t;\mathbf{p}=0)")
    # else
    #     ##for i in 1:simsetup.Nmeas+1
    #     #for (i, time) in enumerate( [measurearray[1].time, measurearray[Int64(floor(end/2))].time, measurearray[end].time])
    #     #    plot!(plots["statprop.png"],[element.lev for element in disc.fftwhelper], measurearray[i].FS, line=(line...,color[i]), marker=(marker...,color[1]), label="scalar: tm=" * string(round(time,digits=2)))
    #     #    plot!(plots["statprop.png"],[element.lev for element in disc.fftwhelper], measurearray[i].FT, line=(line...,color[i]), marker=(marker...,color[2]), label="trans: tm=" * string(round(time,digits=2)))
    #     #    plot!(plots["statprop.png"],[element.lev for element in disc.fftwhelper], measurearray[i].FL, line=(line...,color[i]), marker=(marker...,color[3]), label="long: tm=" * string(round(time,digits=2)))
    #     #end
    #     # t=0 meas
    #     plot!(plots["statprop.png"],[element.lev for element in disc.fftwhelper], measurearray[1].FS, line=(line...,color[1]), marker=(marker...,color[1]), label="scalar: tm=0")
    #     plot!(plots["statprop.png"],[element.lev for element in disc.fftwhelper], measurearray[1].FT, line=(line...,color[2]), marker=(marker...,color[2]), label="trans: tm=0")
    #     plot!(plots["statprop.png"],[element.lev for element in disc.fftwhelper], measurearray[1].FL, line=(line...,color[3]), marker=(marker...,color[3]), label="long: tm=0")
    # end
    

    #############################################
    # plot dispersion relation
    #if mode == "c" plots["dispersionrelation.png"] = plot( xlabel = L"k_L", ylabel = L"\omega_k")
    #else
    #    #plot!(plots["dispersionrelation.png"],[element.lev for element in disc.fftwhelper], measurearray[1].omegaS, line=(line...,color[1]), marker=(marker...,color[1]), label="scalar: tm=" * string(measurearray[end].time))
    #    #plot!(plots["dispersionrelation.png"],[element.lev for element in disc.fftwhelper], measurearray[1].omegaT, line=(line...,color[2]), marker=(marker...,color[2]), label="trans: tm=" * string(measurearray[end].time))
    #    #plot!(plots["dispersionrelation.png"],[element.lev for element in disc.fftwhelper], measurearray[1].omegaL, line=(line...,color[3]), marker=(marker...,color[3]), label="long: tm=" * string(measurearray[end].time))
    #    ##for i in 1:simsetup.Nmeas+1
    #    for (i, time) in enumerate( [measurearray[1].time, measurearray[Int64(floor(end/2))].time, measurearray[end].time])
    #        plot!(plots["dispersionrelation.png"],[element.lev for element in disc.fftwhelper], measurearray[i].omegaS, line=(line...,color[i]), marker=(marker...,color[1]), label="scalar: tm=" * string(round(time,digits=2)))
    #        plot!(plots["dispersionrelation.png"],[element.lev for element in disc.fftwhelper], measurearray[i].omegaT, line=(line...,color[i]), marker=(marker...,color[2]), label="trans: tm=" * string(round(time,digits=2)))
    #        plot!(plots["dispersionrelation.png"],[element.lev for element in disc.fftwhelper], measurearray[i].omegaL, line=(line...,color[i]), marker=(marker...,color[3]), label="long: tm=" * string(round(time,digits=2)))
    #    end
    #end
    # plot Log (1 + 1/np)
    #if mode == "c" plots["BEdistribution.png"] = plot(xlabel = L"\epsilon_k",ylabel = L"Log[1+1/n_k(t)]",)
    #else
    #    #for i in 1:simsetup.Nmeas+1
    #    for (i, time) in enumerate( [measurearray[1].time, measurearray[Int64(floor(end/2))].time, measurearray[end].time])
    #        plot!(plots["BEdistribution.png"],measurearray[i].omega, log.(1 ./ measurearray[i].n .+ 1 ), line=(line...,color[i]), marker=(marker...,color[i]), label= label * ": tm=" * string(time))
    #    end
    #end
    # plot np(omega)
    #if mode == "c" plots["particledistribution.png"] = plot(xlabel = L"\epsilon_k",ylabel = L"n_k(t)",)
    #else
    #    #for i in 1:simsetup.Nmeas+1
    #    for (i, time) in enumerate( [measurearray[1].time, measurearray[Int64(floor(end/2))].time, measurearray[end].time])
    #        plot!(plots["particledistribution.png"],measurearray[i].omega, measurearray[i].n, line=(line...,color[i]), marker=(marker...,color[i]), label= label * ": tm=" * string(time))
    #    end
    #end

    # plot Memory
    #if mode == "c" plots["SigmaF.png"] = plot(xlabel = L"t-t^\prime" ,ylabel = L"\Sigma^F(t,t';k=0)")
    #else
    #    plot!(plots["SigmaF.png"], measurearray[simsetup.Nmeas+1].SigFS, line=(line...,color[1]), marker=(marker...,color[1]), label="scalar")
    #    plot!(plots["SigmaF.png"], measurearray[simsetup.Nmeas+1].SigFT, line=(line...,color[2]), marker=(marker...,color[2]), label="trans")
    #    plot!(plots["SigmaF.png"], measurearray[simsetup.Nmeas+1].SigFL, line=(line...,color[3]), marker=(marker...,color[3]), label="long")
    #end
    #if mode == "c" plots["Sigmar.png"] = plot(xlabel = L"t-t^\prime" ,ylabel = L"\Sigma^{\rho}(t,t';k=0)")
    #else
    #    plot!(plots["Sigmar.png"], measurearray[simsetup.Nmeas+1].SigrS, line=(line...,color[1]), marker=(marker...,color[1]), label= "scalar" )
    #    plot!(plots["Sigmar.png"], measurearray[simsetup.Nmeas+1].SigrT, line=(line...,color[2]), marker=(marker...,color[2]), label= "trans" )
    #    plot!(plots["Sigmar.png"], measurearray[simsetup.Nmeas+1].SigrL, line=(line...,color[3]), marker=(marker...,color[3]), label= "long" )
    #end
    ### plot F_0(t,1)
    #@show [(i-1) * disc.dt * disc.Mass for (i,j) in enumerate(simdata.indices[1]:(simdata.indices[2]-1)) ]
    ## plot F_0(t,1)
    # if mode == "c" plots["offdiagF.png"] = plot(xlabel = "PAST - index - FUTURE",ylabel = L"F(t,1;k=0)")
    # else
    #    plot!(plots["offdiagF.png"],[element.FSonet for element in measurearray], line=(line...,color[1]), marker=(marker...,color[1]), label="scalar")
    #    plot!(plots["offdiagF.png"],[element.FTonet for element in measurearray], line=(line...,color[2]), marker=(marker...,color[2]), label="trans")
    #    plot!(plots["offdiagF.png"],[element.FLonet for element in measurearray], line=(line...,color[3]), marker=(marker...,color[3]), label="long")
    # end
    ## plot r_0(t,1)
    #if mode == "c" plots["offdiagr.png"] = plot(xlabel = "PAST - index - FUTURE",ylabel = L"\rho(t,1,k=0)")
    #else
    #    plot!(plots["offdiagr.png"],[element.rSonet for element in measurearray], line=(line...,color[1]), marker=(marker...,color[1]), label="scalar")
    #    plot!(plots["offdiagr.png"],[element.rTonet for element in measurearray], line=(line...,color[2]), marker=(marker...,color[2]), label="trans")
    #    plot!(plots["offdiagr.png"],[element.rLonet for element in measurearray], line=(line...,color[3]), marker=(marker...,color[3]), label="long")
    #end

end
