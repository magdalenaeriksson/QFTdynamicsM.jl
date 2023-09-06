using Plots
using LaTeXStrings
using ColorSchemes

export plotdata
function plotdata(plots::Dict, sol::QFTdynamicsSolutionGenericModel, mode::String, style::Int64=1, label::String="test")
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

    # plot measured quantity
    if mode == "c" plots["someMeasurement.png"] = plot(xlabel = L"tm",ylabel = L"sum(x_i(t))")
    else
        plot!(plots["someMeasurement.png"],[x.time for x in measurearray], [x.quantity for x in measurearray], line=(line...,color[1]), marker=(marker...,color[1]), label=label)
    end

end
