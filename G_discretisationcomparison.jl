using FFTWhelper
using Plots

# lattice
Nx = 12
sdim = 3
masslat = 0.5

fftwhelper = getfftwhelper(Nx, sdim)
nrmomentalat = length(fftwhelper)
momentalat = [element.lev for element in fftwhelper] * masslat
deltakmax = fftwhelper[2].lev * masslat
deltakmin = (fftwhelper[end].lev - fftwhelper[end-1].lev) * masslat 

# continuum
masscon = (2/pi) * masslat

# cont1 - same nr of momenta
momentacont1 = LinRange(0.,maximum(momentalat),nrmomentalat)
# cont2 - match minimal spacing of lattice
nrmomentacon2 = maximum(momentalat) / deltakmin
momentacont2 = LinRange(0., maximum(momentalat), ceil(Int,nrmomentacon2))
# cont3 - match maximal spacing of lattice
nrmomentacon3 = maximum(momentalat) / deltakmax
momentacont3 = LinRange(0.,maximum(momentalat), ceil(Int, nrmomentacon3))

# plot
#theplot = plot(xlabel="|k|", yaxis=false, yticks=nothing)
fontsize = 10
theplot = plot(xlabel="|k|", yticks=nothing, xticks=nothing, xaxis=false, yaxis=false, draw_arrow = true, xlim=(0,maximum(momentalat)+0.05), ylim=(0.3,4.5), size=(800,250),margin=7Plots.mm , xtickfont=font(fontsize), legend = false, grid=false)
plot!(theplot, momentalat,      4 .+ zero(momentalat),   linewidth=2, markerstrokewidth=3, markersize=6, markershape=:vline, annotations = ([maximum(momentalat)], [4.14], text("3D lattice, #modes=" * string(length(momentalat)   ), fontsize, :right, :bottom  )    )  )
plot!(theplot, momentacont1,    3 .+ zero(momentacont1), linewidth=2, markerstrokewidth=3, markersize=6, markershape=:vline, annotations = ([maximum(momentalat)], [3.14], text("Continuum 1, #modes=" * string(length(momentacont1)), fontsize, :right, :bottom  )    )  )
plot!(theplot, momentacont2,    2 .+ zero(momentacont2), linewidth=2, markerstrokewidth=3, markersize=6, markershape=:vline, annotations = ([maximum(momentalat)], [2.14], text("Continuum 2, #modes=" * string(length(momentacont2)), fontsize, :right, :bottom  )    )  )
plot!(theplot, momentacont3,    1 .+ zero(momentacont3), linewidth=2, markerstrokewidth=3, markersize=6, markershape=:vline, annotations = ([maximum(momentalat)], [1.14], text("Continuum 3, #modes=" * string(length(momentacont3)), fontsize, :right, :bottom  )    )  )
plot!(theplot, [0,maximum(momentalat)+0.05], [0.5,0.5], arrow=true, linewidth=2, color=:black  )  
plot!(theplot, [0,0], [0.5,4.5], linewidth=2, color=:black  )  

savefig(theplot, "discretisationplot.pdf")
