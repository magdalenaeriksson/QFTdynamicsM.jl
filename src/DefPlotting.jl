using Plots
# common GR plotting error (use env variable) with .png files: https://github.com/JuliaPlots/Plots.jl/issues/1612
# anyways: use .pdf, its better
using LaTeXStrings
using ColorSchemes

export displayplots
function displayplots(plots::Dict)
    for plotname in keys(plots)
        display(plots[plotname])
    end
end

export saveplots
function saveplots(plots::Dict, plotdir)
    for plotname in keys(plots)
        savefig( plots[plotname], joinpath(plotdir, plotname) ) # .png version for webview
        savefig( plots[plotname], joinpath(plotdir * "/pdf", plotname[1:end-4] * ".pdf") ) #.pdf version for publications
        println("Saved as " * joinpath(plotdir, plotname) ); flush(stdout)
    end
    # permissions
	cp("etc/index.php", joinpath(plotdir,"index.php"), force=true )
    println("Recusevly changing permissions in " * string(plotdir)); flush(stdout)
    chmod(plotdir, 0o777, recursive=true)
end