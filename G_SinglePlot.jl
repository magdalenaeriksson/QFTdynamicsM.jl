using QFTdynamics
using Serialization
using ArgParse

###################################################################################################################################
# Set up
###################################################################################################################################
lines = readlines("etc/paths.txt")
datapath = chop(lines[2],head=5,tail=0)
postpath = chop(lines[3],head=5,tail=0)
plotpath = chop(lines[4],head=5,tail=0)
   
###################################################################################################################################
# Argparser
###################################################################################################################################
ap = ArgParseSettings()
@add_arg_table ap begin
    "--file"
end
parsed_args = parse_args(ARGS, ap)

###################################################################################################################################
# Do it
###################################################################################################################################
plots = Dict()

# make plots
#datatag1 = "phi4LambdaNLO_L600_ON1_RNone_Nx16_sdim2_Nsteps1000_dt100_M50_NstepsinMemory1002_testing"
datatag1 = parsed_args["file"]
thesolution = deserialize( datapath * "/" * datatag1 * "/SimSolution.jld"  )
plotdata(plots, thesolution, "c")  #create canvas
plotdata(plots, thesolution, "h", 1, "myplot") # fill canvas

# save plots
dirname = datatag1 * "_pdf"
mkpath(plotpath * "/" * dirname)
saveplots(plots, plotpath * "/" * dirname)

