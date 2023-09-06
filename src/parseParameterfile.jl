using ArgParse

#"""
#this file allows to use Simulation.jl with an argument --Parameterfile=<theparameterfile.jl>
#"""

#
# Argparser
#
ap = ArgParseSettings()
@add_arg_table ap begin
    # Discretization & Measurement
    "--Parameterfile"
    help = ".jl file specifying the parameters of the simulation"
	arg_type = String
end
parsed_args = parse_args(ARGS, ap)

include(parsed_args["Parameterfile"])