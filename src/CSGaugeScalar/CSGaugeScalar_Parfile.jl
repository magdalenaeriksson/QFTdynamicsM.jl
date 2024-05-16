modelfile = CSGaugeScalarFile()
const REPL = true

    parameters = Dict(
        # Discretisation
        "Mass"  =>  1,
        "Nx"    =>  32, 
        "sdim"  =>  3, 
        "dt"    =>  0.05,
        "Nsteps"=>  50,
        "NstepsinMemory"  =>  0, # 3 - min for LO
        "Nmeas" =>  0, # either 0 (all timesteps are measured) or >=3 (because the plot script wants to draw 3 graphs)
        # Model (options: CS_SUNgaugeScalar, CS_SUNgaugeScalarTachyonic, CS_U1gaugeScalar)
        "Mod"  =>  "CS_SUNgaugeScalar",
        "Lambda" =>  0.125,
        "g" =>  1,#0.667,
        "N" =>  2,
        #"ONgroup"=> 1,
        # Initialisation (options: MCPnr, MCThermal, Pnr, Thermal, TopHatT1, TopHatT2, TopHatT3)
        "init"  =>  "MCPnr",
        "T"     =>  5,
        "n"     =>  1,
        # Renormalisation (not applicable for CS)
        "Reno"  =>  "RNone",
        # PertExp (options: CS, LOloop, LambdaLO, LambdaNLO)
        "Pexp"  =>  "CS",
        # Numerics
        #"Num"   =>  "GPUreduced",
        "Num"   =>  "CPU",
        "Runs"   =>  8, # should be >= num.Threads() = 8 by default
        "seed"   =>  1,
        "Nchunks"   =>  0, # if 0 -> set to @Threads.nthreads
        # Extras
        "tag"   =>  "newInitalCondsTest", #dont use "_" in the tag
        "ov"    =>  true)