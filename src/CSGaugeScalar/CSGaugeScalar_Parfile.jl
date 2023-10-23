modelfile = CSGaugeScalarFile()
const REPL = true

    parameters = Dict(
        # Discretisation
        "Mass"  =>  1,
        "Nx"    =>  32, 
        "sdim"  =>  3, 
        "dt"    =>  0.05,
        "Nsteps"=>  200,
        "NstepsinMemory"  =>  0, # 3 - min for LO
        "Nmeas" =>  0, # either 0 (all timesteps are measured) or >=3 (because the plot script wants to draw 3 graphs)
        #model
        "Mod"  =>  "CS_SUNgaugeScalar",
        #"Mod"  =>  "CS_U1gaugeScalar",
        "Lambda" =>  0.125,
        "g" =>  0.002,#0.667,
        "N" =>  2,
        #"ONgroup"=> 1,
        #initialisation
        "init"  =>  "Pnr",
        #"init"  =>  "Thermal",
        #"init"  =>  "TopHatT1",
        "T"     =>  5,
        "n"     =>  1,
        #Reno
        "Reno"  =>  "RNone",
        #PertExp
        #"Pexp"  =>  "LambdaLO",
        #"Pexp"  =>  "LambdaNLO",
        #"Pexp"  =>  "LOloop",
        "Pexp"  =>  "CS",
        # Numerics
        #"Num"   =>  "GPUreduced",
        "Num"   =>  "CPU",
        "Runs"   =>  90,#90,
        "seed"   =>  1,
        "Nchunks"   =>  0, # if 0 -> set to @Threads.nthreads
        #extra
        "tag"   =>  "newInitalCondsTest", #dont use "_" in the tag
        "ov"    =>  true)