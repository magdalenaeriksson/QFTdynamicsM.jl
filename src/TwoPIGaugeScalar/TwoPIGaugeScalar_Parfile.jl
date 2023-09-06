# specify which src/<Model> to use
modelfile = TwoPIGaugeScalarfile()
const REPL = true

parameters = Dict(
    # Discretisation
    "Mass"  =>  1,
    "Nx"    =>  32, 
    "sdim"  =>  3, 
    "dt"    =>  0.05,
    "Nsteps"=>  10000,
    "NstepsinMemory"  =>  3, # 3 - min for LO (0 for all steps)
    "Nmeas" =>  5000, # either 0 (all timesteps are measured) or >=3 (because the plot script wants to draw 3 graphs)
    #model
    "Mod"  =>  "SUNgaugeScalar",
    #"Mod"  =>  "U1gaugeScalar",
    "Lambda" =>  1,#0.125,
    "g" =>  1,#0.667,
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
    "Pexp"  =>  "LOloop",
    #"Pexp"  =>  "NLOloop",
    # Numerics
    #"Num"   =>  "GPUreduced",
    "Num"   =>  "CPUfull",
    "Nchunks"   =>  0, # if 0 -> set to @Threads.nthreads
    #extra
    "tag"   =>  "test", #dont use "_" in the tag
    "ov"    =>  true)
