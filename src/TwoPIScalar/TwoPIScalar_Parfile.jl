# specify which src/<Model> to use
modelfile = TwoPIScalarfile()

parameters = Dict(
    # Discretisation
    "Mass"  =>  0.7,
    "Nx"    =>  16, 
    "sdim"  =>  2, 
    "dt"    =>  0.1, #0.05
    "Nsteps"=>  300,
    "NstepsinMemory"  =>  0, # 3 - min for LO
    "Nmeas" =>  0, # either 0 (all timesteps are measured) or >=3 (because the plot script wants to draw 3 graphs)
    #model
    "Mod"  =>  "phi4",
    "Lambda"=>  6,
    "ONgroup"=> 8,
    #initialisation
    "init"  =>  "Pnr",
    #"init"  =>  "Thermal",
    #"init"  =>  "TopHatT1",
    #"init"  =>  "Quench",
    #"init"  =>  "Gauss",
    "xi"    => 0.4,
    "eta"   => 5.0,
    "sig"   => 1.0,
    "T"     =>  5,
    "n"     =>  2,
    #Reno
    "Reno"  =>  "RNone",
    #"Reno"  =>  "RMass",
    #PertExp
    #"Pexp"  =>  "LambdaLO",
    #"Pexp"  =>  "LambdaNLO",
    #"Pexp"  =>  "NinverseLO",
    "Pexp"  =>  "NinverseNLO",
    # Numerics
    #"Num"   =>  "CPUfull",
    "Num"   =>  "CPUred",
    #"Num"   =>  "CPUred2",
    "Nchunks"   =>  0, # if 0 -> set to @Threads.nthreads
    #extra
    "tag"   =>  "testing", #dont use "_" in the tag
    "ov"    =>  true)
