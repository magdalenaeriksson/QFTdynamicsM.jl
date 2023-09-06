# specify which src/<Model> to use
modelfile = GenericModelFile()

parameters = Dict(
    # Discretisation
    "Mass"  =>  0.7,
    "Nx"    =>  4, 
    "sdim"  =>  1, 
    "dt"    =>  0.1, #0.05
    "Nsteps"=>  100,
    "NstepsinMemory"  =>  0, # 3 - min for LO
    "Nmeas" =>  0,
    #model
    "Mod"  =>  "GenericModelA",
    "Omega"=>  1.5,
    #initialisation
    "init"  =>  "aninitialisation",
    #Reno
    "Reno"  =>  "arenormalisation",
    #PertExp
    "Pexp"  =>  "anexpansion",
    #extra
    "tag"   =>  "test",
    "ov"    =>  true)
