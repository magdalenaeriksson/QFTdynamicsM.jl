using QFTdynamics
using Serialization
using Parameters

function createrandomF(disc)
    if disc.sdim==1 A = ones(disc.Nx) end
    if disc.sdim==2 A = ones(disc.Nx,disc.Nx) end
    if disc.sdim==3 A = ones(disc.Nx,disc.Nx,disc.Nx) end
    for element in disc.fftwhelper
        aval = rand() - 0.5
        for idx in element.ind
            A[idx] = aval
        end
    end
    return A
end

function TestTwoPIScalar_LambdaNLO_full(Nx, sdim, Nchunks)
    #
    # compares 
    #
    # set up parameters
    #
    modelfile = TwoPIScalarfile()
    parameters = Dict(
        # Discretisation
        "Mass"  =>  1,
        "Nx"    =>  Nx, 
        "sdim"  =>  sdim, 
        "dt"    =>  0.01, #0.05
        "Nsteps"=>  12,
        "NstepsinMemory"  =>  0, # 3 - min for LO
        "Nmeas" =>  0, # either 0 (all timesteps are measured) or >=3 (because the plot script wants to draw 3 graphs)
        #model
        "Mod"  =>  "phi4",
        "Lambda"=>  10,
        "ONgroup"=> 8,
        #initialisation
        #"init"  =>  "Pnr",
        #"init"  =>  "Thermal",
        "init"  =>  "TopHatT1",
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
        "Pexp"  =>  "LambdaNLO",
        #"Pexp"  =>  "NinverseLO",
        #"Pexp"  =>  "NinverseNLO",
        # Numerics
        #"Num"   =>  "GPUreduced",
        "Num"   =>  "CPUfull",
        "Nchunks"   =>  Nchunks, # if 0 -> set to @Threads.nthreads
        #extra
        "tag"   =>  "testing", #dont use "_" in the tag
        "ov"    =>  true)

    #
    # create Problem
    #
    theproblem = QFTdynamicsProblem(modelfile, parameters)
    thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)
    initialize!(thesolution)

    #
    # evolve to some point
    #
    evto = 10
    println("Evolving to timestep:", evto)
    for i in 3:evto
    evolve!(thesolution, tmpdata, i)
    end

    #
    # check I functions at evto+1
    #
    t=evto+1
    println("Preparing for evolving to timestep:", t)
    @unpack problem, simdata, measurearray = thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem
    expandSimData!(simdata) 

    # calcSelfEnergies
    tmone = t-1
    println("Set random to F and r at timestep: ", tmone)
    # add random data to simdata and make deepcopy
    for tp in simdata.indices[1]:simdata.indices[2] # set it on every tp value, even if the last one actually overwrites other ones
        simdata.F[tmone,tp] .= createrandomF(disc) 
        simdata.r[tmone,tp] .= createrandomF(disc)
    end

    simdata2 = deepcopy(simdata)

    println("Calculate Self Energies at timestep: ", tmone, ", tp range: ", 1, " to ", tmone)
    println("Calculating...")
    # update tmpdata.threadranges for calcSelfEnergies! (it is parallized and works with threadranges inside)
    nrofindices = (simdata.indices[2]-1) - simdata.indices[1] + 1
    tmpdata.threadranges .= splitter(simdata.indices[1], nrofindices, num.nchunks)
    calcSelfEnergies!(model, pexp, disc, simdata, tmpdata, tmone)
    calcSelfEnergiessimple!(model, pexp, disc, simdata2, tmpdata, tmone)

    println("Comparing...")
    for tp in 1:(tmone)
        if sum(isapprox.(simdata.Sigr[tp],simdata2.Sigr[tp])) != disc.vol 
            @show tp
            @show simdata.Sigr[tp][1]
            @show simdata2.Sigr[tp][1]
            return false 
        end
        if sum(isapprox.(simdata.SigF[tp],simdata2.SigF[tp])) != disc.vol 
            @show simdata.SigF[tp]
            @show simdata2.SigF[tp]
            return false
        end
    end

    return true
end