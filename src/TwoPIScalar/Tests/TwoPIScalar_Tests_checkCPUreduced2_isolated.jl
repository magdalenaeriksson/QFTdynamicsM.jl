using QFTdynamics
using Serialization
using Parameters

export createrandomF
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

function createrandomFred(disc)
    A = rand(length(disc.fftwhelper)) .- 0.5
    return A
end

Nx = 4
sdim = 2
Nchunks = 1
#function TestTwoPIScalar_NinverseNLO_red2(Nx, sdim, Nchunks)
    # elementwise comparison with isapprox

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
        #"Num"   =>  "GPUreduced",
        "Num"   =>  "CPUred2",
        "Nchunks"   =>  Nchunks, # if 0 -> set to @Threads.nthreads
        #extra
        "tag"   =>  "testing", #dont use "_" in the tag
        "ov"    =>  true)

    #
    # create Problem
    #
    theproblem = QFTdynamicsProblem(modelfile, parameters)
    thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)
    initialize!(thesolution, tmpdata)

    #
    # evolve to some point
    #
    evto = 10
    #evto = 3
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

    # add random data to simdata and make deepcopy
    println("Set random to F and r everywhere")
    for t in 3:evto
        for tp in simdata.indices[1]:t
            simdata.F[t,tp] .= createrandomFred(disc)
            simdata.r[t,tp] .= createrandomFred(disc)
        end
    end

    ichunk = 1
    for t in 1:evto
        for tp in 1:t
            calc_Fr_F2kr2!(t,tp, simdata, tmpdata, disc, ichunk)
            #@show (tp,t)
            #@show thesolution2.simdata.Fr[tp,t]
        end
    end
 
    simdata2 = deepcopy(simdata)
    tmpdata2 = deepcopy(tmpdata)
    # calcSelfEnergies
    tmone = t-1

    println("Calculate Ir/Is at timestep: ", tmone, ", tp range: ", 1, " to ", tmone)
    println("Calculating...")
    println("")
    # update tmpdata.threadranges for calcIrs! and calcIFs! (it is parallized and works with threadranges inside)
    nrofindices = (simdata.indices[2]-1) - simdata.indices[1] + 1
    tmpdata.threadranges .= splitter(simdata.indices[1], nrofindices, num.nchunks)

    calcIrs!(   model, pexp, disc, simdata,  tmpdata, tmone) # for all tp!
    calcIrsnew!(model, pexp, disc, simdata2, tmpdata2, tmone) # for all tp!

    
    println("CalcIFs...")
    calcIFs!(   model, pexp, disc, simdata,  tmpdata, tmone) # for all tp!
    println("CalcIFsnew...")
    calcIFsnew!(model, pexp, disc, simdata2, tmpdata2, tmone) # for all tp!

    println("Comparing...")
    println("")
    for tp in 1:tmone
        if sum(isapprox.(simdata.Irx[tp],simdata2.Irx[tp])) != length(simdata.Irx[1]) # number of elements in Irx is not volume!
            println("ERROR IN Ir")
            @show isapprox.(simdata.Irx[tp],simdata2.Irx[tp])
            @show tp
            @show simdata.Irx[tp]
            @show simdata2.Irx[tp]
#           return false
        end
        if sum(isapprox.(simdata.IFx[tp],simdata2.IFx[tp])) != length(simdata.IFx[1]) # number of elements in Irx is not volume!
            println("ERROR IN IF")
            @show isapprox.(simdata.IFx[tp],simdata2.IFx[tp])
            @show tp
            @show simdata.IFx[tp][1]
            @show simdata2.IFx[tp][1]
#           return false
        end
    end

#    return true
#end
