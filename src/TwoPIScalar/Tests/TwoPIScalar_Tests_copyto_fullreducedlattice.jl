using QFTdynamics
using Serialization
using Parameters
using BenchmarkTools

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


function TestTwoPIScalar_copy_fullreduced(Nx, sdim)
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
        #"Pexp"  =>  "LambdaNLO",
        #"Pexp"  =>  "NinverseLO",
        "Pexp"  =>  "NinverseNLO",
        # Numerics
        #"Num"   =>  "GPUreduced",
        "Num"   =>  "CPUred",
        "Nchunks"   =>  1, # if 0 -> set to @Threads.nthreads
        #extra
        "tag"   =>  "testing", #dont use "_" in the tag
        "ov"    =>  true)

    #
    # create Problem
    #
    theproblem = QFTdynamicsProblem(modelfile, parameters)
    thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)
    initialize!(thesolution, tmpdata)

    @unpack problem, simdata, measurearray = thesolution
    @unpack model, pexp, disc, init, reno, simsetup = problem
    #
    # test copy to (full)lattice on F[2,2]
    #
    tmone = 2
    tp = 2
    ichunk = 1

    # copy to (full)lattice
    tmpdata2 = deepcopy(tmpdata)
    copytolattice!( tmpdata.F[ichunk], simdata.F[tmone,tp], disc.fftwhelper)
    copytofulllattice!( tmpdata2.F[ichunk], simdata.F[tmone,tp], tmpdata2)
    # benchmarking
    #Afull = deepcopy(tmpdata2.F[ichunk])
    #Ared = deepcopy(simdata.F[tmone,tp]) 
    #@benchmark copytolattice!( $Afull, $Ared, disc.fftwhelper) seconds=30
    #@benchmark copytofulllattice!(      $Afull, $Ared, $tmpdata2) seconds=30

    if sum(isapprox.(tmpdata.F[ichunk],tmpdata2.F[ichunk])) != disc.vol 
        println("Error in copy to lattice/copy to full lattice")
        #@show tmpdata2.F[ichunk]
        #@show tmpdata.F[ichunk]
        return false
    end

    # copy to reducedlattice
    simdata2 = deepcopy(simdata)
    copytoreducedlattice!( simdata.SigF[tp], tmpdata.F[ichunk], disc.fftwhelper)
    copytoreducedlattice!( simdata2.SigF[tp], tmpdata2.F[ichunk], tmpdata2)
    # benchmarking
    #Afull = deepcopy(tmpdata2.F[ichunk])
    #Ared = deepcopy(simdata.F[tmone,tp]) 
    #@benchmark copytoreducedlattice!(    $Ared, $Afull, $disc.fftwhelper)
    #@benchmark copytoreducedlattice!(    $Ared, $Afull, $tmpdata2)

    if sum(isapprox.(simdata.SigF[tp], simdata2.SigF[tp])) != length(disc.fftwhelper)
        println("Error in copy to reduced lattice")
        #@show simdata.SigF[tp]
        #@show simdata2.SigF[tp]
        return false
    end

    return true
end
