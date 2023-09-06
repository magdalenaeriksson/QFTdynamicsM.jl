using QFTdynamics
using Serialization
using Parameters
using FFTWhelper
using LinearAlgebra

include("../TwoPIScalar_Testhelpers.jl")

function TestTwoPIScalar_CPUred2(Nx, sdim, Nchunks, Pexpstring, evolveto)
    #
    # set up parameters
    #
    parameters = Dict(
        # Discretisation
        "Mass"  =>  1,
        "Nx"    =>  Nx, 
        "sdim"  =>  sdim, 
        "dt"    =>  0.01, #0.05
        "Nsteps"=>  evolveto+10,
        "NstepsinMemory"  =>  0, # 3 - min for LO
        "Nmeas" =>  0, # either 0 (all timesteps are measured) or >=3 (because the plot script wants to draw 3 graphs)
        #model
        "Mod"  =>  "phi4",
        "Lambda"=>  10,
        "ONgroup"=> 8,
        #initialisation
        "init"  =>  "Thermal",
        "xi"    => 0.4,
        "eta"   => 5.0,
        "sig"   => 1.0,
        "T"     =>  5,
        "n"     =>  2,
        #Reno
        "Reno"  =>  "RNone",
        #"Reno"  =>  "RMass",
        #PertExp
        "Pexp"  =>  Pexpstring,
        # Numerics
        "Num"   =>  "CPUred",
        "Nchunks"   =>  Nchunks, # if 0 -> set to @Threads.nthreads
        #extra
        "tag"   =>  "testing", #dont use "_" in the tag
        "ov"    =>  true)

    parameters["Num"] = "CPUred"
    parameters["NstepsinMemory"] = 0
    parameters["Nmeas"] = 0
    evolveto = 2 # dont change!!
    (thesolution, tmpdata) = TwoPIScalar_testevolution_returnall(TwoPIScalarfile(), parameters, evolveto)
    #@show thesolution.simdata.indices
    #@show thesolution.simdata.F[1,1]
    #@show thesolution.simdata.F[2,2]
    #@show thesolution.simdata.F[3,3]
    #@show thesolution.simdata.F[4,4]

    parameters["Num"] = "CPUred2"
    parameters["NstepsinMemory"] = 0
    parameters["Nmeas"] = 0
    (thesolution2, tmpdata2) = TwoPIScalar_testevolution_returnall(TwoPIScalarfile(), parameters, 0)
    for i in 1:2 expandSimData!(thesolution2.simdata) end

    if sum( isapprox.(thesolution.simdata.indices, thesolution2.simdata.indices) ) != 2 println("simdatas have different indeices - adjust expand") end 

    # copy r and F to thesolution2 
    for t in 1:(evolveto+2)
        for tp in 1:t
            @show (tp,t)
            thesolution2.simdata.F[t,tp] .= thesolution.simdata.F[t,tp]
            thesolution2.simdata.r[t,tp] .= thesolution.simdata.r[t,tp]

            if sum(isapprox.(thesolution.simdata.r[t,tp],thesolution2.simdata.r[t,tp])) != length(thesolution.problem.disc.fftwhelper) println("init wrong") end
            if sum(isapprox.(thesolution.simdata.F[t,tp],thesolution2.simdata.F[t,tp])) != length(thesolution.problem.disc.fftwhelper) println("init wrong") end
        end
    end

    # now i need to calc (Fr) and (F2kr2) for all t tp values in thesolution2  
    @show thesolution2.simdata.Fr[1,3]
    ichunk = 1
    for t in 1:(evolveto+2)
        for tp in 1:t
            calc_Fr_F2kr2!(t,tp,thesolution2.simdata, tmpdata2, thesolution2.problem.disc, ichunk)
            @show (tp,t)
            @show thesolution2.simdata.Fr[tp,t]
        end
    end
    println("**************************************************")
    println("*********** Preparation is complete!!! ***********")
    println("**************************************************")

    #
    # now calc Irs for CPUred and CPUred2
    #

    println("CALC IRS")
    tmone = 3
    calcIrs!( thesolution.problem.model,  thesolution.problem.pexp,  thesolution.problem.disc,  thesolution.simdata,  tmpdata, tmone) # for all tp!
    calcIrs!(thesolution2.problem.model, thesolution2.problem.pexp, thesolution2.problem.disc, thesolution2.simdata, tmpdata2, tmone) # for all tp!

    for i in 1:tmone
        @show i
        if sum( isapprox.(real.(thesolution.simdata.Irx[i]),real.(thesolution2.simdata.Irx[i]) ) ) != length(thesolution.simdata.Irx[i])
            println("Error in Irx calc")
            for j in 1:length(thesolution.simdata.Irx[i])
                println(real.(thesolution.simdata.Irx[i][j])," vs. " , real.(thesolution2.simdata.Irx[i][j] ) )
            end
            #@show real.(thesolution.simdata.Irx[i])
            #@show real.(thesolution2.simdata.Irx[i])
        end
    end

    #
    # now calc IFs for CPUred and CPUred2
    #
    println("CALC IFS")

    # copy Irx to simdata2 to make sure we have the same 
    for i in 1:tmone
        thesolution2.simdata.Irx[i] .= thesolution.simdata.Irx[i]
    end

    # check once more F and r
    for t in 1:(evolveto+2)
        for tp in 1:t
            if sum(isapprox.(thesolution.simdata.r[t,tp],thesolution2.simdata.r[t,tp])) != length(thesolution.problem.disc.fftwhelper) println("init wrong AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA") end
            if sum(isapprox.(thesolution.simdata.F[t,tp],thesolution2.simdata.F[t,tp])) != length(thesolution.problem.disc.fftwhelper) println("init wrong AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA") end
        end
    end

    tmone = 3
    # calculates at (tmone, 1) ... (tmone, tmone)
    calcIFs!( thesolution.problem.model,  thesolution.problem.pexp,  thesolution.problem.disc,  thesolution.simdata,  tmpdata, tmone) # for all tp!
    calcIFs!(thesolution2.problem.model, thesolution2.problem.pexp, thesolution2.problem.disc, thesolution2.simdata, tmpdata2, tmone) # for all tp!

    for i in 1:tmone
        @show i
        if sum( isapprox.(real.(thesolution.simdata.IFx[i]),real.(thesolution2.simdata.IFx[i]) ) ) != length(thesolution.simdata.IFx[i])
            println("Error in IFx calc")
            for j in 1:length(thesolution.simdata.IFx[i])
                println(real.(thesolution.simdata.IFx[i][j])," vs. " , real.(thesolution2.simdata.IFx[i][j] ) )
            end
            #@show real.(thesolution.simdata.Irx[i])
            #@show real.(thesolution2.simdata.Irx[i])
        end
    end

#    println("results in F2kr2")
#     for t in 1:tmone
#        for tp in 1:t
##            @show thesolution2.simdata.F2kr2[tp,t] # in k space
#            A =  thesolution2.simdata.F2kr2[tp,t]
#            copytolattice!( tmpdata.F[1], A, thesolution.problem.disc.fftwhelper)
#            mul!(tmpdata.Fx[1], tmpdata.ftplan, tmpdata.F[1]) # i need Fr in x space here
#            @show (tp,t)
#            @show tmpdata.Fx[1]./4
#        end
#    end
#
#    for tp in 1:tmone
#        A =  thesolution2.simdata.F2kr2[tp,tmone]
#        copytolattice!( tmpdata.F[1], A, thesolution.problem.disc.fftwhelper)
#        mul!(tmpdata.Fx[1], tmpdata.ftplan, tmpdata.F[1]) # i need Fr in x space here
#        @show tmpdata.Fx[1]
#    end

#    for i in 1:tmone
#        @show i
#        if sum( isapprox.(real.(thesolution.simdata.IFx[i]),real.(thesolution2.simdata.IFx[i]) ) ) != length(thesolution.simdata.Irx[i])
#            println("Error in IFx calc")
#            for j in 1:length(thesolution.simdata.IFx[i])
#                println(real.(thesolution.simdata.IFx[i][j])," vs. " , real.(thesolution2.simdata.IFx[i][j] ) )
#            end
#            #@show real.(thesolution.simdata.Irx[i])
#            #@show real.(thesolution2.simdata.Irx[i])
#        end
#    end

    return true
end

# no evolution at all
TestTwoPIScalar_CPUred2(4, 1, 1, "NinverseNLO", 0)