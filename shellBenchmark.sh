#!/bin/bash

part=cpu28
threads=28

# LambdaNLO
#for i in 1 2 3 4 5 6 7 8 9 10
#    do
#    for Nx in 16 32 42 54 64 112 128 144 162 180
#    #for Nx in 16 32
#        do
#            for NstepsinMemory in 400 600
#            do
#                sbatch -J BM -p $part -n $threads -c 1 --nodes 1 --wrap="julia -t $threads --project=. G_SingleBenchmarkTwoPIScalar.jl --Mass 0.7 --Nx $Nx --sdim 3 --dt 0.01 --Nsteps 1000 --NstepsinMemory $NstepsinMemory --Nmeas 200 --Mod phi4 --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag BMv2 --ov"
#                sbatch -J BM -p $part -n $threads -c 1 --nodes 1 --wrap="julia -t $threads --project=. G_SingleBenchmarkTwoPIScalar.jl --Mass 0.6 --Nx $Nx --sdim 3 --dt 0.01 --Nsteps 1000 --NstepsinMemory $NstepsinMemory --Nmeas 200 --Mod phi4 --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUred --Nchunks 0 --tag BMv2 --ov"
#            done
#    done
#done

#NinverseNLO
for i in 1
    do
    #for Nx in 16 32 42 54 64 112 128 144 162 180
    for Nx in 16 #32
        do
            #for NstepsinMemory in 400 600
            for NstepsinMemory in 400
            do
                sbatch -J BMNinv -p $part -n $threads -c 1 --nodes 1 --wrap="julia -t $threads --project=. G_SingleBenchmarkTwoPIScalar.jl --Mass 0.7 --Nx $Nx --sdim 3 --dt 0.01 --Nsteps 1000 --NstepsinMemory $NstepsinMemory --Nmeas 200 --Mod phi4 --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUfull --Nchunks 0 --tag BMv2 --ov"
                sbatch -J BMNinv -p $part -n $threads -c 1 --nodes 1 --wrap="julia -t $threads --project=. G_SingleBenchmarkTwoPIScalar.jl --Mass 0.6 --Nx $Nx --sdim 3 --dt 0.01 --Nsteps 1000 --NstepsinMemory $NstepsinMemory --Nmeas 200 --Mod phi4 --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUred --Nchunks 0 --tag BMv2 --ov"
            done
    done
done
