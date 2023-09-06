#!/bin/bash
#SBATCH --partition=cpu28
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=1
#SBATCH --job-name=BMNLgcfree
#SBATCH --output=SLURM-%x.%j.out
#SBATCH --mem=245G

#
## BM Script -> not working
#
#julia -t 20 --project=. G_BenchmarkTwoPIScalar.jl --Pexp LambdaNLO
#julia -t 20 --project=. G_BenchmarkTwoPIScalar.jl --Pexp NinverseNLO

#julia -t 28 --project=. G_BenchmarkTwoPIScalar.jl --Pexp LambdaNLO
#julia -t 28 --project=. G_BenchmarkTwoPIScalar.jl --Pexp NinverseNLO

#julia -t 36 --project=. G_BenchmarkTwoPIScalar.jl --Pexp LambdaNLO
#julia -t 36 --project=. G_BenchmarkTwoPIScalar.jl --Pexp NinverseNLO

#
# BM individual ones
#
#julia -t 28 --project=. G_SingleBenchmarkTwoPIScalar.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.01 --Nsteps 1000 --NstepsinMemory 600 --Nmeas 200 --Mod phi4 --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag benchi --ov

#for Nx in 16 32
#    do
#        for NstepsinMemory in 400 600
#        do
#            echo "Working on $Nx / $NstepsinMemory ..."
#            julia -t 28 --project=. G_SingleBenchmarkTwoPIScalar.jl --Mass 0.7 --Nx $Nx --sdim 3 --dt 0.01 --Nsteps 1000 --NstepsinMemory $NstepsinMemory --Nmeas 200 --Mod phi4 --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag BMv1 --ov
#        done
#done

#
# old, was about looking for bugs
#

# try the one at which it fucks up
#julia -t 1 --project=. G_SingleBenchmarkTwoPIScalar.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.01 --Nsteps 1000 --NstepsinMemory 600 --Nmeas 200 --Mod phi4 --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag benchi --ov

# try one of the biggest ones of CPUred
#julia -t 28 --project=. G_SingleBenchmarkTwoPIScalar.jl --Mass 0.7 --Nx 128 --sdim 3 --dt 0.01 --Nsteps 1000 --NstepsinMemory 600 --Nmeas 200 --Mod phi4 --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RNone --Pexp LambdaNLO --Num CPUred --Nchunks 0 --tag benchi --ov
#julia -t 28 --project=. G_SingleBenchmarkTwoPIScalar.jl --Mass 0.7 --Nx 108 --sdim 3 --dt 0.01 --Nsteps 1000 --NstepsinMemory 600 --Nmeas 200 --Mod phi4 --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RNone --Pexp LambdaNLO --Num CPUred --Nchunks 0 --tag benchi --ov

# try one of the biggest ones of CPUfull
#julia -t 28 --project=. G_SingleBenchmarkTwoPIScalar.jl --Mass 0.7 --Nx 42 --sdim 3 --dt 0.01 --Nsteps 1000 --NstepsinMemory 600 --Nmeas 200 --Mod phi4 --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RNone --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag benchi --ov
#julia -t 28 --project=. G_SingleBenchmarkTwoPIScalar.jl --Mass 0.7 --Nx 36 --sdim 3 --dt 0.01 --Nsteps 1000 --NstepsinMemory 600 --Nmeas 200 --Mod phi4 --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RNone --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag benchi --ov

## print memory consumption
#/usr/bin/time -f "%M" julia -t 1 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.01 --Nsteps 100 --NstepsinMemory 0 --Nmeas 0 --Mod phi4 --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag benchi --ov
#/usr/bin/time -f "%M" julia -t 28 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 36 --sdim 3 --dt 0.01 --Nsteps 10 --NstepsinMemory 600 --Nmeas 0 --Mod phi4 --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag benchi --ov
#
