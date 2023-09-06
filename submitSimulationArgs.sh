#!/bin/bash
#SBATCH --partition=cpu28
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=1
#SBATCH --job-name=CPUredtest
#SBATCH --output=SLURM-%x.%j.out
#SBATCH --mail-user=gerhard.ungersback@uis.no
#SBATCH --mail-type=ALL

####################################################################################################################################
# 1/N Expansion
####################################################################################################################################
# NNLO paper - Figure 7
#julia -t 10 --project=. G_Simulation.jl --Mass 1 --Nx 1 --sdim 1 --dt 0.01 --Nsteps 600 --NstepsinMemory 0 --Nmeas 0 --Mod phi4 --Lambda 1 --ONgroup 8 --init Gauss --T 5 --n 2 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RNone --Pexp NinverseNLO --Num CPUfull --Nchunks 0 --tag AndersNNLOpaper --ov
#julia -t 10 --project=. G_Simulation.jl --Mass 1 --Nx 1 --sdim 1 --dt 0.01 --Nsteps 600 --NstepsinMemory 0 --Nmeas 0 --Mod phi4 --Lambda 1 --ONgroup 8 --init Gauss --T 5 --n 2 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RNone --Pexp LambdaLO    --Num CPUfull --Nchunks 0 --tag AndersNNLOpaper --ov
#julia -t 10 --project=. G_Simulation.jl --Mass 1 --Nx 1 --sdim 1 --dt 0.01 --Nsteps 600 --NstepsinMemory 0 --Nmeas 0 --Mod phi4 --Lambda 1 --ONgroup 8 --init Gauss --T 5 --n 2 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RNone --Pexp NinverseLO  --Num CPUfull --Nchunks 0 --tag AndersNNLOpaper --ov

# increase Nx, use Thermal
#julia -t 10 --project=. G_Simulation.jl --Mass 1 --Nx 512 --sdim 1 --dt 0.01 --Nsteps 600 --NstepsinMemory 0 --Nmeas 0 --Mod phi4 --Lambda 6 --ONgroup 1 --init Thermal --T 1 --n 2 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RNone --Pexp NinverseNLO --Num CPUfull --Nchunks 0 --tag testhigherNx --ov
#julia -t 10 --project=. G_Simulation.jl --Mass 1 --Nx 512 --sdim 1 --dt 0.01 --Nsteps 600 --NstepsinMemory 0 --Nmeas 0 --Mod phi4 --Lambda 6 --ONgroup 1 --init Thermal --T 1 --n 2 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RNone --Pexp LambdaLO    --Num CPUfull --Nchunks 0 --tag testhigherNx --ov
#julia -t 10 --project=. G_Simulation.jl --Mass 1 --Nx 512 --sdim 1 --dt 0.01 --Nsteps 600 --NstepsinMemory 0 --Nmeas 0 --Mod phi4 --Lambda 6 --ONgroup 1 --init Thermal --T 1 --n 2 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RNone --Pexp NinverseLO  --Num CPUfull --Nchunks 0 --tag testhigherNx --ov

#Tachionic
#julia -t 20 --project=. G_Simulation.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.1 --Nsteps 250 --NstepsinMemory 0 --Nmeas 0 --Mod phi4 --Lambda 1 --ONgroup 4 --init Pnr --T 5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RNone --Pexp NinverseLO --Num CPUfull --Nchunks 0 --tag Tachyonic --ov
#julia -t 20 --project=. G_Simulation.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.05 --Nsteps 500 --NstepsinMemory 0 --Nmeas 0 --Mod phi4 --Lambda 1 --ONgroup 4 --init Pnr --T 5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseLO --Num CPUfull --Nchunks 0 --tag Tachyonic --ov

####################################################################################################################################
# Lambda Expansion, Equilibration in 3+1 -> reproduction of Anders paper
####################################################################################################################################
# Test #
#julia -t 1 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 100 --NstepsinMemory 0 --Nmeas 0 --Mod phi4 --Lambda 6 --ONgroup 1 --init Thermal --T 2.7 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag Equilibration --ov

#julia -t 10 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 1400 --NstepsinMemory 400 --Nmeas 200 --Mod phi4 --Lambda 6 --ONgroup 1 --init Thermal --T 2.7 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag EquilibrationTincor --ov
#julia -t 10 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 1400 --NstepsinMemory 400 --Nmeas 200 --Mod phi4 --Lambda 6 --ONgroup 1 --init Thermal --T 2.2 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag EquilibrationTincor --ov
#julia -t 10 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 1400 --NstepsinMemory 400 --Nmeas 200 --Mod phi4 --Lambda 6 --ONgroup 1 --init Thermal --T 1.8 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag EquilibrationTincor --ov

#julia -t 10 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 1400 --NstepsinMemory 800 --Nmeas 200 --Mod phi4 --Lambda 6 --ONgroup 1 --init Thermal --T 1.4 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag EquilibrationTincor --ov
#julia -t 10 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 1400 --NstepsinMemory 800 --Nmeas 200 --Mod phi4 --Lambda 6 --ONgroup 1 --init Thermal --T 1.1 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag EquilibrationTincor --ov

#julia -t 10 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 1400 --NstepsinMemory 800 --Nmeas 200 --Mod phi4 --Lambda 6 --ONgroup 1 --init Thermal --T 0.7 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag EquilibrationTincor --ov
#julia -t 10 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 1400 --NstepsinMemory 800 --Nmeas 200 --Mod phi4 --Lambda 6 --ONgroup 1 --init Thermal --T 0.3 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag EquilibrationTincor --ov

####################################################################################################################################
# make some runs on lambda expanison Equilibration
####################################################################################################################################
#julia -t 2 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 10 --NstepsinMemory 0 --Nmeas 0 --Mod phi4 --Lambda 6 --ONgroup 1 --init Thermal --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag EquSize --ov

#julia -t 36 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 1400 --NstepsinMemory 400 --Nmeas 200 --Mod phi4 --Lambda 6 --ONgroup 1 --init Thermal --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag EquSize --ov

#julia -t 36 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.1 --Nsteps 1400 --NstepsinMemory 400 --Nmeas 200 --Mod phi4 --Lambda 6 --ONgroup 1 --init Thermal --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag EquSize --ov

#julia -t 28 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 56 --sdim 3 --dt 0.1 --Nsteps 1400 --NstepsinMemory 400 --Nmeas 200 --Mod phi4 --Lambda 6 --ONgroup 1 --init Thermal --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag EquSize28 --ov

####################################################################################################################################
# make some runs on Ninverse expanison Tachyonic
####################################################################################################################################
# LO
#julia -t 2 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.1 --Nsteps 200 --NstepsinMemory 0 --Nmeas 0 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseLO --Num CPUfull --Nchunks 0 --tag Tachyonic --ov
#julia -t 2 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.1 --Nsteps 200 --NstepsinMemory 0 --Nmeas 0 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaLO --Num CPUfull --Nchunks 0 --tag Tachyonic --ov

# NLO
#julia -t 10 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.1 --Nsteps 200 --NstepsinMemory 0 --Nmeas 0 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUfull --Nchunks 0 --tag Tachyonic --ov
#julia -t 14 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.1 --Nsteps 200 --NstepsinMemory 0 --Nmeas 0 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUfull --Nchunks 0 --tag Tachyonicfast --ov


#julia -t 28 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.01 --Nsteps 2000 --NstepsinMemory 400 --Nmeas 200 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUfull --Nchunks 0 --tag Tachyonic --ov
#julia -t 28 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.01 --Nsteps 2000 --NstepsinMemory 800 --Nmeas 200 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUfull --Nchunks 0 --tag Tachyonic --ov

####################################################################################################################################
# compare CPUfull and CPUred
####################################################################################################################################
#julia -t 28 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.1 --Nsteps 1400 --NstepsinMemory 400 --Nmeas 200 --Mod phi4 --Lambda 6 --ONgroup 1 --init Thermal --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUfull --Nchunks 0 --tag tryCPUfull2 --ov
#julia -t 28 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.1 --Nsteps 1400 --NstepsinMemory 400 --Nmeas 200 --Mod phi4 --Lambda 6 --ONgroup 1 --init Thermal --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp LambdaNLO --Num CPUred --Nchunks 0 --tag tryCPUred2 --ov

#julia -t 28 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.1 --Nsteps 200 --NstepsinMemory 0 --Nmeas 0 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUfull --Nchunks 0 --tag tryCPUfull2 --ov
#julia -t 28 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.1 --Nsteps 200 --NstepsinMemory 0 --Nmeas 0 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUred --Nchunks 0 --tag tryCPUred2 --ov
