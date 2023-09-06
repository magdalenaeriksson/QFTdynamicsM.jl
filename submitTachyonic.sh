#!/bin/bash
#SBATCH --partition=cpu20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --job-name=TachyonicNLOcpu20
#SBATCH --output=SLURM-%x.%j.out
#SBATCH --mail-user=gerhard.ungersback@uis.no
#SBATCH --mail-type=ALL

####################################################################################################################################
# make some runs on Ninverse expanison Tachyonic
####################################################################################################################################

################################################# LO ################################################# 
#julia -t 1 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16  --sdim 3 --dt 0.1 --Nsteps 200 --NstepsinMemory 0 --Nmeas 0 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseLO --Num CPUred --Nchunks 0 --tag Tachyonic --ov
#julia -t 1 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 32  --sdim 3 --dt 0.1 --Nsteps 200 --NstepsinMemory 0 --Nmeas 0 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseLO --Num CPUred --Nchunks 0 --tag Tachyonic --ov
#julia -t 1 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 64  --sdim 3 --dt 0.1 --Nsteps 200 --NstepsinMemory 0 --Nmeas 0 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseLO --Num CPUred --Nchunks 0 --tag Tachyonic --ov
#julia -t 1 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 128 --sdim 3 --dt 0.1 --Nsteps 200 --NstepsinMemory 0 --Nmeas 0 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseLO --Num CPUred --Nchunks 0 --tag Tachyonic --ov

################################################# NLO ################################################ 
#julia -t 36 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 10000 --NstepsinMemory 550 --Nmeas 200 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUred --Nchunks 0 --tag TachyonicEXHAUST --ov
#julia -t 36 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.1 --Nsteps 10000 --NstepsinMemory 550 --Nmeas 200 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUred --Nchunks 0 --tag TachyonicEXHAUST --ov


#julia -t auto threadtest.jl

#julia -t 20 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 10000 --NstepsinMemory 550 --Nmeas 200 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUred --Nchunks 0 --tag TachyonicOPTIMcpu20 --ov
#julia -t 28 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 10000 --NstepsinMemory 550 --Nmeas 200 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUred --Nchunks 0 --tag TachyonicOPTIMcpu28 --ov
#julia -t 36 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 10000 --NstepsinMemory 550 --Nmeas 200 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUred --Nchunks 0 --tag TachyonicOPTIM72 --ov
#julia -t 72 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 10000 --NstepsinMemory 550 --Nmeas 200 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUred --Nchunks 0 --tag TachyonicOPTIM72 --ov
#julia -t 36 --check-bounds=no --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 10000 --NstepsinMemory 550 --Nmeas 200 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUred --Nchunks 0 --tag TachyonicOPTIMnobounds --ov
#julia -t 36 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 10000 --NstepsinMemory 550 --Nmeas 200 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUred --Nchunks 0 --tag TachyonicOPTIM --ov
#julia -t 36 --optimize=3 --math-mode=fast --check-bounds=no --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 10000 --NstepsinMemory 550 --Nmeas 200 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUred --Nchunks 0 --tag TachyonicOPTIMnoboundsfasto3 --ov
#julia -t 36 --check-bounds=no --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.1 --Nsteps 10000 --NstepsinMemory 550 --Nmeas 200 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUred --Nchunks 0 --tag TachyonicOPTIMnobounds --ov

#CPUred2
#julia -t 1 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 100 --NstepsinMemory 0 --Nmeas 0 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUred2 --Nchunks 0 --tag TachyonicCPUred2 --ov

julia -t 20 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 10000 --NstepsinMemory 550 --Nmeas 200 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUred2 --Nchunks 0 --tag TachyonicCPUred2cpu20 --ov

#julia -t 36 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 16 --sdim 3 --dt 0.1 --Nsteps 10000 --NstepsinMemory 550 --Nmeas 200 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUred2 --Nchunks 0 --tag TachyonicCPUred2 --ov
#julia -t 36 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.1 --Nsteps 10000 --NstepsinMemory 550 --Nmeas 200 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUred2 --Nchunks 0 --tag TachyonicCPUred2 --ov
#julia -t 36 --project=. G_SimulationArgs.jl --Mass 0.7 --Nx 64 --sdim 3 --dt 0.1 --Nsteps 10000 --NstepsinMemory 550 --Nmeas 200 --Mod phi4tachyonic --Lambda 1 --ONgroup 4 --init Pnr --T 1.5 --n 0 --xi 0.4 --eta 5.0 --sig 1.0 --Reno RMass --Pexp NinverseNLO --Num CPUred2 --Nchunks 0 --tag TachyonicCPUred2 --ov

################################################# CS ################################################# 
#julia -t 1 --project=. G_CSSimulationArgs.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.1 --Nsteps 200 --Nmeas 0 --Mod CSphi4tachyonicUNS --Lambda 1 --ONgroup 4 --init Pnr --n 0 --Reno RNone --Num CPU --Runs 10 --tag TachyonicRuns10Thread1 --ov
#julia -t 10 --project=. G_CSSimulationArgs.jl --Mass 0.7 --Nx 32 --sdim 3 --dt 0.1 --Nsteps 200 --Nmeas 0 --Mod CSphi4tachyonicUNS --Lambda 1 --ONgroup 4 --init Pnr --n 0 --Reno RNone --Num CPU --Runs 10 --tag TachyonicRuns10Thread10 --ov
