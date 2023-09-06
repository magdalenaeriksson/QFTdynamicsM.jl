#!/bin/bash
#SBATCH --partition=cpu20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#julia -t 1 --project=. G_SinglePlot.jl --file phi4LambdaNLO_L600_ON1_RNone_Nx16_sdim2_Nsteps1000_dt100_M50_NstepsinMemory1002_testing
#julia -t 1 --project=. G_SinglePlot.jl --file phi4LambdaNLO_L600_ON1_RNone_Nx16_sdim2_Nsteps2000_dt50_M50_NstepsinMemory800_testing
#julia -t 1 --project=. G_SinglePlot.jl --file phi4LambdaNLO_L600_ON1_RMass_Nx16_sdim3_Nsteps1400_dt100_M70_NstepsinMemory400_Anders
#julia -t 1 --project=. G_SinglePlot.jl --file phi4LambdaNLO_L600_ON1_RNone_Nx16_sdim3_Nsteps1400_dt100_M70_NstepsinMemory400_Anders
#julia -t 1 --project=. G_SinglePlot.jl --file phi4LambdaNLO_L1800_ON1_RNone_Nx32_sdim2_Nsteps10000_dt100_M100_NstepsinMemory400_twoplusone

julia -t 1 --project=. G_SinglePlot.jl --file phi4LambdaLO_L0_ON1_RNone_Nx16_sdim2_Nsteps1000_dt100_M50_NstepsinMemory400_testing
