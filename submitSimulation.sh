#!/bin/bash
#SBATCH --partition=cpu28
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --job-name=NNLO

#julia -t 1 --project=. G_Simulation.jl
#julia -t 10 --project=. src/TwoPIScalar/Simulationfiles/Anders_NNLOpaper.jl
#julia -t 10 --project=. src/TwoPIScalar/Simulationfiles/Anders_NNLOpaperbiggerNx.jl

julia -t 10 --project=. src/TwoPIScalar/Simulationfiles/Anders_Tachyonic.jl
