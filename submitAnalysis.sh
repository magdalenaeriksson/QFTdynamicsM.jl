#!/bin/bash
#SBATCH --partition=cpu28
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=Analysis

julia -t 1 --project=. G_Analysis.jl
