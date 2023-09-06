#!/bin/bash -l
#SBATCH --partition=cpu20
#SBATCH --nodes=2
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00 # time=HH:MM:SS
#SBATCH -J "thejob"   # job name

# usage: sbatch BMsubmit_ClusterManager.sh
#julia --project=~/gorina11 BM_ClusterManager.jl
#julia --project=. /home/stud/gerhard/bhome/gerhard/coding/juliatesting/BM_ClusterManager.jl
julia /home/stud/gerhard/bhome/gerhard/coding/juliatesting/BM_ClusterManager.jl
