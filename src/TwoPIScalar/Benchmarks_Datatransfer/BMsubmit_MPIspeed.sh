#!/bin/bash -l
#SBATCH --partition=cpu20
#SBATCH --nodes=2
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00 # time=HH:MM:SS
#SBATCH -J "thejob"   # job name

# usage: sbatch BMsubmit_ClusterManager.sh
mpiexecjl -np 40 julia --project=/home/stud/gerhard/bhome/gerhard/coding/juliatesting /home/stud/gerhard/bhome/gerhard/coding/juliatesting/BM_MPIspeed.jl
