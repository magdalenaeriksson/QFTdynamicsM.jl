#!/bin/bash
# usage: sbatch -p gpuA100 --gres=gpu:1 BMsubmit_GPUspeed.sh

julia --project=~/gorina11 BM_GPUspeed.jl
