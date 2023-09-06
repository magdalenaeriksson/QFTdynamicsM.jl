#!/bin/bash
# usage: sbatch -p gpuA100 --gres=gpu:1 BMsubmit_TwoPIScalar.sh

julia --project=. BM_TwoPIScalar.jl
