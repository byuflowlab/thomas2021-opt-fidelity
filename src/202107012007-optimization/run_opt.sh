#!/bin/bash

#SBATCH --time=00:00:05   # walltime
#SBATCH --ntasks=3   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=2048M   # memory per CPU core
#SBATCH -J "LES low-ti"   # job name
#SBATCH --qos=test

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load julia

julia run_opt.jl
