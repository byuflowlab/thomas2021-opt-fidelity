#!/bin/bash

#SBATCH --time=02:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=2048M   # memory per CPU core
#SBATCH -J "LES low-ti"   # job name
#SBATCH --array=0-16     # job array of size 200

echo ${SLURM_ARRAY_TASK_ID}

# load julia module
module load julia

# set up
dirbins=(5 10 15 20 30 40 50 70 90 110 140 170 200 240 280 320 360)

# run julia
julia run_opt.jl --firstrun 1 --nruns 100 --case "low-ti" --out-dir "./test/" --dir-bins ${dirbins[${SLURM_ARRAY_TASK_ID}]}