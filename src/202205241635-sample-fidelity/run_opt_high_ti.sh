#!/bin/bash

#SBATCH --time=48:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH -J "LES high-ti dir fidelity"   # job name
#SBATCH --array=0-7     # job array of size 200

echo ${SLURM_ARRAY_TASK_ID}

# load julia module
module load julia

# set up
dirbins=(1 5 10 20 40 60 80 100)

# run julia
julia run_opt.jl --firstrun 1 --nruns 100 --case "high-ti" --out-dir "./high-ti/" --dir-bins ${dirbins[${SLURM_ARRAY_TASK_ID}]}