#!/bin/bash

#SBATCH --time=02:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=2048M   # memory per CPU core
#SBATCH -J "LES low-ti"   # job name
#SBATCH --array=1-10     # job array of size 200

echo ${SLURM_ARRAY_TASK_ID}

# load julia module
module load julia

# set up
nruns=40
((firstrun=$nruns*${SLURM_ARRAY_TASK_ID}+1))
echo $firstrun

# run julia
julia run_opt.jl --firstrun $firstrun --nruns $nruns --case "low-ti" --out-dir "./test/"