#!/bin/bash
#
#SBATCH --job-name=carbon_rodm
#
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G

SCIPT_DIR=/home/users/degleris/CarbonNetworks.jl/experiments/rodm/

srun hostname
module load julia

srun julia -t 8 ${SCIPT_DIR}${1}