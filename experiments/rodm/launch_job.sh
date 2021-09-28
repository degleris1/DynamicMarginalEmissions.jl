#!/bin/bash
#
#SBATCH --job-name=carbon_rodm
#
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G

SCIPT_DIR=/home/users/degleris/CarbonNetworks.jl/experiments/rodm/

srun hostname
module load julia

srun julia ${SCIPT_DIR}${1}