#!/bin/bash
#SBATCH -a 225-230
#SBATCH --partition=xeon-p8
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2-00:00
#SBATCH -o /home/gridsan/aschmid/relaylogistics/outerr/sense_%a.out
#SBATCH -e /home/gridsan/aschmid/relaylogistics/outerr/sense_%a.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=aschmid@mit.edu

module load julia/1.9.2
module load gurobi/gurobi-1102

julia run_static.jl $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT