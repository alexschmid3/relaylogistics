#!/bin/bash
#SBATCH -a 49-52
#SBATCH --partition=xeon-p8
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2-00:00
#SBATCH -o /home/gridsan/aschmid/relaylogistics/outerr/timedelta_%a.out
#SBATCH -e /home/gridsan/aschmid/relaylogistics/outerr/timedelta_%a.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=aschmid@mit.edu

module load julia/1.9.2
module load gurobi/gurobi-1102

julia run_online.jl $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT