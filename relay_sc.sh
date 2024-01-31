#!/bin/bash
#SBATCH -a 1557-1568
#SBATCH --partition=xeon-p8
#SBATCH --mem=64G
#SBATCH --time=2-00:00
#SBATCH -o /home/gridsan/aschmid/relaylogistics/outerr/resub_%a.out
#SBATCH -e /home/gridsan/aschmid/relaylogistics/outerr/resub_%a.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=aschmid@mit.edu

module load julia/1.9.2
module load gurobi/gurobi-1000

julia run_static.jl $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT