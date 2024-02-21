#!/bin/bash
#SBATCH -a 9-10
#SBATCH --partition=xeon-p8
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2-00:00
#SBATCH -o /home/gridsan/aschmid/relaylogistics/outerr/table3_%a.out
#SBATCH -e /home/gridsan/aschmid/relaylogistics/outerr/table3_%a.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=aschmid@mit.edu

module load julia/1.9.2
module load gurobi/gurobi-1000

julia run_static.jl $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT