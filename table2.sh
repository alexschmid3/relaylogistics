#!/bin/bash
#SBATCH -a 5,21,29,42,66,58
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=sched_mit_sloan_batch
#SBATCH --time=2-00:00
#SBATCH -o /home/aschmid/relaylogistics/outerr/table2_conv_\%a.out
#SBATCH -e /home/aschmid/relaylogistics/outerr/table2_conv_\%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aschmid@mit.edu

module load julia/1.5.2
module load sloan/gurobi/10.0.2

julia run_static.jl $SLURM_ARRAY_TASK_ID