#!/bin/bash
#SBATCH -a 2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=sched_mit_sloan_batch
#SBATCH --time=0-09:00
#SBATCH -o /home/aschmid/relaylogistics/outerr/corridors_\%a.out
#SBATCH -e /home/aschmid/relaylogistics/outerr/corridors_\%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aschmid@mit.edu

module load julia/1.5.2
module load sloan/gurobi/10.0.2

julia run_online.jl $SLURM_ARRAY_TASK_ID