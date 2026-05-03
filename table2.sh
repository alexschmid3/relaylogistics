#!/bin/bash
#SBATCH -a 541,542,543,544,561,562,563,564,961,963,965,967,969,971,973,975,977,979,981,983,985,987,989,991,1473,1474,1475,1476,1477,1478,1479,1480
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=sched_mit_sloan_batch
#SBATCH --time=2-00:00
#SBATCH -o /home/aschmid/relaylogistics/outerr/table2_msomrev_\%a.out
#SBATCH -e /home/aschmid/relaylogistics/outerr/table2_msomrev_\%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aschmid@mit.edu

module load julia/1.5.2
module load sloan/gurobi/10.0.2

julia run_static.jl $SLURM_ARRAY_TASK_ID