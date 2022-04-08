#!/bin/bash

# sbatch --cpus-per-task 32 --mem=64g --gres=lscratch:500 --time=6:00:00 run_sbatch.sh

module load nextflow/21.10.5 
module load singularity
nextflow run test.nf -resume -profile slurm,singularity
