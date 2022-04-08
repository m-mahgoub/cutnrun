#!/bin/bash
module load nextflow/21.10.5 
module load singularity
nextflow run test.nf -resume -profile slurm,singularity