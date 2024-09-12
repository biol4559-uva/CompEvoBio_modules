#!/usr/bin/env bash

#SBATCH -J runFASTQC # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:00:15 ### 15 seconds
#SBATCH --mem 2G
#SBATCH -o /scratch/COMPUTE_ID/logs/demo_1.%A_%a.out # Standard output
#SBATCH -e /scratch/COMPUTE_ID/logs/demo_1.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol4020-aob2x

### run as: sbatch PATH_TO_THIS_FILE
### sacct -j XXXXXXXXX
### cat /scratch/COMPUTE_ID/logs/demo_1.*.err

echo "Hello world"
