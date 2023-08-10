#!/usr/bin/env bash
#
#SBATCH -J makeSingularity # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 4:00:00 ### 6 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/aob2x/compBio/logs/makeSingularity.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/compBio/logs/makeSingularity.%A_%a.err # Standard error
#SBATCH -p instructional
#SBATCH --account biol4559-aob2x


module load singularity
cd /project/biol4559-aob2x/singularity
singularity pull docker://jcbn/dest_v2.5
