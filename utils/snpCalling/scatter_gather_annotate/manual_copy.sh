#!/usr/bin/env bash

#SBATCH -J manual_gather # A single job name for the array
#SBATCH --ntasks-per-node=8 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 4:00:00 ### 1 hours
#SBATCH --mem 64G
#SBATCH -o /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_gather.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_gather.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard


# ijob -A berglandlab -c10 -p standard --mem=50G
# sbatch --array=1-177 ~/CompEvoBio_modules/utils/snpCalling/scatter_gather_annotate/manual_scatter.sh
# sacct -j 4242244

# SLURM_ARRAY_TASK_ID=1
cd /standard/BerglandTeach/mapping_output/
samp=$( ls -d * | tr '\t' '\n' | sed "${SLURM_ARRAY_TASK_ID}q;d" )
echo $samp

if [[ ! -d "/project/berglandlab/DEST/dest_mapped/ExpEvo/$samp" ]]; then
   mkdir /project/berglandlab/DEST/dest_mapped/ExpEvo/$samp
fi

ls /standard/BerglandTeach/mapping_output/${samp}/*.masked.sync.gz*
cp /standard/BerglandTeach/mapping_output/${samp}/*.masked.sync.gz* /project/berglandlab/DEST/dest_mapped/ExpEvo/$samp
