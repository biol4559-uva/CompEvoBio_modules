#!/usr/bin/env bash
#
#SBATCH -J bamqc # A single job name for the array
#SBATCH -c 5 ### 5 cores
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00
#SBATCH --mem 10G
#SBATCH -o /project/biol4559-aob2x/mapping_scripts/COMPUTEID/logs/bamqc.%A_%a.out # Standard output
#SBATCH -e /project/biol4559-aob2x/mapping_scripts/COMPUTEID/logs/bamqc.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol4559-aob2x

### run as: sbatch PATH_TO_THIS_FILE

### modules
  module load qualimap

### which bioproject
  bioproj=SRP002024
  computeId=aob2x
  # SLURM_ARRAY_TASK_ID=1

### make output directory

  if [ ! -d /project/biol4559-aob2x/mapping_scripts/${computeId}/bamqc_output/ ]; then
    mkdir   /project/biol4559-aob2x/mapping_scripts/${computeId}/bamqc_output/
  fi

### run bamqc
  ls /project/biol4559-aob2x/mapping_output/*${bioproj}*/*.original.bam | awk -F/ '{print $5"\t"$0}' > /project/biol4559-aob2x/mapping_scripts/${computeId}/bamqc_output/bamSamples.txt

  qualimap multi-bamqc -r -d /project/biol4559-aob2x/mapping_scripts/${computeId}/bamqc_output/bamSamples.txt -outdir /project/biol4559-aob2x/mapping_scripts/${computeId}/bamqc_output
