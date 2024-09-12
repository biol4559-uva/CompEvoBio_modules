# **Short Read Mapping & Pipelines**

## Overview
In order to turn your short-read data into something useable, you need to map it to a reference genome. Once the data is mapped back to the reference genome you can identify SNPs because of mismatches between the short-reads and the reference genome.

<p align="center">
  <img src="/Module_3/images/tview.png" width="1000"/>
</p>

The process of doing this mapping takes many steps. This includes trimming adapter sequences and low quality sequences from reads, mapping to a reference genome, removing PCR duplicates, filtering, and quality control. Each of these steps is performed by a different tool, and the steps need to be arranged in a pipeline. To facilitate the mapping process, we are using a containerized (Docker) version of this pipeline. A containerized pipeline is a framework to preserve pipeline steps and software versions for complicated pipelines. The Docker container is basically the operating system and software that get loaded on the CPUs that you request for your SLURM job.
<p align="center">
  <img src="/Module_3/images/mapping_pipeline.png" width="1000"/>
</p>

The Docker pipeline that we will be using is part of the DEST project. This is a project that aims to collect allele frequency data for fly samples colleced around the world. The benefit of using the Docker pipeline that we used for DEST in this class is that we can combine our data with it in order to situate your samples in a broader evolutionary context.
<p align="center">
  <img src="/Module_3/images/Figure1.png" width="1000"/>
</p>


Goals:
1. learn about mapping pipelines

Objective:
1. Map your data using the DEST docker. Send the output of your samples to our shared project folder

## Instructions
1. Using the [csv table](/Module_3/sras.txt)that you and your group constructed, divide up the samples roughly evenly. Copy your lines to a new file and save in your scratch folder.

2. Make a new directory in your scratch directory called `/scratch/COMPUTEID/mapping_output`

3. Using Globus, copy `/standard/BerglandTeach/dest_freeze2.6.1_latest.sif` to your scratch directory

2. Copy the `mapping_script.sh` to a new file using Rstudio. Modify that file to have the proper paths and decide if you need PE or SE. Copy the lines below to a dummy text file, and modify `PATH_to_your_small_file` to point it to the small csv file you made in the previous step. Also modify the `PROJECTID` to include the proper SRA project ID. Record the job number.

```
sbatch --array=1-$( cat PATH_to_your_small_file  | wc -l  ) \
PATH_to_your_version_of_mapping_script.sh \
/scratch/COMPUTEID/dest_v2.6.1_latest.sif \
/scratch/COMPUTEID/PROJECTID \
/standard/BerglandTeach/mapping_output \
PATH_to_your_small_file
```

my example
```
sbatch --array=1-$( cat /scratch/aob2x/small_file.csv  | wc -l  ) \
~/CompEvoBio_modules/Module_3/old/mapping_script.aob2x.sh \
/scratch/aob2x/dest_freeze2.6.1_latest.sif \
/scratch/aob2x/compBio/fastq/SRP002024 \
/standard/BerglandTeach/mapping_output \
/scratch/aob2x/small_file.csv
```

3. Confirm your job is successfully running by running the command `sacct -j JOBID`

4. To complete this assignment, upload your JOBID to Cavnas.
