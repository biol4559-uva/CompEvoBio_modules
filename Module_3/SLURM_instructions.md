# **SLURM**

## Overview
Up to this point, we have been using Rivanna in an interactive way. That is great when we are analyzing data and learning to code. But, sometimes, we need to let Rivanna work for hours and we want it to run in the background. To submit jobs to Rivanna, we use the job handler SLURM.

Goals:
1. learn about the structure of a SLURM file
2. learn about modules on Rivanna
3. learn about batch processing

Objective:
1. Run a quality control tool on FASTQ files for your paper using a slurm script

## Instructions -
### Basics of a SLURM script
1. Open a new Rstudio job. For this exercise, we are only going to RStudio to work as a script editor and to use GitHub. Once you have Rstudio open, load your GitHub Repo project, and create a new script. Copy the contents of the contents of `SLURM_demo_1.sh` into the script. Save the file in your GitHub folder in a new folder called "SLURM_assignment". Make sure to save the file name as "SLURM_demo_1.sh"

2. Open a terminal window and make a new "log" directory in your scratch directory. To do this, nagivate to your scratch directory (`/scratch/COMPUTE_ID`). Then, make a new directory (`mkdir logs`)

3. In RStudio, edit the script text to point to your scratch directory (lines 8 & 9).

4. Modify line 13 to point to the SLURM_demo_1.sh script. Save everything.

5. Copy the bit on line 13 after the "run as: " bit. It should be something like `sbatch ~/biol4559/aob2x/SLURM_assignment/SLURM_demo_1.sh`. Go back to your terminal, and paste that in the shell.

6. Copy the job job number and paste it back on line 14, where the XXXXXXX is. Now, copy the `sacct -j XXXXXXXXX` line and paste it back in the terminal. What is the status of your job?

7. Look at the contents of the output files. What do they say?

### Modules & fastq quality control
1. Before we begin analyzing FASTQ data, we want to understand some basics of its quality. For instance, we can ask how base quality scores change along the sequence read. Or, we can ask what the proportions of A, C, T, and G are in the sample. These sorts of assessments can help us identify if the sequencing run "failed", or if it is contaminated with somethign unexpected. (different species have different CG proportions).

2. We will perform this task using two tools. The first is called `fastqc` and the second is `multiqc`. Fastqc takes a single fastq file and runs some analyses. Multiqc takes the output of a bunch of fastqc files and produces a summary report.

3. To use these tools, we need to load them using the `module` command. On Rivanna, most programs need to be loaded first, much like you need to load libraries in R. We also need to make an output folder. Copy these lines to a blank Rscript file; modify the COMPUTEID value for your computeID:
```
module load fastqc  ### this load the module

proj=SRP002024 ### defines the (Short Read Archive) project identifier. (the paper)
srr=SRR036932 ### defines the SRA sample identifier (the unique sample)

if [ ! -d /scratch/COMPUTEID/fastq_QC_Out/ ]; then ### these three lines will make an output folder
  mkdir /scratch/COMPUTEID/fastq_QC_Out/
fi

fastqc -o /scratch/COMPUTEID/fastq_QC_Out --noextract /project/biol4559-aob2x/data/fastq/${proj}/${srr}.fastq.gz

cp -R /scratch/COMPUTEID/fastq_QC_Out ~/fastq_QC_Out

```

4. Down-load the contents of `~/fastq_QC_Out` to your computer using the OpenOnDemand file transfer tool and open the html file. What do you see?

### Basics of batch processing
1. We want to be able to run fastqc on all of the samples associated with our paper. To do this, we will utilize batch processing which is a bit like a foreach loop in R. We have a list of samples we want to process, and each sample is processed independently. For SLURM, this is called an array.

2. To run a slurm script as an array, we will type:

```
sbatch --array=1-5 SLURM_demo_2.sh
```

In this way, our slurm script is a bit like a function that only takes one parameter, a number between X and Y in this case, between 1-5. Within the slurm script, that parameter gets caputred by the `${SLURM_ARRAY_TASK_ID}` variable. To run this script, you will need to copy it into R, and modify the paths, and save on Rivanna. When you run it, what files do you look at to ensure it works?

3. Let's say that we have a list of sample names, and we want to run `fastqc` on each of them. In order to use the `${SLURM_ARRAY_TASK_ID}` variable, we need to link it to a premade table of possible values. Then, we can look up the ith line and use that information as our parameter values that we actually need. Copy `SLURM_demo_3.sh` and modify as appropriate.

### Run FASTQC on all of the samples from your paper
1. Copy and modify `SLURM_demo_4.sh` so that it fits with your path names, as well as your paper `proj` variable. Run it, and check that it is running.

### Run MULTIQC.
1. Combine your FASTQC output files into a single comparative file using MULTIQC:
```
module load multiqc
multiqc -o ~/multiQC_output/ /scratch/COMPUTEID/fastq_QC_Out/
```

2. Download the two contents of that foler and open the html file.

### To complete this assignment
6. To complete this assignment take screen shots of the General Stats table, and the Sequence Quality Histograms figures and upload them to Canvas. Write 3-4 sentences to describe these figures. Watch the insstructional video. Are there any samples that you are concerned about? Why?
