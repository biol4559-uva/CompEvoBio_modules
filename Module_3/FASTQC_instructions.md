

### Modules & fastq quality control
1. Before we begin analyzing FASTQ data, we want to understand some basics of its quality. For instance, we can ask how base quality scores change along the sequence read. Or, we can ask what the proportions of A, C, T, and G are in the sample. These sorts of assessments can help us identify if the sequencing run "failed", or if it is contaminated with somethign unexpected. (different species have different CG proportions).

2. We will perform this task using two tools. The first is called `fastqc` and the second is `multiqc`. Fastqc takes a single fastq file and runs some analyses. Multiqc takes the output of a bunch of fastqc files and produces a summary report. We will use `fastqc` now as a demonstration, and `multiqc` for your final paper's data

3. To use these tools, we need to load them using the `module` command. On Rivanna, most programs need to be loaded first, much like you need to load libraries in R. We also need to make an output folder. Copy these lines to a blank Rscript file; modify the COMPUTEID value for your computeID:
```
module load fastqc  ### this load the module

proj=SRP002024 ### defines the (Short Read Archive) project identifier. (the paper)
srr=SRR036932 ### defines the SRA sample identifier (the unique sample)

if [ ! -d /scratch/COMPUTEID/fastq_QC_Out/ ]; then ### these three lines will make an output folder
  mkdir /scratch/COMPUTEID/fastq_QC_Out/
fi

fastqc -o /scratch/COMPUTEID/fastq_QC_Out --noextract /standard/BerglandTeach/data/fastq/${proj}/${srr}.fastq.gz

```

4. Download the contents of `/scratch/COMPUTEID/fastq_QC_Out` to your computer using the OpenOnDemand file transfer tool and open the html file. What do you see?
