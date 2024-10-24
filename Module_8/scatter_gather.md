# **Scatter-gather**

## Overview

Sometimes we are impatient and want to run our sliding window analysis faster. The easiest way to do this is to use a 'scatter-gather' approach where we use SLURM to run our window script on a subset of the genome. We can run dozens of jobs simultaneously and quickly speed up the analsyis. This is especially important the longer it takes to run whatever analysis you are interested in for your particular window. The scatter phase is using SLURM to launch a bunch of jobs using the 'sbatch --array=1-X' flag in we went over in Module 3, [SLURM_demo_3.sh](https://github.com/biol4559-uva/CompEvoBio_modules/blob/main/Module_3/SLURM_demo_3.sh). We'll want to divide the analysis into about 100 jobs. The windows file is about 2.1K windows long, so we'll aim for about 20 windows per job. The gather phase is collecting the results into a single object and can be run as an interactive job.

To pass the ${SLURM_ARRAY_TASK_ID} to R using this command in your slurm script:
```
Rscript --vanilla FULL_PATH_TO_R_SW_SCRIPT ${SLURM_ARRAY_TASK_ID}
```

To read the ${SLURM_ARRAY_TASK_ID} task ID into your R script, put this at the top of your script:
```
args = commandArgs(trailingOnly=TRUE)
jobId=as.numeric(args[1])
```

In your scatter phases, you'll have to save the window slices to a file using this command:
```
save(output, file=paste("/scratch/USERNAME/sliding_window_output/", jobId, ".Rdata", sep=""))
```

Goals:
1. Learn how to speed up computation with embarrassingly simple parallelization.

Objective:
1. Figure out how to divide up your `windows.to.use` object so that 100 jobs run across ~2.1K windows.

Here is an example [gather script](/CompEvoBio_modules/Module_8/gather.R)
