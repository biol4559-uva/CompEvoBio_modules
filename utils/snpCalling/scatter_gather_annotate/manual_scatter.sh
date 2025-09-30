#!/usr/bin/env bash

#SBATCH -J manual_scatter # A single job name for the array
#SBATCH --ntasks-per-node=8 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 4:00:00 ### 1 hours
#SBATCH --mem 64G
#SBATCH -o /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_gather.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_gather.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard


# ijob -A berglandlab -c10 -p standard --mem=50G
# sbatch --array=1-1002 ~/CompEvoBio_modules/utils/snpCalling/scatter_gather_annotate/manual_scatter.sh
# sacct -j 4244858 | grep -v "COMPLE"
# cat /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_gather.4243100_3.out




module purge
trap 'rm -rf ${tmpdir}' EXIT

#module load htslib bcftools parallel intel/18.0 intelmpi/18.0 mvapich2/2.3.1 R/3.6.3 python/3.6.6 vcftools/0.1.16
#module load htslib/1.10.2 bcftools/1.9 parallel/20200322 intel/18.0 intelmpi/18.0 R/3.6.3 python/3.6.6 vcftools/0.1.16
module load htslib/1.17  bcftools/1.17 parallel/20200322 gcc/11.4.0 openmpi/4.1.4 python/3.11.4 vcftools/0.1.16 R/4.3.1
module load bedtools/2.30.0

### r, mvapch, parallel


## Run params
  popSet=all
  method=PoolSNP
  maf=001
  mac=50
  version=29Sept2025_ExpEvo
  wd=/scratch/aob2x/compBio_SNP_29Sept2025
  script_dir=~/CompEvoBio_modules/utils/snpCalling/
  pipeline_output=/project/berglandlab/DEST/dest_mapped/
  job=${SLURM_ARRAY_TASK_ID}    # job=1



  #### popSet="all"; method="poolSNP"; maf="001"; mac=5; jobs="jobs.csv"; script_dir="/scratch/aob2x/DESTv2/snpCalling"; wd="/scratch/aob2x/DESTv2_output"; SLURM_JOB_ID=1;
  #### pipeline_output="/project/berglandlab/DEST/dest_mapped/"; job="2L,1,10000"

## working & temp directory
  outdir="${wd}/sub_vcfs" #### outdir=${wd}"/sub_vcfs"
    if [ ! -d $outdir ]; then
        mkdir $outdir
    fi

## get list of SNYC files based on popSet & method
### full list
  echo "foo: "${job}
  echo "pipeline_output: "${pipeline_output}
  #echo $( ls -l ${pipeline_output}/*/*/*.sync.gz )
  ls -d ${pipeline_output}/*/*.sync.gz


## get job
  #cat ${script_dir}/scatter_gather_annotate/jobs_genome.csv
  job=$( cat ${script_dir}/scatter_gather_annotate/jobs_genome.csv | sed "${job}q;d" )
  jobid=$( echo ${job} | sed 's/,/_/g' )
  echo "jobid is " $jobid
  echo "job is " $job

## set up RAM disk
  [ ! -d /dev/shm/$USER/ ] && mkdir -p /dev/shm/$USER/  # added -p flag to avoid error when it isn't needed
  [ ! -d /dev/shm/$USER/${SLURM_JOB_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}
  tmpdir=/dev/shm/$USER/${SLURM_JOB_ID}

  echo "Temp dir is $tmpdir"

## get sub section
  subsection () {
    #set +e
    syncFile=${1}
    job=${2}
    jobid=$( echo ${job} | sed 's/,/_/g' )
    tmpdir=${3}

    pop=$( echo ${syncFile} | rev | cut -f1 -d'/' | rev | sed 's/.masked.sync.gz//g' )

    chr=$( echo $job | cut -f1 -d',' )
    start=$( echo $job | cut -f2 -d',' )
    stop=$( echo $job | cut -f3 -d',' )

    echo ${pop}_${jobid}

    #touch ${syncFile}.tbi

    tabix -b 2 -s 1 -e 2 \
    ${syncFile} \
    ${chr}:${start}-${stop} > ${tmpdir}/${pop}_${jobid}
    #set -e

  }
  export -f subsection

  echo "subset"
  if [[ "${method}" == "SNAPE" && "${popSet}" == "PoolSeq" ]]; then
    echo "SNAPE" ${method}
    parallel -j 4 subsection ::: $( ls ${pipeline_output}/*/*/*.masked.sync.gz | tr '  ' '\n' | grep "SNAPE" | grep "monomorphic" ) ::: ${job} ::: ${tmpdir}
  elif [[ "${method}" == "PoolSNP" && "${popSet}" == "all" ]]; then
    echo "PoolSNP" ${method}
    parallel -j 4 subsection ::: $( ls ${pipeline_output}/*/*/*.masked.sync.gz | tr '  ' '\n' | grep -v "SNAPE" ) ::: ${job} ::: ${tmpdir}
  elif [[ "${method}" == "PoolSNP" && "${popSet}" == "PoolSeq" ]]; then
    echo "PoolSNP" ${method}
    parallel -j 4 subsection ::: $( ls ${pipeline_output}/*/*/*.masked.sync.gz | tr '  ' '\n' | grep -v "SNAPE" | grep -v "DGN" ) ::: ${job} ::: ${tmpdir}
  fi




### paste function
  echo "paste"
  Rscript --no-save --no-restore ${script_dir}/scatter_gather_annotate/paste.R ${job} ${tmpdir} ${method}

### run through SNP calling
  echo "SNP calling"

  if [[ "${method}" == "SNAPE" ]]; then
    echo $method
    cat ${tmpdir}/allpops.${method}.sites | python ${script_dir}/PoolSNP/PoolSnp.py \
    --sync - \
    --min-cov 4 \
    --max-cov 0.95 \
    --miss-frac 0.5 \
    --min-count 0 \
    --min-freq 0 \
    --posterior-prob 0.9 \
    --SNAPE \
    --names $( cat ${tmpdir}/allpops.${method}.names |  tr '\n' ',' | sed 's/,$//g' )  > ${tmpdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf

  elif [[ "${method}"=="PoolSNP" ]]; then
    echo $method

    cat ${tmpdir}/allpops.${method}.sites | python ${script_dir}/PoolSNP/PoolSnp.py \
    --sync - \
    --min-cov 4 \
    --max-cov 0.95 \
    --min-count ${mac} \
    --min-freq 0.${maf} \
    --miss-frac 0.5 \
    --names $( cat ${tmpdir}/allpops.${method}.names |  tr '\n' ',' | sed 's/,$//g' )  > ${tmpdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf
  fi

### compress and clean up
  echo "compress and clean"
  # cp ${tmpdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf ${outdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf
  # cat ${tmpdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf | bgzip -c > ${outdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz
   cat ${tmpdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf | vcf-sort | bgzip -c > ${outdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz

  tabix -p vcf ${outdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz


  bedtools intersect -sorted -v -header \
  -b ${script_dir}/scatter_gather_annotate/repeat_bed/repeats.sort.bed.gz \
  -a ${outdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz |
  bgzip -c > \
  ${outdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz

  tabix -p vcf ${outdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz

  #echo "vcf -> bcf "
  #bcftools view -Ou ${tmpdir}/${jobid}.vcf.gz > ${outdir}/${jobid}.bcf

  rm -fr ${tmpdir} # used bash exit trap instead (line 6)

### done
  echo "done"
