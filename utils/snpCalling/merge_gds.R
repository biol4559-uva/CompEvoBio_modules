# ijob -A biol4559-aob2x -c10 -p standard --mem=20G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(SeqArray)
  library(data.table)

### first, convert ExpEvo file
  seqParallelSetup(cluster=10, verbose=TRUE)

  vcf.fn="/scratch/aob2x/compBio_SNP_25Sept2023/dest.PoolSeq.PoolSNP.001.50.25Sept2023.norep.ann.vcf.gz"
  gds.fn="/scratch/aob2x/compBio_SNP_25Sept2023/dest.PoolSeq.PoolSNP.001.50.25Sept2023.norep.ann.gds"
  seqVCF2GDS(vcf.fn, gds.fn, storage.option="ZIP_RA", parallel=10, verbose=T, optimize=T)

### merge
  seqMerge(gds.fn=c("/scratch/aob2x/compBio_SNP_25Sept2023/dest.PoolSeq.PoolSNP.001.50.25Sept2023.norep.ann.gds",
                    "/scratch/aob2x/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds"),
           out.fn="/scratch/aob2x/dest.expEvo.8Jun2023.25Sept2023.norep.ann.gds")


### export
seqGDS2VCF("/scratch/aob2x/compBio_SNP_25Sept2023/dest.PoolSeq.PoolSNP.001.50.25Sept2023.norep.ann.gds",
            "/scratch/aob2x/expEvo.PoolSeq.PoolSNP.001.50.25Sept2023.norep.ann.vcf", info.var=NULL, fmt.var=NULL, chr_prefix="",
  use_Rsamtools=F, verbose=TRUE)

  seqGDS2VCF("/scratch/aob2x/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds",
              "/scratch/aob2x/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.vcf", info.var=NULL, fmt.var=NULL, chr_prefix="",
    use_Rsamtools=F, verbose=TRUE)

samtools

module load htslib/1.10.2 bcftools/1.9 intel/18.0 intelmpi/18.0 parallel/20200322 R/3.6.3 samtools vcftools

bgzip /scratch/aob2x/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.vcf
bgzip -@10 /scratch/aob2x/expEvo.PoolSeq.PoolSNP.001.50.25Sept2023.norep.ann.vcf
tabix -p vcf /scratch/aob2x/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.vcf.gz
tabix -p vcf /scratch/aob2x/expEvo.PoolSeq.PoolSNP.001.50.25Sept2023.norep.ann.vcf.gz

tabix -p vcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz



bcftools merge -0 \
-o /scratch/aob2x/compBio_SNP_25Sept2023/dest.expevo.PoolSNP.001.50.25Sept2023.norep.ann.vcf.gz \
/scratch/aob2x/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.vcf.gz \
/scratch/aob2x/expEvo.PoolSeq.PoolSNP.001.50.25Sept2023.norep.ann.vcf.gz

# ijob -A biol4559-aob2x -c10 -p standard --mem=20G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(SeqArray)
  library(data.table)

### first, convert ExpEvo file
  seqParallelSetup(cluster=10, verbose=TRUE)

  vcf.fn="/scratch/aob2x/compBio_SNP_25Sept2023/dest.expevo.PoolSNP.001.50.11Oct2023.norep.ann.vcf"
  gds.fn="/scratch/aob2x/compBio_SNP_25Sept2023/dest.expevo.PoolSNP.001.50.11Oct2023.norep.ann.gds"
  seqVCF2GDS(vcf.fn, gds.fn, storage.option="ZIP_RA", parallel=10, verbose=T, optimize=T)
