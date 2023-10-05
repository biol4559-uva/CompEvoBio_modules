# ijob -A berglandlab_standard -c10 -p standard --mem=40G
#module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3; R

library(SeqArray)


args = commandArgs(trailingOnly=TRUE)
vcf.fn=args[[1]]
gds.fn=gsub(".vcf", ".gds", vcf.fn)

#vcf.fn=paste(vcf.fn, ".gz", sep="")
#vcf.fn="dest.all.PoolSNP.001.5.test.ann.vcf"
seqParallelSetup(cluster=10, verbose=TRUE)

seqVCF2GDS(vcf.fn, gds.fn, storage.option="ZIP_RA", parallel=10, verbose=T, optimize=T)



#seqVCF2GDS("/project/berglandlab/DEST/vcf/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.vcf.gz",
#            "/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.vgds", storage.option="ZIP_RA", verbose=T, parallel=10, optimize=T)