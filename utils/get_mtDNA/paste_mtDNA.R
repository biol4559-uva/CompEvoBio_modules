# ijob -A berglandlab_standard -c20 -p standard --mem=80G
### module load gcc/11.4.0  openmpi/4.1.4 R/4.3.1; R

### libraries
library(data.table)
library(foreach)
library(doMC)
registerDoMC(5)


tmpdir <- "/standard/BerglandTeach/mtDNA/sync"
### get input files
files <- list.files("/standard/BerglandTeach/mtDNA/sync", pattern=".sync.gz")
length(files)
setwd(tmpdir)

chrom_name = rev(tstrsplit(jobId, "_"))[[3]]
if (chrom_name == "genome") chrom_name = "mitochondrion_genome"

start = rev(tstrsplit(jobId, "_"))[[2]]
end = rev(tstrsplit(jobId, "_"))[[1]]

#files <- files[-1]

### import
o <- foreach(files.i=files, .errorhandling="pass")%dopar%{
  #files.i=files[10]
  tmp <- fread(files.i)
  if(dim(tmp)[1]==0) {
    tmp <- data.table(V1=chrom_name,
                      V2=start:end,
                      V3="N",
                      V4=".:.:.:.:.:.")
  }
  pop <- tstrsplit(files.i, "\\.")[[1]]
  tmp[,pop:=pop]
  tmp
}
o <- rbindlist(o, use.names=T, fill=T)
o[,pop:=gsub(".SNAPE.monomorphic", "", pop)]

dim(o)
o[,.N,pop]

### long to wide
ow <- dcast(o, V1+V2~pop, value.var="V4")

## get reference
popu <- unique(o$pop)
ow.ref <- o[pop==popu[1], c("V1", "V2", "V3"), with=F]

setkey(ow, V1, V2)
setkey(ow.ref, V1, V2)

owr <- merge(ow.ref, ow)

### output
write.table(owr, quote=F, row.names=F, col.names=F, sep="\t", file=paste(tmpdir, "/allpops.", "mtDNA", ".sites", sep=""))
write.table(names(owr)[-c(1,2,3)], quote=F, row.names=F, col.names=F, sep="\t", file=paste(tmpdir, "/allpops.", "mtDNA", ".sites", sep=""))
