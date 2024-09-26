# ijob -A berglandlab_standard -c1 -p standard --mem=6G
### module load gcc/11.4.0  openmpi/4.1.4 R/4.3.1; R

### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(googlesheets4)

### get files
  dirs <- system("ls -d /standard/BerglandTeach/mapping_output/*", intern=T)

### loop
  output <- foreach(i=dirs)%do%{
    # i <- dirs[1]
    full.files <- list.files(i, full.name=T)
    file.names <- list.files(i, full.name=F)
    file.size <- round(file.size(full.files))
    data.table(file.name=file.names, file.size=file.size, proj=tstrsplit(file.names, "_")[[2]],
              sample=paste(tstrsplit(file.names, "_")[[3]], tstrsplit(file.names, "_")[[4]], sep="_"))
  }
  output <- rbindlist(output)
  output.ag <- output[,list(nFiles=length(file.name), origBamSize=file.size[grepl("original.bam", file.name)], syncSize=file.size[grepl("sync.gz", file.name) & !grepl("tbi", file.name)]), list(proj, sample)]
  output.ag.ag <- output.ag[,list(propSuccessful=mean(nFiles==14), minBamSize=min(origBamSize, na.rm=T), minSyncSize=min(syncSize, na.rm=T),
                                  maxBamSize=max(origBamSize, na.rm=T), maxSyncSize=max(syncSize, na.rm=T)), list(proj)]
  output.ag.ag


  foreach(proj.i=output.ag.ag$proj)%do%{

    sheet_write(ss="https://docs.google.com/spreadsheets/d/1S6YrUMkOXBIqAMmgriGd-WVS0OcnEcz2SHT3L3Ps4Xc/edit?usp=sharing",
                sheet=proj.i,
                output[proj==proj.i])
  }




  cat -v file.txt
  sed -i 's/[^[:print:]\t]//g' file.txt
  cat -v file.txt
  nano file.txt
