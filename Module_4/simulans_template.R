### libraries
  .libPaths(c("/project/biol4559-aob2x/biol4559-R-packages/", .libPaths()))
  library(ggplot2)
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(2)
  library(Rsamtools)

### specify the bam file
  fl <- system("ls -d /project/biol4559-aob2x/data/Contamination_test/mapping_output/*/*.srt.flt.bam", intern=T)

### How many reads map to each chromosome?
  rd <- foreach(bamFile=fl, .combine="rbind")%dopar%{
    ### tell me what file we are working on
      message(bamFile)
    ### get the information about number of reads for each chromosome
      stats <- as.data.table(idxstatsBam(bamFile))

    ### we need to subset to the main autosomal arms of melanogaster and simulans
      setkey(stats, seqnames)
      stats_small <- stats[J(c("2L", "2R", "3L", "3R", "X", "sim_2L", "sim_2R", "sim_3L", "sim_3R", "sim_X"))]
      stats_small[grepl("sim", seqnames),species:="sim"]
      stats_small[!grepl("sim", seqnames),species:="mel"]

    ### this command aggregates the data. For each unique value of species, we calculate the total number of mapped reads.
      stats_small.ag <- stats_small[,list(mapped=sum(mapped), seqlength=sum(seqlength)), list(species)]

    ### format output
      out <- data.table(propSim= stats_small.ag[species=="sim"]$mapped/sum( stats_small.ag$mapped))

      out[,nSim:=as.numeric(tstrsplit(bamFile, "\\.")[[8]])]
      out[,nMel:=as.numeric(tstrsplit(bamFile, "\\.")[[10]])]
      out[,samp:=last(tstrsplit(bamFile, "/"))]

    ### return output
      return(out)
  }
  rd[,realProp:=nSim/(nSim+nMel)]

### estimated vs real
  ggplot(data=rd, aes(x=realProp, y=propSim)) + geom_point()
