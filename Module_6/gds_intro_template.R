### run once and change directory to your scratch.
  system("mkdir /scratch/aob2x/biol4559-R-packages-newer")
  .libPaths(c("/scratch/aob2x/biol4559-R-packages-newer")); .libPaths()
  install.packages(c("foreach", "doMC", "data.table", "ggplot2", "patchwork", "BiocManager"), lib="/scratch/aob2x/biol4559-R-packages-newer")
  BiocManager::install("SeqArray", lib="/scratch/aob2x/biol4559-R-packages-newer", force=T)

### libraries
  .libPaths(c("/scratch/aob2x/biol4559-R-packages-newer")); .libPaths()

  library(SeqArray)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(foreach)
  library(doMC)
  registerDoMC(10)

### load this function
  getData <- function(snps=snp.dt[pos==14617051 & chr=="2L"], samples=samps) {
    # snps=snp.dt[pos==14617051 & chr=="2L"]; samples=samps

    ### filter to target
    seqSetFilter(genofile, variant.id=snps$id)

    ### get annotations
    message("Annotations")
    tmp <- seqGetData(genofile, "annotation/info/ANN")
    len1 <- tmp$length
    len2 <- tmp$data

    snp.dt1 <- data.table(len=rep(len1, times=len1),
                          ann=len2,
                          id=rep(snps$id, times=len1))

    # Extract data between the 2nd and third | symbol
    snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
    snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]

    # Collapse additional annotations to original SNP vector length
    snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=","), gene=paste(gene, collapse=",")),
                          list(variant.id=id)]

    snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
    snp.dt1.an[,gene:=tstrsplit(snp.dt1.an$gene,"\\,")[[1]]]

    ### get frequencies
    message("Allele Freqs")

    ad <- seqGetData(genofile, "annotation/format/AD")
    dp <- seqGetData(genofile, "annotation/format/DP")

    if(class(dp)[1]!="SeqVarDataList") {
      dp <- list(data=dp)
    }


    af <- data.table(ad=expand.grid(ad$data)[,1],
                     dp=expand.grid(dp$data)[,1],
                     sampleId=rep(seqGetData(genofile, "sample.id"), dim(ad$data)[2]),
                     variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(ad$data)[1]))

    ### tack them together
    message("merge")
    afi <- merge(af, snp.dt1.an, by="variant.id")
    afi <- merge(afi, snp.dt, by.x="variant.id", by.y="id")

    afi[,af:=ad/dp]

    ### calculate effective read-depth
    afis <- merge(afi, samples[,c("sampleId", "nFlies", "locality",
                                  "lat", "long", "continent", "country", "province", "city",
                                  "min_day", "max_day", "min_month", "max_month", "year", "jday",
                                  "bio_rep", "tech_rep", "exp_rep", "loc_rep", "subsample", "sampling_strategy",
                                  "SRA_Accession"), with=F], by="sampleId")

    afis[chr=="X|Y", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
    afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
    afis[,af_nEff:=round(af*nEff)/nEff]
    setnames(afis, "col", "annotation")
    ### return
    afis[,-c("n"), with=F]
  }

### open GDS file
  genofile <- seqOpen("/scratch/aob2x/GDS/dest.expevo.PoolSNP.001.50.25Sept2023.norep.ann.gds")
  genofile

### load meta-data file
  samps <- fread("/scratch/aob2x/GDS/biol4559_sampleMetadata.csv")

### subset samps to only show records associated with your paper. One easy way to do this is to use the `grepl` command.
### This command returns a vector of TRUE or FALSE if the search pattern is found (or not)
### For example:
  samps[grepl("SRP002024", sampleId)]

### First, we need to extract from our GDS file a dictionary that contains the position information for every SNP in the dataset
### common SNP.dt
  seqResetFilter(genofile)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       nAlleles=seqGetData(genofile, "$num_allele"),
                       id=seqGetData(genofile, "variant.id"))

  snp.dt <- snp.dt[nAlleles==2] ### subset to sites with only two alleles
  seqSetFilter(genofile, snp.dt$id)
  snp.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data] ### this is the global average frequency

### YOUR TURN:
### make a plot of the global average Site Frequency Spectrum that partitions on chromosome. Use geom_density() & facet_grid()


### extract allele frequency data for specific samples and sites. What columns are there?
### The columns that we care about for frequencies are "af_nEff" and "nEff".
### "nEff" is the effective read depth, assuming pools of males. Is this a correct assumption for you? Do you need to modify the function?
### "af_nEff" is the allele frequency estimate rounded such that af_nEff * nEff is an integer.
  dat <- getData(snps=snp.dt[pos==14617051 & chr=="2L"], samples=samps[grepl("SRP002024", sampleId)])

### YOUR TURN:
### use FlyBase to find the full coordinate range of Adh, and extract all SNPs in the region
### Go to https://flybase.org. Type "AdH" into the "Jump to Gene" / "J2G" field.
### make a plot of allele frequencies across your samples. Hint: to connect a line for each SNP use the `group=variant.id` parameter in the `aes` command
### Make sure to rotate the x-axis annotations. Unsure how to do this? Ask Google: "rotate the x-axis annotations ggplot2"\
### Hint: to extract a range of SNPs from snp.dt, you can use `pos>X & pos<Y` inside the hard-brackets


### One basic "rule" in genetics is that mutations that cause amino-acid changes (non-synonymous) are likely to be deleterious.
### On the other hand, mutations in coding sequence of genes that do not cause amino acid changes (synonymous) are likely to be neutral.
### As a result, the average allele frequency of Non-Synonymous mutations should be lower than for Synonymous mutations.
### There should also be fewer mutations at NS sites compared to Syn sites.
### Do you see evidence of this?
###

  ### first, we are going to subset the SNP data to make our life easier. We will re-evaluate this later
    subsamp <- sort(sample(snp.dt$id, 10000)) ### to begin, we are going to work with a subset of data. This gets
    setkey(snp.dt, id) ### to quickly subset on these IDs, we are going to use an indexed join command. First, set the index (key)
    snp.dt[J(subsamp)] ### double check that it works.

  ### next, we are going to extract out allele frequency data
    dat <- getData(snps=snp.dt[J(subsamp)], samples=samps[grepl("SRP002024", sampleId)])
    dat ### The dimensions of the output should be 1000*(num samples). Are they?

  ### YOUR TURN:
  ### take the average frequency across the samples and mutation annotations (the "annotation" column).
  ### Use an aggregation command and calculate the `mean(af_nEff, na.rm=T)` frequency


  ### kinda slow, isn't it? And, 1000 SNPs is not that many out of 4.3M. That is going to make it hard to scale up to a larger number of SNPs.
  ### Does parallelization to make things go faster. Use a foreach() loop to iterate through your samples. use `%dopar%`

    system.time(dat <- getData(snps=snp.dt[J(subsamp)], samples=samps[grepl("SRP002024", sampleId)]))

    system.time(dat <- foreach(samp.i=samps[grepl("SRP002024", sampleId)]$sampleId)%dopar%{
        # samp.i <- "ExpEvo_SRP002024_ACO_1_1975-MM-DD"
        dat <- getData(snps=snp.dt[J(subsamp)], samples=samps[sampleId==samp.i])
        dat
      })
    dat <- rbindlist(dat)

  ### YOUR TURN. Rewrite the foreach loop to loop over each SNP for all samples in your data.



### YOUR TURN:
### make a composite plot
### make sure that axis labels are legible, that subplots are labeled with "A", "B", "C", etc. Arrange the plot in a way that looks good to you
