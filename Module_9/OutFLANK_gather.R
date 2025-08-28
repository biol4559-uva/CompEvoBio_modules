
### libraries
  library(curl)
  library(SeqArray)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(foreach)
  library(doMC)
  registerDoMC(10)

### source a few functions
  source("https://raw.githubusercontent.com/j-a-thia/genomalicious/master/R/outflank_input_freqs.R")
  source("https://raw.githubusercontent.com/j-a-thia/genomalicious/master/R/outflank_mod_fst_correct.R")
  source("https://raw.githubusercontent.com/j-a-thia/genomalicious/master/R/outflank_mod_fst_nocorrect.R")
  #source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/OutFLANK.R")
  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/Likelihood%20functions%20for%20OutFLANK.R")
  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/Fst%20Diploids.R")
  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/FST%20functions.R")
  source("https://raw.githubusercontent.com/biol4559-uva/CompEvoBio_modules/refs/heads/main/Module_9/OutFLANK.R")

### collect
  fl <- list.files("/scratch/aob2x/outflank_output", "window.Rdata", full.names=T)
  out <- foreach(fl.i=fl)%dopar%{
    load(fl.i)
    return(out)
  }
  out <- rbindlist(out, fill=T)
  setnames(out, "meanAF", "meanAlleleFreq")
  setnames(out, "variant.id", "LocusName")

### save big object
  save(out, file="/scratch/aob2x/outFLANK.Rdata")
  load("/scratch/aob2x/outFLANK.Rdata")

### outlier detection
  tmp <- out[!is.na(FST)][FST>0 & FST<1][!is.na(FSTNoCorr)]

  of <- OutFLANK(as.data.frame(tmp), LeftTrimFraction=0.05, RightTrimFraction=0.05,
                  Hmin=0.01, NumberOfSamples=2, qthreshold=0.05)
  of <- as.data.table(of)

### summarize
  ggplot(data=of, aes(results.pvaluesRightTail)) + geom_histogram()

  table(of$results.OutlierFlag)
  of[which.min(results.pvaluesRightTail)]

  ggplot(data=of, aes(x=results.FST, y=log10(results.pvaluesRightTail))) + geom_point()

### some basic plots
  ggplot(data=of, aes(x=results.pos, y=results.FST, color=results.chr)) + geom_point() + facet_grid(~results.chr)

### a basic enrichment test for non-synonymoyus SNPs. What other types of categories are enriched?
  tab.ns <- table(of$results.pvaluesRightTail<.1, of$results.annotation=="missense_variant")
  fisher.test(tab.ns)
  tab.ns

### what candidate genes come up as top candidates? Check them out on Flybase.
  of[results.pvaluesRightTail<.06][results.annotation=="missense_variant"]$results.gene

### are top variants clinal? This uses a clinal analysis from Machado et al 2021 eLife
  clinal <- fread(file="/standard/BerglandTeach/clinal_table.csv")
  setnames(of, c("results.chr", "results.pos"), c("chr", "pos"))
  of.clinal <- merge(of, clinal, by=c("chr", "pos"))
  tab.clinal <- table(of.clinal$results.pvaluesRightTail<.05, of.clinal$clinal.p<.05)
  fisher.test(tab.clinal)

### are top variants seasonal? This uses a seasonal analysis from Machado et al 2021 eLife
  seasonal <- fread(file="/standard/BerglandTeach/seasonal_table.csv")
  of.seas <- merge(of, seasonal, by=c("chr", "pos"))
  tab.seas <- table(of.seas$results.pvaluesRightTail<.05, of.seas$seas.p<.05)
  fisher.test(tab.seas)

### maybe it will be important to look world-wide allele frequency change of individual SNPs that are candidates?

  ### load this function
    source("https://raw.githubusercontent.com/biol4559-uva/CompEvoBio_modules/refs/heads/main/utils/misc/getData_function.R")

  ### open GDS file
    genofile <- seqOpen("/scratch/aob2x/dest.expevo.PoolSeq.PoolSNP.001.50.28Sept2024_ExpEvo.norep.gds")

  ### load meta-data file
    samps <- fread("https://raw.githubusercontent.com/biol4559-uva/CompEvoBio_modules/refs/heads/main/data/full_sample_metadata.28Sept2024_ExpEvo.csv")

  ### common SNP.dt
    seqResetFilter(genofile)
    snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                         pos=seqGetData(genofile, "position"),
                         nAlleles=seqGetData(genofile, "$num_allele"),
                         id=seqGetData(genofile, "variant.id"))

    snp.dt <- snp.dt[nAlleles==2] ### subset to sites with only two alleles
    seqSetFilter(genofile, snp.dt$id)

  ### extract out some candidates
    setnames(of, "results.LocusName", "id")

    worldwide <- getData(snps=of[results.pvaluesRightTail<.05][1:10], samples=samps)

    ggplot(data=worldwide, aes(x=lat, y=af_nEff, color=continent)) + geom_point() + facet_wrap(~variant.id)
    ggplot(data=worldwide, aes(x=jday, y=af_nEff)) + geom_point()  + facet_wrap(~variant.id)
    ggplot(data=worldwide, aes(x=jday, y=af_nEff)) + geom_point()  + facet_wrap(~variant.id)


### you can also pretty up the genome-wide Fst result by doing a sliding window summarization of the Fst & p-values.
### There are a few ways to go about this: You can take the average Fst value in a window. You can summarize the p-values by combining them together.
### one way to combine p-values is to take their product, which after a quick transformation follows a chi-sq distribution.

  ### build a dictionary of windows
    win.bp <- 1e5
    step.bp <- 5e4

    setkey(of, "chr")

    wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),.combine="rbind", .errorhandling="remove")%dopar%{
      # chr.i <- "2L"
      tmp <- of[J(chr.i)]
      data.table(chr=chr.i,
                 start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                 end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
    }

    wins[,window:=1:dim(wins)[1]]
    dim(wins)

  ### wZa function for combining p-values
    wZa.fun <- function(p, h, ret="wZa") {
      Z <- qnorm(p, 0, 1)
      wZa=sum(h*Z, na.rm=T)/(sqrt(sum(h^2, na.rm=T)))
      wZa.p=pnorm(sum(h*Z, na.rm=T)/(sqrt(sum(h^2, na.rm=T))), lower.tail=T)
      if(ret=="wZa") return(wZa)
      if(ret=="wZa.p") return(wZa.p)
    }


  ### run windows
    setkey(of, chr, pos)

    window.out <- foreach(i=1:dim(wins)[1], .errorhandling="remove", .combine="rbind")%dopar%{
      #i <- 100
      message(i)
      ### get data for your sample for this window
        focalSNPs <- of[J(wins[i]$chr)][pos>=wins[i]$start & pos<=wins[i]$end]

      ### Now we summarize.
        window.ag <- focalSNPs[,list(wZa.p=wZa.fun(p=results.pvaluesRightTail, h=results.He, ret="wZa.p"),
                                        wZa=wZa.fun(p=results.pvaluesRightTail, h=results.He, ret="wZa"),
                                        Fst_average=mean(results.FST),
                                        start=min(pos), end=max(pos), nSNPs=length(results.FST)),
                                    list(chr)]

      ### return
        return(window.ag)
    }

  ### your turn - plot it
    -log10(wZa.p)

    library(rgbif)
    coords <- data.frame(decimalLatitude = 54.481084,
                         decimalLongitude = -3.220625)

    elevation(coords, username = "biol4020")

    ll <- samps[,c("long", "lat", "sampleId")]
    setnames(ll, c("long", "lat"), c("x", "y"))
    ll <- na.omit(ll)
    df_elev_epqs <- get_elev_point(as.data.frame(ll), prj = 4326, src = "aws")
    sampse <- merge(samps, df_elev_epqs, by="sampleId")
    sampse[elevation<0, elevation:=0]

    write.csv(sampse, quote=F, row.names=F, file="~/samps_with_elevation.csv")


    set.seed(65.7)
examp_df <- data.frame(x = runif(3, min = -73, max = -72.5), y = runif(3, min = 42,
    max = 43))
crs_dd <- 4326

# Create and example data.frame with additional columns
cats <- data.frame(category = c("H", "M", "L"))

examp_df2 <- data.frame(examp_df, cats)

# Create an example
examp_sf <- sf::st_as_sf(examp_df2, coords = c("x", "y"), crs = crs_dd)
df_elev_epqs <- get_elev_point(examp_df, prj = crs_dd, src = "epqs")
