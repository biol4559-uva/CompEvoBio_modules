getData <- function(snps=snp.dt[pos==14617051 & chr=="2L"], samples=samps) {
  # snps=snp.dt[pos==14617051 & chr=="2L"]; samples=samps[grepl("SRP002024", sampleId)]$sampleId

  ### filter to target
  seqSetFilter(genofile, variant.id=snps$id, sample.id=samples$sampleId)

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
  afi <- merge(afi, snps, by.x="variant.id", by.y="id")

  afi[,af:=ad/dp]

  ### calculate effective read-depth
  afis <- merge(afi, samples[,c("sampleId", "set", "nFlies", "locality",
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
