### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)
  library(dplyr)

### paths
  folder <- "/scratch/aob2x/DEST2_analysis/seasonality/glm_test_SEPT_29_2023/"
  fileStem <- "Rdata"

### get files
  fl <- list.files(folder, fileStem, full.names=T)

### collect completed jobs
  o <- foreach(fl.i=fl, .errorhandling="pass")%dopar%{
    # fl.i <- fl[1]
    message(paste(which(fl.i==fl), length(fl), sep=" / "))

    load(fl.i) ### what did you save your objects as? the line below uses the placeholder of oo

    oo[,invName:=case_when(
          chr=="2L" & pos >	2225744	 & pos < 13154180	 ~ "2Lt",
          chr=="2R" & pos >	15391154 & pos < 	20276334 ~ 	"2RNS",
          chr=="3R" & pos >	11750567 & pos < 	26140370 ~ 	"3RK",
          chr=="3R" & pos >	21406917 & pos < 	29031297 ~ 	"3RMo",
          chr=="3R" & pos >	16432209 & pos < 	24744010 ~ 	"3RP",
          chr=="3L" & pos >	3173046	 & pos < 16308841	 ~ "3LP",
          TRUE ~ "noInv")]
    return(oo)
    #oo[perm==0]

  }


### plot your sliding window
