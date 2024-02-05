# A R function to Import IDAT file without SampleSheet
# Author: Tian

# manifest_directory <- "../../ChAMP-Script/Illumina_Methylation_Manifests/MethylationEPIC_v-1-0_B4.csv"

directory <- "~/Downloads/TCGA-PAAD/"

library("ChAMP")
library("illuminaio")
library("tidyverse")

champ.ImportIDAT <- function(directory = getwd(), offset = 100, arraytype="450K")
{
  message("[=================================]")
  message("[<<<< ChAMP.IMPORTIDAT START >>>>>]")
  message("----------------------------------")
  
  idat_files <- dir(directory, pattern = "*.idat$", recursive = TRUE, full.names = TRUE)
  
  GrnPath <- idat_files[str_detect(idat_files, "_Grn.idat")] %>% sort %>% .[1:10]
  RedPath <- idat_files[str_detect(idat_files, "_Red.idat")] %>% sort %>% .[1:10]
  
  G.idats <- lapply(GrnPath, function(x){ message("  Loading:",x," ---- (",which(GrnPath == x),"/",length(GrnPath),")");readIDAT(x)})
  R.idats <- lapply(RedPath, function(x){ message("  Loading:",x," ---- (",which(RedPath == x),"/",length(RedPath),")");readIDAT(x)})
  
  names(G.idats) <- basename(GrnPath) %>% gsub("_Grn.idat", "", .)
  names(R.idats) <- basename(RedPath) %>% gsub("_Red.idat", "", .)
  
  checkunique <- unique(c(sapply(G.idats, function(x) nrow(x$Quants)),sapply(R.idats, function(x) nrow(x$Quants))))
  
  if(length(checkunique) > 1) 
  {
    message("\n  !!! Important !!! ")
    message("  Seems your IDAT files not from one Array, because they have different numbers of probe.")
    message("  ChAMP wil continue analysis with only COMMON CpGs exist across all your IDAt files. However we still suggest you to check your source of data.\n")
  }
  
  CombineIDAT <- append(G.idats, R.idats)
  commonAddresses <- as.character(Reduce("intersect", lapply(CombineIDAT, function(x) rownames(x$Quants))))
  
  message("\n  Extract Mean value for Green and Red Channel Success")
  GreenMean <- do.call(cbind, lapply(G.idats, function(xx) xx$Quants[commonAddresses, "Mean"]))
  RedMean <- do.call(cbind, lapply(R.idats, function(xx) xx$Quants[commonAddresses, "Mean"]))
  message("    Your Red Green Channel contains ",nrow(GreenMean)," probes.")
  
  message("\n  Extract Mean value for Green and Red Channel Success")
  GreenMean <- do.call(cbind, lapply(G.idats, function(xx) xx$Quants[commonAddresses, "Mean"]))
  RedMean <- do.call(cbind, lapply(R.idats, function(xx) xx$Quants[commonAddresses, "Mean"]))
  message("    Your Red Green Channel contains ",nrow(GreenMean)," probes.")
  
  G.Load <- do.call(cbind,lapply(G.idats,function(x) x$Quants[commonAddresses,"Mean"]))
  R.Load <- do.call(cbind,lapply(R.idats,function(x) x$Quants[commonAddresses,"Mean"]))
  
  message("\n  Reading ", arraytype, " Annotation >>")
  if(arraytype %in% c("EPIC", "EPICv2")) {
    message("    !!! Important, since version 2.29.1, ChAMP set default `EPIC` arraytype as EPIC version 2. ",
            "\n        You can set 'EPIC' or 'EPICv2' to use version 2 EPIC annotation",
            "\n        If you want to use the old version (v1), please specify arraytype parameter as `EPICv1`. ",
            "\n        For 450K array, still use `450K`")
    data("AnnoEPICv2")
  } else if(arraytype == "EPICv1") {
    data("AnnoEPICv1")
  } else { 
    data(Anno450K)
  }
  
  message("\n  Fetching NEGATIVE ControlProbe.")
  control_probe <- rownames(Anno$ControlProbe)[which(Anno$ControlProbe[,1]=="NEGATIVE")]
  message("    Totally, there are ",length(control_probe)," control probes in Annotation.")
  control_probe <- control_probe[control_probe %in% rownames(R.Load)]
  message("    Your data set contains ",length(control_probe)," control probes.")
  rMu <- matrixStats::colMedians(R.Load[control_probe,])
  rSd <- matrixStats::colMads(R.Load[control_probe,])
  gMu <- matrixStats::colMedians(G.Load[control_probe,])
  gSd <- matrixStats::colMads(G.Load[control_probe,])
  
  rownames(G.Load) <- paste("G",rownames(G.Load),sep="-")
  rownames(R.Load) <- paste("R",rownames(R.Load),sep="-")
  
  IDAT <- rbind(G.Load,R.Load)
  
  message("\n  Generating Meth and UnMeth Matrix")
  
  message("    Extracting Meth Matrix...")
  M.check <- Anno$Annotation[,"M.index"] %in% rownames(IDAT)
  message("      Totally there are ",nrow(Anno$Annotation)," Meth probes in ",arraytype," Annotation.")
  message("      Your data set contains ",length(M.check), " Meth probes.")
  M <- IDAT[Anno$Annotation[,"M.index"][M.check],]
  
  message("    Extracting UnMeth Matrix...")
  U.check <- Anno$Annotation[,"U.index"] %in% rownames(IDAT)
  message("      Totally there are ",nrow(Anno$Annotation)," UnMeth probes in ",arraytype," Annotation.")
  message("      Your data set contains ",length(U.check), " UnMeth probes.")
  U <- IDAT[Anno$Annotation[,"U.index"][U.check],]
  
  if(!identical(M.check,U.check))
  {
    stop("  Meth Matrix and UnMeth Matrix seems not paried correctly.")
  } else {
    CpG.index <- Anno$Annotation[,"CpG"][M.check]
  }
  
  rownames(M) <- CpG.index
  rownames(U) <- CpG.index
  
  message("\n  Generating beta Matrix")
  BetaValue <- M / (M + U + offset)
  message("  Generating M Matrix")
  MValue <- log2(M/U)
  message("  Generating intensity Matrix")
  intensity <-  M + U
  
  message("  Calculating Detect P value")
  detP <- matrix(NA,nrow=nrow(intensity),ncol=ncol(intensity))
  rownames(detP) <- rownames(intensity)
  colnames(detP) <- colnames(intensity)
  
  type_II <- rownames(Anno$Annotation)[Anno$Annotation[,"Channel"] == "g+r"]
  type_II <- type_II[type_II %in% rownames(detP)]
  type_I.red <- rownames(Anno$Annotation)[Anno$Annotation[,"Channel"] == "r"]
  type_I.red <- type_I.red[type_I.red %in% rownames(detP)]
  type_I.grn <- rownames(Anno$Annotation)[Anno$Annotation[,"Channel"] == "g"]
  type_I.grn <- type_I.grn[type_I.grn %in% rownames(detP)]
  for(i in 1:ncol(detP))
  {
    detP[type_II,i] <- 1 - pnorm(intensity[type_II,i], mean=rMu[i]+gMu[i], sd=rSd[i]+gSd[i])
    detP[type_I.red,i] <- 1 - pnorm(intensity[type_I.red,i], mean=rMu[i]*2, sd=rSd[i]*2)
    detP[type_I.grn,i] <- 1 - pnorm(intensity[type_I.grn,i], mean=gMu[i]*2, sd=gSd[i]*2)
  }
  if(sum(is.na(detP))) message("    !!! There are NA values in your detP matrix.\n")
  
  
  message("  Counting Beads")
  NBeads <- do.call(cbind, lapply(R.idats, function(x) x$Quants[commonAddresses, "NBeads"]))
  Mbead <- NBeads[substr(Anno$Annotation$M.index[M.check],3,100),]
  Ubead <- NBeads[substr(Anno$Annotation$U.index[U.check],3,100),]
  Ubead[Ubead < 3 | Mbead < 3] <- NA
  rownames(Ubead) <- rownames(intensity)

  message("\n[<<<<< ChAMP.IMPORTIDAT END >>>>>>]")
  message("[===========================]")
  message("[You may want to process champ.filter() next.]\n")
  
  return(list("beta"=BetaValue,"M"=MValue,"intensity"=intensity,"detP"=detP,"beadcount"=Ubead,"Meth"=M,"UnMeth"=U))
}