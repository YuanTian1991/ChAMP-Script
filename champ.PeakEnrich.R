# This is the script to find gene-body or promoter related genes. Data is prepared from MAnorm2 software.
# Currently, this function only works for MAnorm2 preprocessed ChIP-seq data (merely developed for (h)MeDIP-seq). In the future, it may extend to other data types.
# Author: Tian

library("GenomicRanges")
library("data.table")

# load("../65.MAnormSmallBin/M.RData")

# CNT <- M[, c(1:3, grep("cnt", colnames(M)))]
# Occupancy <- M[, c(1:3, grep("occupancy", colnames(M)))]
# phenos <- c("NC", "TC", "LT", "NL")

PrepareAnno <- function(CNT, Occupancy, phenos) {
    
    message("Preparing Annotation Data")
    
    Anno <- CNT[, 1:3]
    colnames(Anno) <- c("chr", "start", "end")
    Anno$pos <- (Anno$start + Anno$end) / 2
    
    message("Calculate Average Intensity, Sum Occupancy, and Fraction of Occupancy for each 100bp bin.")
    for (i in phenos) {
        message("Calculate Mean value for ", i, "...")
        meanCNT <- rowMeans(CNT[, substr(colnames(CNT), 1, 2) == i])
        sumOccupancy <- rowSums(Occupancy[, substr(colnames(Occupancy), 1, 2) == i])
        tmpDF <- data.frame(meanCNT, sumOccupancy, sumOccupancy / sum(substr(colnames(Occupancy), 1, 2) == i))
        colnames(tmpDF) <- paste0(i, c("_cnt", "_ocp", "_fraction"))
        Anno <- cbind(Anno, tmpDF)
    }
    return(Anno)
}

champ.PeakEnrich <- function(SigPeaks, AllPeaks, CandidateRegion, SampleTime = 1000) {
  message("Overlap Sig Peak Set with Gene Set")
  ov <- suppressWarnings(GenomicRanges::findOverlaps(makeGRangesFromDataFrame(SigPeaks[, c("chr", "start", "end")]),
                                                     makeGRangesFromDataFrame(CandidateRegion[, c("chr", "start", "end")])))
  ov <- as.data.frame(ov)
  message("There are ", length(unique(ov$queryHit)), " peaks mapped on ", length(unique(ov$subjectHits)), " regions.")

  SigPeakOV.dt <- setDT(SigPeaks[ov$queryHit, c("chr", "start", "end", "avgValue", "avgFraction")])
  SigPeakOV.dt[, cl := ov$subjectHits]
  message("Calculated sum for CandidateRegions")
  Collapse.dt <- SigPeakOV.dt[, .(chr[1], min(start), max(end), sum(avgValue), .N, sum(avgFraction >= 0.75)), by = .(cl)]
  colnames(Collapse.dt) <- c("symbol", "chr", "start", "end", "avgSum", "count", "countValid")
  Collapse.dt[, symbol := CandidateRegion[symbol, "symbol"]]
  Collapse <- as.data.frame(Collapse.dt)

  RandomMatrix <- list()
  for (i in 1:(SampleTime / 50)) {
    message("Generating Null Distribution Data: ", i * 50, "/", SampleTime)
    RandomMatrix[[i]] <- rowsum(matrix(sample(AllPeaks, nrow(SigPeakOV.dt) * 50, replace = TRUE), ncol = 50), SigPeakOV.dt$cl)
  }

  nulldistribution <- do.call("cbind", RandomMatrix)
  pValGreater <- 1 - (rowSums(matrix(rep(Collapse$avgSum, SampleTime), nrow = nrow(Collapse)) > nulldistribution) / SampleTime)
  pValLess <- 1 - (rowSums(matrix(rep(Collapse$avgSum, SampleTime), nrow = nrow(Collapse)) < nulldistribution) / SampleTime)
  adjPValGreater <- p.adjust(pValGreater, "BH")
  adjPValLess <- p.adjust(pValLess, "BH")
  Collapse <- data.frame(Collapse, data.frame(pValGreater = pValGreater, adjPValGreater = adjPValGreater, pValLess = pValLess, adjPValLess = adjPValLess))

  OverlapStatus <- ov
  OverlapStatus$subjectHits <- CandidateRegion[ov$subjectHits, "symbol"]

  return(list(Collapse = Collapse, OverlapStatus = OverlapStatus))
}


