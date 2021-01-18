# This is a script I wrote to do TFEA between TFregulomeR and DhMR from Mouse Intestine

library(GenomicRanges)

champ.TFEA <- function(myROI, myRandomRegion, myTFPeaks)
{
    message("Number of myROI: ", nrow(myROI))
    message("Number of myRandomRegion: ", nrow(myRandomRegion))

    FisherExactTest <- function(myROI, myRandomRegion, TFPeak)
    {
        if(nrow(myROI) != nrow(myRandomRegion)) stop("Something is wrong, the row number should be same for ROI and RandomRegion.")

        ROI <- makeGRangesFromDataFrame(myROI, keep.extra.columns=TRUE)
        myRandomRegion <- makeGRangesFromDataFrame(myRandomRegion, keep.extra.columns=TRUE)
        
        Regulome.gr <- makeGRangesFromDataFrame(TFPeak, keep.extra.columns=TRUE)
        Regulome.gr <- GenomicRanges::promoters(Regulome.gr, 250, 250)
        
        ovROI <- data.frame(findOverlaps(Regulome.gr, ROI))
        ovRandom <- data.frame(findOverlaps(Regulome.gr, myRandomRegion))
        
        joinTable <- c(length(unique(ovROI$queryHits)), length(unique(ovRandom$queryHits)))
        joinTable <- rbind(joinTable, nrow(TFPeak) - joinTable)
        pValue <- fisher.test(joinTable)$p.value
        colnames(ovROI) <- c("TFPeakHits", "ROIHits")
        return(list(ovROI=ovROI, pValue=pValue))
    }
    
    message("Calculating Fisher Exact Test for each TF.")
    myTFEA <- lapply(myTFPeaks$TFPeaks, function(x) FisherExactTest(myROI, myRandomRegion[sample(1:nrow(myRandomRegion),nrow(myROI)), ], x))


    message("Organise TF List (TFList) and PeakOverlap (PeakOV) for return.")
    TFList <- myTFPeaks$TFList

    keyColumns <- do.call("rbind", lapply(myTFEA, function(x) c(length(unique(x$ovROI$TFPeakHits)), length(unique(x$ovROI$ROIHits)), x$pValue)))
    colnames(keyColumns) <- c("TFPeakHits", "ROIHits", "pValue")
    TFList <- cbind(TFList, keyColumns)
    TFList <- TFList[order(TFList$pValue), ]

    PeakOV <- myTFEA[TFList$ID]

    return(list(TFList=TFList, PeakOV=PeakOV))
}
