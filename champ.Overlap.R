# A script to overlap two Genome-Based segmenets. For example, peaks (CpGs) with genome annotation (CpGIslands, Histone Position)
# Author: Tian

library("GenomicRanges")
library("data.table")

ConcatCharacterArray <- function(x) {
    paste(x, collapse=";")
}

champ.Overlap <- function(x,
                          y,
                          yAttributes=c("feature", "transcriptClass", "geneSymbol"),
                          collapseFunctions=list(ConcatCharacterArray, ConcatCharacterArray, ConcatCharacterArray))
{
    if(!all(yAttributes %in% colnames(y))) stop("At least one yAttributes does not exist in y data frame.")

    x.gr <- makeGRangesFromDataFrame(x)
    y.gr <- makeGRangesFromDataFrame(y)

    message("Finding overlap between x and y.")
    ov <- GenomicRanges::findOverlaps(x.gr, y.gr)

    ovY.dt <- setDT(y[subjectHits(ov), ])
    ovY.dt[, xIndex := queryHits(ov)]
    ovY.dt[, yIndex := subjectHits(ov)]

    message("Collapse y for each x.")
    ovXResult.dt <- ovY.dt[, .(.N, yIndexList = paste(yIndex, collapse=";")), by = .(xIndex)]

    if(length(yAttributes) > 0) {
        for(i in 1:length(yAttributes)) {
            message("Collapsing ", yAttributes[i], " ...")
            ovY.dt[, "tmpColumn" := ovY.dt[, get(yAttributes[i])]]
            tmp.dt <- ovY.dt[, .(collapseFunctions[[i]](tmpColumn)), by = .(xIndex)]
            ovXResult.dt <- merge(ovXResult.dt, tmp.dt)
        }
    }

    message("Prepare data frame for return.")
    ovXResult <- as.data.frame(ovXResult.dt)
    colnames(ovXResult) <- c("xIndex", "yOvNumber", "yIndexList", yAttributes)
    ovXResult <- cbind(x[ovXResult$xIndex, ], ovXResult)

    return(ovXResult)
}
