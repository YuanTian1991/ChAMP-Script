# A script to convert UCSC MySQL cpgIslandsExt.
# Author: Tian

library(RMySQL)
library(GenomicRanges)
library(stringr)

champ.CpGIslands <- function()
{
    message("Extract hg38 cpgIslandsExt from UCSC MySQL database.")
    con_ucsc <- dbConnect(RMySQL::MySQL(), db = "hg38", user = "genome", host = "genome-mysql.soe.ucsc.edu")
    cpgIslandExt <- suppressWarnings(dbGetQuery(con_ucsc, stringr::str_interp("SELECT * FROM cpgIslandExt")))
    dbDisconnect(con_ucsc)

    message("Prepare CpG Islands")
    cpgIslands <- cpgIslandExt[,2:7]
    cpgIslands$index <- paste0("cgi", 1:nrow(cpgIslands))
    cpgIslands <- makeGRangesFromDataFrame(cpgIslands, keep.extra.columns=TRUE)
    cpgIslands$feature <- "cpgIsland"

    message("Prepare CpG Shores")
    cpgShoresNorth <- cpgShoresSouth <- cpgIslands
    end(cpgShoresSouth) <- start(cpgShoresSouth)
    start(cpgShoresSouth) <- start(cpgShoresSouth) - 1999
    cpgShoresSouth$feature <- "southShore"

    start(cpgShoresNorth) <- end(cpgShoresNorth)
    end(cpgShoresNorth) <- end(cpgShoresNorth) + 1999
    cpgShoresNorth$feature <- "northShore"

    message("Prepare CpG Shelves")
    cpgShelvesSouth <- cpgShoresSouth
    cpgShelvesNorth <- cpgShoresNorth
    end(cpgShelvesSouth) <- start(cpgShelvesSouth)
    start(cpgShelvesSouth) <- start(cpgShelvesSouth) - 1999
    cpgShelvesSouth$feature <- "southShelf"

    start(cpgShelvesNorth) <- end(cpgShelvesNorth)
    end(cpgShelvesNorth) <- end(cpgShelvesNorth) + 1999
    cpgShelvesNorth$feature <- "northShelf"

    cgiAnnotation <- rbind(as.data.frame(cpgIslands), as.data.frame(cpgShoresNorth), as.data.frame(cpgShoresSouth), as.data.frame(cpgShelvesNorth), as.data.frame(cpgShelvesSouth))
    cgiAnnotation <- cgiAnnotation[order(cgiAnnotation$seqnames, cgiAnnotation$start),]

    return(cgiAnnotation)
}
