# A script to get TFregulomeR peaks I need.
# Author: Tian

library(TFregulomeR)

champ.getTFPeaks <- function(species, organ)
{
    
    message("Preparing TFregulomeR TF Peaks for ", organ, " from ", species, " database.")
    
    all_record <- suppressMessages(dataBrowser())
    index <- which(all_record$species == species & all_record$organ == organ)
    if(length(index) == 0) stop("Not even one TF found in TFregulomeR for your species/organ")
    
    TFList <- all_record[index, ]
    TFPeaks <- list()
    for(i in TFList$ID) {
        TFPeaks[[i]] <- suppressMessages(loadPeaks(id = i, includeMotifOnly = FALSE))
    }
    return(list(TFList=TFList, TFPeaks=TFPeaks))
}
