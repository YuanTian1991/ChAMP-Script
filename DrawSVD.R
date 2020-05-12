# A script to better show SVD result from ChAMP
# Author: Yuan Tian

library("pheatmap")
library("RColorBrewer")

DrawSVDPlot <- function(mySVD) {
    
    mySVD_log <- log(mySVD, 10)
    breaks.v <- c(-10000, -10, -5, -2, log10(0.05), 0)
    tags <- c(4, 3, 2, 1, 0)
    
    group_tags <- cut(as.numeric(mySVD_log), breaks = breaks.v, include.lowest = FALSE, 
        right = FALSE, labels = tags)
    
    SVD.df <- matrix(as.numeric(as.character(group_tags)), nrow = nrow(mySVD))
    colnames(SVD.df) <- colnames(mySVD)
    rownames(SVD.df) <- paste("pc", 1:nrow(mySVD), sep = "-")
    
    myPalette <- c("snow", brewer.pal(4, "Reds"))
    myLegendLabel <- c("p < 1*10^(-10)", "p < 1*10^(-5)", "p < 0.01", "p < 0.05", 
        "p > 0.05")
    
    pheatmap(t(SVD.df[nrow(SVD.df):1, ncol(SVD.df):1]), cluster_cols = F, cluster_rows = F, 
        color = myPalette, legend_breaks = c(4, 3, 2, 1, 0), legend_labels = myLegendLabel)
}
