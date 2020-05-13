# locFDR outlier detection
library("locfdr")
library("plyr")
library("ggplot2")

# meta <- annot data <- mval

locFDR <- function(data, meta) {
    
    pca <- as.data.frame(unclass(princomp(data)$loadings))
    zstat_comp1 <- (pca$Comp.1 - mean(pca$Comp.1))/sd(pca$Comp.1)
    w <- locfdr(as.matrix(zstat_comp1))
    fdr <- as.data.frame(w$fdr)
    rownames(fdr) <- rownames(pca)
    outmeta <- meta[which(fdr <= 0.1), ]
    outmeta$pca <- "FDR<=0.1"
    inmeta <- meta[which(!(fdr <= 0.1)), ]
    inmeta$pca <- "FDR>=0.1"
    outliers.locfdr <- as.character(outmeta$Sample_Name)
    meta_out <- rbind(outmeta, inmeta)
    meta_out <- meta_out[order(match(colnames(data), meta_out$Sample_Name)), ]

    print(outliers.locfdr)

    pca_plot <- pca[, 1:5]
    pca_plot$Sample_Name <- rownames(pca_plot)
    
    pca_plot <- join(pca_plot, meta_out, by = "Sample_Name")
    p <- ggplot(pca_plot, aes(Comp.1, Comp.2, color = pca)) + geom_point(shape = 19, size = 4, 
        alpha = 0.6) + theme_bw() + xlab("PC1") + ylab("PC2") + scale_color_manual(name = "Outlier\nThreshold", 
        values = c("#005C66", "#733556"))
    
    return(list(w=w,p=p))
}
# Remove samples that were detected as outliers from data and repeat until no
# more outliers are detected
