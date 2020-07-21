# This is a script to unfold DMR result from ChAMP, which means it would return
# DMPs related in it and corresponding genes.  Author: Yuan Tian

library("ChAMP")

champ.unfoldDMR <- function(DMR = myDMR, beta = myNorm, pheno = myLoad$pd$Sample_Group, 
    compare.group = NULL, arraytype = "450K") {

    tmpbeta <- beta
    tmppheno <- pheno

    message("Calculating DMP")
    DMP <- champ.DMP(beta = tmpbeta, pheno = tmppheno, adjPVal = 1, adjust.method = "BH", 
        compare.group = compare.group, arraytype = arraytype)
    DMP <- DMP[[1]]
    
    
    if (arraytype == "EPIC") data(probe.features.epic) else data(probe.features)
    probe.features <- probe.features[rownames(beta), ]
    
    
    message("  Generating Annotation File")
    DMR[[1]]$seqnames <- as.factor(substr(DMR[[1]]$seqnames, 4, 100))
    
    index <- apply(DMR[[1]], 1, function(x) which(probe.features$CHR == x[1] & probe.features$MAPINFO >= 
        as.numeric(x[2]) & probe.features$MAPINFO <= as.numeric(x[3])))
    
    Anno <- data.frame(DMRindex = unname(unlist(sapply(names(index), function(x) rep(x, 
        length(index[[x]]))))), probe.features[do.call(c, index), 1:8])
    message("  Generating Annotation File Success")

    message("  Combind DMP with DMR index.")
    Anno <- data.frame(Anno, DMP[rownames(Anno),1:10])
    
    Anno <- split(Anno, Anno$DMRindex)

    message("\n--------------------------------")
    message("------- DMR Unfolded -----------")
    message("Assume your result is named as myUnfoldedDMR:")
    message("  \n1. Returned Result is a List, for example use myUnfoldedDMR[[1]] to see each result.")
    message("  \n2. If you want to join all list into a large data frame, try: result <- do.call(rbind,myUnfoldedDMR)")
    message("  \n3. If you only want DMP, instead of all CpGs, try: result <- lapply(myUnfoldedDMR,function(x) x[x$adj.P.Val <= 0.05,])")
    message("  \n4. If you want DMR-related genes, try: genes <- unique(do.call(rbind,myUnfoldedDMR)[,'gene'])")

    return(Anno)
}
