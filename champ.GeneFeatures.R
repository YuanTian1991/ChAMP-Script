# A script to convert UCSC MySQL refGene into gene features for mapping.  Author:
# Tian

library(RMySQL)
library(GenomicRanges)
library(stringr)

champ.GeneFeatures <- function(db = "hg19", promoterRange = c(2000, 2000)) {
    
    message("Fetch ", db, " refGene from UCSC mySQL database")
    con_ucsc <- dbConnect(RMySQL::MySQL(), db = db, user = "genome", host = "genome-mysql.soe.ucsc.edu")
    refGene <- suppressWarnings(dbGetQuery(con_ucsc, stringr::str_interp("SELECT * FROM refGene")))
    dbDisconnect(con_ucsc)
    
    geneRank <- rank(refGene$chrom, refGene$txStart)
    
    message("Generate Promoter: Upstream ", promoterRange[1], " to Downstream ", 
        promoterRange[2], " around TSS.")
    RefInfo <- refGene[, c("chrom", "txStart", "txEnd", "strand", "name", "name2")]
    colnames(RefInfo) <- c("seqnames", "start", "end", "strand", "id", "symbol")
    Ref.gr <- makeGRangesFromDataFrame(RefInfo, keep.extra.columns = TRUE)
    Promoter <- as.data.frame(promoters(Ref.gr, upstream = promoterRange[1], downstream = promoterRange[2]))
    Promoter$feature <- "promoter"
    
    
    message("Generate 5'UTR: TSS to CDS Start Site")
    UTR5 <- data.frame(seqnames = refGene$chrom, start = refGene$txStart, end = refGene$cdsStart, 
        width = refGene$cdsStart - refGene$txStart, strand = refGene$strand, id = refGene$name, 
        symbol = refGene$name2, feature = "UTR5")
    
    
    message("Generate 3'UTR: CDS End Site to Gene End")
    UTR3 <- data.frame(seqnames = refGene$chrom, start = refGene$cdsEnd, end = refGene$txEnd, 
        width = refGene$txEnd - refGene$cdsEnd, strand = refGene$strand, id = refGene$name, 
        symbol = refGene$name2, feature = "UTR3")
    
    message("Generate Exons: From each exonStarts to each exonEnds")
    Exons <- suppressWarnings(do.call("rbind", apply(refGene, 1, function(x) {
        names(x) <- colnames(refGene)
        data.frame(seqnames = x["chrom"], start = as.numeric(strsplit(x["exonStarts"], 
            split = ",")[[1]]), end = as.numeric(strsplit(x["exonEnds"], split = ",")[[1]]), 
            width = as.numeric(strsplit(x["exonEnds"], split = ",")[[1]]) - as.numeric(strsplit(x["exonStarts"], 
                split = ",")[[1]]), strand = x["strand"], id = x["name"], symbol = x["name2"], 
            feature = paste0("Exon_", 1:x["exonCount"]))
    })))
    
    message("Generate Intros: From each exonEnds to next exonStarts")
    Introns <- suppressWarnings(do.call("rbind", apply(refGene, 1, function(x) {
        
        names(x) <- colnames(refGene)
        exonEnds <- as.numeric(strsplit(x["exonEnds"], split = ",")[[1]])
        exonStarts <- as.numeric(strsplit(x["exonStarts"], split = ",")[[1]])
        exonCount <- as.numeric(x["exonCount"])
        
        if (exonCount <= 1) 
            return(NULL)
        
        data.frame(seqnames = x["chrom"], start = exonEnds[-length(exonEnds)], end = exonStarts[-1], 
            width = exonStarts[-1] - exonEnds[-length(exonEnds)], strand = x["strand"], 
            id = x["name"], symbol = x["name2"], feature = paste0("Intron_", 1:(exonCount - 
                1)))
    })))
    
    message("Merge all above Gene Features togather.")
    geneFeature <- rbind(Promoter, UTR5, Exons, Introns, UTR3)
    message("Sorting...")
    geneFeature <- geneFeature[order(geneFeature$id, geneFeature$seqnames, geneFeature$start), 
        ]
    
    return(geneFeature)
}
