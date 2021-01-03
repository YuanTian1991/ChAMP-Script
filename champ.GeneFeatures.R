# A script to convert UCSC MySQL refGene into gene features for mapping.
# Author: Tian

library(RMySQL)
library(GenomicRanges)
library(stringr)

champ.GeneFeatures <- function(db = 'hg19', 
                               promoterRange=c(2000, 2000), 
                               features=c("Enhancer","Promoter", "TSS200_1500", "UTR5", "UTR3", "Exons", "Introns"))
{
    
    message("Fetch ",db , " refGene from UCSC mySQL database.")
    con_ucsc <- dbConnect(RMySQL::MySQL(), db = db, user = "genome", host = "genome-mysql.soe.ucsc.edu")
    refGene <- suppressWarnings(dbGetQuery(con_ucsc, stringr::str_interp("SELECT * FROM refGene")))
    dbDisconnect(con_ucsc)
    
    refGene <- refGene[order(refGene$chrom, refGene$txStart), ]
    refGene$rank <- 1:nrow(refGene)

    result <- list()
    
    if(any(c("Promoter", "Enhancer") %in% features)){
        message("Generate Promoter: Upstream ", promoterRange[1], " to Downstream ", promoterRange[2], " around TSS.")
        RefInfo <- refGene[,c("chrom", "txStart", "txEnd", "strand", "name", "name2", "rank")]
        colnames(RefInfo) <- c("seqnames", "start", "end", "strand", "id", "symbol", "rank")
        Ref.gr <- makeGRangesFromDataFrame(RefInfo, keep.extra.columns=TRUE)
        Promoter <- as.data.frame(promoters(Ref.gr,upstream=promoterRange[1], downstream=promoterRange[2]))
        Promoter$feature <- "Promoter"
        result$Promoter = Promoter
    }
    
    if("Enhancer" %in% features) {
        message("Generate Enhancer: Upstream 2000 to Promoter.")
        promoter.gr <- makeGRangesFromDataFrame(Promoter, keep.extra.columns=TRUE)
        Enhancer <- as.data.frame(promoters(promoter.gr, upstream=2000, downstream=0))
        Enhancer$feature <- "Enhancer"
        result$Enhancer = Enhancer
    }
    
    if("TSS200_1500" %in% features) {
        message("Generate TSS200: Upstream 200 of TSS.")
        RefInfo <- refGene[,c("chrom", "txStart", "txEnd", "strand", "name", "name2", "rank")]
        colnames(RefInfo) <- c("seqnames", "start", "end", "strand", "id", "symbol", "rank")
        Ref.gr <- makeGRangesFromDataFrame(RefInfo, keep.extra.columns=TRUE)
        TSS200 <- as.data.frame(promoters(Ref.gr,upstream=200, downstream=0))
        TSS200$feature <- "TSS200"
        result$TSS200 <- TSS200
        
        message("Generate TSS1500: Upstream 200 to 1500 of TSS.")
        TSS1500 <- as.data.frame(promoters(Ref.gr,upstream=1500, downstream=0))
        TSS1500$feature <- "TSS1500"
        TSS1500[TSS1500$strand == '+', "end"] <- TSS200[TSS200$strand == '+', "start"]
        TSS1500[TSS1500$strand == '-', "start"] <- TSS200[TSS200$strand == '-', "end"]
        TSS1500$width <- 1300
        result$TSS1500 <- TSS1500
    }
    
    if("UTR5" %in% features) {
        message("Generate 5'UTR: TSS to CDS Start Site.")
        UTR5 <- data.frame(seqnames=refGene$chrom,
                           start=refGene$txStart,
                           end=refGene$cdsStart,
                           width=refGene$cdsStart - refGene$txStart,
                           strand=refGene$strand,
                           id=refGene$name,
                           symbol=refGene$name2,
                           rank=refGene$rank,
                           feature="UTR5"
        )
        result$UTR5 <- UTR5
    }

    if("UTR3" %in% features) {
        message("Generate 3'UTR: CDS End Site to Gene End.")
        UTR3 <- data.frame(seqnames=refGene$chrom,
                           start=refGene$cdsEnd,
                           end=refGene$txEnd,
                           width=refGene$txEnd - refGene$cdsEnd,
                           strand=refGene$strand,
                           id=refGene$name,
                           symbol=refGene$name2,
                           rank=refGene$rank,
                           feature="UTR3"
        )
        result$UTR3 <- UTR3
    }

    if("Exons" %in% features)
    {
        message("Generate Exons: From each exonStarts to each exonEnds.")
        Exons <- suppressWarnings(do.call('rbind', apply(refGene, 1, function(x) {
                                                       names(x) <- colnames(refGene)
                                                       data.frame(seqnames=x["chrom"],
                                                                  start=as.numeric(strsplit(x["exonStarts"], split=",")[[1]]),
                                                                  end=as.numeric(strsplit(x["exonEnds"], split=",")[[1]]),
                                                                  width=as.numeric(strsplit(x["exonEnds"], split=",")[[1]]) - as.numeric(strsplit(x["exonStarts"], split=",")[[1]]),
                                                                  strand=x["strand"],
                                                                  id=x["name"],
                                                                  symbol=x["name2"],
                                                                  rank=as.numeric(x["rank"]),
                                                                  feature=paste0("Exon_", 1:x["exonCount"])
                                                       )})))
        result$Exons <- Exons
    }

    if("Introns" %in% features){
        message("Generate Introns: From each exonEnds to next exonStarts.")
        Introns <-  suppressWarnings(do.call('rbind', apply(refGene, 1, function(x) {
                                                          names(x) <- colnames(refGene)
                                                          exonEnds <- as.numeric(strsplit(x["exonEnds"], split=",")[[1]])
                                                          exonStarts <- as.numeric(strsplit(x["exonStarts"], split=",")[[1]])
                                                          exonCount <- as.numeric(x["exonCount"])
                                                          
                                                          if(exonCount <= 1) return(NULL)

                                                          data.frame(seqnames=x["chrom"],
                                                                     start=exonEnds[-length(exonEnds)],
                                                                     end=exonStarts[-1],
                                                                     width=exonStarts[-1] - exonEnds[-length(exonEnds)],
                                                                     strand=x["strand"],
                                                                     id=x["name"],
                                                                     symbol=x["name2"],
                                                                     rank=as.numeric(x["rank"]),
                                                                     feature=paste0("Intron_", 1:(exonCount - 1))
                                                          )})))
        result$Introns <- Introns
    }
    
    message("Merge all above Gene Features togather.")
    geneFeature <- do.call("rbind", result)

    if(length(result) > 0){
        message("Sorting...")
        geneFeature <- geneFeature[order(geneFeature$rank, geneFeature$start), ]
        geneFeature <- geneFeature[, setdiff(colnames(geneFeature), c("rank"))]
        rownames(geneFeature) <- c()
    }

    result <- lapply(result, function(x) {
                         tmp <- x[, setdiff(colnames(x), c("rank"))]
                         rownames(tmp) <- c()
                         return(tmp)
        })
    
    return(list(geneFeature=geneFeature, geneFeatureList=result))
}
