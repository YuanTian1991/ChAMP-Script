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
        Promoter$InternalRank <- 2
        Promoter$InternalRank[Promoter$strand == '-'] <- 10004
        result$Promoter = Promoter
    }
    
    if("Enhancer" %in% features) {
        message("Generate Enhancer: Upstream 2000 to Promoter.")
        promoter.gr <- makeGRangesFromDataFrame(Promoter, keep.extra.columns=TRUE)
        Enhancer <- as.data.frame(promoters(promoter.gr, upstream=2000, downstream=0))
        Enhancer$feature <- "Enhancer"
        Enhancer$InternalRank <- 1
        Enhancer$InternalRank[Enhancer$strand == '-'] <- 10005
        result$Enhancer = Enhancer
    }
    
    if("TSS200_1500" %in% features) {
        message("Generate TSS200: Upstream 200 of TSS.")
        RefInfo <- refGene[,c("chrom", "txStart", "txEnd", "strand", "name", "name2", "rank")]
        colnames(RefInfo) <- c("seqnames", "start", "end", "strand", "id", "symbol", "rank")
        Ref.gr <- makeGRangesFromDataFrame(RefInfo, keep.extra.columns=TRUE)
        TSS200 <- as.data.frame(promoters(Ref.gr,upstream=200, downstream=0))
        TSS200$feature <- "TSS200"
        TSS200$InternalRank <- 4
        TSS200$InternalRank[TSS200$strand == '-'] <- 10002
        result$TSS200 <- TSS200
        
        message("Generate TSS1500: Upstream 200 to 1500 of TSS.")
        TSS1500 <- as.data.frame(promoters(Ref.gr,upstream=1500, downstream=0))
        TSS1500$feature <- "TSS1500"
        TSS1500[TSS1500$strand == '+', "end"] <- TSS200[TSS200$strand == '+', "start"]
        TSS1500[TSS1500$strand == '-', "start"] <- TSS200[TSS200$strand == '-', "end"]
        TSS1500$width <- 1300
        TSS1500$InternalRank <- 3
        TSS1500$InternalRank[TSS1500$strand == '-'] <- 10003
        result$TSS1500 <- TSS1500
    }
    
    if("UTR5" %in% features) {
        message("Generate 5'UTR: TSS to CDS Start Site.")

        strandPositive <- refGene[refGene$strand == "+",]

        UTR5Positive <- data.frame(seqnames=strandPositive$chrom,
                           start=strandPositive$txStart,
                           end=strandPositive$cdsStart,
                           width=strandPositive$cdsStart - strandPositive$txStart,
                           strand=strandPositive$strand,
                           id=strandPositive$name,
                           symbol=strandPositive$name2,
                           rank=strandPositive$rank,
                           feature="UTR5",
                           InternalRank=5
        )

        strandNegative <- refGene[refGene$strand == "-",]

        UTR5Negative <- data.frame(seqnames=strandNegative$chrom,
                           start=strandNegative$cdsEnd,
                           end=strandNegative$txEnd,
                           width=strandNegative$txEnd - strandNegative$cdsEnd,
                           strand=strandNegative$strand,
                           id=strandNegative$name,
                           symbol=strandNegative$name2,
                           rank=strandNegative$rank,
                           feature="UTR5",
                           InternalRank=10001
        )
        UTR5 <- rbind(UTR5Positive, UTR5Negative)
        UTR5 <- UTR5[order(UTR5$rank), ]
        result$UTR5 <- UTR5
    }

    if("UTR3" %in% features) {
        message("Generate 3'UTR: CDS End Site to Gene End.")

        strandPositive <- refGene[refGene$strand == "+",]

        UTR3Positive <- data.frame(seqnames=strandPositive$chrom,
                           start=strandPositive$cdsEnd,
                           end=strandPositive$txEnd,
                           width=strandPositive$txEnd - strandPositive$cdsEnd,
                           strand=strandPositive$strand,
                           id=strandPositive$name,
                           symbol=strandPositive$name2,
                           rank=strandPositive$rank,
                           feature="UTR3",
                           InternalRank=9999
        )

        strandNegative <- refGene[refGene$strand == "-",]

        UTR3Negative <- data.frame(seqnames=strandNegative$chrom,
                           start=strandNegative$txStart,
                           end=strandNegative$cdsStart,
                           width=strandNegative$cdsStart - strandNegative$txStart,
                           strand=strandNegative$strand,
                           id=strandNegative$name,
                           symbol=strandNegative$name2,
                           rank=strandNegative$rank,
                           feature="UTR3",
                           InternalRank=9
        )
        UTR3 <- rbind(UTR3Positive, UTR3Negative)
        UTR3 <- UTR3[order(UTR3$rank), ]
        result$UTR3 <- UTR3
    }

    if("Exons" %in% features)
    {
        message("Generate Exons: From each exonStarts to each exonEnds.")
        Exons <- suppressWarnings(do.call('rbind', apply(refGene, 1, function(x) {

                                                       names(x) <- colnames(refGene)

                                                       exonStarts <- as.numeric(strsplit(x["exonStarts"], split=",")[[1]])
                                                       exonEnds <- as.numeric(strsplit(x["exonEnds"], split=",")[[1]])

                                                       if(x["strand"] == "+") {
                                                           exonFeatures <- paste0("Exon_", 1:x["exonCount"])
                                                       } else {
                                                           exonFeatures <- paste0("Exon_", x["exonCount"]:1)
                                                       }

                                                       data.frame(seqnames=x["chrom"],
                                                                  start=exonStarts,
                                                                  end=exonEnds,
                                                                  width=exonEnds - exonStarts,
                                                                  strand=x["strand"],
                                                                  id=x["name"],
                                                                  symbol=x["name2"],
                                                                  rank=as.numeric(x["rank"]),
                                                                  feature=exonFeatures,
                                                                  InternalRank=100
                                                       )})))
        result$Exons <- Exons
    }

    if("Introns" %in% features){
        message("Generate Introns: From each exonEnds to next exonStarts, and gaps between cds and exons.")
        Introns <-  suppressWarnings(do.call('rbind', apply(refGene, 1, function(x) {
                                                          names(x) <- colnames(refGene)

                                                          exonEnds <- as.numeric(strsplit(x["exonEnds"], split=",")[[1]])
                                                          exonStarts <- as.numeric(strsplit(x["exonStarts"], split=",")[[1]])

                                                          IntronStarts <- exonEnds[-length(exonEnds)]
                                                          IntronEnds <- exonStarts[-1]

                                                          cdsStart <- as.numeric(x['cdsStart'])
                                                          cdsEnd <- as.numeric(x['cdsEnd'])

                                                          # Judge if there is a gap between cdsStart and first exon
                                                          if(cdsStart < exonStarts[1]) {
                                                              IntronStarts <- c(cdsStart, IntronStarts)
                                                              IntronEnds <- c(exonStarts[1], IntronEnds)
                                                          }

                                                          # Judge if there is a gap between last exon and cdsEnd
                                                          if(cdsEnd > tail(exonEnds,1)) {
                                                              IntronStarts <- c(IntronStarts,  tail(exonEnds,1))
                                                              IntronEnds <- c(IntronEnds, cdsEnd)
                                                          }

                                                          if(length(IntronStarts) <= 1) return(NULL)
                                                          
                                                          if(x["strand"] == "+") {
                                                              IntronFeatures <- paste0("Intron_", 1:length(IntronStarts))
                                                          } else {
                                                              IntronFeatures <- paste0("Intron_", length(IntronStarts): 1)
                                                          }
                                                          

                                                          data.frame(seqnames=x["chrom"],
                                                                     start=IntronStarts,
                                                                     end=IntronEnds,
                                                                     width=IntronEnds - IntronStarts,
                                                                     strand=x["strand"],
                                                                     id=x["name"],
                                                                     symbol=x["name2"],
                                                                     rank=as.numeric(x["rank"]),
                                                                     feature=IntronFeatures,
                                                                     InternalRank=100
                                                          )})))
        result$Introns <- Introns
    }
    
    message("Merge all above Gene Features togather.")
    geneFeature <- do.call("rbind", result)

    if(length(result) > 0){
        message("Sorting...")
        geneFeature <- geneFeature[order(geneFeature$rank, geneFeature$InternalRank, geneFeature$start), ]
        geneFeature <- geneFeature[, setdiff(colnames(geneFeature), c("rank", "InternalRank"))]
        rownames(geneFeature) <- c()
    }

    result <- lapply(result, function(x) {
                         tmp <- x[, setdiff(colnames(x), c("rank", "InternalRank"))]
                         rownames(tmp) <- c()
                         return(tmp)
        })
    
    return(list(geneFeature=geneFeature, geneFeatureList=result))
}
