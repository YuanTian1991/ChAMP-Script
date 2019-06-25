##############################################################
# Author: Yuan Tian
# Description: This is a script to run paired analysis on DNA Methylation array, which could generate DMR result in accord with other ChAMP pipeline. What this script is doing is hacked into bumphunter algorithm, replaced the non-paired analysis coef with paried result.
##############################################################
# Usage: After prepared myNorm, myLoad .e.g, just source(champ.PairedDMR) Then use this function in the R session, as ChAMP.PairedDMR().
##############################################################

if(getRversion() >= "3.1.0") utils::globalVariables(c("myNorm","myLoad","probe.features.epic","probe.features"))

champ.PairedDMR <- function(beta = myNorm,
                            pair = NULL,
                            pheno = myLoad$pd$Sample_Group,
                            cutoff=NULL,
                            pickCutoff=TRUE,
                            B=250,
                            cores = 3,
                            maxGap = 300,
                            minProbes = 7,
                            bpspan=1000,
                            adjPvalDmr=0.05,
                            arraytype = "450K")
{
    message("[===========================]")
    message("[<< ChAMP.PairedDMR START >>]")
    message("-----------------------------")

    if(is.null(pair)){
        stop("Pair information must be prepared for Paired DMR Analysis.")
    }

    if(is.null(pheno) | length(unique(pheno))<=1)
    {
        stop("pheno parameter is invalid. Please check the input, pheno MUST contain at least two phenotypes.")
    }
	
    message("\n<< Checking Pairs >>")
    if(!all(table(pair)==2))
        stop("Valid Pairs for compare sampels are required. Odd numbers have been detected in your compared data. But in paired information, each patient's name should appear exactly twice.")
    if(!all(table(pair,pheno)==1))
        stop("Pheno and Pairs must corrsponding to each other. The match between your paired information and pheno is not correct.")

    message("\n<<Pairs Checking Success, Running PairedDMR now...>>")

    if(arraytype=="EPIC"){
        RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b2.hg19"))
    }else{
        RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylation450k",annotation = "ilmn12.hg19"))
    }
    message("[Paired DMR] Calculation clusters.")
    probe.features <- getAnnotation(RSobject)
    cpg.idx <- intersect(rownames(beta),rownames(probe.features))
    Anno <- probe.features[cpg.idx,]
    Anno <- Anno[order(Anno$chr,Anno$pos),]
    cpg.idx <- rownames(Anno)
    cl <- clusterMaker(Anno$chr,Anno$pos,maxGap=maxGap)
    names(cl) <- cpg.idx
    bumphunter.idx <- cpg.idx[which(cl %in% names(which(table(cl) > minProbes)))]
    beta <- beta[bumphunter.idx,]
    chr <- Anno[bumphunter.idx,]$chr
    pos <- Anno[bumphunter.idx,]$pos
    cluster <- cl[bumphunter.idx]

    beta <- beta[bumphunter.idx,]
    beta<- replace(beta,which(beta <= 0.001),0.001)
    beta <- replace(beta,which(beta >= 0.999),0.999)
    Y <- log((beta/(1-beta)),2)

    # Get all t statistic 
    message("[Paired DMR] Computing t statistic.")
    rawBeta <- champ.PairedDMP(Y, pair, pheno,adjPVal = 1,arraytype=arraytype)
    rawBeta <- rawBeta[rownames(Y),"t"]
    names(rawBeta) <- rownames(Y)

    mySmooth <- function(y,x,mycluster,bpspan=bpspan,cores=cores)
    {
        tmpPos <- split(x,mycluster)
        tmpLength <- unlist(lapply(tmpPos,function(x) length(x)))
        names(tmpLength) <- names(tmpPos)
        tmpSpan <- unlist(lapply(tmpPos,function(x) bpspan/median(diff(x))/length(x)))
        names(tmpSpan) <- names(tmpPos)
        tmpSpan[which(tmpSpan>1)] <- 1

        if(class(y)=="numeric")
        {
            message("Only one vector.")
            tmpCpG <- split(y,mycluster)
            tmpSmooth <- unlist(sapply(names(tmpCpG),function(m) limma::loessFit(tmpCpG[[m]],tmpPos[[m]],span = tmpSpan[m])$fitted))
        }
        else
        {
            smallfunction <- function(tmp_y,tmp_cluster,tmp_pos,tmp_span)
            {
                tmpCpG <- split(tmp_y,tmp_cluster)
                smallsmooth <- unlist(sapply(names(tmpCpG),function(m) limma::loessFit(tmpCpG[[m]],tmp_pos[[m]],span = tmp_span[m])$fitted))
                return(smallsmooth)
            }
            message("Parallel on Matrix.")
            parallelclusters <- makeCluster(cores)
            registerDoParallel(parallelclusters)
            getDoParWorkers()
            tmpSmooth <- foreach(i = 1:ncol(y), .combine = cbind) %dopar% smallfunction(y[,i],mycluster,tmpPos,tmpSpan)
        }
        return(tmpSmooth)
    }

        # Smooth 
    message("[Paired DMR] Smoothing t statistic.")
    rawBeta_fitted <- mySmooth(y = rawBeta, x = pos, mycluster = cluster,bpspan=bpspan)

    # bootstrap or permuation 
    message("[Paired DMR] Applying Permutations.")
    probe_permutation <- replicate(B,sample(rawBeta))

    # Computing Marginals
    greaterOrEqual <- function(x,y) {
        precision <- sqrt(.Machine$double.eps)
        (x >= y) | (abs(x-y) <= precision)
    }

    # Smoothing on Permutation Beta, This process if very slow...
    message("[Paired DMR] Smoothing t statistic for permutation.")
    permBeta <- mySmooth(probe_permutation,pos,cluster,bpspan=bpspan,cores=cores)

    # Set cutoff.
    message("[Paired DMR] Setting Cutoff.")
    if(pickCutoff == FALSE & is.null(cutoff))
    {
        message("User setted cutoff value here.")
        cutoff <- cutoff
    }else{
        message("User did not set cutoff, 0.99 quantile of all permutation smooth value will be used as cutoff.")
        cutoff <- quantile(abs(permBeta), 0.99, na.rm = TRUE)
    }
    message(sprintf("[Paired DMR] cutoff: %s",round(cutoff, 3)))

    # Start to find bumps
    message("[Paired DMR] Finding regions.")
    tab <- regionFinder(x = rawBeta_fitted, chr = chr, pos = pos, cluster = cluster,cutoff = cutoff, verbose = FALSE)
    message(sprintf("[Paired DMR] Found %s bumps.",nrow(tab)))

    message("[Paired DMR] Computing regions for each Permutation.")
    chunksize <- ceiling(B/cores)
    subMat <- NULL
    nulltabs <- foreach(subMat = iter(permBeta, by = "col", chunksize = chunksize),.combine = "c", .packages = "bumphunter") %dorng% {
    apply(subMat, 2, regionFinder, chr = chr, pos = pos,cluster = cluster, cutoff = cutoff,verbose = FALSE)}

        # Calulcate p value and FWER.
    message("[Paired DMR] Estimating p-values and FWER.")
    L <- V <- A <- as.list(rep(0, B))
    for (i in 1:B) {
        nulltab <- nulltabs[[i]]
        if (nrow(nulltab) > 0) {
            L[[i]] <- nulltab$L
            V[[i]] <- nulltab$value
            A[[i]] <- nulltab$area
        }
    }
    mytots <- apply(tab[,c("L","value")],1,function(i) sapply(c(1:B),function(x) sum(greaterOrEqual(L[[x]],i[1]) & greaterOrEqual(abs(V[[x]]),abs(i[2])))))
    mytots2 <- t(sapply(tab$area,function(i) sapply(c(1:B),function(x) sum(greaterOrEqual(A[[x]],i)))))
    rate1 <- apply(mytots,2,function(x) mean(x>0))
    pvalues1 <- apply(mytots,2,function(x) sum(x))/sum(sapply(nulltabs, nrow))
    rate2 <- rowMeans(mytots2 > 0)
    pvalues2 <- rowSums(mytots2)/sum(sapply(nulltabs, nrow))

    # Finally
    tab$p.value <- pvalues1
    tab$fwer <- rate1
    tab$p.valueArea <- pvalues2
    tab$fwerArea <- rate2
    tab <- tab[order(tab$fwer, -tab$area), ]

    message("<< Calculate Paired DMR success. >>")
    PairedDMR <- tab[which(tab$p.valueArea <= adjPvalDmr),]
    message("ChAMP.PairedDMR detected ",nrow(PairedDMR)," Paired DMRs with P value <= ",adjPvalDmr,".")

    if(nrow(PairedDMR) == 0) stop("No Paired DMR detected.")

    rownames(PairedDMR) <- paste("PairedDMR",1:nrow(PairedDMR),sep="_")
    PairedDMR <- data.frame(PairedDMR[,1:3],width=PairedDMR[,3]-PairedDMR[,2],strand="*",PairedDMR[,4:14])
    colnames(PairedDMR)[1:3] <- c("seqnames","start","end")

    message("[<<< ChAMP.PairedDMR END >>>]")
    message("[===========================]")
    message("[You may want to process PairedDMR.GUI() or champ.GSEA() next.]\n")
    return(PairedDMR)
}
