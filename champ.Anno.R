# A R function to generate ChAMP annotation from original Illumina CSV file.
# Author: Tian

# manifest_directory <- "../../ChAMP-Script/Illumina_Methylation_Manifests/MethylationEPIC_v-1-0_B4.csv"

champ.Anno <- function(manifest_directory)
{
    message("Extract Annotation for ChAMP from original Illumina Human or Mouse Methylation Array")
    message("Read CSV manifest into R.")
    skipline <- which(substr(readLines(manifest_directory),1,7) == "[Assay]")
    Manifest <- read.csv(manifest_directory, head=T, sep=",", skip=skipline, as.is=T)
    
    Anno <- Manifest[,c("IlmnID", "AddressA_ID", "AddressB_ID", "Infinium_Design_Type", "Color_Channel")]
    message("Total Rows in this manifest: ", nrow(Anno))
    
    message("There are totally 3 types of Probes, Type-II, Type-I-Green and Type-I-Red")
    Type_II <- Anno[which(Anno$Infinium_Design_Type=="II"),]
    Type_I.Red <- Anno[which(Anno$Infinium_Design_Type == "I" & Anno$Color_Channel=="Red"),]
    Type_I.Grn <- Anno[which(Anno$Infinium_Design_Type == "I" & Anno$Color_Channel=="Grn"),]
    
    CpG.ID <- c(as.character(Type_II$IlmnID),as.character(Type_I.Red$IlmnID),as.character(Type_I.Grn$IlmnID))
    
    message(paste0("\033[0;43m", ' Key Step ' ,"\033[0m", " Get mapping index between Manifest and Illuminano pacakge output to generate Methylated Signal"))
    M.index <- c(paste("G",Type_II$AddressA_ID, sep="-"),
                 paste("R",Type_I.Red$AddressB_ID, sep="-"),
                 paste("G",Type_I.Grn$AddressB_ID, sep="-"))
    
    message(paste0("\033[0;43m", ' Key Step ' ,"\033[0m", " Get mapping index between Manifest and Illuminano pacakge output to generate UnMethylated Signal"))
    U.index <- c(paste("R",Type_II$AddressA_ID,sep="-"),
                 paste("R",Type_I.Red$AddressA_ID,sep="-"),
                 paste("G",Type_I.Grn$AddressA_ID,sep="-"))
    
    
    message("generate Colour Channel distribution")
    P.Channal <- c(rep("g+r",nrow(Type_II)),rep("r",nrow(Type_I.Red)),rep("g",nrow(Type_I.Grn)))
    
    message("Merge them into one Data.Frame for future usage.")
    ProbeInfo <- data.frame(CpG=CpG.ID, M.index=M.index, U.index=U.index, Channel=P.Channal)
    rownames(ProbeInfo) <- ProbeInfo$CpG
    
    message("Generate Control Probes.")
    message("Reading Control Probes from CSV.")
    skipline <- which(substr(readLines(manifest_directory),1,10) == "[Controls]")
    ControlProbe <- read.csv(manifest_directory, head=FALSE, sep=",", skip=skipline, as.is=TRUE)
    
    message("Generate CpG Annotation.")
    message("Read CSV manifest into R.")
    skipline <- which(substr(readLines(manifest_directory),1,7) == "[Assay]")
    Manifest <- read.csv(manifest_directory, head=T, sep=",", skip=skipline, as.is=T)
    annoCol <- c("IlmnID",
                 "CHR",
                 "MAPINFO",
                 "Strand",
                 "Infinium_Design_Type",
                 "UCSC_RefGene_Name",
                 "UCSC_RefGene_Group",
                 "Relation_to_UCSC_CpG_Island",
                 "DHS",
                 "Enhancer",
                 "Phantom",
                 "Probe_SNPs",
                 "Probe_SNPs_10",
                 "SNP_ID",
                 "SNP_DISTANCE",
                 "SNP_MinorAlleleFrequency",
                 "DNase_Hypersensitivity_NAME",
                 "OpenChromatin_NAME",
                 "TFBS_NAME",
                 "TFBS_Evidence_Count",
                 "CHR_hg38",
                 "Start_hg38",
                 "End_hg38",
                 "Strand_hg38")
    
    CpGAnno <- Manifest[, intersect(annoCol,colnames(Manifest))]
    rownames(CpGAnno) <- CpGAnno$IlmnID
    CpGAnno <- CpGAnno[rownames(ProbeInfo), ]
    
    return(list(ProbeInfo=ProbeInfo, CpGAnno=CpGAnno, ControlProbe=ControlProbe))
}
