
manifest_dir <- "~/Desktop/MouseMethylation-12v1-0_A2.csv"

skipline <- which(substr(readLines(manifest_dir),1,7) == "[Assay]")
A <- read.csv(manifest_dir, head=T, sep=",", stringsAsFactors=F, skip=skipline)
AnnoLoader450K <- A[,c("IlmnID", "AddressA_ID", "AddressB_ID", "Infinium_Design_Type", "Color_Channel")]

Type_II <- AnnoLoader450K[which(AnnoLoader450K$Infinium_Design_Type==2),]
Type_I.Red <- AnnoLoader450K[which(AnnoLoader450K$Infinium_Design_Type == 1 & AnnoLoader450K$Color_Channel=="Red"),]
Type_I.Grn <- AnnoLoader450K[which(AnnoLoader450K$Infinium_Design_Type == 1 & AnnoLoader450K$Color_Channel=="Grn"),]

M.name <- c(as.character(Type_II[,1]),as.character(Type_I.Red[,1]),as.character(Type_I.Grn[,1]))
M.index <- c(paste("G",Type_II[,2],sep="-"),paste("R",Type_I.Red[,3],sep="-"),paste("G",Type_I.Grn[,3],sep="-"))
U.index <- c(paste("R",Type_II[,2],sep="-"),paste("R",Type_I.Red[,2],sep="-"),paste("G",Type_I.Grn[,2],sep="-"))
P.Channal <- c(rep("g+r",nrow(Type_II)),rep("r",nrow(Type_I.Red)),rep("g",nrow(Type_I.Grn)))

Annotation <- data.frame(CpG=M.name,M.index=M.index,U.index=U.index,Channel=P.Channal)

rownames(Annotation) <- Annotation$CpG

Annotation$CpG <- as.character(Annotation$CpG)
Annotation$M.index <- as.character(Annotation$M.index)
Annotation$U.index <- as.character(Annotation$U.index)
Annotation$Channel <- as.character(Annotation$Channel)

