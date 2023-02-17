
library("data.table")
library("tidyverse")
library("glue")
library("magrittr")


manifest_dir <- "./MethylationEPIC v2.0 Files/EPIC-8v2-0_A1.csv"
assay_line <- which(substr(readLines(manifest_dir),1,7) == "[Assay]")
control_line <- which(substr(readLines(manifest_dir),1,10) == "[Controls]")

message("CpG to Probe matching...", appendLF = FALSE)
cpg <- fread(manifest_dir, skip=7, nrows =control_line-assay_line-2) %>%
  .[,c("IlmnID", "AddressA_ID", "AddressB_ID", "Infinium_Design_Type", "Color_Channel"), with=FALSE]

type_II <- cpg[Infinium_Design_Type == "II",]
type_I_Red <- cpg[Infinium_Design_Type == "I" & Color_Channel=="Red",]
type_I_Grn <- cpg[Infinium_Design_Type == "I" & Color_Channel=="Grn",]

M_name <- c(type_II$IlmnID, type_I_Red$IlmnID, type_I_Grn$IlmnID)

M_index <- c(glue("G-{type_II$AddressA_ID}"), 
             glue("R-{type_I_Red$AddressB_ID}"),
             glue("G-{type_I_Grn$AddressB_ID}"))

U_index <- c(glue("R-{type_II$AddressA_ID}"), 
             glue("R-{type_I_Red$AddressA_ID}"),
             glue("G-{type_I_Grn$AddressA_ID}"))

P_Channal <- c(rep("g+r",nrow(type_II)),
               rep("r",nrow(type_I_Red)),
               rep("g",nrow(type_I_Grn)))
message("Done")


message("Extract Probes...", appendLF = FALSE)
cpg_anno <- data.frame(CpG=M_name, M.index=M_index, U.index=U_index, Channel=P_Channal)
rownames(cpg_anno) <- cpg_anno$CpG

controls_anno <- fread(manifest_dir, skip=control_line + 1, header=FALSE) %>% .[,1:4] %>% as.data.frame
rownames(controls_anno) <- controls_anno$V1
controls_anno <- controls_anno[,2:4]
message("Done")


Anno <- list(Annotation=cpg_anno, ControlProbe=controls_anno)
save(Anno, file="AnnoEPICv2.rda")

