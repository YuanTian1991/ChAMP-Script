##############################################################
# Author: Yuan Tian
# Description: This is a script to check the mismatch between pd csv file and IDAT file in user's folder. This script would based on csv file, only find samples shows on the csv, but can not be find in the folder. Users may run this script IN THE IDAT FILE FOLDER, and find out the missing samples.
#
# In most case, there are couple reasons your csv file can not match the IDAT file: 1. Excel automatically changed Sentrix_ID format; 2. You have duplicated copied IDAT in this folder, means you may have two files with the same name in different folder; 3. You have missing samples on your csv files.
##############################################################
# Usage: Open R session, setwd() to your IDAT folder. Then copy and paste this script in the R session.
##############################################################


directory = getwd()
csvfile <- list.files(directory,recursive=TRUE,pattern="csv$",full.names=TRUE)

skipline <- which(substr(readLines(csvfile),1,6) == "[Data]")
suppressWarnings(pd <- read.csv(csvfile,skip=skipline,stringsAsFactor=FALSE,        header=TRUE))

colnames(pd)[which(colnames(pd)=="Sentrix_Position")] <- "Array"
colnames(pd)[which(colnames(pd)=="Sentrix_ID")] <- "Slide"

GrnPath <- unlist(sapply(paste(pd$Slide,pd$Array,"Grn.idat",sep="_"), function(x)       grep(x,list.files(directory,recursive=T,full.names=TRUE), value = TRUE)))
RedPath <- unlist(sapply(paste(pd$Slide,pd$Array,"Red.idat",sep="_"), function(x)       grep(x,list.files(directory,recursive=T,full.names=TRUE), value = TRUE)))

csvGrnName <- paste(pd$Slide,pd$Array,"Grn.idat",sep="_")
csvRedName <- paste(pd$Slide,pd$Array,"Red.idat",sep="_")

GrnFiles <- names(GrnPath)
RedFiles <- names(RedPath)

     message("Below are Green IDAT in csv file, but not in your folder.")
     message("-----------------------------------------------------")
     print(setdiff(csvGrnName,GrnFiles))
     message("-----------------------------------------------------\n\n")

     message("Below are Red IDAT in csv file, but not in your folder.")
     message("-----------------------------------------------------")
     print(setdiff(csvRedName,RedFiles))
     message("-----------------------------------------------------")
