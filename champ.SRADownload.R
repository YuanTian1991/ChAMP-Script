# A very good script I wrote to download SRA in minutes.

champ.SRADownload <- function(file, threads)
{
    sra <- read.csv(file, stringsAsFactors=FALSE)
    
    SRADownload <- function(x)
    {
        cmd = paste("parallel-fastq-dump --sra-id ", x ," --threads ",threads," --outdir ./Data/ --split-files --gzip")
        cat(cmd,"\n")
        system(cmd)
    }
    
    for(f in sra$Run) {
        SRADownload(f)
    }
    return()
}
