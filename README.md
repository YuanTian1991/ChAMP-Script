# Scripts might be used along with ChAMP

During runing ChAMP, some users required some extra features like paired calculation, sample matching .e.g. So I wrote some scripts here for users to use, in most case, users just need to source these script, then run them along with normal ChAMP pipeline.

**champ.PairedDMP.R**: To find paired Differential Methylation Probes in your Beta matrix (M matrix).

**champ.PairedDMR.R**: Using bumphunter algorithm to get paired DMR result.

**MatchSample.R**: To help user to find mismatched sample between pd csv file and IDAT folder.

**champ.Anno.R**: Convert Illumina origin CSV file to ChAMP Annotation.

**champ.unfoldDMR.R**: Unfold generate DMR result.

**champ.SlimCombat.R**: A quick and short function to do Combat, need to be polished.

**champ.GeneFeatures.R**: A function to generate gene features like promoter, 1 exon, intron, 5UTR, 3UTR from raw UCSC RefGene.

**champ.CpGIslands.R**: A function to generate CpG Islands, Shore and Shelf from UCSC table cpgIslandsExt.

**champ.getTFPeaks.R**: A function to get Transcript Factor peaks from TFregulomeR package. The result should be feed to champ.TFEA.R function.

**champ.TFEA.R**: A function to do Transcript Factor Binding Site Enrichment Analysis (TFEA) for a list of Region of Interest (ROI), like genes, promoters, differential methylated regions, differential hydroxty-methylated regions .etc.

**champ.PeakEnrich.R**: A function to do peak set enrichment, currently only support MAnorm2 preprocessed ChIP-seq, actually only validated on MeDIP-seq and hMeDIP-seq.

**champ.Overlap.R**: A function to do overlap collapsing between two genome segments(peaks). Overlap x and y, and collapse y records to each matched x record. Functions for collapse can be self-defined.
