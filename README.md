# Scripts might be used along with ChAMP

During runing ChAMP, some users required some extra features like paired calculation, sample matching .e.g. So I wrote some scripts here for users to use, in most case, users just need to source these script, then run them along with normal ChAMP pipeline.

**champ.PairedDMP.R**: To find paired Differential Methylation Probes in your Beta matrix (M matrix).

**champ.PairedDMR.R**: Using bumphunter algorithm to get paired DMR result.

**MatchSample.R**: To help user to find mismatched sample between pd csv file and IDAT folder.
