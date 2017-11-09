# lincs-pathprint
New pipeline for pathprinting LINCS data

## Scripts
1. **lincs_pathprint_download.R** : Script to download the Level 3 GSE from LINCS (GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx) and split them in a number of chunks (adjustable by number_of_chunks variable). Before saving to an equal number of matrices, the gene IDs are updated using the entrezUpdate.R script from the package GMAfunctions. 

2. **lincs_pathprint_runner.R** : Script that initiates a number of cores, each one utilizing a number of experiment expression matrices to run the *gctx2fingerprint.R* script.

3. **gctx2fingerprint.R** : Script that ... and produce the experiment-fingerprints.