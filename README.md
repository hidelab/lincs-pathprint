# lincs-pathprint
New pipeline for pathprinting LINCS data

## Scripts
1. **gctx_to_SCE.R** : Script to download the Level 3 GSE from LINCS (GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx) and split them in a number of chunks (adjustable by number_of_chunks variable). Before saving to an equal number of matrices, the gene IDs are updated using the entrezUpdate.R script from the package GMAfunctions. Finally, the script initiates a number of cores each one handling running the gctx2fingerprint.R script to produce the experiment fingerprints.

2. **gctx2fingerprint.R** : Script that based on a gctx partial matrix produces the expreriment fingerprints
