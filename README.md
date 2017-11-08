# lincs-pathprint
New pipeline for pathprinting LINCS data

## Scripts
1. **gctx_to_SCE.R** : Download/read the Level 3 GSE from LINCS (GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx) and create the sampleSet.RData file that holds all expreriments (each column of .gctx matrix) to be used in the next step. The downloaded experiments are split into 15 chunks to increase download nad updating speed as well as be ready for the multi-thread fingerprinting in the following scripts.
