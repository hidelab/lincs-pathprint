# lincs-pathprint
New pipeline for pathprinting LINCS data

#### Status
Operation | Status | Note
--- | --- | ---
Experiment download | `Imperfect` | Less experiments than expected, see [Issue 2](https://github.com/hidelab/lincs-pathprint/issues/2)
Experiment Fingerprinting | `Complete` | -
Post-processing step 1 | `Complete` | -
Post-processing step 2 | `Running on Iceberg` | -

## Scripts

### Pre-processing

1. **lincs_pathprint_download.R** : Script to download the Level 3 GSE from LINCS (GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx) and split them in a number of chunks (adjustable by number_of_chunks variable). Before saving to an equal number of matrices, the gene IDs are updated using the entrezUpdate.R script from the package GMAfunctions. 

2. **lincs_pathprint_runner.R** : Script that initiates a number of cores, each one utilizing a number of experiment expression files to run a *lincs_parallel_x.R* script.

3. **lincs_parallel_x.R** : Parallel scripts(14). Each initiates the experiment-fingerprinting process of multiple experiments by running the *gctx2fingerprint.R* script.


3. **gctx2fingerprint.R** : Script that accepts a experiment name, loads it's matrix and produces the experiment fingerprints. It's a merge of geo2fingerprint.R and exprs2fingerprint.R. It utilizes *custom.single.chip.enrichment.R* to produce the experiment fingerprints.

4. **custom.single.chip.enrichment.R** : Slightly tweaked version of *single.chip.enrichment.R* that omits some unrelated parameters.

### Post-processing

1. **Fingerprint_post_processing_step1_lincs.R** :  Creates the *sq_.frame.xxxx-xx-xx.RData* file (SCE dataframe) and the *platform.frame.xxxx-xx-xx.RData* file (platform dataframe). It also creates the *sq_.pathway.SCE.x.RData* files for each pathway.

2. **Fingerprint_post_processing_step2_lincs.R** :  Creates the POE (probability of expression) file for each pathway. 

3. **Fingerprint_post_processing_step3_lincs.R** :  Creates the *sq_.matrix.xxxx-xx-xx.RData* POE matrix file (one file for all pathways). 

4. **matrix_metadata_lincs.R** : Creates the *LINCS.metadata.matrix.RData* file.
 
5. **constructing POE thresholds.R** : Obtains POE threshold values and creates *LINCS.fingerprint.matrix.RData* file.


## Notes

1. Optimize speed in pre-processing 
