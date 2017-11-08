# My Job: Get the .gctx matrix and start threads which create fingerprint files per experiment 
#
# Author: Sokratis Kariotis
# Status: v1
# Timestamp: 07.11.2017
# Script to download the Level 3 GSE from LINCS 
# (GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx) and 
# split them in a number of chunks (adjustable by number_of_chunks variable).
# Before saving to an equal number of matrices, the gene IDs are updated using
# the entrezUpdate.R script from the package GMAfunctions. Finally, the script 
# initiates a number of cores each one handling running the gctx2fingerprint.R
# script to produce the experiment fingerprints.


# For this script we need:
# > (optional) The .gctx matrix downloaded here "~/LINCS project/LINCS specific scripts/gctx_Folder"

# Script start
# load library tha handles .gctx files
library(cmapR)
library(GMAfunctions)

# Select the number of chunks
number_of_chunks <- 15


# Check if the .gctx matrix is in the gctx_Folder. If not download manually and add it.
if(!file.exists("~/LINCS project/LINCS specific scripts/gctx_Folder/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328.gctx",showWarnings = TRUE)[1])
    stop("Cannot find .gctx file.")

# column/experiment and gene count
experiments_count <- length(parse.gctx("~/LINCS project/LINCS specific scripts/gctx_Folder/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328.gctx", 
                   rid=1:1, matrix_only = TRUE)@cid)
gene_count <- length(parse.gctx("~/LINCS project/LINCS specific scripts/gctx_Folder/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328.gctx", 
                   cid=1:1, matrix_only = TRUE)@rid)

# The experiment files are here
gctx_exp_path <- "~/LINCS project/LINCS specific scripts/gctx_experiments/"

# Locate the Entrez history file to be used for updating
historyFile <- "~/LINCS project/LINCS specific scripts/Entrez Gene History/EntrezHistoryList.RData"

# We have to split the matrix column-wise/experiment-wise since its 12.6 Gb altogether
chunk <- as.integer(experiments_count / number_of_chunks)

# gctx matrices are located here
exp_path <- "~/LINCS project/LINCS specific scripts/gctx_Folder/gctx_matrix_"

sampleset <- vector("list", number_of_chunks)
sampleset_count <- 1

for (matrix_chunk in seq(1,experiments_count, by = chunk)){
    current_start <- matrix_chunk

    # Check if its last loop
    if (matrix_chunk + chunk > experiments_count) current_end <- experiments_count
    else current_end <- matrix_chunk + chunk - 1

    partial_gctx_matrix <- parse.gctx("~/LINCS project/LINCS specific scripts/gctx_Folder/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328.gctx",
                     cid=current_start:current_end, matrix_only = TRUE)@mat
    
    # Updating the Entrez IDs
    rownames(partial_gctx_matrix)<-entrezUpdate(rownames(partial_gctx_matrix), historyFile = historyFile)
    print("Entrez IDs updated.")
    
    # Fill the sample set
    sampleset[[sampleset_count]] <- colnames(partial_gctx_matrix)
    
    save_name <- paste(exp_path,current_start,"-",current_end,".RData", sep="")
    save(partial_gctx_matrix, file = save_name)
    print(paste("Downloaded and saved experiments ", current_start, " to ", current_end, " in ", save_name))
    sampleset_count <- sampleset_count + 1
    
    for (e in 1:length(colnames(partial_gctx_matrix))) {

      expression_column <- matrix(partial_gctx_matrix[,e])
      gene_column <- matrix(rownames(partial_gctx_matrix))
      exp_matrix <- cbind(gene_column,expression_column)
      
      experiment_name <- colnames(partial_gctx_matrix)[e]
      experiment_name <- gsub(":", "_", experiment_name)
      full_experiment_name <- paste(gctx_exp_path,experiment_name,".RData", sep = "")
      save(exp_matrix, file = full_experiment_name)
    }
    
}

save(sampleset, file = "sampleset.RData")
print("Sample set saved.")
