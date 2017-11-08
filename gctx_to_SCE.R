# My Job: Get the .gctx matrix and start threads which create fingerprint files per experiment 
#
# Author: Sokratis Kariotis
# Status: v1
# Timestamp: 07.11.2017
# Script to download or read the Level 3 GSE from LINCS 
# (GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx) and create
# the sampleSet.RData file that holds the gene expression of all genes in one expreriment
# (per column of .gctx matrix). Then the sampleSet object  gets split up and passed to 
# the fingerprinting scripts. Currently split into 14 parts as the fingerprinting script 
# is passed to 14 cores.

# For this script we need:
# > (optional) The .gctx matrix downloaded here "~/LINCS project/LINCS specific scripts/gctx_Folder"

# Script start
# load library tha handles .gctx files
library(cmapR)
library(GMAfunctions)

# Check if the .gctx matrix is in the gctx_Folder. If not download manually and add it.
if(!file.exists("~/LINCS project/LINCS specific scripts/gctx_Folder/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328.gctx",showWarnings = TRUE)[1])
    stop("Cannot find .gctx file.")

# column/experiment and gene count
experiments_count <- length(parse.gctx("~/LINCS project/LINCS specific scripts/gctx_Folder/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328.gctx", 
                   rid=1:1, matrix_only = TRUE)@cid)
gene_count <- length(parse.gctx("~/LINCS project/LINCS specific scripts/gctx_Folder/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328.gctx", 
                   cid=1:1, matrix_only = TRUE)@rid)

# The experiment files are here
exp_path <- "~/LINCS project/LINCS specific scripts/gctx_experiments/"

# Locate the Entrez history file to be used for updating
historyFile <- "~/LINCS project/LINCS specific scripts/Entrez Gene History/EntrezHistoryList.RData"

# We have to split the matrix column-wise/experiment-wise since its 12.6 Gb altogether
complete_matrix <- matrix(nrow = gene_count, ncol = experiments_count/10)
chunk <- as.integer(experiments_count / 1  )

for (exp in 1:experiments_count) {
        temp_exp <- parse.gctx("~/LINCS project/LINCS specific scripts/gctx_Folder/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328.gctx",
                         cid=exp:exp)@mat
        
        # Update Entrez Ids
        new <- entrezUpdate(rownames(temp_exp), historyFile = historyFile)
        old <- rownames(temp_exp)
        print(Reduce(setdiff, strsplit(c(new, old), split = "")))
        
        # rownames(temp_exp)<-entrezUpdate(rownames(temp_exp), historyFile = historyFile)
        
        exp_name <- gsub(":", "_", colnames(temp_exp))
        save_name <- paste(exp_path,exp_name,".RData", sep="")
        save(temp_exp, file = save_name)
        print(paste("Downloaded experiment ",exp_name))
        print(paste(exp," experiments complete."))
}


for (matrix_chunk in seq(1,experiments_count/10, by = chunk)){
    current_start <- matrix_chunk

    # Check if its last loop
    if (matrix_chunk + chunk > experiments_count) current_end <- experiments_count
    else current_end <- matrix_chunk + chunk - 1

    print(paste("Select ", current_start, " to ", current_end ))
    temp_mat <- parse.gctx("~/LINCS project/LINCS specific scripts/gctx_Folder/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328.gctx",
                     cid=current_start:current_end, matrix_only = TRUE)@mat

    # Add to the complete matrix
    complete_matrix <- cbind(complete_matrix,temp_mat)
    print(paste("Added experiments ", current_start, " to ", current_end , " to the complete matrix"))
}

# save(complete_matrix, file = "~/LINCS project/LINCS specific scripts/complete_matrix.RData")


# load("/shared/hidelab2/shared/Sokratis/pathprint/fingerprint/data/chipframe36.rda")
# 
# 
# # Use GEOmetaDB to retrieve the arrays
# wd<-getwd()
# library(GEOmetadb)
# 
# # save GEOmetadb database locally on Iceberg
# setwd("/shared/hidelab2/shared/Sokratis/pathprint/Databases/GEOmetadb")
# 
# # download and extract database
# #getSQLiteFile()
# con <- dbConnect(SQLite(), "GEOmetadb.sqlite")
# 
# # Retrieve information from GEO
# #chipnames<-names(chipframe)
# chipnames<-"GPL7202"
# chipvector<-vector("character", 0)
# print("Querying GEOmetaDB")
# for (i in 1:length(chipnames)){
#   query<-paste("SELECT gsm FROM gsm WHERE gpl = '", chipnames[i], "'", sep = "")
#   print(query)
#   gsmList<-dbGetQuery(con, query)
#   chipvector<-c(chipvector, unlist(gsmList, use.names = FALSE))
#   }
# 
# save(gsmList, file = "gsmList.RData")
# print("gsm list")
# 
#   
# GEOsamples.new<-chipvector
# dbDisconnect(con)
# 
# # reset working directory
# setwd(wd)
# 
# # scramble sample names to avoid bias in processing time asscociated with the larger chips
# # apportion into 12 or 14
# print(paste(length(GEOsamples.new), "to be fingerprinted", sep = " ")) 
# GEOsamples.rand <- sample(GEOsamples.new)
# len <- length(GEOsamples.rand) %/% 14
# 
# sampleset<-vector("list", 14)
# for (i in 1:length(sampleset)){
# 	sampleset[[i]] <- GEOsamples.rand[(1+((i-1)*len)):(i*len)]
# 	}
# sampleset[[14]] <- GEOsamples.rand[(1+(13*len)):length(GEOsamples.rand)]
# save(sampleset, file = "sampleset.RData")
# print("Sampleset saved!")
