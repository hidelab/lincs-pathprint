# Script's Job: Produces experiment-fingerprints
#
# Author: Sokratis Kariotis
# Status: v1
# Timestamp: 07.11.2017
# Script that initiates a number of cores, each one utilizing a number of 
# experiments to run the gctx2fingerprint.R script and produce the 
# experiment-fingerprints.

# Script start
# load libraries
library(GMAfunctions)
library(parallel)
library(foreach)
library(doMC)

# Go to the parallel path
setwd("~/LINCS project/LINCS specific scripts/lincs_parallel")

# Select the number of core
number_of_cores <- 14
registerDoMC(cores = number_of_cores)

# Load the sample set
load("~/LINCS project/LINCS specific scripts/sampleset.RData")

# The experiment, Entrez history and gctx matrices file paths
gctx_exp_path <- "~/LINCS project/LINCS specific scripts/gctx_experiments/"
history_file_path <- "~/LINCS project/LINCS specific scripts/Entrez Gene History/EntrezHistoryList.RData"
matrices_path <- "~/LINCS project/LINCS specific scripts/gctx_Folder/"

# These scripts run for single chip enrichment

a = mclapply(source("lincs_parallel_1.R"))
b = mclapply(source("lincs_parallel_2.R"))
c = mclapply(source("lincs_parallel_3.R"))
d = mclapply(source("lincs_parallel_4.R"))
e = mclapply(source("lincs_parallel_5.R"))
f = mclapply(source("lincs_parallel_6.R"))
g = mclapply(source("lincs_parallel_7.R"))
h = mclapply(source("lincs_parallel_8.R"))
i = mclapply(source("lincs_parallel_9.R"))
j = mclapply(source("lincs_parallel_10.R"))
k = mclapply(source("lincs_parallel_11.R"))
l = mclapply(source("lincs_parallel_12.R"))
m = mclapply(source("lincs_parallel_13.R"))
n = mclapply(source("lincs_parallel_14.R"))

# Collecting processes/threads
print("Collecting threads...")
Sys.sleep(10)
x=do.call(rbind,list(a,b,c,d,e,f,g,h,i,j,k,l,m,n))
print("Finished creating experiment fingerprints successfully!")
