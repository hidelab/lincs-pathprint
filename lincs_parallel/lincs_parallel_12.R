# Fingerprinting worker script - running single chip enrichment for a group of
# experiments. Delay allows time to monitor process starting up

# Loading the neccessary scripts and files
source("/shared/hidelab2/shared/Sokratis/pathprint/fingerprint/package/GMAfunctions/R/GSMtable2exprs.R")
source("/shared/hidelab2/shared/Sokratis/pathprint_lincs/gctx2fingerprint.R")
load("/shared/hidelab2/shared/Sokratis/pathprint_lincs/genesets.rda")
load("/shared/hidelab2/shared/Sokratis/pathprint_lincs/pathprint.Hs.gs.rda")
load("/shared/hidelab2/shared/Sokratis/pathprint_lincs/sampleset.RData")

print("preparing for script launch, 5s synch delay")
Sys.sleep(5)

# Load the sample set in
current_samples <- sampleset[[12]];
total <- length(current_samples)

# load contents of directory once
# this used to be within the loop
e_path = "~/LINCS project/LINCS specific scripts/sq_gctx_experiments"
dir<-dir(path = e_path)

for (i in 1:length(current_samples)){
	print(paste(i, "of", total, sep = " "))
    filename<-paste("sq", current_samples[i], sep = "_")
    
    # Only fingerprint experiments when you can't find their sq_file (experiment
    # fingerprint)
	if (filename %in% dir){
		print ("File already analyzed")
	}
	else if (!(filename %in% dir(path = path))){
	    #try-catching to skip erroneous experiments to avoid crushing
	    temp<-tryCatch(gctx2fingerprint(
    	    GSM = current_samples[i], 
    	    GEOthreshold = FALSE,
    	    GEOpath = GEOpath,
    	    geneset = "KEGG and Wikipathways and static",
    	    enrichmentMethod = "SCE",
    	    transformation = "squared.rank",
    	    statistic = "mean",
    	    normalizedScore = FALSE,
    	    progressBar = FALSE
    	# Do nothing if an error occurs [add handlers later]*
        ),error=function(e)NULL)
	    
	    temp1<-list("name" = temp)
	    names(temp1)<-current_samples[i]
	    save(temp1, file = paste(e_path, filename, sep =""))
		  }
		}
