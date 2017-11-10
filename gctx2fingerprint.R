# Script's Job: Produces experiment-fingerprints
#
# Author: Sokratis Kariotis
# Status: v1
# Timestamp: 08.11.2017
# Script that accepts a experiment name, loads it's matrix and produces the
# 

# Script start - start editing here 10-11-2017: I need to merge 
# geo2fingerprint-exprs2fingerprint-single.chip.enrichment.R in this one 
# script (as lots of code that handles probes etc is useless in this project)
# Remember to code and test on HPC. Upload evrything on my pathprint-lincs folder

# Load the genesets
load("/shared/hidelab2/shared/Sokratis/pathprint/fingerprint/package/pathprint/data/genesets.rda")

# This is a function to produce a normalized fingerprint directly from an expression 
# matrix. EnrichmentMethod - algorithm used to calculate gene enrichment score and
# it can only be  single.chip.enrichment (SCE). Transformation, statistic, 
# normalizedScore and progressBar are relevant to SCE.
gctx2fingerprint<-function(	GSM, GEOpath = tempdir(), GEOthreshold = TRUE,
							geneset = "KEGG and Wikipathways and static",
							enrichmentMethod = "SCE", transformation = "squared.rank",
							statistic = "mean", normalizedScore = FALSE,
							progressBar = FALSE
							)
	{

	# load list of GEO files - for speed only do this if object not found already
	if (!(exists("GEOfiles"))){	
	assign("GEOfiles", gsub(".soft", "", dir(path = GEOpath)), .GlobalEnv)
		}
	
	# retrieve GEO file and strip out required meta data
	library(GEOquery)
	print("Downloading/loading GEO file")
	# Only downloaded file if not already in file path
	if (GSM %in% GEOfiles){
		print("using locally saved GEO file")
		try(geo <- getGEO(GSM, file = (paste(GEOpath, GSM, ".soft", sep = ""))))
		}
	else if (!(GSM %in% GEOfiles)){
		try(geo<-getGEO(GSM, destdir = GEOpath))
		}
	

	return(list(SCG = geo.SCG, species = species, platform = platform))
	}

