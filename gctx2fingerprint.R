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

this.dir <- dirname(parent.frame(2)$ofile)
source(paste(this.dir, "exprs2fingerprint_options.R", sep="/"))
load("/shared/hidelab2/shared/Sokratis/pathprint/fingerprint/package/pathprint/data/genesets.rda")

geo2fingerprint<-function(	GSM, GEOpath = tempdir(), GEOthreshold = TRUE,
							geneset = "KEGG and Wikipathways and static",
							enrichmentMethod = "SCE", transformation = "rank",
							statistic = "mean", normalizedScore = FALSE,
							progressBar = FALSE
							)
# this is a function to produce a normalized fingerprint directly from a GEO accession number
# it is designed to be run on the server_updated so that normalized dataframe is not included
# enrichmentMethod - algorithm used to calculate gene enrichment score
# Algorithm can be single.chip.GSEA (SCG) or single.chip.enrichment (SCE)
# Note that transformation, statistic, normalizedScore and progressBar are only relevant to SCE

	{

#	for hpc111 previously used GEOpath <- "/home/galtschu2/Documents/Databases/GEOfiles/"

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
	
	# There now follows a series of checks to ensure that the correct data and meta data is present
	
	if (exists("geo")){					
		if (!(is.null(Meta(geo)$organism_ch1))){organism <- Meta(geo)$organism_ch1}
		else if (!(is.null(Meta(geo)$organism_ch2))){organism <- Meta(geo)$organism_ch2}
		else {organism <- "species field empty"}
# keep conversion to "human" and "mouse" to retain backwards compatibility with old dataframe
		if (organism == "Mus musculus"){
			species <- "mouse"
			}
		else if (organism == "Homo sapiens"){
			species <- "human"
			} else {
				species <- organism
				}
		platform<-Meta(geo)$platform_id
		if (!(exists("genesets"))) stop("Need to load genesets list")
		if ((length(Table(geo)) > 0) & (organism %in% names(genesets))){
		exprs<-GSMtable2exprs(Table(geo))
			geo.SCG <- exprs2fingerprint_options(											exprs = exprs, platform = platform, species = species, GEOthreshold = GEOthreshold, geneset = geneset,
											transformation = transformation,
											statistic = statistic, normalizedScore = normalizedScore, progressBar = progressBar
										)	

			}

		else if (!(platform %in% names(chipframe))){
			print (paste("Returning null values as no available chip annotation for", platform, sep = " "))
			geo.SCG = 0
			geo.norm = 0
			if(!(exists("missingchips"))){ 
		 		missingchips <- vector("list", 0)
		 		}
		 	if (!(platform %in% missingchips)){
		 		missingchips <- append(missingchips, platform)
		 		assign("missingchips", missingchips, envir = .GlobalEnv)
		 		}
			}
		else if (((platform %in% names(chipframe)) & !(organism %in% names(genesets)))){
			print (paste("Returning null values as analysis does not include", organism, sep = " "))
			geo.SCG = 0
			geo.norm = 0
			if(!(exists("missingspecies"))){ 
		 		missingspecies <- vector("list", 0)
		 		}
		 	if (!(organism %in% missingspecies)){
		 		missingspecies <- append(missingspecies, organism)
		 		assign("missingspecies", missingspecies, envir = .GlobalEnv)
		 		}
		 	}
		 else if ((length(Table(geo)) < 1)){
			print (paste("Returning null values as no array data for", GSM, sep = " "))
			geo.SCG = 0
			geo.norm = 0
			if(!(exists("badfiles"))){ 
		 		badfiles <- vector("list", 0)
		 		}
		 	if (!(GSM %in% badfiles)){
		 		badfiles <- append(badfiles, GSM)
		 		assign("badfiles", badfiles, envir = .GlobalEnv)
		 		}
	 		
			}
		}
		else if (!(exists("geo"))){
			print ("Error in GEO file")	
			geo.SCG = 0
			geo.norm = 0
			species = "unknown"
			platform = "unknown"
			if(!(exists("badfiles"))){ 
		 		badfiles <- vector("list", 0)
		 		}
		 	if (!(GSM %in% badfiles)){
		 		badfiles <- append(badfiles, GSM)
		 		assign("badfiles", badfiles, envir = .GlobalEnv)
		 		}
		 	}	

	return(list(SCG = geo.SCG, species = species, platform = platform))
	}

