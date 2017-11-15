# Compile metadata for POE matrix
# Author: Sokratis Kariotis
# Script to compile a table of metadata for the fingerprint matrix



# I AM HERE: must find metadata etc


# define path of repository
definePath<-function(){
readline("define pathprint repository, or blank for default (~/LINCS project/LINCS specific scripts/post_processing_general_files/) : ")
}
pathprintRepository<-definePath()
if (pathprintRepository == ""){
	pathprintRepository<-"~/LINCS project/LINCS specific scripts/post_processing_general_files/"
	}

# load fingerprint matrix
defineFile<-function(){
readline("enter filename: ")
}
dataPath<-"~/LINCS project/LINCS specific scripts/post_processing_general_files/"

dir(path = dataPath, pattern = "platform.frame.20")
print("select platform file")
platformfile<-defineFile()
load(paste(dataPath, platformfile, sep = ""))

# Read metadata table
metadata_table = read.table("~/LINCS project/LINCS specific scripts/gctx_Folder/GSE70138_Broad_LINCS_inst_info.txt", fill = TRUE)

gsm.char = dbGetQuery(con, "select gsm, characteristics_ch1 from gsm")
gsm.source = dbGetQuery(con, "select gsm, source_name_ch1 from gsm")
gsm.species = dbGetQuery(con, "select gsm, organism_ch1 from gsm")
gsm.gse<-geoConvert(rownames(platform.frame), c("gse"))$gse
gsm.gpl<-geoConvert(rownames(platform.frame), c("gpl"))$gpl
dbDisconnect(con)

# construct metadata table
GEO.metadata.matrix<-as.data.frame(matrix(nrow = length(rownames(platform.frame)), ncol = 7))
colnames(GEO.metadata.matrix)<-c( "Title", "GSE", "GPL", "Species", "Source", "Characteristics")
GEO.metadata.matrix$Title<-as.character(rownames(platform.frame))
GEO.metadata.matrix$GSE<-gsm.gse$to_acc[match(GEO.metadata.matrix$GSM, gsm.gse$from_acc)]
GEO.metadata.matrix$GPL<-gsm.gpl$to_acc[match(GEO.metadata.matrix$GSM, gsm.gpl$from_acc)]
GEO.metadata.matrix$Species<-gsm.species$organism_ch1[match(GEO.metadata.matrix$GSM, gsm.species$gsm)]
GEO.metadata.matrix$Source<-gsm.source$source_name_ch1[match(GEO.metadata.matrix$GSM, gsm.source$gsm)]
GEO.metadata.matrix$Characteristics<-gsm.char$characteristics_ch1[match(GEO.metadata.matrix$GSM, gsm.char$gsm)]

#####
# Convert to ascii strings for better package compatibility

GEO.metadata.matrix.temp<-GEO.metadata.matrix

# first convert "ï" to "i"
GEO.metadata.matrix.temp[,5] <- gsub("ï", "i", GEO.metadata.matrix.temp[,5])
GEO.metadata.matrix.temp[,6] <- gsub("ï", "i", GEO.metadata.matrix.temp[,6])
GEO.metadata.matrix.temp[,7] <- gsub("ï", "i", GEO.metadata.matrix.temp[,7])

GEO.metadata.matrix.temp[,5] <- iconv(GEO.metadata.matrix.temp[,5], "latin1", "ASCII", "byte")
GEO.metadata.matrix.temp[,6] <- iconv(GEO.metadata.matrix.temp[,6], "latin1", "ASCII", "byte")
GEO.metadata.matrix.temp[,7] <- iconv(GEO.metadata.matrix.temp[,7], "latin1", "ASCII", "byte")



# show changes
all.equal(GEO.metadata.matrix.temp, GEO.metadata.matrix)

GEO.metadata.matrix<-GEO.metadata.matrix.temp

# save file into fingerprint package
try(system(paste("cp ",
             pathprintRepository,
             "GEO.metadata.matrix.RData ",
             dataPath,
             "GEO.metadata.matrix.RData.old",
             sep = "")))

save(GEO.metadata.matrix, file = paste(pathprintRepository, "GEO.metadata.matrix.RData", sep = ""), compress = TRUE)

# end
