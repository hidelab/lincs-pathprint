# Post processing metadata matrix creation step 4
# Author: Sokratis Kariotis
# Description: Script to compile a table of metadata for the fingerprint matrix.
# Due to the nature of the metadata found in GSE70138_Broad_LINCS_inst_info.txt, 
# the LINCS.metadata.matrix will differ from the correspodning GEO file. In addition
# to the present columns, we have GPL and GSE columns. Subject to change!

# Define path of repository
definePath<-function(){
readline("define pathprint repository, or blank for default (~/LINCS project/LINCS specific scripts/post_processing_general_files/) : ")
}
pathprintRepository<-definePath()

if (pathprintRepository == ""){
	pathprintRepository<-"~/LINCS project/LINCS specific scripts/post_processing_general_files/"
	}

# Load fingerprint matrix
defineFile<-function(){
readline("enter filename: ")
}
dataPath<-"~/LINCS project/LINCS specific scripts/post_processing_general_files/"

dir(path = dataPath, pattern = "platform.frame.20")
print("select platform file")
platformfile<-defineFile()
load(paste(dataPath, platformfile, sep = ""))

# Read metadata table
LINCS.metadata.matrix = read.table("~/LINCS project/LINCS specific scripts/gctx_Folder/GSE70138_Broad_LINCS_inst_info.txt", fill = TRUE, header = TRUE)

# Construct metadata table
LINCS.metadata.matrix <- as.data.frame(append(LINCS.metadata.matrix, list(Species = "Homo sapiens"), after = 1))
LINCS.metadata.matrix <- as.data.frame(append(LINCS.metadata.matrix, list(GPL = "GPL20573"), after = 1))
LINCS.metadata.matrix <- as.data.frame(append(LINCS.metadata.matrix, list(GSE = "GSE70138"), after = 1))

# Filtering out incomplete lines
LINCS.metadata.matrix <- LINCS.metadata.matrix[!(is.na(LINCS.metadata.matrix$inst_id == "") | LINCS.metadata.matrix$cell_id == "" | LINCS.metadata.matrix$det_plate == ""), ]

# Convert to ascii strings for better package compatibility

LINCS.metadata.matrix.temp<-LINCS.metadata.matrix

# first convert "ï" to "i"
LINCS.metadata.matrix.temp[,5] <- gsub("ï", "i", LINCS.metadata.matrix.temp[,5])
LINCS.metadata.matrix.temp[,6] <- gsub("ï", "i", LINCS.metadata.matrix.temp[,6])
LINCS.metadata.matrix.temp[,7] <- gsub("ï", "i", LINCS.metadata.matrix.temp[,7])

LINCS.metadata.matrix.temp[,5] <- iconv(LINCS.metadata.matrix.temp[,5], "latin1", "ASCII", "byte")
LINCS.metadata.matrix.temp[,6] <- iconv(LINCS.metadata.matrix.temp[,6], "latin1", "ASCII", "byte")
LINCS.metadata.matrix.temp[,7] <- iconv(LINCS.metadata.matrix.temp[,7], "latin1", "ASCII", "byte")

# Show changes
all.equal(LINCS.metadata.matrix.temp, LINCS.metadata.matrix)

LINCS.metadata.matrix<-LINCS.metadata.matrix.temp

# Save file into fingerprint package
try(system(paste("cp ",
             pathprintRepository,
             "LINCS.metadata.matrix.RData ",
             dataPath,
             "LINCS.metadata.matrix.RData.old",
             sep = "")))

save(LINCS.metadata.matrix, file = paste(pathprintRepository, "LINCS.metadata.matrix.RData", sep = ""), compress = TRUE)

# end
