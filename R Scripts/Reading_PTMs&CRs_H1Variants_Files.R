############################################################################################################################################
#                                       Reading PTMs&CRs and Histone 1 variants files 
############################################################################################################################################
# 1| PTMs and CRs data normalized by length: Reading the files -----------------------------------------------------------------------------
# Setting Working Directory & Loading the Scripts 
setwd(PTMsDir)
# Storing all the PTMs & CRs filenames of interest
print("2| LOADING OF THE BED FILES")
myFiles <- list.files(pattern="*.bed")
FileNames <- myFiles

List <- data.frame()
Bands_Names <- c()

# For each Cytoband type (5), for each PTM (12), do:
print(paste("   PTMs Analyzed: ", paste(PTMs,collapse=",")))
for (band in 1:length(CytoBands)){
  for (PTM in 1:length(PTMs)){
      if (substr(CytoBands[band],1,4)=="Gneg") {
        FileName <- myFiles [grepl(as.character(PTMs[PTM]),myFiles) & grepl(substr(CytoBands[band],1,4),myFiles)] 
      } 
      if (substr(CytoBands[band],1,4)=="Gpos") {
        FileName <- myFiles [grepl(as.character(PTMs[PTM]),myFiles) & grepl(CytoBands[band],myFiles)] 
      }
      
      # Creating a discrete variable for each PTM within each Cytoband type
      PTM_band <- read.table (FileName, header=FALSE)
      PTM_band$V8 <- paste(PTM_band$V4, PTM_band$V5, PTM_band$V1, sep="_")
      PTM_band$V8 <- gsub("gneg","Gneg",PTM_band$V8); PTM_band$V8 <- gsub("gpos","Gpos",PTM_band$V8) 
      PTM_band$V9 <- paste(PTMs[PTM])
      PTM_band <- PTM_band[match(rownames(GnegGposGroups)[which(GnegGposGroups$CytobandGroup==CytoBands[band])], PTM_band$V8) ,]
      
      PTM_band$V8 <- gsub("Gneg",CytoBands[band],PTM_band$V8)
      do.call("<-",list(paste(PTMs[PTM],CytoBands[band], sep="_"), PTM_band))
      
      # These variables will store 
      List <- rbind(List, PTM_band)
      
  }
  Bands_Names <- c(Bands_Names, PTM_band$V8)
  
  # Creating a discrete variable storing all the PTMs within each Cytoband groups
  do.call("<-",list(paste("Total",CytoBands[band],"PTMsCRs", sep="_"), 
                    as.data.frame(matrix(data=List[ ,7], nrow = length(PTM_band$V8), ncol=length(PTMs), dimnames = list(c(PTM_band$V8),c(PTMs)) )))) 
  
  # Cleaking up the variables in order to start a new cycle of the loop
  List <- data.frame()
}
save(Bands_Names, file= paste(ResultsDir,"/Bands_Names.Rdata",sep=""))

# Creating a discrete variable storing all the PTMs within all the Cytoband groups and storing the number of bands per group 
TOTAL_PTMsCRs <- data.frame()
for (band in CytoBands){
  TOTAL_PTMsCRs <- rbind(TOTAL_PTMsCRs, eval(as.name(paste("Total_",band,"_PTMsCRs",sep="") )) )
}

# Removing undesired bands
Discarded_Bands <- c()
for (band in DiscardedBands){ # 1.Searching ALL the bands that present the undesired characteristics
  Discarded_Bands <- c(Discarded_Bands, c(grep(band,rownames(TOTAL_PTMsCRs), value=TRUE)))
}
for (band in Discarded_Bands){ # 2.Removing each of the bands from "TOTAL_PTMsCRs"
  if (band %in% rownames(TOTAL_PTMsCRs) == TRUE ){
    TOTAL_PTMsCRs <- TOTAL_PTMsCRs[-which(rownames(TOTAL_PTMsCRs) %in% band), ]
  }
}

# Storing the number of bands per group
Bands_Length <- c()
for (band in CytoBands){
  Bands_Length <- c(Bands_Length, length(grep(band,rownames(TOTAL_PTMsCRs))))
}

# Checking the content and the dimensions of the variables "TOTAL_PTMsCRs", which stores all the information
# Scaling the data  
TOTAL_PTMsCRs_scaled <- as.data.frame(apply(TOTAL_PTMsCRs, 2, scale))
rownames(TOTAL_PTMsCRs_scaled) <- rownames(TOTAL_PTMsCRs)

# In this case, the variable present a manageable name but it is going to be stored including the cell line type 
save(TOTAL_PTMsCRs, file = paste(ResultsDir,"/TOTAL_PTMsCRs_",CellLine,".RData",sep=""))


# 2| H1 variants data normalized: Reading the files ----------------------------------------------------------------------------------------
# Setting Working Directory & Loading the Scripts 
setwd(HistonesDir)

# Storing all the histone variants' filenames of interest
myFiles <- list.files(pattern="*.bed")
FileNames <- myFiles

List <- data.frame()

# For each Cytoband type (5), for each H1 variant (8), do:
for (band in 1:length(CytoBands)){
  for (histone in 1:length(HistonesFiles)){
    if (substr(CytoBands[band],1,4)=="Gneg") {
      FileName <- myFiles [grepl(HistonesFiles[histone],myFiles) & grepl(substr(CytoBands[band],1,4),myFiles)][1]
    } 
    if (substr(CytoBands[band],1,4)=="Gpos") {
      FileName <- myFiles [grepl(HistonesFiles[histone],myFiles) & grepl(CytoBands[band],myFiles)][1]
    }
  
    # Creating a discrete variable for each PTM within each Cytoband type
    H1_band <- read.table (FileName, header=FALSE)
    H1_band$V8 <- paste(H1_band$V4, H1_band$V5, H1_band$V1, sep="_")
    H1_band$V8 <- gsub("gneg","Gneg",H1_band$V8); H1_band$V8 <- gsub("gpos","Gpos",H1_band$V8) 
    H1_band$V9 <- paste(Histones[histone])
    H1_band <- H1_band[match(rownames(GnegGposGroups)[which(GnegGposGroups$CytobandGroup==CytoBands[band])], H1_band$V8) ,]
    
    H1_band$V8 <- gsub("Gneg",CytoBands[band],H1_band$V8)
    do.call("<-",list(paste(Histones[histone],CytoBands[band], sep="_"), H1_band))
    
    # These variables will store 
    List <- rbind(List, H1_band)
  }
  # Creating a discrete variable storing all the H1s within each Cytoband groups
  do.call("<-",list(paste("Total",CytoBands[band],"H1", sep="_"), 
                    as.data.frame(matrix(data=List[ ,6], nrow = length(H1_band$V8), ncol=length(Histones), dimnames = list(c(H1_band$V8),c(Histones)) )))) 
  
  # Cleaking up the variables in order to start a new cycle of the loop
  List <- data.frame()
}

# Creating a discrete variable storing all the H1 variants within all the Cytoband groups
TOTAL_H1Variants <- data.frame()
for (band in CytoBands){
  TOTAL_H1Variants <- rbind(TOTAL_H1Variants, eval(as.name(paste("Total_",band,"_H1",sep="") )) )
}

# Increasing values (they tend to be little) and removing undesired bands from H1 variants data.frame
for (band in Discarded_Bands){ # 2.Removing each of the bands from "TOTAL_PTMsCRs"
  if (band %in% rownames(TOTAL_H1Variants) == TRUE ){
    TOTAL_H1Variants <- TOTAL_H1Variants[-which(rownames(TOTAL_H1Variants) %in% band), ]
  }
}
TOTAL_H1Variants <- TOTAL_H1Variants[complete.cases(TOTAL_H1Variants),]
if (CellLine == "HeLa1"){
  TOTAL_H1Variants <- TOTAL_H1Variants*100
}

# Scaling the data  
H1_variants_scaled <- as.data.frame(apply(TOTAL_H1Variants, 2, scale))
rownames(H1_variants_scaled) <- rownames(TOTAL_H1Variants)

# Checking the content and the dimensions of the variables "Total_GposGneg_PTMsCRs", which stores all the information
save(TOTAL_H1Variants, file = paste(ResultsDir,"/TOTAL_H1Variants_",CellLine,".RData",sep=""))


