##############################################################################################################################################
#                                                              AB COMPARTMENTS ANALYSIS
##############################################################################################################################################
# This R script has been designed to analyze A/B compartments obtained from a Hi-C experiment.It is composed of four main sections: 
# H1 Abundance in A/B Compartments, PTMs&CRs Abundance in A/B Compartments, A/B Compartments vs CytoBands and A/B Compartments vs TADs Groups

##############################################################################################################################################
# H1 Abundance in AB Compartments
##############################################################################################################################################

#--- Setting Working Directory and required variables ----------------------------------------------------------------------------------------
setwd("/home/andrea/Bioinformatica_Núria/HiC_TADs/H1abundance_TADs/")
AB_WT2 <- read.table("H1.2endo_r2_hg19.bdg_AB_01_WT.bed", header = F, sep = "\t")
AB_KD2 <- read.table("H1.2endo_r2_hg19.bdg_AB_02_KD.bed",header = F, sep = "\t")
AB_WTx <- read.table("H1Xendo_r3_hg19.bdg_AB_01_WT.bed", header = F, sep = "\t")
AB_KDx <- read.table("H1Xendo_r3_hg19.bdg_AB_02_KD.bed",header = F, sep = "\t")

#--- BoxPlot representing the H1 abundance in A/B Compartments (WT) --------------------------------------------------------------------------
pdf("/home/andrea/Bioinformatica_Núria/HiC_TADs/T47D_H1abundance_AB.pdf")
boxplot (AB_WT2$V5[AB_WT2$V4=="A"], AB_WT2$V5[AB_WT2$V4=="B"], 
         AB_WTx$V5[AB_WTx$V4=="A"], AB_WTx$V5[AB_WTx$V4=="B"],outline = F, at = c(1,2,3.5,4.5), 
         ylab = "ChIP-Seq Signal", main = "T47D - H1 abundance within AB compartments", xaxt="n",
         col = c("lightblue","sandybrown","lightblue","sandybrown"), cex.axis = 0.8,
         names = c("H1.2_A_wt","H1.2_B_wt","H1X_A_wt","H1X_B_wt"))
abline(h = -0.002, col = "lightblue", lwd = 1, lty = 1)
abline(h = 0, col = "gray", lwd = 1,lty = 1)
abline(h = 0.002, col = "forestgreen", lwd = 1, lty = 1)
abline(v = 2.75)
axis(1, c(1.5,4),c("H1.2","H1X") , cex.axis=0.8)
title(sub="Histone H1 (Wild-Type)", adj=0.5, line=2.5)
legend("topright",legend=c("A","B"), fill = c("lightblue","sandybrown"))
dev.off()

#--- BoxPlot representing the H1 abundance in A/B Compartments (WT+KD) -----------------------------------------------------------------------
pdf("/home/andrea/Bioinformatica_Núria/HiC_TADs/T47D_H1abundance_AB.pdf")
boxplot (AB_WT2$V5[AB_WT2$V4=="A"], AB_WT2$V5[AB_WT2$V4=="B"], 
         AB_KD2$V5[AB_WT2$V4=="A"], AB_KD2$V5[AB_WT2$V4=="B"],
         AB_WTx$V5[AB_WTx$V4=="A"], AB_WTx$V5[AB_WTx$V4=="B"],
         AB_KDx$V5[AB_WTx$V4=="A"], AB_KDx$V5[AB_WTx$V4=="B"], outline = F, at = c(1,2,3,4,5.5,6.5,7.5,8.5), 
         ylab = "ChIP-Seq Signal", main = "T47D - H1 abundance within AB compartments", xaxt="n",
         col = c("lightblue","sandybrown","lightblue","sandybrown","lightblue","sandybrown","lightblue","sandybrown"),
         names = c("H1.2_A_wt","H1.2_B_wt","H1.2_A_kd","H1.2_B_kd","H1X_A_wt","H1X_B_wt","H1X_A_kd","H1X_B_kd"))
abline(h = -0.002, col = "lightblue", lwd = 1, lty = 1)
abline(h = 0, col = "gray", lwd = 1,lty = 1)
abline(h = 0.002, col = "forestgreen", lwd = 1, lty = 1)
abline(v = 4.75)
axis(1, c(1.5,3.5,6,8),c("WT","KD","WT","KD") , cex.axis=0.75)
title(sub="H1.2", adj=0.24, line=2.5); title(sub="H1X", adj=0.76, line=2.5)
legend("topright",legend=c("A","B"), fill = c("lightblue","sandybrown"))
dev.off()


##############################################################################################################################################
# PTMs&CRs Abundance at A/B Compartments 
##############################################################################################################################################

#--- Setting Working Directory and required variables ----------------------------------------------------------------------------------------
setwd("/home/andrea/Bioinformatica_Núria/HiC_TADs/AB_Compartments/AB_vs_PTMs/")
myFiles <- list.files()
PTMs <- c("CHD1","CTCF","EZH2","H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9ac","H3K9me3","P300","RNAPOLII")  
ColorsCode <- list("CTCF"="lightblue1","H3K4me3"="lightgreen","H3K27ac"="aquamarine4", 
                   "H3K9ac"="darkolivegreen1","RNAPOLII"="royalblue","CHD1"="steelblue1", 
                   "H3K4me1"="cyan3","H3K36me3"="peachpuff","P300"="gold",
                   "H3K27me3"="coral2","EZH2"="sandybrown","H3K9me3"="brown")

TOTAL <- data.frame()

for (PTM in PTMs){
  #--- Reading the BED files with the PTMs abundance in A/B Compartments and normalizing by the length of the compartments -------------------
  FileName <- myFiles[grepl(as.character(PTM),myFiles)][1]
  File <- read.table(FileName, stringsAsFactors = F, sep = "\t")
  File$V5 <- (File$V5 / (File$V3-File$V2))*1000000
  #--- Separating A from B Compartments and including both in a global variable called "TOTAL" -----------------------------------------------
  FileA <- File[File$V4=="A",]; FileB <- File[File$V4=="B",]
  TOTAL <- rbind(TOTAL,FileA$V5,FileB$V5)
}

#--- Scaling "TOTAL" to compare the different PTMs and distributing the normalized values into independent variables -------------------------
TOTAL <- t(scale(TOTAL))
A_CHD1 <- TOTAL[,1]; B_CHD1 <- TOTAL[,2]; A_CTCF <- TOTAL[,3]; B_CTCF <- TOTAL[,4]; A_EZH2 <- TOTAL[,5]; 
B_EZH2 <- TOTAL[,6]; A_H3K27ac <- TOTAL[,7]; B_H3K27ac <- TOTAL[,8]; A_H3K27me3 <- TOTAL[,9]; B_H3K27me3 <- TOTAL[,10];
A_H3K36me3 <- TOTAL[,11]; B_H3K36me3 <- TOTAL[,12]; A_H3K4me1 <- TOTAL[,13]; B_H3K4me1 <- TOTAL[,14]; A_H3K4me3 <- TOTAL[,15];
B_H3K4me3 <- TOTAL[,16]; A_H3K9ac <- TOTAL[,17]; B_H3K9ac <- TOTAL[,18]; A_H3K9me3 <- TOTAL[,19]; B_H3K9me3 <- TOTAL[,20]; 
A_P300 <- TOTAL[,21]; B_P300 <- TOTAL[,22]; A_RNAPOLII <- TOTAL[,23]; B_RNAPOLII <- TOTAL[,24]

# Preparing loop's required variables
Plots_Names <- c("TFs_5","PTMs_5","PTMs_7")
TFs_5 <- c("RNAPOLII","CHD1","CTCF","EZH2","P300") 
PTMs_5 <- c("H3K27ac","H3K4me3","H3K36me3", "H3K27me3", "H3K9me3") 
PTMs_7 <- c("H3K4me1","H3K27ac","H3K4me3","H3K9ac","H3K36me3","H3K27me3","H3K9me3") 

# Opening a connection with a new PDF file
pdf("/home/andrea/Bioinformatica_Núria/HiC_TADs/AB_Compartments/AB_vs_PTMs/AB_vs_PTMs.pdf")
for (plot in Plots_Names){
  Plot_Data <- c()
  SpacesVector <- c()
  Count <- 1
  for (comp in c("A","B")){
    for (mark in eval(as.name(plot))){
      # Generating the data which will be used when performing each plot (for each band, for each mark: add their information to "Plot_Data")
      Plot_Data <- c(Plot_Data, list(eval(as.name(paste(comp,mark,sep = "_")))) )
      SpacesVector <- c(SpacesVector, Count)
      Count <- Count+1
    }
    Count <- Count+1
  }
  Colors <- c()
  for (PTM in eval(as.name(plot))){
    Colors <- c(Colors, ColorsCode[[PTM]])
  }
  boxplot(Plot_Data, col = Colors, outline = F, ylab = paste(plot," abundance",sep=""), xaxt="n", cex.lab=0.75,
          names = rep(eval(as.name(plot)),2), xaxt="n", cex = 0.75, cex.axis = 0.75, at = SpacesVector,
          main = paste(plot," abundance at AB Compartments (normalyzed to length, scaled)",sep=""), xlab=NULL)
  abline(h = 0, col = "grey", lwd = 1)
  abline(h = 0.5, col = "forestgreen", lwd = 1)
  abline(h = -0.5, col = "lightblue", lwd = 1)
  
  #-- Adding to X axis the name of the groups and the lines separating grous of boxes of PTMs -------------------------------------------
  VectorPositions <- c()
  BandPosition <- length(eval(as.name(plot)))-length(eval(as.name(plot)))/2.5
  for (band in 1:3){
    VectorPositions <- c(VectorPositions, BandPosition) 
    BandPosition <- BandPosition +1+ length(eval(as.name(plot)))
  }
  axis(1, VectorPositions[1:2], c("A Compartment","B Compartment"),cex.axis=0.75)
  abline(v=length(eval(as.name(plot)))+1 , col = "grey", lwd = 1)
  legend("topright", inset=.01, eval(as.name(plot)), fill = Colors, cex=0.7, bg ="white" )
}
dev.off()


################################################################################################################################################
# A/B Compartments vs CytoBands
################################################################################################################################################

#--- Setting Working Directory and required variables ------------------------------------------------------------------------------------------
WorkingDir <- "/home/andrea/Bioinformatica_Núria/HiC_TADs/AB_Compartments/AB_vs_CytoBands/"
setwd(WorkingDir)
CytoBands <- c("Gneg1","Gneg2","Gneg3","Gneg4","Gpos25","Gpos50","Gpos75","Gpos100")
myFiles <- list.files()

AB_bed <- read.table("/home/andrea/Bioinformatica_Núria/HiC_TADs/ABcompartments_01_WT.bed", stringsAsFactors = F, sep = "\t")
A_Compartments <- AB_bed[AB_bed$V4=="A",]; B_Compartments <- AB_bed[AB_bed$V4=="B",]
A_Compartments[,"Length"] <- A_Compartments$V3- A_Compartments$V2
B_Compartments[,"Length"] <- B_Compartments$V3- B_Compartments$V2
A_Total <- sum(A_Compartments$Length); B_Total <- sum(B_Compartments$Length)

#--- Calculating the Intersection between A/B Compartments and CytoBands -----------------------------------------------------------------------
for (band in CytoBands){
  #--- Reading the BED files with the intersection and calculating the total amount of overlapping base bairs ----------------------------------
  FileName <- myFiles[grepl(as.character(band),myFiles)][1]     
  File <- read.table(FileName, stringsAsFactors = F, sep = "\t")
  File$V5 <- (File$V3-File$V2)
  #--- Separating A from B Compartments, calculating the percentage of overlapping and declaring independent variables -------------------------
  FileA <- File[File$V4=="A",]; FileB <- File[File$V4=="B",]
  FileA$V5 <- FileA$V5/A_Total*100; FileB$V5 <- FileB$V5/B_Total*100
  do.call("<-", list(paste("A",band,sep="_"), FileA$V5))
  do.call("<-", list(paste("B",band,sep="_"), FileB$V5))
}

#--- Plotting a BoxPlot with the average overlapping between A/B Compartments and CytoBands in a PDF file --------------------------------------
pdf(paste(WorkingDir,"AB_CytoBands.pdf",sep = ""))
boxplot(A_Gneg1,A_Gneg2,A_Gneg3,A_Gneg4,A_Gpos25,A_Gpos50,A_Gpos75,A_Gpos100,
        B_Gneg1,B_Gneg2,B_Gneg3,B_Gneg4,B_Gpos25,B_Gpos50,B_Gpos75,B_Gpos100, 
        col = c("Gneg1"="#b3d9ff","Gneg2"="#99e699","Gneg3"="#ffccb3","Gneg4"="#e59acc",
                "Gpos25"="#66b3ff","Gpos50"="#47d147","Gpos75"="#ffaa80","Gpos100"="#d147a3"), outline = F,
        at=c(1,2,3,4,5.5,6.5,7.5,8.5,10,11,12,13,14.5,15.5,16.5,17.5), ylab = "Overlapping (%)", xaxt = "n",
        main = "T47D - Overlapping between CytoBands & AB Compartments", cex.axis=0.9)
abline(v=9.25, col="gray",lwd=1)
axis(1, c(4.75,13.75), c("A Compartment","B Compartment"),cex.axis=0.9)
legend("topright",c("Gneg1","Gneg2","Gneg3","Gneg4","Gpos25","Gpos50","Gpos75","Gpos100"), bg="white",
       fill = c("#b3d9ff","#99e699","#ffccb3","#e59acc","#66b3ff","#47d147","#ffaa80","#d147a3"), cex = 0.85)
dev.off()


################################################################################################################################################
# A/B Compartments vs TADs Groups
################################################################################################################################################

#--- Setting Working Directory and required variables ------------------------------------------------------------------------------------------
WorkingDir <- "/home/andrea/Bioinformatica_Núria/HiC_TADs/AB_Compartments/AB_vs_TADsGroups/"
setwd(WorkingDir)
TADsGroups <- c("H1.2_enriched","LIKE_H1.2","H1X_enriched","LIKE_H1X")
myFiles <- list.files()

#--- Calculating the Intersection between A/B Compartments and TADs Groups ---------------------------------------------------------------------
AB_TADs <- c(); TOTAL <- data.frame()
for (group in TADsGroups){
  FileName <- myFiles[grepl(as.character(group),myFiles)][1]
  File <- read.table(FileName, stringsAsFactors = F, sep = "\t")
  File$V5 <- (File$V3-File$V2)
  FileA <- File[File$V4=="A",]; FileB <- File[File$V4=="B",]
  FileA$V5 <- FileA$V5/A_Total*100 ; FileB$V5 <- FileB$V5/B_Total*100
  do.call("<-", list(paste("A",group,sep="_"), FileA$V5))
  do.call("<-", list(paste("B",group,sep="_"), FileB$V5))
  AB_TADs <- c(AB_TADs, paste("A",group,sep="_"), paste("B",group,sep="_"))
  TOTAL <- rbind(TOTAL,FileA$V5,FileB$V5)
}

#--- Plotting a BoxPlot with the average overlapping between A/B Compartments and TADs Groups in a PDF file ------------------------------------
pdf(paste(WorkingDir,"AB_TADsGroups.pdf",sep = ""), height = 9, width = 8)
par(mar=c(6,4,3,2))
boxplot(A_H1.2_enriched,A_LIKE_H1.2,A_LIKE_H1X,A_H1X_enriched,
        B_H1.2_enriched,B_LIKE_H1.2,B_LIKE_H1X,B_H1X_enriched,
        col= c("lightblue2","cyan3","sandybrown","coral2"), cex = 0.9,
        main = "T47D - Overlapping between TADs Groups & AB Compartments", cex.axis=0.9,
        outline = F, at=c(1,2,3,4,5.5,6.5,7.5,8.5), ylab = "Overlapping (%)", xaxt = "n")
abline(v=4.75, col="gray",lwd=1)
axis(1, c(2.5,7), c("A Compartment","B Compartment"),cex.axis=0.9)
legend("topright",c("H1.2","H1.2_Like","H1X_Like","H1X"),fill = c("lightblue2","cyan3","sandybrown","coral2"), cex = 0.85)
dev.off()


################################################################################################################################################
# A/B Compartments: Changing Bins
################################################################################################################################################

#--- Setting Working Directory and required variables ------------------------------------------------------------------------------------------
setwd("/home/andrea/Bioinformatica_Núria/HiC_TADs/TADs/TADs_Groups/TADsGroups_H1abundance_X2Ratio/")
H12_TADs <- read.table("H1.2_enriched.bed",sep="\t",stringsAsFactors=F) 
H12like_TADs <- read.table("LIKE_H1.2_enriched.bed",sep="\t",stringsAsFactors=F)
H1Xlike_TADs <- read.table("LIKE_H1X_enriched.bed",sep="\t",stringsAsFactors=F)
H1X_TADs <- read.table("H1X_enriched.bed",sep="\t",stringsAsFactors=F)

WorkingDir <- "/home/andrea/Bioinformatica_Núria/HiC_TADs/BINS/"; setwd(WorkingDir)
H12_Bins_AtoB <- read.table("H1.2endo_r2_hg19.bdg_changeA_B.bed",sep="\t",stringsAsFactors=F)
H12_Bins_BtoA <- read.table("H1.2endo_r2_hg19.bdg_changeB_A.bed",sep="\t",stringsAsFactors=F)
H1X_Bins_AtoB <- read.table("H1Xendo_r3_hg19.bdg_changeA_B.bed",sep="\t",stringsAsFactors=F)
H1X_Bins_BtoA <- read.table("H1Xendo_r3_hg19.bdg_changeB_A.bed",sep="\t",stringsAsFactors=F)

#--- Ratio 2/X Comparison: TADsGroups vs Changing A/B Bins -------------------------------------------------------------------------------------
H12_Bins_AtoB <- H12_Bins_AtoB[H12_Bins_AtoB$V1!="chr9",]; H12_Bins_AtoB <- H12_Bins_BtoA[H12_Bins_BtoA$V1!="chr9",] # Removing Chr9
Ratio_AtoB <- (H12_Bins_AtoB$V4+10)/(H1X_Bins_AtoB$V4+10)
Ratio_BtoA <- (H12_Bins_BtoA$V4+10)/(H1X_Bins_BtoA$V4+10)

#--- Plotting a BoxPlot to compare the ratio 2/X within TADs Groups vs changins A/B & B/A bins in a PDF file -----------------------------------
pdf(paste(WorkingDir,"2XRatio_TADsGroups_ChangingBinsAB.pdf",sep=""),width = 8)
boxplot(list(H12_TADs$V6, H12like_TADs$V6, H1Xlike_TADs$V6, H1X_TADs$V6, Ratio_AtoB, Ratio_BtoA), outline = F, cex.main=1.5, cex.lab=1.3, cex.axis=1.2,
        col = c("lightblue","cyan3","coral2","sandybrown","#4dd6e6","#4dd6e6"), at=c(1,2,3,4,5.5,6.5),
        names=c("H1.2","H1.2-like","H1X-like","H1X","A-B","B-A"), xlab="       TADs Groups                                        Bins",
        main="TADs Groups vs Changing A/B Bins", ylab="2/X Ratio")
abline(v=4.75)
dev.off()

#--- H1 Abundance Comparison: H1.2, H1.4 & H1X within A/B Bins ---------------------------------------------------------------------------------
HistonesAnalysis <- c("H1.2","H1.4","H1X"); Bins <- c("A_B","B_A")
myFiles <- list.files(WorkingDir, pattern="*.bed")
HistonesBins <- c()
for (histone in HistonesAnalysis){
  for (bin in Bins){
    #--- Iterating over H1s and changing bins, loading the corresponding files and storing their content in independent variables --------------
    FileName <- myFiles [grepl(histone,myFiles) & grepl(bin,myFiles)][1]
    print(FileName)
    Histone_Bin <- read.table (FileName, header=FALSE, stringsAsFactors = F)
    do.call("<-",list(paste(histone,bin, sep="_"), na.omit(Histone_Bin)))
    HistonesBins <- c(HistonesBins, paste(histone,bin, sep="_"))
  }
}

#--- Plotting a BoxPlot representing H1 abundance within changing A/B & B/A bins in a PDF file -------------------------------------------------
pdf(paste(WorkingDir,"T47D_H1Abundance_ChangingBinsAB.pdf",sep=""),width = 7.5)
boxplot(list(H1.2_A_B$V4,H1.2_B_A$V4,H1.4_A_B$V4,H1.4_B_A$V4,H1X_A_B$V4,H1X_B_A$V4), at=c(1,2,3.5,4.5,6,7),
        outline = F, cex.main=1.5, cex.lab=1.3, cex.axis=1.2, col = c("lightblue","sandybrown"), 
        xlab="H1.2                          H1.4                          H1X", names=c("A-B","B-A","A-B","B-A","A-B","B-A"),
        main="H1 Abundance in Changing A/B Bins", ylab="ChIP-Seq Signal") 
abline(v=c(2.75,5.25))
abline(h=0,col="grey")
boxplot(list(H1.2_A_B$V4,H1.2_B_A$V4,H1.4_A_B$V4,H1.4_B_A$V4,H1X_A_B$V4,H1X_B_A$V4),at=c(1,2,3.5,4.5,6,7),
        outline = F, col = c("lightblue","sandybrown"), add=TRUE, xaxt="n", yaxt="n") 
dev.off()
################################################################################################################################################
