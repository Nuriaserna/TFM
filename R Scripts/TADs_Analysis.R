##############################################################################################################################################
#                                                              TADs ANALYSIS
##############################################################################################################################################
# This R script has been designed to analyze TADs obtained from a Hi-C experiment.It is composed of four main sections: 
# H1 Abundance in A/B Compartments, PTMs&CRs Abundance in A/B Compartments, A/B Compartments vs CytoBands and A/B Compartments vs TADs Groups



########################################################################################################################################################
#                                             Calculating the average number of TADs edges within each cytoband type
########################################################################################################################################################
# Opening the files containing the Gneg&Gpos bands information and the TADs coordinates
for (num in c(25,50,75,100)){
  do.call("<-", list(paste("Gpos",num,sep=""),
                     read.table(paste("/home/andrea/Bioinformatica_Núria/Cytobands/CytoBands_new/cytoBand_Gpos",num,".bed",sep=""),
                                header = F, sep="\t", stringsAsFactors = F)
  ))
}
for (num in c(1,2,3,4)){
  do.call("<-", list(paste("Gneg",num,sep=""),
                     read.table(paste("/home/andrea/Bioinformatica_Núria/Cytobands/CytoBands_new/Gneg",num,".bed",sep=""),
                                header = F, sep="\t", stringsAsFactors = F)
  ))
}
TADs <- read.table("/home/andrea/Bioinformatica_Núria/HiC_TADs/TADs_01_WT.bed",header=F,sep="\t")

# Modifying the CytoBands' variables to just use the columns with the coordinates
Gpos25 <- Gpos25[,1:3]; Gpos50 <- Gpos50[,1:3]; Gpos75 <- Gpos75[,1:3]; Gpos100 <- Gpos100[,1:3]
Gneg1 <- Gneg1[,1:3]; Gneg2 <- Gneg2[,1:3]; Gneg3 <- Gneg3[,1:3]; Gneg4 <- Gneg4[,1:3] 

########################################################################################################################################################
# For each CytoBand:
Total_Counts <- c()
for (CytoBand in c("Gneg1","Gneg2","Gneg3","Gneg4","Gpos25","Gpos50","Gpos75","Gpos100")){
  Bands_Counts <- c()
  print(CytoBand)
  # For each band within a specific CytoBand:
  for (band in 1:length(rownames(eval(as.name(CytoBand)))) ){
    Band_Count <- 0
    
    # For each TAD:
    for (TAD in 1:length(rownames(TADs)) ){
      # If both the evaluated band and TAD are located in the same chromosome and one of the edges is included in the band, add 1 to the "Band_Count" counter
      if (TADs[TAD,1] == eval(as.name(CytoBand))[band,1] & TADs[TAD,2] > eval(as.name(CytoBand))[band,2] & TADs[TAD,2] < eval(as.name(CytoBand))[band,3]){
        Band_Count <- Band_Count + 1
      }
      if (TAD == length(rownames(TADs)) & TADs[TAD,1] == eval(as.name(CytoBand))[band,1] & TADs[TAD,3] > eval(as.name(CytoBand))[band,2]){
        Band_Count <- Band_Count + 1
      }
    }
    # If the specific band has at least a count, add its value to "Bands_Count"
    if (Band_Count > 0){
      Bands_Counts <- c(Bands_Counts, Band_Count/(eval(as.name(CytoBand))[band,3]-eval(as.name(CytoBand))[band,2])  )
    }
  }
  Bands_Counts <- Bands_Counts*1000000
  Total_Counts <- c(Total_Counts, list(Bands_Counts)) # Add all the bands' counts of a specific CytoBand group to the variable "Total_Counts"
}
save(Total_Counts, file = "/home/andrea/Bioinformatica_Núria/HiC_TADs//Total_TADedges_Counts.RData")
########################################################################################################################################################
# Opening a connection with a PDF output file and printing the BoxPlot of interest
pdf("/home/andrea/Bioinformatica_Núria/HiC_TADs/Average_TADs_edges_CytoBands_no0.pdf")
boxplot(Total_Counts, cex = 0.85, cex.axis = 0.85, outline =FALSE,
        ylab = "TADs edges", main = "Average number of TADs edges per cytoband type",
        names = c("Gneg1","Gneg2","Gneg3","Gneg4","Gpos25","Gpos50","Gpos75","Gpos100"),
        col = c("Gneg1"="#b3d9ff","Gneg2"="#99e699","Gneg3"="#ffccb3","Gneg4"="#e59acc",
                "Gpos25"="#66b3ff","Gpos50"="#47d147","Gpos75"="#ffaa80","Gpos100"="#d147a3"))
dev.off()
########################################################################################################################################################
########################################################################################################################################################







#--------------------------------------------TADs Groups-------------------------------------------------------------------------------------
setwd("/home/andrea/Bioinformatica_Núria/HiC_TADs/H1abundance_TADs/")
H1.2_TADs <- read.table("H1.2endo_r2_hg19.bdg_TADs_01_WT.bed", header = F, sep = "\t")
H1X_TADs <- read.table("H1Xendo_r3_hg19.bdg_TADs_01_WT.bed", header = F, sep = "\t")
######################################################################################################################################################
# 1b| CONSTRUCTING TWO GROUPS OF TADs ACCORDING TO THEIR H1.2b OR H1X ENRICHMENT (RATIO)
# Merging H1.2 and H1X data and defining two different groups according to their H1.2 and H1X content
H1variants_TADs <- merge(H1.2_TADs, H1X_TADs, by="row.names")
head(H1variants_TADs); dim(H1variants_TADs)

rownames(H1variants_TADs) <- H1variants_TADs$Row.names
H1variants_TADs <- H1variants_TADs[,c(2,3,4,5,9)]
colnames(H1variants_TADs) <- c("Chr","From","To","H1.2","H1X")

H1variants_TADs_Ratio <- c()

for (TAD in 1:length(rownames(H1variants_TADs))){
  Ratio <- (H1variants_TADs[TAD,"H1.2"]+10)/(H1variants_TADs[TAD,"H1X"]+10)
  print(Ratio)
  H1variants_TADs_Ratio <- rbind(H1variants_TADs_Ratio, data.frame(H1variants_TADs[TAD,],Ratio))
}

Sorted <- H1variants_TADs_Ratio[ match(sort(H1variants_TADs_Ratio[,"Ratio"],decreasing = TRUE), H1variants_TADs_Ratio[,"Ratio"]),]
length(rownames(Sorted))/4
Group_H12 <- Sorted[1:802,]
Group_H12like <- Sorted[803:1604,]
Group_H1Xlike <- Sorted[1605:2406,]
Group_H1X <- Sorted[2407:3209,]

# BoxPlot showing the average length of both H1.2 and H1X enriched TADs
pdf(file="/home/andrea/Bioinformatica_Núria/HiC_TADs/H1abundance_TADs/H12_H1X_Ratios_enrichedTADs_length.pdf", width = 8, height = 9)
Data <- c(c(Group_H12[3]-Group_H12[2]), c(Group_H12like[3]-Group_H12like[2]), 
          c(Group_H1Xlike[3]-Group_H1Xlike[2]), c(Group_H1X[3]-Group_H1X[2]))
boxplot(Data, outline = F, names = c("H1.2 enriched","H1.2 enriched-like","H1X enriched-like","H1X enriched"), main="Average TADs length",
        ylab = "TADs length (bp)", col=c("lightblue2","cyan3","coral2","sandybrown"), cex=0.9, cex.axis=0.9)
dev.off()



Scatter_TADs <- H1variants_TADs[, c(4,5)]  # Remove the first column, which just stored the rownames
colnames(Scatter_TADs) <- c("H1.2","H1X")
head(Scatter_TADs)

a <- cor.test(Scatter_TADs$H1X, Scatter_TADs$H1.2, method= "pearson")
Correlation_Table <- data.frame(x= c(a$p.value),
                                y= c(a$estimate),
                                z= c("H1X") )
Correlation_Table <- as.data.frame(t(Correlation_Table))
row.names(Correlation_Table) <- c("p.value", "corr.coef", "VS.Variant")
colnames(Correlation_Table)  <- c("H1.2")

TADs_H1.2_H1X <- ggplot(Scatter_TADs, aes(x=H1X, y=H1.2)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour="ligthblue"),  size=0.2, show.legend = F) + theme_classic() + 
  theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
                          plot.title = element_text(hjust = 0.5, color = "#666666"),plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
  ggtitle(label = "H1.2 vs H1X (TADs)", subtitle = paste("p.value = ",format(a$p.value,digits=3),"\ncorr.coef = ",round(a$estimate,3),sep="")) + scale_color_manual(values="steelblue1")
TADs_H1.2_H1X



pdf("/home/andrea/Bioinformatica_Núria/HiC_TADs/Rplot.pdf")
library(ggplot2)
ggplot()+
  geom_point(data=Group_H12, aes(x=H1X, y=H1.2, color="H1.2"), size=0.2, show.legend = F)+
  geom_point(data=Group_H12like, aes(x=H1X, y=H1.2, color="H1.2_Like"), size=0.2, show.legend = F)+
  geom_point(data=Group_H1Xlike, aes(x=H1X, y=H1.2, color="H1X_Like"), size=0.2, show.legend = F)+
  geom_point(data=Group_H1X, aes(x=H1X, y=H1.2, color="H1X"), size=0.2, show.legend = F) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5, color = "#666666"),
                          plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
  scale_color_manual(values=c("lightblue2","cyan3","sandybrown","coral2"))+
  ggtitle(label = "H1.2 vs H1X (TADs)", subtitle = "p.value = 6.94e-54\ncorr.coef = -0.268")
dev.off()

Group_H12 <- Group_H12[,c(1,2,3)]; Group_H12 <- Group_H12[ order( as.numeric(rownames(Group_H12))) ,] 
Group_H12[,"From"] <- as.integer(Group_H12[,"From"]); Group_H12[,"To"] <- as.integer(Group_H12[,"To"]) 

Group_H12like <- Group_H12like[,c(1,2,3)]; Group_H12like <- Group_H12like[ order( as.numeric(rownames(Group_H12like))) ,] 
Group_H12like[,"From"] <- as.integer(Group_H12like[,"From"]); Group_H12like[,"To"] <- as.integer(Group_H12like[,"To"]) 

Group_H1Xlike <- Group_H1Xlike[,c(1,2,3)]; Group_H1Xlike <- Group_H1Xlike[ order( as.numeric(rownames(Group_H1Xlike))) ,] 
Group_H1Xlike[,"From"] <- as.integer(Group_H1Xlike[,"From"]); Group_H1Xlike[,"To"] <- as.integer(Group_H1Xlike[,"To"]) 

Group_H1X <- Group_H1X[,c(1,2,3)]; Group_H1X <- Group_H1X[ order( as.numeric(rownames(Group_H1X))) ,] 
Group_H1X[,"From"] <- as.integer(Group_H1X[,"From"]); Group_H1X[,"To"] <- as.integer(Group_H1X[,"To"]) 

write.table(Group_H12, file = "/home/andrea/Bioinformatica_Núria/HiC_TADs/TADs/TADs_Groups/H1.2_enriched_TADs.bed", quote=F,col.names = F,row.names = F,sep = "\t")
write.table(Group_H1X, file = "/home/andrea/Bioinformatica_Núria/HiC_TADs/TADs/TADs_Groups/H1X_enriched_TADs.bed", quote=F,col.names = F,row.names = F,sep = "\t")
write.table(Group_H12like, file = "/home/andrea/Bioinformatica_Núria/HiC_TADs/TADs/TADs_Groups/H1.2_enrichedLIKE_TADs.bed", quote=F,col.names = F,row.names = F,sep = "\t")
write.table(Group_H1Xlike, file = "/home/andrea/Bioinformatica_Núria/HiC_TADs/TADs/TADs_Groups/H1X_enrichedLIKE_TADs.bed", quote=F,col.names = F,row.names = F,sep = "\t")
#---------------------------------------------------------------------------------------------------------------------------------


#########################################################################################################################
#                                         TADs vs. CytoBands
#########################################################################################################################
#-------------------------------------------------------------------------------------------
WorkingDir <- "/home/andrea/Bioinformatica_Núria/HiC_TADs/TADs/TADs_vs_CytoBands/"
setwd(WorkingDir)
CytoBands <- c("Gneg1","Gneg2","Gneg3","Gneg4","Gpos25","Gpos50","Gpos75","Gpos100")
TADsGroups <-  c("H1.2_enriched","LIKE_H1.2_enriched","H1X_enriched","LIKE_H1X_enriched")
myFiles <- list.files()

TADs_CytoBands <- c(); TOTAL <- data.frame()
for (group in TADsGroups){
  TADs_bed <- read.table(paste("/home/andrea/Bioinformatica_Núria/HiC_TADs/TADs/TADs_Groups/",group,".bed",sep=""), stringsAsFactors = F, sep = "\t")
  TADs_bed[,"Length"] <- as.numeric(TADs_bed$V3 - TADs_bed$V2)
  TADs_Total <- sum(TADs_bed$Length)
  for (band in CytoBands){
    FileName <- myFiles[grepl(as.character(group),myFiles) & grepl(as.character(band),myFiles)][1]
    File <- read.table(FileName, stringsAsFactors = F, sep = "\t")
    File$V5 <- (File$V3-File$V2)
    File$V5 <- File$V5
    do.call("<-", list(paste("New",group,band,sep="_"), File$V5))
    TADs_CytoBands <- c(TADs_CytoBands, paste("New",group,band,sep="_"))
    TOTAL <- rbind(TOTAL,File$V5)
  }
}

pdf(paste(WorkingDir,"TADsGroups_CytoBands.pdf",sep = ""))
boxplot(New_H1X_enriched_Gneg1,New_H1X_enriched_Gneg2,New_H1X_enriched_Gneg3,New_H1X_enriched_Gneg4,
        New_H1X_enriched_Gpos25,New_H1X_enriched_Gpos50,New_H1X_enriched_Gpos75,New_H1X_enriched_Gpos100,
        New_LIKE_H1X_enriched_Gneg1,New_LIKE_H1X_enriched_Gneg2,New_LIKE_H1X_enriched_Gneg3,New_LIKE_H1X_enriched_Gneg4,
        New_LIKE_H1X_enriched_Gpos25,New_LIKE_H1X_enriched_Gpos50,New_LIKE_H1X_enriched_Gpos75,New_LIKE_H1X_enriched_Gpos100,
        
        New_LIKE_H1.2_enriched_Gneg1,New_LIKE_H1.2_enriched_Gneg2,New_LIKE_H1.2_enriched_Gneg3,New_LIKE_H1.2_enriched_Gneg4,
        New_LIKE_H1.2_enriched_Gpos25,New_LIKE_H1.2_enriched_Gpos50,New_LIKE_H1.2_enriched_Gpos75,New_LIKE_H1.2_enriched_Gpos100,
        New_H1.2_enriched_Gneg1,New_H1.2_enriched_Gneg2,New_H1.2_enriched_Gneg3,New_H1.2_enriched_Gneg4,
        New_H1.2_enriched_Gpos25,New_H1.2_enriched_Gpos50,New_H1.2_enriched_Gpos75,New_H1.2_enriched_Gpos100,
        
        col = c("Gneg1"="#b3d9ff","Gneg2"="#99e699","Gneg3"="#ffccb3","Gneg4"="#e59acc",
                "Gpos25"="#66b3ff","Gpos50"="#47d147","Gpos75"="#ffaa80","Gpos100"="#d147a3"), outline = F,
        at=c(1,2,3,4,5,6,7,8,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,
             18,19,20,21,22,23,24,25,26.5,27.5,28.5,29.5,30.5,31.5,32.5,33.5), ylab = "Overlapping (%)", xaxt = "n",
        main = "T47D - Overlapping between CytoBands & TADs Groups", cex.axis=0.9)
abline(v=c(8.75,17.25,25.75), col="gray",lwd=1)
axis(1, c(4.5,13,21.5,30), c("H1X enriched","H1X_Like enriched","H1.2_Like enriched","H1.2 enriched"),cex.axis=0.8)

legend("topright",c("Gneg1","Gneg2","Gneg3","Gneg4","Gpos25","Gpos50","Gpos75","Gpos100"), bg="white",
       fill = c("#b3d9ff","#99e699","#ffccb3","#e59acc","#66b3ff","#47d147","#ffaa80","#d147a3"), cex = 0.85)

dev.off()


#--- Per Chromosome -------------------------------------------------------------------------------------------
Chromosomes <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                 "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")

TADs_CytoBands_Sorted_CyBa <- c( "New_H1X_enriched_Gneg1","New_H1X_enriched_Gneg2","New_H1X_enriched_Gneg3","New_H1X_enriched_Gneg4",
                                 "New_H1X_enriched_Gpos25","New_H1X_enriched_Gpos50","New_H1X_enriched_Gpos75","New_H1X_enriched_Gpos100",
                                 "New_LIKE_H1X_enriched_Gneg1","New_LIKE_H1X_enriched_Gneg2","New_LIKE_H1X_enriched_Gneg3","New_LIKE_H1X_enriched_Gneg4",
                                 "New_LIKE_H1X_enriched_Gpos25","New_LIKE_H1X_enriched_Gpos50","New_LIKE_H1X_enriched_Gpos75","New_LIKE_H1X_enriched_Gpos100",
                                 "New_LIKE_H1.2_enriched_Gneg1","New_LIKE_H1.2_enriched_Gneg2","New_LIKE_H1.2_enriched_Gneg3","New_LIKE_H1.2_enriched_Gneg4",  
                                 "New_LIKE_H1.2_enriched_Gpos25","New_LIKE_H1.2_enriched_Gpos50","New_LIKE_H1.2_enriched_Gpos75","New_LIKE_H1.2_enriched_Gpos100",
                                 "New_H1.2_enriched_Gneg1","New_H1.2_enriched_Gneg2","New_H1.2_enriched_Gneg3","New_H1.2_enriched_Gneg4",
                                 "New_H1.2_enriched_Gpos25","New_H1.2_enriched_Gpos50","New_H1.2_enriched_Gpos75","New_H1.2_enriched_Gpos100")

PerChromosome <- data.frame(row.names = Chromosomes)
for (var in TADs_CytoBands_Sorted_CyBa){
  PerChromosome[,var] <- rep(0,23)
  for (chr in Chromosomes){
    Chr_Count <- 0
    Chr_Count <- Chr_Count + sum(eval(as.name(var))[eval(as.name(var))$V1==chr, "V5"])
    PerChromosome[chr,var] <- Chr_Count}
}

PerChromosome_No0 <- c()
for (col in 1:length(colnames(PerChromosome))){
  PerChromosome_No0 <- c(PerChromosome_No0, list(PerChromosome[,col][PerChromosome[,col]>0]) )
}

pdf(paste(WorkingDir,"TADsGroups_CytoBands_CytoBands.pdf",sep = ""),height = 8, width = 10)
par(mar=c(6,4,3,6))
boxplot(PerChromosome_No0,
        col = c("Gneg1"="#b3d9ff","Gneg2"="#99e699","Gneg3"="#ffccb3","Gneg4"="#e59acc",
                "Gpos25"="#66b3ff","Gpos50"="#47d147","Gpos75"="#ffaa80","Gpos100"="#d147a3"), outline = F,
        at=c(1,2,3,4,5,6,7,8,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,
             18,19,20,21,22,23,24,25,26.5,27.5,28.5,29.5,30.5,31.5,32.5,33.5), ylab = "Overlapping (bp)",
        main = "T47D - Overlapping between CytoBands & TADs Groups", cex.axis=0.7, xaxt="n",ylim = c(0, 36000000))
abline(v=c(8.75,17.25,25.75), col="gray",lwd=1)
axis(1, c(4.5,13,21.5,30), c("H1X enriched","H1X_Like enriched","H1.2_Like enriched","H1.2 enriched"),cex.axis=0.8)
text(y=-7000000, x=18, labels="TADs Group", cex = 1, xpd=TRUE)
legend("topleft",c("Gneg1","Gneg2","Gneg3","Gneg4","Gpos25","Gpos50","Gpos75","Gpos100"), bg="white",
       fill = c("#b3d9ff","#99e699","#ffccb3","#e59acc","#66b3ff","#47d147","#ffaa80","#d147a3"), cex = 0.85)
dev.off()

#---- separated graphs ----
pdf(paste(WorkingDir,"TADs_Groups_CytoBands.pdf",sep = ""), width = 8.5)
par(mar=c(6,5,3,6))
boxplot(PerChromosome_No0[1:8],
        col = c("Gneg1"="#b3d9ff","Gneg2"="#99e699","Gneg3"="#ffccb3","Gneg4"="#e59acc",
                "Gpos25"="#66b3ff","Gpos50"="#47d147","Gpos75"="#ffaa80","Gpos100"="#d147a3"), outline = F,
        at=c(1,2,3,4,5.5,6.5,7.5,8.5), cex.axis=1.6, xaxt="n")
axis(1, c(1,2,3,4,5.5,6.5,7.5,8.5), c("1","2","3","4","25","50","75","100"),cex.axis=1.8)
text(y=-5400000, x=c(2.5,7), labels=c("Gneg","Gpos"), cex = 1.8, xpd=TRUE)
title(main=list("H1X enriched",cex=2), ylab = list("Overlapping (bp)",cex=1.7))
#legend("topright",c("Gneg1","Gneg2","Gneg3","Gneg4","Gpos25","Gpos50","Gpos75","Gpos100"), bg="white",
#       fill = c("#b3d9ff","#99e699","#ffccb3","#e59acc","#66b3ff","#47d147","#ffaa80","#d147a3"), cex = 1.5)

par(mar=c(6,5,3,6))
boxplot(PerChromosome_No0[9:16],
        col = c("Gneg1"="#b3d9ff","Gneg2"="#99e699","Gneg3"="#ffccb3","Gneg4"="#e59acc",
                "Gpos25"="#66b3ff","Gpos50"="#47d147","Gpos75"="#ffaa80","Gpos100"="#d147a3"), outline = F,
        at=c(1,2,3,4,5.5,6.5,7.5,8.5), cex.axis=1.5, xaxt="n")
title(main=list("H1X enriched-like",cex=2), ylab = list("Overlapping (bp)",cex=1.7))
axis(1, c(1,2,3,4,5.5,6.5,7.5,8.5), c("1","2","3","4","25","50","75","100"),cex.axis=1.8)
text(y=-2700000, x=c(2.5,7), labels=c("Gneg","Gpos"), cex = 1.8, xpd=TRUE)

par(mar=c(6,5,3,6))
boxplot(PerChromosome_No0[17:24],
        col = c("Gneg1"="#b3d9ff","Gneg2"="#99e699","Gneg3"="#ffccb3","Gneg4"="#e59acc",
                "Gpos25"="#66b3ff","Gpos50"="#47d147","Gpos75"="#ffaa80","Gpos100"="#d147a3"), outline = F,
        at=c(1,2,3,4,5.5,6.5,7.5,8.5), cex.axis=1.5, xaxt="n")
title(main=list("H1.2 enriched-like",cex=2), ylab = list("Overlapping (bp)",cex=1.7))
axis(1, c(1,2,3,4,5.5,6.5,7.5,8.5), c("1","2","3","4","25","50","75","100"),cex.axis=1.8)
text(y=-3200000, x=c(2.5,7), labels=c("Gneg","Gpos"), cex = 1.8, xpd=TRUE)

par(mar=c(6,5,3,6))
boxplot(PerChromosome_No0[25:32],
        col = c("Gneg1"="#b3d9ff","Gneg2"="#99e699","Gneg3"="#ffccb3","Gneg4"="#e59acc",
                "Gpos25"="#66b3ff","Gpos50"="#47d147","Gpos75"="#ffaa80","Gpos100"="#d147a3"), outline = F,
        at=c(1,2,3,4,5.5,6.5,7.5,8.5), cex.axis=1.5,cex=1.5, xaxt="n")
title(main=list("H1.2 enriched",cex=2), ylab = list("Overlapping (bp)",cex=1.7))
axis(1, c(1,2,3,4,5.5,6.5,7.5,8.5), c("1","2","3","4","25","50","75","100"),cex.axis=1.8)
text(y=-7800000, x=c(2.5,7), labels=c("Gneg","Gpos"), cex = 1.8, xpd=TRUE)
dev.off()






###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

# FIRST OF ALL, WE DETERMINED THE GROUPS OF TADs, THEIR RELATION WITH CYTOBANDS etc (BECAUSE WE DO NOT HAVE INFO ABOUT KD)
# THEN WE INCLUDE THE DOX EXPERIMENT INFO:

###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###









#--------------------------------------------WT_KD_TADs_Pre/Post DOX-------------------------------------------------------------------------------------
# Loading the files of interest
TADs_WT <- read.table("/home/andrea/Bioinformatica_Núria/HiC_TADs/Original_Files/TADs_01_WT.bed")
TADs_KD <- read.table("/home/andrea/Bioinformatica_Núria/HiC_TADs/Original_Files/TADs_02_multiH1-KD.bed")

# Preparing the required variables
colnames(TADs_WT) <- c("Chr","From","To","Str","Int"); colnames(TADs_KD) <- c("Chr","From","To","Str","Int")
Groups_N_WT <- data.frame(rep(0,10)); Groups_N_KD <- data.frame(rep(0,10))
colnames(Groups_N_WT) <- "Times"; colnames(Groups_N_KD) <- "Times" 

#-----------------------------------------------------------------------------------------------------------
# Counting 
for (TAD in 1:length(rownames(TADs_WT)) ){
  for (N in 1:10){
    if (TADs_WT[TAD,"Str"]==N){Groups_N_WT[as.character(N),"Times"] <- Groups_N_WT[as.character(N),"Times"]+1}
  }
}

for (TAD in 1:length(rownames(TADs_KD)) ){
  for (N in 1:10){
    if (TADs_KD[TAD,"Str"]==N){Groups_N_KD[as.character(N),"Times"] <- Groups_N_KD[as.character(N),"Times"]+1}
  }
}

Groups_Times_WT <- Groups_N_WT$Times
Groups_Times_KD <- Groups_N_KD$Times

#--------------------------------------------------------------------------------------------------
# Representing the frequency of TADs with a specific border strength with a BarPlot

pdf("/home/andrea/Bioinformatica_Núria/HiC_TADs/T47D_TADs_BorderStrength_Freq.pdf")
Wt <- barplot(Groups_Times_WT, main = "T47D - TADs' Groups by Border Strength", cex.axis = 1.2, names = 1:10, col = "lightblue")
title(xlab = "TADs' Border Strength (Wild-Type)",cex.lab=1.2, line = 3.5)
text(x=Wt, y =Groups_Times_WT+80, labels=as.character(Groups_N_WT$Times), cex = 1, xpd=TRUE)
title(ylab ="TADs Borders",cex.lab=1.2, line = 3)

Kd <- barplot(Groups_Times_KD, main = "T47D - TADs' Groups by Border Strength", cex.axis = 1.2, names = 1:10, col = "lightblue")
text(x=Kd, y =Groups_Times_KD+80, labels=as.character(Groups_N_KD$Times), cex = 1, xpd=TRUE)
title(xlab = "TADs' Border Strength (Knock-Down)", cex.lab=1.2, line = 3.5)
title(ylab ="TADs Borders",cex.lab=1.2, line = 3)

B <- barplot(Groups_Times_KD-Groups_Times_WT, main = "T47D - TADs Border Strength Re-Distribution after Dox Treatment", cex.axis = 1.2,
             col = c("lightblue","lightblue","indianred","indianred","indianred","lightblue","indianred","indianred","lightblue","lightblue"))
text(x=B, y =c(5,5,-11,-16,-44,5,-33,-13,9,102), labels=as.character(Groups_Times_KD-Groups_Times_WT), cex = 1, xpd=TRUE)
text(x=B, y=-53, labels = 1:10, xpd = TRUE, cex = 1)
title(xlab="TADs' Border Strength (Subtraction)",cex.lab=1.2, line = 3.8)
title(ylab ="TADs Borders",cex.lab=1.2, line = 3)

dev.off()
##################################################################################################################









###############################################################################################################################################################
#                                           TADs & AB Compartments Analysis
#######################################################################################################################
# 1| TADs Comparison: AJV vs François, Dox+ vs Dox-, A vs B Compartments, High vs Low 2/X 
library(VennDiagram)
library(gridExtra)
WorkingDir <- "/home/andrea/Bioinformatica_Núria/HiC_TADs/TADs/FRA_vs_AJV/"
FRA_TADs <- read.table("/home/andrea/Bioinformatica_Núria/TADs/TADs_OriginalFiles/T47D_TADs.bed", stringsAsFactors = F, sep="\t")
AJV_TADs_WT <- read.table("TADs_01_WT_NoCen.bed", stringsAsFactors = F, sep="\t")
AJV_TADs_KD <- read.table("TADs_02_multiH1-KD_NoCen.bed", stringsAsFactors = F, sep="\t")
AB_WT <- read.table("/home/andrea/Bioinformatica_Núria/HiC_TADs/AB_Compartments/ABcompartments_01_WT.bed", stringsAsFactors = F, sep="\t")
AB_KD <- read.table("/home/andrea/Bioinformatica_Núria/HiC_TADs/AB_Compartments/ABcompartments_02_multiH1-KD.bed", stringsAsFactors = F, sep="\t")
TADs_H12 <- read.table("H1.2_enriched_TADs.bed", stringsAsFactors = F, sep="\t")
TADs_H12like <- read.table("H1.2_enrichedLIKE_TADs.bed", stringsAsFactors = F, sep="\t")
TADs_H1Xlike <- read.table("H1X_enrichedLIKE_TADs.bed", stringsAsFactors = F, sep="\t")
TADs_H1X <- read.table("H1X_enriched_TADs.bed", stringsAsFactors = F, sep="\t")

# In order to optimize the script, the comparison between TADs are done per chromosome
Chromosomes <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                 "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")

#--- AJV vs François --------------------------------------------------------------------------------------------------
FRA_Total <- length(rownames(FRA_TADs)); AJV_WT_Total <- length(rownames(AJV_TADs_WT))
for (bin in c(100,1000,10000,100000)){
  Overlapping_Borders <- data.frame()
  for (chr in Chromosomes){
    FRABorders <- FRA_TADs[FRA_TADs$V1==chr,]$V3
    AJVBorders <- AJV_TADs_WT[AJV_TADs_WT$V1==chr,]$V3
    for (FRABorder in FRABorders){
      for (AJVBorder in AJVBorders){
        if (AJVBorder < FRABorder+bin & AJVBorder > FRABorder-bin){
          print(FRABorder); print(AJVBorder)
          Overlapping_Borders <- rbind(Overlapping_Borders, data.frame(chr, AJVBorder, FRABorder))
        }
      }
    }
  }
  do.call("<-", list(paste("Overlapping_Borders",as.integer(bin),sep="_"), Overlapping_Borders))
  write.table(Overlapping_Borders, file=paste(WorkingDir,"Overlapping_Borders_",bin,".bed",sep=""), quote=F, row.names=F, col.names=F, sep="\t")
}


FraAjv1 <- draw.pairwise.venn(FRA_Total, AJV_WT_Total, length(rownames(Overlapping_Borders_100)), category = c("François (100 bp)", "Albert"), lty=1, 
                              fill = c("lightblue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = c(0.08, 0.04), lwd=1)
FraAjv2 <- draw.pairwise.venn(FRA_Total, AJV_WT_Total, length(rownames(Overlapping_Borders_1000)), category = c("François (1,000 bp)", "Albert"), lty=1, 
                              fill = c("lightblue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = c(0.08, 0.04), lwd=1)
FraAjv3 <- draw.pairwise.venn(FRA_Total, AJV_WT_Total, length(rownames(Overlapping_Borders_10000)), category = c("François (10,000 bp)", "Albert"), lty=1,
                              fill = c("lightblue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = c(0.08, 0.04), lwd=1)
FraAjv4 <- draw.pairwise.venn(FRA_Total, AJV_WT_Total, length(rownames(Overlapping_Borders_100000)), category = c("      François (100,000 bp)", "Albert      "), lty=1,
                              fill = c("lightblue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = c(0.06, 0.02), lwd=1)

pdf(paste(WorkingDir,"Venn_FraAjv.pdf"))
grid.arrange(gTree(children=FraAjv1),gTree(children=FraAjv2),gTree(children=FraAjv3),gTree(children=FraAjv4),ncol=2)
dev.off()

############################################################################################################################################################
#                                                          TADs Groups in Dox+/Dox- Conditions Venn Di
############################################################################################################################################################
#--- Dox+ vs Dox- --------------------------------------------------------------------------------------------------
for (bin in c(100,1000,10000,100000)){
  Overlapping_Borders <- data.frame()
  for (chr in Chromosomes){
    WTBorders <- AJV_TADs_WT[AJV_TADs_WT$V1==chr,]
    KDBorders <- AJV_TADs_KD[AJV_TADs_KD$V1==chr,]
    for (WTB in WTBorders){
      for (KDB in KDBorders){
        if (WTBorders[WTB,3] < KDBorders[KDB,3]+bin & WTBorders[WTB,3] > KDBorders[KDB,3]-bin &
            WTBorders[WTB,2] < KDBorders[KDB,2]+bin & WTBorders[WTB,2] > KDBorders[KDB,2]-bin){
          Overlapping_Borders <- rbind(Overlapping_Borders, data.frame(chr, WTBorder, KDBorder))
        }
      }
    }
  }
  do.call("<-", list(paste("Overlapping_Dox_Borders",as.integer(bin),sep="_"), Overlapping_Borders))
  write.table(Overlapping_Borders, file=paste(WorkingDir,"Overlapping_Dox_Borders_",bin,".bed",sep=""), quote=F, row.names=F, col.names=F, sep="\t")
}

WT_Total <- length(rownames(AJV_TADs_WT)); KD_Total <- length(rownames(AJV_TADs_KD))

Dox1 <- draw.pairwise.venn(WT_Total, KD_Total, length(rownames(Overlapping_Dox_Borders_100)), category = c("      Dox- (100 bp)", "Dox+      "), lty=1, 
                           fill = c("lightblue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = c(0.05, 0.05), lwd=1)
Dox2 <- draw.pairwise.venn(WT_Total, KD_Total, length(rownames(Overlapping_Dox_Borders_1000)), category = c("      Dox- (1,000 bp)", "Dox+      "), lty=1, 
                           fill = c("lightblue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = c(0.05, 0.05), lwd=1)
Dox3 <- draw.pairwise.venn(WT_Total, KD_Total, length(rownames(Overlapping_Dox_Borders_10000)), category = c("      Dox-(10,000 bp)", "Dox+      "), lty=1,
                           fill = c("lightblue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = c(0.05, 0.05), lwd=1)
Dox4 <- draw.pairwise.venn(WT_Total, KD_Total, length(rownames(Overlapping_TADs)), category = c("                Dox- (100,000 bp)", "Dox+               "), lty=1,
                           fill = c("lightblue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = c(0.05, 0.05), lwd=1)

pdf(paste(WorkingDir,"Venn_Dox.pdf"))
grid.arrange(gTree(children=Dox1),gTree(children=Dox2),gTree(children=Dox3),gTree(children=Dox4),ncol=2)
dev.off()
############################################################################################################################################################
#                                                         TADs Groups and Dox+/Dox- Conditions
############################################################################################################################################################

#------ BarPlot TADs_A_WT,. TADs_B_WT, TADs_A_KD and TADs_B_KD ---------
pdf(paste(WorkingDir,"BarPlot_TADs_in_AB.pdf",sep=""), width = 6)
barplot(c(WT_A_Borders,WT_B_Borders,WT_N_Borders,KD_A_Borders,KD_B_Borders,KD_N_Borders), col = c("lightblue", "sandybrown"))
text(y=-150, x=c(1.3,3.7), labels=c("Wild-Type","Knock-Down"), cex = 1.1, xpd=TRUE)
text(y=c(WT_A_Borders,WT_B_Borders,KD_A_Borders,KD_B_Borders)+70, x= c(0.7,1.9,3.1,4.3), 
     labels = c(WT_A_Borders,WT_B_Borders,KD_A_Borders,KD_B_Borders), xpd=TRUE)
title(main="TADs in A/B Compartments", ylab = "TADs", xlab = "Condition", cex.lab=1.2, cex.main=1.2) 
legend("topright", c("A","B"), fill = c("lightblue", "sandybrown"))
dev.off()
############################################################################################################################################################



############################################################################################################################################################
#                                                                    Border Strength
############################################################################################################################################################
#--- Border Strength vs. TADs Groups (BarPlot) -------------------------------------------------------------------------------
Strength_vs_Groups <- data.frame(rep(0,10),rep(0,10),rep(0,10),rep(0,10))
colnames(Strength_vs_Groups) <- c("H1.2","H1.2like","H1Xlike","H1X")

for (strength in 1:10){
  if (strength %in% AJV_TADs_WT$V4){
    TADs_Strength <- AJV_TADs_WT[AJV_TADs_WT$V4==strength,]
    for (TAD1 in 1:length(rownames(TADs_Strength))){
      for (TAD2 in 1:length(rownames(TADs_H12))){
        if (TADs_Strength[TAD1,1] == TADs_H12[TAD2,1] & TADs_Strength[TAD1,2] == TADs_H12[TAD2,2] & TADs_Strength[TAD1,3] == TADs_H12[TAD2,3]){
          Strength_vs_Groups[strength,1] <- Strength_vs_Groups[strength,1]+1
        }
      }
      for (TAD3 in 1:length(rownames(TADs_H12like))){
        if (TADs_Strength[TAD1,1] == TADs_H12like[TAD3,1] & TADs_Strength[TAD1,2] == TADs_H12like[TAD3,2] & TADs_Strength[TAD1,3] == TADs_H12like[TAD3,3]){
          Strength_vs_Groups[strength,2] <- Strength_vs_Groups[strength,2]+1
        }
      }
      for (TAD4 in 1:length(rownames(TADs_H1Xlike))){
        if (TADs_Strength[TAD1,1] == TADs_H1Xlike[TAD4,1] & TADs_Strength[TAD1,2] == TADs_H1Xlike[TAD4,2] & TADs_Strength[TAD1,3] == TADs_H1Xlike[TAD4,3]){
          Strength_vs_Groups[strength,3] <- Strength_vs_Groups[strength,3]+1
        }
      }
      for (TAD5 in 1:length(rownames(TADs_H1X))){
        if (TADs_Strength[TAD1,1] == TADs_H1X[TAD5,1] & TADs_Strength[TAD1,2] == TADs_H1X[TAD5,2] & TADs_Strength[TAD1,3] == TADs_H1X[TAD5,3]){
          Strength_vs_Groups[strength,4] <- Strength_vs_Groups[strength,4]+1
        }
      }
    }
  }
}

BarPlot_Data <- rev(as.matrix(Strength_vs_Groups))

pdf("/home/andrea/Bioinformatica_Núria/HiC_TADs/TADs/T47D_TADsGroups_BorderStrength.pdf", width = 10)
barplot(c(BarPlot_Data[31:40],BarPlot_Data[21:30],BarPlot_Data[11:20],BarPlot_Data[1:10]), ylab = "TADs Borders", main="T47D - TADs Groups vs. Border Strengh",
        xlab="TADs Group", col=c("#3ca1c3","#63b4cf","#76bdd5","#8ac7db","#9dd0e1","lightblue","#b1dae7","#c4e3ed","#d8ecf3","#ebf6f9"), cex.axis = 1)
abline(h=c(50,400), col="grey15",lty=2, lwd=0.8)
text(x=c(0.8,2,3,4.3,5.5,6.8,8,9,10.3,11.5,12.7,13.9,15,16.3,17.5,18.7,19.8,21,22.3,23.5,24.7,25.8,27,28.3,29.5,30.7,31.9,33,34.3,35.5,36.7,38,39.1,40.3,41.5,42.8,43.9,45.1,46.3,47.6), 
     y =c(BarPlot_Data[31:40],BarPlot_Data[21:30],BarPlot_Data[11:20],BarPlot_Data[1:10])+12, cex = 0.9, xpd=TRUE,
     labels=as.character(c(BarPlot_Data[31:40],BarPlot_Data[21:30],BarPlot_Data[11:20],BarPlot_Data[1:10])))
text(y=-20, x=c(6,18,30,42), labels=c("H1.2","H1.2-like","H1X-like","H1X"), cex = 1, xpd=TRUE)
dev.off()


############################################################################################################################################################
#                                                                    TADs Interactions
############################################################################################################################################################
#--- Defining Interactions Groups ------------------------------------------------------
head(AJV_TADs_WT)
TADs_Interactions <- AJV_TADs_WT[,c(1,2,3,5)]; colnames(TADs_Interactions) <- c("Chr","From","To","Int")

Sorted <- TADs_Interactions[ match(sort(TADs_Interactions[,"Int"],decreasing = TRUE), TADs_Interactions[,"Int"]),]
length(rownames(Sorted))/4
Group_Int1 <- Sorted[1:802,]
Group_Int2 <- Sorted[803:1604,]
Group_Int3 <- Sorted[1605:2406,]
Group_Int4 <- Sorted[2407:3209,]

# BoxPlot showing the average length of both H1.2 and H1X enriched TADs
pdf(file="/home/andrea/Bioinformatica_Núria/HiC_TADs/TADs_InteractionsGroups_length.pdf", width = 8, height = 9)
Data <- c(c(Group_Int1[3]-Group_Int1[2]), c(Group_Int2[3]-Group_Int2[2]), 
          c(Group_Int3[3]-Group_Int3[2]), c(Group_Int4[3]-Group_Int4[2]))
boxplot(Data, outline = F, names = c("Int 1","Int 2","Int 3","Int 4"), main="Average TADs length",
        ylab = "TADs length (bp)", col=c("lightblue2","cyan3","coral2","sandybrown"), cex=0.9, cex.axis=0.9)
dev.off()

Group_Int1 <- Group_Int1[,c(1,2,3)]; Group_Int1 <- Group_Int1[ order( as.numeric(rownames(Group_Int1))) ,] 
Group_Int1[,"From"] <- as.integer(Group_Int1[,"From"]); Group_Int1[,"To"] <- as.integer(Group_Int1[,"To"]) 

Group_Int2 <- Group_Int2[,c(1,2,3)]; Group_Int2 <- Group_Int2[ order( as.numeric(rownames(Group_Int2))) ,] 
Group_Int2[,"From"] <- as.integer(Group_Int2[,"From"]); Group_Int2[,"To"] <- as.integer(Group_Int2[,"To"]) 

Group_Int3 <- Group_Int3[,c(1,2,3)]; Group_Int3 <- Group_Int3[ order( as.numeric(rownames(Group_Int3))) ,] 
Group_Int3[,"From"] <- as.integer(Group_Int3[,"From"]); Group_Int3[,"To"] <- as.integer(Group_Int3[,"To"]) 

Group_Int4 <- Group_Int4[,c(1,2,3)]; Group_Int4 <- Group_Int4[ order( as.numeric(rownames(Group_Int4))) ,] 
Group_Int4[,"From"] <- as.integer(Group_Int4[,"From"]); Group_Int4[,"To"] <- as.integer(Group_Int4[,"To"]) 

write.table(Group_Int1, file = "/home/andrea/Bioinformatica_Núria/HiC_TADs/TADs/TADs_Groups/Interactions1_TADs.bed", quote=F,col.names = F,row.names = F,sep = "\t")
write.table(Group_Int2, file = "/home/andrea/Bioinformatica_Núria/HiC_TADs/TADs/TADs_Groups/Interactions2_TADs.bed", quote=F,col.names = F,row.names = F,sep = "\t")
write.table(Group_Int3, file = "/home/andrea/Bioinformatica_Núria/HiC_TADs/TADs/TADs_Groups/Interactions3_TADs.bed", quote=F,col.names = F,row.names = F,sep = "\t")
write.table(Group_Int4, file = "/home/andrea/Bioinformatica_Núria/HiC_TADs/TADs/TADs_Groups/Interactions4_TADs.bed", quote=F,col.names = F,row.names = F,sep = "\t")

Group_Int1 <- Sorted[1:802,]; Group_Int1[,"Group"] <- 1
Group_Int2 <- Sorted[803:1604,]; Group_Int2[,"Group"] <- 2
Group_Int3 <- Sorted[1605:2406,]; Group_Int3[,"Group"] <- 3
Group_Int4 <- Sorted[2407:3209,]; Group_Int4[,"Group"] <- 4

TADs_Interactions_Groups <- rbind(Group_Int1,Group_Int2,Group_Int3,Group_Int4)

#--- Interactions vs. TADs Groups -------------------------------------------------------------------------------------------------------------------------
Interactions_vs_Groups <- data.frame(rep(0,4),rep(0,4),rep(0,4),rep(0,4))
colnames(Interactions_vs_Groups) <- c("H1.2","H1.2like","H1Xlike","H1X")

for (interaction in 1:4){
  if (interaction %in% TADs_Interactions_Groups$Group){
    TADs_Interactions <- TADs_Interactions_Groups[TADs_Interactions_Groups$Group==interaction,]
    for (TAD1 in 1:length(rownames(TADs_Interactions))){
      for (TAD2 in 1:length(rownames(TADs_H12))){
        if (TADs_Interactions[TAD1,1] == TADs_H12[TAD2,1] & TADs_Interactions[TAD1,2] == TADs_H12[TAD2,2] & TADs_Interactions[TAD1,3] == TADs_H12[TAD2,3]){
          Interactions_vs_Groups[interaction,1] <- Interactions_vs_Groups[interaction,1]+1
        }
      }
      for (TAD3 in 1:length(rownames(TADs_H12like))){
        if (TADs_Interactions[TAD1,1] == TADs_H12like[TAD3,1] & TADs_Interactions[TAD1,2] == TADs_H12like[TAD3,2] & TADs_Interactions[TAD1,3] == TADs_H12like[TAD3,3]){
          Interactions_vs_Groups[interaction,2] <- Interactions_vs_Groups[interaction,2]+1
        }
      }
      for (TAD4 in 1:length(rownames(TADs_H1Xlike))){
        if (TADs_Interactions[TAD1,1] == TADs_H1Xlike[TAD4,1] & TADs_Interactions[TAD1,2] == TADs_H1Xlike[TAD4,2] & TADs_Interactions[TAD1,3] == TADs_H1Xlike[TAD4,3]){
          Interactions_vs_Groups[interaction,3] <- Interactions_vs_Groups[interaction,3]+1
        }
      }
      for (TAD5 in 1:length(rownames(TADs_H1X))){
        if (TADs_Interactions[TAD1,1] == TADs_H1X[TAD5,1] & TADs_Interactions[TAD1,2] == TADs_H1X[TAD5,2] & TADs_Interactions[TAD1,3] == TADs_H1X[TAD5,3]){
          Interactions_vs_Groups[interaction,4] <- Interactions_vs_Groups[interaction,4]+1
        }
      }
    }
  }
}

BarPlot_Data_Int <- rev(as.matrix(Interactions_vs_Groups))

pdf("/home/andrea/Bioinformatica_Núria/HiC_TADs/TADs/T47D_TADsGroups_Interactions.pdf")
barplot(c(BarPlot_Data_Int[c(16,15,14,13)],BarPlot_Data_Int[c(12,11,10,9)],BarPlot_Data_Int[c(8,7,6,5)],BarPlot_Data_Int[c(4,3,2,1)]), 
        col=c("#3ca1c3","#9dd0e1","#c4e3ed","#ebf6f9"), cex.axis = 1.2,space=c(0.1,0.1,0.1,0.1,0.6,0.1,0.1,0.1,0.6,0.1,0.1,0.1,0.6,0.1,0.1,0.1) )
abline(h=c(150), col="grey15",lty=2, lwd=0.8)
title(ylab = "TADs Borders",cex.lab=1.2, line = 3)
title(main="T47D - TADs Groups vs. Interaction Density",cex.lab=1.3, line = 2.8)
title(xlab="TADs Group",cex.lab=1.2, line = 3.5)
text(x=c(0.55,1.65,2.85,3.95,5.45,6.6,7.7,8.8,10.35,11.5,12.55,13.7,15.3,16.4,17.45,18.6), 
     y =c(BarPlot_Data_Int[c(16,15,14,13)],BarPlot_Data_Int[c(12,11,10,9)],BarPlot_Data_Int[c(8,7,6,5)],BarPlot_Data_Int[c(4,3,2,1)])+10, cex = 1, xpd=TRUE,
     labels=as.character(c(BarPlot_Data_Int[c(16,15,14,13)],BarPlot_Data_Int[c(12,11,10,9)],BarPlot_Data_Int[c(8,7,6,5)],BarPlot_Data_Int[c(4,3,2,1)])))
text(y=-25, x=c(2.3,7.2,12,17), labels=c("H1.2","H1.2-like","H1X-like","H1X"), cex = 1, xpd=TRUE)
text(y=-10, x=c(0.55,1.65,2.85,3.95,5.45,6.6,7.7,8.8,10.35,11.5,12.55,13.7,15.3,16.4,17.45,18.6), labels=rep(1:4,4), cex = 1, xpd=TRUE)
dev.off()

#--- Border Strength / Interactions Density vs. TADs Groups (BoxPlots) ----------------------------------------------------------------------------------------------------

for (condition in c("WT","KD")){
  if(condition=="WT"){AnalyzedTADs <- AJV_TADs_WT}; if(condition=="KD"){AnalyzedTADs <- AJV_TADs_KD}
  for (propertie in c("BoxPlot_Interactions","BoxPlot_Strength")){
    Str_H12 <- c(); Str_H12like <- c(); Str_H1Xlike <- c(); Str_H1X <- c()
    Col <- ifelse (propertie =="BoxPlot_Strength",4,5)
    for (chr in Chromosomes){ 
      TADs_Chr <- AnalyzedTADs[AnalyzedTADs$V1==chr,]
      for (TAD1 in 1:length(rownames(TADs_Chr))){
        for (TAD2 in 1:length(rownames(TADs_H12))){
          if (TADs_Chr[TAD1,1] == TADs_H12[TAD2,1] & TADs_Chr[TAD1,2] == TADs_H12[TAD2,2] & TADs_Chr[TAD1,3] == TADs_H12[TAD2,3]){
            Str_H12 <- c(Str_H12, TADs_Chr[TAD1,Col])
          }
        }
        for (TAD2 in 1:length(rownames(TADs_H12like))){
          if (TADs_Chr[TAD1,1] == TADs_H12like[TAD2,1] & TADs_Chr[TAD1,2] == TADs_H12like[TAD2,2] & TADs_Chr[TAD1,3] == TADs_H12like[TAD2,3]){
            Str_H12like <- c(Str_H12like, TADs_Chr[TAD1,Col])
          }
        }
        for (TAD2 in 1:length(rownames(TADs_H1Xlike))){
          if (TADs_Chr[TAD1,1] == TADs_H1Xlike[TAD2,1] & TADs_Chr[TAD1,2] == TADs_H1Xlike[TAD2,2] & TADs_Chr[TAD1,3] == TADs_H1Xlike[TAD2,3]){
            Str_H1Xlike <- c(Str_H1Xlike, TADs_Chr[TAD1,Col])
          }
        }
        for (TAD2 in 1:length(rownames(TADs_H1X))){
          if (TADs_Chr[TAD1,1] == TADs_H1X[TAD2,1] & TADs_Chr[TAD1,2] == TADs_H1X[TAD2,2] & TADs_Chr[TAD1,3] == TADs_H1X[TAD2,3]){
            Str_H1X <- c(Str_H1X, TADs_Chr[TAD1,Col])
          }
        }
      }
    }
    do.call("<-", list(paste(propertie,condition,sep = "_"), list(Str_H12,Str_H12like,Str_H1Xlike,Str_H1X)))
  }
}

pdf("Interactions_WTvsKD_TADsGroups.pdf")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
boxplot(list(BoxPlot_Interactions_WT[[1]],BoxPlot_Interactions_KD[[1]],BoxPlot_Interactions_WT[[2]],BoxPlot_Interactions_KD[[2]],
             BoxPlot_Interactions_WT[[3]],BoxPlot_Interactions_KD[[3]],BoxPlot_Interactions_WT[[4]],BoxPlot_Interactions_KD[[4]]),
        outline = F, names=c("WT","KD","WT","KD","WT","KD","WT","KD"), main="Interactions Density in TADs Groups", cex.main=1.4,cex.axis=1.1, cex.lab=1.3,
        ylab="Interactions Density", xlab="Condition",col = c("lightblue2","lightblue2","cyan3","cyan3","coral2","coral2","sandybrown","sandybrown"))
        legend("topright", inset=c(-0.28,0), legend=c("H1.2","H1.2-like","H1X-like","H1X"), title="TADs Group", fill = c("lightblue2","cyan3","coral2","sandybrown"))
dev.off()

pdf("BorderStrength_WTvsKD_TADsGroups.pdf") 
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
boxplot(list(BoxPlot_Strength_WT[[1]],BoxPlot_Strength_KD[[1]],BoxPlot_Strength_WT[[2]],BoxPlot_Strength_KD[[2]],
             BoxPlot_Strength_WT[[3]],BoxPlot_Strength_KD[[3]],BoxPlot_Strength_WT[[4]],BoxPlot_Strength_KD[[4]]), ylim=c(0.5,10.5),
        outline = F, names=c("WT","KD","WT","KD","WT","KD","WT","KD"), main="Border Strength in TADs Groups", cex.main=1.4,cex.axis=1.1, cex.lab=1.3,
        ylab="Border Strength", xlab="Condition",col = c("lightblue2","lightblue2","cyan3","cyan3","coral2","coral2","sandybrown","sandybrown"))
        legend("topright", inset=c(-0.28,0), legend=c("H1.2","H1.2-like","H1X-like","H1X"), title="TADs Group", fill = c("lightblue2","cyan3","coral2","sandybrown"))
dev.off()

#--- Border Strength vs. Interactions Density (ScatterPLot) ----------------------------------------------------------------------------------------------------

StrInt_WT <- data.frame(AJV_TADs_WT$V4, AJV_TADs_WT$V5); StrInt_KD <- data.frame(AJV_TADs_KD$V4, AJV_TADs_KD$V5)
colnames(StrInt_WT) <- c("Strength", "Interactions"); colnames(StrInt_KD) <- c("Strength", "Interactions")
wt <- cor.test(StrInt_WT$Strength, StrInt_WT$Interactions, method = "pearson")
kd <- cor.test(StrInt_KD$Strength, StrInt_KD$Interactions, method = "pearson")

library(ggplot2)
pdf("/home/andrea/Bioinformatica_Núria/HiC_TADs/TADs/Correlation_Strength_vs_Interactions.pdf")
Strengths_Interactions_WT <- ggplot(StrInt_WT, aes(x=Strength, y=Interactions)) +
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour="ligthblue"),  size=0.2, show.legend = F) + theme_classic() + 
  theme_classic() + theme(plot.title = element_text(hjust = 0.5, color = "#666666"),plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
  ggtitle(label = "Border Strength vs Interactions Density (WT)", subtitle = paste("p.value = ",format(wt$p.value,digits=3),"\ncorr.coef = ",round(wt$estimate,3),sep="")) + scale_color_manual(values="steelblue1")
Strengths_Interactions_WT

Strengths_Interactions_WT <- ggplot(StrInt_KD, aes(x=Strength, y=Interactions)) +
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour="ligthblue"),  size=0.2, show.legend = F) + theme_classic() + 
  theme_classic() + theme(plot.title = element_text(hjust = 0.5, color = "#666666"),plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
  ggtitle(label = "Border Strength vs Interactions Density (KD)", subtitle = paste("p.value = ",format(kd$p.value,digits=3),"\ncorr.coef = ",round(kd$estimate,3),sep="")) + scale_color_manual(values="steelblue1")
Strengths_Interactions_WT
dev.off()



