##############################################################################################################################################
#                                                      H1 Abundance at CytoBands 
##############################################################################################################################################
print("4| H1 VARIANTS ANALYSIS")
# Opening a connection with a new PDF file
pdf(paste(ResultsDir,"/H1_",CellLine,"_Abundance_Cytobands_",Bands_Studied,".pdf",sep=""), width = 8, height = 9)
Colors <- c()
for (histone in Histones){
  Colors <- c(Colors, ColorCodeH1s[[histone]])
}
Plot_Data <- c()
SpacesVector <- c()
Count <- 1
for (band in CytoBands){
  for (histone in colnames(TOTAL_H1Variants)){
    # Generating the data which will be used when performing each plot (for each band, for each histone: add their information to "Plot_Data")
    Plot_Data <- c(Plot_Data, list(TOTAL_H1Variants[grep(band, rownames(TOTAL_H1Variants)),histone]) )
    SpacesVector <- c(SpacesVector, Count)
    Count <- Count+1
  }
  Count <- Count+1
}
boxplot(Plot_Data, col = Colors, outline = F, ylab = "ChIP-Seq Signal", cex.lab=1.2,
        xaxt="n", cex = 1, cex.axis = 1, at = SpacesVector, cex.main = 1.3, xlab = "CytoBand",
        main = paste(CellLine," H1 variants abundance at ",Bands_Studied," cytobands",sep=""), xlab=NULL)
  
# Adding the name of the groups to the x axis and lines separating grous of boxes of H1 variants belonging to a specific group
# Both "GroupPosition" and "LinePosition" variables will depend on the number of histones analyzed
VectorPositions <- c()
BandPosition <- length(colnames(TOTAL_H1Variants))-length(colnames(TOTAL_H1Variants))/2.5
for (band in 1:length(CytoBands) ){
  VectorPositions <- c(VectorPositions, BandPosition) 
  BandPosition <- BandPosition +1+ length(colnames(TOTAL_H1Variants))
}
axis(1, VectorPositions, CytoBands, cex.axis=1)
    
LinePosition <- length(colnames(TOTAL_H1Variants))+1
for (band in 1:(length(CytoBands)-1)){
  abline(v=LinePosition, col = "grey", lwd = 1)
  LinePosition <- LinePosition +1+ length(colnames(TOTAL_H1Variants))
}
abline(h = 0, col = "grey", lwd = 1)
abline(h = 0.005, col = "#66b3ff", lwd = 1)
abline(h = -0.005, col = "#d147a3", lwd = 1)
boxplot(Plot_Data,at = SpacesVector, col = Colors, outline = F, add=TRUE, xaxt="n", yaxt="n") 
legend("topright", inset=.01, colnames(TOTAL_H1Variants), fill = Colors, cex=0.95, bg ="white" )
dev.off()
  
#############################################################################################################################################################