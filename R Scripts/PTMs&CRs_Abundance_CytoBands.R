##############################################################################################################################################
#                                                 PTMs&CRs Abundance at CytoBands 
##############################################################################################################################################
print("3| PTMs&TFs ANALYSIS")
if (length(PTMs)==12){
      # Preparing loop's required variables
      Plots_Names <- c("TFs_5","PTMs_5","PTMs_7")
      TFs_5 <- c("RNAPOLII","CHD1","CTCF","EZH2","P300")
      PTMs_5 <- c("H3K27ac","H3K4me3","H3K36me3", "H3K27me3", "H3K9me3")
      PTMs_7 <- c("H3K4me1","H3K27ac","H3K4me3","H3K9ac","H3K36me3","H3K27me3","H3K9me3")
      
      # Opening a connection with a new PDF file
      pdf(paste(ResultsDir,"/H1_",CellLine,"_PTMs&CRs_Cytobands_",Bands_Studied,".pdf",sep=""), width = 9)
      for (plot in Plots_Names){
        Plot_Data <- c()
        SpacesVector <- c()
        Count <- 1
        for (band in CytoBands){
          for (mark in eval(as.name(plot))){
            # Generating the data which will be used when performing each plot (for each band, for each mark: add their information to "Plot_Data")
            Plot_Data <- c(Plot_Data, list(TOTAL_PTMsCRs_scaled[grep(band, rownames(TOTAL_PTMsCRs_scaled)),mark]) )
            SpacesVector <- c(SpacesVector, Count)
            Count <- Count+1
          }
          Count <- Count+1
        }
        Colors <- c()
        for (PTM in eval(as.name(plot))){
          Colors <- c(Colors, ColorCodePTMs[[PTM]])
        }
        boxplot(Plot_Data, col = Colors, outline = F, ylab = paste(plot," abundance",sep=""), xaxt="n", cex.lab=1.2, xlab = "CytoBand",
                names = rep(eval(as.name(plot)),length(CytoBands)), xaxt="n", cex = 1, cex.axis = 1, at = SpacesVector, cex.main=1.3,
                main = paste(plot," abundance at ",Bands_Studied," cytobands (normalyzed and scaled)",sep=""), xlab=NULL)
        abline(h = 0, col = "grey", lwd = 1)
        abline(h = 1, col = "#66b3ff", lwd = 1)
        abline(h = -1, col = "#d147a3", lwd = 1)
        boxplot(Plot_Data, col = Colors, outline = F, xaxt="n", yaxt="n",at = SpacesVector, add = TRUE)
        
        # Adding the name of the groups to the x axis and lines separating grous of boxes of H1 variants belonging to a specific group
        # Both "GroupPosition" and "LinePosition" variables will depend on the number of histones analyzed
        VectorPositions <- c()
        BandPosition <- length(eval(as.name(plot)))-length(eval(as.name(plot)))/2.5
        for (band in 1:length(CytoBands) ){
          VectorPositions <- c(VectorPositions, BandPosition) 
          BandPosition <- BandPosition +1+ length(eval(as.name(plot)))
        }
        axis(1, VectorPositions, CytoBands, cex.axis=1)
        
        LinePosition <- length(eval(as.name(plot)))+1
        for (band in 1:(length(CytoBands)-1)){
          abline(v=LinePosition, col = "grey", lwd = 1)
          LinePosition <- LinePosition +1+ length(eval(as.name(plot)))
        }
        
      legend("topright", inset=.01, eval(as.name(plot)), fill = Colors, cex=1, bg ="white" )
      } 
      
      dev.off()
}

