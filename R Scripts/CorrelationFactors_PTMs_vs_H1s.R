#############################################################################################################################################################
#                                          Correlation Factors (HeLa1 vs PTMs) - HeLa 1 experiment
#############################################################################################################################################################
# Creation of scatter plots in order to see correlations between Hela1 H1 variants and PTMs (by grouping with both Gpos and Gpos+Gneg)
print("5| H1 VARIANTS vs PTMs&TFs ANALYSIS")
if (length(PTMs)==12){
  
      #--- Merging, Sorting and Scaling the Data ------------------------------------------------------------------------------------------------------------
      H1variants_PTMs_CRs <- merge(TOTAL_PTMsCRs, TOTAL_H1Variants, by="row.names")
      head(H1variants_PTMs_CRs); dim(H1variants_PTMs_CRs)
      
      rownames(H1variants_PTMs_CRs) <- H1variants_PTMs_CRs$Row.names
      H1variants_PTMs_CRs <- H1variants_PTMs_CRs[, -1]  # Removes the first column, which just stored the rownames
      head(H1variants_PTMs_CRs)
      
      # Re-ordering (necessary because when merging the data, values are sorted numerically according to the chromosome band)
      H1variants_PTMs_CRs$V1 <- rownames(H1variants_PTMs_CRs)
      order_H1abundance_script <- rownames(TOTAL_H1Variants) 
      
      # The final data frame will contain information about the analyzed PTMs&TFs and histone H1
      H1variants_PTMs_CRs <- H1variants_PTMs_CRs[match(order_H1abundance_script, H1variants_PTMs_CRs$V1),]
      H1variants_PTMs_CRs <- H1variants_PTMs_CRs[, c(PTMs,Histones)] 
      head(H1variants_PTMs_CRs); dim(H1variants_PTMs_CRs)
      
      # Scaling the data (required because both PTMs and H1 variants are going to be merged)
      H1variants_PTMs_CRs_scaled <-  as.data.frame(apply(H1variants_PTMs_CRs, 2, scale))
      rownames(H1variants_PTMs_CRs_scaled) <- rownames(H1variants_PTMs_CRs)
      head(H1variants_PTMs_CRs_scaled); dim(H1variants_PTMs_CRs_scaled); range(H1variants_PTMs_CRs_scaled)
      ScatterData <- H1variants_PTMs_CRs_scaled
      
      # This step is required in order to not taking into account extrem values by considering 1st and 3rd quartiles
      UpLim <- c(); DownLim <- c()
      for (histone in Histones){
        DownLim <- c(DownLim, summary(ScatterData[,histone])[[2]]- 4)
        UpLim <- c(UpLim, summary(ScatterData[,histone])[[5]]+ 4)
      }
      UpLim <- mean(UpLim); DownLim <- mean(DownLim)
      
      # Removing rows with extrem values
      dim(ScatterData)
      for (histone in Histones){
        ScatterData <- ScatterData[(ScatterData[,histone]<UpLim),]
        ScatterData <- ScatterData[(ScatterData[,histone]>DownLim),]
      }
      dim(ScatterData)
      
      Number_of_Bands <- c()
      for (band in CytoBands){
        Number_of_Bands <- c(Number_of_Bands, length(grep(band,rownames(ScatterData))))
      }
      
      CytoBands_Scatter <- c()
      for (band in 1:length(CytoBands)){
        CytoBands_Scatter <- c(CytoBands_Scatter, rep(CytoBands[band], Number_of_Bands[band]))
      }
      ScatterData$Cytoband_Group <- CytoBands_Scatter
      head(ScatterData)
      
      #--- Correlation Tests --------------------------------------------------------------------------------------------------------------------------------
      x <- c(); y <- c(); z <- c()
      Cols <- c()
      
      # Storing in x, y and z the p.values, estimates and histone names, respectively
      for (histone in Histones){
        for (PTM in PTMs){
          Result <- cor.test(ScatterData[ ,histone] , ScatterData[ ,PTM], method = "pearson") 
          x <- c(x, format(Result$p.value,digits=3))
          y <- c(y, format(Result$estimate,digits=3))
          Cols <- c(Cols, PTM)
        }
        z <- c(z, rep(histone, length(PTMs)) )
        
      }
      
      for (histone1 in Histones){
        for (histone2 in Histones){
          if (histone1 != histone2){
            Result <- cor.test(ScatterData[ ,histone1] , ScatterData[ ,histone2], method = "pearson") 
            if (!(Result$estimate %in% y)){
              x <- c(x, Result$p.value)
              y <- c(y, Result$estimate)
              z <- c(z, histone2)
              Cols <- c(Cols, histone1)
            }
          }
        }
      } 
      Correlation_Table <- data.frame(x, y, z)
      Correlation_Table <- as.data.frame(t(Correlation_Table),stringsAsFactors = FALSE)
      row.names(Correlation_Table) <- c("p.value", "corr.coef", "VS.Variant")
      colnames(Correlation_Table)  <- Cols
    
      
      #--- Scatter Plots -----------------------------------------------------------------------------------------------------------------------------------
      # The "ggplot" function is going to be used to represent the calculated correlation factors
      Scatter_Colors <- c("Gneg1"="#b3d9ff","Gneg2"="#99e699","Gneg3"="#ffccb3","Gneg4"="#e59acc",
                          "Gpos25"="#66b3ff","Gpos50"="#47d147","Gpos75"="#ffaa80","Gpos100"="#d147a3")
      Scatters_Names <- c()
      
      for (histone in Histones){
        Histone_Table <- Correlation_Table[ ,Correlation_Table["VS.Variant",]==histone]
        for (element in colnames(Histone_Table)){
          do.call("<-", list(paste(histone,element,sep="_"),
                             ggplot(ScatterData, aes_string(x=element, y=histone)) +
                               geom_smooth(method = "lm", se = FALSE, color="black") + 
                               geom_point(aes(colour=Cytoband_Group),  show.legend = F, size=1) + theme_classic() + 
                               theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
                                                       plot.title = element_text(hjust = 0.5, color = "#666666"), 
                                                       plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
                               scale_color_manual(values=Scatter_Colors) +
                               ggtitle(label=paste(element,"vs",histone,sep=" "), 
                                       subtitle = paste("p.value = ",Histone_Table["p.value", element],"\ncorr.coef = ",Histone_Table["corr.coef", element],sep="")) 
          ))
         Scatters_Names <- c(Scatters_Names, paste(histone,element,sep="_")) 
        }
      }
      
      #  H1.A vs H1.B  ######################################################################################
      H1.A_H1.B <- ggplot(ScatterData, aes_string(x=Histones[2], y=Histones[1])) + 
        geom_smooth(method = "lm", se = FALSE, color="black") + 
        geom_point(aes(colour=Cytoband_Group),  show.legend = T, size=1) + theme_classic() + 
        theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
                                plot.title = element_text(hjust = 0.5, color = "#666666"),
                                plot.subtitle = element_text(hjust = 0.5, color = "#666666")) + 
        ggtitle(label = "H1.2 vs H1.X", subtitle = paste("p.value = ",format(as.numeric(Correlation_Table["p.value", Histones[1]]),digit=3),
                                                         "\ncorr.coef = ",format(as.numeric(Correlation_Table["corr.coef", Histones[1]]),digit=3),sep=""))  
      H1.A_H1.B <- H1.A_H1.B + scale_color_manual(values=Scatter_Colors)
      
      # Opening the PDF file which will store all the graphs and executing "Command" strings as code
      pdf(paste(ResultsDir,"/H1_",CellLine,"_H1vsPTMs_Correlations_",Bands_Studied,".pdf",sep=""), width = 12)
      
      Marks <- list(c("H3K9me3", "H3K27me3", "H3K4me3", "H3K4me1", "H3K27ac", "H3K9ac"),
                    c("H3K36me3", "EZH2", "CTCF", "P300", "CHD1", "RNAPOLII"))
      
      # This step is required because there is no way to use "grid.arrange" function without manually specifying ggplot2 originated objects
      # A solution is to construct the required query in a string format and then evaluate it as code with "eval(parse(text= STRING))"
      for (histone in Histones){
        Command <- ""
        for (set in Marks){
          for (mark in set){
            Command <- paste(Command,paste("eval(as.name('",histone,"_",mark,"')),",sep=""),sep="") 
            }
          }
        Command <- paste("grid.arrange(",Command,"nrow = 2)",sep="")
        eval(parse(text=Command))
      }
      
      print(H1.A_H1.B)
      dev.off() # Closing the connection with the PDF 
}

#############################################################################################################################################################