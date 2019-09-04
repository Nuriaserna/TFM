#############################################################################################################################################################
#                                                           H1 Correlation Tests: Scatter Plots
#############################################################################################################################################################

#--- Preparing the Data ---------------------------------------------------------------------------------------------------------------------------
ScatterData <- H1_variants_scaled[,c(1:length(Histones))]

# This step is required in order to not taking into account extrem values by considering 1st and 3rd quartiles
UpLim <- c(); DownLim <- c()
for (histone in Histones){
  DownLim <- c(DownLim, summary(ScatterData[,histone])[[2]]- 0.5)
  UpLim <- c(UpLim, summary(ScatterData[,histone])[[5]]+ 0.3)
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
ScatterData$CytoBand <- CytoBands_Scatter
head(ScatterData)


#--- Correlation Tests ----------------------------------------------------------------------------------------------------------------------------
x <- c(); y <- c(); z <- c()

# By using a "for" loop vectors x, y and z will incorporate p.values, estimates and histone names, respectively
# A second for loop will obtain the values for the correlations between histones
Cols <- c()
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

#--- Scatter Plots --------------------------------------------------------------------------------------------------------------------------------
# Generating the Scatter Plot showing the correlation between H1.2b vs. H1X within TADs
Scatter_Colors <- c("Gneg1"="#b3d9ff","Gneg2"="#99e699","Gneg3"="#ffccb3","Gneg4"="#e59acc",
                    "Gpos25"="#66b3ff","Gpos50"="#47d147","Gpos75"="#ffaa80","Gpos100"="#d147a3")

Scatters_Names <- c()

for (X in 1:length(Correlation_Table)){
  Leg <- ifelse (X == length(Correlation_Table), TRUE, FALSE)
  histoneY <- colnames(Correlation_Table)[X]
  histoneX <- Correlation_Table["VS.Variant",X]
  do.call("<-", list(paste(histoneX,histoneY,sep="_"),
                      ggplot(ScatterData, aes_string(x=histoneX, y=histoneY)) +
                        geom_smooth(method = "lm", se = FALSE, color="black") + 
                        geom_point(aes(colour=CytoBand),  show.legend = F, size=1) + theme_classic() + 
                        theme_classic() + theme(plot.title = element_text(hjust = 0.5, color = "#666666"), 
                                                plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
                        scale_color_manual(values=Scatter_Colors) +
                        ggtitle(label=paste(histoneX,"vs",histoneY,sep=" "), 
                                subtitle = paste("p.value = ",format(as.numeric(Correlation_Table["p.value", X]),digits=4),"\ncorr.coef = ",
                                                 format(as.numeric(Correlation_Table["corr.coef", X]),digit=3),sep="")) 
  )) 
  Scatters_Names <- c(Scatters_Names, paste(histoneX,histoneY,sep="_")) 
  
}

Command <- ""
for (corr in Scatters_Names){
  Command <- paste(Command,paste("eval(as.name('",corr,"')),",sep=""),sep="") 
}

pdf(paste(ResultsDir,"/H1_",CellLine,"_H1s_Correlations_",Bands_Studied,".pdf",sep=""),width = 20, height = 12)
Command <- paste("grid.arrange(",Command,"nrow = 2)",sep="")
eval(parse(text=Command))
dev.off()

#############################################################################################################################################################