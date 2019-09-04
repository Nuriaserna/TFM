#!/usr/bin/env Rscript
#############################################################################################################################################################
#                                                           Donstream Analysis Pipeline
#############################################################################################################################################################
#--- Loading the required Libraries -------------------------------------------------------------------------------------------------------------------------
options(warn=-1)
library("optparse", warn.conflicts = FALSE) 
library(ggplot2, warn.conflicts = FALSE)
library(gridExtra, warn.conflicts = FALSE)

#--- Arguments Storage --------------------------------------------------------------------------------------------------------------------------------------
option_list = list(
  make_option(c("-l", "--line"), action="store", type="character", help="Cell Line", metavar="character"),
  make_option(c("-y", "--pathptm"), action="store", type="character", help="PTMs&TFs Path", metavar="character"),
  make_option(c("-x", "--pathh1"), action="store", type="character", help="Histone H1 Variants Path", metavar="character"),
  make_option(c("-e", "--exchromosomes"), action="store", type="character", help="Chromosomes to be excluded", metavar="character"),
  make_option(c("-f", "--fileh1"), action="store", type="character", help="Histone H1 Variants to be analyzed (Filename Identifier)", metavar="comma-separated list"),
  make_option(c("-n", "--nameh1"), action="store", type="character", help="Histone H1 Variants to be analyzed (H1 Variant Common Name)", metavar="comma-separated list")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (length(opt)<6){
  stop("At least one argument must be supplied\nUsage: Rscript Main_Script.R -a <celline> -b <pathH1> -c <pathPTMs> -d <chr> -e <H1variants>", call.=FALSE)
}

CellLine <- c(opt$line)
PTMsDir <- c(opt$pathptm)
HistonesDir <- c(opt$pathh1)
DiscardedBands <- c(opt$exchromosomes)
HistonesFiles <- c(); Histones <- c()
for (i in strsplit(opt$fileh1,",")){HistonesFiles <- c(HistonesFiles, i)}
for (i in strsplit(opt$nameh1,",")){Histones <- c(Histones, i)}

print("1| ARGUMENTS STORAGE")
print(paste("   Cell Line: ",CellLine)) 
print(paste("   PTMs&TFs Directory: ",PTMsDir))
print(paste("   H1 Variants Directory: ",HistonesDir))
print(paste("   Chromosomes Excluded: ",DiscardedBands))
print(paste("   H1 Variants Filenames: ",opt$fileh1))
print(paste("   H1 Variants Names: ",opt$nameh1))

#--- Constants Storage --------------------------------------------------------------------------------------------------------------------------------------
load("/home/andrea/Bioinformatica_NÃºria/R_Scripts/GnegGposGroups.Rdata")
Bands_Studied <- "Gpos&Gneg"
CytoBands <- c("Gneg1","Gneg2","Gneg3","Gneg4","Gpos25","Gpos50","Gpos75","Gpos100")
Bands_Colors <- c("Gneg1"="#b3d9ff","Gneg2"="#99e699","Gneg3"="#ffccb3","Gneg4"="#e59acc",
                  "Gpos25"="#66b3ff","Gpos50"="#47d147","Gpos75"="#ffaa80","Gpos100"="#d147a3")
ColorCodeH1s <- list("H1.2"="lightblue","H1X"="sandybrown","H1Xb"="sienna1","H1.4"="gold","H1.5"="darkolivegreen1","H1.5b"="darkolivegreen3",
                     "H1.2_WT"="lightblue","H1X_WT"="sandybrown","H1.4_WT"="#ffe866","H1.5_WT"="darkolivegreen1",
                     "H1.2_KD"="lightblue3","H1X_KD"="sienna1","H1.4_KD"="gold","H1.5_KD"="darkolivegreen3")
PTMs <- c("RNAPOLII","CHD1","CTCF","H3K4me1","H3K27ac","H3K4me3","H3K9ac","P300","H3K36me3","EZH2","H3K27me3","H3K9me3")
ColorCodePTMs <- list("CTCF"="lightblue1","H3K4me3"="lightgreen","H3K27ac"="aquamarine4","H3K9ac"="darkolivegreen1", 
                   "RNAPOLII"="royalblue","CHD1"="steelblue1","H3K4me1"="cyan3","H3K36me3"="peachpuff",
                   "P300"="gold","H3K27me3"="coral2","EZH2"="sandybrown","H3K9me3"="brown")

#--- Results Directory Creation -----------------------------------------------------------------------------------------------------------------------------
dir.create(paste(HistonesDir, CellLine,"_Results",sep=""), showWarnings = FALSE)
ResultsDir <-paste(HistonesDir, CellLine,"_Results",sep="")

#--- Scripts Calling ----------------------------------------------------------------------------------------------------------------------------------------
source("/home/andrea/Bioinformatica_NÃºria/R_Scripts/TFM/TFM_Reading_PTMs&CRs_H1Variants_Files.R")
source("/home/andrea/Bioinformatica_NÃºria/R_Scripts/TFM/TFM_PTMs&CRs_Abundance_CytoBands.R")
source("/home/andrea/Bioinformatica_NÃºria/R_Scripts/TFM/TFM_H1_Abundance_CytoBands.R")
source("/home/andrea/Bioinformatica_NÃºria/R_Scripts/TFM/TFM_H1s_Correlations.R")
source("/home/andrea/Bioinformatica_NÃºria/R_Scripts/TFM/TFM_CorrelationFactors_PTMs_vs_H1s.R")

#############################################################################################################################################################
