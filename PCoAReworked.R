# PCoA with adjustment for working with OTUtables

## Author Matthew C. B. Tsilimigras
## 6-14-2017
## Generates PCoA Plotting and MDS axes values and eigenvales
## Only works for qiime currently
## ROLE IN THE MANUSCRIPT AND ITS SUPPLEMENTS:
## Used to generate MDS axes required for MouseStressModeling.R to run correctly
## Otherwise not used

rm(list=ls())

library("vegan")
library("calibrate")
#library("psych")

otuTable <- TRUE
dataType <- "closedQIIMER1"
## dataType <- "denovoabundantOTU"

if (.Platform$OS.type == "windows") {
  if (dataType == "closedQIIMER1") {
    baseDir <- "F:/Google Drive/Datasets/MarkExperiment/closedQIIMER1/"
  } else {
    baseDir <- "/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/ArgonneSequencing/JanuaryResequencing/rawData/rg_results/abundantOTU"
  }
} else {
  if (dataType == "closedQIIMER1") {
    baseDir <- "/Users/mbrown67/Google Drive/Datasets/MarkExperiment/closedQIIMER1/"
    ##baseDir <- "/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/ArgonneSequencing/JanuaryResequencing/rawData/rg_results/closedQIIMER1"
    
  } else {
    baseDir <- "/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/ArgonneSequencing/JanuaryResequencing/rawData/rg_results/abundantOTU"
  }
}

metadataDir <- paste(baseDir, "metadata", sep="/")
## setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/ArgonneSequencing/JanuaryResequencing/metadata/")
dataDir <- paste(baseDir, "data", sep="/")
processedDir <- paste(baseDir, "processed", sep="/")
analysisDir <- paste(baseDir, "analysis", sep="/")
if (otuTable == TRUE) {
  taxonomicLevels <- 1
} else {
  taxonomicLevels <- c(2:7)
}
for (taxa in taxonomicLevels){
  setwd(processedDir)
  if( otuTable == TRUE) {
    if ( dataType == "closedQIIMER1") {
      inFileName <- "qiimeOTULogNormwithMetadata.txt"
    }
    else {
      inFileName <- "abundantOTUR1LogNormwithMetadata.txt"
    }
  } else {
    
    inFileName <- paste0(dataType, "_L", taxa, "_LogNormwithMetadata.txt")
  }
  myT <- read.table(inFileName, header=TRUE, sep="\t", comment.char="@")
  ## Some concerns about column classes
  setwd(analysisDir)
  if(otuTable == TRUE) {
    endMetadataIndex <- which(colnames(myT) == "counts")
  } else {
    endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
  }
  
  ## Display all 5 tissue sample types
  myMDS <- capscale(myT[,(endMetadataIndex +1):ncol(myT)]~1,distance="bray")
  
  write.table(myMDS$CA$u, sep="\t", file=paste0(dataType, "_L_", taxa, "_pcoa.txt"))
  write.table(myMDS$CA$eig, file=paste0(dataType, "_L_", taxa, "_eigenValues.txt"), sep="\t")
  
  percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
  
  pdf( paste0(dataType, "_L_", taxa, "_topMDS.pdf"))
  ## Color by tissue
  for (xrun in 1:5) {
    for (yrun in 2:5) {
      if(xrun == yrun){
        break
      }
      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main=paste0("PCoA at level ", taxa, " \n(Sample Tissue)"), cex=2.0, ##bty="L",
           pch = 16,
           col=ifelse(myT$Source=="Cecal Content", "green", ifelse(myT$Source == "duo", "black", ifelse(myT$Source == "feces", "brown", ifelse(myT$Source == "ileum", "yellow", "red")))))
      par(xpd=TRUE)
      
      legend("topright", inset=c(-0.9,0),
             c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
             pch = 16, cex = 1.1,
             col=c("green", "blue", "brown", "yellow", "red"))
      ## Color by Sex
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main=paste0("PCoA at level", taxa," \n(Mouse Sex)"), cex=2.0, ##bty="L",
           pch = 16,
           col=ifelse(myT$Sex == "male", "blue", "pink"))
      par(xpd=TRUE)
      
      legend("topright", inset=c(-0.9,0),
             c("Male", "Female"),
             pch = 16, cex = 1.1,
             col=c("blue", "pink"))
      ## Color by Stress
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main=paste0("PCoA at level", taxa," \n(Stress/Unstressed)"), cex=2.0, ##bty="L",
           pch = 16,
           col=ifelse(myT$Exp.or.Ctrl == "Exp", "red", "black"))
      par(xpd=TRUE)
      
      legend("topright", inset=c(-0.9,0),
             c("Stress", "Control"),
             pch = 16, cex = 1.1,
             col=c("red", "black"))
      
      ## Color by cage
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main=paste0("PCoA at level ", taxa," \n(Cage)"), cex=2.0, ##bty="L",
           pch = 16,
           col=ifelse(myT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(myT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(myT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(myT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(myT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(myT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(myT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
      par(xpd=TRUE)
      
      legend("topright", inset=c(-0.9,0),
             c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
             pch = 16, cex = 1.1,
             col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
    }
  }
  
  dev.off()
  
  ## Display only 2 tissue sample types
  tissuemyT <- rbind(myT[myT$Source == "feces",], myT[myT$Source == "Cecal Content",])
  
  myMDS <- capscale(tissuemyT[,(endMetadataIndex +1):ncol(tissuemyT)]~1,distance="bray")
  
  ## This likely won't be used for anything
  write.table(myMDS$CA$u, sep="\t", file=paste(dataType, "_L_", taxa, "_cecalfecal_pcoa.txt",sep=""))
  write.table(myMDS$CA$eig, file=paste(dataType, "_L_", taxa, "_cecalfecal_eigenValues.txt", sep=""), sep="\t")
  
  percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
  
  pdf( paste(dataType, "_L_", taxa, "_cecalfecal_topMDS.pdf",sep=""))
  for (xrun in 1:5) {
    for (yrun in 2:5) {
      if(xrun == yrun){
        break
      }
      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
      ## color by tissue
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main=paste0("PCoA at level ",taxa," \n(Sample Tissue)"), cex=2.0, ##bty="L",
           pch = 16,
           col=ifelse(tissuemyT$Source=="Cecal Content", "green", ifelse(tissuemyT$Source == "duo", "black", ifelse(tissuemyT$Source == "feces", "brown", ifelse(tissuemyT$Source == "ileum", "yellow", "red")))))
      par(xpd=TRUE)
      
      legend("topright", inset=c(-0.9,0),
             c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
             pch = 16, cex = 1.1,
             col=c("green", "blue", "brown", "yellow", "red"))
      
      ## Color by Sex
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main=paste0("PCoA at level ",taxa," \n(Mouse Sex)"), cex=2.0, ##bty="L",
           pch = 16,
           col=ifelse(tissuemyT$Sex == "male", "blue", "pink"))
      par(xpd=TRUE)
      
      legend("topright", inset=c(-0.9,0),
             c("Male", "Female"),
             pch = 16, cex = 1.1,
             col=c("blue", "pink"))
      
      ## Color by stress
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main=paste0("PCoA at level ",taxa," \n(Stress/Unstressed)"), cex=2.0, ##bty="L",
           pch = 16,
           col=ifelse(tissuemyT$Exp.or.Ctrl == "Exp", "red", "black"))
      par(xpd=TRUE)
      
      legend("topright", inset=c(-0.9,0),
             c("Stress", "Control"),
             pch = 16, cex = 1.1,
             col=c("red", "black"))
      
      ## Color by cage
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main=paste0("PCoA at level ",taxa," \n(Cage)"), cex=2.0, ##bty="L",
           pch = 16,
           col=ifelse(tissuemyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(tissuemyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(tissuemyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(tissuemyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(tissuemyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(tissuemyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(tissuemyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
      par(xpd=TRUE)
      
      legend("topright", inset=c(-0.9,0),
             c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
             pch = 16, cex = 1.1,
             col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
    }
  }
  
  dev.off()
  
  tissuemyT <- myT[myT$Source == "Cecal Content",]
  ## Display only cecal tissue sample types
  myMDS <- capscale(tissuemyT[,(endMetadataIndex +1):ncol(tissuemyT)]~1,distance="bray")
  
  write.table(myMDS$CA$u, sep="\t", file=paste(dataType, "_L_", taxa, "_cecal_pcoa.txt",sep=""))
  write.table(myMDS$CA$eig,file=paste(dataType, "_L_", taxa, "_cecal_eigenValues.txt", sep=""), sep="\t")
  
  percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
  
  pdf( paste(dataType, "_L_", taxa, "cecal_topMDS.pdf",sep=""))
  for (xrun in 1:5) {
    for (yrun in 2:5) {
      if(xrun == yrun){
        break
      }
      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
      ## Color by tissue (this is useless at this point)
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main=paste0("PCoA at level ", taxa, " \n(Sample Tissue)"), cex=2.0, ##bty="L",
           pch = ifelse(tissuemyT$Exp.or.Ctrl == "Exp", 17, 16),
           col=ifelse(tissuemyT$Sex == "male", "blue", "red"))
      par(xpd=TRUE)
      
      legend("topright", inset=c(-0.9,0),
             c("Male\nControl", "Female\nControl", "Male\nStressed", "Female\nStressed"),
             pch = c(16, 16, 17, 17), cex = 1.1,
             col=c("blue", "red", "blue", "red"))
      
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main=paste0("PCoA at level ", taxa," \n(Mouse Sex)"), cex=2.0, ##bty="L",
           pch = 16,
           col=ifelse(tissuemyT$Sex == "male", "blue", "pink"))
      par(xpd=TRUE)
      
      legend("topright", inset=c(-0.9,0),
             c("Male", "Female"),
             pch = 16, cex = 1.1,
             col=c("blue", "pink"))
      
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main=paste0("PCoA at level ",taxa," \n(Stress/Unstressed)"), cex=2.0, ##bty="L",
           pch = 16,
           col=ifelse(tissuemyT$Exp.or.Ctrl == "Exp", "red", "black"))
      par(xpd=TRUE)
      
      legend("topright", inset=c(-0.9,0),
             c("Stress", "Control"),
             pch = 16, cex = 1.1,
             col=c("red", "black"))
      
      
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main=paste0("PCoA at level ",taxa," \n(Cage)"), cex=2.0, ##bty="L",
           pch = 16,
           col=ifelse(tissuemyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(tissuemyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(tissuemyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(tissuemyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(tissuemyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(tissuemyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(tissuemyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
      par(xpd=TRUE)
      
      legend("topright", inset=c(-0.9,0),
             c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
             pch = 16, cex = 1.1,
             col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
    }
  }
  
  dev.off()
  
  ## This is the one that is utimately following in the analysis
  tissuemyT <- myT[myT$Source == "feces",]
  
  myMDS <- capscale(tissuemyT[,(endMetadataIndex +1):ncol(tissuemyT)]~1,distance="bray")
  
  write.table(myMDS$CA$u, sep="\t", file=paste(dataType, "_L", taxa, "_fecal_pcoa.txt",sep=""))
  write.table(myMDS$CA$eig,file=paste(dataType, "_L", taxa, "_fecal_eigenValues.txt", sep=""), sep="\t")
  
  percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
  
  pdf( paste(dataType, "_L", taxa, "_fecal_topMDS.pdf",sep=""))
  for (xrun in 1:5) {
    for (yrun in 2:5) {
      if(xrun == yrun){
        break
      }
      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
      ## This plot is useless
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main=paste0("PCoA at level ",taxa," \n(Sample Tissue)"), cex=2.0, ##bty="L",
           pch = ifelse(tissuemyT$Exp.or.Ctrl == "Exp", 17, 16),
           col= ifelse(tissuemyT$Sex=="male", rgb(0,0,1,0.75), rgb(1,0,0,0.75)))
      par(xpd=TRUE)
      
      legend("topright", inset=c(-0.9,0),
             c("Male\nControl", "Female\nControl", "Male\nStressed", "Female\nStressed"),
             pch = c(16, 16, 17, 17), cex = 1.1,
             col=c(rgb(0,0,1,0.75), rgb(1,0,0,0.75), rgb(0,0,1,0.75), rgb(1,0,0,0.75)))
      ## Color by Sex
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main=paste0("PCoA at level ",taxa," \n(Mouse Sex)"), cex=2.0, ##bty="L",
           pch = 16,
           col=ifelse(tissuemyT$Sex == "male", "blue", "pink"))
      par(xpd=TRUE)
      
      legend("topright", inset=c(-0.9,0),
             c("Male", "Female"),
             pch = 16, cex = 1.1,
             col=c("blue", "pink"))
      ## Color by stressed
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main=paste0("PCoA at level ",taxa," \n(Stress/Unstressed)"), cex=2.0, ##bty="L",
           pch = 16,
           col=ifelse(tissuemyT$Exp.or.Ctrl == "Exp", "red", "black"))
      par(xpd=TRUE)
      
      legend("topright", inset=c(-0.9,0),
             c("Stress", "Control"),
             pch = 16, cex = 1.1,
             col=c("red", "black"))
      
      ## Color by cage
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main=paste0("PCoA at level ",taxa," \n(Cage)"), cex=2.0, ##bty="L",
           pch = 16,
           col=ifelse(tissuemyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(tissuemyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(tissuemyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(tissuemyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(tissuemyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(tissuemyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(tissuemyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
      par(xpd=TRUE)
      
      legend("topright", inset=c(-0.9,0),
             c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
             pch = 16, cex = 1.1,
             col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
    }
  }
  
  dev.off()
  
  ## Plotting in line with figure from originally submitted manuscript
  
  
  ## Over all three behavioral tests
  if (taxa == 2){
    ## It doesn't matter which tissue here
    behaviormyT <- myT[myT$Source == "feces",]
    myMDS <- capscale(behaviormyT[,22:60]~1,distance="bray")
    
    write.table(myMDS$CA$u, sep="\t", file=paste(dataType, "_L", taxa, "_behavior_pcoa.txt",sep=""))
    write.table(myMDS$CA$eig,file=paste(dataType, "_L", taxa, "_behavior_eigenValues.txt", sep=""), sep="\t")
    
    percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
    
    pdf("behavior_topMDS.pdf")
    for (xrun in 1:4) {
      for (yrun in 2:4) {
        if(xrun == yrun){
          break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main=paste0("PCoA at otu level \n(Sample Tissue)"), cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Source=="Cecal Content", "green", ifelse(behaviormyT$Source == "duo", "black", ifelse(behaviormyT$Source == "feces", "brown", ifelse(behaviormyT$Source == "ileum", "yellow", "red")))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main=paste0("PCoA at otu level \n(Mouse Sex)"), cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Sex == "male", "blue", "pink"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Male", "Female"),
               pch = 16, cex = 1.1,
               col=c("blue", "pink"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Stress/Unstressed)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Exp.or.Ctrl == "Exp", "red", "black"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Stress", "Control"),
               pch = 16, cex = 1.1,
               col=c("red", "black"))
        
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Cage)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(behaviormyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(behaviormyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(behaviormyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(behaviormyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(behaviormyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(behaviormyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
      }
    }
    
    dev.off()
    
    ## elevated.plus.maze
    behaviormyT <- myT[myT$Source == "feces",]
    ## behaviormyT <- behaviormyT[behaviormyT$Sex == "male",]
    myMDS <- capscale(behaviormyT[,22:32]~1,distance="bray")
    
    write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_epm_pcoa.txt",sep=""))
    write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_epm_eigenValues.txt", sep=""), sep="\t")
    
    percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
    
    pdf("elevatedplusmaze_topMDS.pdf")
    for (xrun in 1:4) {
      for (yrun in 2:4) {
        if(xrun == yrun){
          break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Sample Tissue)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Source=="Cecal Content", "green", ifelse(behaviormyT$Source == "duo", "black", ifelse(behaviormyT$Source == "feces", "brown", ifelse(behaviormyT$Source == "ileum", "yellow", "red")))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Mouse Sex)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Sex == "male", "blue", "pink"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Male", "Female"),
               pch = 16, cex = 1.1,
               col=c("blue", "pink"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Stress/Unstressed)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Exp.or.Ctrl == "Exp", "red", "black"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Stress", "Control"),
               pch = 16, cex = 1.1,
               col=c("red", "black"))
        
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Cage)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(behaviormyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(behaviormyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(behaviormyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(behaviormyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(behaviormyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(behaviormyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
      }
    }
    
    dev.off()
    behaviormyT <- myT[myT$Source == "feces",]
    behaviormyT <- behaviormyT[behaviormyT$Sex == "male",]
    myMDS <- capscale(behaviormyT[,22:32]~1,distance="bray")
    
    write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_epm_pcoa_MALE.txt",sep=""))
    write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_epm_eigenValues_MALE.txt", sep=""), sep="\t")
    
    percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
    
    pdf("elevatedplusmaze_topMDS_MALE.pdf")
    for (xrun in 1:4) {
      for (yrun in 2:4) {
        if(xrun == yrun){
          break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Sample Tissue)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Source=="Cecal Content", "green", ifelse(behaviormyT$Source == "duo", "black", ifelse(behaviormyT$Source == "feces", "brown", ifelse(behaviormyT$Source == "ileum", "yellow", "red")))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Mouse Sex)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Sex == "male", "blue", "pink"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Male", "Female"),
               pch = 16, cex = 1.1,
               col=c("blue", "pink"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Stress/Unstressed)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Exp.or.Ctrl == "Exp", "red", "black"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Stress", "Control"),
               pch = 16, cex = 1.1,
               col=c("red", "black"))
        
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Cage)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(behaviormyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(behaviormyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(behaviormyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(behaviormyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(behaviormyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(behaviormyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
      }
    }
    
    dev.off()
    
    
    behaviormyT <- myT[myT$Source == "feces",]
    behaviormyT <- behaviormyT[behaviormyT$Sex == "female",]
    myMDS <- capscale(behaviormyT[, 22:32]~1,distance="bray")
    write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_epm_pcoa_FEMALE.txt",sep=""))
    write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_epm_eigenValues_FEMALE.txt", sep=""), sep="\t")
    
    percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
    
    pdf("elevatedplusmaze_FEMALE_topMDS.pdf")
    for (xrun in 1:4) {
      for (yrun in 2:4) {
        if(xrun == yrun){
          break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Sample Tissue)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Source=="Cecal Content", "green", ifelse(behaviormyT$Source == "duo", "black", ifelse(behaviormyT$Source == "feces", "brown", ifelse(behaviormyT$Source == "ileum", "yellow", "red")))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Mouse Sex)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Sex == "male", "blue", "pink"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Male", "Female"),
               pch = 16, cex = 1.1,
               col=c("blue", "pink"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Stress/Unstressed)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Exp.or.Ctrl == "Exp", "red", "black"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Stress", "Control"),
               pch = 16, cex = 1.1,
               col=c("red", "black"))
        
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Cage)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(behaviormyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(behaviormyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(behaviormyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(behaviormyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(behaviormyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(behaviormyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
      }
    }
    
    dev.off()
    
    
    ## ldmyT <- behaviormyT[,12:21]
    ## Display all 5 tissue sample types
    behaviormyT <- myT[myT$Source == "feces",]
    myMDS <- capscale(behaviormyT[,33:42]~1,distance="bray")
    
    write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_lightdark_pcoa.txt",sep=""))
    write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_lightdark_eigenValues.txt", sep=""), sep="\t")
    
    percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
    
    pdf("lightdark_topMDS.pdf")
    for (xrun in 1:4) {
      for (yrun in 2:4) {
        if(xrun == yrun){
          break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Sample Tissue)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Source=="Cecal Content", "green", ifelse(behaviormyT$Source == "duo", "black", ifelse(behaviormyT$Source == "feces", "brown", ifelse(behaviormyT$Source == "ileum", "yellow", "red")))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Mouse Sex)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Sex == "male", "blue", "pink"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Male", "Female"),
               pch = 16, cex = 1.1,
               col=c("blue", "pink"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Stress/Unstressed)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Exp.or.Ctrl == "Exp", "red", "black"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Stress", "Control"),
               pch = 16, cex = 1.1,
               col=c("red", "black"))
        
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Cage)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(behaviormyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(behaviormyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(behaviormyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(behaviormyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(behaviormyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(behaviormyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
      }
    }
    
    dev.off()
    
    behaviormyT <- myT[myT$Source == "feces",]
    behaviormyT <- behaviormyT[behaviormyT$Sex == "male",]
    myMDS <- capscale(behaviormyT[,33:42]~1,distance="bray")
    
    write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_lightdark_pcoa_MALE.txt",sep=""))
    write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_lightdark_eigenValues_MALE.txt", sep=""), sep="\t")
    
    percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
    
    pdf("lightdark_topMDS_MALE.pdf")
    for (xrun in 1:4) {
      for (yrun in 2:4) {
        if(xrun == yrun){
          break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Sample Tissue)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Source=="Cecal Content", "green", ifelse(behaviormyT$Source == "duo", "black", ifelse(behaviormyT$Source == "feces", "brown", ifelse(behaviormyT$Source == "ileum", "yellow", "red")))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Mouse Sex)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Sex == "male", "blue", "pink"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Male", "Female"),
               pch = 16, cex = 1.1,
               col=c("blue", "pink"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Stress/Unstressed)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Exp.or.Ctrl == "Exp", "red", "black"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Stress", "Control"),
               pch = 16, cex = 1.1,
               col=c("red", "black"))
        
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Cage)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(behaviormyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(behaviormyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(behaviormyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(behaviormyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(behaviormyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(behaviormyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
      }
    }
    
    dev.off()
    
    behaviormyT <- myT[myT$Source == "feces",]
    behaviormyT <- behaviormyT[behaviormyT$Sex == "female",]
    myMDS <- capscale(behaviormyT[,33:42]~1,distance="bray")
    
    write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_lightdark_pcoa_FEMALE.txt",sep=""))
    write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_lightdark_eigenValues_FEMALE.txt", sep=""), sep="\t")
    
    percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
    
    pdf("lightdark_topMDS_FEMALE.pdf")
    for (xrun in 1:4) {
      for (yrun in 2:4) {
        if(xrun == yrun){
          break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Sample Tissue)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Source=="Cecal Content", "green", ifelse(behaviormyT$Source == "duo", "black", ifelse(behaviormyT$Source == "feces", "brown", ifelse(behaviormyT$Source == "ileum", "yellow", "red")))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Mouse Sex)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Sex == "male", "blue", "pink"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Male", "Female"),
               pch = 16, cex = 1.1,
               col=c("blue", "pink"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Stress/Unstressed)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Exp.or.Ctrl == "Exp", "red", "black"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Stress", "Control"),
               pch = 16, cex = 1.1,
               col=c("red", "black"))
        
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Cage)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(behaviormyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(behaviormyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(behaviormyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(behaviormyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(behaviormyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(behaviormyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
      }
    }
    
    dev.off()
    
    ## ofmyT <- behaviormyT[,22:39]
    ## Display all 5 tissue sample types
    behaviormyT <- myT[myT$Source == "feces",]
    myMDS <- capscale(behaviormyT[,43:60]~1,distance="bray")
    
    write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_openfield_pcoa.txt",sep=""))
    write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_openfield_eigenValues.txt", sep=""), sep="\t")
    
    percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
    
    pdf("openfield_topMDS.pdf")
    for (xrun in 1:4) {
      for (yrun in 2:4) {
        if(xrun == yrun){
          break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Sample Tissue)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Source=="Cecal Content", "green", ifelse(behaviormyT$Source == "duo", "black", ifelse(behaviormyT$Source == "feces", "brown", ifelse(behaviormyT$Source == "ileum", "yellow", "red")))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Mouse Sex)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Sex == "male", "blue", "pink"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Male", "Female"),
               pch = 16, cex = 1.1,
               col=c("blue", "pink"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Stress/Unstressed)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Exp.or.Ctrl == "Exp", "red", "black"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Stress", "Control"),
               pch = 16, cex = 1.1,
               col=c("red", "black"))
        
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Cage)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(behaviormyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(behaviormyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(behaviormyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(behaviormyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(behaviormyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(behaviormyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
      }
    }
    
    dev.off()
    behaviormyT <- myT[myT$Source == "feces",]
    behaviormyT <- behaviormyT[behaviormyT$Sex == "male",]
    myMDS <- capscale(behaviormyT[,43:60]~1,distance="bray")
    
    write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_openfield_pcoa_MALE.txt",sep=""))
    write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_openfield_eigenValues_MALE.txt", sep=""), sep="\t")
    
    percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
    
    pdf("openfield_topMDS_MALE.pdf")
    for (xrun in 1:4) {
      for (yrun in 2:4) {
        if(xrun == yrun){
          break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Sample Tissue)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Source=="Cecal Content", "green", ifelse(behaviormyT$Source == "duo", "black", ifelse(behaviormyT$Source == "feces", "brown", ifelse(behaviormyT$Source == "ileum", "yellow", "red")))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Mouse Sex)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Sex == "male", "blue", "pink"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Male", "Female"),
               pch = 16, cex = 1.1,
               col=c("blue", "pink"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Stress/Unstressed)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Exp.or.Ctrl == "Exp", "red", "black"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Stress", "Control"),
               pch = 16, cex = 1.1,
               col=c("red", "black"))
        
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Cage)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(behaviormyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(behaviormyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(behaviormyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(behaviormyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(behaviormyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(behaviormyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
      }
    }
    
    dev.off()
    
    behaviormyT <- myT[myT$Source == "feces",]
    behaviormyT <- behaviormyT[behaviormyT$Sex == "female",]
    myMDS <- capscale(behaviormyT[,43:60]~1,distance="bray")
    
    write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_openfield_pcoa_FEMALE.txt",sep=""))
    write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_openfield_eigenValues_FEMALE.txt", sep=""), sep="\t")
    
    percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
    
    pdf("openfield_topMDS_FEMALE.pdf")
    for (xrun in 1:4) {
      for (yrun in 2:4) {
        if(xrun == yrun){
          break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Sample Tissue)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Source=="Cecal Content", "green", ifelse(behaviormyT$Source == "duo", "black", ifelse(behaviormyT$Source == "feces", "brown", ifelse(behaviormyT$Source == "ileum", "yellow", "red")))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Mouse Sex)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Sex == "male", "blue", "pink"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Male", "Female"),
               pch = 16, cex = 1.1,
               col=c("blue", "pink"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Stress/Unstressed)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Exp.or.Ctrl == "Exp", "red", "black"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Stress", "Control"),
               pch = 16, cex = 1.1,
               col=c("red", "black"))
        
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Cage)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(behaviormyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(behaviormyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(behaviormyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(behaviormyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(behaviormyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(behaviormyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
      }
    }
    
    dev.off()
    
  }
  ## Compute MDS on behavioral data and output
  
  ## Plotting for presentation
  tissuemyT <- myT[myT$Source == "feces",]
  
  myMDS <- capscale(tissuemyT[,(endMetadataIndex +1):ncol(tissuemyT)]~1,distance="bray")
  
  percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
  
  pdf( paste(dataType, "_L", taxa, "_PRESENTATION_topMDS.pdf",sep=""), height=9, width = 16)
  
  for (xrun in 1:4) {
    for (yrun in 2:4) {
      if(xrun == yrun){
        break
      }
      par(mar=c(2, 2, 3, 1) + 0.1,
          oma=c(4, 4, 1, 1) + 0.1,
          mgp=c(3, 1, 0),
          xpd=TRUE, mfrow = c(1,3))
      
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun], ylim = c(-0.4, 0.4),
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main="Coloration by Sex of Mouse", cex=3.0, cex.main = 2.5, ##bty="L",
           pch = 16,
           col=ifelse(tissuemyT$Sex == "male", "blue", "pink"), cex.axis = 2)
      ##             par(xpd=TRUE)
      
      legend("topright",
             c("Male", "Female"),
             pch = 16, cex = 2,
             col=c("blue", "pink"))
      mtext(paste("MDS Axis", yrun, ": ", format(percentVariance[yrun],digits=3),sep=""), side = 2, line = 1.5, outer=TRUE, cex=2)
      
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun], ylim = c(-0.4, 0.4),
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main="Coloration by Stress Versus Control", cex=3.0, cex.main = 2.5, ##bty="L",
           pch = 16,
           col=ifelse(tissuemyT$Exp.or.Ctrl == "Exp", "red", "black"), cex.axis = 2)
      ##                     par(xpd=TRUE)
      
      legend("topright",
             c("Stress", "Control"),
             pch = 16, cex = 2,
             col=c("red", "black"))
      
      mtext(paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""), line = 1.5, side = 1, outer=TRUE, cex=2)
      
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun], ylim = c(-0.4, 0.4),
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main="Coloration by Cage", cex.main = 2.5, cex = 3, ##bty="L",
           pch = 16,
           col=ifelse(tissuemyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(tissuemyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(tissuemyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(tissuemyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(tissuemyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(tissuemyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(tissuemyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))), cex.axis=2)
      par(xpd=TRUE)
      
      legend("topright",
             c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
             pch = 16, cex = 1.7,
             col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
      
    }
  }
  
  dev.off()
  
  if (taxa == 2) {
    ## Plotting for presentation
    tissuemyT <- myT[myT$Source == "feces",]
    ## For behavioral abundance data
    myMDS <- capscale(tissuemyT[,22:60]~1,distance="bray")
    
    percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
    
    pdf( paste("PRESENTATION_BEHAVIOR_topMDS.pdf",sep=""), height=9, width = 16)
    
    for (xrun in 1:4) {
      for (yrun in 2:4) {
        if(xrun == yrun){
          break
        }
        par(mar=c(2, 2, 3, 1) + 0.1,
            oma=c(4, 4, 1, 1) + 0.1,
            mgp=c(3, 1, 0),
            xpd=TRUE, mfrow = c(1,3))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="Coloration by Sex of Mouse", cex=3.0, cex.main = 2.5, ##bty="L",
             pch = 16,
             col=ifelse(tissuemyT$Sex == "male", "blue", "pink"), cex.axis = 2)
        ##             par(xpd=TRUE)
        
        legend("topright",
               c("Male", "Female"),
               pch = 16, cex = 2,
               col=c("blue", "pink"))
        mtext(paste("MDS Axis", yrun, ": ", format(percentVariance[yrun],digits=3),sep=""), side = 2, line = 1.5, outer=TRUE, cex=2)
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="Coloration by Stress Versus Control", cex=3.0, cex.main = 2.5, ##bty="L",
             pch = 16,
             col=ifelse(tissuemyT$Exp.or.Ctrl == "Exp", "red", "black"), cex.axis = 2)
        ##                     par(xpd=TRUE)
        
        legend("topright",
               c("Stress", "Control"),
               pch = 16, cex = 2,
               col=c("red", "black"))
        
        mtext(paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""), line = 1.5, side = 1, outer=TRUE, cex=2)
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="Coloration by Cage", cex.main = 2.5, cex = 3, ##bty="L",
             pch = 16,
             col=ifelse(tissuemyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(tissuemyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(tissuemyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(tissuemyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(tissuemyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(tissuemyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(tissuemyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))), cex.axis=2)
        par(xpd=TRUE)
        
        legend("topright",
               c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
               pch = 16, cex = 1.7,
               col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
        
      }
    }
    
    dev.off()
  }
  ## Plotting for presentation
  tissuemyT <- myT[myT$Source == "feces",]
  tissuemyT <- tissuemyT[tissuemyT$Sex == "male",]
  ## Display only 2 tissue sample types
  myMDS <- capscale(tissuemyT[,(endMetadataIndex +1):ncol(tissuemyT)]~1,distance="bray")
  
  percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
  
  write.table(myMDS$CA$u, sep="\t", file=paste0(dataType, "_L_", taxa, "_MALEONLYpcoa.txt"))
  write.table(myMDS$CA$eig, file=paste0(dataType, "_L_", taxa, "_MALEONLYeigenValues.txt"), sep="\t")
  
  pdf( paste(dataType, "_L", taxa, "_MALEONLY_topMDS.pdf",sep=""), height=9, width = 16)
  
  for (xrun in 1:4) {
    for (yrun in 2:4) {
      if(xrun == yrun){
        break
      }
      par(mar=c(2, 2, 3, 1) + 0.1,
          oma=c(4, 4, 1, 1) + 0.1,
          mgp=c(3, 1, 0),
          xpd=TRUE, mfrow = c(1,3))
      
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main="Coloration by Sex of Mouse", cex=3.0, cex.main = 2.5, ##bty="L",
           pch = 16,
           col=ifelse(myT$Sex == "male", "blue", "pink"), cex.axis = 2)
      ##             par(xpd=TRUE)
      
      legend("topright",
             c("Male", "Female"),
             pch = 16, cex = 2,
             col=c("blue", "pink"))
      mtext(paste("MDS Axis", yrun, ": ", format(percentVariance[yrun],digits=3),sep=""), side = 2, line = 1.5, outer=TRUE, cex=2)
      
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main="Coloration by Stress Versus Control", cex=3.0, cex.main = 2.5, ##bty="L",
           pch = 16,
           col=ifelse(myT$Exp.or.Ctrl == "Exp", "red", "black"), cex.axis = 2)
      ##                     par(xpd=TRUE)
      
      legend("topright",
             c("Stress", "Control"),
             pch = 16, cex = 2,
             col=c("red", "black"))
      
      mtext(paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""), line = 1.5, side = 1, outer=TRUE, cex=2)
      
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main="Coloration by Cage", cex.main = 2.5, cex = 3, ##bty="L",
           pch = 16,
           col=ifelse(myT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(myT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(myT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(myT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(myT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(myT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(myT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))), cex.axis=2)
      par(xpd=TRUE)
      
      legend("topright",
             c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
             pch = 16, cex = 1.7,
             col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
      
    }
  }
  
  dev.off()
  
  ## Plotting for presentation
  tissuemyT <- myT[myT$Source == "feces",]
  tissuemyT <- tissuemyT[tissuemyT$Sex == "female",]
  ## Display only 2 tissue sample types
  myMDS <- capscale(tissuemyT[,(endMetadataIndex +1):ncol(tissuemyT)]~1,distance="bray")
  
  percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
  write.table(myMDS$CA$u, sep="\t", file=paste0(dataType, "_L_", taxa, "_FEMALEONLYpcoa.txt"))
  write.table(myMDS$CA$eig, file=paste0(dataType, "_L_", taxa, "_FEMALEONLYeigenValues.txt"), sep="\t")
  
  
  pdf( paste(dataType, "_L", taxa, "_FEMALEONLY_topMDS.pdf",sep=""), height=9, width = 16)
  
  for (xrun in 1:4) {
    for (yrun in 2:4) {
      if(xrun == yrun){
        break
      }
      par(mar=c(2, 2, 3, 1) + 0.1,
          oma=c(4, 4, 1, 1) + 0.1,
          mgp=c(3, 1, 0),
          xpd=TRUE, mfrow = c(1,3))
      
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main="Coloration by Sex of Mouse", cex=3.0, cex.main = 2.5, ##bty="L",
           pch = 16,
           col=ifelse(tissuemyT$Sex == "male", "blue", "pink"), cex.axis = 2)
      ##             par(xpd=TRUE)
      
      legend("topright",
             c("Male", "Female"),
             pch = 16, cex = 2,
             col=c("blue", "pink"))
      mtext(paste("MDS Axis", yrun, ": ", format(percentVariance[yrun],digits=3),sep=""), side = 2, line = 1.5, outer=TRUE, cex=2)
      
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main="Coloration by Stress Versus Control", cex=3.0, cex.main = 2.5, ##bty="L",
           pch = 16,
           col=ifelse(tissuemyT$Exp.or.Ctrl == "Exp", "red", "black"), cex.axis = 2)
      ##                     par(xpd=TRUE)
      
      legend("topright",
             c("Stress", "Control"),
             pch = 16, cex = 2,
             col=c("red", "black"))
      
      mtext(paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""), line = 1.5, side = 1, outer=TRUE, cex=2)
      
      plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
           xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
           ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
           main="Coloration by Cage", cex.main = 2.5, cex = 3, ##bty="L",
           pch = 16,
           col=ifelse(tissuemyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(tissuemyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(tissuemyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(tissuemyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(tissuemyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(tissuemyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(tissuemyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))), cex.axis=2)
      par(xpd=TRUE)
      
      legend("topright",
             c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
             pch = 16, cex = 1.7,
             col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
      
    }
  }
  
  dev.off()
  
  if (taxa == 2) {
    ## Plotting for presentation
    tissuemyT <- myT[myT$Source == "feces",]
    tissuemyT <- tissuemyT[tissuemyT$Sex == "male",]
    ## Display only 2 tissue sample types
    myMDS <- capscale(tissuemyT[,22:60]~1,distance="bray")
    
    percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
    
    write.table(myMDS$CA$u, sep="\t", file=paste0(dataType, "_L_", taxa, "_MALEONLYBEHAVIORpcoa.txt"))
    write.table(myMDS$CA$eig, file=paste0(dataType, "_L_", taxa, "_MALEONLYBEHAVIOReigenValues.txt"), sep="\t")
    
    
    pdf( paste("MALEONLYBEHAVIOR_topMDS.pdf",sep=""), height=9, width = 16)
    
    for (xrun in 1:4) {
      for (yrun in 2:4) {
        if(xrun == yrun){
          break
        }
        par(mar=c(2, 2, 3, 1) + 0.1,
            oma=c(4, 4, 1, 1) + 0.1,
            mgp=c(3, 1, 0),
            xpd=TRUE, mfrow = c(1,3))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="Coloration by Sex of Mouse", cex=3.0, cex.main = 2.5, ##bty="L",
             pch = 16,
             col=ifelse(myT$Sex == "male", "blue", "pink"), cex.axis = 2)
        ##             par(xpd=TRUE)
        
        legend("topright",
               c("Male", "Female"),
               pch = 16, cex = 2,
               col=c("blue", "pink"))
        mtext(paste("MDS Axis", yrun, ": ", format(percentVariance[yrun],digits=3),sep=""), side = 2, line = 1.5, outer=TRUE, cex=2)
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="Coloration by Stress Versus Control", cex=3.0, cex.main = 2.5, ##bty="L",
             pch = 16,
             col=ifelse(myT$Exp.or.Ctrl == "Exp", "red", "black"), cex.axis = 2)
        ##                     par(xpd=TRUE)
        
        legend("topright",
               c("Stress", "Control"),
               pch = 16, cex = 2,
               col=c("red", "black"))
        
        mtext(paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""), line = 1.5, side = 1, outer=TRUE, cex=2)
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="Coloration by Cage", cex.main = 2.5, cex = 3, ##bty="L",
             pch = 16,
             col=ifelse(myT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(myT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(myT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(myT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(myT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(myT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(myT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))), cex.axis=2)
        par(xpd=TRUE)
        
        legend("topright",
               c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
               pch = 16, cex = 1.7,
               col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
        
      }
    }
    
    dev.off()
    
    ## Plotting for presentation
    tissuemyT <- myT[myT$Source == "feces",]
    tissuemyT <- tissuemyT[tissuemyT$Sex == "female",]
    ## Display only 2 tissue sample types
    myMDS <- capscale(tissuemyT[,22:60]~1,distance="bray")
    
    percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
    
    write.table(myMDS$CA$u, sep="\t", file=paste0(dataType, "_L_", taxa, "_FEMALEONLYBEHAVIORpcoa.txt"))
    write.table(myMDS$CA$eig, file=paste0(dataType, "_L_", taxa, "_FEMALEONLYBEHAVIOReigenValues.txt"), sep="\t")
    
    pdf( paste("FEMALEONLYBEHAVIOR_topMDS.pdf",sep=""), height=9, width = 16)
    
    for (xrun in 1:4) {
      for (yrun in 2:4) {
        if(xrun == yrun){
          break
        }
        par(mar=c(2, 2, 3, 1) + 0.1,
            oma=c(4, 4, 1, 1) + 0.1,
            mgp=c(3, 1, 0),
            xpd=TRUE, mfrow = c(1,3))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="Coloration by Sex of Mouse", cex=3.0, cex.main = 2.5, ##bty="L",
             pch = 16,
             col=ifelse(tissuemyT$Sex == "male", "blue", "pink"), cex.axis = 2)
        ##             par(xpd=TRUE)
        
        legend("topright",
               c("Male", "Female"),
               pch = 16, cex = 2,
               col=c("blue", "pink"))
        mtext(paste("MDS Axis", yrun, ": ", format(percentVariance[yrun],digits=3),sep=""), side = 2, line = 1.5, outer=TRUE, cex=2)
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="Coloration by Stress Versus Control", cex=3.0, cex.main = 2.5, ##bty="L",
             pch = 16,
             col=ifelse(tissuemyT$Exp.or.Ctrl == "Exp", "red", "black"), cex.axis = 2)
        ##                     par(xpd=TRUE)
        
        legend("topright",
               c("Stress", "Control"),
               pch = 16, cex = 2,
               col=c("red", "black"))
        
        mtext(paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""), line = 1.5, side = 1, outer=TRUE, cex=2)
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="Coloration by Cage", cex.main = 2.5, cex = 3, ##bty="L",
             pch = 16,
             col=ifelse(tissuemyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(tissuemyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(tissuemyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(tissuemyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(tissuemyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(tissuemyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(tissuemyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))), cex.axis=2)
        par(xpd=TRUE)
        
        legend("topright",
               c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
               pch = 16, cex = 1.7,
               col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
        
      }
    }
    
    dev.off()
  }
  
  if (taxa == 2){
    
    behaviormyT <- myT[myT$Source == "feces",]
    
    myMDS <- capscale(behaviormyT[,22:60]~1,distance="bray")
    
    write.table(myMDS$CA$u, sep="\t", file=paste("BehaviorBIPLOT_pcoa.txt",sep=""))
    write.table(myMDS$CA$eig,file=paste("BehaviorBIPLOT_eigenValues.txt", sep=""), sep="\t")
    
    percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
    
    plot(myMDS, display=c("wa", "bp"))
    plot(myMDS, display=c("sp", "bp"))
    plot(myMDS, display=c("wa", "sp", "bp"))
    
    pdf("BehaviorBIPLOT_topMDS.pdf")
    for (xrun in 1:4) {
      for (yrun in 2:4) {
        if(xrun == yrun){
          break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
        ## This one is not useful
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA colored by \n(Sample Tissue)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Source=="Cecal Content", "green", ifelse(behaviormyT$Source == "duo", "black", ifelse(behaviormyT$Source == "feces", "brown", ifelse(behaviormyT$Source == "ileum", "yellow", "red")))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA Colored by \n(Mouse Sex)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Sex == "male", "blue", "pink"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Male", "Female"),
               pch = 16, cex = 1.1,
               col=c("blue", "pink"))
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Stress/Unstressed)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Exp.or.Ctrl == "Exp", "red", "black"))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("Stress", "Control"),
               pch = 16, cex = 1.1,
               col=c("red", "black"))
        
        
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA Colored by \n(Cage)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(behaviormyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(behaviormyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(behaviormyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(behaviormyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(behaviormyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(behaviormyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
        par(xpd=TRUE)
        
        legend("topright", inset=c(-0.9,0),
               c("SF1", "CF2", "CF1", "SF2", "SM1", "CM2", "CM1", "SM2"),
               pch = 16, cex = 1.1,
               col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
      }
    }
    
    dev.off()
  }
  
}

tissuemyT <- myT[myT$Source == "Cecal Content",]

myMDS <- capscale(tissuemyT[,(endMetadataIndex +1):ncol(tissuemyT)]~1,distance="bray")

write.table(myMDS$CA$u, sep="\t", file=paste(dataType, "_L", taxa, "_cecal_pcoa.txt",sep=""))
write.table(myMDS$CA$eig,file=paste(dataType, "_L", taxa, "_cecal_eigenValues.txt", sep=""), sep="\t")

percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))

pdf( paste(dataType, "_L", taxa, "_OnePerPagececal_topMDS.pdf",sep=""))
for (xrun in 1:5) {
  for (yrun in 2:5) {
    if(xrun == yrun){
      break
    }
    #par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
    plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
         xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
         ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
         main=paste0("PCoA at level ",taxa," \n(Sample Tissue)"), cex=2.0, ##bty="L",
         pch = ifelse(tissuemyT$Exp.or.Ctrl == "Exp", 17, 16),
         col= ifelse(tissuemyT$Sex=="male", rgb(0,0,1,0.75), rgb(1,0,0,0.75)))
    par(xpd=TRUE)
    
    legend("topright", inset=c(-0.9,0),
           c("Male\nControl", "Female\nControl", "Male\nStressed", "Female\nStressed"),
           pch = c(16, 16, 17, 17), cex = 1.1,
           col=c(rgb(0,0,1,0.75), rgb(1,0,0,0.75), rgb(0,0,1,0.75), rgb(1,0,0,0.75)))
  }
}

dev.off()
