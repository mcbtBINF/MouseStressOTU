# Corticosterone Modeling

# Attach the data to the existing OTUTable
# Later add in code for the taxa-by-taxa level modeling
rm(list=ls())

library("Kendall")
library("pscl")
library("lmtest")
library("nlme")
library("psych")
library("ggplot2")
library("vegan")
library("coin")

otuTable <- FALSE
dataType <- "closedQIIMER1"
tissues <- c("feces")
taxa <- 2
## dataType <- "denovoabundantOTU"
if (.Platform$OS.type == "windows") {
  if (dataType == "closedQIIMER1") {
    baseDir <- "F:/Google Drive/Datasets/MarkExperiment/closedQIIMER1/"
  } else {
    baseDir <- "/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/ArgonneSequencing/JanuaryResequencing/rawData/rg_results/abundantOTU"
  }
} else {
  if (dataType == "closedQIIMER1") {
    baseDir <- "/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/ArgonneSequencing/JanuaryResequencing/rawData/rg_results/closedQIIMER1"
  } else {
    baseDir <- "/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/ArgonneSequencing/JanuaryResequencing/rawData/rg_results/abundantOTU"
  }
}

metadataDir <- paste(baseDir, "metadata", sep="/")
dataDir <- paste(baseDir, "data", sep="/")
processedDir <- paste(baseDir, "processed", sep="/")
analysisDir <- paste(baseDir, "analysis", sep="/")

if( otuTable == TRUE) {
  if ( dataType == "closedQIIMER1") {
    otuMapping <- read.table(paste(processedDir,"LyteSharon_r01_cr_MAPPING.txt",sep="/"), sep="\t", header=TRUE)
  } else {
    otuMapping <- read.table(paste(processedDir,"LyteSharon_R01_PL_wTaxaRDP80_MAPPING.txt", sep="/"), sep="\t", header=TRUE)
  }
  otuMapping[,1] <- paste0("X",otuMapping[,1])
  rownames(otuMapping)<-otuMapping[,1]
}



## Add in the corticosterone data
## TODO Make this nicer later
divider <- 4
if (otuTable == TRUE) {
  taxonomicLevels <- 1
} else {
  taxonomicLevels <- c(2:7)
}

bestAt <- data.frame()
index <- 1
setwd(processedDir)
if( otuTable == TRUE) {
  if ( dataType == "closedQIIMER1") {
    inFileName <- "qiimeOTULogNormwithMetadata.txt"
  } else {
    inFileName <- "abundantOTUR1LogNormwithMetadata.txt"
  }
} else {
  inFileName <- paste0(dataType, "_L", taxa, "_LogNormwithMetadata.txt")
}

myT <- read.table(inFileName, header=TRUE, sep="\t", comment.char="@")

# Add in values for the new measurements
setwd(metadataDir)
corticosteroneData <- read.table("Sample_corticosterone.txt", header=TRUE, sep = "\t", comment.char="@")
testMerge <- merge(x=myT, y=corticosteroneData, by = 1, all=TRUE)
myT <- testMerge
## Some concerns about column classes
setwd(analysisDir)
if(otuTable == TRUE) {
  endMetadataIndex <- which(colnames(myT) == "counts")
} else {
  endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
}

myT <- myT[, c(1:endMetadataIndex, ncol(myT), (endMetadataIndex+1):(ncol(myT) - 1))]


subT <- myT[myT$Source == tissues,]
t <- tissues
names <- vector()
pValuesGroup <- vector()
pValuesSex <- vector()
pValuesCorticosterone <- vector()
pVals <- list()
names <- vector()

names[index] <- "Corticosterone"
corticosterone <- subT$AvgCorticosterone.pg.ml.
sex <- factor(subT$Sex)

cage <-  factor(subT$Cage, c("female #5,6,7,8 in same cage", "female #13,14,15,16 in same cage", "female #1,2,3,4 in same cage", "female #9,10,11,12 in same cage",
                             "male #5,6,7,8 in same cage", "male #13,14,15,16 in same cage", "male #1,2,3,4 in same cage", "male #9,10,11,12 in same cage"))

group <- factor(subT$Group)
## Depth was found to not be significant
## This may need to be taken out for the otuTable option
depth <- as.numeric(subT$counts)

myFrame <- data.frame(corticosterone, sex, group, cage)

pValuesSex[index] <- pvalue(wilcox_test(corticosterone~sex, data=myFrame))
pValuesGroup[index] <- pvalue(wilcox_test(corticosterone~group, data=myFrame))


fullModel <-
  gls( corticosterone~  group * sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),data = myFrame )

glsModelNoCage <- gls(corticosterone~  group * sex, data = myFrame)
#print( rownames(anova(glsModelNoCage))[-1])
## Merges the nonparametric with the terms of the MLM, along with a cage effect p-value calculation
pVals[[index]] <- c(pValuesSex[index], pValuesGroup[index], ## pValuesBehavior[index],
                    anova(fullModel)$"p-value"[-1],
                    ##summary(fullModelLME)$tTable[-1,5],
                    anova(fullModel, glsModelNoCage)$"p-value"[2])

# This may have issues due to the 1 measurement
dFrame <- data.frame(names, do.call(rbind, pVals))
colnames(dFrame) <- c("taxonName",
                      "wilcoxonSex", "wilcoxonStress", ## "kendallBehavior",
                      rownames(anova(glsModelNoCage))[-1], "cage")
adjOffset <- 2
numCols <- ncol(dFrame)

for (i in adjOffset:numCols){
  dFrame[,i + numCols - (adjOffset - 1)] <- p.adjust(dFrame[,i], method = "BH")
  colnames(dFrame)[i + numCols - (adjOffset - 1)] <- paste0("adjusted_",colnames(dFrame)[i])
}

write.table(dFrame, file=paste(dataType,"_analysis_CorticosteroneStressSexCage_gls_L", taxa, "_tissue_",t, "_divider_", divider, ".txt",sep=""), sep="\t",row.names=FALSE)

pdf( paste(dataType, "_CorticosteroneStressSexCage_gls_L", taxa, "_tissue_", t, "_divider_", divider, ".pdf", sep=""))

## This might break for the otuTable selection
if(otuTable == TRUE) {
  toIter <- rownames(dFrame)
} else {
  toIter <- 1:dim(dFrame)[1]
}
for (nameIter in toIter) {
  par(mfrow=c(2,2), mgp=c(3, 1.25, 0), oma = c(0, 0, 2, 0))
  #corticosterone <- subT[,as.character(dFrame[nameIter,]$taxonName)]
  corticosterone <- subT$AvgCorticosterone.pg.ml.
  
  sex <- factor(subT$Sex)
  if(otuTable == TRUE) {
    depth <- as.numeric(subT$counts)
  } else {
    depth <- as.numeric(subT$sampleDepth)
  }
  group <- factor(subT$Group)
  
  cage <-  factor( paste( subT$Housing, subT$Sex, sep=""), c("#5,6,7,8 in same cagefemale", "#13,14,15,16 in same cagefemale", "#1,2,3,4 in same cagefemale", "#9,10,11,12 in same cagefemale",
                                                             "#5,6,7,8 in same cagemale", "#13,14,15,16 in same cagemale", "#1,2,3,4 in same cagemale", "#9,10,11,12 in same cagemale"))
  
  myFrame <- data.frame(corticosterone, sex, group, depth, cage)
  
  boxplot( corticosterone ~ group, data = myFrame, xaxt="n", las=2,
           main = paste("Wilcoxon Stress\np-value", format(dFrame[nameIter,]$"adjusted_wilcoxonStress",digits=3) ) )
  axis(1, at=c(1, 2), labels = c("Control", "Stress"), cex.axis=0.9)
  stripchart(corticosterone ~ group,
             data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
  
  boxplot( corticosterone ~ sex, xaxt="n", las=2,
           main = paste("Wilcoxon Sex\np-value", format(dFrame[nameIter,]$"adjusted_wilcoxonSex",digits=3) ) )
  axis(1, at=c(1, 2), labels = c("Female", "Male"), cex.axis=0.9)
  stripchart(corticosterone ~ sex,
             data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
  
  plot.new()
  plot.new()
  
  mtext(paste0(dFrame[nameIter,]$taxonName), outer=TRUE, cex = 0.7)
  
  boxplot( corticosterone ~ group, las=2, xaxt="n",
           main = paste("Stress\np-value", format(dFrame[nameIter,]$"adjusted_group",digits=3) ) )
  axis(1, at=c(1, 2), labels = c("Control", "Stress"), cex.axis=0.9)
  stripchart(corticosterone ~ group,
             data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
  
  boxplot( corticosterone ~ sex, las=2, xaxt = "n",
           main = paste("Sex\np-value", format(dFrame[nameIter,]$"adjusted_sex",digits=3) ) )
  axis(1, at=c(1, 2), labels = c("Female", "Male"), cex.axis=0.9)
  stripchart(corticosterone ~ sex,
             data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
  
  mtext(paste0(dFrame[nameIter,]$taxonName), outer=TRUE, cex = 0.7)
  
  boxplot( corticosterone ~ group*sex, las=2,
           main=paste("Stress*Sex interaction\np-value", format(dFrame[nameIter,]$"adjusted_group:sex",digits=3)),
           xaxt='n')
  axis(1, at=c(1, 2, 3, 4), labels = c("Control\nFemale", "Stress\nFemale", "Control\nMale", "Stress\nMale"), cex.axis=0.9)
  stripchart(corticosterone ~ group*sex,
             data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE )
  
  boxplot( corticosterone ~ cage, las=2,
           main=paste("cage effect\np-value", format(dFrame[nameIter,]$"adjusted_cage",digits=3)), xaxt='n')
  axis(1, at=c(1, 3, 5, 7),
       labels=c("CF1", "SF1", "CM1", "SM1"), cex.axis=0.8)
  axis(1, at=c(2, 4, 6, 8),
       labels=c("CF2", "SF2", "CM2", "SM2"), cex.axis=0.8)
  stripchart(corticosterone ~ cage,
             data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE)
  # This might have to change as well.
  mtext(paste0(dFrame[nameIter,]$taxonName), outer=TRUE, cex = 0.7)
  
  ##         ##             boxplot( bug ~ behave*sex, ylab="Log normalized abundance",
  ##         ##         main=paste("sex*behave interaction p-value", format(dFrame$adjustedpValuesSexBehaveFromMixedInteraction[index],digits=3)),
  ##         ##         xaxt='n')
  ##         ## axis(1, at=c(1, 2, 3, 4), labels = c("Control\nFemale", "Exp.\nFemale", "Control\nMale", "Exp.\nMale"), cex.axis=0.9)
  
  ##         ## stripchart(bug ~ behave*sex,
  ##         ##            data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index] )
  
  ##             mtext(paste0(otuMapping[nameIter,2],"\n",nameIter), outer=TRUE, cex = 0.7)
  
  ##             index = index + 1
  ##         }
}
## This is the selection from the list for however many
## This should be made nicer
# if (otuTable == TRUE) {
#   iterOver <- c(3:8)
# } else {
#   iterOver <- c(2:7)
# }
# for(histIter in iterOver) {
#   hist(dFrame[,histIter], breaks=20, main=colnames(dFrame)[histIter])
# }

dev.off()

