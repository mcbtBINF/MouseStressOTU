## Only works for QIIME right now
## For processing taxonomic levels or OTUtables
## denovo via abundantOTU on R1
## closed via qiime on R1
## abundance models of taxon ~ stress*sex + (1|cage)
## Shannon diversity
## MDS on abundance
## MDS on behavior (4)
## L1 is the simple name for the OTUtable results

## TODO Problems with grabbing the right p-values after the switchover

rm(list=ls())

library("Kendall")
library("pscl")
library("lmtest")
library("nlme")
library("psych")
library("ggplot2")
library("vegan")
library("coin")

otuTable <- TRUE
dataType <- "closedQIIMER1"
## dataType <- "denovoabundantOTU"
#tissues <- c("feces") 
tissues <- c("Cecal Content")

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

## tissues <- unique(myT$Source)
## We are only concerned with fecal (or cecal ) samples.
## divider explores the rarity threshold
## divider <- 8
divider <- 4
if (otuTable == TRUE) {
  taxonomicLevels <- 1
} else {
  taxonomicLevels <- c(2:7)
}
##tissues <- c("feces")

ShannonDiversityList <- list()
for( t in tissues ) {
  ## pValues for Shannon diversity
  pValuesSexS <- vector()
  pValuesGroupS <- vector()
  taxonomicLevel <- vector()
  pValsS <- list()
  indexS <- 1
  ## pValues for MDS Axes Models
  pValuesSexM <- vector()
  pValuesGroupM <- vector()
  pValsM <- list()
  for( taxa in taxonomicLevels) {
    bestAt <- data.frame()
    setwd(processedDir)
    ## inFileName <- "qiimeOTULogNormwithMetadata.txt"
    ## inFileName <- "abundantOTUR1LogNormwithMetadata.txt"
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
    ## Some concerns about column classes
    setwd(analysisDir)
    if(otuTable == TRUE) {
      endMetadataIndex <- which(colnames(myT) == "counts")
    } else {
      endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
    }
    
    subT <- myT[ myT$Source == t, ]
    ## Raw is 22-60
    ## Columns 59 and 60 cause problems with the models
    ## for taxa <-2, c(22:58)
    ## for taxa <-3, c(22:58)
    ## for taxa <-4, c(22:58)
    index <- 1
    names <- vector()
    pValuesGroup <- vector()
    pValuesSex <- vector()
    pValuesBehavior <- vector()
    
    pValuesCage <- vector()
    behavz <- vector()
    pVals <- list()
    
    for ( j in (endMetadataIndex + 1):ncol(subT)) {
      ## DEBUGGING: 66 is the first true case here, then 68
      bug <- subT[,j]
      ## Skip taxa causing numerical instabilities
      ## Are these still true?
      if( sum(bug != 0 ) > nrow(subT) / divider ) {
        if (dataType == "closedQIIMER1") {
          if( (taxa == 2) | ## & j != 84) |
              (taxa == 3) | ## & j != 88) |
              (taxa == 4) | ## & j != 90) |
              (taxa == 5) | ## & j != 129) |
              (taxa == 6) | ## & j != 161) |
              (taxa == 7) ## & j != 144 & j != 183)
          ){
            ## for ( behavior in c(22:58)){
            
            ## taxa 5: j != 255
            ## taxa 6: j != 410
            ## taxa 7: j != 490
            ##  | (taxa == 4 & j != 90)
            ## For a previous model where the behavior was taken into account
            ## pValuesBehave[index] <- wilcox.test( subT[subT$Group=="Experimental", j],
            ##                                   subT[subT$Group=="Control", j])$p.value
            names[index] <- names(subT)[j]
            ## behavz[index] <- colnames(subT)[behavior]
            ## behave <- as.numeric(subT[,behavior])
            sex <- factor(subT$Sex)
            
            cage <- factor(subT$Cage, c("female #5,6,7,8 in same cage", "female #13,14,15,16 in same cage", "female #1,2,3,4 in same cage", "female #9,10,11,12 in same cage",
                                        "male #5,6,7,8 in same cage", "male #13,14,15,16 in same cage", "male #1,2,3,4 in same cage", "male #9,10,11,12 in same cage"))
            
            group <- factor(subT$Group)
            ## Depth was found to not be significant
            ## This may need to be taken out for the otuTable option
            depth <- as.numeric(subT$sampleDepth)
            
            myFrame <- data.frame(bug, sex, group, cage, depth)
            
            ## pValuesSex[index] <- wilcox.test( subT[subT$Sex=="male", j],
            ## subT[subT$Sex=="female", j])$p.value
            pValuesSex[index] <- pvalue(wilcox_test(bug~sex, data=myFrame))
            pValuesGroup[index] <- pvalue(wilcox_test(bug~group, data=myFrame))
            ## pValuesGroup[index] <- wilcox.test( subT[subT$Group=="Control", j],
            ##  subT[subT$Group=="Experimental", j])$p.value
            
            ## Currently not used in this modeling
            ## pValuesBehavior[index] <- cor.test( bug, behave, method="kendall")$p.value[1]
            
            fullModel <-
              gls( bug~  group * sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),data = myFrame )
            fullModelNoc <- gls(bug~  group * sex, data = myFrame)
            
            ##fullModelLME <- lme(bug~ group*sex, method="REML", random = ~1|cage, data = myFrame)
            ## fullModelwithDepth <-
            ##  gls( bug~  group * sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),data = myFrame )
            
            glsModelNoCage <- gls(bug~  group * sex, data = myFrame)
            #print( rownames(anova(glsModelNoCage))[-1])
            ## Merges the nonparametric with the terms of the MLM, along with a cage effect p-value calculation
            pVals[[index]] <- c(pValuesSex[index], pValuesGroup[index], ## pValuesBehavior[index],
                                anova(fullModel)$"p-value"[-1],
                                ##summary(fullModelLME)$tTable[-1,5],
                                anova(fullModel, glsModelNoCage)$"p-value"[2])
            ## pVals[[index]] <- c(pValuesSex[index], pValuesGroup[index], ## pValuesBehavior[index],
            ##                     anova(fullModelwithDepth)$"p-value"[-1],
            ##                     anova(fullModelwithDepth, fullModelNocwithDepth)$"p-value"[2])
            
            ## print(c(taxa, index))
            index <- index + 1
          }
          if (taxa == 1) {
            names[index] <- names(subT)[j]
            ## behavz[index] <- colnames(subT)[behavior]
            ## behave <- as.numeric(subT[,behavior])
            sex <- factor(subT$Sex)
            
            cage <-  factor(subT$Cage, c("female #5,6,7,8 in same cage", "female #13,14,15,16 in same cage", "female #1,2,3,4 in same cage", "female #9,10,11,12 in same cage",
                                         "male #5,6,7,8 in same cage", "male #13,14,15,16 in same cage", "male #1,2,3,4 in same cage", "male #9,10,11,12 in same cage"))
            
            group <- factor(subT$Group)
            ## Depth was found to not be significant
            ## This may need to be taken out for the otuTable option
            depth <- as.numeric(subT$counts)
            
            myFrame <- data.frame(bug, sex, group, cage, depth)
            
            ##pValuesSex[index] <- wilcox.test( subT[subT$Sex=="male", j],
            ##                                  subT[subT$Sex=="female", j])$p.value
            ##pValuesGroup[index] <- wilcox.test( subT[subT$Group=="Control", j],
            ##                                    subT[subT$Group=="Experimental", j])$p.value
            pValuesSex[index] <- pvalue(wilcox_test(bug~sex, data=myFrame))
            pValuesGroup[index] <- pvalue(wilcox_test(bug~group, data=myFrame))
            ## Currently not used in this modeling
            ## pValuesBehavior[index] <- cor.test( bug, behave, method="kendall")$p.value[1]
            
            fullModel <-
              gls( bug~  group * sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),data = myFrame )
            ## fullModelNoc <- gls(bug~  behave * group * sex, data = myFrame)
            
            fullModelLME <- lme(bug~ group*sex, method="REML", random = ~1|cage, data = myFrame)
            
            ## fullModelwithDepth <-
            ##  gls( bug~  group * sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),data = myFrame )
            
            glsModelNoCage <- gls(bug~  group * sex, data = myFrame)
            
            ## Merges the nonparametric with the terms of the MLM, along with a cage effect p-value calculation
            pVals[[index]] <- c(pValuesSex[index], pValuesGroup[index], ## pValuesBehavior[index],
                                anova(fullModel)$"p-value"[-1],
                                ##summary(fullModelLME)$tTable[-1,5],
                                anova(fullModel, glsModelNoCage)$"p-value"[2])
            
            ## pVals[[index]] <- c(pValuesSex[index], pValuesGroup[index], ## pValuesBehavior[index],
            ##                     anova(fullModelwithDepth)$"p-value"[-1],
            ##                     anova(fullModelwithDepth, fullModelNocwithDepth)$"p-value"[2])
            
            #print(c(taxa, index))
            index <- index + 1
            
          }
        }
      }
    }
    
    ## This will not work for the otuTable right now
    if (otuTable == TRUE) {
      namesMapping <- otuMapping[names,2]
      dFrame <- data.frame(names, namesMapping, do.call(rbind, pVals))
      
      colnames(dFrame) <- c("taxonName", "namesMapping",
                            "wilcoxonSex", "wilcoxonStress", ## "kendallBehavior",
                            rownames(anova(glsModelNoCage))[-1], "cage")
      ## print(colnames(dFrame))
      adjOffset <- 3
      rownames(dFrame) <- dFrame$taxonName
    } else {
      dFrame <- data.frame(names, do.call(rbind, pVals))
      colnames(dFrame) <- c("taxonName",
                            "wilcoxonSex", "wilcoxonStress", ## "kendallBehavior",
                            rownames(anova(glsModelNoCage))[-1], "cage")
      adjOffset <- 2
    }
    
    numCols <- ncol(dFrame)
    # Creates MHC columns and names them.
    for (i in adjOffset:numCols){
      dFrame[,i + numCols - (adjOffset - 1)] <- p.adjust(dFrame[,i], method = "BH")
      colnames(dFrame)[i + numCols - (adjOffset - 1)] <- paste0("adjusted_",colnames(dFrame)[i])
    }
    
    bestAt <- rbind(bestAt, dFrame[which(dFrame[grep("adjust", colnames(dFrame))] < 0.05, arr.ind=TRUE)[,1],])

    write.table(dFrame, file=paste(dataType,"_analysis_StressSexCage_gls_L", taxa, "_tissue_",t, "_divider_", divider, ".txt",sep=""), sep="\t",row.names=FALSE)
    
    ## This might be overcounting things.
    bestAt <- unique(bestAt)
    write.table(bestAt, file=paste(dataType,"_SigHitsat_StressSexCage_gls_L", taxa, "_tissue_", t, "_divider_", divider, ".txt",sep=""), sep="\t",row.names=FALSE)
    
    
    pdf( paste(dataType, "_StressSexCage_gls_L", taxa, "_tissue_", t, "_divider_", divider, ".pdf", sep=""))
    
    ## This might break for the otuTable selection
    if(otuTable == TRUE) {
      toIter <- rownames(dFrame)
    } else {
      toIter <- 1:dim(dFrame)[1]
    }
    for (nameIter in toIter) {
      par(mfrow=c(2,1), mgp=c(3, 1.25, 0), oma = c(0, 0, 2, 0))
      bug <- subT[,as.character(dFrame[nameIter,]$taxonName)]
      
      sex <- factor(subT$Sex)
      if(otuTable == TRUE) {
        depth <- as.numeric(subT$counts)
      } else {
        depth <- as.numeric(subT$sampleDepth)
      }
      group <- factor(subT$Group)
      
      cage <-  factor( paste( subT$Housing, subT$Sex, sep=""), c("#5,6,7,8 in same cagefemale", "#13,14,15,16 in same cagefemale", "#1,2,3,4 in same cagefemale", "#9,10,11,12 in same cagefemale",
                                                                 "#5,6,7,8 in same cagemale", "#13,14,15,16 in same cagemale", "#1,2,3,4 in same cagemale", "#9,10,11,12 in same cagemale"))
      
      myFrame <- data.frame(bug, sex, group, depth, cage)
      
      boxplot(outline=FALSE, bug ~ group, data = myFrame,
               main = paste("Wilcoxon Stress\np-value", format(dFrame[nameIter,]$"adjusted_wilcoxonStress",digits=3) ), ylab="Log normalized abundance" )
      stripchart(bug ~ group,
                 data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
      
      boxplot(outline=FALSE, bug ~ sex,
               main = paste("Wilcoxon Sex\np-value", format(dFrame[nameIter,]$"adjusted_wilcoxonSex",digits=3) ), ylab="Log normalized abundance" )
      stripchart(bug ~ sex,
                 data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
      
      mtext(paste0(dFrame[nameIter,]$taxonName,"\n",dFrame[nameIter,]$namesMapping), outer=TRUE, cex = 0.7)
      
      boxplot(outline=FALSE, bug ~ group,
               main = paste("Stress\np-value", format(dFrame[nameIter,]$"adjusted_group",digits=3) ), ylab="Log normalized abundance" )
      stripchart(bug ~ group,
                 data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
      
      boxplot(outline=FALSE, bug ~ sex,
               main = paste("Sex\np-value", format(dFrame[nameIter,]$"adjusted_sex",digits=3) ), ylab="Log normalized abundance" )
      stripchart(bug ~ sex,
                 data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
      
      mtext(paste0(dFrame[nameIter,]$taxonName,"\n",dFrame[nameIter,]$namesMapping), outer=TRUE, cex = 0.7)
      
      boxplot(outline=FALSE, bug ~ group*sex, ylab="Log normalized abundance",
               main=paste("Stress*Sex interaction\np-value", format(dFrame[nameIter,]$"adjusted_group:sex",digits=3)),
               xaxt='n')
      axis(1, at=c(1, 2, 3, 4), labels = c("Control\nFemale", "Stress\nFemale", "Control\nMale", "Stress\nMale"), cex.axis=0.9)
      stripchart(bug ~ group*sex,
                 data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE )
      
      boxplot(outline=FALSE, bug ~ cage ,las=2, ylab="Log normalized abundance",
               main=paste("Cage effect\np-value", format(dFrame[nameIter,]$"adjusted_cage",digits=3)), xaxt='n')
      axis(1, at=c(1, 3, 5, 7),
           labels=c("CF1", "SF1", "CM1", "SM1"), cex.axis=0.8)
      axis(1, at=c(2, 4, 6, 8),
           labels=c("CF2", "SF2", "CM2", "SM2"), cex.axis=0.8)
      stripchart(bug ~ cage,
                 data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
      # This might have to change as well.
      mtext(paste0(dFrame[nameIter,]$taxonName,"\n",dFrame[nameIter,]$namesMapping), outer=TRUE, cex = 0.7)
      
      ##         ##             boxplot(outline=FALSE, bug ~ behave*sex, ylab="Log normalized abundance",
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
    if (otuTable == TRUE) {
      iterOver <- c(3:8)
    } else {
      iterOver <- c(2:7)
    }
    for(histIter in iterOver) {
      hist(dFrame[,histIter], breaks=20, main=colnames(dFrame)[histIter], xlab="Unadjusted p-Values")
    }
    
    dev.off()
    
    ## Just print the best
    pdf( paste(dataType, "_SigHits_StressSexCage_gls_L", taxa, "_tissue_", t, "_divider_", divider, ".pdf", sep=""))
    if(otuTable == TRUE) {
      toIter <- rownames(bestAt)
    } else {
      toIter <- 1:dim(bestAt)[1]
    }
    for (nameIter in toIter) {
      par(mfrow=c(2,1), mgp=c(3, 1.25, 0), oma = c(0, 0, 2, 0))
      bug <- subT[,as.character(bestAt[nameIter,]$taxonName)]
      
      sex <- factor(subT$Sex)
      if (otuTable == TRUE) {
        depth <- as.numeric(subT$counts)
      } else {
        depth <- as.numeric(subT$sampleDepth)
      }
      group <- factor(subT$Group)
      
      cage <-  factor( paste( subT$Housing, subT$Sex, sep=""), c("#5,6,7,8 in same cagefemale", "#13,14,15,16 in same cagefemale", "#1,2,3,4 in same cagefemale", "#9,10,11,12 in same cagefemale",
                                                                 "#5,6,7,8 in same cagemale", "#13,14,15,16 in same cagemale", "#1,2,3,4 in same cagemale", "#9,10,11,12 in same cagemale"))
      
      myFrame <- data.frame(bug, sex, group, depth, cage)
      
      boxplot(outline=FALSE, bug ~ group, data = myFrame,
               main = paste("Wilcoxon Stress\np-value", format(bestAt[nameIter,]$"adjusted_wilcoxonStress",digits=3) ), ylab="Log normalized abundance" )
      stripchart(bug ~ group,
                 data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
      
      boxplot(outline=FALSE, bug ~ sex,
               main = paste("Wilcoxon Sex\np-value", format(bestAt[nameIter,]$"adjusted_wilcoxonSex",digits=3) ), ylab="Log normalized abundance" )
      stripchart(bug ~ sex,
                 data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
      
      mtext(paste0(bestAt[nameIter,]$taxonName,"\n",bestAt[nameIter,]$namesMapping), outer=TRUE, cex = 0.7)
      
      boxplot(outline=FALSE, bug ~ group,
               main = paste("Stress\np-value", format(bestAt[nameIter,]$"adjusted_group",digits=3) ), ylab="Log normalized abundance" )
      stripchart(bug ~ group,
                 data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
      
      boxplot(outline=FALSE, bug ~ sex,
               main = paste("Sex\np-value", format(bestAt[nameIter,]$"adjusted_sex",digits=3) ), ylab="Log normalized abundance" )
      stripchart(bug ~ sex,
                 data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
      
      mtext(paste0(bestAt[nameIter,]$taxonName,"\n",bestAt[nameIter,]$namesMapping), outer=TRUE, cex = 0.7)
      
      boxplot(outline=FALSE, bug ~ group*sex, ylab="Log normalized abundance",
               main=paste("Stress*Sex interaction\np-value", format(bestAt[nameIter,]$"adjusted_group:sex",digits=3)),
               xaxt='n')
      axis(1, at=c(1, 2, 3, 4), labels = c("Control\nFemale", "Stress\nFemale", "Control\nMale", "Stress\nMale"), cex.axis=0.9)
      stripchart(bug ~ group*sex,
                 data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE )
      
      boxplot(outline=FALSE, bug ~ cage ,las=2, ylab="Log normalized abundance",
               main=paste("Cage effect\np-value", format(bestAt[nameIter,]$"adjusted_cage",digits=3)), xaxt='n')
      axis(1, at=c(1, 3, 5, 7),
           labels=c("CF1", "SF1", "CM1", "SM1"), cex.axis=0.8)
      axis(1, at=c(2, 4, 6, 8),
           labels=c("CF2", "SF2", "CM2", "SM2"), cex.axis=0.8)
      stripchart(bug ~ cage,
                 data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
      
      mtext(paste0(bestAt[nameIter,]$taxonName,"\n",bestAt[nameIter,]$namesMapping), outer=TRUE, cex = 0.7)
      
      ##         ##             boxplot(outline=FALSE, bug ~ behave*sex, ylab="Log normalized abundance",
      ##         ##         main=paste("sex*behave interaction p-value", format(dFrame$adjustedpValuesSexBehaveFromMixedInteraction[index],digits=3)),
      ##         ##         xaxt='n')
      ##         ## axis(1, at=c(1, 2, 3, 4), labels = c("Control\nFemale", "Exp.\nFemale", "Control\nMale", "Exp.\nMale"), cex.axis=0.9)
      
      ##         ## stripchart(bug ~ behave*sex,
      ##         ##            data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index] )
      
      ##             mtext(paste0(otuMapping[nameIter,2],"\n",nameIter), outer=TRUE, cex = 0.7)
      
      ##             index = index + 1
      ##         }
    }
    dev.off()
    
    ## Start of Shannon diversity Section
    pdf( paste(dataType, "_StressSexCage_ShannonDiversity_gls_L", taxa, "_tissue_", t, ".pdf", sep=""))
    par(mfrow=c(2,1), mgp=c(3, 1.25, 0), oma = c(0, 0, 2, 0))
    ##     ## Remember that this can cause problems
    bug <- apply(subT[,(endMetadataIndex+1):ncol(subT)],1,diversity)
    ## bug <- subT$ShannonDiversity
    ## ShannonDiversityList[[index]] <- bug
    sex <- factor(subT$Sex)
    
    cage <-  factor(subT$Cage, c("female #5,6,7,8 in same cage", "female #13,14,15,16 in same cage", "female #1,2,3,4 in same cage", "female #9,10,11,12 in same cage",
                                 "male #5,6,7,8 in same cage", "male #13,14,15,16 in same cage", "male #1,2,3,4 in same cage", "male #9,10,11,12 in same cage"))
    
    group <- factor(subT$Group)
    
    myFrame <- data.frame(bug, sex, group, cage)
    
    ##pValuesSexS[indexS] <- wilcox.test( bug[subT$Sex=="male"],
    ##                                    bug[subT$Sex=="female"])$p.value
    ##pValuesGroupS[indexS] <- wilcox.test( bug[subT$Group=="Control"],
    ##                                      bug[subT$Group=="Experimental"])$p.value
    pValuesSexS[indexS] <- pvalue(wilcox_test(bug~sex, data=myFrame))
    pValuesGroupS[indexS] <- pvalue(wilcox_test(bug~group, data=myFrame))
    fullModel <- gls( bug~  group * sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),data = myFrame )
    ## fullModel <- lme(bug~ group*sex, method="REML", random = ~1|cage, data = myFrame)
    
    fullModelNoc <- gls(bug~ group * sex, data = myFrame)
    pValsS[[indexS]] <- c(pValuesSexS[indexS], pValuesGroupS[indexS],
                          ##summary(fullModel)$tTable[-1,5],
                          anova(fullModel)$"p-value"[-1],
                          anova(fullModel, fullModelNoc)$"p-value"[2])
    names(pValsS[[indexS]]) <- c("wilcoxonSex", "wilcoxonStress", "Stress", "Sex", "Stress:Sex", "Cage")
    
    boxplot(outline=FALSE, bug ~ group, data = myFrame,
             main = paste("Wilcoxon Stress\np-value", format(pValsS[[indexS]]["wilcoxonStress"],digits=3) ), ylab="Shannon diversity" )
    stripchart(bug ~ group,
               data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
    
    boxplot(outline=FALSE, bug ~ sex,
             main = paste("Wilcoxon Sex\np-value", format(pValsS[[indexS]]["wilcoxonSex"],digits=3) ), ylab="Shannon diversity" )
    stripchart(bug ~ sex,
               data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
    
    mtext(paste0("Shannon diversity at level ", taxa), outer=TRUE, cex = 0.7)
    
    boxplot(outline=FALSE, bug ~ group,
             main = paste("Stress\np-value", format(pValsS[[indexS]]["Stress"],digits=3) ), ylab="Shannon diversity" )
    stripchart(bug ~ group,
               data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
    
    boxplot(outline=FALSE, bug ~ sex,
             main = paste("Sex\np-value", format(pValsS[[indexS]]["Sex"],digits=3) ), ylab="Shannon diversity" )
    stripchart(bug ~ sex,
               data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
    
    mtext(paste0("Shannon diversity at level ", taxa), outer=TRUE, cex = 0.7)
    
    boxplot(outline=FALSE, bug ~ group*sex, ylab="Shannon diversity",
             main=paste("Stress*Sex interaction\np-value", format(pValsS[[indexS]]["Stress:Sex"],digits=3)),
             xaxt='n')
    axis(1, at=c(1, 2, 3, 4), labels = c("Control\nFemale", "Stress\nFemale", "Control\nMale", "Stress\nMale"), cex.axis=0.9)
    stripchart(bug ~ group*sex,
               data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE )
    
    ## mtext(paste0(dFrame[nameIter,]$taxonName), outer=TRUE, cex = 0.7)
    
    boxplot(outline=FALSE, bug ~ cage ,las=2, ylab="Shannon diversity",
             main=paste("Cage effect\np-value", format(pValsS[[indexS]]["Cage"],digits=3)), xaxt='n')
    axis(1, at=c(1, 3, 5, 7),
         labels=c("CF1", "SF1", "CM1", "SM1"), cex.axis=0.8)
    axis(1, at=c(2, 4, 6, 8),
         labels=c("CF2", "SF2", "CM2", "SM2"), cex.axis=0.8)
    stripchart(bug ~ cage,
               data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
    
    mtext(paste0("Shannon diversity at level ", taxa), outer=TRUE, cex = 0.7)
    
    taxonomicLevel[[indexS]] <- taxa
    ## dSFrame <- data.frame(taxa, do.call(rbind, pValsS))
    ## colnames(dSFrame) <- c("taxonomicLevel", names(pValsS[[indexS]]))
    indexS <- indexS + 1
    
    dev.off()
    
    ## Start of MDS Axes section
    ## Abundance
    
    ## Will need to create this for the otuTable data :important:
    if (t == "Cecal Content") {
      mdsData <- read.table(file=paste(dataType, "_L_", taxa, "_cecal_pcoa.txt",sep=""), header=TRUE)
      mdsEigen <- read.table(file=paste(dataType, "_L_", taxa, "_cecal_eigenValues.txt", sep=""), header=TRUE)
    } else {
      mdsData <- read.table(file=paste(dataType, "_L", taxa, "_fecal_pcoa.txt",sep=""), header=TRUE)
      mdsEigen <- read.table(file=paste(dataType, "_L", taxa, "_fecal_eigenValues.txt", sep=""), header=TRUE)
    }
    percentVariance <- mdsEigen/sum(mdsEigen)
    sex <- factor(subT$Sex)
    
    cage <-  factor(subT$Cage, c("female #5,6,7,8 in same cage", "female #13,14,15,16 in same cage", "female #1,2,3,4 in same cage", "female #9,10,11,12 in same cage",
                                 "male #5,6,7,8 in same cage", "male #13,14,15,16 in same cage", "male #1,2,3,4 in same cage", "male #9,10,11,12 in same cage"))
    
    group <- factor(subT$Group)
    
    for (mdsI in 1:10){
      bug <- mdsData[,mdsI]
      myFrame <- data.frame(bug, sex, group, cage)
      ##pValuesSexM[mdsI] <- wilcox.test( bug[subT$Sex=="male"],
      ##                                  bug[subT$Sex=="female"])$p.value
      ##pValuesGroupM[mdsI] <- wilcox.test( bug[subT$Group=="Control"],
      ##                                    bug[subT$Group=="Experimental"])$p.value
      pValuesSexM[mdsI] <- pvalue(wilcox_test(bug~sex, data=myFrame))
      pValuesGroupM[mdsI] <- pvalue(wilcox_test(bug~group, data=myFrame))
      fullModel <-
        gls( bug~  group * sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),data = myFrame )
      ## fullModel <- lme(bug~ group*sex, method="REML", random = ~1|cage, data = myFrame)
      
      fullModelNoc <- gls(bug~ group * sex, data = myFrame)
      pValsM[[mdsI]] <- c(pValuesSexM[mdsI], pValuesGroupM[mdsI],
                          anova(fullModel)$"p-value"[-1],
                          ##summary(fullModel)$tTable[-1,5],
                          anova(fullModel, fullModelNoc)$"p-value"[2])
      
    }
    dMFrame <- data.frame(1:10, percentVariance[1:10,], do.call(rbind, pValsM))
    colnames(dMFrame) <- c("MDSAxis", "percentVariance",
                           "wilcoxonSex", "wilcoxonStress", ## "kendallBehavior",
                           rownames(anova(fullModelNoc))[-1], "cage")
    
    numCols <- ncol(dMFrame)
    for (i in 3:numCols){
      dMFrame[,i + numCols - 2] <- p.adjust(dMFrame[,i], method = "BH")
      colnames(dMFrame)[i + numCols - 2] <- paste0("adjusted_",colnames(dMFrame)[i])
    }
    
    pdf( paste(dataType, "_AbundanceMDS_StressSexCage_gls_L", taxa, "_tissue_", t, ".pdf", sep=""))
    par(mfrow=c(2,1), mgp=c(3, 1.25, 0), oma = c(0, 0, 2, 0))
    for (mdsI in 1:10){
      bug <- mdsData[,mdsI]
      myFrame <- data.frame(bug, sex, group, cage)
      boxplot(outline=FALSE, bug ~ group, data = myFrame,
               main = paste("Wilcoxon Stress\np-value", format(dMFrame[mdsI,]$adjusted_wilcoxonStress,digits=3) ), ylab=paste0("MDS Axis ", mdsI) )
      stripchart(bug ~ group,
                 data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
      
      boxplot(outline=FALSE, bug ~ sex,
               main = paste("Wilcoxon Sex\np-value", format(dMFrame[mdsI,]$adjusted_wilcoxonSex,digits=3) ), ylab=paste0("MDS Axis ", mdsI))
      stripchart(bug ~ sex,
                 data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
      
      mtext(paste0("MDS Axes at taxonomic level ", taxa), outer=TRUE, cex = 0.7)
      
      boxplot(outline=FALSE, bug ~ group,
               main = paste("Stress\np-value", format(dMFrame[mdsI,]$adjusted_group,digits=3) ), ylab=paste0("MDS Axis ", mdsI) )
      stripchart(bug ~ group,
                 data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
      
      boxplot(outline=FALSE, bug ~ sex,
               main = paste("Sex\np-value", format(dMFrame[mdsI,]$adjusted_sex,digits=3) ), ylab=paste0("MDS Axis ", mdsI) )
      stripchart(bug ~ sex,
                 data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
      
      mtext(paste0("MDS Axes at taxonomic level ", taxa), outer=TRUE, cex = 0.7)
      
      boxplot(outline=FALSE, bug ~ group*sex, ylab=paste0("MDS Axis ", mdsI),
               main=paste("Stress*Sex interaction p-value", format(dMFrame[mdsI,]$"adjusted_group:sex",digits=3)),
               xaxt='n')
      axis(1, at=c(1, 2, 3, 4), labels = c("Control\nFemale", "Stress\nFemale", "Control\nMale", "Stress\nMale"), cex.axis=0.9)
      stripchart(bug ~ group*sex,
                 data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE )
      
      ## mtext(paste0(dFrame[nameIter,]$taxonName), outer=TRUE, cex = 0.7)
      
      boxplot(outline=FALSE, bug ~ cage ,las=2, ylab=paste0("MDS Axis ", mdsI),
               main=paste("cage effect p-value", format(dMFrame[mdsI,]$adjusted_cage,digits=3)), xaxt='n')
      axis(1, at=c(1, 3, 5, 7),
           labels=c("CF1", "SF1", "CM1", "SM1"), cex.axis=0.8)
      axis(1, at=c(2, 4, 6, 8),
           labels=c("CF2", "SF2", "CM2", "SM2"), cex.axis=0.8)
      stripchart(bug ~ cage,
                 data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
      
      mtext(paste0("MDS Axes at taxonomic level ", taxa), outer=TRUE, cex = 0.7)
      
    }
    dev.off()
    write.table(dMFrame, file=paste(dataType,"_MDSAxes_StressSexCage_gls_L", taxa, "_tissue_", t, ".txt",sep=""), sep="\t",row.names=FALSE)
    
    ## Not on the abundance data, only on the behavioral data
    ## Needs to be done exactly once then
    if (taxa == 2 & t == "feces") {
      mdsData <- read.table(file=paste(dataType, "_L2", "_behavior_pcoa.txt",sep=""), header=TRUE)
      mdsEigen <- read.table(file=paste(dataType, "_L2", "_behavior_eigenValues.txt", sep=""), header=TRUE)
      percentVariance <- mdsEigen/sum(mdsEigen)
      
      for (mdsI in 1:10){
        bug <- mdsData[,mdsI]
        myFrame <- data.frame(bug, sex, group, cage)
        ##pValuesSexM[mdsI] <- wilcox.test( bug[subT$Sex=="male"],
        ##                                  bug[subT$Sex=="female"])$p.value
        ##pValuesGroupM[mdsI] <- wilcox.test( bug[subT$Group=="Control"],
        ##                                    bug[subT$Group=="Experimental"])$p.value
        pValuesSexM[mdsI] <- pvalue(wilcox_test(bug~sex, data=myFrame))
        pValuesGroupM[mdsI] <- pvalue(wilcox_test(bug~group, data=myFrame))
        fullModel <-
          gls( bug~  group * sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),data = myFrame )
        ## fullModel <- lme(bug~ group*sex, method="REML", random = ~1|cage, data = myFrame)
        
        fullModelNoc <- gls(bug~ group * sex, data = myFrame)
        pValsM[[mdsI]] <- c(pValuesSexM[mdsI], pValuesGroupM[mdsI],
                            ## summary(fullModel)$tTable[-1,5],
                            anova(fullModel)$"p-value"[-1],
                            anova(fullModel, fullModelNoc)$"p-value"[2])
        
      }
      dMFrame <- data.frame(1:10, percentVariance[1:10,], do.call(rbind, pValsM))
      colnames(dMFrame) <- c("MDSAxis", "percentVariance",
                             "wilcoxonSex", "wilcoxonStress", ## "kendallBehavior",
                             rownames(anova(fullModelNoc))[-1], "cage")
      
      numCols <- ncol(dMFrame)
      for (i in 3:numCols){
        dMFrame[,i + numCols - 2] <- p.adjust(dMFrame[,i], method = "BH")
        colnames(dMFrame)[i + numCols - 2] <- paste0("adjusted_",colnames(dMFrame)[i])
      }
      
      pdf( paste(dataType, "_AllBehaviorMDS_StressSexCage_gls_",t, ".pdf", sep=""))
      par(mfrow=c(2,1), mgp=c(3, 1.25, 0), oma = c(0, 0, 2, 0))
      for (mdsI in 1:10){
        bug <- mdsData[,mdsI]
        myFrame <- data.frame(bug, sex, group, cage)
        boxplot(outline=FALSE, bug ~ group, data = myFrame,
                 main = paste("Wilcoxon Stress\np-value", format(dMFrame[mdsI,]$adjusted_wilcoxonStress,digits=3) ), ylab=paste0("MDS Axis ", mdsI) )
        stripchart(bug ~ group,
                   data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
        
        boxplot(outline=FALSE, bug ~ sex,
                 main = paste("Wilcoxon Sex\np-value", format(dMFrame[mdsI,]$adjusted_wilcoxonSex,digits=3) ), ylab=paste0("MDS Axis ", mdsI))
        stripchart(bug ~ sex,
                   data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
        
        mtext(paste0("MDS Axes at taxonomic level ", taxa), outer=TRUE, cex = 0.7)
        
        boxplot(outline=FALSE, bug ~ group,
                 main = paste("Stress\np-value", format(dMFrame[mdsI,]$adjusted_group,digits=3) ), ylab=paste0("MDS Axis ", mdsI) )
        stripchart(bug ~ group,
                   data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
        
        boxplot(outline=FALSE, bug ~ sex,
                 main = paste("Sex\np-value", format(dMFrame[mdsI,]$adjusted_sex,digits=3) ), ylab=paste0("MDS Axis ", mdsI) )
        stripchart(bug ~ sex,
                   data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
        
        mtext(paste0("MDS Axes at taxonomic level ", taxa), outer=TRUE, cex = 0.7)
        
        boxplot(outline=FALSE, bug ~ group*sex, ylab=paste0("MDS Axis ", mdsI),
                 main=paste("Stress*Sex interaction\np-value", format(dMFrame[mdsI,]$"adjusted_group:sex",digits=3)),
                 xaxt='n')
        axis(1, at=c(1, 2, 3, 4), labels = c("Control\nFemale", "Stress\nFemale", "Control\nMale", "Stress\nMale"), cex.axis=0.9)
        stripchart(bug ~ group*sex,
                   data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE )
        
        ## mtext(paste0(dFrame[nameIter,]$taxonName), outer=TRUE, cex = 0.7)
        
        boxplot(outline=FALSE, bug ~ cage ,las=2, ylab=paste0("MDS Axis ", mdsI),
                 main=paste("Cage effect\np-value", format(dMFrame[mdsI,]$adjusted_cage,digits=3)), xaxt='n')
        axis(1, at=c(1, 3, 5, 7),
             labels=c("CF1", "SF1", "CM1", "SM1"), cex.axis=0.8)
        axis(1, at=c(2, 4, 6, 8),
             labels=c("CF2", "SF2", "CM2", "SM2"), cex.axis=0.8)
        stripchart(bug ~ cage,
                   data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
        
        mtext(paste0("MDS Axes at taxonomic level ", taxa), outer=TRUE, cex = 0.7)
        
      }
      dev.off()
      
      write.table(dMFrame, file=paste(dataType,"_MDSAxesAllBehavior_StressSexCage_gls_",t, ".txt",sep=""), sep="\t",row.names=FALSE)
      
      ## Start of EPM Data
      mdsData <- read.table(file=paste("LyteSharon_r01_epm_pcoa.txt",sep=""), header=TRUE)
      mdsEigen <- read.table(file=paste("LyteSharon_r01_epm_eigenValues.txt", sep=""), header=TRUE)
      percentVariance <- mdsEigen/sum(mdsEigen)
      pValuesSexM <- vector()
      pValuesGroupM <- vector()
      pValsM <- list()
      
      for (mdsI in 1:5){
        bug <- mdsData[,mdsI]
        myFrame <- data.frame(bug, sex, group, cage)
        ##pValuesSexM[mdsI] <- wilcox.test( bug[subT$Sex=="male"],
        ##                                  bug[subT$Sex=="female"])$p.value
        ##pValuesGroupM[mdsI] <- wilcox.test( bug[subT$Group=="Control"],
        ##                                    bug[subT$Group=="Experimental"])$p.value
        pValuesSexM[mdsI] <- pvalue(wilcox_test(bug~sex, data=myFrame))
        pValuesGroupM[mdsI] <- pvalue(wilcox_test(bug~group, data=myFrame))
        fullModel <-
          gls( bug~  group * sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),data = myFrame )
        ## fullModel <- lme(bug~ group*sex, method="REML", random = ~1|cage, data = myFrame)
        
        fullModelNoc <- gls(bug~ group * sex, data = myFrame)
        pValsM[[mdsI]] <- c(pValuesSexM[mdsI], pValuesGroupM[mdsI],
                            ##summary(fullModel)$tTable[-1,5],
                            anova(fullModel)$"p-value"[-1],
                            anova(fullModel, fullModelNoc)$"p-value"[2])
        
      }
      dMFrame <- data.frame(1:5, percentVariance[1:5,], do.call(rbind, pValsM))
      colnames(dMFrame) <- c("MDSAxis", "percentVariance",
                             "wilcoxonSex", "wilcoxonStress", ## "kendallBehavior",
                             rownames(anova(fullModelNoc))[-1], "cage")
      
      numCols <- ncol(dMFrame)
      for (i in 3:numCols){
        dMFrame[,i + numCols - 2] <- p.adjust(dMFrame[,i], method = "BH")
        colnames(dMFrame)[i + numCols - 2] <- paste0("adjusted_",colnames(dMFrame)[i])
      }
      
      pdf( paste(dataType, "_EPMBehaviorMDS_StressSexCage_gls_",t, ".pdf", sep=""))
      par(mfrow=c(2,1), mgp=c(3, 1.25, 0), oma = c(0, 0, 2, 0))
      for (mdsI in 1:5){
        bug <- mdsData[,mdsI]
        myFrame <- data.frame(bug, sex, group, cage)
        boxplot(outline=FALSE, bug ~ group, data = myFrame,
                 main = paste("Wilcoxon Stress\np-value", format(dMFrame[mdsI,]$adjusted_wilcoxonStress,digits=3) ), ylab=paste0("MDS Axis ", mdsI) )
        stripchart(bug ~ group,
                   data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
        
        boxplot(outline=FALSE, bug ~ sex,
                 main = paste("Wilcoxon Sex\np-value", format(dMFrame[mdsI,]$adjusted_wilcoxonSex,digits=3) ), ylab=paste0("MDS Axis ", mdsI))
        stripchart(bug ~ sex,
                   data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
        
        mtext(paste0("MDS Axes at taxonomic level ", taxa), outer=TRUE, cex = 0.7)
        
        boxplot(outline=FALSE, bug ~ group,
                 main = paste("Stress\np-value", format(dMFrame[mdsI,]$adjusted_group,digits=3) ), ylab=paste0("MDS Axis ", mdsI) )
        stripchart(bug ~ group,
                   data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
        
        boxplot(outline=FALSE, bug ~ sex,
                 main = paste("Sex\np-value", format(dMFrame[mdsI,]$adjusted_sex,digits=3) ), ylab=paste0("MDS Axis ", mdsI) )
        stripchart(bug ~ sex,
                   data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
        
        mtext(paste0("MDS Axes at taxonomic level ", taxa), outer=TRUE, cex = 0.7)
        
        boxplot(outline=FALSE, bug ~ group*sex, ylab=paste0("MDS Axis ", mdsI),
                 main=paste("Stress*Sex interaction\np-value", format(dMFrame[mdsI,]$"adjusted_group:sex",digits=3)),
                 xaxt='n')
        axis(1, at=c(1, 2, 3, 4), labels = c("Control\nFemale", "Stress\nFemale", "Control\nMale", "Stress\nMale"), cex.axis=0.9)
        stripchart(bug ~ group*sex,
                   data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE )
        
        ## mtext(paste0(dFrame[nameIter,]$taxonName), outer=TRUE, cex = 0.7)
        
        boxplot(outline=FALSE, bug ~ cage ,las=2, ylab=paste0("MDS Axis ", mdsI),
                 main=paste("Cage effect\np-value", format(dMFrame[mdsI,]$adjusted_cage,digits=3)), xaxt='n')
        axis(1, at=c(1, 3, 5, 7),
             labels=c("CF1", "SF1", "CM1", "SM1"), cex.axis=0.8)
        axis(1, at=c(2, 4, 6, 8),
             labels=c("CF2", "SF2", "CM2", "SM2"), cex.axis=0.8)
        stripchart(bug ~ cage,
                   data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
        
        mtext(paste0("MDS Axes at taxonomic level ", taxa), outer=TRUE, cex = 0.7)
      }
      dev.off()
      
      write.table(dMFrame, file=paste(dataType,"_MDSAxesEPMBehavior_StressSexCage_gls_",t, ".txt",sep=""), sep="\t",row.names=FALSE)
      
      ## Start of lightdark data
      mdsData <- read.table(file=paste("LyteSharon_r01_lightdark_pcoa.txt",sep=""), header=TRUE)
      mdsEigen <- read.table(file=paste("LyteSharon_r01_lightdark_eigenValues.txt", sep=""), header=TRUE)
      percentVariance <- mdsEigen/sum(mdsEigen)
      pValuesSexM <- vector()
      pValuesGroupM <- vector()
      pValsM <- list()
      
      for (mdsI in 1:5){
        bug <- mdsData[,mdsI]
        myFrame <- data.frame(bug, sex, group, cage)
        ##pValuesSexM[mdsI] <- wilcox.test( bug[subT$Sex=="male"],
        ##                                  bug[subT$Sex=="female"])$p.value
        ##pValuesGroupM[mdsI] <- wilcox.test( bug[subT$Group=="Control"],
        ##                                    bug[subT$Group=="Experimental"])$p.value
        pValuesSexM[mdsI] <- pvalue(wilcox_test(bug~sex, data=myFrame))
        pValuesGroupM[mdsI] <- pvalue(wilcox_test(bug~group, data=myFrame))
        fullModel <-
          gls( bug~  group * sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),data = myFrame )
        ## fullModel <- lme(bug~ group*sex, method="REML", random = ~1|cage, data = myFrame)
        
        fullModelNoc <- gls(bug~ group * sex, data = myFrame)
        pValsM[[mdsI]] <- c(pValuesSexM[mdsI], pValuesGroupM[mdsI],
                            ##summary(fullModel)$tTable[-1,5],
                            anova(fullModel)$"p-value"[-1],
                            anova(fullModel, fullModelNoc)$"p-value"[2])
        
      }
      dMFrame <- data.frame(1:5, percentVariance[1:5,], do.call(rbind, pValsM))
      colnames(dMFrame) <- c("MDSAxis", "percentVariance",
                             "wilcoxonSex", "wilcoxonStress", ## "kendallBehavior",
                             rownames(anova(fullModelNoc))[-1], "cage")
      
      numCols <- ncol(dMFrame)
      for (i in 3:numCols){
        dMFrame[,i + numCols - 2] <- p.adjust(dMFrame[,i], method = "BH")
        colnames(dMFrame)[i + numCols - 2] <- paste0("adjusted_",colnames(dMFrame)[i])
      }
      
      pdf( paste(dataType, "_LightDarkBehaviorMDS_StressSexCage_gls_",t, ".pdf", sep=""))
      par(mfrow=c(2,1), mgp=c(3, 1.25, 0), oma = c(0, 0, 2, 0))
      for (mdsI in 1:5){
        bug <- mdsData[,mdsI]
        myFrame <- data.frame(bug, sex, group, cage)
        boxplot(outline=FALSE, bug ~ group, data = myFrame,
                 main = paste("Wilcoxon Stress\np-value", format(dMFrame[mdsI,]$adjusted_wilcoxonStress,digits=3) ), ylab=paste0("MDS Axis ", mdsI) )
        
        stripchart(bug ~ group,
                   data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
        
        boxplot(outline=FALSE, bug ~ sex,
                 main = paste("Wilcoxon Sex\np-value", format(dMFrame[mdsI,]$adjusted_wilcoxonSex,digits=3) ), ylab=paste0("MDS Axis ", mdsI))
        stripchart(bug ~ sex,
                   data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
        
        mtext(paste0("MDS Axes at taxonomic level ", taxa), outer=TRUE, cex = 0.7)
        
        boxplot(outline=FALSE, bug ~ group,
                 main = paste("Stress\np-value", format(dMFrame[mdsI,]$adjusted_group,digits=3) ), ylab=paste0("MDS Axis ", mdsI) )
        stripchart(bug ~ group,
                   data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
        
        boxplot(outline=FALSE, bug ~ sex,
                 main = paste("Sex\np-value", format(dMFrame[mdsI,]$adjusted_sex,digits=3) ), ylab=paste0("MDS Axis ", mdsI) )
        stripchart(bug ~ sex,
                   data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
        
        mtext(paste0("MDS Axes at taxonomic level ", taxa), outer=TRUE, cex = 0.7)
        
        boxplot(outline=FALSE, bug ~ group*sex, ylab=paste0("MDS Axis ", mdsI),
                 main=paste("Stress*Sex interaction\np-value", format(dMFrame[mdsI,]$"adjusted_group:sex",digits=3)),
                 xaxt='n')
        axis(1, at=c(1, 2, 3, 4), labels = c("Control\nFemale", "Stress\nFemale", "Control\nMale", "Stress\nMale"), cex.axis=0.9)
        stripchart(bug ~ group*sex,
                   data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE )
        
        ## mtext(paste0(dFrame[nameIter,]$taxonName), outer=TRUE, cex = 0.7)
        
        boxplot(outline=FALSE, bug ~ cage ,las=2, ylab=paste0("MDS Axis ", mdsI),
                 main=paste("Cage effect\np-value", format(dMFrame[mdsI,]$adjusted_cage,digits=3)), xaxt='n')
        axis(1, at=c(1, 3, 5, 7),
             labels=c("CF1", "SF1", "CM1", "SM1"), cex.axis=0.8)
        axis(1, at=c(2, 4, 6, 8),
             labels=c("CF2", "SF2", "CM2", "SM2"), cex.axis=0.8)
        stripchart(bug ~ cage,
                   data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
        
        mtext(paste0("MDS Axes at taxonomic level ", taxa), outer=TRUE, cex = 0.7)
      }
      dev.off()
      
      write.table(dMFrame, file=paste(dataType,"_MDSAxesLightDarkBehavior_StressSexCage_gls_",t, ".txt",sep=""), sep="\t",row.names=FALSE)
      
      ## Start of openfield data
      mdsData <- read.table(file=paste("LyteSharon_r01_openfield_pcoa.txt",sep=""), header=TRUE)
      mdsEigen <- read.table(file=paste("LyteSharon_r01_openfield_eigenValues.txt", sep=""), header=TRUE)
      percentVariance <- mdsEigen/sum(mdsEigen)
      pValuesSexM <- vector()
      pValuesGroupM <- vector()
      pValsM <- list()
      
      for (mdsI in 1:5){
        bug <- mdsData[,mdsI]
        myFrame <- data.frame(bug, sex, group, cage)
        ##pValuesSexM[mdsI] <- wilcox.test( bug[subT$Sex=="male"],
        ##                                  bug[subT$Sex=="female"])$p.value
        ##pValuesGroupM[mdsI] <- wilcox.test( bug[subT$Group=="Control"],
        ##                                    bug[subT$Group=="Experimental"])$p.value
        pValuesSexM[mdsI] <- pvalue(wilcox_test(bug~sex, data=myFrame))
        pValuesGroupM[mdsI] <- pvalue(wilcox_test(bug~group, data=myFrame))
        fullModel <-
          gls( bug~  group * sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),data = myFrame )
        ## fullModel <- lme(bug~ group*sex, method="REML", random = ~1|cage, data = myFrame)
        
        fullModelNoc <- gls(bug~ group * sex, data = myFrame)
        pValsM[[mdsI]] <- c(pValuesSexM[mdsI], pValuesGroupM[mdsI],
                            anova(fullModel)$"p-value"[-1],
                            anova(fullModel, fullModelNoc)$"p-value"[2])
        
      }
      dMFrame <- data.frame(1:5, percentVariance[1:5,], do.call(rbind, pValsM))
      colnames(dMFrame) <- c("MDSAxis", "percentVariance",
                             "wilcoxonSex", "wilcoxonStress", ## "kendallBehavior",
                             rownames(anova(fullModelNoc))[-1], "cage")
      
      numCols <- ncol(dMFrame)
      for (i in 3:numCols){
        dMFrame[,i + numCols - 2] <- p.adjust(dMFrame[,i], method = "BH")
        colnames(dMFrame)[i + numCols - 2] <- paste0("adjusted_",colnames(dMFrame)[i])
      }
      
      pdf( paste(dataType, "_OpenFieldBehaviorMDS_StressSexCage_gls_L", taxa, "_", t,".pdf", sep=""))
      par(mfrow=c(2,1), mgp=c(3, 1.25, 0), oma = c(0, 0, 2, 0))
      for (mdsI in 1:5){
        bug <- mdsData[,mdsI]
        myFrame <- data.frame(bug, sex, group, cage)
        boxplot(outline=FALSE, bug ~ group, data = myFrame,
                 main = paste("Wilcoxon Stress\np-value", format(dMFrame[mdsI,]$adjusted_wilcoxonStress,digits=3) ), ylab=paste0("MDS Axis ", mdsI) )
        stripchart(bug ~ group,
                   data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
        
        boxplot(outline=FALSE, bug ~ sex,
                 main = paste("Wilcoxon Sex\np-value", format(dMFrame[mdsI,]$adjusted_wilcoxonSex,digits=3) ), ylab=paste0("MDS Axis ", mdsI))
        stripchart(bug ~ sex,
                   data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
        
        mtext(paste0("MDS Axes at taxonomic level ", taxa), outer=TRUE, cex = 0.7)
        
        boxplot(outline=FALSE, bug ~ group,
                 main = paste("Stress\np-value", format(dMFrame[mdsI,]$adjusted_group,digits=3) ), ylab=paste0("MDS Axis ", mdsI) )
        stripchart(bug ~ group,
                   data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
        
        boxplot(outline=FALSE, bug ~ sex,
                 main = paste("Sex\np-value", format(dMFrame[mdsI,]$adjusted_sex,digits=3) ), ylab=paste0("MDS Axis ", mdsI) )
        stripchart(bug ~ sex,
                   data = myFrame, vertical = TRUE, pch = 21, method="jitter", add=TRUE, ylab = names[index])
        
        mtext(paste0("MDS Axes at taxonomic level ", taxa), outer=TRUE, cex = 0.7)
        
        boxplot(outline=FALSE, bug ~ group*sex, ylab=paste0("MDS Axis ", mdsI),
                 main=paste("Stress*Sex interaction\np-value", format(dMFrame[mdsI,]$"adjusted_group:sex",digits=3)),
                 xaxt='n')
        axis(1, at=c(1, 2, 3, 4), labels = c("Control\nFemale", "Stress\nFemale", "Control\nMale", "Stress\nMale"), cex.axis=0.9)
        stripchart(bug ~ group*sex,
                   data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE )
        
        ## mtext(paste0(dFrame[nameIter,]$taxonName), outer=TRUE, cex = 0.7)
        
        boxplot(outline=FALSE, bug ~ cage ,las=2, ylab=paste0("MDS Axis ", mdsI),
                 main=paste("Cage effect\np-value", format(dMFrame[mdsI,]$adjusted_cage,digits=3)), xaxt='n')
        axis(1, at=c(1, 3, 5, 7),
             labels=c("CF1", "SF1", "CM1", "SM1"), cex.axis=0.8)
        axis(1, at=c(2, 4, 6, 8),
             labels=c("CF2", "SF2", "CM2", "SM2"), cex.axis=0.8)
        stripchart(bug ~ cage,
                   data = myFrame,vertical = TRUE, pch = 21, method="jitter", add=TRUE)
        
        mtext(paste0("MDS Axes at taxonomic level ", taxa), outer=TRUE, cex = 0.7)
      }
      dev.off()
      
      write.table(dMFrame, file=paste(dataType,"_MDSAxesOpenFieldBehavior_StressSexCage_gls_",t, ".txt",sep=""), sep="\t",row.names=FALSE)
    }
  }
  ## Finish up the Shannon diversity calculations
  dSFrame <- data.frame(taxonomicLevel, do.call(rbind, pValsS))
  colnames(dSFrame) <- c("taxonomicLevel", names(pValsS[[1]]))
  write.table(dSFrame, file=paste(dataType,"_ShannonDiversity_StressSexCage_gls_",t, ".txt",sep=""), sep="\t",row.names=FALSE)
}