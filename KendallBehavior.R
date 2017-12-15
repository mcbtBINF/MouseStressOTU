rm(list=ls())

library("Kendall")
library("pscl")
library("lmtest")
library("nlme")
library("psych")
library("ggplot2")
library("coin")

## setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/ArgonneSequencing/JanuaryResequencing/rawData/rg_results/")

## inFileName <- "qiimeOTULogNormwithMetadata.txt"
## myT <- read.table(inFileName,header=TRUE,sep="\t")

## otuMapping <- read.table("LyteSharon_r01_cr_MAPPING.txt", sep="\t", header=TRUE)
## otuMapping[,1] <- paste0("X",otuMapping[,1])
## rownames(otuMapping)<-otuMapping[,1]

## The otuTable route is not currently working
otuTable <- TRUE
behaviorColumns <- c(43:60) #EPM c(22:32) #LD c(33:42) #OF c(43:60)
dataType <- "closedQIIMER1"
tissues <- c("Cecal Content") #c("feces") # c("Cecal Content")

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

##baseDir <- "/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/ArgonneSequencing/JanuaryResequencing/rawData/rg_results/closedQIIMER1"
metadataDir <- paste(baseDir, "metadata", sep="/")
dataDir <- paste(baseDir, "data", sep="/")
processedDir <- paste(baseDir, "processed", sep="/")
analysisDir <- paste(baseDir, "analysis", sep="/")

if(otuTable == TRUE) {
  if ( dataType == "closedQIIMER1") {
    otuMapping <- read.table(paste(processedDir,"LyteSharon_r01_cr_MAPPING.txt",sep="/"), sep="\t", header=TRUE)
  } else {
    otuMapping <- read.table(paste(processedDir,"LyteSharon_R01_PL_wTaxaRDP80_MAPPING.txt", sep="/"), sep="\t", header=TRUE)
  }
  otuMapping[,1] <- paste0("X",otuMapping[,1])
  rownames(otuMapping)<-otuMapping[,1]
}

divider <- 4
if (otuTable == TRUE) {
  taxonomicLevels <- 1
} else {
  taxonomicLevels <- c(2:7)
}

##tissues <- unique(myT$Source)
##taxonomicLevels <- c(2, 3, 4, 5, 6, 7)
## TODO Need to fix this on the other levels

for( t in tissues ) {
  for( taxa in taxonomicLevels) {
    bestAt <- data.frame()
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
    ## Some concerns about column classes
    setwd(analysisDir)
    if(otuTable == TRUE) {
      endMetadataIndex <- which(colnames(myT) == "counts")
    } else {
      endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
    }
    #endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
    colnames(myT)[behaviorColumns] <- gsub("\\.\\.\\.", " # ", colnames(myT)[behaviorColumns])
    colnames(myT)[behaviorColumns] <- gsub("\\.", " ", colnames(myT)[behaviorColumns])
    
    subT <- myT[ myT$Source == t, ]
    # this is the split
    #subT <- subT[subT$Sex == "male",]
    ## Raw is 22-60
    ## for taxa <-2, c(22:58)
    ## for taxa <-3, c(22:58)
    ## for taxa <-4, c(22:58)
    index <- 1
    ## for ( behavior in c(22:58)){
    
    ## pdf( paste(dataType, "_KendallBehavior_", t, colnames(myT)[behavior], ".pdf", sep=""))
    ## TODO Will also need to save out a version over all the behaviors
    ## index <- 1
    names <- vector()
    pValuesGroup <- vector()
    pValuesSex <- vector()
    pValuesBehavior <- vector()
    
    pValuesBehaveFromMixed <- vector()
    pValuesSexFromMixed <- vector()
    pValuesBehaveGroupFromMixedInteraction <- vector()
    pValuesBehaveSexFromMixedInteraction <- vector()
    pValuesGroupSexFromMixedInteraction <- vector()
    pValues3WayInteraction <- vector()
    pValuesCage <- vector()
    pValuesCagenoBehave <- vector()
    behavz <- vector()
    pVals <- list()
    
    ## Ignores the metadata
    for ( j in (endMetadataIndex + 1):ncol(subT)) {
      ## 66 is the first true case here, then 68
      ## par(mfrow=c(2,2), mgp=c(3, 1.25, 0), oma = c(0, 0, 2, 0))
      bug <- subT[,j]
      ## Check to see if these are problematic with the switchover to nlme
      if( sum(bug != 0 ) > nrow(subT) / 4 ) {
        if( (taxa == 1) | 
            (taxa == 2) | # & j != 84) |
            (taxa == 3) | # & j != 88) |
            (taxa == 4) | # & j != 90) |
            (taxa == 5) | # & j != 129) |
            (taxa == 6) | # & j != 161) |
            (taxa == 7) # & j != 144 & j != 183)
        ){
          for ( behavior in behaviorColumns){
            
            ## taxa 5: j != 255
            ## taxa 6: j != 410
            ## taxa 7: j != 490
            ##  | (taxa == 4 & j != 90)
            ## pValuesBehave[index] <- wilcox.test( subT[subT$Group=="Experimental", j],
            ##                                   subT[subT$Group=="Control", j])$p.value
            names[index] <- names(subT)[j]
            behavz[index] <- colnames(subT)[behavior]
            behave <- as.numeric(subT[,behavior])
            ## group <- factor(subT$Group)
            sex <- factor(subT$Sex)
            
            cage <-  factor(subT$Cage, c("female #5,6,7,8 in same cage", "female #13,14,15,16 in same cage", "female #1,2,3,4 in same cage", "female #9,10,11,12 in same cage",
                                         "male #5,6,7,8 in same cage", "male #13,14,15,16 in same cage", "male #1,2,3,4 in same cage", "male #9,10,11,12 in same cage"))
            
            group <- factor(subT$Group)
            #depth <- as.numeric(subT$sampleDepth)
            
            myFrame <- data.frame(bug, sex, group, behave, cage)
            ## These should be replaced with whatever is current in BLJ
            ##pValuesSex[index] <- wilcox.test( subT[subT$Sex=="male", j],
            ##                                  subT[subT$Sex=="female", j])$p.value
            ##pValuesGroup[index] <- wilcox.test( subT[subT$Group=="Control", j],
            ##                                    subT[subT$Group=="Experimental", j])$p.value
            pValuesSex[index] <- pvalue(wilcox_test(bug~sex, data=myFrame))
            pValuesGroup[index] <- pvalue(wilcox_test(bug~group, data=myFrame))
            pValuesBehavior[index] <- cor.test( bug, behave, method="kendall")$p.value[1]
            
            # fullModelwithDepth <-
            #   gls( bug~  behave * group * sex + depth, method="REML",correlation=corCompSymm(form=~1|factor(cage),),data = myFrame )
            # 
            # fullModelNocwithDepth <- gls(bug~  behave * group * sex + depth, data = myFrame)
            # fullModelLME <- lme(bug~ behave*group*sex + depth, method="REML", random = ~1|cage, data = myFrame)
            # 
            ## TODO Did I use depth in the paper itself?
            ## I was able to drop it from consideration.
            # noBehaveModel <- lme(bug ~ group*sex, method="REML", random=~1|cage, data = myFrame)
            # noBehavenoCageModel <- gls(bug ~ group*sex, data = myFrame)
            # pValuesCage[index] <-  anova(fullModelLME, fullModelNocwithDepth)$"p-value"[2]
            # pValuesCagenoBehave[index] <- anova(noBehaveModel, noBehavenoCageModel)$"p-value"[2]
            ## DONE Add in the model from the paper itself, but as nlme
            pVals[[index]] <- c(pValuesSex[index], pValuesGroup[index], pValuesBehavior[index])       
            index <- index + 1
          }
        }
      }
    }
    
    if (otuTable == TRUE) {
      namesMapping <- otuMapping[names, 2]
      dFrame <- data.frame(names, namesMapping, behavz, do.call(rbind, pVals))
      
      colnames(dFrame) <- c("taxonName", "namesMapping", "behavior",
                            "wilcoxonSex", "wilcoxonGroup", "kendallBehavior")
      ## print(colnames(dFrame))
      adjOffset <- 4
      rownames(dFrame) <- paste0(dFrame$taxonName, "_", dFrame$behavior)
    } else {
      dFrame <- data.frame(names, behavz, do.call(rbind, pVals))
      colnames(dFrame) <- c("taxonName", "behavior",
                            "wilcoxonSex", "wilcoxonGroup", "kendallBehavior")
      adjOffset <- 3
    }
    
    #dFrame <- data.frame(names, behavz, do.call(rbind, pVals))
    #colnames(dFrame) <- c("taxonName", "behaviorName",
    #                      "wilcoxonSex", "wilcoxonGroup", "kendallBehavior")
    ## dFrame <- dFrame [order(dFrame$pValuesBehaveFromMixed),]
    numCols <- ncol(dFrame)
    for (i in adjOffset:numCols){
      dFrame[,i + numCols - (adjOffset - 1)] <- p.adjust(dFrame[,i], method = "BH")
      colnames(dFrame)[i + numCols - (adjOffset - 1)] <- paste0("adjusted_",colnames(dFrame)[i])
    }
    
    bestAt <- rbind(bestAt, dFrame[which(dFrame[grep("adjust", colnames(dFrame))] < 0.05, arr.ind=TRUE)[,1],])
    
    ## rownames(dFrame) <- dFrame$taxonName
    
    write.table(dFrame, file=paste(dataType, "_", tissues, "_analysisKendallBehaviorSTRICT_", taxa, "_", colnames(myT)[behavior], ".txt",sep=""), sep="\t",row.names=FALSE)
    
    bestAt <- unique(bestAt)
    write.table(bestAt, file=paste(dataType, "_", tissues, "_SigHitsatKendallSTRICTL_", taxa, ".txt",sep=""), sep="\t",row.names=FALSE)
    
    pdf( paste(dataType, "_KendallBehavior_L", taxa, "_tissue_", t, "_divider_", divider, ".pdf", sep=""))
    
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
      
      boxplot( bug ~ group, data = myFrame,
               main = paste("Wilcoxon group p-value", format(dFrame[nameIter,]$"adjusted_wilcoxonGroup",digits=3) ), ylab="Log normalized abundance" )
      stripchart(bug ~ group,
                 data = myFrame,vertical = TRUE, pch = 21, add=TRUE)
      
      boxplot( bug ~ sex,
               main = paste("Wilcoxon sex p-value", format(dFrame[nameIter,]$"adjusted_wilcoxonSex",digits=3) ), ylab="Log normalized abundance" )
      stripchart(bug ~ sex,
                 data = myFrame, vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])
      
      mtext(paste0(dFrame[nameIter,]$taxonName), outer=TRUE, cex = 0.7)
      
      boxplot( bug ~ group,
               main = paste("MLM group p-value", format(dFrame[nameIter,]$"adjusted_group",digits=3) ), ylab="Log normalized abundance" )
      stripchart(bug ~ group,
                 data = myFrame,vertical = TRUE, pch = 21, add=TRUE)
      
      boxplot( bug ~ sex,
               main = paste("MLM sex p-value", format(dFrame[nameIter,]$"adjusted_sex",digits=3) ), ylab="Log normalized abundance" )
      stripchart(bug ~ sex,
                 data = myFrame, vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])
      
      mtext(paste0(dFrame[nameIter,]$taxonName), outer=TRUE, cex = 0.7)
      
      boxplot( bug ~ group*sex, ylab="Log normalized abundance",
               main=paste("group*sex interaction p-value", format(dFrame[nameIter,]$"adjusted_group:sex",digits=3)),
               xaxt='n')
      axis(1, at=c(1, 2, 3, 4), labels = c("Control\nFemale", "Exp.\nFemale", "Control\nMale", "Exp.\nMale"), cex.axis=0.9)
      stripchart(bug ~ group*sex,
                 data = myFrame,vertical = TRUE, pch = 21, add=TRUE )
      
      boxplot( bug ~ cage ,las=2, ylab="Log normalized abundance",
               main=paste("cage effect p-value", format(dFrame[nameIter,]$"adjusted_cage",digits=3)), xaxt='n')
      axis(1, at=c(1, 3, 5, 7),
           labels=c("FC1", "FX1", "MC1", "MX1"), cex.axis=0.8)
      axis(1, at=c(2, 4, 6, 8),
           labels=c("FC2", "FX2", "MC2", "MX2"), cex.axis=0.8)
      stripchart(bug ~ cage,
                 data = myFrame,vertical = TRUE, pch = 21, add=TRUE)
      # This might have to change as well.
      mtext(paste0(dFrame[nameIter,]$taxonName), outer=TRUE, cex = 0.7)
      
      ##         ##             boxplot( bug ~ behave*sex, ylab="Log normalized abundance",
      ##         ##         main=paste("sex*behave interaction p-value", format(dFrame$adjustedpValuesSexBehaveFromMixedInteraction[index],digits=3)),
      ##         ##         xaxt='n')
      ##         ## axis(1, at=c(1, 2, 3, 4), labels = c("Control\nFemale", "Exp.\nFemale", "Control\nMale", "Exp.\nMale"), cex.axis=0.9)
      
      ##         ## stripchart(bug ~ behave*sex,
      ##         ##            data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index] )
      
      ##             mtext(paste0(otuMapping[nameIter,2],"\n",nameIter), outer=TRUE, cex = 0.7)
      
      ##             index = index + 1
      ##         }
    }
    ## This is the selection from the list for however many
    ## This should be made nicer
    if (otuTable == TRUE) {
      iterOver <- c(4:9)
    } else {
      iterOver <- c(3:8)
    }
    for(histIter in iterOver) {
      hist(dFrame[,histIter], breaks=20, main=colnames(dFrame)[histIter])
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
      
      boxplot( bug ~ group, data = myFrame,
               main = paste("Wilcoxon group p-value", format(bestAt[nameIter,]$"adjusted_wilcoxonGroup",digits=3) ), ylab="Log normalized abundance" )
      stripchart(bug ~ group,
                 data = myFrame,vertical = TRUE, pch = 21, add=TRUE)
      
      boxplot( bug ~ sex,
               main = paste("Wilcoxon sex p-value", format(bestAt[nameIter,]$"adjusted_wilcoxonSex",digits=3) ), ylab="Log normalized abundance" )
      stripchart(bug ~ sex,
                 data = myFrame, vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])
      
      mtext(paste0(bestAt[nameIter,]$taxonName), outer=TRUE, cex = 0.7)
      
      boxplot( bug ~ group,
               main = paste("MLM group p-value", format(bestAt[nameIter,]$"adjusted_group",digits=3) ), ylab="Log normalized abundance" )
      stripchart(bug ~ group,
                 data = myFrame,vertical = TRUE, pch = 21, add=TRUE)
      
      boxplot( bug ~ sex,
               main = paste("MLM sex p-value", format(bestAt[nameIter,]$"adjusted_sex",digits=3) ), ylab="Log normalized abundance" )
      stripchart(bug ~ sex,
                 data = myFrame, vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])
      
      mtext(paste0(bestAt[nameIter,]$taxonName), outer=TRUE, cex = 0.7)
      
      boxplot( bug ~ group*sex, ylab="Log normalized abundance",
               main=paste("group*sex interaction p-value", format(bestAt[nameIter,]$"adjusted_group:sex",digits=3)),
               xaxt='n')
      axis(1, at=c(1, 2, 3, 4), labels = c("Control\nFemale", "Exp.\nFemale", "Control\nMale", "Exp.\nMale"), cex.axis=0.9)
      stripchart(bug ~ group*sex,
                 data = myFrame,vertical = TRUE, pch = 21, add=TRUE )
      
      boxplot( bug ~ cage ,las=2, ylab="Log normalized abundance",
               main=paste("cage effect p-value", format(bestAt[nameIter,]$"adjusted_cage",digits=3)), xaxt='n')
      axis(1, at=c(1, 3, 5, 7),
           labels=c("FC1", "FX1", "MC1", "MX1"), cex.axis=0.8)
      axis(1, at=c(2, 4, 6, 8),
           labels=c("FC2", "FX2", "MC2", "MX2"), cex.axis=0.8)
      stripchart(bug ~ cage,
                 data = myFrame,vertical = TRUE, pch = 21, add=TRUE)
      
      mtext(paste0(bestAt[nameIter,]$taxonName), outer=TRUE, cex = 0.7)
      
      ##         ##             boxplot( bug ~ behave*sex, ylab="Log normalized abundance",
      ##         ##         main=paste("sex*behave interaction p-value", format(dFrame$adjustedpValuesSexBehaveFromMixedInteraction[index],digits=3)),
      ##         ##         xaxt='n')
      ##         ## axis(1, at=c(1, 2, 3, 4), labels = c("Control\nFemale", "Exp.\nFemale", "Control\nMale", "Exp.\nMale"), cex.axis=0.9)
      
      ##         ## stripchart(bug ~ behave*sex,
      ##         ##            data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index] )
      
      ##             mtext(paste0(otuMapping[nameIter,2],"\n",nameIter), outer=TRUE, cex = 0.7)
      
      ##             index = index + 1
      ##         }
    }
    dev.off()
    
    
    pdf( paste(dataType, "_", tissues, "_KendallBehavior_", taxa, ".pdf", sep=""))
    
    hist(dFrame$kendallBehavior, breaks=20, main="")
    
    dev.off()
    }
}