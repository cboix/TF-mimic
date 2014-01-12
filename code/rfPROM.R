#!/usr/bin/R
fc <- read.table('/home/carles/QCB301/promoters/FINALcounts_2.tsv')
FC <- fc[,c(-1,-6169)]
FC <- t(FC)
FC <- data.frame(FC)
names(FC) <- fc$TF

stor.Nam <- rownames(FC)
# TURN it all into classification vars!
FC <- apply(FC,FUN=function(x){if(x != 0){return(1)};return(0)},MARGIN=c(1,2))
FC <- data.frame(FC)
# clean up very non informative rows:
cutOff <- 25
cont <- colSums(FC)
FC <- droplevels(FC[ ,cont >= cutOff])
cont <- colSums(FC)
FC <- droplevels(FC[,cont <= nrow(FC) - cutOff])
# make them all factors!!
# FROM SO: @Thierry:
convert.magic <- function(obj, type){
  FUN1 <- switch(type, character = as.character, numeric = as.numeric, factor = as.factor)
  out <- lapply(obj, FUN1)
  as.data.frame(out)
}
FC <- convert.magic(FC,'factor')
FC$gene <- stor.Nam

# Interchangeable:
gt <- read.table('/home/carles/QCB301/MIG1/iZ0-90i',header=T,sep="\t")
names(gt)[6] <- 'LogFC'
if (intersect(names(gt),names(FC)) == 'gene'){
    FCL <- merge(FC,gt[,c(1,6)])
    } else { 
        print("ERROR in column names")
}
# Bin the data! Lets do this by SD!
FCL <- FCL[FCL$LogFC != Inf,]
FCL <- FCL[FCL$LogFC != -Inf,]

meanLG <- mean(FCL$LogFC)
sdLG <- sd(FCL$LogFC)
maxLG <- max(FCL$LogFC)
minLG <- min(FCL$LogFC)
# Equal width cuts
segments <- 10
# 10 sd on either side is probably enough? -- 20 bins?
br <- seq(-10, 10, by=1)
br <- br * sdLG + meanLG
br <- c(-Inf,br,Inf)
cutLG <- cut(FCL$LogFC, breaks=br , include.lowest=T)
FCL$LogFC <- cutLG

# subset correctly by dropping levels.
smFCL <- droplevels(FCL[,c(2:50,ncol(FCL))])


print(levels(FCL$LogFC))

# Trying w/ simple RF:
library(randomForest)
rf1sm <- randomForest( LogFC ~ ., data=smFCL)

# lets try this? dropping the gene column...
rf1 <- randomForest( LogFC ~ ., data=droplevels(FCL[,-1]))

dtf.rf <- data.frame(TF=rownames(importance(rf1)),importance=c(rf1$importance))
dtf.rf <- dtf.rf[order(dtf.rf$importance),]

# library(party)
# # NEED TO TRY TO GET RID OF HIGHLY CORRELATED VARIABLE BIAS.
# data.controls <- cforest_unbiased(ntree=1000, mtry=round(sqrt(ncol(smFCL))) )
# data.cforest <- cforest(LogFC ~ .,data=smFCL,controls=data.controls)
# data.cforest.varimp <- varimp(data.cforest, conditional = TRUE)
# data.cforest.varimp
