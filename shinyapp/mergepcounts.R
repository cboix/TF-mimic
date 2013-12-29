#!/usr/bin/R
DIR='/home/carles/QCB301/promoters/'

OTG <- read.table('/home/carles/QCB301/ReferenceGenomesandAnnotations/ORFtoGENE',header=T)
OTG<- OTG[order(OTG$ORFname),]
names(OTG) <- c('pwm','TF')

TFS <- read.table('/home/carles/QCB301/ReferenceGenomesandAnnotations/AllMotifs.pcounts',header=F)
names(TFS) <- c('pwm','All')
TFS <- merge(OTG,TFS)


pclist <- list.files(path=DIR,pattern='*.pcounts')
genes <- gsub(".pcounts","",pclist)
dtf <- data.frame(Gene=genes,counts=pclist)

for (i in 1:nrow(dtf)){
    filename2 <- paste(DIR,dtf$counts[i],sep="")
    if(file.exists(filename2)){
            gL <- read.table(filename2,header=FALSE)
            names(gL) <- c(as.character(dtf$Gene[i]),'pwm')
            gL2 <- data.frame(0,as.character(TFS$pwm[!TFS$pwm %in% gL$pwm]))
            names(gL2) <- c(as.character(dtf$Gene[i]),'pwm')
            gL <- rbind(gL,gL2)
            TFS <-  merge(gL[,c(2,1)],TFS,all=TRUE)
            print(as.character(dtf$Gene[i]))
    }
}

write.table(x=TFS,file=paste(DIR,'FINALcounts.tsv',sep="")
