read.table('/home/carles/QCB301/MIG1prom',sep=" ",header=F) -> mp
read.table('/home/carles/QCB301/MSN2prom',sep=" ",header=F) -> np

MP <- mp[mp$V1 == 'YGL035C',c(2,3,4,5,9)]
names(MP) <- c('gene','start','end','occ','val')

NP <- np[np$V5 >= 0.2,c(2,3,4,5,9)]
names(NP) <- c('gene','start','end','occ','val')

library(ggplot2)
ggplot(MP,aes(start)) + geom_histogram(binwidth=10) + labs(x='MIG1 Binding Site Position',y='Count')
ggsave('~/QCB301/Images/MIGpl1.png',dpi=600,units='in',width=5,height=3)
dev.new()
ggplot(NP,aes(start)) + geom_histogram(binwidth=10)
# MSN2 seems to be picking up too much noise in motifs. I may need a better pwm... 


# Now look at exp. levels! in ATF!
exp.mig <- read.table('/home/carles/QCB301/MIG1/iFA0-90i',header=T,sep="\t")
exploc.mig <- merge(MP,exp.mig)
# only get ~ 5500 genes.
ggplot(exploc.mig,aes(start,log2.fold_change.)) + geom_point(alpha=.5)

# lets look at the rolling mean?
library(zoo)
findwindowEXP <- function(dtf,K){
    meanO <- rollmean(dtf$log2.fold_change.,k=K)
    sdO <- rollapply(dtf$log2.fold_change.,FUN=sd,width=K)
    out <- data.frame(position=dtf[K:nrow(dtf),1],mean=meanO,min=meanO - sdO,max=meanO + sdO)
}
findEXP2 <- function(dtf,K){
    dtf <- dtf[,c(2,10)]
    dtf <- dtf[order(dtf),]
    dtf <- na.omit(dtf)
    dtf <- dtf[dtf$log2.fold_change. != Inf,]
    dtf <- dtf[dtf$log2.fold_change. != -Inf,]
    return(findwindowEXP(dtf,K))
}

elmout <- findEXP2(exploc.mig,500)
ggplot(elmout,aes(position,mean))  + geom_point(alpha=.5,color='black') + geom_ribbon(aes(x=position,ymax=max,ymin=min),alpha=.1) + labs(x='Away from Gene\t\t\t\t \t\t\t\t Motif Position  \t\t\t\t \t\t\t\t Towards Gene',y = 'Mean Fold Change')
# ggplot(elm500,aes(position,mean)) + geom_point(alpha=.5) + geom_point(data=elm100,aes(position,mean),color="red",alpha=.5) + geom_point(data=elm1000,aes(position,mean),color="blue",alpha=.5)
# Both show a pattern of going up:
ggsave('~/QCB301/Images/MIG1Amotifplacement500.png',dpi=600,units='in',width=5,height=3)


msn2t <- read.table('QCB301/Lists/MSN2_targets'); names(msn2t) <- 'gene'
# -MSN2=-----------------------Same work:
# Now look at exp. levels! in ATF!
exp.mig <- read.table('/home/carles/QCB301/MSN2/mFA0-90i',header=T,sep="\t")
exp.mig <- merge(exp.mig,msn2t)
exploc.mig <- merge(NP,exp.mig)
# only get ~ 5500 genes.
# ggplot(exploc.mig,aes(start,log2.fold_change.)) + geom_point(alpha=.5)
elmout <- findEXP2(exploc.mig,500)
ggplot(elmout,aes(position,mean))  + geom_point(alpha=.5,color='blue') + geom_ribbon(aes(x=position,ymax=max,ymin=min),alpha=.1) + labs(x='Away from Gene\t\t\t\t \t\t\t\t Motif Position  \t\t\t\t \t\t\t\t Towards Gene',y = 'Mean Fold Change')

# Now look at exp. levels! in ZEV!
exp.mig <- read.table('/home/carles/QCB301/MIG1/iZ0-90i',header=T,sep="\t")
exploc.mig <- merge(MP,exp.mig)
elmout <- findEXP2(exploc.mig,500)
ggplot(elmout,aes(position,mean))  + geom_point(alpha=.5,color='black') + geom_ribbon(aes(x=position,ymax=max,ymin=min),alpha=.1) + labs(x='Away from Gene\t\t\t\t \t\t\t\t Motif Position  \t\t\t\t \t\t\t\t Towards Gene',y = 'Mean Fold Change')
ggsave('~/QCB301/Images/MIG1Zmotifplacement500.png',dpi=600,units='in',width=5,height=3)

# Now look at exp. levels! in ZEV!
exp.mig <- read.table('/home/carles/QCB301/MSN2/mZ0-90i',header=T,sep="\t")
#exp.mig <- merge(exp.mig,msn2t)
exploc.mig <- merge(MP,exp.mig)
elmout <- findEXP2(exploc.mig,1000)
ggplot(elmout,aes(position,mean))  + geom_point(alpha=.5,color='blue') + geom_ribbon(aes(x=position,ymax=max,ymin=min),alpha=.1)

# ggplot(elm500,aes(position,mean)) + geom_point(alpha=.5) + geom_point(data=elm100,aes(position,mean),color="red",alpha=.5) + geom_point(data=elm1000,aes(position,mean),color="blue",alpha=.5)

mf <- fc[fc$TF == 'MIG1',]
mf <- mf[c(2:(length(mf)-2))]
t(mf) -> dtf
dtf <- data.frame(gene=rownames(dtf),mc = as.numeric(as.character(dtf)))
m <- read.table('QCB301/MIG1/iFA0-90i',sep="\t",header=T)
merge(m,dtf) -> a
a <- a[ a$log2.fold_change. != -Inf,]
a <- a[ a$log2.fold_change. != Inf,]
a <- a[ order(abs(a$log2.fold_change.),decreasing=T),]
a <- data.frame(count=a$mc[1:5092])
a$rank <- 1:nrow(a)
b <- a[1:1500,]

findwindowMOT <- function(dtf,K){
    meanO <- rollmean(dtf$count,k=K)
    sdO <- rollapply(dtf$count,FUN=sd,width=K)
    out <- data.frame(position=dtf[K:nrow(dtf),2],mean=meanO,min=meanO - sdO,max=meanO + sdO)
}
aout <- findwindowMOT(b,100)
ggplot(aout,aes(position,mean))  + geom_point(alpha=.5,color='blue') + geom_ribbon(aes(x=position,ymax=max,ymin=min),alpha=.1) + labs(x='Rank: Absolute value of fold change from 0 to 90 minutes',y = 'Average number of Mig1 sites')
ggsave('QCB301/Images/MIGsitesRM.png',dpi=600,height=3.5,width=5,units='in')

#-------------------

mf <- fc[fc$TF == 'MSN2',]
mf <- mf[c(2:(length(mf)-2))]
t(mf) -> dtf
dtf <- data.frame(gene=rownames(dtf),mc = as.numeric(as.character(dtf)))
m <- read.table('QCB301/MSN2/mFA0-90i',sep="\t",header=T)
merge(m,dtf) -> a
a <- a[ a$log2.fold_change. != -Inf,]
a <- a[ a$log2.fold_change. != Inf,]
a <- a[ order(abs(a$log2.fold_change.),decreasing=T),]
a <- data.frame(count=a$mc[1:5092])
a$rank <- 1:nrow(a)
b <- a[1:1500,]

findwindowMOT <- function(dtf,K){
    meanO <- rollmean(dtf$count,k=K)
    sdO <- rollapply(dtf$count,FUN=sd,width=K)
    out <- data.frame(position=dtf[K:nrow(dtf),2],mean=meanO,min=meanO - sdO,max=meanO + sdO)
}
aout <- findwindowMOT(b,200)
ggplot(aout,aes(position,mean))  + geom_point(alpha=.5,color='blue') + geom_ribbon(aes(x=position,ymax=max,ymin=min),alpha=.1) + labs(x='Rank: Absolute value of fold change from 0 to 90 minutes',y = 'Average number of Msn2 sites')
ggsave('QCB301/Images/MSNsitesRMA.png',dpi=600,height=3.5,width=5,units='in')

