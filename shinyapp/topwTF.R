#!/usr/bin/R
OTG <- read.table('/home/carles/QCB301/ReferenceGenomesandAnnotations/ORFtoGENE',sep=" ",header=T)
pDIR='/home/carles/QCB301/promoters/'

# Want to check all combinations, but lets start w/ MIG1 + iZ0-90i
motif <- 'MIG1'
ft <- read.table('/home/carles/QCB301/MIG1/iZ0-90i',sep="\t",header=T)
number <- 55

process_logdtf <- function(dtf){
    dtf <- dtf[dtf$log2.fold_change != Inf,]
    dtf <- dtf[dtf$log2.fold_change != -Inf,]
    dtf <- dtf[order(dtf$significant,abs(dtf$log2.fold_change.),decreasing=T),]
    return(dtf)
}


counts <- read.table('/home/carles/QCB301/promoters/FINALcounts_2.tsv')
a <- counts[counts$TF == motif,c(-1,-6169,-6170)]
dtf <- data.frame(gene = names(a),count = as.numeric(as.character(unlist(a))))

# Now take only the genes that have a MIG1 target. merge them w/ FT:
dtf <- dtf[dtf$count > 0,]
topMot <- process_logdtf(merge(dtf,ft))

out <- as.character(head(topMot$gene,number))
outFILE <- paste('/home/carles/QCB301/MEME_iZ090i_',number,'_',motif,sep="")
system(paste('rm -f',outFILE))
for(n in out){ 
    x <- paste('cat ',pDIR,n,'.fasta >> ',outFILE,sep="")
    system(x)
}


# The same, but with the MIG1 targets instead!
mg <- read.table('/home/carles/QCB301/Lists/MIG1_targets')
ft <- read.table('/home/carles/QCB301/MIG1/iZ0-90i',sep="\t",header=T)
# Take top 25 no matter what?
number <- 25
# Now take only the genes that have a MIG1 target. merge them w/ FT:
names(mg) <- 'gene'
topMot <- process_logdtf(merge(mg,ft))

out <- as.character(head(topMot$gene,number))
outFILE <- paste('/home/carles/QCB301/MEME_iZ090i_',number,'_MIGtargets',sep="")
system(paste('rm -f',outFILE))
for(n in out){ 
    fname <- paste(pDIR,n,'.fasta',sep="")
    if(!file.exists(fname)) {
        O_n <- OTG[OTG$gene == n,]$ORFname
        fname <- paste(pDIR,O_n,'.fasta',sep="")
    }
        x <- paste('cat',fname,'>>',outFILE)
    system(x)
}

# The same, but with the MSN2 targets instead!
mg <- read.table('/home/carles/QCB301/Lists/MSN2_targets')
ft <- read.table('/home/carles/QCB301/MSN2/mZ0-90i',sep="\t",header=T)
# Take top 25 no matter what?
number <- 50
# Now take only the genes that have a MSN2 target. merge them w/ FT:
names(mg) <- 'gene'
topMot <- process_logdtf(merge(mg,ft))

out <- as.character(head(topMot$gene,number))
outFILE <- paste('/home/carles/QCB301/MEME_mZ090i_',number,'_MSNtargets',sep="")
system(paste('rm -f',outFILE))
for(n in out){ 
    fname <- paste(pDIR,n,'.fasta',sep="")
    if(!file.exists(fname)) {
        O_n <- OTG[OTG$gene == n,]$ORFname
        fname <- paste(pDIR,O_n,'.fasta',sep="")
    }
        x <- paste('cat',fname,'>>',outFILE)
    system(x)
}


# The same, but with the MIG1 targets in ATF!!
mg <- read.table('/home/carles/QCB301/Lists/MIG1_targets')
ft <- read.table('/home/carles/QCB301/MIG1/iA0-90i',sep="\t",header=T)
# Take top 25 no matter what?
number <- 50
# Now take only the genes that have a MIG1 target. merge them w/ FT:
names(mg) <- 'gene'
topMot <- process_logdtf(merge(mg,ft))

out <- as.character(head(topMot$gene,number))
outFILE <- paste('/home/carles/QCB301/MEME_iA090i_',number,'_MIGtargets',sep="")
system(paste('rm -f',outFILE))
for(n in out){ 
    fname <- paste(pDIR,n,'.fasta',sep="")
    if(!file.exists(fname)) {
        O_n <- OTG[OTG$gene == n,]$ORFname
        fname <- paste(pDIR,O_n,'.fasta',sep="")
    }
        x <- paste('cat',fname,'>>',outFILE)
    system(x)
}


# The same, but with the MSN2 targets in ATF!
mg <- read.table('/home/carles/QCB301/Lists/MSN2_targets')
ft <- read.table('/home/carles/QCB301/MSN2/mA0-90i',sep="\t",header=T)
# Take top 25 no matter what?
number <- 50
# Now take only the genes that have a MSN2 target. merge them w/ FT:
names(mg) <- 'gene'
topMot <- process_logdtf(merge(mg,ft))

out <- as.character(head(topMot$gene,number))
outFILE <- paste('/home/carles/QCB301/MEME_mA090i_',number,'_MSNtargets',sep="")
system(paste('rm -f',outFILE))
for(n in out){ 
    fname <- paste(pDIR,n,'.fasta',sep="")
    if(!file.exists(fname)) {
        O_n <- OTG[OTG$gene == n,]$ORFname
        fname <- paste(pDIR,O_n,'.fasta',sep="")
    }
        x <- paste('cat',fname,'>>',outFILE)
    system(x)
}

