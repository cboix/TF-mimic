#!/usr/bin/R
# Want to check all combinations, but lets start w/ MIG1 + iZ0-90i
motif <- 'MIG1'
ft <- read.table('/home/carles/QCB301/MIG1/iZ0-90i',sep="\t",header=T)
number <- 55


counts <- read.table('/home/carles/QCB301/promoters/FINALcounts_2.tsv')
a <- counts[counts$TF == motif,c(-1,-6169,-6170)]
dtf <- data.frame(gene = names(a),count = as.numeric(as.character(unlist(a))))

# Now take only the genes that have a MIG1 target. merge them w/ FT:
dtf <- dtf[dtf$count > 0,]
topMot <- merge(dtf,ft)
topMot <- topMot[topMot$log2.fold_change != Inf,]
topMot <- topMot[topMot$log2.fold_change != -Inf,]
topMot <- topMot[order(topMot$significant,abs(topMot$log2.fold_change.),decreasing=T),]

out <- as.character(head(topMot$gene,number))

system(paste('rm -f /home/carles/QCB301/MEME_iZ090i_',number,'_',motif,sep=""))
for(n in out){ 
    x <- paste('cat /home/carles/QCB301/promoters/',n,'.fasta >> /home/carles/QCB301/MEME_iZ090i_',number,'_',motif,sep="")
    system(x)
}
