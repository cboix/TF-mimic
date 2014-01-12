library(shiny)
library(ggplot2)
library(reshape)
library(seqLogo)

# Define server logic required to plot various variables against TFs
shinyServer(function(input, output) {
    output$fullPlot <- renderPlot({
        tfIN <- paste('~/QCB301/',input$tf,input$protein,input$time,sep="")
        if (!file.exists(tfIN))
            return(NULL)
        tfIN <- read.table(tfIN,sep="\t",header=T)
        tfIN <- tfIN[order(tfIN$gene),]
        p <- ggplot(tfIN,aes(log2.fold_change.)) + geom_density(fill='grey',alpha=.5) + labs(x="Log2 - Fold Change",y="Density")
        gLC <- c(tfIN[tfIN$gene== paste(input$gene),]$log2.fold_change.);
        if ( length(gLC) == 1 ){ p <- p + geom_vline(xintercept=gLC,color='red',linetype='dashed')}
        inFile <- input$file1
        if (!is.null(inFile)){
            gn <- read.table(inFile$datapath,sep="\t",header=F)
            gn <- data.frame(gene=sort(gn[,1]))
            gn <- merge(gn,tfIN)
            wxc <- wilcox.test(gn$log2.fold_change.,tfIN$log2.fold_change.)
            p <- p + geom_density(data=gn,aes(log2.fold_change.),color="blue",fill="lightblue",alpha=.4) 
        if(input$wxc){ p <- p + geom_text(data=NULL,x=5,y=0.4,label="") + labs(title=paste("P-value: ",formatC(wxc$p.value, width=3, flag="0")))}}
        print(p)
    })
        
    output$motifPlot <- renderPlot({
        tfIN <- paste('~/QCB301/',input$tf,input$protein,input$time,sep="")
        if (!file.exists(tfIN))
            return(NULL)
        tfIN <- read.table(tfIN,sep="\t",header=T)
        tfIN <- tfIN[order(tfIN$gene),]
        tfMin <- tfIN[,c(1,6)]
        dtf <- cbind(tfMin,'All')
        names(dtf) <- c('Gene','Log2FC','Set')

        tsv <- read.table('/home/carles/QCB301/promoters/FINALcounts_2.tsv')
        if (input$highTol)
            tsv <- read.table('/home/carles/QCB301/promoters/FINALcounts.tsv')
        # List of Motifs to look at:
        ex <- gsub(",","','",input$motifComp); ex <- gsub(" ","",ex)
        ex <- paste("motifs <- c('",ex,"')",sep="")
        eval(parse(text=ex))
        # all genes vs. motifs:
        mot <- tsv[tsv$TF %in% motifs,]
        stor.dtf <- dtf$Log2FC
        if(nrow(mot) > 0){
            lMot <- nrow(mot)
            mot2 <- data.frame(cbind(Gene=names(mot[2:6166]),t(mot[,2:6166])))
            names(mot2) <- c('Gene',as.character(mot$TF))
            for(i in 2:(lMot+1)){mot2[,i] <-  as.numeric(as.character(mot2[,i]))}
            # Get genes with the motif
            getGenesM <- function(dtf,i,w){
                hits <- dtf[,i] > 0
                g <- as.character(dtf[hits,1])
                if(!w){ 
                    return(g)
                } else {
                    c <- dtf[hits,i]
                    return(rep(g,times=c))
            } }
            num.topmot <- c()
            if (lMot > 1){ 
                motGenes <- sapply(2:(lMot+1),dtf=mot2,w=input$weightedMotif,getGenesM) 
                for(i in 1:lMot){
                    mdt <- tfMin[tfMin$gene %in% unlist(motGenes[[i]]),]
                    mdt <- cbind(mdt,names(mot2)[i+1])
                    names(mdt) <- c('Gene','Log2FC','Set')
                    print(head(mdt))
                    dtf <- rbind(dtf,mdt)
                    if (length(num.topmot) == 0){ stor.mdt <- mdt$Log2FC; num.topmot <- nrow(mdt) }
                }
            } else {
                motGenes <- getGenesM(mot2,2,input$weightedMotif) 
                mdt <- tfMin[tfMin$gene %in% unlist(motGenes),]
                mdt <- cbind(mdt,names(mot2)[2])
                names(mdt) <- c('Gene','Log2FC','Set')
                print(head(mdt))
                dtf <- rbind(dtf,mdt)
                num.topmot <- nrow(mdt)
                stor.mdt <- mdt$Log2FC;
            }
        } else {mot2 <- c('','no motifs'); stor.mdt <- 0; num.topmot= 0}
        wxcM <- wilcox.test(stor.mdt,stor.dtf)
        # Selection of data + plotting:
        p <- ggplot(dtf,aes(Log2FC,fill=Set,color=Set)) + geom_density(alpha=.5) + labs(x="Log2 - Fold Change",y="Density") + labs(title=paste(names(mot2)[2],' has ',num.topmot,' instances. Wxc p-value = ',wxcM$p.value,sep=""))
        # Accessory genes/list:
            gLC <- c(tfIN[tfIN$gene== paste(input$geneMot),]$log2.fold_change.);
            if ( length(gLC) == 1 ){ p <- p + geom_vline(xintercept=gLC,color='red',linetype='dashed')}
            inFile <- input$fileMot
            if (!is.null(inFile)){
                gn <- read.table(inFile$datapath,sep="\t",header=F)
                gn <- data.frame(gene=sort(gn[,1]))
                gn <- merge(gn,tfIN)
                wxc <- wilcox.test(gn$log2.fold_change.,tfIN$log2.fold_change.)
                p <- p + geom_density(data=gn,aes(log2.fold_change.),color="blue",fill="lightblue",alpha=.4)
            #if(input$wxcMotif){ p <- p + geom_text(data=NULL,x=5,y=0.4,label="") + labs(title=paste("P-value: ",formatC(wxc$p.value, width=3, flag="0")))}
            }
        print(p)
    })

    output$againstPlot <- renderPlot({
        # Gene set 1
        tfIN <- paste('~/QCB301/',input$tf,input$protein,input$time,sep="")
        if (!file.exists(tfIN))
            return(NULL)
        tfIN <- read.table(tfIN,sep="\t",header=T)
        tfIN <- tfIN[order(tfIN$gene),]
        names(tfIN)[6] <- 'log2_tf'

        # Gene set 2
        tfIN2 <- paste('~/QCB301/',input$tf2,input$protein2,input$time2,sep="")
        if (!file.exists(tfIN2))
            return(NULL)
        tfIN2 <- read.table(tfIN2,sep="\t",header=T)
        tfIN2 <- tfIN2[order(tfIN2$gene),]
        names(tfIN2)[6] <- 'log2_tf2'

        tt <- merge(tfIN[,c(1,6)],tfIN2[,c(1,6)])

        p <- ggplot(tt,aes(log2_tf,log2_tf2)) + geom_point(alpha=.4) + labs(x=paste(input$time,'for',strsplit(input$tf,"/")[[1]][1]),y=paste(input$time2,'for',strsplit(input$tf2,"/")[[1]][1])) + scale_x_continuous(lim=c(-15,15)) + scale_y_continuous(lim=c(-15,15))

        gLC <- tt[tt$gene== paste(input$geneDiff),];
        if (length(gLC) == 3){ p <- p + geom_point(data=gLC,aes(x=log2_tf,y=log2_tf2),color='red',cex=5)}

        # Labels a gene set
        inFile2 <- input$file2
        if (!is.null(inFile2)){
            gn <- read.table(inFile2$datapath,sep="\t",header=F)
            gn <- data.frame(gene=sort(gn[,1]))
            gn <- merge(gn,tt)
            p <- p + geom_point(data=gn,aes(x=log2_tf,y=log2_tf2),color='blue',cex=3.5,alpha=.4)
        }

        print(p)
    })

    output$filenameTF <- renderPrint({
        tf <- strsplit(input$tf,"/")
        if(input$protein=='A'){ pr <- 'ATF' } else { pr <- 'ZEV' }
        txt <- paste(tf[[1]][1],pr,'with time range',input$time)
        inFile <- input$file1
        if (!is.null(inFile)){ txt <- paste(inFile[1],'against',txt) }
        print(txt)
    })


    output$promIndTable <- renderDataTable({
        filename <- paste('/home/carles/QCB301/promoters/',input$genePromoter,'.tfbs',sep="")
        if (file.exists(filename)){
            OTG <- read.table('/home/carles/QCB301/ReferenceGenomesandAnnotations/ORFtoGENE',header=T)
            OTG<- OTG[order(OTG$ORFname),]
            names(OTG) <- c('pwm','TF')

            tfbs <- read.table(filename,header=T)
            tfbs <- tfbs[order(tfbs$pwm),]
            tfbs2 <- merge(OTG,tfbs)
            tfbs3 <- tfbs2[tfbs2$occ >= as.numeric(input$threshold),]
            p <- tfbs3[,c(2,4,5,6,10)]

            inFile <- input$fileTFlist
            if (!is.null(inFile)){
                tfList <- as.numeric(as.character(read.table(inFile$datapath,sep="\t",header=F)))
                p <- p[p$TF %in% tfList,]
            }
            # This supercedes the TF file
            if (input$onlyThree){
                p <- p[p$TF %in% c('MIG1b','RAP1b','MIG1','RAP1','MSN2'),]
            }
            names(p) <- c('TF','Start','End','Occupancy','Strand')
            p <- p[order(p$TF),]
        } else { p <- data.frame(Gene=0,Does=0,Not=0,Exist=0) }
        # border case:
        if(nrow(p) == 0){ p <- data.frame(No=0,TFs=0,IN=0,Promoter=0)}
        p
    },options = list(bSortClasses = TRUE,aLengthMenu = c(10,25,40,60,100,1000), iDisplayLength = 10))        


    # PROMOTER ANALYSIS:
    output$promoterHits <- renderPlot({
        OTG <- read.table('/home/carles/QCB301/ReferenceGenomesandAnnotations/ORFtoGENE',header=T)
        OTG<- OTG[order(OTG$ORFname),]
        names(OTG) <- c('pwm','TF')

        # No need to check for existence? Lets just not move this file :)
        TFS <- read.table('/home/carles/QCB301/ReferenceGenomesandAnnotations/AllMotifs.pcounts',header=F)
        TFS <- cbind(TFS,'All')
        names(TFS) <- c('pwm','Count','Type')
        TFS$Prop <- TFS$Count/6167
        TFS <- merge(OTG,TFS)

        # ADD Genelist to TFS
        fileIn <- input$fileTFGenelist
        if (!is.null(fileIn)){
            geneList <- as.character(unlist(read.table(fileIn$datapath,sep="\t",header=F)));
            if (input$onlyMMR){
                working <- data.frame(pwm=c('YGL035C','YMR037C','YNL216W'),Count = 0)
            } else { 
                working <- data.frame(pwm=TFS[,1],Count=0)
            }
            for (name in geneList){
                filename2 <- paste('/home/carles/QCB301/promoters/',name,'.pcounts',sep="")
                if(file.exists(filename2)){
                    gL <- read.table(filename2,header=F)
                    names(gL) <- c('c','pwm')
                    gL2 <-  merge(gL[,c(2,1)],working)
                    gL2$Count <- gL2$Count + gL2$c
                    working <- working[!(working$pwm %in% factor(as.character(gL$pwm))),]
                    working <- rbind(working,gL2[,c(1,3)])
                }
            }
            Gtfs <- cbind(working,fileIn[1])
            names(Gtfs) <- c('pwm','Count','Type')
            Gtfs <- merge(OTG,Gtfs)
            Gtfs$Prop <- Gtfs$Count/length(geneList)
            TFS <- rbind(TFS,Gtfs)
        }

        filename <- paste('/home/carles/QCB301/promoters/',input$geneTFhist,'.pcounts',sep="")
        print(filename)
        if (file.exists(filename)){
            tfs <- read.table(filename,header=F)
            tfs <- cbind(tfs[,c(2,1)],input$geneTFhist)
            names(tfs) <- c('pwm','Count','Type')
            tfs$Prop <- tfs$Count 
            print(tfs)
            tfs <- tfs[order(tfs$pwm),]
            tfs <- merge(OTG,tfs)
            if (input$onlyMMR){
                tfscut <- tfs[tfs$TF %in% c('MIG1b','RAP1b','MIG1','RAP1','MSN2'),]
                TFScut <- TFS[TFS$TF %in% c('MIG1b','RAP1b','MIG1','RAP1','MSN2'),]
                if( nrow(tfscut) == 0){
                    both <- TFScut 
                } else { 
                    both <- rbind(tfscut,TFScut) 
                }
            } else {
                TFScut <- TFS[TFS$TF %in% factor(as.character(tfs$TF)),]
                both <- rbind(tfs,TFScut)
            }
            p <- ggplot(both,aes(TF,Prop,fill=Type)) + geom_histogram(color='black',position="dodge") + labs(x="Transcription Factor",y="Normalized Count") + coord_flip()
        } else { p <- NULL }
        print(p)
    })        


    # PROMOTER ANALYSIS 2 (from table!!)
    output$promoterHits2 <- renderPlot({
        tsv <- read.table('/home/carles/QCB301/promoters/FINALcounts_2.tsv')
        fin <- data.frame(TF=tsv$TF,Count=tsv$All,Type='All',Prop=tsv$All/6165)

        # Add genelist
        fileIn <- input$fileTFGenelist
        if (!is.null(fileIn)){
            geneList <- as.character(unlist(read.table(fileIn$datapath,header=F)));
            idxv <- (names(tsv) %in% geneList) * 1:ncol(tsv)
            gL <- tsv[,idxv]
            dtf <- data.frame(TF=tsv$TF,Count=rowSums(gL),fileIn[1])
            names(dtf) <- c('TF','Count','Type')
            dtf$Prop <- dtf$Count/length(names(gL))
            print(names(dtf))
            print(head(dtf))
            fin <- rbind(fin,dtf)
        }

        # test <- 'RAP1,MSN2, MIG1'
        # ex <- gsub(",","','",test)
        # For the multiple gene list:
        ex <- gsub(",","','",input$geneTFhist)
        ex <- gsub(" ","",ex)
        ex <- paste("genes <- c('",ex,"')",sep="")
        eval(parse(text=ex))
        idx1 <- (names(tsv) %in% genes) * 1:ncol(tsv)
        if (sum(idx1) > 0){
            gL2 <- data.frame(tsv[,idx1])
            if (ncol(gL2) == 1){ names(gL2) <- genes}
            dtf2 <- data.frame(TF=tsv$TF,gL2)
            dtf2 <- melt(dtf2,id.vars='TF')
            names(dtf2) <- c('TF','Type','Count')
            dtf2$Prop <- dtf2$Count
            fin <- rbind(fin,dtf2)
        }

        fin <- fin[fin$Count > 0,]
        # make sure that we use this set of TFs as the minimal set?
        if (input$onlyMMR){
            fin <- fin[fin$TF %in% c('MIG1b','RAP1b','MIG1','RAP1','MSN2'),]

        }
        if (nrow(fin) > 0){
            p <- ggplot(fin,aes(TF,Prop,fill=Type)) + geom_histogram(color='black',position="dodge") + labs(x="Transcription Factor",y="Normalized Count") + coord_flip()
        } else { p <- NULL }
        print(p)
    })        

    output$showLogo <- renderPlot({
        DIR = '/home/carles/QCB301/ReferenceGenomesandAnnotations/1.02/ALIGNED_ENOLOGO_FORMAT_PFMS/'
        OTG <- read.table('/home/carles/QCB301/ReferenceGenomesandAnnotations/ORFtoGENE',header=T)
        OTG<- OTG[order(OTG$ORFname),]
        names(OTG) <- c('pwm','TF')

        motif <- input$motifLogo
        a <- OTG[OTG$TF == motif,]
        pw <- as.character(a$pwm)
        pw <- paste(pw,"_",sep="")
        if (nrow(a) > 0){
            library(seqLogo)
            pflist <- list.files(path=DIR,pattern='*.pfm')
            b <- c()
            for (n in pflist){
                if(length(grep(pw,n)) == 1){
                    b <- n
                    print(n)
                }
            }
            if (length(b) > 0){
                m <- read.table(paste(DIR,b,sep=""))
                p <- seqLogo(m[c(1,4,3,2),-1])
            }
        } else { p <- NULL }
        print(p)
    })

    output$geneTable <- renderDataTable({
        tfIN <- paste('~/QCB301/',input$tf,input$protein,input$time,sep="")
        if (!file.exists(tfIN))
            return(NULL)
        tfIN <- read.table(tfIN,sep="\t",header=T)
        tfOUT <- tfIN[,c(-2,-3)]
        names(tfOUT) <- c('Gene','Reads_1','Reads_2','Log2 Fold Change','P-value','Q-value','Significant?')
        return(tfOUT)
    })
})
