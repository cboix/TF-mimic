library(shiny)
library(ggplot2)

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
        if(input$pval){ p <- p + geom_text(data=NULL,x=5,y=0.4,label="") + labs(title=paste("P-value: ",formatC(wxc$p.value, width=3, flag="0")))}}
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

        p <- ggplot(tt,aes(log2_tf,log2_tf2)) + geom_point(alpha=.4) + labs(x=paste(input$time,'for',strsplit(input$tf,"/")[[1]][1]),y=paste(input$time2,'for',strsplit(input$tf2,"/")[[1]][1]))

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
                p <- p[p$TF %in% c('MIG1','RAP1','MSN2'),]
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

        # No need to check for existence? Lets  just not move this file :)
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
                tfscut <- tfs[tfs$TF %in% c('MIG1','RAP1','MSN2'),]
                TFScut <- TFS[TFS$TF %in% c('MIG1','RAP1','MSN2'),]
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

})
