library(shiny)
library(ggplot2)

# Define UI for dataset viewing
shinyUI(pageWithSidebar(
  # Application title
  headerPanel("RNA-seq data visualization"),


    # Sidebar with controls to select the variable to plot
    sidebarPanel(
        selectInput("time", "Time Range:", list("0-15 induced" = "0-15i", "0-90 induced" = "0-90i", "0-90 uninduced" = "0-90u", "90-90 uninduced" = "90i-90u"),'0-15 induced'),
        radioButtons("tf", "Transcription Factor:", c("MIG1" = "MIG1/n", "MSN2" = "MSN2/m", "RAP1" = "RAP1/r"),'MSN2'),
        radioButtons('protein', 'Dataset', c('ATF'='A', 'ZEV'='Z'), 'ATF')
    ),

  mainPanel(
    tabsetPanel(
        tabPanel("Subsets",
            plotOutput('fullPlot'),
            verbatimTextOutput('filenameTF'),
            wellPanel(
                fileInput('file1', 'Upload gene list', accept=c('text/tsv', 'text/tab-separated-values,text/plain', '.tsv')),
                tags$hr(), 
                checkboxInput('pval', 'Wilcoxon Test', FALSE), 
                textInput('gene','Gene:')
                # TODO Keep adding genes and have a clear button if possible.
                )
            ),
                
        tabPanel("Promoter Analysis",
# SHOW Only tfs that the gene has, against average (or cummulative) TF counts. 
# USE the gene set to look at gene promoters.
                 plotOutput("promoterHits"),
                 wellPanel(
"The threshold occupancy for these TFs is 0.25", textInput('geneTFhist','Gene:'),
checkboxInput('onlyMMR','Only show MIG1, MSN2, and RAP1:',FALSE),
fileInput('fileTFGenelist', 'Upload Gene list', accept=c('text/tsv', 'text/tab-separated-values,text/plain', '.tsv'))
                           )),

        # TODO TF COUNTS TAB with a histogram against the mean...
        tabPanel("Promoter Table", 
                # Individual 
                 dataTableOutput("promIndTable"),
                 wellPanel(
checkboxInput('onlyThree','Only show MIG1, MSN2, and RAP1:',FALSE),
textInput('genePromoter','Gene:'),
radioButtons("threshold", "Threshold:", c("0.025!"="0.025","0.05!"="0.05","0.1!"="0.1","0.25!"="0.25","0.5!"="0.5"),"0.25!"),
fileInput('fileTFlist', 'Upload TF list', accept=c('text/tsv', 'text/tab-separated-values,text/plain', '.tsv'))
                 )),


        # TODO SIDE BY SIDE COMPARISONS OF THE GENE IN DIFF DATASETS, over time.
        tabPanel("By Gene", verbatimTextOutput("summary")),

        tabPanel("Comparison", 
            plotOutput('againstPlot'),
            wellPanel(
            selectInput("time2", "Time Range:", list("0-15 induced" = "0-15i", "0-90 induced" = "0-90i", "0-90 uninduced" = "0-90u", "90-90 uninduced" = "90i-90u"),'0-90 induced'),
            radioButtons("tf2", "Transcription Factor:", c("MIG1" = "MIG1/n", "MSN2" = "MSN2/m", "RAP1" = "RAP1/r"),'MSN2'),
            radioButtons('protein2', 'Dataset', c('ATF'='A', 'ZEV'='Z'), 'ATF'),
            textInput('geneDiff','Gene:'),
            fileInput('file2', 'Upload gene list', accept=c('text/tsv', 'text/tab-separated-values,text/plain', '.tsv')),
                tags$hr())
            )
        )
    )
))

# TODO add a Go term menu where you can pick several GO sets.
