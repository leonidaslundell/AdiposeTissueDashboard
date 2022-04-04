# rm(list=ls())
print(getwd())
load("dat.Rdata")
# source("plotFunctions.R")
# load("Adipose gene expression analysis/interactiveDash/dat.Rdata") |> try()
# source("Adipose gene expression analysis/interactiveDash/plotFunctions.R") |> try()

library(ggpubr)
library(ggplot2)
library(ggiraph)
library(shinydashboard)
library(shiny)
library(data.table)
library(shinycssloaders)

####################################################
#functions

boxPlot <- function(highlight = NULL){
  if(is.null(highlight)){
    g <- ggplot() +
      annotate(geom = "text", label = "Type in you favorite gene in the bar above,\nor click on the plots below to explore a gene", x = .6,y = .5, size = 18) +
      theme_void()
    
    return(g)
  }
  if(length(highlight)>1){
    g <- ggplot() +
      annotate(geom = "text", label = "Select a single gene", x = .6,y = .5, size = 18) +
      theme_void()
    
    return(g)
  }
  g <- ggplot(datBox[Symbol == highlight, ],
              aes(x = Group, y = value, fill = Diagnosis)) +
    geom_boxplot() +
    scale_x_discrete(labels = rep(c("PRE", "POST", "REC"), 2)) +
    scale_y_log10(labels = scales::comma) +
    xlab("") +
    ylab("Expression (counts per million)") +
    theme_bw() +
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 18))
  
  if(nrow(datRes[Symbol == highlight & value != "ns", ])>0){
    sigs <- nrow(datRes[Symbol == highlight & value != "ns", ])
    positions <- datBox[Symbol == highlight, ]$value |> ceiling() |> sort(decreasing = T) |> log10()
    g <- g +
      ggpubr::stat_pvalue_manual(datRes[Symbol == highlight & value != "ns", ],
                                 label = "value",
                                 y.position = if(sigs==1){positions[1]}else{positions[1:sigs]},
                                 label.size = 6, 
                                 bracket.size = 1, bracket.nudge.y = 0.5)
    g
  }
  g <- g + ggtitle(highlight)
  return(g)
}

exp2logFC <- function(aveExp,
                      logFC,
                      sig,
                      Symbol = datMD$Symbol,
                      title = "", 
                      xZoom = c(NA, NA),
                      yZoom = c(NA, NA),
                      highlight = NULL){
  
  y <- data.table(Symbol,
                  aveExp,
                  logFC,
                  sig,
                  stringsAsFactors = F)
  y <- y[!is.na(y$aveExp),]
  
  g <- ggplot(mapping = aes(x = aveExp, y = logFC, tooltip = Symbol, data_id = Symbol)) +
    geom_point(data = y[y$sig == "*",], color = "black", size = 1.4) +
    geom_point_interactive(data = y[y$sig == "ns",], color = "gray", alpha = .5, size = 0.5) +
    geom_point_interactive(data = y[y$sig == "*",], color = "red", alpha = .9, size = 1) +
    geom_hline(yintercept = 0, color = "gray") +
    scale_x_log10() +
    xlab("Expression (counts per million)") + 
    ylab("logFC") +
    # coord_cartesian(xlim = xZoom, ylim = yZoom) +
    theme_bw() +
    theme(text = element_text(size = 14),
          plot.margin = margin(0,0,0,0))
  g
  if(!is.null(highlight)){
    g <- g + geom_point_interactive(data = y[Symbol %in% highlight,], size = 6, color = "blue")
  }
  
  return(g)
}

girafize <- function(g){
  girafe(ggobj = g, 
         options = list(
           opts_sizing(width = .7),
           opts_zoom(max = 5),
           opts_toolbar(saveaspng = F),
           opts_selection(
             type = "single")
         ))
}


####################################################
#UI
ui <- function(){
  
  header <- dashboardHeader(
    title = "Gene expresssion changes in adipose tissue after exercise", titleWidth = 800
  )
  
  side <- dashboardSidebar(disable = T)
  
  body <- dashboardBody({
                          fixedPage(
                            box("Search for you favorite gene",
                                selectizeInput(inputId = "gene",
                                               choices = NULL,
                                               multiple = T,
                                               label = ""),
                                uiOutput("extraInfo"),
                                width = 12),
                            box(width = 12,
                                plotOutput("boxPlot") |> withSpinner()),
                            column(width = 12,
                                   # tags$head(tags$style(HTML('.box {margin: 0px;}'))), 
                                   box("Gene expression immediately after exercise in Non-diabetics", #width = 12,
                                       girafeOutput("NGT_PREvsPOST") |> withSpinner()),
                                   box("Gene expression 3 hours after exercise in Non-diabetics", #width = 12,
                                       girafeOutput("NGT_PREvsREC") |> withSpinner()),
                                   box("Gene expression immediately after exercise in Diabetics", #width = 12,
                                       girafeOutput("T2D_PREvsPOST") |> withSpinner()),
                                   box("Gene expression 3 hours after exercise in Diabetics", #width = 12,
                                       girafeOutput("T2D_PREvsREC") |> withSpinner()) 
                            )
                          )
  })
  
  
  dashboardPage(skin = "black",
    header,
    side,
    body
  )
}

####################################################
#Server

server <- function(input, output, session) {
  
  updateSelectizeInput(inputId = "gene", 
                       choices = datRes$Symbol |> unique() |> sort(decreasing = F), 
                       server = TRUE)

  ####################################################
  #get external info'
  output$extraInfo <- renderUI({
    req(gene$gene)
    stringdb <- fread(paste0("https://string-db.org/api/tsv/get_link?identifiers=",
                             gene$gene,
                             "&species=9606"))[1, url] |> try()
    genecard <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                       gene$gene)
    fluidPage(
      actionButton(inputId = "A", 
                   label = "StringDB", 
                   onclick = paste0("window.open('", stringdb, "', '_blank')")),
      actionButton(inputId = "B",
                   label = "Gene cards",
                   onclick = paste0("window.open('", genecard, "', '_blank')"))
    )
  })
  
  ####################################################
  #main plots
  
  #############
  #NGT pre post
  output$NGT_PREvsPOST <- renderGirafe({
    exp2logFC(aveExp = datMD$NGT_PREvsPOST, 
              logFC = datMD$NGT_PREvsPOST.logFC,
              sig = datMD$NGT_PREvsPOST.sig,
              highlight = gene$gene) |> girafize()
  })
  #############
  #NFT pre rec
  output$NGT_PREvsREC <- renderGirafe({
    exp2logFC(aveExp = datMD$NGT_PREvsREC, 
              logFC = datMD$NGT_PREvsREC.logFC,
              sig = datMD$NGT_PREvsREC.sig,
              highlight = gene$gene) |> girafize()
  })
  #############
  #T2D pre post
  output$T2D_PREvsPOST <- renderGirafe({
    exp2logFC(aveExp = datMD$T2D_PREvsPOST,
              logFC = datMD$T2D_PREvsPOST.logFC,
              sig = datMD$T2D_PREvsPOST.sig,
              highlight = gene$gene) |> girafize()
  })
  #############
  #T2D pre rec
  output$T2D_PREvsREC <- renderGirafe({
    exp2logFC(aveExp = datMD$T2D_PREvsREC,
              logFC = datMD$T2D_PREvsREC.logFC,
              sig = datMD$T2D_PREvsREC.sig,
              highlight = gene$gene) |> girafize()
  })

  #############
  #create a container for the highlight
  gene <- reactiveValues(gene = NULL)
  
  observeEvent(input$NGT_PREvsPOST_selected,{
    gene$gene <- input$NGT_PREvsPOST_selected[length(input$NGT_PREvsPOST_selected)]
    updateSelectizeInput(inputId = "gene",
                         choices = datRes$Symbol |> unique() |> sort(decreasing = F),
                         server = TRUE,
                         selected = NULL)#clears the search bar when clicking
  })
  observeEvent(input$NGT_PREvsREC_selected,{
    gene$gene <- input$NGT_PREvsREC_selected[length(input$NGT_PREvsREC_selected)]
    updateSelectizeInput(inputId = "gene",
                         choices = datRes$Symbol |> unique() |> sort(decreasing = F),
                         server = TRUE,
                         selected = NULL)#clears the search bar when clicking
  })
  observeEvent(input$T2D_PREvsPOST_selected,{
    gene$gene <- input$T2D_PREvsPOST_selected[length(input$T2D_PREvsPOST_selected)]
    updateSelectizeInput(inputId = "gene",
                         choices = datRes$Symbol |> unique() |> sort(decreasing = F),
                         server = TRUE,
                         selected = NULL)#clears the search bar when clicking
  })
  observeEvent(input$T2D_PREvsREC_selected,{
    gene$gene <- input$T2D_PREvsREC_selected[length(input$T2D_PREvsREC_selected)]
    updateSelectizeInput(inputId = "gene",
                         choices = datRes$Symbol |> unique() |> sort(decreasing = F),
                         server = TRUE,
                         selected = NULL)#clears the search bar when clicking
  })
  
  observeEvent(input$gene,{
    gene$gene <- input$gene
  })
  
  #############
  #boxplot
  output$boxPlot <- renderPlot({
    boxPlot(highlight = gene$gene)}
  )
}

####################################################
#Run

shinyApp(ui, server)




