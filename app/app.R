load("dat.Rdata")


library(shinydashboard)
library(shiny)
library(data.table)
library(shinycssloaders)

source("plotFunctions.R")

####################################################
#UI
ui <- function(){
  
  header <- dashboardHeader(
    title = "Gene expresssion changes in adipose tissue after exercise", titleWidth = 800
  )
  
  side <- dashboardSidebar(disable = T)
  
  body <- dashboardBody({
                          fixedPage(
                            box("Search for your favorite gene",
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
              Symbol = datMD$Symbol,
              highlight = gene$gene) |> girafize()
  })
  #############
  #NFT pre rec
  output$NGT_PREvsREC <- renderGirafe({
    exp2logFC(aveExp = datMD$NGT_PREvsREC, 
              logFC = datMD$NGT_PREvsREC.logFC,
              sig = datMD$NGT_PREvsREC.sig,
              Symbol = datMD$Symbol,
              highlight = gene$gene) |> girafize()
  })
  #############
  #T2D pre post
  output$T2D_PREvsPOST <- renderGirafe({
    exp2logFC(aveExp = datMD$T2D_PREvsPOST,
              logFC = datMD$T2D_PREvsPOST.logFC,
              sig = datMD$T2D_PREvsPOST.sig,
              Symbol = datMD$Symbol,
              highlight = gene$gene) |> girafize()
  })
  #############
  #T2D pre rec
  output$T2D_PREvsREC <- renderGirafe({
    exp2logFC(aveExp = datMD$T2D_PREvsREC,
              logFC = datMD$T2D_PREvsREC.logFC,
              sig = datMD$T2D_PREvsREC.sig,
              Symbol = datMD$Symbol,
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
    boxPlot(highlight = gene$gene,
            datBox = datBox,
            datRes = datRes)}
  )
}

####################################################
#Run

shinyApp(ui, server)




