#' @export start_GSEA
#' @import shiny
#' @import tidyverse
#' @import vroom

#library(shiny)
#library(tidyverse)
#library(vroom)

ui.GSEA <- fluidPage(
  titlePanel("GSEA analysis", windowTitle = "Ss"),
  sidebarLayout(
    sidebarPanel(width = 3,
                 h4("Import Data"),
                 fileInput("upload", NULL, accept = c("tsv", "csv", "txt")),

                 h4("Analysis Select"),
                 selectInput("enrich.method", "Method",choices = c("gseGO", "gseKEGG", "Reactome"), selected = "gseKEGG"),

                 h4("GO Ont"),
                 selectInput("ont", "ont select",choices = c("BP", "MF", "CC", "ALL"), selected = "MF"),

                 h4("Download result"),
                 sliderInput("height", "height", min = 100, max = 2000, value = 400, step = 1),
                 sliderInput("width", "width", min = 100, max = 2000, value = 700, step = 1),
                 downloadButton('plotDown',label="Download Plot"),
                 downloadButton("download.df", "Download Detail"),
                 downloadButton("download.FCdf", "Download FC"),
                 tableOutput("summary")),

    mainPanel(
      tabsetPanel(
        tabPanel("Plot",
                 sidebarLayout(
                   sidebarPanel(
                     h4("Plot"),
                     selectInput("pathwayuse", "Pathway", choices = NULL),
                     #h4("Gene Select"),
                     #numericInput("p", "P-value", value = 0.05, min = 0, max = 1, step = 0.001),
                     #numericInput("fc", "FoldChange", value = 0, min = 0, max = 99999, step = 1),
                     h4("Thresold"),
                     numericInput("p.kegg", "pvalueCutoff", value = 0.05, min = 0, max = 1, step = 0.001),
                     numericInput("q.kegg", "qvalueCutoff", value = 0.2, min = 0, max = 1, step = 0.001),
                     selectInput("p.adjmethod", "pAdjustMethod", choices = c( "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), selected = "BH"),
                     numericInput("mingene", "minGSSize", value = 10, min = 1, max = 99999, step = 1),
                     numericInput("maxgene", "maxGSSize", value = 500, min = 1, max = 99999, step = 1)
                   ),
                   mainPanel(plotOutput("kegg.p"))
                 )
        ),
        tabPanel("Table", dataTableOutput("kegg.df")),
        tabPanel("fcTable", dataTableOutput("inputTable"))
      )
    )

  ))


server.GSEA <- function(input, output, session) {

  data <- reactive({
    req(input$upload)
    df <- vroom(input$upload$datapath)
    # filter
    df.use <- df %>% mutate(
      stat =  case_when(
        (log2FoldChange > 0) ~ "up",
        (log2FoldChange < 0) ~ "down",
        TRUE ~ "None"
      )
    )

    if (input$enrich.method == "gseKEGG") {
      enrichobj <- kegg.gsea(
        gene.table = df,
        pvalueCutoff = input$p.kegg,
        pAdjustMethod = input$p.adjmethod,
        minGSSize = input$mingene,
        maxGSSize = input$maxgene
      )
    }else if (input$enrich.method == "gseGO") {
      enrichobj <- go.gsea(
        gene.table = df,
        pvalueCutoff = input$p.kegg,
        pAdjustMethod= input$p.adjmethod,
        minGSSize = input$mingene,
        maxGSSize= input$maxgene,
        ont= input$ont
      )
    }


    return(list(enrichobj, df.use, df.use))
  })

  summ <- reactive({
    data()[[3]] %>%
      group_by(stat) %>%
      summarise(ngenes = n())
  })

  observeEvent(data()[[1]], {
    updateSelectInput(inputId = "pathwayuse", choices = data()[[1]]@result[["ID"]])
  })

  output$download.df <- downloadHandler(
    filename = function() {
      paste0(input$enrich.method, 'detail','.tsv')
    },
    content = function(file) {
      vroom::vroom_write(data()[[1]]@result, file)
    })

  output$download.FCdf <- downloadHandler(
    filename = function() {
      paste0(input$enrich.method, 'Foldchange','.tsv')
    },
    content = function(file) {
      vroom::vroom_write(data()[[2]], file)
    })

  output$plotDown <- downloadHandler(
    filename = function(){
      paste0(input$enrich.method, '.','png')
    },
    content = function(file){

      Cairo::CairoPDF(
        file = file,
        width = input$width/100,
        height = input$height/100)
      print(plot.gsea(
        data()[[1]], pathway_name = input$pathwayuse))
      dev.off()

      #ggsave(filename = file, plot = plot.gsea(
      #  data()[[1]], pathway_name = input$pathwayuse), width = input$width/100, height = input$height/100)
    }
  )

  output$kegg.p <- renderPlot(plot.gsea(
    data()[[1]], pathway_name = input$pathwayuse),width = function()input$width, height = function()input$height)

  output$kegg.df <- renderDataTable(data()[[1]]@result, options = list(pageLength = 10))

  output$inputTable <- renderDataTable(data()[[2]])

  output$summary <- renderTable(summ())

}

start_GSEA <- function(ui = ui.GSEA, server = server.GSEA) {
  shiny::shinyApp(ui, server)
}



