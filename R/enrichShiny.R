#' @export start_KEGG_GO
#' @import shiny
#' @import vroom

ui.keggo <- fluidPage(
  titlePanel("Enrich analysis", windowTitle = "Ss"),
  sidebarLayout(
    sidebarPanel(width = 3,
                 h4("Import Data"),
                 fileInput("upload", NULL, accept = c("tsv", "csv", "txt")),

                 h4("Analysis Select"),
                 selectInput("enrich.method", "Method",choices = c("GO", "KEGG"), selected = "KEGG"),

                 h4("GO Ont"),
                 selectInput("ont", "ont select",choices = c("BP", "MF", "CC", "ALL"), selected = "MF"),

                 h4("Download result"),
                 sliderInput("height", "height", min = 100, max = 2000, value = 400, step = 1),
                 sliderInput("width", "width", min = 100, max = 2000, value = 400, step = 1),
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
                     selectInput("p.type", "type", choices = c("barplot", "dotplot"), selected = "dotplot"),
                     numericInput("show", "showCategory", value = 10, min = 1, max = 99999, step = 1),
                     h4("Gene Select"),
                     numericInput("p", "P-value", value = 0.05, min = 0, max = 1, step = 0.001),
                     numericInput("fc", "FoldChange", value = 0, min = 0, max = 99999, step = 1),
                     selectInput("st", "stat",choices = c("up", "down"), selected = "up"),
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


server.keggo <- function(input, output, session) {
  data <- reactive({
    req(input$upload)
    df <- vroom(input$upload$datapath)
    # filter
    df <- df %>% mutate(
      stat =  case_when(
        (padj <= input$p & log2FoldChange >= input$fc) ~ "up",
        (padj <= input$p & log2FoldChange <= -input$fc) ~ "down",
        TRUE ~ "None"
      )
    )

    df.use <- df %>% filter(stat == input$st)

    enrichobj <- enrich.Select(
      df = df.use,
      enrich.method = input$enrich.method,
      pvalueCutoff = input$p.kegg,
      pAdjustMethod = input$p.adjmethod,
      minGSSize = input$mingene,
      maxGSSize = input$maxgene,
      qvalueCutoff = input$q.kegg,
      ont = input$ont
    )

    return(list(enrichobj, df.use, df))
  })

  summ <- reactive({
    data()[[3]] %>%
      group_by(stat) %>%
      summarise(ngenes = n())
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
      ggsave(filename = file,plot = enrich.p(
        data()[[1]],
        show = input$show, type = input$p.type), width = input$width/100, height = input$height/100)
    }
  )

  output$kegg.p <- renderPlot(enrich.p(
    data()[[1]],
    show = input$show, type = input$p.type),width = function()input$width, height = function()input$height)

  output$kegg.df <- renderDataTable(enrich.df(data()[[1]]) ,options = list(pageLength = 10))
  output$inputTable <- renderDataTable(add.id(data()[[3]]))
  output$summary <- renderTable(summ())

}

start_KEGG_GO <- function(ui = ui.keggo, server = server.keggo) {
  shinyApp(ui, server)}
