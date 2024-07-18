library(shiny)
library(shinyFiles)
library(Seurat)
library(ggplot2)



# Increase the maximum file upload size to 10GB
options(shiny.maxRequestSize = 10 * 1024^3)

ui <- fluidPage(
    titlePanel("Seurat Object Viewer"),
    sidebarLayout(
        sidebarPanel(
            fileInput("seurat_file", "Upload Seurat Object", accept = c(".rds")),
            textInput("gene", "Enter Gene Name", value = ""),
            numericInput("plot_width", "Plot Width (px)", value = 600, min = 100),
            numericInput("plot_height", "Plot Height (px)", value = 400, min = 100)
        ),
        mainPanel(
            plotOutput("dimplot", width = "100%"),
            plotOutput("featureplot", width = "100%"),
            plotOutput("vlnplot", width = "100%")
        )
    )
)

server <- function(input, output, session) {
    seurat_obj <- reactiveVal(NULL)
    
    observeEvent(input$seurat_file, {
        req(input$seurat_file)
        seurat_obj(readRDS(input$seurat_file$datapath))
    })
    
    output$dimplot <- renderPlot({
        req(seurat_obj())
        DimPlot(seurat_obj())
    }, width = reactive({input$plot_width}), height = reactive({input$plot_height}))
    
    output$featureplot <- renderPlot({
        req(seurat_obj())
        FeaturePlot(seurat_obj(), features = input$gene, label = T, repel = T)
    }, width = reactive({input$plot_width}), height = reactive({input$plot_height}))
    
    output$vlnplot <- renderPlot({
        req(seurat_obj())
        VlnPlot(seurat_obj(), features = input$gene)
    }, width = reactive({input$plot_width}), height = reactive({input$plot_height}))
}

shinyApp(ui, server)
