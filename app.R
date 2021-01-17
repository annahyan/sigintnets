library(chrwiseSignatures)
library(ggplot2)
library(ggraph)

load("/home/ahakobyan/ClusterProjects/GitLab/signaturesNetwork/tcga.all.sigs.RData")



ui <- fluidPage(
    
    # Application title
    titlePanel("Survival analysis"),
    
    # Sidebar with a slider input for number of bins 
    # sidebarLayout(
    #     sidebarPanel(
    #         sliderInput("bins",
    #                     "Number of bins:",
    #                     min = 1,
    #                     max = 50,
    #                     value = 30)
    #     ),
    #     
    #     # Show a plot of the generated distribution
    #     mainPanel(
    #         plotOutput("distPlot")
    #     )
    # )
    selectInput(
        inputId = "dataset",
        label = "Choose the dataset",
        choices = list("TCGA", "PCAWG"),
        selected = NULL,
        multiple = FALSE,
        selectize = TRUE,
        width = NULL,
        size = NULL),
    
    selectInput(
        inputId = "metric_type",
        label = "Choose metric type",
        choices = names(tcga.all.sigs),
        selected = NULL,
        multiple = FALSE,
        selectize = TRUE,
        width = NULL,
        size = NULL
    ),
    mainPanel(
        plotOutput("networkPlot")
    )
)




# Define server logic required to draw a histogram
server <- function(input, output) {
    
    # output$distPlot <- renderPlot({
    #     # generate bins based on input$bins from ui.R
    #     x    <- faithful[, 2]
    #     bins <- seq(min(x), max(x), length.out = input$bins + 1)
    #     
    #     # draw the histogram with the specified number of bins
    #     hist(x, breaks = bins, col = 'darkgray', border = 'white')
    # })
    output$networkPlot = renderPlot({
        
        cor.matrix = tcga.all.sigs[[input$metric_type]]
        sig_c = nrow(cor.matrix)
        
        if(input$metric_type == "cooc") {
            threshold = 1
            cor.matrix[cor.matrix > 6] = 6
            
            cor.matrix[ abs(cor.matrix) < 2] = 0
            
        } else {
            threshold = 0.2
        }
        plot_network(cor.matrix[1:sig_c, 1:sig_c], min_threshold = threshold)
    })
}


shinyApp(ui, server)