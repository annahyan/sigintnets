library(chrwiseSignatures)
library(ggplot2)
library(ggraph)
library(igraph)
library(visNetwork)

library(shinythemes)
library(shinycssloaders)
library(ggrepel)

load("/home/ahakobyan/ClusterProjects/GitLab/signaturesNetwork/tcga.all.sigs.RData")
load("/home/ahakobyan/ClusterProjects/GitLab/signaturesNetwork/TCGA.tissue.all.sigs.list.RData")

signature_etiologies = read.delim("signature_annotations.tsv", h = T)
table(signature_etiologies$Class)

ui <- fluidPage(
    navbarPage("SigIntNets", theme = shinytheme("cosmo"),
               ### All samples
               tabPanel("All samples", fluid = TRUE, # tags$style(button_color_css),
                        
                        sidebarLayout (
                            sidebarPanel (
                                ## Application title
                                titlePanel("Signature interactions"),
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
                                    size = NULL )
                            ),
                            ##
                            ## Main Panel
                            ## 
                            mainPanel(
                                visNetworkOutput("intNetworkPlot")
                            )
                                
                        )
                        ),
               ### TCGA samples
               tabPanel("TCGA samples", fluid = TRUE,
                        sidebarLayout(
                            sidebarPanel(
                                titlePanel("Signature interactions in TCGA tissues"),

                                selectInput(
                                    inputId = "tissue",
                                    label = "Choose the tissue",
                                    choices = names(TCGA.tissue.all.sigs.list),
                                    selected = NULL,
                                    multiple = FALSE,
                                    selectize = FALSE,
                                    width = NULL,
                                    size = NULL ),
                                
                                selectInput(
                                    inputId = "metric_type",
                                    label = "Choose metric type",
                                    choices = names(TCGA.tissue.all.sigs.list[[1]]),
                                    selected = NULL,
                                    multiple = FALSE,
                                    selectize = TRUE,
                                    width = NULL,
                                    size = NULL) 
                            ),
                            mainPanel(
                                visNetworkOutput("intNetworkPlot")
                            )
                            
                        )
                        )
               
               )
)



# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$intNetworkPlot <- renderVisNetwork({
        # minimal example
        cor.matrix = tcga.all.sigs[[input$metric_type]]
        sig_c = nrow(cor.matrix)

        if(input$metric_type == "cooc") {
            threshold = 2
            cor.matrix[cor.matrix > 6] = 6

            cor.matrix[ abs(cor.matrix) < 2] = 0

        } else {
            threshold = 0.2
        }
 
        graph.input = cor.matrix
        graph.input[ abs(graph.input) < threshold ] = 0
        
        diag(graph.input) = 0
        
        graph.input[abs(graph.input) >= threshold ] = 1
        
        ## cor.matrix[ abs(cor.matrix) < 0.2 ]
        
        sig.adjacency = graph.adjacency(graph.input)
        
        graph <- graph_from_data_frame(get.edgelist(sig.adjacency) )
        
        data = toVisNetworkData(graph)
        
        nodes = data$nodes
        nodes$group = signature_etiologies[match( 
            nodes$id, signature_etiologies$Signature), "Class"]
        nodes$font.size = 25
        
        edges = data$edges
        
        edges$intensity = sapply( 1:length(E(graph)),
                                  function(x) {
                                      gends = ends(graph, x)
                                      cor.matrix[gends[1], gends[2] ]
                                  } )
        
        edges$width = abs(edges$intensity)
        
        
            col.edges = c( rgb(0, 140, 160, maxColorValue = 255),
                           rgb(210, 50, 60, maxColorValue = 255) )
            
            edges$color =  col.edges[ (edges$intensity > 0) + 1 ] 
       #  }
        
        visNetwork(nodes, edges) %>% visIgraphLayout(layout = "layout_nicely") %>%
            visGroups(groupname = "Ageing", color = "#97c2fc") %>%
            visGroups(groupname ="APOBEC/AID", color = "#ffff00") %>% 
            visGroups(groupname = "BER/HR", color = "#fb7e81") %>%
            visGroups(groupname = "Environment", color = "#7be141") %>%
            visGroups(groupname = "Unknown", color = "#ffc0cb") %>%
            visGroups(groupname = "MMR", color = "#eb7df4") %>%
            visGroups(groupname = "UV", color = "#ffa807") %>%
            visGroups(groupname = "Polymerase", color = "#ad85e4") %>%
            visGroups(groupname = "Technical", color = "gray") %>%
            visGroups(groupname = "Chemo", color = "#ffc0cb") %>%
            visLegend()
        
    }) 
}

shinyApp(ui, server)
