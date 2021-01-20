library(chrwiseSignatures)
library(ggplot2)
library(ggraph)
library(igraph)
library(visNetwork)

library(shinythemes)
library(shinycssloaders)
library(ggrepel)
library(shiny)

load("/home/ahakobyan/ClusterProjects/GitLab/signaturesNetwork/tcga.all.sigs.RData")
load("/home/ahakobyan/ClusterProjects/GitLab/signaturesNetwork/TCGA.tissue.all.sigs.list.RData")

signature_etiologies = read.delim("signature_annotations.tsv", h = T)
table(signature_etiologies$Class)

plotAllInteractiveInput = function(id, paneltitle = "Signature Interactions") {

    ns <- NS(id)
    tagList(
            selectInput(
                inputId = ns("dataset"),
                label = "Choose the dataset",
                choices = list("TCGA", "PCAWG"),
                selected = NULL,
                multiple = FALSE,
                selectize = TRUE,
                width = NULL,
                size = NULL),
            selectInput(
                inputId = ns("metric_type"),
                label = "Choose metric type",
                choices = names(tcga.all.sigs),
                selected = NULL,
                multiple = FALSE,
                selectize = TRUE,
                width = NULL,
                size = NULL ),
            visNetworkOutput(ns("intNetworkPlot")))
}

getNodesEdges = function(cor.matrix, metric) {
    
    sig_c = nrow(cor.matrix)
    
    if(metric == "cooc") {
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
    return(list(nodes = nodes, edges = edges))
}

plotVisNetwork = function(innodes, inedges) {
    visNetwork(innodes, inedges) %>%
        visIgraphLayout(layout = "layout_nicely") %>%
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
}


plotAllInteractiveServer = function(id) {
    moduleServer(
        id,
        function(input, output, session) {
            cor.matrix = reactive({tcga.all.sigs[[input$metric_type]]})
            
            outlist = reactive({ 
                getNodesEdges(cor.matrix(), input$metric_type)
            })
            
            
            # output$out = renderText(names(outlist()))
            output$intNetworkPlot = renderVisNetwork({
                plotVisNetwork(outlist()$nodes, outlist()$edges)
            })
            # output$out = renderText({"blah"})
            # return(list(nodeout = nodes, edgeout = edges))
        }
    )
}

plotTissueInteractiveInput = function(id, paneltitle = "Signature Interactions") {

    ns <- NS(id)
    tagList (
            selectInput(
                inputId = ns("tissue"),
                label = "Choose the tissue",
                choices = names(TCGA.tissue.all.sigs.list),
                selected = NULL,
                multiple = FALSE,
                selectize = FALSE,
                width = NULL,
                size = NULL ),
            
            selectInput(
                inputId = ns("metric_type"),
                label = "Choose metric type",
                choices = names(TCGA.tissue.all.sigs.list[[1]]),
                selected = NULL,
                multiple = FALSE,
                selectize = TRUE,
                width = NULL,
                size = NULL),
            # verbatimTextOutput(ns("out"))
            visNetworkOutput(ns("intNetworkPlot"))   
        )
}

plotTissueInteractiveServer = function(id) {
    moduleServer(
        id,
        function(input, output, session) {

            cor.matrix = reactive({TCGA.tissue.all.sigs.list[[input$tissue]][[input$metric_type]]})

            outlist = reactive({
                getNodesEdges(cor.matrix(), input$metric_type)
            })

            
            output$intNetworkPlot = renderVisNetwork({
                plotVisNetwork(outlist()$nodes, outlist()$edges)
            })
            # output$out = renderText({"blah"})
            # return(list(nodeout = nodes, edgeout = edges))
        }
    )
}
