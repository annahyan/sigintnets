ui <- fluidPage(
    navbarPage("SigIntNets", theme = shinytheme("cosmo"),
               ### All samples
               tabPanel("All samples", fluid = TRUE, # tags$style(button_color_css),
                        plotAllInteractiveInput("all_samples")
                        ),
### TCGA samples
               tabPanel("TCGA samples", fluid = TRUE,
                                plotTissueInteractiveInput("tissues")
                            )
               )
)
