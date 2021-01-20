server <- function(input, output, session) {
    plotAllInteractiveServer("all_samples")
    plotTissueInteractiveServer("tissues")
}
