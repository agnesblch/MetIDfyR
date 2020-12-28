#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Function for module UI
titlePanelUI <- function(id, name){
    
    ns = NS(id)
    fluidPage(
        fluidRow(column(2, imageOutput(ns("logo"), height = "100%", width = "100%")),
                 column(8, titlePanel(name))), hr()
    )
    
}

# Function for module server logic
titlePanelMod <- function(input, output, session) {
    
    output$logo = renderImage({
        list(src = "www/logo_metidfyr.png",
             alt = "Missing logo image",
             width = 100, height = 100)
    }, deleteFile = F)
    
}
