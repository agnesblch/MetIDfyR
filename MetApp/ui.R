#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#
if(!"pacman" %in% installed.packages()) install.packages("pacman")
pacman::p_load("shiny", "shinythemes", "shinyFiles", "DT", "shinyjs", "parallel", "shinyalert", "dplyr", "rclipboard","rmarkdown", "rsvg")

source("modules/titlePanelMod.R")

labelMandatory <- function(label) {
  tagList(
    label,
    span("*", class = "mandatory_star")
  )
}

mandatoryCSS <- ".mandatory_star { color: red; }"

#### UI ####

navbarPage("MetApp for MetIDfyR",
           theme = shinytheme("yeti"),
           
           
# Launcher for MetIDfyR
           tabPanel("Launcher",
                    
                    titlePanelUI(id = "t1", name = "MetIDfyR launcher"),
                    
                    h4("Required :"),
                      p("-a TSV file (see Input tab)"),
                      p("-an output directory and an output name"),
                      p("-a config file (see TEMPLATE_config.R or build it using Configuration tab)"),
                    
                    br(),
                    fluidRow(column(3, tags$b("Path to Rscript if not in PATH"), br(),
                                    shinyFilesButton("Rscript", "Select Rscript", 
                                                     "Path to Rscript", multiple = F), br(),
                                    textOutput("select_Rscript")),
                            column(3, labelMandatory(tags$b("Path to input file")), br(),
                                    shinyFilesButton("input", "Select input", 
                                                     "Path to input file", multiple = F), br(),
                                    textOutput("select_input")),
                             column(3, labelMandatory(tags$b("Path to save the output from the software")),br(),
                                    shinyDirButton("outdir", "Select output", 
                                                   "Path to save the output from the software"),
                                    textOutput("select_outdir"))
                            ),
                    fluidRow(column(3,
                                    textInput("output", labelMandatory("Name of the output directory"), placeholder = "my_output")),
                             column(3, tags$b("Path to config file"), br(),
                                    shinyFilesButton("config", "Select config file", 
                                                     "Path to config file", multiple = F), br(),
                                    textOutput("select_config"))
                             ), br(),
                    
                    textInput("command", label="", width="200%"),
                    rclipboardSetup(),
                    uiOutput("clip")
                    
                    ),
           
           
# Input file creation
           tabPanel("Input",
                    titlePanelUI(id="t2", "Input file generation"),
                    h4("Fill the following table and save it if you don't have an input TSV."),
                    p("Double click to modify a value. The file will be saved in the input directory."),
                    p("Leave the cell EMPTY if there's nothing to put in."),
                    actionButton("add", "Add row"),
                    DT::dataTableOutput('mlcTSV'),
                    br(),
                    textInput("tsvname", "Name of the input TSV"),
                    actionButton("save", "Save tsv")
                    ),
# Configuration file creation
           tabPanel("Configuration",
                    shinyjs::useShinyjs(),
                    useShinyalert(),
                    shinyjs::inlineCSS(mandatoryCSS),
                    
                    titlePanelUI(id="t3", "Configuration"),
                    h4("Please enter your parameters for MetIDfyR"),
                    h4("Responses will be saved as a R file in the input directory."),
                    fluidRow(id = "form",
                             column(6,
                                    textInput("config_name", "Name of the config file to save", placeholder = "config.R"),
                                    
                                    selectInput("bool_phase_2", "Include phase II transformations", choices = c("FALSE", "TRUE")),
                                    
                                    selectInput("cores", "Number of cores to use", choices = c(1:detectCores()-1)),
                                    
                                    tags$b("Path to personal R library if required"), br(),
                                    shinyDirButton("lib_perso", "Select library", "Path to personal R library if required"),
                                    textOutput("select_lib_perso"),
                                    
                                    selectInput("nb_transformation", "Number of transformations to perform", choices = c(1:6))
                                    
                             ),
                             column(6,
                                    textInput("min_peak_intensity", "Minimum intensity to consider a peak", value="5e5"),
                                    
                                    textInput("mz_ppm", "M/Z tolerance in ppm", value = "10"),
                                    
                                    textInput("rt_windows", "Retention time windows to consider a chromatogram peak", value="5"),
                                    
                                    selectInput("nb_scan", "Minimum number of scan", choices = c(1:10)),
                                    
                                    textInput("wdw_mz_ms2", "Quadrupole selection window width (m/z)", value="5"),
                                    
                                    textInput("minimum_mz", "Minimum MZ", value="200")
                             )
                        
                    ),
                    
                    actionButton("submit", "Submit"),
                    actionButton("resetbut", "Reset")
           ),
           
# Visualization of data
           tabPanel("Visualization",
                    # Application title
                    titlePanelUI(id="t4", "Data visualization"),
                    h4("Here you can visualize the metabolites output produced by MetIDfyR"),
                    
                    fluidRow(column(5, tags$b("Comment area"),
                                    textAreaInput("comment", "", height = "100%", width = "200%", resize = "both"))
                             ), 
                    
                    br(), br(),
                    h4(textOutput("header")),
                    shinyDirButton("dir", "Select directory", "Upload"),
                    downloadButton("report", "Save report"), br(),
                    
                    sidebarLayout(
                      sidebarPanel(
                        checkboxGroupInput("columns","Select Columns"),
                        width = 2
                      ),
                      
                      mainPanel(
                        DT::dataTableOutput('table')
                      )
                    )
                    
                    ,
                    
                    hr(),
                    
                    # Add tabsets 
                    tabsetPanel(type = "tabs", id = "polarity",
                                tabPanel(title = "Positive", value = "POS", 
                                         tags$br(),
                                         fluidRow(column(4, selectInput("POS", "Select a metabolite", choices = NULL)),
                                                  column(4, selectInput("POS_peak", "Select a peak", choices = NULL))),
                                         fluidRow(column(8, imageOutput("figure_POS")),
                                                  column(4, tableOutput("table_POS")))
                                ),
                                
                                tabPanel(title = "Negative", value = "NEG",
                                         tags$br(),
                                         fluidRow(column(4, selectInput("NEG", "Select a metabolite", choices = NULL)),
                                                  column(4, selectInput("NEG_peak", "Select a peak", choices = NULL))),
                                         fluidRow(column(8, imageOutput("figure_NEG")),
                                                  column(4, tableOutput("table_NEG")))
                                )
                    )
           )
           
           
           # tabPanel("Contact",
           #          h3("About MetIDfyR"),
           #          h3("GitHub link")
           # )
           
           
)
