#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  # # #
  ### Launcher Panel
  # # # 
  
  # Display template to save input TSV for the compound
  tsv_template = as.data.frame(readr::read_tsv("../input/TEMPLATE_start_mlc.tsv"))
  output$mlcTSV = DT::renderDataTable(tsv_template, selection = 'none', rownames = FALSE, editable=T, server=F)
  proxy = dataTableProxy("mlcTSV")
  
  # Update the table when a cell is edited
  observeEvent(input$mlcTSV_cell_edit, {
    info = input$mlcTSV_cell_edit
    i = info$row
    j = info$col + 1
    v = info$value
    tsv_template[i, j] <<- DT:::coerceValue(v, tsv_template[i, j])
  })
  
  # Add new row in table 
  observeEvent(input$add, {
    proxy %>% addRow(as.data.frame(readr::read_tsv("../input/TEMPLATE_start_mlc.tsv")))
  })
  
  # Save the input TSV in the input directory
  observeEvent(input$save, {
    
    readr::write_tsv(tsv_template, paste0("../input/", input$tsvname, ".tsv"))
    
    shinyalert(
      title = "MetApp",
      text = paste0("You saved ", input$tsvname, ".tsv in the input directory."),
      html = FALSE,
      type = "success",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
      confirmButtonCol = "#AEDEF4",
      timer = 0,
      imageUrl = "",
      animation = TRUE
    )
  })
  
  # Rscript path
  shinyFileChoose(input, "Rscript", roots = c(home = normalizePath("/")))
  Rscript = reactive(input$Rscript)
  observeEvent(input$Rscript, { 
    validate(need( length(unlist(Rscript()[1])) > 1, "" ))
    output$select_Rscript = renderText( paste( unlist(Rscript()[1]), collapse = "/" ) )
  })
  
  # INPUT file button
  shinyFileChoose(input, "input", roots = c(home = normalizePath("~")))
  infile = reactive(input$input)
  observeEvent(input$input, { 
    validate(need( length(unlist(infile()[1])) > 1, "" ))
    output$select_input = renderText( paste( unlist(infile()[1]), collapse = "/" ) )
  })
  
  
  # OUT PATH directory button
  shinyDirChoose(input, "outdir", roots = c(home = normalizePath("~")))
  outdir = reactive(input$outdir)
  observeEvent(input$outdir, { 
    validate(need( length(unlist(outdir()[1])) > 1, "" ))
    output$select_outdir = renderText( paste(unlist(outdir()[1]), collapse = "/") )
  })
  
  # OUT PATH directory button
  shinyFileChoose(input, "config", roots = c(home = normalizePath("~")))
  config = reactive(input$config)
  observeEvent(input$config, { 
    validate(need( length(unlist(config()[1])) > 1, "" ))
    output$select_config = renderText( paste(unlist(config()[1]), collapse = "/") )
  })
  
  observe({
    rscript_path = ifelse(length(unlist(Rscript()[1])) > 1, 
                          paste( unlist(Rscript()[1]), collapse = "/" ),
                          "Rscript") 
    input_path = paste( unlist(infile()[1]), collapse = "/" )
    out_path = paste0(paste(unlist(outdir()[1]), collapse = "/"), "/", input$output)
    updateTextInput(session, "command", 
                    value= paste0(rscript_path, " MetIDfyR.R -i ~", input_path, 
                                  " -o ~", out_path,
                                  ifelse(length(unlist(config()[1])) > 1,
                                         paste0(" -c ~", paste(unlist(config()[1]), collapse = "/")),
                                         "") ) )
    
  })
  
  # Add clipboard buttons
  output$clip <- renderUI({
    rclipButton("clipbtn", "rclipButton Copy", input$command, icon("clipboard"))
  })
  
  # # #
  ### Configuration Panel
  # # #
  
  configFields = c("lib_perso", "bool_phase_2", "cores", "nb_transformation", "min_peak_intensity", "mz_ppm", "rt_windows", 
                   "nb_scan", "wdw_mz_ms2", "minimum_mz")
  
  fieldsAll = c("bool_phase_2", "cores", "nb_transformation", "min_peak_intensity", "mz_ppm", 
                "rt_windows", "nb_scan", "wdw_mz_ms2", "minimum_mz")
  fieldsSelect = c("lib_perso", "ms2_reference")
  
  observe({
    
    # MS2 REFERENCE file button
    shinyFileChoose(input, "ms2_reference", roots = c(home = normalizePath("~")))
    ms2_tsv = reactive(input$ms2_reference)
    observeEvent(input$ms2_reference, { 
      validate(need( length(unlist(ms2_tsv()[1])) > 1, "" ))
      output$select_ms2_ref = renderText( paste( unlist(ms2_tsv()[1]), collapse = "/" ) )
    })
    
    # R LIBRARY directory button
    shinyDirChoose(input, "lib_perso", roots = c(home = normalizePath("~")))
    libperso = reactive(input$lib_perso)
    observeEvent(input$lib_perso, { 
      validate(need( length(unlist(libperso()[1])) > 1, "" ))
      output$select_lib_perso = renderText( paste(unlist(libperso()[1]), collapse = "/") )
    })
    
    # Check mandatory fields are filled, enable submit if TRUE
    mandatoryFilled <- input$config_name != "" 
    shinyjs::toggleState(id = "submit", condition = mandatoryFilled)
    
    # Store data from the form
    formData = reactive({
      config = sapply(fieldsSelect, function(x) paste(unlist(input[[x]][1]), collapse = "/"))  # Generate path
      config = c(config, sapply(fieldsAll, function(x) input[[x]]))  # Extract parameters values
      t(config)
    }) 
    
    # When submit is pressed
    observeEvent(input$submit, {
      req(mandatoryFilled)
      
      tmp_form = tibble::as_tibble(formData())
      
      file = paste0("../input/", input$config_name, ".R")
      
      myForm = as.data.frame(cbind(c(configFields, "list_transfo"), 
                                   c(tmp_form[configFields], "\"input/list_transfo.tsv\"")))
      
      write(x="#Config file for MetIDfyR", file)
      for(i in 1:nrow(myForm)){
        if( myForm[i,1] == "lib_perso" & length(unlist(libperso()[1])) > 1 ){
          content = paste0("\"", myForm[i,2], "\"")
        }else content = myForm[i,2]
        
        write(paste0(myForm[i,1], " = ", ifelse(content=="", NA, content)), file , append=TRUE)
      }
      
      updateTextInput(session, "config_name", value = "", placeholder = "config.R")
      
      shinyalert(
        title = "MetApp",
        text = "The file is saved in the input directory.",
        html = FALSE,
        type = "success",
        showConfirmButton = TRUE,
        confirmButtonText = "OK",
        confirmButtonCol = "#AEDEF4",
        timer = 0,
        imageUrl = "",
        animation = TRUE
      )
      
    })
    
  }) 
  
  observeEvent(input$resetbut, {
    shinyjs::reset("form")
    output$select_input = renderText({})
    output$select_outdir = renderText({})
    output$select_lib_perso = renderText({})
    output$select_ms2_ref = renderText({})
  })
  
  
  
### Visualization panel
  
  shinyDirChoose(input, "dir", roots = c(home = normalizePath("~")))
  dir = reactive(input$dir)
  path = reactive({ 
    home = normalizePath("~")
    if(length(unlist(dir()[1])) > 1){
      
      file.path(home, paste(unlist(dir()$path[-1]), collapse = "/"))
    }else{
      file.path()
    }
  })
  
  #Update the metabolites list when the output folder change
  observeEvent(path(), {
    Product = reactive({
      validate(need( length(path()) > 0, "Choose a directory"))
      total_table = readr::read_tsv(list.files(path(), pattern = "out.*\\.tsv$", full.names = T))
      num_col = names(Filter(is.numeric, total_table))
      total_table[num_col] = as_tibble(apply(total_table[num_col], 2, function(col) round(col, 3)))
      total_table$intensity = scales::scientific(total_table$intensity)
      total_table = plyr::rename(total_table, c("name"="Name", "formula"="Formula", "polarity"="Polarity", "adduct"="Adduct", "mz"="m/z", 
                                           "diff"="Change", "rt"="RT", "intensity"="Intensity", "abscore"="iAScore", 
                                           "dotp_ms2"="MS2 dotp", "common_ms2_peak"="MS2p", "mono_ppm"="Delta m/z (ppm)",
                                           "index_peak"="Index", "nb_transfo"="Nb Transfo", "score"="Score"))
      total_table %>%
        arrange(desc(Score))
    })
    
    if(length(path()) > 0){
      output$header = renderText({ paste0("Currently open: ", unique(Product()$Name) ) })
    }else{
      output$header = renderText({"Select a directory which contains the output figures and table for a parent drug in one mzML:"})
    }
    
    
    
    observe({ shinyjs::toggleState(id = "report", condition = length(path()) & length(input$table_rows_selected)) 
      
      output$report = downloadHandler(
        filename = "report.pdf",
        
        content = function(file){
          rmarkdown::render("report_template.Rmd", output_file = file , quiet = T)
        }
      )
    })
    
    updateCheckboxGroupInput(session, "columns", choices = names(Product()),
                             selected =c("Formula", "Polarity", "Change", "m/z", "RT", "MS2p", "iAScore",
                                         "MS2 dotp", "Delta m/z (ppm)", "Intensity", "Score"))
    #Print table
    output$table <- DT::renderDataTable({
      datatable(Product()[, input$columns], rownames = F, selection = "multiple", 
                options = list(scrollX = TRUE, searchHighlight = T))
    })
    
    #Observe click event on table rows
    observeEvent(input$table_cell_clicked, {
      if(length(input$table_cell_clicked) > 0){
        data = Product()[input$table_cell_clicked$row, ]
        
        form = data$Formula
        peak = data$Index
        find_figure = list.files(path(), pattern = paste0(form, "_", peak), recursive = T, full.names = T)
        
        polar = ifelse(data$Polarity=="Positive mode", "POS", "NEG")
        updateTabsetPanel(session, "polarity", selected = polar)
        
        updateSelectInput(session, input$polarity, selected = form )
        observeEvent(input[[input$polarity]], {
          
          updateSelectInput(session, paste0(input$polarity, "_peak"), selected = peak)
        })
        
      }
    })
    
    
    ### POS ###
    
    filesname_pos = list.files(path(), recursive = T, pattern = ".*_POS.svg|.*plus.*svg|.*_POS.png|.*plus.*png")
    metabolites = unique(lapply(basename(filesname_pos), function(x) strsplit(x, '_')[[1]][1]))
    
    updateSelectInput(session, "POS", choices = metabolites )
    
    #Update the chromatogram peaks list
    observeEvent(input$POS, {
      
      filesname_pos = list.files(path(), recursive = T, pattern = paste0("^", input$POS, "_.*_POS|^", input$POS,"_.*plus"))
      peaks = unique(lapply(basename(filesname_pos), function(x) strsplit(x, '_')[[1]][2]))
      
      updateSelectInput(session, "POS_peak", choices = peaks )
      
    })
    
    #Update the POS page image
    observeEvent(input$POS_peak, {
      
      output$figure_POS = renderImage({
        validate(need( length(path()) > 0, "Select directory"))
        filename = grep(paste0(input$POS, "_", input$POS_peak, "_"), filesname_pos, value = T)
        
        list(src = file.path(path(), filename),
             alt = "Missing figure in pos",
             width = 800,
             height = 800)
      }, deleteFile = F)
      
      
      output$table_POS <- renderTable({
        tabsname_pos = list.files(path(), recursive = T, 
                                  pattern = paste0(input$POS, "_", input$POS_peak, ".*_POS.tsv|",
                                                   input$POS,"_", input$POS_peak, "_.*plus.*tsv"))
        validate(need( length(path()) > 0, "Select a directory"),
                 need( length(tabsname_pos) == 1, "No MS2 table to show"))
        table=readr::read_tsv(file.path(path(), tabsname_pos))
        table
      })
      
    })
    
    
    ### NEG ###
    
    filesname_neg = list.files(path(), recursive = T, pattern = ".*_NEG.png|.*minus.*png|.*_NEG.svg|.*minus.*svg")
    metabolites = unique(lapply(basename(filesname_neg), function(x) strsplit(x, '_')[[1]][1]))
    
    updateSelectInput(session, "NEG", choices = metabolites )
    
    
    #Update the chromatogram peaks list
    observeEvent(input$NEG, {
      
      filesname_neg = list.files(path(), recursive = T, pattern = paste0("^", input$NEG,".*_NEG|^", input$NEG, "_.*minus"))
      peaks = unique(lapply(basename(filesname_neg), function(x) strsplit(x, '_')[[1]][2]))
      
      updateSelectInput(session, "NEG_peak", choices = peaks )
      
    })
    
    #Update the NEG page image
    observeEvent(input$NEG_peak, {
      
      output$figure_NEG = renderImage({
        validate(need( length(path()) > 0, ""))
        filename = grep(paste0(input$NEG, "_", input$NEG_peak, "_"), filesname_neg, value=T)
        
        list(src = file.path(path(), filename),
             alt = "Missing figure in NEG",
             width = 800,
             height = 800)
      }, deleteFile = F)
      
      
      output$table_NEG <- renderTable({
        tabsname_neg = list.files(path(), recursive = T, 
                                  pattern = paste0(input$NEG, "_", input$NEG_peak, ".*_NEG.tsv|",
                                                   input$NEG,"_", input$NEG_peak, "_.*minus.*tsv"))
        validate(need( length(path()) > 0, "Select a directory"),
                 need( length(tabsname_neg) == 1, "No MS2 table to show"))
        table=readr::read_tsv(file.path(path(), tabsname_neg))
        table
      })
      
      
    })
    
  })
  
})
