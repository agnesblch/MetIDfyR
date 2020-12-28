#### DEPENDENCIES ####

if(!"pacman" %in% installed.packages()) install.packages("pacman")
pacman::p_load("BiocManager", "optparse")

#### PARAMETERS ####

option_list = list(
  make_option(c("-i", "--input"), type="character",
              help="input informations (see TEMPLATE_start_mlc.tsv)", metavar="character"),
  make_option(c("-o", "--output"), type="character",
              help="output directory", metavar="character"),
  make_option(c("-c", "--config"), type="character", default="input/config.R",
              help="config path", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Test option
if (is.null(opt$input )){
  print_help(opt_parser)
  stop("Input tsv required.n", call.=FALSE)
}
if (is.null(opt$output )){
  print_help(opt_parser)
  stop("Output directory required.n", call.=FALSE)
}

source(opt$config) #load config informations
source("util.R") #load functions file
source("plot_functions.R") #load plot functions file


.libPaths( c( .libPaths(), lib_perso))
suppressPackageStartupMessages(pacman::p_load("dplyr", "plyr", "doParallel", "stringr", "MSnbase", "Rdisop", 
                                              "readr", "tibble", "arrangements", "snow", "svglite"))

options(digits=10) #allow more digits
E_MASS = c(minus = 5.489E-4, plus = -5.489E-4) 

dir.create(opt$output, showWarnings = FALSE)



#### Initialisation : Data loading and preparation ####



cat("### Initialisation ###\n")
data_tsv = read_tsv(opt$input, col_types = cols()) #table containing the molecules
data_tsv = as_tibble(cbind(data_tsv, fill_table(data_tsv)))

transfo_init = getShift(list_transfo) #table containing the potential transformation
transfo_init$phase = as.factor(transfo_init$phase)
transfo_init$type = getType(transfo_init)

# print config
cat("### CONFIGURATION ###\n")
writeLines(readLines(opt$config))
cat("\n")

# Loop on each row of the input table
for(row in 1:nrow(data_tsv)){

  start=Sys.time()

  ##Launch cluster with n-1 cores (n = number of available cores)
  cores = ifelse(is.na(cores) | cores == "", detectCores()-1, cores)
  cl = makeCluster(cores)
  clusterApply(cl, 1:cores, opt, fun= function(x,opt){
    source(opt$config)
    source("util.R")
    source("plot_functions.R")
    library(pacman)
    p_load("tibble", "MSnbase", "plyr", "ggplot2", "ggpubr", "Rdisop", "dplyr", "ggrepel",
           "stringr", "htmlwidgets", "readr")
  })

  transfo = transfo_init ; transfo$possible = TRUE
  data = data_tsv[row, ] ; do = c()
  do$plus = !is.na(data$adduct_plus) ; do$minus = !is.na(data$adduct_minus)
  if(!do$plus) data$adduct_plus = FALSE ; if(!do$minus) data$adduct_minus = FALSE
  if(is.na(data$ms2_reference_tsv)) data$ms2_reference_tsv = FALSE
  

  #### Possible transformation ####

  cat(paste0("### Select combination for ", data$name, " ###\n"))

  #check if phase 2 is performed
  if(!bool_phase_2) transfo$possible[which(transfo$phase==2)] = FALSE

  #check the presence of specifics atoms in the parent drug formula
  for(atom in c("N", "F", "Cl", "Br")){
    # if the atom is not present in the parent drug formula, then it can't be remove
    if(!atom %in% colnames(data)){
      row_atom=grep(atom, transfo$remove)
      transfo$possible[row_atom] = FALSE
    }
  }

  #Get all combinations with replacement
  #Check the feasibility of each combination
  combin_transfo = arrangements::combinations(which(transfo$possible), nb_transformation, replace=TRUE)

  #Build groups by cores for foreach loop
  group = rep(1:cores, each=floor(nrow(combin_transfo)/cores))
  group = append(group, rep(1, nrow(combin_transfo)%%cores))



  #### Transformation ####



  cat(paste0("### Do transformation for ", data$name, " ###\n"))
  clusterExport(cl, ls())

  info_all_combi = parLapply (cl, unique(group), function(current_grp){

    if(ncol(combin_transfo) >= 1) {
      grp_index = which(group == current_grp)
      if(ncol(combin_transfo) == 1) {current_cmbn = as.data.frame(combin_transfo[grp_index, ])
      }else{current_cmbn = combin_transfo[grp_index, ]}
      bool = check_combn(transfo, data, current_cmbn)

      list_cmbn = current_cmbn[which(bool),]

      info_all_combi = getCombiFormula(data, transfo, list_cmbn)

      as_tibble(do.call(rbind, info_all_combi))
    }else{
      tibble(Molecule = data$name, Transformation = data$name, Formula = gsub(" ", "", data$formula),
                 Diff = "", Nb_Transfo = 0)
    }
  })
  info_all_combi = unique(do.call(rbind, info_all_combi))
  
  # Close cluster
  stopCluster(cl)
  closeAllConnections()


  #### Chromatogram : check signal for each metabolites ####

  cat(paste0("### Start chromatogram for ", data$name, " ###\n"))

  cat("### Create Directories ###\n")

  #Create output directories in the working directory for the current molecule
  out_current_mlc = paste0(opt$output, "/out_", data$name, "_", sub('\\..*$', '', basename(data$ms_file)), 
                           "_", format(Sys.time(), "%d%m%y_%H%M"))

  lapply( c("/", "/POS", "/NEG", "/POSNEG", "/input"),
          function(newdir) dir.create(paste0(out_current_mlc, newdir), showWarnings = FALSE) )



  #### SEARCH REFERENCE MS2 ####



  cat("### Get reference MS2 ###\n")

  # Load mzML files and adduct info
  ms_file = c() ; adduct = c()
  if(do$minus){
    #get adduct in neg
    adduct$minus = list(formula = data$adduct_minus,
                        mz = ifelse(data$adduct_minus=="H",
                                    -getMolecule(data$adduct_minus)$isotopes[[1]][1,1],
                                    getMolecule(data$adduct_minus)$isotopes[[1]][1,1]))
    ms_file$minus = readMSData(data$ms_file, mode="onDisk") %>% filterPol(polarity = 0)
  }

  if(do$plus){
    #get adduct in pos
    adduct$plus = list(formula = data$adduct_plus,
                       mz = getMolecule(data$adduct_plus)$isotopes[[1]][1,1])
    ms_file$plus = readMSData(data$ms_file,mode="onDisk") %>% filterPol(polarity = 1)
  }

   # If the parameter is not empty
  if(data$ms2_reference_tsv != FALSE ) {
    if(file.exists(data$ms2_reference_tsv)){
      ref_ms2_tsv = read_tsv(data$ms2_reference_tsv, col_types = cols())
      ref_ms2_tsv$polarity = ifelse(ref_ms2_tsv$intensity > 0, "plus", "minus")
    }else{
      message("MS2 tsv provided doesn't exists, please check")
    }
  }
  ref_ms2 = getMS2Reference(names(do)[which(do==T)])

  cat(paste0("### Start foreach loop : ", length(unique(info_all_combi$Formula)), " metabolites ###\n"))
  
  # Save input for the current molecule
  file.copy(opt$config, paste0(out_current_mlc, "/input"))
  write_tsv(data_tsv[row, c(1:6)], paste0(out_current_mlc, "/input/molecule.tsv"))
  if(exists("ref_ms2_tsv")) file.copy(data$ms2_reference_tsv, paste0(out_current_mlc, "/input"))


  #### PARALLEL LOOP ####
  

  # OPEN CLUSTER
  cl_big = makeCluster(cores)
  clusterApply(cl_big, 1:cores, opt, fun= function(x,opt){
    source(opt$config)
    source("util.R")
    source("plot_functions.R")
    library(pacman)
    
    p_load("tibble", "MSnbase", "plyr", "ggplot2", "ggpubr", "Rdisop", "dplyr", "ggrepel",
           "stringr", "htmlwidgets", "readr")
  })
  
  # Initialisation of cluster environnement : export variables and packages
  clusterExport(cl_big, ls())
  
  registerDoParallel(cl_big) # do parallel analysis
  
  #Do in parallel chromatogram for the current molecule for each formula obtained
  
  BIG_TABLE = foreach::foreach(current_formula = unique(info_all_combi$Formula), .combine = rbind ) %dopar% {
    
    current_mlc = getMolecule(current_formula)
    
    #Check if the current molecule is valid based on the formula
    if(getValid(current_mlc)=="Valid"){
      
      
      ## Init plot legend
      
      # Select the shortest transformation to obtain the current formula
      current_transfo = info_all_combi$Transformation[which(info_all_combi$Formula == current_formula)]
      split_transfo = lapply(strsplit(current_transfo, ";"), length)
      optim_transfo = current_transfo[which(split_transfo == min(unlist(split_transfo)))]
      
      optim_transfo = ifelse(length(optim_transfo) > 1, optim_transfo[which.min(nchar(optim_transfo))], optim_transfo)
      nb_transfo = info_all_combi$Nb_Transfo[which(info_all_combi$Transformation == optim_transfo)]
      
      current_diff = unique(info_all_combi$Diff[which(info_all_combi$Formula == current_formula)])
      
      plot_subtitle = paste("Shorter transformation :", optim_transfo,
                            ifelse(current_diff != "", paste("\nFormula differences :", current_diff), ""))
      
      
      combine_plot_chromato = data.frame() ; combine_ms_data = data.frame()
      for(pol in c("minus", "plus")){
        # MINUS : chromatogram and mass spectrum data
        if(do[[pol]]){
          plot_chromato = getChromato(ms_file[[pol]], current_mlc, adduct[[pol]], E_MASS[[pol]], min_intensity = min_peak_intensity,
                                            mz_precision = mz_ppm, min_scan = nb_scan, min_mz = minimum_mz)
          ms_data = getMassSpectrum(ms_file[[pol]], current_mlc, adduct[[pol]], plot_chromato, E_MASS[[pol]], peak_wdw = rt_windows,
                                          min_intensity = min_peak_intensity)
          
          # If there is a mass spectrum then give index to peak
          if(!empty(ms_data)){
            plot_chromato$index = 0 ; plot_chromato$polarity = pol
            plot_chromato$index[which((plot_chromato$rt/60) %in% unique(ms_data$rtime) &
                                              plot_chromato$isotope == 1)] = unique(ms_data$index)
            ms_data$polarity = pol
            combine_plot_chromato = rbind(combine_plot_chromato, plot_chromato)
            combine_ms_data = rbind(combine_ms_data, ms_data)
          }
        }
      }
      

      ### BUILD AND SAVE FIGURES ###

      

      # If there is a signal
      if(!empty(combine_plot_chromato)){

        # Get the polarity of the signal
        if(length(unique(combine_plot_chromato$polarity)) == 2){
          out_polarity = "POSNEG"
        }else if(unique(combine_plot_chromato$polarity) == "minus"){
          out_polarity = "NEG"
        }else if(unique(combine_plot_chromato$polarity) == "plus") {
          out_polarity = "POS"
        }

        plot_path = paste0(out_current_mlc, "/", out_polarity)

        tmp_table = data.frame()

        # for each polarity with signal
        for(pol in unique(combine_plot_chromato$polarity)){
          perc_common_peak = c() ; dotp_ms2 = c() ; mono_ppm = c()

          ms2_ref = ref_ms2[[pol]][2]  ;  mz_ref = ref_ms2[[pol]][1]

          plot_chromato = combine_plot_chromato[combine_plot_chromato$polarity == pol, ]
          ms_data = combine_ms_data[which(combine_ms_data$polarity == pol), ]

          ## Build plots
          group_ms = split(ms_data, ms_data$rtime)
          
          # Contain the SVG figures name to remove false positive
          filenames = c()

          # For each chromatogram peak / retention time
          for(current_ms in group_ms) {

            ggp = doChromato(plot_chromato[which(plot_chromato$rt/60 < unique(current_ms$rtime) +2 &
                                                   plot_chromato$rt/60 > unique(current_ms$rtime) -2), ]) +
              labs(title = "Metabolite chromatogram",
                   subtitle = paste(plot_subtitle, "\nm/z : ", round(min(plot_chromato$mz),5)))

            ggp_ms = doMassSpectrum(current_ms)
            mono_ppm = c(mono_ppm, current_ms$exp_ppm[!is.na(current_ms$exp_ppm)][1])

            ggparr_ms = ggarrange(ggp, ggp_ms)
            
            if( 2 %in% unique(msLevel(ms_file[[pol]])) ){
              # If there is a reference MS2
              if(!is.null(ms2_ref)){
                data_ms2 = compareMS2(ms_file[[pol]], optim_transfo, ms2_ref[[1]], mz_ref[[1]], unique(plot_chromato$mz),
                                      unique(current_ms$rtime), wdw_mz = wdw_mz_ms2, tol_mz = mz_ppm)
                # If there is a signal in MS2
                ggp_ms2 = doMS2(data_ms2)
                dotp_ms2 = c(dotp_ms2, ifelse(is.null(data_ms2$dotp_ms2), 0, data_ms2$dotp_ms2))
                
                if(length(data_ms2) > 0){
                  #write table with match m/z
                  filename_tsv = paste0("table_", current_formula, "_", unique(current_ms$index), "_",
                                        ifelse(out_polarity=="POSNEG", pol, ""), "_", out_polarity,".tsv")
                  write_tsv(data_ms2$match, file.path(plot_path, filename_tsv))
                }
              }else{
                data_ms2 = compareMS2(ms_file[[pol]], optim_transfo, mz_exp = unique(plot_chromato$mz)[1],
                                      rt_exp = unique(current_ms$rtime), exists_ref_ms2 = F, wdw_mz = wdw_mz_ms2)
                ggp_ms2 = doMS2(data_ms2, exists_ref = F)
                dotp_ms2 = c(dotp_ms2, 0)
              }
              
              #final ggplot with chromatogram and mass spectrum
              ggp_tot = ggarrange(ggparr_ms, ggp_ms2, nrow=2)
            }else{
              data_ms2 = c()
              dotp_ms2 = c(dotp_ms2, 0)
              ggp_tot = ggparr_ms
            }
            
            # Annote figure with MS2 score based on common peaks
            perc_peak = 0
            # If there is ms2 data and there is at least 1 match with reference
            if(length(data_ms2) > 1 && nrow(data_ms2$match) > 0) {
              nb_ref_peaks = length(grep("ref", data_ms2$data$type)) #number of reference peak 
              perc_peak = round(max(data_ms2$match$index_peak) / nb_ref_peaks*100,2)

              ggp_tot = annotate_figure(ggp_tot, bottom = text_grob(
                paste0("Percentage of common peaks : ", perc_peak, "%", " / Dotp : ", round(dotp_ms2[length(dotp_ms2)],3)), size = 8 ))
            }

            # Plot title
            ggp_tot = annotate_figure(ggp_tot, top = text_grob(
              paste(data$name, ":", current_formula,
                    ifelse(pol == "plus", "+", "-"), 
                    adduct[[pol]]$formula,
                    "/ peak number", unique(current_ms$index))
              , size = 9))

            perc_common_peak = c(perc_common_peak, perc_peak)

            filename = paste0(current_formula, "_", unique(current_ms$index),"_",
                              ifelse(out_polarity=="POSNEG", pol, ""), "_", out_polarity,".svg")
            ggsave(file.path(plot_path, filename), ggp_tot)
            
            filenames = append(filenames, filename)

          } # end for retention time

          peak_info = unique(ms_data[, c("rtime", "dotp", "rscore", "abscore", "peak_intensity", "index")])
          
          tmp_table = rbind.data.frame(tmp_table,
                                       cbind.data.frame(name = data$name, formula = current_formula, polarity = pol,
                                                        adduct = unique(plot_chromato$adduct), mz = round(min(plot_chromato$mz),5),
                                                        transfo = optim_transfo, diff = current_diff, rt = round(peak_info$rtime, 3),
                                                        abscore = round(peak_info$abscore, 3), dotp_ms2 = round(dotp_ms2, 3), 
                                                        common_ms2_peak = perc_common_peak, mono_ppm = mono_ppm,
                                                        intensity = peak_info$peak_intensity, index_peak = peak_info$index,
                                                        nb_transfo = nb_transfo,
                                                        filepath = plot_path, filename = filenames))

        } #end for polarity

        tmp_table

      }
    }

  } #end foreach
  
  # Close cluster
  stopCluster(cl)
  closeAllConnections()
  
  
  if(!empty(BIG_TABLE)){
    BIG_TABLE$nb_transfo = as.numeric(as.character(BIG_TABLE$nb_transfo))

    ## Compute score

    # Search for parent signal in both polarities
    parent_intensity = list(plus="", minus="")
    for(pol in c("plus", "minus")){
      if(do[[pol]]){
        parent = which(BIG_TABLE$nb_transfo==0 & BIG_TABLE$polarity==pol)
        # If needed fix the parent retention time
        if(length(parent) > 1){
          if( length(ref_ms2[[pol]]) == 3 && !is.na(ref_ms2[[pol]][3]) ){
            good_parent = which.min( abs( BIG_TABLE$rt[parent] - ref_ms2[[pol]][[3]] ) )
            parent_intensity[[pol]] = BIG_TABLE$intensity[parent[good_parent]]
          }else{
            good_parent = which.max(BIG_TABLE$intensity[parent])
            parent_intensity[[pol]] = BIG_TABLE$intensity[parent[good_parent]]
          }
        }else if(length(parent) == 0){
          parent_intensity[[pol]] = 0
          message("No parent found for this molecule in ", pol)
        }else{
          parent_intensity[[pol]] = BIG_TABLE$intensity[parent]
        }
      }
    }

    # Fix limit value to intensity ratio
    ratio_intensity = ifelse(BIG_TABLE$polarity == "plus",
                             ifelse(BIG_TABLE$intensity/parent_intensity$plus == Inf,
                                    0, BIG_TABLE$intensity/parent_intensity$plus),
                             ifelse(BIG_TABLE$intensity/parent_intensity$minus == Inf,
                                    0, BIG_TABLE$intensity/parent_intensity$minus))
    ratio_intensity[ratio_intensity>1] = 1

    # Compute Score and penalty based on the number of transformations
    penalty = ifelse(BIG_TABLE$nb_transfo == 0, 1, 1/sqrt(BIG_TABLE$nb_transfo))
    BIG_TABLE$rintensity = ratio_intensity
    BIG_TABLE$score = round(penalty * ( 1/2 * BIG_TABLE$common_ms2_peak / 100 * BIG_TABLE$dotp_ms2 +
                                          1/2 * BIG_TABLE$abscore * ratio_intensity )
                            , 3 )
    
    BIG_TABLE$score[BIG_TABLE$score == Inf] = 1
    BIG_TABLE = arrange(BIG_TABLE, desc(score))
    
    # Remove isobare (close M/Z)
    for(tab_row in 1:nrow(BIG_TABLE)){
      isobare = close_match(x= BIG_TABLE$mz, target=BIG_TABLE$mz[tab_row], tolerance = BIG_TABLE$mz[tab_row]/1e6*mz_ppm)
      true_isobare = which(BIG_TABLE$formula[isobare] != BIG_TABLE$formula[tab_row] & 
                             BIG_TABLE$polarity[isobare] == BIG_TABLE$polarity[tab_row])
      if(length(true_isobare) > 0) BIG_TABLE = BIG_TABLE[-isobare[true_isobare],]
    }
    
    BIG_TABLE[,c("formula", "filepath", "filename")] = apply(BIG_TABLE[,c("formula", "filepath", "filename")], 2, as.character)
    # Group by RT and remove isotopes false positive
    BIG_TABLE_FINAL = split(BIG_TABLE, BIG_TABLE$rt) %>%
      lapply(function(current_rt){
        
        index = 1
        while (index <= nrow(current_rt)) {
          # Get the mono isotopic information
          mol = getMolecule(current_rt$formula[index])
          mz_iso = mol$isotopes[[1]][1, 2:6] + adduct[[current_rt$polarity[index]]]$mz + E_MASS[[current_rt$polarity[index]]]
          # Search for false positive with MZ close to isotopes m/z
          false_pos = unlist( lapply(mz_iso, function(iso) which(between(current_rt$mz, iso-iso/1e6*mz_ppm, iso+iso/1e6*mz_ppm))) )
          
          # Remove figures and tables for false positive
          names = unlist(strsplit(current_rt$filename[false_pos], ".svg"))
          file.remove(
            list.files(path = unique(current_rt$filepath[false_pos]), pattern = paste0(names, collapse="|"), full.names = T)
          )
          # Save the results
          if(length(false_pos) > 0){
            current_rt = current_rt[-false_pos,]
          }
          index = index +1
        }
        
        return(current_rt[,!names(current_rt) %in% c("filepath", "filename")])
        
      }) 
    
    BIG_TABLE_FINAL = do.call(rbind, BIG_TABLE_FINAL)
    
    
    BIG_TABLE_FINAL$polarity = ifelse(BIG_TABLE_FINAL$polarity == "plus", "Positive mode", "Negative mode")

    # Write output table
    BIG_TABLE_FINAL %>%
      write_tsv(paste0(out_current_mlc, "/out_", sub('\\..*$', '', basename(data$ms_file)), ".tsv"))
  }

  stop=Sys.time()

  cat(paste("### Execution time :", data$name, round(difftime(stop, start, units="mins"),2), "mins ###\n"))

} #end for loop
