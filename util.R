# # #
# Functions for MetIDfyR
# author : Agnes Barnabe
# # #

###
# Get available cores in the system
###
availableCore = function(){
  # If on Windows, simply search for the number of cores on the device
  if(.Platform$OS.type == "windows"){
    parallel::detectCores()
  # Else, use top to get free cores
  }else{
    usedCores = system("top -b -n 1 1 | grep Cpu | cut -c 9-13", intern=T)
    usedCores = as.numeric(gsub(",", ".", usedCores))
    return(detectCores() - sum(usedCores > 80))
  }
}

###
# Estimate mass shift for each transformation
# Need a list of transformation
###
getShift = function(transfo_path){
  transfos = read_tsv(transfo_path, col_types = cols())
  
  removeshift<-rep(0, nrow(transfos))
  addshift<-rep(0, nrow(transfos))
  
  for(n in 1:nrow(transfos)){
    if( !is.na(transfos$remove[n]) ){
      removeshift[n] = getMolecules(transfos$remove[n])$isotopes[[1]][1,1]
    }
    if( !is.na(transfos$add[n]) ){
      addshift[n] = getMolecules(transfos$add[n])$isotopes[[1]][1,1]
    }
  }
  
  transfos$mass_change = addshift - removeshift
  
  return(transfos)
  
}

###
# Filter of the polarity  
# Need a mzML object 
###

filterPol = function(object, polarity) {
  if (missing(polarity)) return(object)
  object <- object[polarity(object) %in% polarity]
  object
}

###
# Search closest value from target in a vector
# Need : vector x, target value to match and tolerance 
# Return : sorted value depending of the target distance if tolerance = NULL,
# if tolerance != NULL, return the vector with 1 if the value is close to the target else NA
###

close_match = function(x, target, tolerance=NULL){
  if(is.null(tolerance)){
    x[order(x-target)]
  }else{
    match = c()
    for(tar in target){
      match = append(match, which(ifelse(x <= tar+tolerance & x >= tar-tolerance, 1, NA) == 1))
    }
    return(sort(match))
  }
}

###
# Get the molecule composition based on the formula
# Need : the molecular formula with space between each atom
# Return : a list with the count of each atom in the formula
###

mol_compo = function(formula){
  #extract number of each atoms in the formula
  tmp = gsub("\\b([A-Za-z]+)\\b","\\11",formula)
  value = as.integer(str_extract_all(tmp,"\\d+")[[1]])
  names(value) = str_extract_all(tmp,"[A-Za-z]+")[[1]]
  
  #create a list with atom as name and count as value
  compo = aggregate(values ~ ind, stack(value), sum)
  compo$ind = as.character(compo$ind)
  return(compo)
}


###
# Fill the table with the atom count of the molecule
# Need : a table containing the molecules informations
# Return : a table completed with the count of each atom
###

fill_table = function(tab){
  #for each molecule, get formula composition
  tmp = data.frame()
  for(row in 1:nrow(tab)) {
    compo = mol_compo(tab$formula[row])
    tmp[row, compo$ind] = compo$values
  }
  
  tmp[is.na(tmp)]=0
  return(tmp)
}


###
# Replacement transformation
# Need : the formula, what to add and to remove
# Return : the final molecule
###

replaceMolecules = function(formula, add, remove){
  tmp = subMolecules(formula, remove)
  tmp = addMolecules(tmp$formula, add)
  return(tmp)
}


###
# Get formula after each transformation in a combination
# Need : molecule informations, list of transformations and the combination index
# Return : list of formula with, as names, the combined name of transformations used
###

getCombiFormula = function(data, transfo, list_cmbn){
  
  info_all_combi = plyr::alply(list_cmbn, 1, function(cmbn) {
    
    molecule = getMolecules(data$formula)
    
    # Init table with 0 transformation 
    i=1 ; info_transfo = tibble(Molecule = data$name, Transformation = data$name, Formula = molecule$formula, 
                                Diff = "", Nb_Transfo = 0)
    
    
    continue_cmbn = check_combn(transfo, data, matrix(cmbn[1], nrow=1)) 
    # for each transformation of the combination
    # while i <= number of transformations, combinaison is valid and the number of phase 2 transformation is <= 2
    while(i <= length(cmbn) & continue_cmbn & sum(transfo$phase[cmbn[1:i]]==2) <= 2){
      index=cmbn[i]
      #add, remove or replace atoms depending on the transformation type
      res_transfo = switch(transfo$type[index],
                           'add'= addMolecules(info_transfo$Formula[i], getMolecules(transfo$add[index])$formula),
                           'remove' = subMolecules(info_transfo$Formula[i], getMolecules(transfo$remove[index])$formula),
                           'replace' = replaceMolecules(info_transfo$Formula[i], getMolecules(transfo$add[index])$formula, 
                                                        getMolecules(transfo$remove[index])$formula)
      )
      
      info_transfo = info_transfo %>% add_row(Molecule = data$name, Transformation = paste0(transfo$name[cmbn[1:i]], collapse = ";"),
                                              Formula = res_transfo$formula, 
                                              Diff = getCount(transfo, data, cmbn[1:i], "difference"), 
                                              Nb_Transfo = i )
      
      i=i+1
      # Check if the current combination is possible
      continue_cmbn = check_combn(transfo, data, matrix(cmbn[1:i], nrow=1))
    }
    
    if(nrow(info_transfo)){
      return(info_transfo)
    }
    
  })
  return(info_all_combi)
}


###
# getMolecule using the formula with spaces 
# Need : the formula
# Return : the molecule informations
###

getMolecules = function(formula){
  return(getMolecule(gsub(" ", "", formula)))
}


###
# getMass for a specific isotope
# Need : the molecule and the isotope number
# Return : the isotope exact mass
###

getMasse = function(molecule, isotope){
  return(molecule$isotopes[[1]][1,isotope])
}


###
# Get the type of transformation
# Need : table of transformations informations
# Return: : the type of transformation for each one
###

getType = function(transfo){
  type = rep(NA, nrow(transfo))
  type[is.na(transfo$remove)] = "add"
  type[is.na(transfo$add)] = "remove"
  type[is.na(type)] = "replace"
  return(type)
}

###
# Get the formula after transformation
# Need : table with transformations informations, table with current molecule composition 
#   the combination selected and the response type
# Return : the new composition after combination 
#   OR a string with the composition difference between formula
###

getCount = function(transfo, data_count, cmbn, response = c("compo", "difference")){
  
  # Get number of atoms to add and remove
  remove_trans = paste0(transfo[cmbn,]$remove, collapse=" ") #remove transformation
  remove_count = mol_compo(remove_trans)
  
  add_trans = paste0(transfo[cmbn,]$add, collapse=" ") #adding transformation
  add_count = mol_compo(add_trans)
  
  # Add columns if doesn't exists
  tmp = unique(c(add_count$ind, remove_count$ind))
  data_count[setdiff(tmp, names(data_count))]=0
  
  initial_count = data_count
  nums = unlist(lapply(initial_count, is.numeric)) #get numeric columns
  nums_name = names(which(nums))
  
  # Get count after transformations
  data_count[remove_count$ind] = data_count[remove_count$ind] - remove_count$values
  data_count[add_count$ind] = data_count[add_count$ind] + add_count$values
  
  if(response == "compo"){
    
    if("NA" %in% names(data_count)) data_count = subset(data_count, select=-c(`NA`))
    return( data_count )
    
  }else if(response == "difference"){
    
    # Get difference btw initial formula and post-combination formula 
    diff = data_count[nums_name] - initial_count[nums_name]
    if("NA" %in% names(diff)) diff = subset(diff, select=-c(`NA`)) #remove column NA
    
    sup = paste0( names(diff)[diff<0], ifelse(abs(diff[diff < 0]) > 1, abs(diff[diff < 0]), "" ) )
    add = paste0( names(diff)[diff>0], ifelse(abs(diff[diff > 0]) > 1, abs(diff[diff > 0]), "" ) )
    
    diff_exp = paste0( ifelse(length(sup) > 0,
                              paste0(" - ", do.call(paste0, as.list(sup))), "")
                       ,ifelse(length(add) > 0,
                               paste0(" + ", do.call(paste0, as.list(add))), "")
    )
    
    return(diff_exp)
  }
  
}


###
# List all possible combination of transformations with the resulting formula, 
# the composition difference and the number of transformation.
###
createCombiTable = function(current_grp, current_data){
  # If there is at least one transformation
  if(ncol(combin_transfo) >= 1) {
    grp_index = which(group == current_grp)
    if(ncol(combin_transfo) == 1) {
      current_cmbn = as.data.frame(combin_transfo[grp_index, ])
    }else if(nrow(combin_transfo) == 1){
      current_cmbn = rbind(combin_transfo[grp_index, ])
    }else{
      current_cmbn = combin_transfo[grp_index, ]
    }
    bool = check_combn(transfo, current_data, current_cmbn)
    
    list_cmbn = current_cmbn[which(bool),]
    
    info_all_combi = getCombiFormula(current_data, transfo, current_cmbn)
    
    as_tibble(do.call(rbind, info_all_combi))
  }else{
    tmp = tibble(Molecule = current_data$name, Transformation = current_data$name, Formula = gsub(" ", "", current_data$formula),
                 Diff = "", Nb_Transfo = 0)
    
    tmp
  }
  
}

###
# Check if the combinations of transformations are possible
# Need : table with transformations informations, table with current molecule composition 
#   and the list of all existing combinations
# Return : table of boolean, TRUE if the combination is possible
###

check_combn = function(transfo, data, combin){
  #check if the combination is possible according to the molecule composition
  do_combn = apply(combin, 1, function(cmbn){
    
    #return boolean, TRUE if there is no negative count so the combination is possible
    tab_count = getCount(transfo, data, cmbn, "compo")
    return( sum(Filter(is.numeric, tab_count) < 0) == 0 )
    
  })
  
  return( do_combn )
}

###
# Calculate dotp value between theoretical and experimental mass spectrum
# Need : theoretical and experimental mass spectrum
# Return : dotp values and mz ppm difference between matched peaks OR False if there is less than two isotopes
###
doDotp = function(th, exp, mz_ppm = 10){
  
  # Keep ratio >= 5% or at least two values
  if(length(which(th$perctot >= 0.05)) > 2){
    th = th[which(th$perctot >= 0.05),]
  }else{ th = th[1:2,] }
  
  exp_peak = tibble(th_inten = integer(), exp_inten = integer(), ppm = integer(), mz = integer()) ; x = 1
  tmp = close_match(exp$mz, th$mz[x], tolerance = th$mz[x]/1e6*mz_ppm)
  
  while(length(tmp)){
    
    # Keep the max intensity as the real peak
    tmp = tmp[which.max(exp$intensity[tmp])]
    # Return the relative abundance
    exp_peak = add_row(exp_peak, th_inten = th$perctot[x], exp_inten = exp$perctot[tmp], 
                       ppm = (exp$mz[tmp]-th$mz[x])/th$mz[x]*1e6, mz = exp$mz[tmp] )
    
    x=x+1
    # Get the closest experimental mz from the theoretical mz
    tmp = close_match(exp$mz, th$mz[x], tolerance = th$mz[x]/1e6*mz_ppm)  
  }
  
  #If there is only the monoisotopic, can't calculate the dotp so skip the RT
  if(nrow(exp_peak) <= 1) return(FALSE)
  
  # a = theoretical ratio ; b = experimental ratio
  a = exp_peak[,1] ; b = exp_peak[,2]/max(exp_peak[,2])
  
  return( list(dotp = sum(a*b)/sqrt(sum(a^2)*sum(b^2)),
               rscore = 1-sum(abs(a-b))/sum(a),
               abscore = 1-sum(abs(a-b)),
               exp_ppm = exp_peak[,c(3,4)]) )
}

### 
# Get data to build Mass Spectrum if there is signal
# Need : ms file, molecule informations, adduct and chromatogram data
# Return : a table with data to build mass spectrum for each local maximum in the chromatogram
###
getMassSpectrum = function(ms_file, current_mlc, adduct, chrom_dataframe, masse_electron, 
                           nb_iso_th=3, min_intensity=5e5, peak_wdw = 5){
  
  # If there is a chromatogram, search for mass spectrum
  if(!empty(chrom_dataframe)) {
    
    chromato_1 = new("Spectrum1", mz = chrom_dataframe$rt[which(chrom_dataframe$isotope==1)]/60, 
                     intensity = chrom_dataframe$intensity[which(chrom_dataframe$isotope==1)], centroided = F)
    
    ### LOCAL MAXIMUM in chromatogram for isotope 1
    p = pickPeaks(chromato_1, halfWindowSize = peak_wdw)
    #keep only maximum above the minimal intensity
    realpeak = which(p@intensity > min_intensity)
    
    ### Get mass spectrum for selected chromatogram peak
    spectrum_data = data.frame() 
    index_peak = 1
    
    # For each retention time / each chromatogram peak search for the mass spectrum
    for(current_peak in realpeak){
      
      # Extract rtime and intensity of the chromatogram peak
      rtime = p@mz[current_peak] ; inten_peak = p@intensity[current_peak]
      
      ### THEORETICAL RELATIVE ABUNDANCE
      nb_iso = 1 
      abund_th = matrix(ncol=2, nrow=0) 
      colnames(abund_th)=c("mz", "intensity")  
      
      # get isotopes informations for nb_iso_th isotopes 
      while(nb_iso <= nb_iso_th){
        abund_th = rbind(abund_th, getIsotope(current_mlc, nb_iso))
        nb_iso = nb_iso +1
      }
      
      ### EXPERIMENTAL RELATIVE ABUNDANCE
      if(!empty(abund_th)){ 
        # add adduct and electron mass to isotopes mz and convert total to relative abundance
        abund_th = as_tibble(abund_th)
        abund_th$mz = abund_th$mz + adduct$mz + masse_electron
        abund_th$perctot = abund_th$intensity/max(abund_th$intensity)
        
        mass_spectrum = which(round(rtime(ms_file),5) == round(rtime*60,5))
        
        current_spectra = ms_file[[ mass_spectrum ]]
        mass_wd = which( current_spectra@mz > min(abund_th$mz)-1.5 & current_spectra@mz < max(abund_th$mz)+4.5 )
        
        # if there is signal in the mz interval
        if(!length(mass_wd) == 0 ){
          abund_exp = cbind(mz = current_spectra@mz[mass_wd], intensity = current_spectra@intensity[mass_wd])
          abund_exp = as_tibble(abund_exp)
          
          # Get the mono isotopic pic to compute relative abundance
          mono = close_match(abund_exp$mz, abund_th$mz[1], tolerance = abund_th$mz[1]/1e6*mz_ppm)
          # Keep the max intensity as the real peak
          mono = mono[which.max(abund_exp$intensity[mono])]
          
          abund_exp$perctot = abund_exp$intensity/abund_exp$intensity[mono]
          
          # Compare experimental and theorical spectrum via dotproduct 
          scores = doDotp(abund_th, abund_exp, mz_ppm)
          
          if(length(scores) > 1){
            
            #Combine output table for isotopic pattern 
            spectrum_data = rbind(spectrum_data, 
                                  cbind(abund_exp, rtime = rtime, type = "exp", dotp = scores$dotp, 
                                        rscore = scores$rscore, abscore = scores$abscore, exp_ppm = NA,
                                        index=index_peak, peak_intensity = inten_peak, row.names=NULL),
                                  cbind(abund_th, rtime = rtime, type = "th", dotp = scores$dotp, 
                                        rscore = scores$rscore, abscore = scores$abscore, exp_ppm = NA,
                                        index=index_peak, peak_intensity = inten_peak, row.names=NULL))
            #Add mz ppm difference to match peak
            spectrum_data$exp_ppm[spectrum_data$mz %in% scores$exp_ppm$mz & 
                                    spectrum_data$index == index_peak] = scores$exp_ppm$ppm
            
            index_peak = index_peak + 1
          }
          
        }
      }
    } # end for 
    
    return(spectrum_data)
    
  }
  #Else return an empty dataframe
  return(data.frame())
  
}

###
# Search current molecule in the sample (ms_file)
# Need : ms file, molecule, adduct information, electron mass
# Return : chromatogram data
###

getChromato = function(ms_file, molecule, adduct, masse_electron, 
                       max_nb_iso = 5, min_mz = 200, min_intensity = 1e5, min_scan = 5, mz_precision = 5){
  
  #get mz for each isotopes 
  tmp_iso = which(getIsotope(molecule)[2,] >= 0.1)
  list_iso = matrix(ncol=3, nrow=0) 
  colnames(list_iso)=c("mz", "intensity", "nb_iso")  
  for(i in tmp_iso) list_iso = rbind(list_iso, c(getIsotope(molecule, i), i))
  
  list_iso = as_tibble(list_iso)
  list_iso$perctot = list_iso$intensity/list_iso$intensity[1]
  #add adduct and electron mass to molecule mz
  mz_adduct = list_iso$mz + adduct$mz + masse_electron 
  
  chrom_dataframe = tibble(rt = double(), intensity = double(), nb_scan = double(),
                           isotope = integer(), adduct = character(), mz = integer())
  
  #for each isotope
  for(iso in 1:length(mz_adduct)){
    signal=F
    #if isotope mass superior to the threshold
    if(mz_adduct[iso] > min_mz) {
      
      mz_confidence = c(mz_adduct[iso]-mz_adduct[iso]/1e6*mz_precision,
                        mz_adduct[iso]+mz_adduct[iso]/1e6*mz_precision) #get mz interval
      
      chrs <- chromatogram(ms_file, mz = mz_confidence, missing = 1) #extract chromatogram info
      #if chromatogram for the polarity
      if(length(chrs) > 0){
        intensity_iso = chrs[[1]]@intensity
        adapt_intensity = min_intensity*list_iso$perctot[iso]
        
        nb_above_thr = length(which(intensity_iso > adapt_intensity - adapt_intensity*0.2))
        
        #if enough signal intensity and number of scans
        if(nb_above_thr >= min_scan){
          chrom_dataframe = add_row(chrom_dataframe, rt = chrs[[1]]@rtime, intensity = intensity_iso, nb_scan = nb_above_thr,
                                    isotope = list_iso$nb_iso[iso], adduct = adduct$formula, mz = mz_adduct[iso])
          signal = T
        }
        
      }
    }
    if(!signal) break()
    
  }
  chrom_dataframe$isotope = factor(chrom_dataframe$isotope, levels=c(1,2,3))
  
  return(chrom_dataframe)
  
  
}

###
# Search reference MS2 in the mzML
# Need : current molecule information
# Return : list with the ms2 and mz of the parent drug in both polarities
###
getMS2Reference = function(data, polarities = c("plus", "minus")){
  ref_mlc = getMolecules(data$formula)
  plot_chromato = list() ; ms_data = list() ; RT = list()
  
  for(pol in polarities){
    if(do[[pol]]){
      plot_chromato[[pol]] = getChromato(ms_file[[pol]], ref_mlc, adduct[[pol]], E_MASS[[pol]], min_intensity = min_peak_intensity,
                                        mz_precision = mz_ppm, min_scan = nb_scan)
      ms_data[[pol]] = getMassSpectrum(ms_file[[pol]], ref_mlc, adduct[[pol]], plot_chromato[[pol]], E_MASS[[pol]], peak_wdw = rt_windows,
                                      min_intensity = min_peak_intensity)
    }
    
    if( !empty(ms_data[[pol]]) ) RT[[pol]] = unique(ms_data[[pol]]$rtime)
    # Use given rtime if provided, else keep the max intensity peak
    if(length(RT[[pol]]) > 1){
      RT[[pol]] = RT[[pol]][which.max(unique(ms_data[[pol]]$peak_intensity))]
    }
  }
  
  # If there is pos and neg data && RT are different, keep the max TIC
  if(length(RT) == 2 && round(RT$plus,2) != round(RT$minus, 2)){
    ind_max = which.max(c(unique(ms_data$plus$peak_intensity[which(ms_data$plus$rtime == RT$plus)]),
                        unique(ms_data$minus$peak_intensity[which(ms_data$minus$rtime == RT$minus)])))
    RT[[polarities[-ind_max]]] = RT[[polarities[ind_max]]]
  }
  
  # Search for MS2 reference in both polarities 
  ref_data = list(plus = c(), minus = c())
  for(pol in polarities){
    if( do[[pol]] && 2 %in% unique(msLevel(ms_file[[pol]])) ){
      rt = RT[[pol]]
      
      tmp_mz = unique(plot_chromato[[pol]]$mz)
      mz_ref = tmp_mz[1]
      # If there is a ms2 spectra provided by the user in the current polarity
      if(exists("ref_ms2_tsv") && sum(ref_ms2_tsv$polarity == pol) > 0){
        ms2_ref = new("Spectrum2", mz = ref_ms2_tsv$mz[ref_ms2_tsv$polarity == pol],
                            intensity = abs(ref_ms2_tsv$intensity[ref_ms2_tsv$polarity == pol]))
        ref_data[[pol]] = c(mz = mz_ref, ms2 = ms2_ref)
        
      }else if( length(rt) >=1 ){
        
        ms2_ref = close_match(precursorMz(ms_file[[pol]]), mz_ref, tolerance = wdw_mz_ms2/2)
        # Keep as reference the MS2 with the closest rtime
        
        ms2_ref = ms2_ref[close_match(rtime(ms_file[[pol]])[ms2_ref], rt*60, tolerance = 5)]

        ms2_ref = ms2_ref[ which.max(unlist(lapply(ms2_ref, function(x) ms_file[[pol]][[x]]@tic))) ]
        
        if(length(ms2_ref) > 0){
          ms2_ref = ms_file[[pol]][[ms2_ref]]
          ref_data[[pol]] = c(mz = list(tmp_mz), ms2 = ms2_ref, rt = ms2_ref@rt/60)
        }
      }else{
        message("No reference MS2 spectra found in ", pol)
      }
      
    }else if( do[[pol]] && ! 2 %in% unique(msLevel(ms_file[[pol]])) ){
      message("No reference MS2 spectra found in ", pol)
    }
  }
  
  return(ref_data)
  
}


###
# Compare MS2 
# Need : ms file, reference MS2 & mz, experimental mz & rt
# Return : MS2 data reference vs experimental
###

compareMS2 = function(ms_file, optim_transfo, ref_ms2, mz_ref, mz_exp, rt_exp, wdsize = 5, exists_ref_ms2 = T, wdw_mz = 1,
                      tol_mz = 5, tol_rt = 5, noise_fraction = 0.05){
  
  # Search the closest M/Z to the precursor M/Z
  exp_ms2_num = close_match(precursorMz(ms_file), mz_exp[1], tolerance = wdw_mz/2)
  # Keep as metabolite MS2 the one with the closest rtime
  exp_ms2_num = exp_ms2_num[close_match(rtime(ms_file)[exp_ms2_num], rt_exp*60, tolerance = tol_rt)]
  
  if(length(exp_ms2_num) > 0){
    
    # Select the best MS2 if there are many
    exp_ms2_num = exp_ms2_num[ which.max(unlist(lapply(exp_ms2_num, function(x) ms_file[[x]]@tic))) ]
    # Get MS2, select only maximun scan in peaks and remove noise
    exp_ms2 = ms_file[[exp_ms2_num]]
    exp_ms2 = pickPeaks( removePeaks(exp_ms2, max(exp_ms2@intensity)*noise_fraction), halfWindowSize = wdsize ) 
    
    # If there is a reference ms2 search for COMMON PEAKS
    if(exists_ref_ms2){
      # If this is an extracted MS2, remove noise and pick peaks, else this is a given MS2 so do nothing
      if(length(ref_ms2@rt) > 0){
        ref_ms2 = pickPeaks( removePeaks(ref_ms2, max(ref_ms2@intensity)*noise_fraction), halfWindowSize = wdsize )
      }
      
      # Select the max intensity peaks except precursor
      precursor_ref = close_match( ref_ms2@mz, mz_ref, tolerance = mz_ref[1]/1e6*tol_mz ) 
      precursor_exp = close_match(exp_ms2@mz, mz_exp, tolerance = mz_exp[1]/1e6*tol_mz )
      
      ref_max_peaks = which(ref_ms2@intensity >= max(ref_ms2@intensity)*noise_fraction )
      
      plot_ref = normalize(ref_ms2) ; plot_exp = normalize(exp_ms2)
      type_ref = rep("ref", length(plot_ref@mz)) ; type_ref[precursor_ref] = "precursor_ref"
      type_exp = rep("exp", length(plot_exp@mz)) ; type_exp[precursor_exp] = "precursor_exp"
      
      nb_common_peak = 0 ; i = 1
      mz_exp_frag = exp_ms2@mz[exp_ms2@intensity > 0]
      
      label_peak = tibble(ref_mz = integer(), exp_mz = integer(), diff_mz = integer(), 
                          transfo = character(), index_peak = integer())
      index_ref = rep("", length(ref_ms2@intensity)) ; index_exp = rep("", length(mz_exp_frag))
      
      # For each reference peak selected, search for common peak in the experimental MS2
      for(ref_peak in ref_max_peaks) {
        tmp_mz_peak = ref_ms2@mz[ref_peak]
        
        # Get each transformation possible for the current combinaison of transformation
        split_transfo = strsplit(optim_transfo, ";")[[1]]
        cmbn = c()
        
        for(nb in 1:length(split_transfo)){
          tmp = arrangements::combinations(split_transfo, nb)
          cmbn = append( cmbn, split(tmp, row(tmp)))
        }
        
        cmbn = unique(cmbn)
        # Get m/z difference for each transformation
        mz_diff = lapply(cmbn, function(x){ 
          tmp_index = match(x, transfo$name)
          sum(transfo$mass_change[tmp_index])
        })
        names(mz_diff) = lapply(cmbn, function(x) paste0(x, collapse = ";"))
        mz_diff = mz_diff[!duplicated(mz_diff)]
        # add no transformation to list of m/z shift
        mz_diff = append(mz_diff, list(None = 0))
        
        # Search for matching peak with each m/z difference
        for(id_diff in 1:length(mz_diff)){
          match_peak = close_match(mz_exp_frag - mz_diff[[id_diff]], tmp_mz_peak, tolerance = tmp_mz_peak/1e6*tol_mz)
          # If there is a match
          if(length(match_peak) > 0) { 
            tmp_diff = round( mz_exp_frag[match_peak] - ref_ms2@mz[ref_peak], 4 )
            
            label_peak = add_row(label_peak, ref_mz = tmp_mz_peak, exp_mz = mz_exp_frag[match_peak],
                                 diff_mz = tmp_diff, transfo = names(mz_diff)[id_diff], index_peak = i)
            
            index_ref[ref_peak] = i ; index_exp[match_peak] = i
            
          }
        }
        
        #If there is no match for the current reference fragment
        if(!tmp_mz_peak %in% label_peak$ref_mz){
          label_peak = add_row(label_peak, ref_mz = tmp_mz_peak, exp_mz = 0,
                               diff_mz = 0, transfo = "No match", index_peak = 0)
        }
        
        if(i %in% label_peak$index_peak){ i=i+1}
        
      } #end for reference peak
      
      
      # Combine MS2 plot data
      
      combine_plot_data = as_tibble(
        rbind.data.frame(
          cbind.data.frame(mz = round(plot_ref@mz,5), intensity = plot_ref@intensity, type=type_ref, index = index_ref), 
          cbind.data.frame(mz = round(plot_exp@mz[plot_exp@intensity > 0],5), 
                           intensity = -plot_exp@intensity[plot_exp@intensity > 0], 
                           type=type_exp[plot_exp@intensity > 0], index = index_exp)))
      
      # Compute MS2 dotp and recalibrate shifted peaks
      exp_data = combine_plot_data[grep("exp", combine_plot_data$type),]
      exp_data$mz_calibrate = combine_plot_data$mz[grep("exp", combine_plot_data$type)]
      
      match_peak = match(round(label_peak$exp_mz, 5), exp_data$mz)
      match_peak = match_peak[!is.na(match_peak)]
      
      exp_int = exp_data$intensity
      # If there are matched peaks, unshift them
      if(length(match_peak)){
        # Unshift shifted peaks and experimental precursor 
        exp_data$mz_calibrate[match_peak] = exp_data$mz[match_peak] - label_peak$diff_mz[label_peak$index_peak > 0]
        
        # Sum intensity for double match peaks
        exp_int = aggregate(intensity ~ mz_calibrate, data=exp_data, FUN=sum)
        
        # Build new MS2 spectrum with calibrate data
        spec_exp = new("Spectrum2", mz = exp_int$mz_calibrate, intensity = -exp_int$intensity, centroided=T)
      }else{
        # Build new MS2 spectrum with calibrate data
        spec_exp = new("Spectrum2", mz = exp_data$mz_calibrate, intensity = -exp_int, centroided=T)
      }
      
      
      return(list( data = combine_plot_data[combine_plot_data$intensity != 0, ], 
                   match = label_peak, 
                   dotp_ms2 = compareSpectra(plot_ref, spec_exp, "dotp")))
      
    }else{ #else plot the experimental ms2
      return(normalize(exp_ms2))
    }
    
  }else{
    return(list())
  }
  
}
