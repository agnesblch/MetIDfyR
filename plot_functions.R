# # #
# Functions to build plot based on data produced with MetIDfyR
# author : Agnes Barnabe
# # #

###
# Build chromatogram plot
# Need : chromatogram data
# Return : a ggplot
###

doChromato = function(chrom_dataframe, line_size = 0.5, alpha = 0.4){
  
  myColors = RColorBrewer::brewer.pal(length(levels(chrom_dataframe$isotope)), "Set1")
  
  gpl = ggplot(chrom_dataframe, aes(x=rt/60, y=intensity, alpha=alpha,
                                    group=isotope, color=isotope))+
    geom_line(size=line_size) +
    geom_text(aes(label = ifelse(index > 0, index, '')),
              hjust = -0.2, vjust= 1, col="black", size=3) +
    theme_classic() +
    guides(alpha = FALSE) +  #don't show legend for alpha
    scale_y_continuous(expand=c(0,0),
                       labels = function(x) format(x, scientific = TRUE)) +
    scale_color_manual(values = myColors,
                       labels = levels(chrom_dataframe$isotope),
                       drop=FALSE)+
    labs(x="Retention time (min)",y="Intensity", color = "Isotope") +
    theme(legend.position = "bottom",
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-20,-10,0,-10),
          axis.text=element_text(size = 6, colour="black"),
          axis.ticks = element_line(colour="black"),
          plot.subtitle = element_text(size = 6), 
          plot.title = element_text(size = 7),
          legend.text = element_text(size=6),
          legend.title = element_text(size=6),
          axis.title = element_text(size=6),
          strip.text = element_text(size=6)) 
  
  return(gpl)
  
}

###
# Build mass spectrum plot
# Need : mas spectrum data
# Return : a ggarrange of ggplot
###

doMassSpectrum = function(data_ms, polarity = NULL){ 
  
  abund_exp = data_ms[which(data_ms$type == "exp"),]
  abund_th = data_ms[which(data_ms$type == "th"),]
  subtit = paste0( abund_exp$index, " / rtime : ", round(abund_exp$rtime, 2), " / dotproduct : ", round(abund_exp$abscore, 4),
                   " ", polarity)
  
  ggp = ggplot(data=abund_exp, aes(mz, perctot)) +
    geom_segment(data=abund_th, aes(xend=mz, yend=0, col = "red", alpha=0.8, linetype="Theoric"), size=2) +
    geom_segment(aes(xend=mz, yend=0, linetype = "Experimental"), col = "blue", alpha = 0.8) +
    theme_classic() + guides(alpha = F, colour  = F, size=F) +
    geom_label_repel(
      aes(mz, perctot, label = ifelse(!is.na(exp_ppm), 
                                      paste0(round(mz,5), "\n(", sprintf("%+g", round(exp_ppm, 3)), " ppm)"), '')),
      size = 3,
      point.padding = unit(0.1, "lines"),
      box.padding = unit(1, "lines"),
      segment.color = 'grey50'
    ) +
    labs(title = "Mass Spectrum theoric versus experimental", subtitle = subtit ) +
    xlab("") + ylab("") +
    scale_linetype_manual(name ="", values = c(1,1),
                          guide = guide_legend(override.aes = list(color=c("blue", "red")))) +
    scale_x_continuous(limits = c(abund_th$mz[1]-1.5, abund_th$mz[1]+4.5 )) +
    coord_cartesian(ylim=c(0, 1)) +
    theme(plot.subtitle = element_text(size=6),
          axis.text = element_text(size=6, colour="black"),
          legend.text = element_text(size=6), 
          plot.title = element_text(size = 7),
          legend.position = "bottom",
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-20,-10,0,-10),
          axis.ticks= element_line(colour="black"))
  
  
  return(ggp)
}


###
# Build MS2 plot
# Need : ms2 data 
# Return : MS2 plot reference vs experimental
###

doMS2 = function(data_ms2, exists_ref = T){
  if(length(data_ms2) > 0){
    if(exists_ref){
      label_peak = data_ms2$match
      data_ms2 = data_ms2$data
      index_ref = data_ms2$index[which(data_ms2$type=="ref")] ; index_exp = data_ms2$index[which(data_ms2$type=="exp")]
      
      # Color in black precursors of the drug and the isotopes detected in the chromatogram
      # Highlight matched peaks
      MyColor = ifelse(data_ms2$type=="ref",
                       ifelse(data_ms2$index != "", "#313695", "#74ADD1"), 
                       ifelse(data_ms2$index != "", "#A50026", "#F46D43"))
      MyColor[grep("precursor", data_ms2$type)] = "black"
      
      ggplot(data_ms2, aes(mz, intensity, color = type)) +
        theme_classic() + geom_hline(yintercept = 0, size=0.5, alpha = 0.8) +
        geom_segment(aes(xend=mz, yend=0),
                     color = MyColor) +
        geom_point(aes(mz, intensity), size = 1,
                   color = MyColor) +
        geom_text(data=data_ms2 %>% filter(type=="ref"), aes(label=index_ref),
                  hjust = -0.5, col="black", size=2) +
        geom_text(data=data_ms2 %>% filter(type=="exp"), aes(label=index_exp),
                  hjust = -0.5, col="black", size=2) +
        theme(plot.margin = unit(c(1,1,0,1), "lines"), 
              axis.text = element_text(size=6,colour="black"), 
              axis.title = element_text(size = 6), 
              plot.title = element_text(size = 7),
              axis.ticks= element_line(colour="black")) +
        xlab("m/z") + ylab("Relative abundance (%)") + ggtitle("MS2 reference (A) versus metabolite (B)") + 
        annotate("text", x = max(data_ms2$mz)+5, y = c(1, -1), label = c("A", "B"), size = 4)
      
    }else{
      ggplot(as.data.frame(cbind(data_ms2@mz, data_ms2@intensity)), aes(data_ms2@mz, data_ms2@intensity)) +
        theme_classic() + theme(axis.text = element_text(size=6,colour="black"), 
                                axis.title = element_text(size = 6), 
                                plot.title = element_text(size = 7),
                                axis.ticks=element_line(colour="black")) +
        geom_segment(aes(xend=data_ms2@mz, yend=0), col="#F46D43") +
        geom_point(col="#F46D43") +
        scale_y_continuous(expand=c(0,0)) +
        ggtitle("Metabolite experimental MS2") + xlab("m/z") + ylab("Relative abundance (%)")
    }
  }else{
    return(NULL)
  }
  
}
