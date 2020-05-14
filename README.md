# MetIDfyR
This tool allows prediction and detection of metabolites in a sample (mzML file). The pipeline returns as output spectrum for each metabolites if there is a signal in positive and/or negative analysis.

## Need : 
- a configuration file with the parameter (see TEMPLATE_config.R). 
- a mzML file containing the sample informations
- a tsv file containing the parent drug informations (see TEMPLATE_start_mlc.tsv)

Launch in the software directory in command line (or add the directory to your PATH) :
```
Rscript MetIDfyR.R -i path-to-input-file -o output-directory -c config-file
```
Or from Rstudio :
```
system("Rscript MetIDfyR.R -i path-to-input-file -o output-directory -c config-file")
```

## Example :
To test MetIDfyR, you can perform analyses using provided LGD-4033 and Cocaine datasets. 

You can unzip the datasets from Rstudio using :
```
unzip("input/LGD_DIA_peak-picked.mzML.zip",exdir = "input")
unzip("input/U_H_COCA_peak-picked.mzML.zip",exdir = "input")
```

Then you can launch MetIDfyR from Rstudio using :
```
system("Rscript MetIDfyR.R -i input/lgd_DIA_peak-picked.tsv -o LGD4033_results -c input/config_LGD.R")
system("Rscript MetIDfyR.R -i input/cocaine_peak-picked.tsv -o Cocaine_results -c input/config_cocaine.R")
```

## Dependencies : 
- pacman
- BiocManager
- optparse
- ggplot2
- dplyr / plyr
- foreach / doParallel
- stringr
- MSnbase
- Rdisop
- readr
- tibble
- ggpubr
- ggrepel
- arrangements
- snow
- RColorBrewer

## MetApp
Shiny application to display resulting metabolites figures.
This app allow to save PDF report for selected metabolites thanks to the packages "rmarkdown" and "rsvg". 
The SVG figures are printed in the report using the R package "magick" (https://cran.rstudio.com/bin/windows/contrib/3.6/magick_2.3.zip). 

