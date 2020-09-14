# `MetIDfyR`
`MetIDfyR` is an open-source, cross-platform and versatile R script to predict and detect metabolites in mass spectrometry data (mzML) based on the raw formula of the drug of interest.

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
To test `MetIDfyR`, you can perform analyses using provided LGD-4033 and Cocaine datasets. 

1. Clone the `MetIDfyR` repository.
2. Create a new R project in the `MetIDfyR` directory.
3. Unzip the example datasets using the following command in R console :
```
unzip("input/LGD_DIA_peak-picked.mzML.zip",exdir = "input")
unzip("input/U_H_COCA_peak-picked.mzML.zip",exdir = "input")
unzip("input/S9_Diclofenac_PhI.zip",exdir = "input")
unzip("input/S9_Diclofenac_PhII.zip",exdir = "input")
unzip("input/S9_Diclofenac_std.zip",exdir = "input")
```
4. Launch `MetIDfyR` from R console using :
```
system("Rscript MetIDfyR.R -i input/lgd_DIA_peak-picked.tsv -o LGD4033_results -c input/config_LGD.R")
system("Rscript MetIDfyR.R -i input/cocaine_peak-picked.tsv -o Cocaine_results -c input/config_cocaine.R")
system("Rscript MetIDfyR.R -i input/diclofenac-picked.tsv -o Diclofenac_results -c input/config_diclofenac.R")
```
5. When the run is done, open the file "ui.R" in MetApp directory and click "run app" to launch MetApp.
6. Visualize and generate a PDF report using the "Visualization" tab in the app.

## Dependencies : 
MetIDfyR was developed and tested on <a href="https://www.r-project.org/" title="More about R">R</a> version 3.6.1. 
Pacman package is used to install all dependencies needed to run `MetIDfyR`.
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
This app saves PDF report for selected metabolites thanks to the packages "rmarkdown" and "rsvg". 
The SVG figures are printed in the report using the R package "magick" (https://cran.rstudio.com/bin/windows/contrib/3.6/magick_2.3.zip). 

## Citation
>Delcourt V, Barnabé A, Loup B, Garcia P, André F, Chabot B, Trévisiol S, Moulard Y, Popot M-A & Bailly-Chouriberry L -
> *`MetIDfyR`, an Open-Source R Package to Decipher Small-Molecule Drugs Metabolism Through High Resolution Mass Spectrometry*. 2020. AnalChem. https://doi.org/10.1021/acs.analchem.0c02281
