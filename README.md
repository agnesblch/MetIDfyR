![metidfyr_300px](https://user-images.githubusercontent.com/44233550/202466223-ca5d2c08-2e8c-4452-afda-e69e9bdf1480.png)


# `MetIDfyR`
`MetIDfyR` is an open-source, cross-platform and versatile R script to predict and detect metabolites in mass spectrometry data (mzML) based on the raw formula of the drug of interest.

## How to get `MetIDfyR` ?
To use `MetIDfyR` on your computer you need to clone the project from GitHub or by using Rstudio version control (preferred). You can also download it from GitHub or by using the following commands in R:
```
download.file(url = "https://github.com/agnesblch/MetIDfyR/archive/master.zip",
              destfile = "MetIDfyR-master.zip")

unzip('MetIDfyR-master.zip')
setwd('MetIDfyR-master/')
```

## Need
- a configuration file with the parameter (see TEMPLATE_config.R). 
- a mzML file containing the sample informations
- a tsv file containing the parent drug informations (see TEMPLATE_start_mlc.tsv)

Launch in the software directory in command line (or add the directory to your PATH):
```
Rscript MetIDfyR.R -i path-to-input-file -o output-directory -c config-file
```
Or from Rstudio:
```
system("Rscript MetIDfyR.R -i path-to-input-file -o output-directory -c config-file")
```

## Example
To test `MetIDfyR`, you can perform analysis using provided LGD-4033, Cocaine and Diclofenac datasets. 

1. Clone the `MetIDfyR` repository or download it from GitHub.
2. Create a new R project in the `MetIDfyR` directory.
3. Unzip the example datasets using the following command in R console:
```
unzip("input/LGD_DIA_peak-picked.mzML.zip",exdir = "input")
unzip("input/U_H_COCA_peak-picked.mzML.zip",exdir = "input")
unzip("input/S9_Diclofenac_PhI.zip",exdir = "input")
unzip("input/S9_Diclofenac_PhII.zip",exdir = "input")
unzip("input/S9_Diclofenac_std.zip",exdir = "input")
```
4. Launch `MetIDfyR` from R console using:
```
system("Rscript MetIDfyR.R -i input/lgd_DIA_peak-picked.tsv -o LGD4033_results -c input/config_LGD.R")
system("Rscript MetIDfyR.R -i input/cocaine_peak-picked.tsv -o Cocaine_results -c input/config_cocaine.R")
system("Rscript MetIDfyR.R -i input/diclofenac-picked.tsv -o Diclofenac_results -c input/config_diclofenac.R")
```
5. When the run is done, open the file "ui.R" in MetApp directory and click "run app" to launch MetApp.
6. In the "Visualization" panel open the output folder containing the file "out_....tsv".
7. You can generate a PDF report containing the selected metabolites.

## Dependencies
MetIDfyR was developed and tested on <a href="https://www.r-project.org/" title="More about R">R</a> version 3.6.1. It works normally up to R version 4.0. 
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

You can find this application in the repository: https://github.com/GIELCH/MetApp

## Citation
>Delcourt V, Barnabé A, Loup B, Garcia P, André F, Chabot B, Trévisiol S, Moulard Y, Popot M-A & Bailly-Chouriberry L -
> *`MetIDfyR`, an Open-Source R Package to Decipher Small-Molecule Drugs Metabolism Through High Resolution Mass Spectrometry*. 2020. AnalChem. https://doi.org/10.1021/acs.analchem.0c02281
