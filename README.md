# penda

An R package that performs personalized differential analysis of omics data.
The __penda__ package provide different methods to perform differential analysis.

## Installation

To get the current development version from github:
git clone https://github.com/bcm-uga/penda.git
cd penda 
R

## Build vignettes

```R
setwd("vignette")
rmarkdown::render("vignette_simulation.Rmd")
# source(knitr::purl(("vignette_simulation.Rmd")))

