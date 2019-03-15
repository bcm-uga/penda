# penda

An R package that performs personalized differential analysis of omics data.
The __penda__ package provide different methods to perform differential analysis.

## Installation

To get the current development version from github:

```R
install.packages("devtools")
devtools::install_github("bcm-uga/penda")

## Build vignettes

```R
setwd("vignette")
rmarkdown::render("vignette_simulation.Rmd")
# source(knitr::purl(("vignette_simulation.Rmd")))