# penda

An R package that performs personalized differential analysis of omics data.
The __penda__ package provide different methods to perform differential analysis.

## Installation

To get the current development version from github:

```
git clone https://github.com/bcm-uga/penda.git
cd penda 
R
```


## Build package

```R
install.packages("devtools")
# mixtools htmltools scales yaml lazyeval plyr rlang ggplot2 gtools caTools KernSmooth
devtools::load_all(); devtools::document(); devtools::install()
```
## Build vignettes

```R
setwd("vignettes")
rmarkdown::render("vignette_simulation.Rmd")
# source(knitr::purl(("vignette_simulation.Rmd")))

#Depends of your number of controls:
rmarkdown::render("vignette_penda.Rmd")
rmarkdown::render("vignette_penda_1ctrl.Rmd")

```

