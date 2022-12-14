# penda

An R package that performs personalized differential analysis of omics data.
The __penda__ package provide different methods to perform differential analysis.

*PenDA, a rank-based method for personalized differential analysis: Application to lung cancer*
Magali Richard,Clémentine Decamps, Florent Chuffart, Elisabeth Brambilla, Sophie Rousseaux, Saadi Khochbin, Daniel Jost.
**PLoS Comput Biol**, 2020 May 11;16(5):e1007869. 
https://doi.org/10.1371/journal.pcbi.1007869.


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

