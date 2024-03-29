
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Reproduction-Material-Spatial-Monte-Carlo-Experiments

This repository provides all materials for reproducing the results of
“Spatial Regression Models: A Systematic Comparison of Different Model
Specifications using Monte Carlo Experiments”.

## Requirements

The codes needs the following folders: “01\_Script”, “02\_Data”,
“03\_Output”. All R Scripts are required in folder “01\_Script”.

The following packages are necessary for reproduction of main results:

``` r
install.packages("ggmap")
install.packages("RgoogleMaps") 
install.packages("foreign")
install.packages("sp") 
install.packages("GISTools")
install.packages("sp") 
install.packages("spdep")
install.packages("texreg") 
install.packages("xtable")
install.packages("rgdal") 
install.packages("raster")
install.packages("extrafont")
```

## Explanations

All R scripts named “01~” contain the functions required for the
simulations. To reproduce the results of the final paper, you need to
run the R script “02\_Monte Carlo Simulation Spatial\_Main Analyses”.
Results of the Supplement are produced by running “03\_Monte Carlo
Simulation Spatial\_Supplementary Analyses”. You only need to change the
superordinate directory “dir” (line 12), containing the three
sub-folders. Generated data will be saved in “02\_Data”, output will be
generated in the “03\_Output” directory.

Note that reproducing the results with R=1000 simulation runs may take a
lot of time on a conventional computer. Thus, all simulation programs
are also available as “~\_hpc” version for usage on a high performance
cluster with the package “doParallel”. Depending on your system, you
need to set up additional batch files to pass the code to a high
performance cluster. However, it may decrease simulation time
dramatically.

## System and version information

Platform: Windows 7 (x86\_64-w64-mingw32)

Version: R version 3.5.0
