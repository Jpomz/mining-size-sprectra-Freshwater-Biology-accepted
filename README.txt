README

This R project contains all of the code necessary to reproduce the results and figures in:
Pomeranz, Warburton and Harding. Anthropogenic mining alters macroinvertebrate size spectra in streams. Freshwater Biology. in press. 

Data are available at: https://doi.org/10.5061/dryad.v6g985s

In order to conduct the analysis, download the data from Data Dryad, and place the two .csv files in the "data" file. 

All results and figures will be saved in the "results" file. 

The scripts should be run in numerical order (e.g. 1_calculate...R, 2_size...R etc.)

Each script lists the libraries necessary for that script to properly execute. You must first install the packages before you can load them. 
e.g.
> install.packages("plyr")
> library(plyr)
You only need to install the package once, but need to load it each time you start R. 
Some of these packages have naming discrepencies (e.g. plyr and dplyr both have functions with the same name), therefore it is recommended to restart your R session after running each script. (default shortcut: Ctrl + Shift + F10)

Any questions can be referred to the corresponding author, Justin Pomeranz, at jfpomeranz@gmail.com