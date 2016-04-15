# CCC method
The package uses the principles of the Compound Characteristics Comparison method to reorder candidate structures given by the in-silico fragmenter MetFrag and return the best candidates. To perform the analysis with your own data, you need to follow the instructions in the order reported here.
#### First step: 
Build the standard dataset: The dataset.building function transforms a dataset of standards into a dataset readable by the CCC_method. the original dataset must contain the compound name, the molecular formula, the monoisotopic mass, the retention time and the smiles code. NB: If you use the dataset given by the package, you can jump this step.
#### Second step:
Test the dataset: you need to test the your dataset to see if it is appropriate for the method. The model.testing function performs standard test and gives the results in a easy readable data.frame. Jump this step if you use the dataset given in the package (which is already validated)
#### Third step:
Get your predictions: The function apply.model uses the xsaFA object from CAMERA package, clusters the m/z according to their isotopic pattern, measures their isotopic ratio and gives the most probable amount of carbon atoms in the formula. From this informations, it gives the CCC_method features.
#### Forth step:
Get the best candidates: The functions reorder.csv and reorder.sdf use the CCC_method features to reorder the candidate structures given by MetFrag. The functions do exactly the same, but one can be used with csv files, and the other one with the sdf files (the two possible outcomes of Metfrag). We suggest to use one or the other format depending on the fact you want to keep working on R or not. Indeed, the reorder.csv function gives a data.frame that can be saved as csv file. On the other hand, the reorder.sdf gives a SDFset file readable from the package ChemmineR. 

# Installation
First install all dependencies

```R
install.packages("devtools")
install.packages("stringr")
install.packages("rcdk")
install.packages("pls")
install.packages("ppls")
install.package("glmnet")
source("http://bioconductor.org/biocLite.R")
biocLite()
library(BiocInstaller)
biocLite("xcms")
biocLite("CAMERA")
biocLite("ChemmineR")
biocLite("fmcsR")
library(devtools)
install_github("lucanard/CCC")
```
