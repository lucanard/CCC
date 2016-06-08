# CCC method
The package uses the principles of the Compound Characteristics Comparison method to rank candidate structures given by the in-silico fragmenter MetFrag and returns the best candidates. Multiple function are available inside the package: 
1) It is possoble to create your own CCC models, to test them and to use them in your predictions. 
2) It is possible to use the native CCC method predictions to obtain your best candidates from the in-silico fragmentation simulator MetFrag.

#### 1) Building your own standards' dataset
The dataset.building function transforms a dataset of standards into a dataset readable by the CCC_method. the original dataset must contain the compound name, the molecular formula, the monoisotopic mass, the retention time and the smiles code.
#### 1) Build your own CCC models
The function model.building gives you all the tools to build a similar model to the CCC method one using your own standards dataset. 
#### 1) Test your own CCC models
Once you create your own models you can confirm their predictivity using the model.testing function. Loading the dataset built with the dataset.building function, it will evaluate the perfomance of your models.

#### Third step:
Get your predictions: The function apply.model uses the xsaFA object from CAMERA package, clusters the m/z according to their isotopic pattern, measures their isotopic ratio and gives the most probable amount of carbon atoms in the formula. From this informations, it gives the CCC_method features.
#### Fourth step:
Get the best candidates: The functions reorder.csv and reorder.sdf use the CCC_method features to reorder the candidate structures given by MetFrag. The functions do exactly the same, but one can be used with csv files, and the other one with the sdf files (the two possible outcomes of Metfrag). We suggest to use one or the other format depending on the fact you want to keep working on R or not. Indeed, the reorder.csv function gives a data.frame that can be saved as csv file. On the other hand, the reorder.sdf gives a SDFset file readable from the package ChemmineR. 

![Sample image](https://github.com/lucanard/CCC/blob/master/CCC%20flowchart%20-%20Standard.png "CCC workflow")

# Installation
First install dependencies from CRAN

```R
install.packages("devtools")
install.packages("stringr")
install.packages("rcdk")
install.packages("rinchi")
install.packages("ppls")
install.packages("glmnet")
```
Then install dependencies from bioconductor

```R
source("http://bioconductor.org/biocLite.R")
biocLite()
library(BiocInstaller)
biocLite("ChemmineR")
biocLite("fmcsR")
```
Last, install metfRag and the CCC_method from github

```R
library(devtools)
install_github("c-ruttkies/MetFragR/metfRag")
install_github("lucanard/CCC", subdir="CCC_method")
```
