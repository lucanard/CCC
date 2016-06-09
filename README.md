# CCC method
The package uses the principles of the Compound Characteristics Comparison method to rank candidate structures given by the in-silico fragmenter MetFrag and returns the best candidates. Multiple function are available inside the package: 
1) It is possoble to create your own CCC models, to test them and to use them in your predictions. 
2) It is possible to use the native CCC method predictions to obtain your best candidates from the in-silico fragmentation simulator MetFrag.

### 1) Your own CCC model
#### Building your own standards' dataset
The dataset.building function transforms a dataset of standards into a dataset readable by the CCC_method. the original dataset must contain the compound name, the molecular formula, the monoisotopic mass, the retention time and the smiles code.
#### Build your own CCC models
The function model.building gives you all the tools to build a similar model to the CCC method one using your own standards dataset. 
#### Test your own CCC models
Once you create your own models you can confirm their predictivity using the model.testing function. Loading the dataset built with the dataset.building function, it will evaluate the perfomance of your models.

### 2) Rerank the structures proposed by MetFrag
#### Obtain the CCC features
The function CCC.peakTable and CCC.xsaFA use the peaktable/xsaFA objects, to cluster the features according to their isotopic pattern, measuring their isotopic ratio and they give the most probable amount of carbon atoms in the formula and the CCC_method features.
#### Obtain the related pseudo-spectra. It returns a clusterTable
The function ps_spec gives the pseudospectra related to the clusterTable clusters. 
#### Rerank the metFrag candidates 
Get the best candidates: The functions rerank2010, rerank2016 use the clustertable features to rerank the candidate structures given by MetFrag2010 or MetfragBeta (2016). The functions do exactly the same, but one can be used with csv and sdf files on the 2010 metFrag files, while the other one only with the csv files of the MetFragBeta version (2016). 

![Sample image](https://github.com/lucanard/CCC/blob/master/CCC%20flowchart%20-%20Standard%20(1).png?raw=true "CCC workflow")

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
