# CCC
CCC method
The package uses the principles of the Compound Characteristics Comparison method to reorder candidate structures given by the in-silico fragmenter MetFrag. To perform the analysis with your own data, you need to follow the instructions in the order reported here.
# First step: 
The dataset.building function transforms a dataset of standards containing the compound name, molecular formula, monoisotopic mass, retention time and smiles code into a dataset readable by the CCC_method. NB: If you use the dataset given by the package, you can jump this step.
# Second step:
Test the dataset: you need to test the dataset you just build to see if it is appropriate for your purpose. The model.testing function performs standard test and gives the results in a easy readable data.frame. Jump this step if you use the dataset given in the package
#Third step

