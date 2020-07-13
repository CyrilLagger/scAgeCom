# scAgeCom
Analysis of Age-related Changes in Intercellular Communication from scRNA-seq

# Folder/file description
### Top level
* analysis - main folder with detailed analysis
* src - folder with utility functions
* test - folder to test and benchmark specific parts of codes
* config.yml - config file with different sets of parameters 

### analysis

### test
Note: each script can be run independently, as long as some data files are stored.
* test_1_scDiffcom.R - checks that scDiffCom returns correct values
    * requires a Seurat file for testing (e.g: "../data_scAgeCom/data_seurat_example.rds")
    * shows typical usage of scDiffCom
    * checks scores and p-values (from independent functions)
    * create plots to compare distributions from permutations
* test_2_compare_to_cpdb.R - check that scDiffCom returns similar values compared to CellPhoneDB
    * requires preprocessed CPDB output
    * shows (deprecated) usage of run_cpdb_from_seurat
    * shows that LR scores are exactly the same between the methods
    * shows that LR specificity p-values are similar (up to randomness)
    * shows high similarity in classifying significant/non-significant CCIs