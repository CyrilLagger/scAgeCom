# scAgeCom
Analysis of Age-related Changes in Intercellular Communication from scRNA-seq

# Folder/file description
### Top level
* analysis - main folder with detailed analysis
* src - folder with utility functions
* test - folder to test and benchmark specific parts of codes
* config.yml - config file with different sets of parameters 

### analysis
Note: it should not be necessary to run analysis 1-4 for downstream analysis. Data are already stored.
* analysis_1_Seurat_processing.R - verifies the integrity of Seurat files
    * controls the content of the files
    * checks that QC has been done properly
    * produces some statistical information and plots about each dataset/tissue/cell-type
* analysis_2_LRdb_selection.R - compares the 4 LR databases
    * summary and Venn diagram
    * selection of relevant LR-pairs
    * GO enrichment analysis on the LR pairs
* analysis_3_run_scDiffCom.R - does scDiffCom analysis on all datasets
    * run in parallel over all tissues (one dataset at a time)
    * need to be called each time for each dataset
* analysis_4_filtering - binds and filters scDiffCom results
    * depends on src_1_filtering.R
    * need to improve the cutoff exploration

### test
Note: each script can be run independently, as long as some data files are stored. Additional files can also be added if more tests are needed.
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
* test_3_scDiffCom_parameters - compare scDiffCom results for different methods/parameters
    * compares normalization, log-scale and x-sided test
    * checks the classification of each CCI for each choice of parameters
    * checks the shape of the distributions from the permutation tests
* TODO
    * do test with another file, e.g. from Calico
    * add more test related to LR_score cutoff
    * clarify how test_3 can help us to justify our choice of parameters
    * clarify the origin of the strang cases whith very low distribution means