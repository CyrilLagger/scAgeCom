####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - June 2020
##
## Test CellPhoneDB as run from the package
## scDiffCom
##
####################################################
##

library(scDiffCom)
library(Seurat)

#load a Seurat objects
#seurat_tms_test_file.rds corresponds to the Liver tissue from TMS FACS data.
seurat_test <- readRDS("test_and_comparison/seurat_tms_test_file.rds")
seurat_test <- NormalizeData(seurat_test, assay = "RNA")
seurat_test$age_group <- ifelse(seurat_test$age %in% c('1m', '3m'), 'young', 'old' )
seurat_test$cell_ontology_class <- as.character(seurat_test$cell_ontology_class)

#not sure if the following works well in this version of scDiffCom
#not really needed to run it since an example file already exists
#you also need to have CellPhoneDB installed in your current environment


# run_cpdb_from_seurat(
#   seurat_obj = seurat_test,
#   assay = "RNA",
#   slot = "data",
#   log_scale = FALSE,
#   seurat_cell_type_id = "cell_ontology_class",
#   min_cells = 5,
#   condition_id = "age_group",
#   input_dir = paste0(getwd(),"/cpdb_example"),
#   create_plots = FALSE,
#   method = 'statistical_analysis',
#   iterations = 10,
#   threshold = NULL,
#   result_precision = NULL,
#   counts_data = "gene_name",
#   output_format = NULL,
#   verbose = TRUE,
#   debug_seed = NULL,
#   threads = NULL,
#   subsampling = FALSE,
#   subsampling_log = FALSE,
#   subsampling_num_pc = NULL,
#   subsampling_num_cells = NULL,
#   return_full_dt = TRUE
# )
