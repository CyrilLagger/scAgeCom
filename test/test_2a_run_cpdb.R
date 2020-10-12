####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - September 2020
##
## Script that show how to run CELLPHONEDB from 
## the scDiffCom package.
##
####################################################
##

## Libraries ####

library(Seurat)
library(scDiffCom)
library(data.table)
library(ggplot2)

## Data path ####
dir_data_test <- "../data_scAgeCom/test/"

## Load a Seurat objects ####

#here the file corresponds to the Liver tissue from TMS FACS data.
seurat_t2 <- readRDS(paste0(dir_data_test, "data_seurat_example_facs_liver.rds"))
seurat_t2 <- NormalizeData(seurat_t2, assay = "RNA")
seurat_t2$age_group <- ifelse(seurat_t2$age %in% c('1m', '3m'), 'young', 'old' )
seurat_t2$cell_ontology_class <- as.character(seurat_t1$cell_ontology_class)

## Run CellPhoneDB ####

#Note: the function scDiffCom::run_cpdb_from_seurat is a R wrapper
#      around the original CPDB python package. It only works on
#      Linux and if CPDB is installed in your current python environment.

#typical usage (for the parameters see the CPDB manual: https://github.com/Teichlab/cellphonedb)
#the output is the typical CPDB output with various files. We also create a single data.table that 
#summarize the main results.
run_cpdb_from_seurat(
  seurat_obj = seurat_t2,
  assay = "RNA",
  slot = "data",
  log_scale = FALSE,
  seurat_cell_type_id = "cell_ontology_class",
  input_species = "mouse",
  min_cells = 5,
  condition_id = "age_group",
  input_dir = paste0(dir_data_test,"cpdb_test_results"),
  create_plots = FALSE,
  method = 'statistical_analysis',
  iterations = 1000,
  threshold = NULL,
  result_precision = NULL,
  counts_data = "gene_name",
  output_format = NULL,
  verbose = TRUE,
  debug_seed = NULL,
  threads = NULL,
  subsampling = FALSE,
  subsampling_log = FALSE,
  subsampling_num_pc = NULL,
  subsampling_num_cells = NULL,
  return_full_dt = TRUE
)


