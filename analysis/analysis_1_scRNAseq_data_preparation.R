####################################################
##
## Project: scAgeCom
##
## Last update - April 2021
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## scRNA-seq data preparation
##
####################################################
##

## Dataset Overview ####

# Tabula Muris Senis scRNA-seq data have been retrieved from Amazon S3 
# following the link given here https://github.com/czbiohub/tabula-muris-senis

# Calico2019 scRNA-seq data have been downloaded from here
# https://mca.research.calicolabs.com/data

# There were both stored as h5ad scanpy files, that we converted to 
# Seurat objects and stored on your lab server.

## Libraries ####

library(Seurat)
library(data.table)
## Dataset paths on Server ####

server_path <- "/home/nis/lcyril/work/lcyril_data/scRNA_seq/seurat_processed/"

dataset_path <- c(
  tms_facs = paste0(server_path, "seurat_final_tms_facs.rds"),
  tms_droplet = paste0(server_path, "seurat_final_tms_droplet.rds"),
  calico_kidney = paste0(server_path, "seurat_final_calico_kidney.rds"),
  calico_lung = paste0(server_path, "seurat_final_calico_lung.rds"),
  calico_spleen = paste0(server_path, "seurat_final_calico_spleen.rds")
)

## Load datasets ####

# requires several GB of RAM

seurat_objects <- list(
  tms_facs = readRDS(dataset_path[["tms_facs"]]),
  tms_droplet = readRDS(dataset_path[["tms_droplet"]]),
  calico_kidney = readRDS(dataset_path[["calico_kidney"]]),
  calico_lung = readRDS(dataset_path[["calico_lung"]]),
  calico_spleen = readRDS(dataset_path[["calico_spleen"]])
)

## Check dataset content ####

## Prepare data for scAgeComShiny ####

## Prepare figures for manuscript ####


