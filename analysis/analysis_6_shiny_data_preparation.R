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
## collect and prepare all results for shiny
##
####################################################
##

## Libraries ####

library(scDiffCom)
library(data.table)

## save all previous results ####

scAgeCom_shiny_data <- c(
  readRDS(
    "../data_scAgeCom/analysis/outputs_data/data_2_LRI_data_preparation.rds"
  ),
  readRDS(
    "../data_scAgeCom/analysis/outputs_data/data_4_tissue_specific_results.rds"
  ),
  readRDS(
    "../data_scAgeCom/analysis/outputs_data/data_5_tissue_shared_results.rds"
  )
)

saveRDS(
  scAgeCom_shiny_data,
  "../data_scAgeCom/analysis/outputs_data/scAgeCom_shiny_data.rds"
)

