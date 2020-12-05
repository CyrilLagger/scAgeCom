##
## Project: scAgeCom
##
## Last update - December 2020
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## Process scDiffCom results and check filtering
## parameters.
##
####################################################
##

## Libraries ####
library(scDiffCom)
library(data.table)
library(ggplot2)
#library(pbapply)

## Specify the directory with all scDiffCom results ####
dir_results <- "../data_scAgeCom/scdiffcom_results_30_11_2020"
dir_data_analysis <- "../data_scAgeCom/analysis/"

## Retrieve all the results for each tissue and dataset ####
RESULT_PATHS <- list.dirs(dir_results, recursive = FALSE)
RESULT_PATHS

RESULT_PATHS <- RESULT_PATHS[c(1, 4, 8)]

DATASETS <- lapply(
  RESULT_PATHS, 
  function(path) {
    tissues <- gsub(".*scdiffcom_(.+)\\.rds.*", "\\1", list.files(path))
    temp <- lapply(
      X = tissues,
      FUN = function(
        tiss
      ) {
        res <- readRDS(paste0(path, "/scdiffcom_", tiss, ".rds"))
      }
    )
    names(temp) <- tissues
    return(temp)
  }
)

## Set the names manually, be careful to check what you are doing ####
RESULT_PATHS
#names(DATASETS) <- c("calico")
names(DATASETS) <- c("calico", "droplet_mixed", "facs_mixed")
#names(DATASETS) <- c("calico", "droplet_female", "droplet_male", "droplet_mixed", "droplet_sex",
#                     "facs_female", "facs_male", "facs_mixed", "facs_sex")

DATASETS_COMBINED <- lapply(
  seq_along(DATASETS),
  function(i) {
    scDiffCom:::Combine_scDiffCom(
      l = DATASETS[[i]],
      object_name = names(DATASETS)[[i]],
      verbose = TRUE
    )
  }
)
names(DATASETS_COMBINED) <- names(DATASETS)

#saveRDS(DATASETS_COMBINED, "../../../../../temp_agecom.rds")
#DATASETS_COMBINED <- readRDS("../../../../../temp_agecom.rds")

## See the parameters ####
lapply(
  DATASETS_COMBINED,
  function(i) {
    parameters(i)
  }
)

## Inspect the impact of threshold logfc and threshold score ####

test <- lapply(
  DATASETS_COMBINED,
  function(data_comb) {
    temp_seq <- seq(1.1, 2, 0.1)
    temp_list <- lapply(
      temp_seq,
      function(thr) {
        temp_obj <- FilterCCI(
          object = data_comb,
          new_threshold_logfc = log(thr),
          skip_ora = TRUE
        )
        get_cci_detected(temp_obj)[, {
          temp = .N
          .SD[,.(pct=.N/temp*100),by=REGULATION_SIMPLE]
        }, by = ID ]
      }
    )
    names(temp_list) <- as.character(temp_seq)
    temp_dt <- rbindlist(
    temp_list,
    use.names = TRUE,
    idcol = "logfc_threshold"
    )
    temp_dt$logfc_threshold <- as.numeric(temp_dt$logfc_threshold)
    return(temp_dt)
  }
)

ggplot(test$calico[REGULATION_SIMPLE == "FLAT"], aes(x = logfc_threshold, y = pct, color = ID)) + geom_point()

ggplot(test$calico[REGULATION_SIMPLE == "FLAT"], aes(x = logfc_threshold, y = pct, group = logfc_threshold )) + geom_boxplot()
ggplot(test$droplet_mixed[REGULATION_SIMPLE == "FLAT"], aes(x = logfc_threshold, y = pct, group = logfc_threshold )) + geom_boxplot()
ggplot(test$facs_mixed[REGULATION_SIMPLE == "FLAT"], aes(x = logfc_threshold, y = pct, group = logfc_threshold )) + geom_boxplot()

## Test with a logfc of log(1.5) ####
DATASETS_COMBINED_log15 <- lapply(
  DATASETS_COMBINED,
  function(i) {
    FilterCCI(
      object = i,
      new_threshold_logfc = log(1.5),
      skip_ora = FALSE,
      verbose = TRUE
    )
  }
)

saveRDS(DATASETS_COMBINED_log15, "../../../../../test15_combined.rds")

DATASETS_COMBINED_log15 <- readRDS("../../../../../test15_combined.rds")

BuildNetwork(
  DATASETS_COMBINED_log15$droplet_mixed,
  network_type = "bipartite",
  subobject_name = "Kidney"
)

test <- lapply(
  DATASETS_COMBINED_log15,
  function(i) {
    i@cci_raw <- list()
    return(i)
  }
)

saveRDS(test, "../../../../../test_light.rds")

test <- readRDS("../../../../../test_light.rds")

test <- lapply(
  test,
  function(i) {
    RunORA(
      object = i,
      stringent = "stringent",
      stringent_logfc_threshold = log(2),
      verbose = TRUE
    )
  }
)

test <- lapply(
  test,
  function(i) {
    RunORA(
      object = i,
      stringent = "stringent",
      stringent_logfc_threshold = log(2),
      verbose = TRUE,
      global = FALSE
    )
  }
)
 

generate_interactive_network(
  DATASETS_COMBINED_log15,
  network_type = "bipartite"
)

test <- DATASETS_COMBINED_log15$facs_mixed@ora_combined_default$LR_GENES
test2 <- DATASETS_COMBINED$facs_mixed@ora_combined_default$LR_GENES

DATASETS_COMBINED$facs_mixed@ora_default$LR_GENES[OR_UP > 1 & BH_P_VALUE_UP <= 0.05][, .N, by = c("VALUE")][order(-N)]

test3 <- DATASETS_COMBINED_log15$facs_mixed@ora_default$LR_GENES[OR_DOWN > 1 & BH_P_VALUE_DOWN <= 0.05][, .N, by = c("VALUE")][order(-N)]

test4 <- DATASETS_COMBINED_log15$facs_mixed@ora_default$GO_TERMS[OR_DOWN > 1 & BH_P_VALUE_DOWN <= 0.05][, .N, by = c("VALUE")][order(-N)]
test5 <- DATASETS_COMBINED_log15$facs_mixed@ora_default$GO_TERMS[OR_UP > 1 & BH_P_VALUE_UP <= 0.05][, .N, by = c("VALUE")][order(-N)]

test6 <- DATASETS_COMBINED_log15$droplet_mixed@ora_default$GO_TERMS[OR_UP > 1 & BH_P_VALUE_UP <= 0.05][, .N, by = c("VALUE")][order(-N)]

test7 <- DATASETS_COMBINED_log15$facs_mixed@ora_default$GO_TERMS[OR_DOWN > 1 & BH_P_VALUE_DOWN <= 0.05][, .N, by = c("VALUE")][order(-N)]
test8 <- DATASETS_COMBINED_log15$droplet_mixed@ora_default$GO_TERMS[OR_DOWN > 1 & BH_P_VALUE_DOWN <= 0.05][, .N, by = c("VALUE")][order(-N)]
test9 <- DATASETS_COMBINED_log15$calico@ora_default$GO_TERMS[OR_DOWN > 1 & BH_P_VALUE_DOWN <= 0.05][, .N, by = c("VALUE")][order(-N)]

intersect(test7[N>3]$VALUE, test5[N>3]$VALUE)

## See how to deal with sex ####


## See how to deal with age ####



## Analyse the influence of the score and logfc thresholds  ####

#function to see the impact of these two parameters on the distribution of REGULATION
get_regulation_behaviour <- function(
  object,
  logfc_cuts,
  score_cuts
) {
  rbindlist(
    lapply(
      logfc_cuts,
      function(logfc_cut) {
        rbindlist(
          lapply(
            score_cuts,
            function(score_cut) {
              dt <- get_cci_table_filtered(
                run_filtering_and_ora(
                  object = object,
                  new_threshold_quantile_score = score_cut,
                  new_threshold_logfc = logfc_cut,
                  skip_ora = TRUE,
                  verbose = FALSE
                )
              )
              res <- dt[, .N, by = c("REGULATION")][, pct := N/sum(N)*100]
            }
          ),
          use.names = TRUE,
          idcol = "score_cut"
        )
      }
    ),
    use.names = TRUE,
    idcol = "logfc_cut"
  )
}

get_regulation_behaviour_bind <- function(
  dataset
) {
  temp <- lapply(
    seq_along(dataset),
    function(i) {
      print(names(dataset)[i])
      temp2 <- lapply(
        seq_along(dataset[[i]]),
        function(j) {
          print(names(dataset[[i]])[j])
          get_regulation_behaviour(
            object = dataset[[i]][[j]],
            logfc_cuts = logfc_cuts,
            score_cuts = score_cuts
          )
        }
      )
      names(temp2) <-names(dataset[[i]])
      rbindlist(
        temp2,
        use.names = TRUE,
        idcol = "Tissue"
      )
    }
  )
  names(temp) <- names(dataset)
  res <- rbindlist(
    temp,
    use.names = TRUE,
    idcol = "Dataset"
  )
  res[, Dataset_Tissue := paste(Dataset, Tissue, sep = "_")]
  return(res)
}

#define a series of cutoffs

logfc_cuts <- c(1.1, 1.2, 1.3, 1.4, 1.5)
names(logfc_cuts) <- logfc_cuts
logfc_cuts <- log(logfc_cuts)

score_cuts <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
names(score_cuts) <- score_cuts

#computation for each dataset and tissue (can take more than one hour)
regulation_behaviour_mixed <- get_regulation_behaviour_bind(
  dataset = DATASETS_mixed
)
regulation_behaviour_sex <- get_regulation_behaviour_bind(
  dataset = DATASETS_sex
)

#save the results

saveRDS(regulation_behaviour_sex, "../data_scAgeCom/analysis/a4_data_regulation_behaviour_sex_2datasets.rds")
saveRDS(regulation_behaviour_mixed, "../data_scAgeCom/analysis/a4_data_regulation_behaviour_mixed_3datasets.rds")

## Plot the behaviour for different tissues and dataset ####

#box plots over all tissues
ggplot(regulation_behaviour_mixed[REGULATION == "FLAT" & score_cut == 0.25],
       aes(x = logfc_cut, y = pct, color = Dataset)) +
  geom_boxplot()

ggplot(regulation_behaviour_sex[REGULATION == "FLAT" & score_cut == 0.25],
       aes(x = logfc_cut, y = pct, color = Dataset)) +
  geom_boxplot()

#box plots over common tissues
ggplot(regulation_behaviour_mixed[REGULATION == "FLAT" & score_cut == 0.25 & Tissue %in% c("Kidney", "Lung", "Spleen")],
       aes(x = logfc_cut, y = pct, color = Dataset)) +
  geom_boxplot()

#line plots over all tissues
ggplot(regulation_behaviour_mixed[Dataset == "facs_mixed" & REGULATION == "FLAT" & score_cut == 0.25],
       aes(x = logfc_cut, y = pct, color = Tissue, group = Tissue)) +
  geom_line()

ggplot(regulation_behaviour[Dataset == "droplet_mixed" & REGULATION == "FLAT" & score_cut == 0.25],
       aes(x = logfc_cut, y = pct, color = Tissue, group = Tissue)) +
  geom_line()

ggplot(regulation_behaviour[Dataset == "calico" & REGULATION == "FLAT" & score_cut == 0.25],
       aes(x = logfc_cut, y = pct, color = Tissue, group = Tissue)) +
  geom_line()

ggplot(regulation_behaviour_sex[Dataset == "facs_sex" & REGULATION == "FLAT" & score_cut == 0.25],
       aes(x = logfc_cut, y = pct, color = Tissue, group = Tissue)) +
  geom_line()

#line plots over common tissues
ggplot(regulation_behaviour[REGULATION == "FLAT" & score_cut == 0.25 & Tissue %in% c("Kidney", "Lung", "Spleen")],
       aes(x = logfc_cut, y = pct, color = Dataset_Tissue, group = Dataset_Tissue)) +
  geom_line()

ggplot(regulation_behaviour[REGULATION == "FLAT" & score_cut == 0.25 & Tissue %in% c("Kidney")],
       aes(x = logfc_cut, y = pct, color = Dataset_Tissue, group = Dataset_Tissue)) +
  geom_line()

ggplot(regulation_behaviour[REGULATION == "FLAT" & score_cut == 0.25 & Tissue %in% c("Spleen")],
       aes(x = logfc_cut, y = pct, color = Dataset_Tissue, group = Dataset_Tissue)) +
  geom_line()

ggplot(regulation_behaviour[REGULATION == "FLAT" & score_cut == 0.25 & Tissue %in% c("Lung")],
       aes(x = logfc_cut, y = pct, color = Dataset_Tissue, group = Dataset_Tissue)) +
  geom_line()


## Compare the logfc ####
full_cci_table_filtered <- rbindlist(
  lapply(
    DATASETS_mixed,
    function(dataset) {
      rbindlist(
        lapply(
          dataset,
          function(tissue) {
            get_cci_table_filtered(tissue)
          }
        ),
        use.names = TRUE,
        idcol = "Tissue",
        fill = TRUE
      )
    }
  ),
  use.names = TRUE,
  idcol = "Dataset"
)

full_cci_table_filtered_sex <- rbindlist(
  lapply(
    DATASETS_sex,
    function(dataset) {
      rbindlist(
        lapply(
          dataset,
          function(tissue) {
            get_cci_table_filtered(tissue)
          }
        ),
        use.names = TRUE,
        idcol = "Tissue",
        fill = TRUE
      )
    }
  ),
  use.names = TRUE,
  idcol = "Dataset"
)

full_cci_table_bind <- rbindlist(
  list(age = full_cci_table_filtered[, c("Dataset", "Tissue", "LOGFC")], sex = full_cci_table_filtered_sex[, c("Dataset", "Tissue", "LOGFC")]),
  use.names = TRUE,
  idcol = TRUE
)

ggplot(full_cci_table_filtered, aes(x = LOGFC, color = Dataset, ..scaled..)) +
  geom_density(position = "identity", alpha = 0.2)
ggplot(full_cci_table_filtered, aes(x = LOGFC, color = Dataset)) +
  geom_density(position = "identity", alpha = 0.2) + xlim(c(-5, 5))
ggplot(full_cci_table_filtered[Tissue %in% c("Kidney", "Lung", "Spleen")], aes(x = LOGFC, color = Dataset)) +
  geom_density(position = "identity", alpha = 0.2) + xlim(c(-5, 5))

ggplot(full_cci_table_filtered_sex, aes(x = LOGFC, color = Dataset)) +
  geom_density(position = "identity", alpha = 0.2)

ggplot(full_cci_table_bind, aes(x = LOGFC, color = Dataset)) +
  geom_density(position = "identity", alpha = 0.2)


## We keep the initial parameters log(1.2) and 0.25 and create a new object without the raw data ####

DATASETS_light <- lapply(
  DATASETS_mixed,
  function(i) {
    lapply(
      i,
      function(tiss) {
        tiss <- set_cci_table_raw(tiss, list()) 
        return(tiss)
      }
    )
  }
)

#saveRDS(DATASETS_light, "../data_scAgeCom/analysis/a4_data_results_all_cases.rds")
saveRDS(DATASETS_light, "../data_scAgeCom/analysis/a4_data_results_3cases.rds")

## Read light dataset ####
DATASETS_light <- readRDS("../data_scAgeCom/analysis/a4_data_results_all_cases.rds")

## Look at ligand vs receptor contributions to scores and logfc ####

coherence_plot <- function(
  dataset,
  axis_lim
) {
  dt_filt <- rbindlist(
    lapply(
      dataset,
      function(i) {
        scDiffCom:::get_cci_table_filtered(i)
      }
    ),
    fill = TRUE
  )
  dt_filt[, LOG2FC_L := log2(
    pmin(
      L1_EXPRESSION_OLD,
      L2_EXPRESSION_OLD,
      na.rm = TRUE
      )
    /
      pmin(
        L1_EXPRESSION_YOUNG,
        L2_EXPRESSION_YOUNG,
        na.rm = TRUE
        )
    ) ]
  dt_filt[, LOG2FC_R := log2(
    pmin(
      R1_EXPRESSION_OLD,
      R2_EXPRESSION_OLD,
      R3_EXPRESSION_OLD,
      na.rm = TRUE
    )
    /
      pmin(
        R1_EXPRESSION_YOUNG,
        R2_EXPRESSION_YOUNG,
        R3_EXPRESSION_YOUNG,
        na.rm = TRUE
      )
  ) ]
  dt_filt[, LOG2FC_L := ifelse(is.infinite(LOG2FC_L) & LOG2FC_L > 0, axis_lim,
                               ifelse(is.infinite(LOG2FC_L) & LOG2FC_L < 0, -axis_lim, LOG2FC_L))]
  dt_filt[, LOG2FC_R := ifelse(is.infinite(LOG2FC_R) & LOG2FC_R > 0, axis_lim,
                               ifelse(is.infinite(LOG2FC_R) & LOG2FC_R < 0, -axis_lim, LOG2FC_R))]
  ggplot(dt_filt, aes(x = LOG2FC_L, y = LOG2FC_R, color = REGULATION_SIMPLE)) +
    geom_point() +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    xlim(c(-axis_lim-1, axis_lim+1)) +
    ylim(c(-axis_lim-1,axis_lim+1))
}

coherence_plot(DATASETS_light$facs_mixed, 12)
coherence_plot(DATASETS_light$droplet_mixed, 9)
coherence_plot(DATASETS_light$calico, 4)

####
####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - November 2020
##
## Perfom ORA analysis and result interpretation
## per tissue.
##
####################################################
##

## Libraries ####
library(scDiffCom)
library(data.table)
library(ggplot2)

## Specify data directory ####
dir_data_analysis <- "../data_scAgeCom/analysis/"

## Load scDiffCom (light) results ####
DATASETS_light <- readRDS(paste0(dir_data_analysis, "a4_data_results_3cases.rds"))

names(DATASETS_light)
#DATASETS_light <- DATASETS_light[c(1,4,8)]

## Add specific terms to cell-types and genes ####

#cell-type families
cell_types_dt <- setDT(read.csv(paste0(dir_data_analysis, "scDiffCom_cell_types.csv"), stringsAsFactors = FALSE))

#genage longevity genes
genage_mouse <- setDT(read.csv(paste0(dir_data_analysis, "genage_mouse.tsv"), sep = "\t", header = TRUE))

#add the terms to the data.table
DATASETS_light <- lapply(
  DATASETS_light,
  function(dataset) {
    lapply(
      dataset,
      function(tiss) {
        dt <- scDiffCom:::get_cci_table_filtered(tiss)
        if(identical(dt, list())) return(tiss)
        dt[cell_types_dt, on = "L_CELLTYPE==scDiffCom.cell.type", L_CELL_FAMILY := i.Family...broad]
        dt[cell_types_dt, on = "R_CELLTYPE==scDiffCom.cell.type", R_CELL_FAMILY := i.Family...broad]
        dt[, LR_CELL_FAMILY := paste(L_CELL_FAMILY, R_CELL_FAMILY, sep = "_")]
        dt[,GENAGE := ifelse(
          LIGAND_1 %in% genage_mouse$Gene.Symbol | LIGAND_2 %in% genage_mouse$Gene.Symbol |
            RECEPTOR_1 %in% genage_mouse$Gene.Symbol | RECEPTOR_2 %in% genage_mouse$Gene.Symbol,
          "longevity_associated",
          "not_longevity_associated"
        )
        ]
        tiss <- scDiffCom:::set_cci_table_filtered(tiss, dt)
      }
    )
  }
)

## Redo ORA with new categories ####

DATASETS_light <- lapply(
  DATASETS_light,
  function(dataset) {
    lapply(
      dataset,
      function(tiss) {
        scDiffCom::run_ora(
          object = tiss,
          categories = c("LR_CELLTYPE", "LR_CELL_FAMILY", "LR_NAME", "GO", "GENAGE"),
          overwrite = TRUE,
          logfc_threshold = log(1.2)
        )
      }
    )
  }
)

## Save results with new ORA ####

saveRDS(DATASETS_light, paste0(dir_data_analysis, "a4_data_results_3cases_with_ora.rds"))

DATASETS_light <- readRDS(paste0(dir_data_analysis, "a4_data_results_all_cases.rds"))

facs_aorta <- get_cci_table_filtered(DATASETS_light$facs_mixed$Lung)

facs_aorta[, score_age := LOGFC*(-log(PVAL_DIFF+1E-4))]

unique(facs_aorta$LIGAND_1)
unique(facs_aorta$RECEPTOR_1)

test <- facs_aorta[, c("LR_CELLTYPE", "LR_NAME", "score_age", "LR_SCORE_OLD", "LR_SCORE_YOUNG", "REGULATION_SIMPLE", "REGULATION")]
test[, avg_score := (LR_SCORE_OLD + LR_SCORE_YOUNG)/2]
test[, max_score := pmax(LR_SCORE_OLD, LR_SCORE_YOUNG)]



scores <- c(test$LR_SCORE_OLD, test$LR_SCORE_YOUNG)
scorse_norm <- (scores-min(scores))/(max(scores) - min(scores))
scorse_z <- (scores - mean(scores))/sd(scores)
hist(scorse_z, breaks = 100)

test$LR_SCORE_OLD_norm <- scorse_norm[1:nrow(test)]
test$LR_SCORE_YOUNG_norm <- scorse_norm[(nrow(test)+1):(2*nrow(test))]

test[, CCI := paste(LR_CELLTYPE, LR_NAME, sep = "_")]


test2 <- test[, c("CCI", "LR_SCORE_YOUNG_norm", "LR_SCORE_OLD_norm")]
ggplot(test2, aes(x = LR_SCORE_YOUNG_norm, xend = LR_SCORE_OLD_norm, y = y1, yend = y2)) + geom_segment() + scale_x_log10()
test2 <- melt.data.table(test2)
test2[, age := ifelse(grepl("YOUNG", variable), "YOUNG", "OLD") ]

test2$y1 <- 0
test2$y2 <- 1

ggplot(test2, aes(x = variable, y = value)) + geom_point() + scale_y_log10() + geom_segment(aes(xend = YOUNG, yend = OLD))



ggplot(test, aes(log10(LR_SCORE_YOUNG_norm), log10(LR_SCORE_OLD_norm), color = REGULATION)) + geom_point()

ggplot(test, aes(log10(LR_SCORE_YOUNG_norm), log10(LR_SCORE_OLD_norm))) +
  geom_point(aes(colour=score_age)) +
  scale_colour_gradient2()

ggplot(test, aes(log10(LR_SCORE_YOUNG), log10(LR_SCORE_OLD))) +
  geom_point(aes(colour=score_age)) +
  scale_colour_gradient(low = "blue", high = "red")

hist(scorse_norm, breaks = 50)

hist(c(test$LR_SCORE_OLD, test$LR_SCORE_YOUNG), breaks = 100)

ggplot(test, aes(x = log10(max_score), y = score_age, color = REGULATION_SIMPLE)) + geom_point()
ggplot(test, aes(x = log10(avg_score), y = score_age, color = REGULATION_SIMPLE)) + geom_point()

test2 <- get_cci_table_filtered(DATASETS$droplet_mixed$Tongue)
test3 <- get_cci_table_raw(DATASETS$droplet_mixed$Tongue)

test4 <- test3[LR_SORTED %in% test2$LR_SORTED]
test4[, LOGFC := log(LR_SCORE_OLD/LR_SCORE_YOUNG)]
test4[, LR_CELLTYPE := paste(L_CELLTYPE, R_CELLTYPE, sep = "_")]



min(test4$LOGFC)


test5 <- dcast.data.table(
  test4[, c("LR_SORTED", "LR_CELLTYPE", "LOGFC")],
  formula = LR_SORTED ~ LR_CELLTYPE,
  value.var = "LOGFC"
)

test6 <- as.matrix(test5, rownames = "LR_SORTED")

pheatmap::pheatmap(test6, scale = "row")
pheatmap::pheatmap(test6)
heatmap(test6, scale = "row")

ComplexHeatmap::Heatmap(t(test6), clustering_distance_rows = "pearson")
ComplexHeatmap::Heatmap(t(test6))


ComplexHeatmap::Heatmap(test6, col = rev(rainbow(10)))
ComplexHeatmap::Heatmap(test6, col = colorRamp2(seq(min(test6), max(test6), length = 3), c("blue", "#EEEEEE", "red")))

densityHeatmap(t(test6))

saveRDS(DATASETS_light, "../data_scAgeCom/analysis/a4_data_results_all_cases.rds")

##########


####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - November 2020
##
## Perfom ORA analysis and result interpretation
## per tissue.
##
####################################################
##

## Libraries ####
library(scDiffCom)
library(data.table)
library(ggplot2)

## Specify data directory ####
dir_data_analysis <- "../data_scAgeCom/analysis/"

## Load Data ####
DATASETS_light <- readRDS(paste0(dir_data_analysis, "a4_data_results_3cases_with_ora.rds"))

## Create full table per dataset #####

DATASETS_full_cci_table <- lapply(
  DATASETS_light,
  function(dataset) {
    rbindlist(
      lapply(
        dataset,
        function(tiss) {
          get_cci_table_filtered(tiss)
        }
      ),
      use.names = TRUE,
      idcol = "Tissue",
      fill = TRUE
    )  
  }
)

## Create new scDiffCom objects that combines all tissues and do ORA analysis ####

DATASETS_light_combined <- lapply(
  seq_along(DATASETS_light),
  function(i) {
    temp <- DATASETS_light[[i]]$Spleen
    temp_param <- parameters(temp)
    temp_param$object_name <- "Combined_Tissues"
    parameters(temp) <- temp_param
    temp <- set_cci_table_filtered(temp, DATASETS_full_cci_table[[i]] )
    temp <- run_ora(temp, overwrite = TRUE, categories = c("LR_CELL_FAMILY", "GO", "LR_CELLTYPE", "LR_NAME", "GENAGE"))
  }
)
names(DATASETS_light_combined) <- names(DATASETS_light)

## Save/read the results ####

saveRDS(DATASETS_light_combined, paste0(dir_data_analysis, "a6_data_combinedTissues_3cases.rds"))
DATASETS_light_combined <- readRDS(paste0(dir_data_analysis, "a6_data_combinedTissues_3cases.rds"))

## We need to filter the ORA for only the LR-interactions that appears in at least 3 tissues ####

ora_test <- get_ora_tables(DATASETS_light_combined$facs_mixed)$LR_NAME

counts_test <- DATASETS_light_combined$facs_mixed@cci_table_filtered[
  REGULATION_SIMPLE == "DOWN", .N, by = c("Tissue", "LR_NAME")
  ][, .N, by = c("LR_NAME")]

dc_test <- merge.data.table(
  ora_test,
  counts_test,
  by.x = "Value",
  by.y = "LR_NAME",
  all.x = TRUE
)
setnafill(dc_test, fill = 0, cols = "N")

dc_test[, SCORE_DOWN := -log10(pval_adjusted_DOWN)*log(OR_DOWN)]

ggplot(dc_test[N>4][order(-SCORE_DOWN)][1:30], aes(x = SCORE_DOWN, y = reorder(Value, SCORE_DOWN))) +
  geom_point(aes(size =N))

## Try to do it for GO terms as well #####

## Gene centric analysis for FACS #####

#retrieve all CCIs
facs_cci <- copy(DATASETS_full_cci_table$facs_mixed)
facs_cci[, N_CC := (uniqueN(L_CELLTYPE))^2, by = "Tissue"]
facs_cci[, N_TISSUES_FOR_LR := uniqueN(Tissue), by = "LR_NAME"]
facs_cci[, N_CCI_PER_TISSUE_FOR_LR := .N, by = c("LR_NAME", "Tissue") ]
facs_cci[, PCT_CCI_PER_TISSUE_FOR_LR := N_CCI_PER_TISSUE_FOR_LR/N_CC]
facs_cci[, AVG_PCT_CCI_PER_TISSUE_FOR_LR := mean(PCT_CCI_PER_TISSUE_FOR_LR), by = "LR_NAME"]
facs_cci[facs_cci[REGULATION_SIMPLE == "UP", uniqueN(Tissue), by = "LR_NAME"], on = "LR_NAME", N_TISSUES_FOR_LR_UP := i.V1 ]
setnafill(facs_cci, fill = 0, cols = "N_TISSUES_FOR_LR_UP")
facs_cci[facs_cci[REGULATION_SIMPLE == "DOWN", uniqueN(Tissue), by = "LR_NAME"], on = "LR_NAME", N_TISSUES_FOR_LR_DOWN := i.V1 ]
setnafill(facs_cci, fill = 0, cols = "N_TISSUES_FOR_LR_DOWN")

#retrive cell-families

facs_family_dt <- facs_cci[,{
  totwt = .N
  .SD[,.(frac=.N/totwt),by=LR_CELL_FAMILY]
},by=LR_NAME]

facs_family_dt[, main_family := max(frac), by = "LR_NAME"]
facs_family_dt[, MAIN_LR_CELL_FAMILY := paste(.SD[frac == max(frac)]$LR_CELL_FAMILY, collapse = "/"), by = "LR_NAME"]


#retrieve all ora for each tissue
facs_ora <- rbindlist(
  lapply(
    DATASETS_light$facs_mixed,
    function(i) {
      dt <- get_ora_tables(i)$LR_NAME
    }
  ),
  use.names = TRUE,
  idcol = "Tissue"
)
facs_ora <- facs_ora[, c(1,2,8,12,13, 18,22, 23)]
facs_ora[, regulation := ifelse(
  pval_adjusted_UP <= 0.05 & OR_UP > 1,
  ifelse(
    pval_adjusted_DOWN <= 0.05 & OR_DOWN > 1,
    "BOTH",
    "UP"
  ),
  ifelse(
    pval_adjusted_DOWN <= 0.05 & OR_DOWN > 1,
    "DOWN",
    "NONE"
  )
)
]

facs_ora_regulation <- dcast.data.table(
  facs_ora[,.N, by = c("Value", "regulation")],
  formula = Value ~ regulation,
  value.var = "N"
)
facs_ora_regulation <- facs_ora_regulation[, c("Value", "UP", "DOWN")]
facs_ora_regulation[is.na(facs_ora_regulation)] <- 0
setnames(facs_ora_regulation, old = c("UP", "DOWN"), new = c("N_TISSUES_ORA_UP", "N_TISSUE_ORA_DOWN"))

# retrieve ora on combined tissues
facs_ora_combined <- DATASETS_light_combined$facs_mixed@ora_tables$LR_NAME
facs_ora_combined <- facs_ora_combined[, c(1,2,7,11,12,17,21,22)]


# get LR_genes summary table
facs_LR_GENES_summary <- unique(facs_cci[, c("LR_NAME", "N_TISSUES_FOR_LR", "AVG_PCT_CCI_PER_TISSUE_FOR_LR",
                                             "N_TISSUES_FOR_LR_UP","N_TISSUES_FOR_LR_DOWN")])
facs_LR_GENES_summary <- merge.data.table(
  facs_LR_GENES_summary,
  facs_ora_regulation,
  by.x = "LR_NAME",
  by.y = "Value",
  all.x = TRUE
)
facs_LR_GENES_summary <- merge.data.table(
  facs_LR_GENES_summary,
  facs_ora_combined,
  by.x = "LR_NAME",
  by.y = "Value",
  all.x = TRUE
)

facs_LR_GENES_summary <- merge.data.table(
  facs_LR_GENES_summary,
  unique(facs_family_dt[,c("LR_NAME", "MAIN_LR_CELL_FAMILY")]),
  by= "LR_NAME"
)

ggplot(facs_LR_GENES_summary, aes(x = N_TISSUES_ORA_UP, y = ORA_score_UP)) + geom_point() + stat_smooth(model = "lm", formula=y~x)

## Gene centric analysis for Droplet #####

#retrieve all CCIs
droplet_cci <- copy(DATASETS_full_cci_table$droplet_mixed)
droplet_cci[, N_CC := (uniqueN(L_CELLTYPE))^2, by = "Tissue"]
droplet_cci[, N_TISSUES_FOR_LR := uniqueN(Tissue), by = "LR_NAME"]
droplet_cci[, N_CCI_PER_TISSUE_FOR_LR := .N, by = c("LR_NAME", "Tissue") ]
droplet_cci[, PCT_CCI_PER_TISSUE_FOR_LR := N_CCI_PER_TISSUE_FOR_LR/N_CC]
droplet_cci[, AVG_PCT_CCI_PER_TISSUE_FOR_LR := mean(PCT_CCI_PER_TISSUE_FOR_LR), by = "LR_NAME"]
droplet_cci[droplet_cci[REGULATION_SIMPLE == "UP", uniqueN(Tissue), by = "LR_NAME"], on = "LR_NAME", N_TISSUES_FOR_LR_UP := i.V1 ]
setnafill(droplet_cci, fill = 0, cols = "N_TISSUES_FOR_LR_UP")
droplet_cci[droplet_cci[REGULATION_SIMPLE == "DOWN", uniqueN(Tissue), by = "LR_NAME"], on = "LR_NAME", N_TISSUES_FOR_LR_DOWN := i.V1 ]
setnafill(droplet_cci, fill = 0, cols = "N_TISSUES_FOR_LR_DOWN")

#retrive cell-families

droplet_family_dt <- droplet_cci[,{
  totwt = .N
  .SD[,.(frac=.N/totwt),by=LR_CELL_FAMILY]
},by=LR_NAME]

droplet_family_dt[, main_family := max(frac), by = "LR_NAME"]
droplet_family_dt[, MAIN_LR_CELL_FAMILY := paste(.SD[frac == max(frac)]$LR_CELL_FAMILY, collapse = "/"), by = "LR_NAME"]


#retrieve all ora for each tissue
droplet_ora <- rbindlist(
  lapply(
    DATASETS_light$droplet_mixed,
    function(i) {
      dt <- get_ora_tables(i)$LR_NAME
    }
  ),
  use.names = TRUE,
  idcol = "Tissue"
)
droplet_ora <- droplet_ora[, c(1,2,8,12,13, 18,22, 23)]
droplet_ora[, regulation := ifelse(
  pval_adjusted_UP <= 0.05 & OR_UP > 1,
  ifelse(
    pval_adjusted_DOWN <= 0.05 & OR_DOWN > 1,
    "BOTH",
    "UP"
  ),
  ifelse(
    pval_adjusted_DOWN <= 0.05 & OR_DOWN > 1,
    "DOWN",
    "NONE"
  )
)
]

droplet_ora_regulation <- dcast.data.table(
  droplet_ora[,.N, by = c("Value", "regulation")],
  formula = Value ~ regulation,
  value.var = "N"
)
droplet_ora_regulation <- droplet_ora_regulation[, c("Value", "UP", "DOWN")]
droplet_ora_regulation[is.na(droplet_ora_regulation)] <- 0
setnames(droplet_ora_regulation, old = c("UP", "DOWN"), new = c("N_TISSUES_ORA_UP", "N_TISSUE_ORA_DOWN"))

# retrieve ora on combined tissues
droplet_ora_combined <- DATASETS_light_combined$droplet_mixed@ora_tables$LR_NAME
droplet_ora_combined <- droplet_ora_combined[, c(1,2,7,11,12,17,21,22)]


# get LR_genes summary table
droplet_LR_GENES_summary <- unique(droplet_cci[, c("LR_NAME", "N_TISSUES_FOR_LR", "AVG_PCT_CCI_PER_TISSUE_FOR_LR",
                                                   "N_TISSUES_FOR_LR_UP","N_TISSUES_FOR_LR_DOWN")])
droplet_LR_GENES_summary <- merge.data.table(
  droplet_LR_GENES_summary,
  droplet_ora_regulation,
  by.x = "LR_NAME",
  by.y = "Value",
  all.x = TRUE
)
droplet_LR_GENES_summary <- merge.data.table(
  droplet_LR_GENES_summary,
  droplet_ora_combined,
  by.x = "LR_NAME",
  by.y = "Value",
  all.x = TRUE
)

droplet_LR_GENES_summary <- merge.data.table(
  droplet_LR_GENES_summary,
  unique(droplet_family_dt[,c("LR_NAME", "MAIN_LR_CELL_FAMILY")]),
  by= "LR_NAME"
)

ggplot(droplet_LR_GENES_summary, aes(x = N_TISSUES_ORA_UP, y = ORA_score_UP)) + geom_point() + stat_smooth(model = "lm", formula=y~x)

## Save genes summary ####

LR_GENES_summary <- list(
  facs_mixed = facs_LR_GENES_summary,
  droplet_mixed = droplet_LR_GENES_summary
)

saveRDS(LR_GENES_summary,"../scAgeCom/shinyApp/data/genes_summary_shiny.rds")

###########
#table(facs_cci$LR_CELL_FAMILY)
#test <- facs_cci[, .N, by = c("LR_NAME", "Tissue")]
#ggplot(test, aes(x=Tissue, y = LR_NAME)) + geom_point(aes(size = N))
#test2 <- test[, list(text = paste(substring(Tissue, 1, 3), collapse = "-")), by = "LR_NAME"]

#test3 <- facs_cci[, .N, by = c("LR_NAME", "LR_CELL_FAMILY")]



#a <- testx[LR_NAME %in% facs_LR_GENES_summary[N_TISSUE_ORA_DOWN > 4]$LR_NAME]
#b <- sapply(
#  unique(testx$LR_NAME),
#  function(i) {
#    testx[LR_NAME == i][order(-frac)]$LR_CELL_FAMILY[[1]]
#  },
#  USE.NAMES = TRUE
#)

testx[, is_max := ifelse()]

facs_LR_GENES_summary[N_TISSUES_ORA_UP > 6]$LR_NAME

ggplot(testx[LR_NAME %in% facs_LR_GENES_summary[N_TISSUE_ORA_DOWN > 4]$LR_NAME], aes(x= LR_CELL_FAMILY, y = LR_NAME)) + geom_point(aes(size = frac))


testxx <- as.matrix(dcast.data.table(
  testx,
  formula = LR_NAME ~ LR_CELL_FAMILY,
  value.var = "frac"
),
rownames = "LR_NAME"
)
testxx[is.na(testxx)] <- 0

heatmap(testxx)

library(circlize)
ComplexHeatmap::Heatmap(
  testxx,
  #clustering_distance_rows = "pearson",
  show_column_names = FALSE,
  show_row_names = FALSE
)

ComplexHeatmap::Heatmap(
  testxx,
  row_km = 10,
  column_km = 10,
  show_column_names = FALSE,
  show_row_names = FALSE
)

cl = kmeans(testxx, centers = 10)
cr = kmeans(t(testxx), centers = 10)

cl_list <- lapply(
  unique(cl$cluster),
  function(i) {
    names(cl$cluster[cl$cluster == i])
  }
)

cr_list <- lapply(
  unique(cr$cluster),
  function(i) {
    names(cr$cluster[cr$cluster == i])
  }
)

cl$cluster

library(seriation)
o = seriate(testxx, method = "BEA_TSP")
Heatmap(testxx, name = "mat", 
        row_order = get_order(o, 1), column_order = get_order(o, 2),
        column_title = "seriation by BEA_TSP method")


pheatmap::pheatmap(testxx, scale = "row")
pheatmap::pheatmap(testxx, scale = "column")

test4 <- test3[, list(text = paste(LR_CELL_FAMILY, collapse = "_")), by = "LR_NAME"]
table(test4$text)



ggplot(test3, aes(x= LR_CELL_FAMILY, y = LR_NAME)) + geom_point(aes(size = N))

facs_LR_GENES_summary
###########



## Tissue specificity#####




##################



#for each LR_NAME we count several values



### divergent signals?

ora_facs <- DATASETS_light_combined$facs_mixed@ora_tables$LR_NAME
ora_facs$ORA_score_UP

ora_facs[pval_adjusted_UP <= 0.05 & pval_adjusted_DOWN <= 0.05 & OR_UP > 1 & OR_DOWN > 1]





test <- dcast.data.table(
  ora_facs_tissues_bind[, c("Tissue", "Value", "regulation")],
  formula = Value ~ Tissue, 
  value.var = "regulation"
)

## Plot some ORA results ####

plot_ora(
  DATASETS_light_combined$facs,
  category = "GO",
  OR_val = "OR_DOWN",
  pval_val = "pval_DOWN",
  ORA_score_val = "ORA_score_DOWN", 
  max_value = 10
)

plot_ora(
  DATASETS_light_combined$droplet_mixed,
  category = "LR_NAME",
  OR_val = "OR_UP",
  pval_val = "pval_UP",
  ORA_score_val = "ORA_score_UP", 
  max_value = 50
)

plot_ora(
  DATASETS_light_combined$facs_mixed,
  category = "LR_NAME",
  OR_val = "OR_UP",
  pval_val = "pval_UP",
  ORA_score_val = "ORA_score_UP", 
  max_value = 10
)




## Counts per tissue ####

DATASETS_light_combined$droplet_mixed@cci_table_filtered[
  REGULATION_SIMPLE == "UP", .N, by = c("Tissue", "LR_NAME")
  ][
    , .N, by = c("LR_NAME")
    ][order(-N)][LR_NAME == "Scgb1a1:Lmbr1l"]

DATASETS_full_cci_table$droplet_mixed[REGULATION_SIMPLE == "UP", .N, by = c("Tissue", "LR_CELLTYPE")][, .N, by = c("LR_CELLTYPE")][order(-N)]

DATASETS_full_cci_table$facs_mixed[REGULATION_SIMPLE == "UP", .N, by = c("Tissue", "LR_NAME")][, .N, by = c("LR_NAME")][order(-N)]
DATASETS_full_cci_table$facs_mixed[REGULATION_SIMPLE == "UP", .N, by = c("Tissue", "LR_CELLTYPE")][, .N, by = c("LR_CELLTYPE")][order(-N)]

## LR_name specific
facs_table <- DATASETS_full_cci_table$facs_mixed

LR_GENES_facs <- sort(unique(facs_table$LR_NAME))

facs_table[LR_NAME == LR_GENES_facs[[1]]][,.N, by = c("Tissue", "REGULATION_SIMPLE")]
hist(facs_table[LR_NAME == LR_GENES_facs[[1]]]$LOGFC, breaks = 50)

hist(facs_table[GENAGE == "longevity_associated"]$LOGFC, breaks = 100)



##

#top 10 LR_pairs in common
top10_facs_LR_GENES <- DATASETS_full_cci_table$facs_mixed[REGULATION_SIMPLE == "UP", .N, by = c("Tissue", "LR_NAME")][, .N, by = c("LR_NAME")][order(-N)][1:10]

full_cci_table_droplet <- DATASETS_full_cci_table$droplet_mixed
full_cci_table_droplet[, Tissue_CELLTYPES := paste0(Tissue, LR_CELLTYPE)]
full_cci_table_droplet[, pct := .N/1315*100, by = LR_NAME]
full_cci_table_droplet[, N_LR_NAME := .N, by = LR_NAME]

full_cci_table_droplet[order(-pct)]

full_cci_table_droplet[REGULATION_SIMPLE == "UP", .N, by = LR_NAME][order(-N)][1:20]
full_cci_table_droplet[REGULATION_SIMPLE == "DOWN", .N, by = LR_NAME][order(-N)][1:20]


full_cci_table_facs <- DATASETS_full_cci_table$facs_mixed
full_cci_table_facs[, Tissue_CELLTYPES := paste0(Tissue, LR_CELLTYPE)]
full_cci_table_facs[, pct := .N/1315*100, by = LR_NAME]
full_cci_table_facs[, N_LR_NAME := .N, by = LR_NAME]

full_cci_table_facs[order(-pct)]

full_cci_table_facs[REGULATION_SIMPLE == "UP", .N, by = LR_NAME][order(-N)][1:20]
full_cci_table_facs[REGULATION_SIMPLE == "DOWN", .N, by = LR_NAME][order(-N)][1:30]

##


full_cci_table_droplet[LR_NAME == "Lgals1:Itgb1"][, .N, by = L_CELLTYPE][order(-N)][1:20]
full_cci_table_droplet[LR_NAME == "Lgals1:Itgb1"][, .N, by = R_CELLTYPE][order(-N)][1:20]
full_cci_table_droplet[LR_NAME == "Lgals1:Itgb1"][, .N, by = LR_CELLTYPE][order(-N)][1:20]
full_cci_table_droplet[LR_NAME == "Lgals1:Itgb1"][, .N, by = REGULATION_SIMPLE]
full_cci_table_droplet[LR_NAME == "Lgals1:Itgb1"][, .N, by = Tissue][order(-N)][1:20]
test <- full_cci_table_droplet[LR_NAME == "Lgals1:Itgb1"][, .N, by = c("Tissue", "REGULATION_SIMPLE")]


##



#metalloprotease adam, Timp2:Itgb1, etc


full_cci_table_facs[LR_NAME == "Hspa8:Adrb2"][, .N, by = L_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Hspa8:Adrb2"][, .N, by = R_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Hspa8:Adrb2"][, .N, by = LR_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Hspa8:Adrb2"][, .N, by = REGULATION_SIMPLE]
full_cci_table_facs[LR_NAME == "Hspa8:Adrb2"][, .N, by = Tissue][order(-N)][1:20]
test <- full_cci_table_facs[LR_NAME == "Hspa8:Adrb2"][, .N, by = c("Tissue", "REGULATION_SIMPLE")]

#heat-shock protein related to aging, what's its role in IC, with adrenergic receptors...
#see https://pubmed.ncbi.nlm.nih.gov/23153586/
#postulate new mechanism? e.g. development, see heart, marrow
#The present data provide new evidence that serum concentration of Hsp70 decreases with age in a normal population.
#Our study also shows that higher levels of Hsp70 are associated with inflammation and frailty in elderly patients.

full_cci_table_facs[LR_NAME == "Psen1:Notch2"][, .N, by = L_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Psen1:Notch2"][, .N, by = R_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Psen1:Notch2"][, .N, by = LR_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Psen1:Notch2"][, .N, by = REGULATION_SIMPLE]
full_cci_table_facs[LR_NAME == "Psen1:Notch2"][, .N, by = Tissue][order(-N)][1:20]
test <- full_cci_table_facs[LR_NAME == "Psen1:Notch2"][, .N, by = c("Tissue", "REGULATION_SIMPLE")]

#role in AD


full_cci_table_facs[LR_NAME == "Col3a1:Itgb1"][, .N, by = L_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Col3a1:Itgb1"][, .N, by = R_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Col3a1:Itgb1"][, .N, by = LR_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Col3a1:Itgb1"][, .N, by = REGULATION_SIMPLE]
full_cci_table_facs[LR_NAME == "Col3a1:Itgb1"][, .N, by = Tissue][order(-N)][1:20]
test <- full_cci_table_facs[LR_NAME == "Col3a1:Itgb1"][, .N, by = c("Tissue", "REGULATION_SIMPLE")]

#colagen, integrin stuff ###
# The bulk of the current literature suggests that collagen-binding integrins only have a limited role
#in adult connective tissue homeostasis, partly due to a limited availability of cell-binding sites in the 
#mature fibrillar collagen matrices. However, some recent data suggest that, instead, they are more crucial for
#dynamic connective tissue remodeling events 
#– such as wound healing – where they might act specifically to remodel and restore the tissue architecture. 

full_cci_table_facs[LR_NAME == "Uba52:Tgfbr2"][, .N, by = L_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Uba52:Tgfbr2"][, .N, by = R_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Uba52:Tgfbr2"][, .N, by = LR_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Uba52:Tgfbr2"][, .N, by = REGULATION_SIMPLE]
full_cci_table_facs[LR_NAME == "Uba52:Tgfbr2"][, .N, by = Tissue][order(-N)][1:20]
test <- full_cci_table_facs[LR_NAME == "Uba52:Tgfbr2"][, .N, by = c("Tissue", "REGULATION_SIMPLE")]


#Ubb up but Uba52 and Ubc down!!!!

full_cci_table_facs[LR_NAME == "Lgals1:Ptprc"][, .N, by = L_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Lgals1:Ptprc"][, .N, by = R_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Lgals1:Ptprc"][, .N, by = LR_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Lgals1:Ptprc"][, .N, by = REGULATION_SIMPLE]
full_cci_table_facs[LR_NAME == "Lgals1:Ptprc"][, .N, by = Tissue][order(-N)][1:20]
test <- full_cci_table_facs[LR_NAME == "Lgals1:Ptprc"][, .N, by = c("Tissue", "REGULATION_SIMPLE")]

#link with apoptosis ## T cell - T cell death is going down!!!!!
#Although not absolutely required for susceptibility to galectin-1,
#CD45 is a major receptor for galectin-1 on T cells, acts as a negative and positive regulator of galectin-1 death,
#and enhances phagocytic clearance of cells killed by galectin-1 (11,–15).


full_cci_table_facs[LR_NAME == "App:Ncstn"][, .N, by = L_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "App:Ncstn"][, .N, by = R_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "App:Ncstn"][, .N, by = LR_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "App:Ncstn"][, .N, by = REGULATION_SIMPLE]
full_cci_table_facs[LR_NAME == "App:Ncstn"][, .N, by = Tissue][order(-N)][1:20]
test <- full_cci_table_facs[LR_NAME == "App:Ncstn"][, .N, by = c("Tissue", "REGULATION_SIMPLE")]

#related to AD, but why in so many other tissues! see related article! new role of App!
#APP function and Aβ pathology in adipose tissue
#APP function and Aβ pathology in muscle

full_cci_table_facs[LR_NAME == "Apoe:Sorl1"][, .N, by = L_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Apoe:Sorl1"][, .N, by = R_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Apoe:Sorl1"][, .N, by = LR_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Apoe:Sorl1"][, .N, by = REGULATION_SIMPLE]
full_cci_table_facs[LR_NAME == "Apoe:Sorl1"][, .N, by = Tissue][order(-N)]

#Link with AD!!

full_cci_table_facs[LR_NAME == "Ltb:Tnfrsf1a"][, .N, by = L_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Ltb:Tnfrsf1a"][, .N, by = R_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Ltb:Tnfrsf1a"][, .N, by = LR_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Ltb:Tnfrsf1a"][, .N, by = REGULATION_SIMPLE]
full_cci_table_facs[LR_NAME == "Ltb:Tnfrsf1a"][, .N, by = Tissue][order(-N)]

full_cci_table_facs[LR_NAME == "Calm1:Ptpra"][, .N, by = L_CELLTYPE][order(-N)]
full_cci_table_facs[LR_NAME == "Calm1:Ptpra"][, .N, by = R_CELLTYPE][order(-N)]
full_cci_table_facs[LR_NAME == "Calm1:Ptpra"][, .N, by = LR_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Calm1:Ptpra"][, .N, by = REGULATION_SIMPLE]

full_cci_table_facs[LR_NAME == "Ubb:Ripk1"][, .N, by = REGULATION_SIMPLE]
full_cci_table_facs[LR_NAME == "Ubb:Ripk1"][, .N, by = L_CELLTYPE][order(-N)]
full_cci_table_facs[LR_NAME == "Ubb:Ripk1"][, .N, by = R_CELLTYPE][order(-N)]
full_cci_table_facs[LR_NAME == "Ubb:Ripk1"][, .N, by = LR_CELLTYPE][order(-N)][1:20]

full_cci_table_facs[LR_NAME == "B2m:Hfe"][, .N, by = REGULATION_SIMPLE]
full_cci_table_facs[LR_NAME == "B2m:Hfe"][, .N, by = L_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "B2m:Hfe"][, .N, by = R_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "B2m:Hfe"][, .N, by = LR_CELLTYPE][order(-N)][1:20]
## look at iron stuff

full_cci_table_facs[LR_NAME == "Calm2:Insr"][, .N, by = REGULATION_SIMPLE]
full_cci_table_facs[LR_NAME == "Calm2:Insr"][, .N, by = L_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Calm2:Insr"][, .N, by = R_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Calm2:Insr"][, .N, by = LR_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Calm2:Insr"][, .N, by = Tissue][order(-N)][1:20]

full_cci_table_facs[LR_NAME == "Hmgb1:Thbd"][, .N, by = REGULATION_SIMPLE]
full_cci_table_facs[LR_NAME == "Hmgb1:Thbd"][, .N, by = L_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Hmgb1:Thbd"][, .N, by = R_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Hmgb1:Thbd"][, .N, by = LR_CELLTYPE][order(-N)][1:20]
full_cci_table_facs[LR_NAME == "Hmgb1:Thbd"][, .N, by = Tissue][order(-N)][1:20]

LR6db$LR6db_GO$LR_GO_intersection[LR_SORTED == "Hmgb1_Thbd"]

#anti-inflammatory signal !!!

full_cci_table_facs[LR_NAME == "Ubb:Ripk1", c("Tissue", "LR_CELLTYPE", "LOGFC", "REGULATION")][order(-LOGFC)]
hist(full_cci_table_facs[LR_NAME == "Ubb:Ripk1"]$LOGFC, breaks = 50)
summary(full_cci_table_facs[LR_NAME == "Ubb:Ripk1"]$LOGFC)

test <- full_cci_table_facs[Tissue == "Aorta"]

test <- dcast(full_cci_table_facs[REGULATION == "DOWN", .N, by = c("LR_NAME", "Tissue")], LR_NAME ~ Tissue, value.var = "N")
test[is.na(test)] <- 0

cols_tissues <- colnames(test)[2:24]

test[, N_sum := rowSums(.SD), .SDcols = cols_tissues]
test[, paste0(cols_tissues, "_specific") := lapply(
  cols_tissues,
  function(tiss) {
    ifelse(get(tiss) == N_sum, TRUE, FALSE)
  }
)]

test2 <- lapply(
  cols_tissues,
  function(tiss) {
    test[get(paste0(tiss, "_specific")) == TRUE][["LR_NAME"]]
  }
)
names(test2) <- cols_tissues

test2$Diaphragm

lapply(
  seq_along(DATASETS_light$facs_mixed),
  function(i) {
    intersect(
      DATASETS_light$facs_mixed[[i]]@ora_tables$LR_NAME[OR_DOWN > 1 & pval_adjusted_DOWN <= 0.05][["Value"]],
      test2[[names(DATASETS_light$facs_mixed)[[i]]]]
    )
  }
)

DATASETS_light$facs_mixed$Aorta@ora_tables$LR_NAME[OR_DOWN >1 & pval_adjusted_DOWN <= 0.05]

test[, aort_test := ifelse(Aorta == N_sum, TRUE, FALSE)]

test <- full_cci_table_facs[LR_NAME %in% top10_facs_LR_GENES$LR_NAME, c("LR_NAME", "Tissue_CELLTYPES", "LOGFC")]
hist(test$LOGFC, breaks = 100)
top10_facs_mat <- as.matrix(dcast(test, LR_NAME ~ Tissue_CELLTYPES, value.var = "LOGFC"), rownames = "LR_NAME")
top10_facs_mat[is.na(top10_facs_mat)] <- 0

ComplexHeatmap::Heatmap(top10_facs_mat, show_column_names = FALSE)
pheatmap::pheatmap(top10_facs_mat, show_colnames = FALSE)

full_cci_table_facs[, .N, by = LR_NAME][order(-N)]



full_cci_table_facs[, .N/1315*100, by = LR_NAME][order(-V1)][V1 > 10]

hist(full_cci_table_facs[, .N/1315*100, by = LR_NAME][order(-V1)]$V1, breaks = 100)

unique(full_cci_table_facs$Tissue_CELLTYPES)
full_cci_table_facs[,.N, by = Tissue_CELLTYPES][, .N]

table(full_cci_table_facs[, .N, by = LR_NAME][order(-N)]$N)

mat_full_facs <- as.matrix(
  dcast(
    full_cci_table_facs[pct > 10, c("LR_NAME", "Tissue_CELLTYPES", "LOGFC")],
    LR_NAME ~ Tissue_CELLTYPES,
    value.var = "LOGFC"
  ),
  rownames = "LR_NAME"
)
mat_full_facs[is.na(mat_full_facs)] <- 0

hist(full_cci_table_facs$LOGFC, breaks = 100)

min(mat_full_facs)

ComplexHeatmap::Heatmap(mat_full_facs, show_column_names = FALSE, show_row_names = FALSE)
pheatmap::pheatmap(mat_full_facs, show_rownames = FALSE, show_colnames = FALSE )




####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - November 2020
##
####################################################
##

## Load libraries ####

library(data.table)
library(scDiffCom)

## Specify data directory ####
dir_data_analysis <- "../data_scAgeCom/analysis/"

## Load scDiffCom (light) results ####
DATASETS_light <- readRDS(paste0(dir_data_analysis, "a4_data_results_3cases_with_ora.rds"))

names(DATASETS_light)
names(DATASETS_light) <- c("Calico Data", "TMS Droplet Data", "TMS FACS Data")

## Process data for tissue-specific analysis ####

# Add columns and change colnames

DATASETS_light <- lapply(
  DATASETS_light,
  function(dataset) {
    lapply(
      dataset,
      function(tiss) {
        cci_table_filtered <- get_cci_table_filtered(tiss)
        cci_table_filtered[, c("LOG2FC") := list(LOGFC*log2(exp(1)))]
        temp_score <- c(cci_table_filtered$LR_SCORE_YOUNG, cci_table_filtered$LR_SCORE_OLD)
        #temp_score <- (temp_score-median(temp_score))/(quantile(temp_score,0.75)-quantile(temp_score, 0.25))
        temp_score <- (temp_score-min(temp_score))/(quantile(temp_score,0.9)-min(temp_score))
        cci_table_filtered$LR_SCORE_SCALED_YOUNG <- temp_score[1:nrow(cci_table_filtered)]
        cci_table_filtered$LR_SCORE_SCALED_OLD <- temp_score[(1+nrow(cci_table_filtered)):(2*nrow(cci_table_filtered))]
        setnames(
          cci_table_filtered,
          old = c("L_CELLTYPE", "R_CELLTYPE", "BH_PVAL_DIFF", "LR_NAME", "LR_SCORE_SCALED_YOUNG", "LR_SCORE_SCALED_OLD"),
          new = c("Emitter Cell Type", "Receiver Cell Type", "Adj. P-Value", "Ligand-Receptor Genes", "SCORE (YOUNG)", "SCORE (OLD)")
        )
        new_object <- set_cci_table_filtered(tiss, cci_table_filtered)
        ora_tables <- get_ora_tables(new_object)
        ora_tables <- ora_tables[c("GO", "LR_NAME", "LR_CELLTYPE", "LR_CELL_FAMILY")]
        names(ora_tables) <- c("GO Terms", "Ligand-Receptor Genes", "Cell Types", "Cell Families")
        ora_tables <- lapply(
          ora_tables,
          function(ORA_dt) {
            setnames(
              ORA_dt,
              old = c("OR_UP", "pval_adjusted_UP", "OR_DOWN", "pval_adjusted_DOWN",
                      "OR_FLAT", "pval_adjusted_FLAT"),
              new = c("Odds Ratio Up", "Adj. P-Value Up", "Odds Ratio Down", "Adj. P-Value Down",
                      "Odds Ratio Stable", "Adj. P-Value Stable")
            )
            if("Value_NAME" %in% colnames(ORA_dt)) {
              setnames(
                ORA_dt,
                old = c("Value", "Value_NAME"),
                new = c("Value_ID", "Value")
              )
            }
            return(ORA_dt)
          }
        )
        new_object <- scDiffCom:::set_ora_tables(new_object, ora_tables)
        return(new_object)
      }
    )
  }
)

saveRDS(DATASETS_light, "shinyApp/data/scdiffcom_objects_shiny.rds")

