####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - November 2020
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
library(pbapply)

## Specify the directory with scDiffCom results ####
dir_results <- "../data_scAgeCom/scdiffcom_results/"
dir_data_analysis <- "../data_scAgeCom/analysis/"

## Retrieve all the results for each tissue and dataset ####
RESULT_PATHS <- list.dirs(dir_results, recursive = FALSE)

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
#names(DATASETS) <- c("calico_nlog")#, "calico_log", "droplet_nlog", "droplet_log", "facs_nlog", "facs_log")
names(DATASETS) <- c("calico", "droplet_female", "droplet_male", "droplet_mixed", "droplet_sex",
                     "facs_female", "facs_male", "facs_mixed", "facs_sex")

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
              temp_obj <- run_filtering_and_ORA(
                object = object,
                new_cutoff_quantile_score = score_cut,
                new_cutoff_logfc = logfc_cut,
                skip_ORA = TRUE,
                verbose = FALSE
              )
              dt <- scDiffCom:::get_cci_table_filtered(temp_obj)
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

#define a series of cutoffs

logfc_cuts <- c(1.1, 1.2, 1.3, 1.4, 1.5)
names(logfc_cuts) <- logfc_cuts
logfc_cuts <- log(logfc_cuts)

score_cuts <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
names(score_cuts) <- score_cuts

#computation for each dataset and tissue (can take more than one hour)
regulation_behaviour <- rbindlist(
  lapply(
    seq_along(DATASETS),
    function(i) {
      print(names(DATASETS)[i])
      rbindlist(
        lapply(
          seq_along(DATASETS[[i]]),
          function(j) {
            print(names(DATASETS[[i]])[j])
            get_regulation_behaviour(
              object = DATASETS[[i]][[j]],
              logfc_cuts = logfc_cuts,
              score_cuts = score_cuts
            )
          }
        ),
        use.names = TRUE,
        idcol = "Tissue"
      )
    }
  ),
  use.names = TRUE,
  idcol = "Dataset"
)

#plot the behaviour for different tissues and dataset
test <- get_regulation_behaviour(
  DATASETS$facs_mixed$`Brain_Non-Myeloid`,
  logfc_cuts = logfc_cuts,
  score_cuts = score_cuts
)
ggplot(test[REGULATION == "FLAT"], aes(x = logfc_cut, y = pct, color = score_cut)) + geom_point()
ggplot(test[REGULATION == "FLAT"], aes(x = logfc_cut, y = N, color = score_cut)) + geom_point()
ggplot(test[REGULATION == "UP"], aes(x = logfc_cut, y = pct, color = score_cut)) + geom_point()
ggplot(test[REGULATION == "UP_APPEARS"], aes(x = logfc_cut, y = pct, color = score_cut)) + geom_point()
ggplot(test[REGULATION == "DOWN"], aes(x = logfc_cut, y = pct, color = score_cut)) + geom_point()
ggplot(test[REGULATION == "DOWN_DISAPPEARS"], aes(x = logfc_cut, y = pct, color = score_cut)) + geom_point()


## We keep the initial parameters log(1.2) and 0.25 and create a new object without the raw data ####

DATASETS_light <- lapply(
  DATASETS,
  function(i) {
    lapply(
      i,
      function(tiss) {
        tiss <- scDiffCom:::set_cci_table_raw(tiss, list()) 
        return(tiss)
      }
    )
  }
)

saveRDS(DATASETS_light, "../data_scAgeCom/analysis/a4_data_results_all_cases.rds")

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
coherence_plot(DATASETS$droplet_mixed, 9)
coherence_plot(DATASETS$calico, 4)




######



reg_counts <- rbindlist(
  lapply(
    DATASETS_light[c(1,4,8)],
    function(i) {
      rbindlist(
        lapply(
          i,
          function(tiss) {
            dt <- scDiffCom:::get_cci_table_filtered(tiss)
            dt[, .N, by = c("LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG",
                            "LR_DETECTED_AND_SIGNIFICANT_IN_OLD",
                            "REGULATION_SIMPLE") ][, V1 := N/sum(N)*100]
          }
        ),
        use.names = TRUE,
        idcol = "Tissue"
      )
    }
  ),
  use.names = TRUE,
  idcol = "Dataset"
)

reg_counts[, dataset := gsub("\\_.*","", Dataset)]
reg_counts[, temp := paste(dataset, REGULATION_SIMPLE, sep = "_")]
reg_counts <- reg_counts[ order(REGULATION_SIMPLE)]

ggplot(reg_counts, aes(x= REGULATION_SIMPLE, y = V1)) + geom_boxplot() + facet_grid(vars(dataset))
ggplot(reg_counts[Tissue %in% c("Spleen", "Kidney", "Lung")], aes(x= REGULATION_SIMPLE, y = V1)) + geom_boxplot() + facet_grid(vars(dataset))
ggplot(reg_counts[Tissue %in% c("Spleen", "Kidney", "Lung")], aes(x= temp, y = V1)) + geom_boxplot()
