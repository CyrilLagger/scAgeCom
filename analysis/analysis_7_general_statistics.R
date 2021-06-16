####################################################
##
## Project: scAgeCom
##
## Last update - June 2021
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## general statistics
##
####################################################
##

## libraries ####
library(data.table)
library(ggplot2)

## load results ####

cci_stat_table <- readRDS("../data_scAgeCom/analysis/outputs_data/data_4_CCI_table_unprocessed.rds")


## create regulation table for supplemental data and stat ####

regulation_table <- cci_stat_table[
  ,
  .N,
  by = c("IS_CCI_EXPRESSED_YOUNG", "IS_CCI_SPECIFIC_YOUNG", "IS_CCI_SCORE_YOUNG",
         "IS_CCI_EXPRESSED_OLD", "IS_CCI_SPECIFIC_OLD", "IS_CCI_SCORE_OLD",
         "IS_DE_LOGFC", "IS_DE_SIGNIFICANT", "DE_DIRECTION",
         "IS_CCI_DETECTED_YOUNG", "IS_CCI_DETECTED_OLD", "IS_CCI_DE",
         "REGULATION")
]

## summary statistics ####

cci_stat_table[, .N]
cci_stat_table[, .N, by = "REGULATION"][, N/sum(N)*100]

regulation_table[
  REGULATION == "NSC",
  sum(N),
  by = c("IS_CCI_DETECTED_YOUNG", "IS_CCI_DETECTED_OLD")
][, V1/sum(V1)*100]

regulation_table[
  REGULATION %in% c("UP", "DOWN"),
  sum(N),
  by = c("IS_CCI_DETECTED_YOUNG", "IS_CCI_DETECTED_OLD")
]

## number of LRI per ER ####

cci_stat_table[, DTER_CELLTYPES := paste(Dataset, Tissue, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, sep = "_")]
NLRI_table <- cci_stat_table[, .N, by = c("Dataset", "DTER_CELLTYPES")]

ggplot(
  data = NLRI_table
) + geom_histogram(
  aes(
    x = N,
    fill = Dataset,
    color = Dataset
  ),
  bins = 100,
  position = "identity",
  alpha = 0.3
)


distr_LRI <- c(NLRI_table$N, rep(0, sum(cci_stat_table[, uniqueN(EMITTER_CELLTYPE), by = c("Dataset", "Tissue")][, V1^2]) - nrow(NLRI_table) ))
mean(distr_LRI)
sd(distr_LRI)
hist(distr_LRI, breaks = 100)

## regulation distribution per tissue ####

regulation_distr <- cci_stat_table[, .N, by = c("Dataset", "Tissue", "REGULATION")]
regulation_distr <- dcast.data.table(
  regulation_distr,
  Dataset + Tissue ~ REGULATION,
  value.var = "N"
)
regulation_distr[is.na(regulation_distr)] <- 0

regulation_distr[
  ,
  total := UP + DOWN + FLAT + NSC
]
regulation_distr[
  ,
  c(
    "UP",
    "DOWN",
    "FLAT",
    "NSC"
  ) :=
    list(
      UP/total,
      DOWN/total,
      FLAT/total,
      NSC/total
    )
]
regulation_distr_long <- melt.data.table(
  regulation_distr,
  id.vars = c("Dataset", "Tissue"),
  measure.vars = c("UP", "DOWN", "FLAT", "NSC"),
  variable.name = "REGULATION",
  value.name = "pct"
)
regulation_distr_long[, id := paste(Dataset, Tissue, sep = "_")]

regulation_distr_long[
  ,
  Dataset := factor(
    Dataset,
    c("TMS FACS (male)", "TMS FACS (female)", "TMS Droplet (male)", "TMS Droplet (female)", "Calico Droplet (male)")
  )
]

## Figure 6 ####

ggplot(
  data = regulation_distr_long,
  aes(
    y = Tissue,
    x = pct,
    fill = REGULATION
  )
) + geom_bar(
  stat = "identity",
  position = "fill",
  alpha = 0.8
) + scale_fill_manual(
  "Age-Regulation",
  values = c("UP" = "red", "DOWN" = "blue", "FLAT" = "green", "NSC" = "grey")
) + facet_wrap(
  ~ Dataset,
  ncol = 5
) + ggplot2::scale_y_discrete(
    limits = sort(
      unique(regulation_distr_long$Tissue),
      decreasing = TRUE
    )
) + xlab(
  "Fraction of CCIs per regulation group"
) + theme(
  text = element_text(size = 34, face = "bold"),
  axis.text.x = element_text(size = 18),
  axis.text.y = element_text(size = 28, face = "bold")
)
#manual save 3000x1400






