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
## scRNA-seq data preparation
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
  text = element_text(size = 30),
  axis.text.x = element_text(size = 16)
)


###########################



stat_table[, N_LRI_per_ER := `Total Cell-Cell Interactions`/(`Total Cell Types`)^2]



test[, frac_up := `UP CCIs`/`Total Cell-Cell Interactions`*100]
test[, frac_down := `Down CCIs`/`Total Cell-Cell Interactions`*100]
test[, frac_flat := `Flat CCIs`/`Total Cell-Cell Interactions`*100]

test2 <- melt.data.table(
  test[, c(1,2,10, 11,12)],
  id.vars = c("Tissue", "Dataset"),
  value.vars = c("frac_up", "frac_down", "fracs_flat")
)

ggplot(test2, aes(x = Dataset, y = value, color = variable)) + geom_boxplot()
ggplot(test2, aes(x = Dataset, y = frac_down)) + geom_boxplot()



ggplot(test, aes(x = frac_up, y = frac_down)) + geom_point()

ggplot(test, aes(x = `Total Cell Types`, y = N_LRI_per_ER)) + geom_point()


sum(test$`Total Cell-Cell Interactions`)
mean(test$N_LRI_per_ER)
sd(test$N_LRI_per_ER)

hist(test$N_LRI_per_ER, breaks = 20)


hist(test$frac_down, breaks = 20)
hist(test$frac_up, breaks = 20)

test3 <- copy(scAgeCom_shiny_data$CCI_table)
test3[, ERI := paste(`Emitter Cell Type`, `Receiver Cell Type`, sep = "_")]
test3[, ID := paste(Dataset, Tissue, sep = "_")]
test3[, ID2 := paste(ID, ERI, sep = "_")]
test3[, .N, by = "ID2"]

setkey(test3, ID, ERI)



test4 <- test3[, .N, by = c("Tissue", "Dataset", "ERI")]
test5 <- test3[, .N, by = c("Tissue", "Dataset", "ERI", "Age Regulation")]

test4[
  test5[`Age Regulation` == "UP"],
  on = c("Tissue", "Dataset", "ERI"),
  N_UP := i.N
]
test4[
  test5[`Age Regulation` == "DOWN"],
  on = c("Tissue", "Dataset", "ERI"),
  N_DOWN := i.N
]
test4[
  test5[`Age Regulation` == "FLAT"],
  on = c("Tissue", "Dataset", "ERI"),
  N_FLAT := i.N
]
test4[
  test5[`Age Regulation` == "NSC"],
  on = c("Tissue", "Dataset", "ERI"),
  N_NSC := i.N
]

test5 <- melt.data.table(
  test4,
  id.vars = c("Tissue", "Dataset", "ERI"),
  measure.vars = c("N_UP", "N_DOWN", "N_FLAT", "N_NSC")
)

test4[is.na(test4)] <- 0

test4[, mean(N)]
test4[, sd(N)]
test4[, median(N)]
hist(test4$N, breaks = 100)

test4[, mean(N_UP)]
test4[, sd(N_UP)]
hist(test4$N_UP, breaks = 100)

test4[, mean(N_DOWN)]
test4[, sd(N_DOWN)]
hist(test4$N_DOWN, breaks = 100)

test4[, mean(N_FLAT)]
test4[, sd(N_FLAT)]
hist(test4$N_FLAT, breaks = 100)





p <- ggplot(
  test4
) + geom_boxplot(
  aes(
    x = N,
    y = Tissue#,
    #color = Dataset
  )
)

ggplot(
  test4
) + geom_boxplot(
  aes(
    x = N_UP,
    y = Tissue#,
    #color = Dataset
  )
)

ggplot(
  test4
) + geom_histogram(
  aes(
    x = N
  ),
  bins = 100
)

ggplot(
  test5
) + geom_boxplot(
  aes(
    y = value,
    x = Tissue,
    color = variable
  )
) + facet_wrap(
  Dataset ~ .
)





