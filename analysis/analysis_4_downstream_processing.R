####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - October 2020
##
## Process scDiffCom results and check filtering
## parameters.
##
####################################################
##

## Libraries ####
library(scDiffCom)
library(data.table)
#library(clusterProfiler)
#library(org.Mm.eg.db)

## Specify the directory with scDiffCom results ####
dir_results <- "../data_scAgeCom/scdiffcom_results/"
dir_data_analysis <- "../data_scAgeCom/analysis/"

## Retrieve all the results for each tissue ####

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
names(DATASETS) <- c("calico_nlog")#, "calico_log", "droplet_nlog", "droplet_log", "facs_nlog", "facs_log")

## Let us focus on the not log-transformed data and change some parameters ####

# Based on our test analysis we choose to focus on the non-log-transformed
DATASETS <- DATASETS[grepl("_nlog", names(DATASETS))]

get_reg_pct_behaviour <- function(
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
                skip_ORA = TRUE
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

logfc_cuts <- c(1.1, 1.2, 1.3, 1.4, 1.5)
names(logfc_cuts) <- logfc_cuts
logfc_cuts <- log(logfc_cuts)

score_cuts <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
names(score_cuts) <- score_cuts

pct_reg_dt <- get_reg_pct_behaviour(DATASETS$calico_nlog$Lung, logfc_cuts = logfc_cuts, score_cuts = score_cuts)


ggplot(pct_reg_dt[REGULATION == "FLAT"], aes(x = logfc_cut, y = pct)) + geom_point() +
  facet_grid(vars(score_cut))

ggplot(pct_reg_dt[REGULATION == "FLAT"], aes(x = logfc_cut, y = pct, color = score_cut)) + geom_point()

ggplot(pct_reg_dt[REGULATION == "FLAT"], aes(x = logfc_cut, y = N, color = score_cut)) + geom_point()

ggplot(pct_reg_dt[REGULATION == "UP"], aes(x = logfc_cut, y = N, color = score_cut)) + geom_point()
## Varies some parameters ####

DATASETS_new <- lapply(
  DATASETS,
  function(i) {
    lapply(
      i,
      function(tiss) {
        run_filtering_and_ORA(
          tiss,
          new_cutoff_quantile_score = 0.1
        )
      }
    )
  }
)

Regulation_1.1 <- dcast(
  rbindlist(
    lapply(
      DATASETS,
      function(i) {
        rbindlist(
          lapply(
            i,
            function(tiss) {
              tiss$scdiffcom_dt_filtered[, .N, by =REGULATION_SIMPLE]
            }
          ),
          use.names = TRUE,
          idcol = "Tissue"
        )
      }
    ),
    use.names = TRUE,
    idcol = "dataset"
  ),
  formula = dataset + Tissue ~ REGULATION_SIMPLE,
  value.var = "N"
)

Regulation_1.2 <- dcast(
  rbindlist(
    lapply(
      DATASETS_new,
      function(i) {
        rbindlist(
          lapply(
            i,
            function(tiss) {
              tiss$scdiffcom_dt_filtered[, .N, by =REGULATION_SIMPLE]
            }
          ),
          use.names = TRUE,
          idcol = "Tissue"
        )
      }
    ),
    use.names = TRUE,
    idcol = "dataset"
  ),
  formula = dataset + Tissue ~ REGULATION_SIMPLE,
  value.var = "N"
)

Regulation_comp <- merge.data.table(
  Regulation_1.1,
  Regulation_1.2
)

temp_logfc_comp <- Regulation_comp
temp_quantile_comp <- Regulation_comp
temp_quantile_comp2 <- Regulation_comp

## We keep the initial parameters and create a new object without the raw data ####

DATASETS_light <- lapply(
  DATASETS,
  function(i) {
    lapply(
      i,
      function(tiss) {
        tiss[["scdiffcom_dt_raw"]] <- NA
        return(tiss)
      }
    )
  }
)

saveRDS(DATASETS_light, "../data_scAgeCom/analysis/a4_data_results_nlog.rds")

## Read light dataset ####
DATASETS_light <- readRDS("../data_scAgeCom/analysis/a4_data_results_nlog.rds")

## Add cell families and redo ORA ####

cell_types_dt <- setDT(read.csv("../data_scAgeCom/analysis/scDiffCom_cell_types.csv", stringsAsFactors = FALSE))

lapply(
  DATASETS_light,
  function(i) {
    lapply(
      i,
      function(tiss) {
        tiss$scdiffcom_dt_filtered[
          cell_types_dt, on = "L_CELLTYPE==scDiffCom.cell.type", L_CELL_FAMILY := i.Family...broad
          ]
        tiss$scdiffcom_dt_filtered[
          cell_types_dt, on = "R_CELLTYPE==scDiffCom.cell.type", R_CELL_FAMILY := i.Family...broad
          ]
        tiss$scdiffcom_dt_filtered[, LR_CELL_FAMILY := paste(L_CELL_FAMILY, R_CELL_FAMILY, sep = "_")]
      }
    )
  }
)

DATASETS_light <- lapply(
  DATASETS_light,
  function(i) {
    lapply(
      i,
      function(tiss) {
        scDiffCom::run_ORA(
          scdiffcom_result = tiss,
          verbose = TRUE,
          categories = c("L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPE", "LR_NAME",
                         "L_CELL_FAMILY", "R_CELL_FAMILY", "LR_CELL_FAMILY")
        )
      }
    )
  }
)


## Add senescence associated genes and longevity genes ####

genage_mouse <- setDT(read.csv("../data_scAgeCom/analysis/genage_mouse.tsv", sep = "\t", header = TRUE))

lapply(
  DATASETS_light,
  function(i) {
    lapply(
      i,
      function(tiss) {
        tiss$scdiffcom_dt_filtered[,
                                   GENAGE := ifelse(
                                     LIGAND_1 %in% genage_mouse$Gene.Symbol | LIGAND_2 %in% genage_mouse$Gene.Symbol |
                                       RECEPTOR_1 %in% genage_mouse$Gene.Symbol | RECEPTOR_2 %in% genage_mouse$Gene.Symbol,
                                     "longevity associated",
                                     "not longevity associated"
                                   )
                                   ]
      }
    )
  }
)

DATASETS_light <- lapply(
  DATASETS_light,
  function(i) {
    lapply(
      i,
      function(tiss) {
        scDiffCom::run_ORA(
          scdiffcom_result = tiss,
          verbose = TRUE,
          categories = c("L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPE", "LR_NAME",
                         "L_CELL_FAMILY", "R_CELL_FAMILY", "LR_CELL_FAMILY", "GENAGE", "GO")
        )
      }
    )
  }
)


##################
table(DATASETS_light$calico_nlog$Kidney$scdiffcom_dt_filtered$LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG)

DATASETS_light$calico_nlog$Kidney$scdiffcom_dt_filtered[, .N, by = c("LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG",
                                                                     "LR_DETECTED_AND_SIGNIFICANT_IN_OLD",
                                                                     "REGULATION") ]


reg_counts <- rbindlist(
  lapply(
    DATASETS_light,
    function(i) {
      rbindlist(
        lapply(
          i,
          function(tiss) {
            tiss$scdiffcom_dt_filtered[, .N, by = c("LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG",
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


######


reg_counts_facs2 <- dcast(
  reg_counts_facs,
  formula = Tissue ~ LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG + LR_DETECTED_AND_SIGNIFICANT_IN_OLD + REGULATION,
  value.var = "V1"
)




reg_counts_droplet <- reg_counts$droplet_nlog
reg_counts_droplet[, temp := paste(LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG, LR_DETECTED_AND_SIGNIFICANT_IN_OLD, REGULATION, sep = "_")]
reg_counts_droplet2 <- dcast(
  reg_counts_droplet,
  formula = Tissue ~ LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG + LR_DETECTED_AND_SIGNIFICANT_IN_OLD + REGULATION,
  value.var = "V1"
)

ggplot(reg_counts_droplet, aes(x= temp, y = V1)) + geom_boxplot()




library(ontoProc)
onto_go <- getGeneOnto()
onto_go_names <- onto_go$name

go_id_name_dt <- data.table(
  name = onto_go_names,
  ID = names(onto_go_names)
)

DATASET_GO <- lapply(
  DATASETS_light,
  function(i) {
    lapply(
      i,
      function(tiss) {
        temp <- scDiffCom::run_ORA(
          scdiffcom_result = tiss,
          verbose = TRUE,
          categories = c("L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPE", "LR_NAME",
                         "GO"),
          logfc_threshold = log(1.1)
        )
      }
    )
  }
)


test <- scDiffCom::run_ORA(
  scdiffcom_result = DATASETS_light$droplet_nlog$Limb_Muscle,
  verbose = TRUE,
  categories = c("L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPE", "LR_NAME",
                 "GO"),
  logfc_threshold = log(1.1)
)



table(test$scdiffcom_dt_filtered$REGULATION_SIMPLE)
hist(test$scdiffcom_dt_filtered$LOGFC, breaks = 100)
ggplot(test$scdiffcom_dt_filtered, aes(x = LOGFC, y = -log10(BH_PVAL_DIFF+1E-4))) +
  geom_point()

test_ora <- test$ORA

test_ora <- merge.data.table(
  test_ora,
  go_id_name_dt,
  by.x = "Value",
  by.y = "ID",
  all.x = TRUE
)

test_ora_goInter <- test_ora[Category == "GO_intersection"]
test_ora_goInter_up <- test_ora_goInter[pval_adjusted_UP <= 0.05 & OR_UP >= 1]
test_ora_goInter_down <- test_ora_goInter[pval_adjusted_DOWN <= 0.05 & OR_DOWN >= 1]

unique(test_ora$Category)

ggplot(test_ora[pval_adjusted_UP <= 0.05 & OR_UP >= 1 & Category == "LR_NAME"], 
       aes(x = log10(OR_UP), y = -log10(pval_adjusted_UP), label = Value)) +
  geom_point() + 
  geom_text(aes(label = Value))

ggplot(test_ora[pval_adjusted_DOWN <= 0.05 & OR_DOWN >= 1 & Category == "LR_NAME"], 
       aes(x = log10(OR_DOWN), y = -log10(pval_adjusted_DOWN), label = Value)) +
  geom_point() + 
  geom_text(aes(label = Value))

ggplot(test_ora[pval_adjusted_UP <= 0.05 & OR_UP >= 1 & Category == "LR_CELLTYPE"], 
       aes(x = log10(OR_UP), y = -log10(pval_adjusted_UP), label = Value)) +
  geom_point() + 
  geom_text(aes(label = Value))

ggplot(test_ora[pval_adjusted_DOWN <= 0.05 & OR_DOWN >= 1 & Category == "LR_CELLTYPE"], 
       aes(x = log10(OR_DOWN), y = -log10(pval_adjusted_DOWN), label = Value)) +
  geom_point() + 
  geom_text(aes(label = Value))

intersect(
  test_ora[pval_adjusted_UP <= 0.05 & OR_UP >= 1 & Category == "LR_CELLTYPE"]$Value,
  test_ora[pval_adjusted_DOWN <= 0.05 & OR_DOWN >= 1 & Category == "LR_CELLTYPE"]$Value
)



ggplot(test_ora_goInter_up, aes(x = log10(OR_UP), y = -log10(pval_adjusted_UP), label = name )) + geom_point() +
  geom_text(aes(label=ifelse(OR_UP >= 1, name, "")), position=position_jitter(width=0.1,height=0.1))

ggplot(test_ora_goInter_down, aes(x = log10(OR_DOWN), y = -log10(pval_adjusted_DOWN), label = name )) + geom_point() +
  geom_text(aes(label=ifelse(OR_DOWN >= 1, name, "")), position=position_jitter(width=0.1,height=0.11))


intersect(test_ora_goInter_up$name, test_ora_goInter_down$name)

test_ora_goInter_up[order(pval_adjusted_UP)][1:50]$name
onto_plot2(onto_go, exclude_descendants(onto_go, c("GO:0003674" ,"GO:0005575") , 
                                        test_ora_goInter_up[order(pval_adjusted_UP)]$Value))


onto_plot2(onto_go, remove_links(onto_go, exclude_descendants(onto_go, c("GO:0003674" ,"GO:0005575") , 
                                                              test_ora_goInter_up[order(pval_adjusted_UP)]$Value), hard = TRUE))




test_ora_goInter_down[order(pval_adjusted_DOWN)][1:50]$Value
onto_plot2(onto_go, exclude_descendants(onto_go, c("GO:0003674" ,"GO:0005575") , 
                                        test_ora_goInter_down[order(pval_adjusted_DOWN)]$Value))

intersect(test_ora_goInter_down$name, test_ora_goInter_up$name)

test_ora_goUnion <- test_ora[Category == "GO_UNION"]

test_ora_goUnion_up <- test_ora_goUnion[pval_adjusted_UP <= 0.05 & OR_UP >= 1]
test_ora_goUnion_up[order(pval_adjusted_UP)][1:50]$name
onto_plot2(onto_go, exclude_descendants(onto_go, c("GO:0003674" ,"GO:0005575") , 
                                        test_ora_goUnion_up[order(pval_adjusted_UP)][1:50]$Value))

test_ora_goUnion_down <- test_ora_goUnion[pval_adjusted_DOWN <= 0.05 & OR_DOWN >= 1]
test_ora_goUnion_down[order(pval_adjusted_DOWN)][1:50]$Value
onto_plot2(onto_go, exclude_descendants(onto_go, c("GO:0003674" ,"GO:0005575") , 
                                        test_ora_goUnion_down[order(pval_adjusted_DOWN)][1:50]$Value))

intersect(test_ora_goUnion_down$name, test_ora_goUnion_up$name)


############


test3_up <- test2[Category == "GO_intersection" & pval_adjusted_UP <= 0.05 & OR_UP >= 4]
test3_down <- test2[Category == "GO_intersection" & pval_adjusted_DOWN <= 0.05 & OR_DOWN >= 4]
test4_up <- as.data.table(sort(onto_go_names[names(onto_go_names) %in% test3_up$Value]))
test4_down <- as.data.table(sort(onto_go_names[names(onto_go_names) %in% test3_down$Value]))

test3_up <- merge.data.table(
  test3_up,
  go_id_name_dt,
  by.x = "Value",
  by.y = "ID",
  all.x = TRUE
)

test3_down <- merge.data.table(
  test3_down,
  go_id_name_dt,
  by.x = "Value",
  by.y = "ID",
  all.x = TRUE
)

intersect(test3_up$Value, test3_down$Value)

test5_up <- names(onto_go_names[names(onto_go_names) %in% test3_up$Value])
test5_down <- names(onto_go_names[names(onto_go_names) %in% test3_down$Value])


onto_plot2(onto_go, exclude_descendants(onto_go, c("GO:0003674" ,"GO:0005575") , test5_up))

exclude_descendants(onto_go, c("GO:0003674" ,"GO:0005575") , test5_up)
intersection_with_descendants(onto_go, "", test5_up)



## Naive GO analysis ####


universe_LR <- unique(unlist(LR6db$LR6db_curated[, c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
universe_LR <- universe_LR[!is.na(universe_LR)]

scdiffcom_all_filtered <- rbindlist(
  lapply(
    DATASETS_light,
    function(i) {
      rbindlist(
        lapply(
          i,
          function(tiss) {
            tiss$scdiffcom_dt_filtered
          }
        ),
        use.names = TRUE,
        idcol = "TISSUE",
        fill = TRUE
      )
    }
  ),
  use.names = TRUE,
  idcol = "DATASET"
)
scdiffcom_all_filtered[, c("L_TCT", "R_TCT", "LOG2FC") := list(
  paste(TISSUE, L_CELLTYPE, sep = ": "),
  paste(TISSUE, R_CELLTYPE, sep = ": "),
  LOGFC*log2(exp(1))
)]


temp <- DATASETS_light$facs_nlog$Lung$scdiffcom_dt_filtered
temp <- scdiffcom_all_filtered[DATASET == "facs_nlog"]
up_genes <- unique(unlist(temp[REGULATION_SIMPLE == "UP" & LOGFC > log(4), c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
up_genes <- up_genes[!is.na(up_genes)]

down_genes <- unique(unlist(temp[REGULATION_SIMPLE == "DOWN" & LOGFC < -log(4), c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
down_genes <- down_genes[!is.na(down_genes)]

intersect(up_genes, down_genes)

ego_up <- enrichGO(
  gene          = up_genes,
  universe      = universe_LR,
  keyType       = "SYMBOL",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = FALSE
)

ego_down <- enrichGO(
  gene          = down_genes,
  universe      = universe_LR,
  keyType       = "SYMBOL",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = FALSE
)

dotplot(ego_up, showCategory = 10)
dotplot(ego_down, showCategory = 10)

ego_MF_strong <- lapply(
  c(1,2,3,4,5,7,8,9),
  function(i) {
    temp <- DATASETS[[i]]$results
    up_genes <- unique(unlist(temp[REGULATION == "UP" & LOGFC > log(2), c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
    up_genes <- up_genes[!is.na(up_genes)]
    down_genes <- unique(unlist(temp[REGULATION == "DOWN" & LOGFC < -log(2), c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
    down_genes <- down_genes[!is.na(down_genes)]
    ego_up <- enrichGO(gene       = up_genes,
                       universe      = universe_LR,
                       keyType       = "SYMBOL",
                       OrgDb         = org.Mm.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = FALSE)
    ego_down <- enrichGO(gene       = down_genes,
                         universe      = universe_LR,
                         keyType       = "SYMBOL",
                         OrgDb         = org.Mm.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = FALSE)
    return(list(ego_up, ego_down))
  }
  
)

plot_go_MF_droplet <- cowplot::plot_grid(
  plotlist = list(
    dotplot(ego_MF_strong[[3]][[1]], showCategory = 10),
    dotplot(ego_MF_strong[[3]][[2]], showCategory = 10),
    dotplot(ego_MF_strong[[4]][[1]], showCategory = 10),
    dotplot(ego_MF_strong[[4]][[2]], showCategory = 10),
    dotplot(ego_MF_strong[[5]][[1]], showCategory = 10),
    dotplot(ego_MF_strong[[5]][[2]], showCategory = 10)
  ),
  ncol = 2,
  align = "v",
  labels = c(paste(names(DATASETS)[3], c("UP", "DOWN")),
             paste(names(DATASETS)[4], c("UP", "DOWN")),
             paste(names(DATASETS)[5], c("UP", "DOWN")))
)
plot_go_MF_droplet
ggsave(filename = paste0(dir_data_analysis, "a4_plot_go_MF_droplet.png"),
       plot = plot_go_MF_droplet, scale = 2)


cowplot::plot_grid(
  plotlist = list(
    dotplot(ego_MF_strong[[6]][[1]], showCategory = 10),
    dotplot(ego_MF_strong[[6]][[2]], showCategory = 10),
    dotplot(ego_MF_strong[[7]][[1]], showCategory = 10),
    dotplot(ego_MF_strong[[7]][[2]], showCategory = 10),
    dotplot(ego_MF_strong[[8]][[1]], showCategory = 10),
    dotplot(ego_MF_strong[[8]][[2]], showCategory = 10)
  ),
  ncol = 2,
  align = "v"
)

## Local vs global analysis ####
facs <- scdiffcom_all_filtered[DATASET == "facs_nlog"]
facs[, LR_TCT := paste(L_TCT, R_TCT, sep = "_")]

length(unique(facs$LR_TCT))
length(unique(facs$TISSUE))

facs_up_all <- facs[REGULATION_SIMPLE == "UP" & LOGFC >= log(1.1), .N, by = c("LR_NAME", "TISSUE")][, .N, by = LR_NAME]
facs_down_all <- facs[REGULATION_SIMPLE == "DOWN" & LOGFC <= -log(1.1), .N, by = c("LR_NAME", "TISSUE")][, .N, by = LR_NAME]

facs_up_all_top <- facs_up_all[N >= 15][order(-N)]
colnames(facs_up_all_top) <- c("LR_GENES", "Tissue_Counts")
facs_down_all_top <- facs_down_all[N>=18][order(-N)]
colnames(facs_down_all_top) <- c("LR_GENES", "Tissue_Counts")


plot_facs_up_all_top <- tableGrob(facs_up_all_top, rows = NULL)
grid.newpage()
grid.draw(plot_facs_up_all_top)

ggsave(filename = paste0(dir_data_analysis, "plot_facs_all_up.png"),
       plot = plot_facs_up_all_top, scale = 1.5)

plot_facs_down_all_top <- tableGrob(facs_down_all_top, rows = NULL)
grid.newpage()
grid.draw(plot_facs_down_all_top)

ggsave(filename = paste0(dir_data_analysis, "plot_facs_all_down.png"),
       plot = plot_facs_down_all_top, scale = 2)


test3 <- facs[REGULATION_SIMPLE == "UP" & LOGFC >= log(1.1), .N, by = c("LR_NAME", "LR_TCT")][, .N, by = LR_NAME]
test5 <- facs[REGULATION_SIMPLE == "DOWN" & LOGFC <= -log(1.1), .N, by = c("LR_NAME", "LR_TCT")][, .N, by = LR_NAME]




## Test of GENAGE ###########

test <- DATASETS_light$facs_nlog$BAT$scdiffcom_dt_filtered

ftable(test$GENAGE, test$REGULATION_SIMPLE)

test <- sort(unique(unlist(LR6db$LR6db_curated[, c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")])))

test2 <- sort(unique(unlist(DATASETS_light$droplet_nlog$Bladder$scdiffcom_dt_filtered[, c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")])))

intersect(test2, genage_mouse$Gene.Symbol)

test <- DATASETS_light$facs_nlog$Lung$ORA[OR_UP >=1 & pval_adjusted_UP <= 0.05 & Category == "GENAGE"]

test2 <- DATASETS_light$facs_nlog$Lung$scdiffcom_dt_filtered[GENAGE == "longevity associated" & REGULATION_SIMPLE == "UP",]

## Test fpm analysis ####

library(arules)
library(glue)
source("src/fpm.R")

test <- analyze_FreqItemSets(
  data = DATASETS_light$facs_nlog$Lung$scdiffcom_dt_filtered,
  support = 0.001,
  confidence = 0.1
)

test2 <- test$sub1
test3 <- test$sub2

test_d <- DATASETS_light$calico_nlog$Kidney$scdiffcom_dt_filtered

DATASETS_light$calico_nlog$Kidney$scdiffcom_dt_filtered$DIFFERENTIALLY_EXPRESSED

## Test heatmaps ####


library(purrr)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(igraph)

source("src/utils.R")
source("src/visualizations.R")

test_heat <- DATASETS_light$facs_nlog$`Brain_Non-Myeloid`
test_heat = add_dummy_tissue(test_heat)

str(test_heat, max.level = 1)

create_heatmaps_dirs = function(base_dir, dir_names){
  # Create paths
  heatmaps_dirs = map_chr(
    dir_names,
    ~ file.path(base_dir, .x)
  )  
  # Create dirs at paths
  map(
    heatmaps_dirs,
    ~ dir.create(.x, showWarnings=FALSE)
  )
  heatmaps_dirs = as.list(heatmaps_dirs)
  names(heatmaps_dirs) = dir_names
  return(heatmaps_dirs)
}

construct_heatmaps = function(ora_dt, dir){
  dt_ctypes = get_celltypes_enrichment(ora_dt)
  build_heatmaps(dt_ctypes, dir)
}

HEATMAPS_BASE_DIR = "../data_scAgeCom/heatmaps"
heatmaps_dirs = create_heatmaps_dirs(HEATMAPS_BASE_DIR, c('test_dataset'))
construct_heatmaps(test_heat$ORA, heatmaps_dirs$test_dataset)


## Test cell bias#####
dt = test_heat$scdiffcom_dt_filtered

(ggplot(data=dt, aes(L_NCELLS_YOUNG, R_NCELLS_YOUNG))
  + geom_point(aes(color=DIFFERENTIALLY_EXPRESSED))
  + ggtitle('N cell bias: young')
)

(ggplot(data=dt, aes(L_NCELLS_OLD, R_NCELLS_OLD))
  + geom_point(aes(color=DIFFERENTIALLY_EXPRESSED))
  + ggtitle('N cell bias: old')
)

(ggplot(data=dt, aes(L_NCELLS_YOUNG, L_NCELLS_OLD))
  + geom_point(aes(color=DIFFERENTIALLY_EXPRESSED))
  + ggtitle('N cell bias: ligand cells young vs old')
)

(ggplot(data=dt, aes(R_NCELLS_YOUNG, R_NCELLS_OLD))
  + geom_point(aes(color=DIFFERENTIALLY_EXPRESSED))
  + ggtitle('N cell bias: receptor cells young vs old')
)


## Test igraph plots #####

dir <- "../data_scAgeCom/analysis/Graphs"
if ( !dir.exists(dir) ) {dir.create(dir)}

analyze_Graph(
  test_heat$ORA,
  test_heat$scdiffcom_dt_filtered,
  "DummyTissue",
  config=GRAPH_CONFIG,
  use_adjpval=TRUE,
  dir=dir,
  analysis_name="DummyTest")


## TOP CCI ####

n_top <- 1000

top_LR_up_FACS <- unique(CCI_tmsFACS_age[order(-LOGFC)][1:n_top]$LR_NAME)
top_LR_up_FACS

top_genes_up_FACS <- sort(unique(unlist(CCI_tmsFACS_age[order(-LOGFC)][1:n_top, c(paste0("LIGAND_", 1:2), paste0("RECEPTOR_", 1:3) )])))
top_genes_up_FACS

top_LR_down_FACS <- unique(CCI_tmsFACS_age[order(LOGFC)][1:n_top]$LR_NAME)
top_LR_down_FACS

top_genes_down_FACS <- sort(unique(unlist(CCI_tmsFACS_age[order(LOGFC)][1:n_top, c(paste0("LIGAND_", 1:2), paste0("RECEPTOR_", 1:3) )])))
top_genes_down_FACS


## global ORA ####

ORA_all <- readRDS("../data_scAgeCom/analysis/analysis_4_data_ora.rds")
ORA_tmsFACS <- ORA_all$tms_facs

top_ORA_up_FACS <- ORA_tmsFACS[Tissue == "All" & Category == "LR_NAME" & pval_UP <= 0.05 & OR_UP >= 1][order(pval_UP)]$Value
top_ORA_up_FACS

top_ORA_down_FACS <- ORA_tmsFACS[Tissue == "All" & Category == "LR_NAME" & pval_DOWN <= 0.05 & OR_DOWN >= 1][order(pval_DOWN)]$Value
top_ORA_down_FACS










## General statistics per (dataset-) tissue ####

CCI_counts <- FULL_RESULTS_AGE[, .N, by = c("DATASET", "TISSUE", "CASE_TYPE")]
CCI_counts_full <- FULL_RESULTS_AGE[, .N, by = c("DATASET", "TISSUE")]
CCI_counts_full2 <- FULL_RESULTS_AGE_2[, .N, by = c("DATASET", "TISSUE")]

CCI_counts_dc <- dcast.data.table(
  CCI_counts,
  formula = DATASET + TISSUE ~ REGULATION,
  value.var = "N"
)

test <- dcast.data.table(
  CCI_counts,
  formula = TISSUE  ~ DATASET + REGULATION,
  value.var = "N"
)

FULL_RESULTS_AGE[, CCIT := paste(TISSUE, LR_CELLTYPE, LR_NAME, sep = "_")]

head(sort(table(FULL_RESULTS_AGE$CCI), decreasing = TRUE))
head(sort(table(FULL_RESULTS_AGE[ REGULATION == "UP"]$CCIT), decreasing = TRUE), 50)
head(sort(table(FULL_RESULTS_AGE[ REGULATION == "DOWN"]$CCIT), decreasing = TRUE), 50)


test <- ORA_RESULTS$tms_facs
test2 <- test[Category == "LR_NAME", c("Tissue", "Category", "Value", "pval_UP", "OR_UP")]
test3<- test[Category == "LR_NAME", c("Tissue", "Category", "Value", "pval_DOWN", "OR_DOWN")]

LRdbcur <- LR6db$LR6db_curated
LRdbcur[LIGAND_1 == "Crlf2"]






