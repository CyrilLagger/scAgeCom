####################################################
##
## Project: scAgeCom
## 
## October 2020
##
## anais.equey@etu.univ-amu.fr 
## cyril.lagger@liverpool.ac.uk 
##
## 
##
####################################################
##

## Load libraries ####

library(scDiffCom)
library(data.table)
library(future)
library(future.apply)

## Analysis data path ####
dir_data_analysis <- "../data_scAgeCom/analysis/"

##Load scDiffCom results ####

scdiffcom_results <- readRDS("../data_scAgeCom/analysis/a4_data_results_nlog.rds")

## 
#compute score and ora for a specific go terms
#needs to compute counts and score for a term that is not in a column of the data.frame

#take a go term for example 

#take go term
#find associated LR interaction
#create a column in data.table with TRUE/FALSE
#perform score and ORA on this column

scDiffCom::run_ORA()

## 11. Compute score and ORA for each GO/family term ####



GO_terms_keep <- GO_terms[GO_terms %in% GO_terms_in_data & GO_terms != ""]
family_terms_keep <- family_terms[family_terms %in% family_terms_in_data & family_terms != ""]

scdiffcom_results_bind <- rbindlist(
  lapply(
    scdiffcom_results,
    function(i) {
      rbindlist(
        lapply(
          i,
          function(tiss) {
            tiss$scdiffcom_dt_filtered
          }
        ),
        use.names = TRUE,
        fill = TRUE,
        idcol = "TISSUE"
      )
    }
  ),
  use.names = TRUE,
  idcol = "DATASET"
)

max_logfc <- max(scdiffcom_results_bind[is.finite(LOGFC)]$LOGFC)
min_logfc <- min(scdiffcom_results_bind[is.finite(LOGFC)]$LOGFC)

scdiffcom_results_bind[, LOGFC_reg := ifelse(
  is.finite(LOGFC),
  LOGFC,
  ifelse(
    LOGFC > 0,
    max_logfc,
    min_logfc
  )
)]
max(scdiffcom_results_bind$LOGFC_reg)

test <- lapply(
  unique(scdiffcom_results_bind$DATASET),
  function(i) {
    print(i)
    l1 <- lapply(
      unique(scdiffcom_results_bind[DATASET == i]$TISSUE),
      function(j) {
        print(j)
        temp <- scdiffcom_results_bind[DATASET == i & TISSUE == j]
        sort(sapply(
          unique(scdiffcom_results_bind$LR_NAME),
          function(k) {
            temp[REGULATION_SIMPLE == "UP"][
              LR_NAME == k
              ][, sum(LOGFC_reg)]/temp[REGULATION_SIMPLE == "UP",sum(LOGFC_reg)]
          }
        ),
        decreasing = TRUE
        )
      }
    )
    names(l1) <- unique(scdiffcom_results_bind[DATASET == i]$TISSUE)
    return(l1)
  }
)

temp <- scdiffcom_results_bind[DATASET == "droplet_nlog"]

all_LR <- unique(temp$LR_NAME)
all_CC <- unique(temp$LR_CELLTYPE)

test_up <- sapply(all_CC, function(i) {
  temp[REGULATION_SIMPLE == "UP" & LR_CELLTYPE == i][,sum(LOGFC_reg)]/temp[REGULATION_SIMPLE == "UP",sum(LOGFC_reg)]
}
)

test_down <- sapply(all_CC, function(i) {
  temp[REGULATION_SIMPLE == "DOWN" & LR_CELLTYPE == i][,sum(LOGFC_reg)]/temp[REGULATION_SIMPLE == "DOWN",sum(LOGFC_reg)]
}
)

head(sort(test_up, decreasing = TRUE))
head(sort(test_down, decreasing = FALSE))

hist(log10(test_up), breaks = 50)
hist(log10(test_down), breaks = 50)

identical(names(test_up), names(test_down))

test_diff <- test_up - test_down

head(sort(test_diff, decreasing = TRUE), 20)

head(sort(test_diff, decreasing = FALSE), 20)

rank_up <- data.table(
  LR_genes = names(sort(test_diff, decreasing = TRUE)),
  Rank = 1:length(test_diff)
)

rank_down <- data.table(
  LR_genes = names(sort(test_diff, decreasing = FALSE)),
  Rank = 1:length(test_diff)
)

rank_up_top <- rank_up[Rank <= 10]
rank_down_top <- rank_down[Rank <= 10]

colnames(rank_up_top) <- c("Cell-types", "Rank")
colnames(rank_down_top) <- c("Cell-types", "Rank")



plot_droplet_up_rank_top <- tableGrob(rank_up_top, rows = NULL)
grid.newpage()
grid.draw(plot_droplet_up_rank_top)

ggsave(filename = paste0(dir_data_analysis, "plot_droplet_rank_up_CC.png"),
       plot = plot_droplet_up_rank_top, scale = 1.5)

plot_droplet_down_rank_top <- tableGrob(rank_down_top, rows = NULL)
grid.newpage()
grid.draw(plot_droplet_down_rank_top)

ggsave(filename = paste0(dir_data_analysis, "plot_droplet_rank_down_CC.png"),
       plot = plot_droplet_down_rank_top, scale = 1.5)




##
future::plan(multiprocess)
options(future.globals.maxSize = 64 * 1000 ^ 3)



start_time <- Sys.time()
scores_go_union_up <- lapply(
  unique(scdiffcom_results_bind$DATASET),
  function(i) {
    print(i)
    l1 <- lapply(
      unique(scdiffcom_results_bind[DATASET == i]$TISSUE),
      function(j) {
        print(j)
        temp <- scdiffcom_results_bind[DATASET == i & TISSUE == j]
        sort(future_sapply(
          GO_terms_keep,
          function(k) {
            temp[REGULATION_SIMPLE == "UP"][
              stri_detect_fixed(GO_union, k)
              #grepl(i, GO_union, fixed = TRUE)
              ][, sum(LOGFC_reg)]/temp[REGULATION_SIMPLE == "UP",sum(LOGFC_reg)]
          }
        ),
        decreasing = TRUE
        )
      }
    )
    names(l1) <- unique(scdiffcom_results_bind[DATASET == i]$TISSUE)
    return(l1)
  }
)
names(scores_go_union_up) <- unique(scdiffcom_results_bind$DATASET)
end_time <- Sys.time()
end_time - start_time

## 12. GO score analysis ######

go_scores_up <- readRDS("../data_scAgeCom/analysis/scores_go_union_up.rds")
go_scores_down <- readRDS("../data_scAgeCom/analysis/scores_go_union_down.rds")

test <- lapply(
  go_scores_up,
  function(i) {
    rn <- names(i[[1]])
    l <- i
    l$GO <- rn
    setDT(l)
  }
)

test2 <- lapply(
  go_scores_down,
  function(i) {
    rn <- names(i[[1]])
    l <- i
    l$GO <- rn
    setDT(l)
  }
)

test3 <- merge.data.table(
  test$facs_nlog,
  test2$facs_nlog,
  by = "GO"
)

test4 <- test3[, c("GO", "Diaphragm.x", "Diaphragm.y")]
test4[, diff := `Diaphragm.x` - `Diaphragm.y`]

test5 <- test3[, c("GO", "GAT.x", "GAT.y")]
test5[, diff := `GAT.x` - `GAT.y`]

####################


go_ora <- lapply(
  GO_terms_keep,
  function(go_t) {
    CCI_all <- lapply(CCI_all, function(dataset) {
      dataset$scdiffcom_dt_filtered[, temp_col := ifelse(grepl(go_t, GO_intersection, fixed = TRUE), TRUE, FALSE)]
      return(dataset)
    })
    test <- scDiffCom::run_ORA(
      CCI_all[[1]],
      verbose = TRUE,
      categories = "temp_col"
    )
    test <- test$ORA
    test <- test[Value == TRUE]
    test[, name := go_t]
  }
)

go_scores <- rbindlist(lapply(
  GO_terms_keep,
  function(go_t) {
    temp <- CCI_all[[1]]$scdiffcom_dt_filtered[grepl(go_t, GO_intersection)]
    res <- data.table(
      GO = go_t,
      sum_up = sum(temp[REGULATION_SIMPLE == "UP"][["LOGFC"]]),
      sum_down = sum(temp[REGULATION_SIMPLE == "DOWN"][["LOGFC"]]),
      sum_flat = sum(temp[REGULATION_SIMPLE == "FLAT"][["LOGFC"]])
    )
    return(res)
  }
))

go_scores[, diff := sum_up + sum_down]

test_ora <- rbindlist(go_ora)

test_ora_up <- test_ora[pval_adjusted_UP <= 0.05 & OR_UP >= 1, c("pval_adjusted_UP", "OR_UP", "name")][order(pval_adjusted_UP)]
test_ora_down <- test_ora[pval_adjusted_DOWN <= 0.05 & OR_DOWN >= 1, c("pval_adjusted_DOWN", "OR_DOWN", "name")][order(pval_adjusted_DOWN)]


ggplot(test_ora, aes(x = log2(OR_UP), y = -log10(pval_adjusted_UP))) + geom_point()

#the part below is very slow and not working properly, I have to think about it

GO_number <- pbsapply(GO_terms[1:500], function(i) {
  sum(grepl(i, CCI_tmsFACS_age$GO_union, fixed = TRUE))
})

hist(GO_number)

head(sort(GO_number), 100)

GO_logfc <- pbsapply(GO_terms, function(i) {
  mean(CCI_tmsFACS_age[grepl(i, CCI_tmsFACS_age$GO_intersection, fixed = TRUE)]$LOGFC)
})

hist(GO_logfc)

head(sort(GO_logfc, decreasing = TRUE))
head(sort(GO_logfc, decreasing = FALSE))


hist(CCI_tmsFACS_age[grepl("zinc ion binding", CCI_tmsFACS_age$GO_union)]$LOGFC, breaks = 100)

test <- CCI_all$calico
test$GO_intersection[[12345]]
test$GO_union[[1234]]

identical(test$GO_intersection[[1]], test$GO_intersection[[25]])

######
library(data.table)
library(scDiffCom)
LR <- LR6db$LR6db_curated

LR_genes <- unique(unlist(LR[,c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
LR_genes <- LR_genes[!is.na(LR_genes)]

genage <- read.csv("../../../../../genage_models_export.tsv", header = TRUE, sep = "\t")

LR_long_genes <- genage$Gene.Symbol[genage$Gene.Symbol %in% LR_genes]


LR_genage <- LR[LIGAND_1 %in% LR_long_genes | LIGAND_2 %in% LR_long_genes |
                  RECEPTOR_1 %in% LR_long_genes | RECEPTOR_2 %in% LR_long_genes | RECEPTOR_3 %in% LR_long_genes]

cellage <- readRDS("../../P2_scAgingNetworks/scSenescence/senescence_2.rds")
setDT(cellage)
cellage <- cellage[database == "cellage"]
cellage_ortho <- get_orthologs(cellage$gene, input_species = "human")

LR_cellage_genes <- cellage_ortho$mouse_symbol[cellage_ortho$mouse_symbol %in% LR_genes]

#####

res_facs <- DATASETS$`TMS FACS Male`$results


test <- res_facs[REGULATION %in% c("UP"), .N, by = c("LR_NAME", "TISSUE", "REGULATION")]

test2 <- test[, .N, by = c("LR_NAME")]

test3 <- dcast.data.table(
  test,
  formula = LR_NAME ~ TISSUE,
  value.var = "N"
)


testd <- res_facs[REGULATION %in% c("DOWN"), .N, by = c("LR_NAME", "TISSUE", "REGULATION")]

testd2 <- testd[, .N, by = c("LR_NAME")]

testd3 <- dcast.data.table(
  testd,
  formula = LR_NAME ~ TISSUE,
  value.var = "N"
)

library(parallel)
library(foreach)

data <- 1:1e9
data_list <- list("1" = data,
                  "2" = data,
                  "3" = data,
                  "4" = data)
detectCores()
cl <- parallel::makeCluster(detectCores())
# Activate cluster for foreach library
doParallel::registerDoParallel(cl)
time_foreach <- system.time({
  r <- foreach::foreach(i = 1:length(data_list),
                        .combine = rbind) %dopar% {
                          mean(data_list[[i]])
                        }
})
time_foreach[3]
# Stop cluster to free up resources
parallel::stopCluster(cl)


## xxx. On how to extract group of genes with specific terms ####

#Example: let's retrieve all the go terms with "immune" and all the genes with "immune"
go_immune <- LR_go_frequency[grepl("immune", name_1006)]
LR_genes_immune_1 <- LR_genes_go[grepl("immune", GO_terms)]
LR_genes_immune_2 <- LR_genes_go[grepl("immune response", GO_terms)]

#inflammation
go_inflam <- LR_go_frequency[grepl("inflam", name_1006)]
LR_genes_inflam_1 <- LR_genes_go[grepl("inflam", GO_terms)]

#collagen
go_colla <- LR_go_frequency[grepl("collagen", name_1006)]
LR_genes_colla_1 <- LR_genes_go[grepl("collagen", GO_terms)]


## Some test
LR_go_frequency_short <- LR_go_frequency[N >= 10]
LR_genes_go_short <- LR_genes_go[name_1006 %in% LR_go_frequency_short$name_1006]

LR_genes_go_combined_short <- as.data.table(
  sapply(unique(LR_genes_go_short$mgi_symbol), function(i) {
    temp <- LR_genes_go_short[mgi_symbol == i]
    res <- paste0(sort(temp$name_1006), collapse = ",")
    return(res)
  },
  USE.NAMES = TRUE),
  keep.rownames = TRUE
)
colnames(LR_genes_go_combined_short) <- c("gene", "GO_terms")

## Some test with omnipathR
library(OmnipathR)
get_annotation_resources()
test <- import_omnipath_annotations(proteins = LR_genes)

LR_genes[!(LR_genes %in% human_LR$mouse_symbol)]

human_LR <- get_orthologs(LR_genes, input_species = "mouse")

annotations <- import_omnipath_annotations(
  proteins = human_LR$human_symbol[1:500],
  resources = c("NetPath")
  #resources = c("SignaLink_pathway")
  #resources = c("GO_Intercell", "SignaLink_pathway")
)

annotations2 <- import_omnipath_annotations(
  proteins = human_LR$human_symbol[501:1000],
  resources = c("NetPath", "SignaLink_pathway")
  #resources = c("SignaLink_pathway")
  #resources = c("GO_Intercell", "SignaLink_pathway")
)

annotations3 <- import_omnipath_annotations(
  proteins = human_LR$human_symbol[1001:1500],
  resources = c("NetPath", "SignaLink_pathway")
  #resources = c("SignaLink_pathway")
  #resources = c("GO_Intercell", "SignaLink_pathway")
)


