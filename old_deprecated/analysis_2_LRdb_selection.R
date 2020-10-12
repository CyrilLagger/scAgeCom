####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - July 2020
##
## Compare the Ligand-Receptor pairs from the various
## databases. Select the most relevant ones for
## the aging analysis.
##
####################################################
##

## Librarires ####

library(scDiffCom)
library(data.table)
library(VennDiagram)
library(ggplot2)
library(ggpmisc)
library(gridExtra)
library(grid)
library(gtable)
library(clusterProfiler)
library(org.Mm.eg.db)

## Data path ####
dir_data <- "../data_scAgeCom/"

## Load data ####
data("LRall")

## Summary data.frame ####

LRdb_summary <- data.frame(
  name = c("SingleCellSignalR", "CellPhoneDB", "NicheNet", "scTensor"),
  original_species = c("Homo Sapiens", "Homo Sapiens", "Homo Sapiens", "Mus Musculus"),
  curated = c("Yes", "Yes", "Partially", "No"),
  type = c("pair","subunit" , "pair", "pair"),
  number_orthologs = c(
    nrow(LRall[scsr == TRUE,]),
    nrow(LRall[cpdb == TRUE,]),
    nrow(LRall[nichenet == TRUE,]),
    nrow(LRall[sctensor == TRUE,])
  ),
  source = c(
    paste(unique(unlist(strsplit(unique(LRall[scsr == TRUE]$source_scsr), ","))), collapse = ", "),
    paste(c("IUPHAR", "DLRP", "iMEX", "Intact", "InnateDB", "DIP", "MINT", "I2D", "MatrixDB", "HPIDB"), collapse = ", "),
    paste(c("Ramilowski", "Kegg", "IUPHAR", "ppi_prediction"), collapse = ", "),
    paste(c("STRING", "TrEMBL", "SWISSPROT"), collapse = ", ")
  )
)

## Venn diagram ####

venn.diagram(
  x = list(
    LRall[scsr == TRUE,]$SYMB_LR,
    LRall[cpdb == TRUE,]$SYMB_LR,
    LRall[nichenet == TRUE,]$SYMB_LR,
    LRall[sctensor == TRUE,]$SYMB_LR
  ),
  category.names = c( "SingleCellSignalR", "CellPhoneDB", "NicheNet", "scTensor"),
  filename = paste0(dir_data, "analysis/analysis_2_plot_LR_venn.png"),
  imagetype = "png",
  col = "black",
  fill =  c("#6b7fff",  "#de4dff", "#ff4059", "#2cff21"),
  cat.col = c("#6b7fff",  "#de4dff", "#ff4059", "#2cff21"),
  cat.cex = 1,
  cat.pos = c(340, 10, 0, 0),
  margin = 0.15
)

## Look at single LR ####
sort(table(LRall$source_scsr), decreasing = TRUE)
sort(table(LRall[sctensor == FALSE &
                        nichenet == FALSE &
                        cpdb == FALSE & 
                        scsr == TRUE,
                      ]$source_scsr), decreasing = TRUE)

sort(table(LRall$source_nichenet), decreasing = TRUE)
sort(table(LRall[sctensor == FALSE &
                        nichenet == TRUE &
                        cpdb == FALSE & 
                        scsr == FALSE,
                      ]$source_nichenet), decreasing = TRUE)

sort(table(LRall$source_sctensor), decreasing = TRUE)
sort(table(LRall[sctensor == TRUE &
                        nichenet == FALSE &
                        cpdb == FALSE & 
                        scsr == FALSE,
                      ]$source_sctensor), decreasing = TRUE)

# look at LR pairs with L==R
LRall[GENESYMB_L == GENESYMB_R,]

## Nice data.frame for presentation ####

#nice data.frame visualization for the paper
df_tr <- LRdb_summary
df_tr<- transpose(df_tr, keep.names = "")
colnames(df_tr) <- df_tr[1,]
df_tr <- df_tr[-c(1),]
df_tr_rownames <- c(
  "Original species",
  "Curated",
  "Interaction type",
  "Interaction number",
  "Source"
)
df_tr$name <- df_tr_rownames
colnames(df_tr)[[1]] <- ""

for(i in 2:5) {
  df_tr[5,i] <- sapply(
    strwrap(
      df_tr[5,i],
      width = 25,
      simplify = FALSE),
    paste,
    collapse = "\n"
  )
}

g <- tableGrob(df_tr, rows = NULL)
grid.newpage()
grid.draw(g)
ggsave(paste0(dir_data, "analysis/analysis_2_plot_LRdb_summary_table.png"),
       plot = g, scale = 1.5)

# example data.frame with a few LR pairs
LRSYMB_forTalk <- c(
  "Gcg_Adora2a",  "Apoe_Sdc4",  "Vegfd_Itga5",  "Thbs1_Itga4",  "Itgb1_Plaur",
  "Hsp90b1_Egfr", "Gip_Pth1r",    "Ccl28_Ackr2",  "Avp_Fshr",     "Apln_Bdkrb1"
)
LR_forTalk <- LRall[SYMB_LR %in% LRSYMB_forTalk,]
g <- tableGrob(LR_forTalk, rows = NULL)
grid.newpage()
grid.draw(g)
ggsave(paste0(dir_data, "analysis/analysis_2_plot_LRdb_example_table.png"),
       plot = g, scale = 2)

## Select relevant LR pairs ####

LR_toKeep <- LRall[scsr == TRUE | cpdb == TRUE, ]
nrow(LR_toKeep)

# look at orthology
table(LR_toKeep$CONF_L, LR_toKeep$CONF_R)
table(LR_toKeep$TYPE_L, LR_toKeep$TYPE_R)

## Enrichment analysis on LR-pairs ####

# look at the genes in more details
Ligand_toKeep <- unique(LR_toKeep$GENESYMB_L)
Receptor_toKeep <- unique(LR_toKeep$GENESYMB_R)
common_toKeep <- intersect(Ligand_toKeep, Receptor_toKeep)
all_toKeep <- unique(c(Ligand_toKeep, Receptor_toKeep))

# list of GO analysis
ego_toEnrich <- list(
  all_MF = list(gene = all_toKeep, ont = "MF"),
  all_BP = list(gene = all_toKeep, ont = "BP"),
  all_CC = list(gene = all_toKeep, ont = "CC"),
  L_MF = list(gene = Ligand_toKeep, ont = "MF"),
  L_BP = list(gene = Ligand_toKeep, ont = "BP"),
  L_CC = list(gene = Ligand_toKeep, ont = "CC"),
  R_MF = list(gene = Receptor_toKeep, ont = "MF"),
  R_BP = list(gene = Receptor_toKeep, ont = "BP"),
  R_CC = list(gene = Receptor_toKeep, ont = "CC")
)

# apply ClusterProfiler
ego_results <- lapply(
  ego_toEnrich,
  function(x) {
    enrichGO(
      gene = x[["gene"]],
      OrgDb = org.Mm.eg.db,
      keyType = 'SYMBOL',
      ont = x[["ont"]],
      pAdjustMethod = "BH",
      pvalueCutoff = 0.01,
      qvalueCutoff  = 0.05
    )
  }
)

ego_plots <- lapply(
  ego_results,
  function(x) {
    dotplot(
      x,
      showCategory = 8,
      font.size = 10
    )
  }
)

g <- cowplot::plot_grid(
  plotlist = ego_plots[c(1,3,4,6,7,9)],
  ncol = 2,
  labels = c("LR - Molecular Function", "LR - Cellular Component",
             "Ligand - Molecular Function", "Ligand - Cellular Component",
             "Receptor - Molecular Function", "Receptor - Cellular Component"),
  align = "v"
)
ggsave(paste0(dir_data, "analysis/analysis_2_plot_enrichGO_LRKeep.png"),
       plot = g, scale = 1.8)




