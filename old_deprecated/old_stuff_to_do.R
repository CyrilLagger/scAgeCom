
go <- ontoProc::getGeneOnto()

go_id <- sample(go$id, 100)


sapply(
  go$id,
  function(i) {
    sapply(
      go$id, 
      function(j) {
        
      }
    )
  }
)


test <- findCommonAncestors()

co <- getCellOnto(useNew=TRUE)

# TODO: wrap in utility function.
parents <- co$parents
self <- rep(names(parents), lengths(parents))
library(igraph)
g <- make_graph(rbind(unlist(parents), self))
g

requireNamespace("Rgraphviz")
requireNamespace("graph")
cl = getCellOnto()
cl3k = c("CL:0000492", "CL:0001054", "CL:0000236", "CL:0000625",
         "CL:0000576", "CL:0000623", "CL:0000451", "CL:0000556")
p3k = ontologyPlot::onto_plot(cl, cl3k)
gnel = make_graphNEL_from_ontology_plot(p3k)
gnel = improveNodes(gnel, cl)
graph::graph.par(list(nodes=list(shape="plaintext", cex=.8)))
gnel = Rgraphviz::layoutGraph(gnel)
Rgraphviz::renderGraph(gnel)



ora_tables <- scDiffCom::get_ora_tables(obj)
names(ora_tables) <- c("GO Terms", "LR_CELLTYPE", "Cell Types", "Cell Families")
ora_tables <- lapply(
  ora_tables,
  function(ORA_dt) {
    dt <- copy(ORA_dt)
    setnames(
      dt,
      new = c("OR_UP", "pval_adjusted_UP", "OR_DOWN", "pval_adjusted_DOWN",
              "OR_FLAT", "pval_adjusted_FLAT"),
      old = c("Odds Ratio Up", "Adj. P-Value Up", "Odds Ratio Down", "Adj. P-Value Down",
              "Odds Ratio Stable", "Adj. P-Value Stable")
    )
    return(dt)
  }
)
obj <- scDiffCom::set_ora_tables(obj, ora_tables)
cci_table_filtered <- copy(scDiffCom:::get_cci_table_filtered(obj))
setnames(
  cci_table_filtered,
  new = c("L_CELLTYPE", "R_CELLTYPE", "BH_PVAL_DIFF", "LR_NAME"),
  old = c("Emitter Cell Type", "Receiver Cell Type", "Adj. P-Value", "Ligand-Receptor Genes")
)
obj <- scDiffCom:::set_cci_table_filtered(obj, cci_table_filtered)

test <- LR6db$LR6db_GO$LR_GO_intersection

#####

library(ontoProc)
onto_go <- getGeneOnto()
onto_go_names <- onto_go$name

go_id_name_dt <- data.table(
  name = onto_go_names,
  ID = names(onto_go_names)
)

test_ora <- merge.data.table(
  test_ora,
  go_id_name_dt,
  by.x = "Value",
  by.y = "ID",
  all.x = TRUE
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

###########



#####

heat_test <- DATASETS_light$facs_mixed$Lung@cci_table_filtered[, c("LR_CELLTYPE", "LR_NAME", "LOGFC")]

heat_test <- dcast.data.table(
  heat_test,
  formula = LR_NAME ~ LR_CELLTYPE,
  value.var = "LOGFC"
)
heat_test[is.na(heat_test)] <- 0

heat_test <- as.matrix(heat_test, rownames = "LR_NAME")

pheatmap::pheatmap(
  heat_test,
  show_rownames = FALSE,
  show_colnames = FALSE,
  scale = "row"
)

pheatmap::pheatmap(
  heat_test,
  show_rownames = FALSE,
  show_colnames = FALSE,
  scale = "row",
  kmeans_k = 10
)

testxx4 <- testxx3
testxx4[testxx4 <= log(1.2)] <- 0

testxx4 <- testxx4[rowSums(testxx4) > 0,]

pheatmap(testxx4, scale = "row")

pheatmap(testxx3, scale = "column")
pheatmap(testxx3, scale = "row")


########################
test_ora_up <- test_ora[, c("Value", "Value_NAME", "pval_adjusted_UP", "OR_UP", "ratio_UP", "my_OR_UP")]

test_ora_down <- test_ora[, c("Value", "Value_NAME", "pval_adjusted_DOWN", "OR_DOWN")]

ggplot(test_ora_up[OR_UP > 1 & pval_adjusted_UP <= 0.05], aes(x = log2(my_OR_UP), y = -log10(pval_adjusted_UP))) +
  geom_point() +
  geom_text(aes(label=ifelse(my_OR_UP >= 5 | pval_adjusted_UP <= 1.E-10, Value_NAME, "")), position=position_jitter(width=0.5,height=0.1))

ggplot(test_ora_up[OR_UP > 1 & pval_adjusted_UP <= 0.05], aes(x = log10(my_OR_UP), y = -log10(pval_adjusted_UP))) +
  geom_point() +
  geom_text(aes(label=ifelse(my_OR_UP >= 5 | pval_adjusted_UP <= 1.E-10, Value_NAME, "")))

test_ora_up[, new_col := -log10(pval_adjusted_UP)*log10(my_OR_UP)]

hist(test_ora_up$new_col, breaks = 50)

ggplot(test_ora_up[OR_UP > 1 & pval_adjusted_UP <= 0.05][order(-new_col)][1:20], aes(Value_NAME, new_col)) +
  geom_bar(stat = "identity") +
  coord_flip()




getLeavesFromTerm("GO:0002682", onto_go)


ggplot(test_ora_down[OR_DOWN > 1], aes(x = log2(OR_DOWN), y = -log10(pval_adjusted_DOWN))) + geom_point()

ggplot(test_ora_up[OR_UP > 1 & pval_adjusted_UP <= 0.05], aes(x = ratio_UP, y = -log10(pval_adjusted_UP))) + geom_point()

ggplot(test_ora_up[OR_UP > 1 & pval_adjusted_UP <= 0.05], aes(x = 1/ratio_UP, y = OR_UP)) + geom_point()
ggplot(test_ora_up[OR_UP > 1 & pval_adjusted_UP <= 0.05], aes(x = my_OR_UP, y = OR_UP)) + geom_point()

test2 <- test_ora_up[OR_UP > 1 & pval_adjusted_UP <= 0.05][order(-OR_UP)]

ggplot(test_ora_up[OR_UP > 1 & pval_adjusted_UP <= 0.05][order(-OR_UP)][1:100], aes(Value_NAME, OR_UP)) +
  geom_bar(stat = "identity") +
  coord_flip()




pheatmap(testxx3, scale = "row", kmeans_k = 50)

heatmap(testxx3)

vec_logfc <- seq(1.1, 2, by = 0.1)
names(vec_logfc) <- vec_logfc
vec_logfc <- log(vec_logfc)

bladder_test <- rbindlist(
  lapply(
    vec_logfc,
    function(logfc) {
      temp <- run_ORA(
        object = test,
        categories = "GO",
        logfc_threshold = logfc,
        overwrite = TRUE
      )
      temp_ORA <- scDiffCom:::get_ora_tables(temp)$GO
      temp_ORA <- temp_ORA[, c("Value", "Value_NAME", "pval_adjusted_UP", "OR_UP")]
    }
  ),
  use.names = TRUE,
  idcol = "logfc_thr"
)

ggplot(bladder_test[Value == "GO:0000165"], aes(x= logfc_thr, y = OR_UP)) + geom_point()
ggplot(bladder_test[Value == "GO:0000165"], aes(x= logfc_thr, y = -log10(pval_adjusted_UP))) + geom_point()

ggplot(bladder_test[Value == "GO:0019221"], aes(x= logfc_thr, y = OR_UP)) + geom_point()
ggplot(bladder_test[Value == "GO:0019221"], aes(x= logfc_thr, y = -log10(pval_adjusted_UP))) + geom_point()

ggplot(bladder_test[Value == "GO:1901143"], aes(x= logfc_thr, y = OR_UP)) + geom_point()
ggplot(bladder_test[Value == "GO:1901143"], aes(x= logfc_thr, y = -log10(pval_adjusted_UP))) + geom_point()



testx <- test@cci_table_filtered[, c("avg", "sd") := list(mean(LOGFC), sd(LOGFC)), by = "LR_NAME"][, c("LR_NAME", "avg", "sd")]
testx <- unique(testx)

ggplot(testx, aes(x= avg, y = sd)) + geom_point()


hist(test@cci_table_filtered[LOGFC > log(1.2),mean(LOGFC), by = "LR_NAME"][order(-V1)]$V1, breaks = 50)

test_ora_LR_name <- test@ORA$LR_NAME

test2 <- test_ora_LR_name[, c("Value", "OR_UP")][test@cci_table_filtered[LOGFC > log(1.5),mean(LOGFC), by = "LR_NAME"], on = "Value==LR_NAME", avg := i.V1]

ggplot(test2, aes(x=log10(OR_UP), y = avg)) + geom_point()

test_ora <- test@ORA$GO



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



ggplot(test_ora_goInter_up, aes(x = log10(OR_UP), y = -log10(pval_adjusted_UP), label = Value_NAME )) + geom_point() +
  geom_text(aes(label=ifelse(OR_UP >= 2, Value_NAME, "")), position=position_jitter(width=0.1,height=0.1))

ggplot(test_ora_goInter_down, aes(x = log10(OR_DOWN), y = -log10(pval_adjusted_DOWN), label = Value_NAME )) + geom_point() +
  geom_text(aes(label=ifelse(OR_DOWN >= 1, Value_NAME, "")), position=position_jitter(width=0.1,height=0.11))


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







## Messy file with some ideas
## will eventually be removed#####

SCAT_facs <- diffcom_results_detected$tms_facs[TISSUE == "SCAT",]

CCI_chord(SCAT_facs, "flat", TRUE, 10)
CCI_chord(SCAT_facs, "Up", TRUE, 5)
CCI_chord(SCAT_facs, "Down", TRUE, 5)

CCI_chord(SCAT_facs, "flat", FALSE, 10)
CCI_chord(SCAT_facs, "Up", FALSE, 5)
CCI_chord(SCAT_facs, "Down", FALSE, 5)

CCI_chord(SCAT_facs, "Down", TRUE, 10)
CCI_chord(SCAT_facs, "Down", FALSE, 10, "../data_scAgeCom/test.png")

gtest <- cowplot::as_grob(CCI_chord(SCAT_facs, "Down", TRUE, 10))

CCI_chord <- function(
  dt,
  direction,
  is_LR,
  filter,
  dir = NULL
) {
  dt_clean <- dt[SIG_TYPE != "FFF", ]
  if(direction == "Up") {
    dt_clean <- dt_clean[SIG_TYPE %in% c("FTT", "TTTU"),  ]
  } else if(direction == "Down") {
    dt_clean <- dt_clean[SIG_TYPE %in% c("TFT", "TTTD"),  ]
  }
  if(is_LR) {
    dt_clean <- merge.data.table(
      x = unique(dt_clean[, c("L_GENE", "R_GENE", "LR_GENES")]),
      y = dt_clean[,.N, by = "LR_GENES"][N >= filter, ],
      by = "LR_GENES",
      all = FALSE
    )[, LR_GENES := NULL]
    dt_clean[, L_GENE := paste0("L-", L_GENE)]
    dt_clean[, R_GENE := paste0("R-", R_GENE)]
  } else{
    dt_clean <- merge.data.table(
    x = unique(dt_clean[, c("L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPES")]),
    y = dt_clean[,.N, by = "LR_CELLTYPES"][N >= filter, ],
    by = "LR_CELLTYPES",
    all = FALSE
  )[, LR_CELLTYPES := NULL]
  dt_clean[, L_CELLTYPE := paste0("L-", L_CELLTYPE)]
  dt_clean[, R_CELLTYPE := paste0("R-", R_CELLTYPE)]
  }
  if(!is.null(dir)) {
    png(dir, width = 2000, height = 2000, res = 200)
  }
  circos.clear()
  chordDiagram(dt_clean, annotationTrack = "grid", 
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(dt_clean))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1.2)
  }, bg.border = NA)

  if(!is.null(dir)) {
    dev.off()
  }
}


B2M <- LRone2one[GENESYMB_L == "B2m",]


#setDT(diffcom_results_detected_strong$tms_facs)
B2M_facs <- diffcom_results_detected_strong$tms_facs[L_GENE == "B2m", ]
B2M_droplet <- diffcom_results_detected_strong$tms_droplet[L_GENE == "B2m", ]
table(B2M_facs$SIG_TYPE)

B2M_facs_up <- B2M_facs[SIG_TYPE %in% c("FTT", "TTTU")]
B2M_facs_up[, TISSUE_LR_CELLTYPES := paste(TISSUE, LR_CELLTYPES, sep = "_")]
B2M_droplet_up <- B2M_droplet[SIG_TYPE %in% c("FTT", "TTTU")]


sort(table(B2M_facs_up$TISSUE), decreasing = TRUE)
sort(table(B2M_droplet_up$TISSUE), decreasing = TRUE)
head(sort(table(B2M_facs_up$TISSUE_LR_CELLTYPES), decreasing = TRUE),10)
sort(table(B2M_facs_up$L_CELLTYPE), decreasing = TRUE)
sort(table(B2M_facs_up$R_GENE), decreasing = TRUE)
ggplot(B2M_facs_up, aes(x= LR_SCORE_old)) + geom_histogram(bins = 50)
ggplot(B2M_facs_up, aes(x= LR_SCORE_young)) + geom_histogram(bins = 50)
ggplot(B2M_facs_up, aes(x= LR_LOGFC)) + geom_histogram(bins = 50)
    


App_nc <- diffcom_results_significant$tms_facs[LR_GENES == "App_Ncstn",]
sort(table(diffcom_results_significant$tms_facs[L_GENE == "App",]$R_GENE))


diffcom_results_detected_strong$tms_facs[L_GENE == "Fgf2",]
test <- diffcom_results_significant$tms_facs[R_GENE == "Fgfr2",]

Itgb1 <- diffcom_results_significant$tms_facs[R_GENE == "Itgb1",]
table(Itgb1$SIG_TYPE)
sort(table(Itgb1$TISSUE), decreasing = TRUE)
sort(table(Itgb1$L_GENE), decreasing = TRUE)
sort(table(Itgb1$L_CELLTYPE), decreasing = TRUE)

Itgb1_drop <- diffcom_results_significant$tms_droplet[R_GENE == "Itgb1",]
table(Itgb1_drop$SIG_TYPE)
sort(table(Itgb1_drop$TISSUE), decreasing = TRUE)
sort(table(Itgb1$L_GENE), decreasing = TRUE)
sort(table(Itgb1$L_CELLTYPE), decreasing = TRUE)


Spp1_facs <- diffcom_results_significant$tms_facs[L_GENE == "Spp1",]
Spp1_drop <- diffcom_results_significant$tms_droplet[L_GENE == "Spp1",]

Itgav_drop <- diffcom_results_significant$tms_droplet[R_GENE == "Itgav",]
table(Itgav_drop$SIG_TYPE)

"Il6" %in% diffcom_results$tms_facs$L_GENE

"Il6" %in% LRone2one$GENESYMB_L

Il6 <- LRone2one[GENESYMB_L == "Il6",]
LRone2one[GENESYMB_R == "Il6ra",]

library(SingleCellSignalR)
data("LRdb")
LRdb
"IL6" %in% LRdb$ligand

LRdb[LRdb$ligand == "IL6",]

data()

IL6sc <- LRdb[LRdb$ligand ==  "IL6",]

LRall$LRall_many2many[LRall$LRall_many2many$GENESYMB_R == "Il6ra", ]


mm2Hs[mm2Hs$`Gene name` == "IL6R",]

test <- diffcom_results_significant$tms_facs
test2 <- diffcom_results_significant$tms_droplet



##volcano plot #####
library(ggrepel)
diffcom_1000iter <- readRDS("test_and_comparison/data_results_diffcom_10000iter_log.rds")
diffcom_1000iter[, LR_DIFF := LR_SCORE_old - LR_SCORE_young]
ggplot(diffcom_1000iter[(LR_DETECTED_young | LR_DETECTED_old), ], 
       aes(x = LR_DIFF, y = -log10(BH_PVAL_DIFF + 1E-4))) +
  geom_point() +
  geom_vline(xintercept = log(1.1)) +
  geom_vline(xintercept = -log(1.1)) +
  geom_hline(yintercept = -log10(0.05))

+
  geom_text(
    data = diffcom_1000iter[(LR_DETECTED_young | LR_DETECTED_old) & BH_PVAL_DIFF < 0.01 & abs(LR_DIFF) > 0.5, ],
    aes(label = LR_GENES),
    size = 5
    #box.padding = unit(0.35, "lines"),
    #point.padding = unit(0.3, "lines")
  )

hist(log10(diffcom_1000iter$LR_SCORE_old), breaks = 100)
hist(log10(diffcom_1000iter$LR_SCORE_young), breaks = 100)

ggplot(diffcom_1000iter[(LR_DETECTED_young | LR_DETECTED_old) & (BH_PVAL_young <= 0.05 | BH_PVAL_old <= 0.05) & 
                          (LR_SCORE_young > 0.1 | LR_SCORE_old > 0.1), ], 
       aes(x = LR_DIFF, y = -log10(BH_PVAL_DIFF + 1E-4))) +
  geom_point() +
  geom_vline(xintercept = log(1.1)) +
  geom_vline(xintercept = -log(1.1)) +
  geom_hline(yintercept = -log10(0.05))


+
  geom_text(
    data = diffcom_1000iter[(LR_DETECTED_young | LR_DETECTED_old) & BH_PVAL_DIFF < 0.01 & abs(LR_DIFF) > 0.3, ],
    aes(label = L_CELLTYPE),
    size = 5
    #box.padding = unit(0.35, "lines"),
    #point.padding = unit(0.3, "lines")
  )






calico_previous <- bind_tissues(data_path$calico, tissue_list$calico)

calico_new <- bind_tissues("../../../../../scDiffCom_results/diffcom_calico_size_factor_log_10000iter_mixed", tissue_list$calico)

table(calico_new$LR_DETECTED_old, calico_new$LR_DETECTED_young)
table(calico_previous$LR_DETECTED_old, calico_previous$LR_DETECTED_young)


diffcom_results <- mapply(bind_tissues, data_path, tissue_list, SIMPLIFY = FALSE)

facs_previous <- bind_tissues(data_path$tms_facs, tissue_list$tms_facs)
facs_news <- bind_tissues("../../../../../scDiffCom_results/diffcom_tms_facs_size_factor_log_10000iter_mixed", tissue_list$tms_facs)

facs_previous[, id_col := paste(LR_GENES, L_CELLTYPE, R_CELLTYPE, TISSUE, sep = "_")]
facs_news[, id_col := paste(LR_GENES, L_CELLTYPE, R_CELLTYPE, TISSUE, sep = "_")]


table(facs_news$LR_DETECTED_old, facs_news$LR_DETECTED_young)
table(facs_previous$LR_DETECTED_old, facs_previous$LR_DETECTED_young)


facs_comp <- facs_news[id_col %in% facs_previous$id_col,]
facs_previous <- unique(facs_previous)


table(facs_comp$LR_DETECTED_old, facs_comp$LR_DETECTED_young)
table(facs_previous$LR_DETECTED_old, facs_previous$LR_DETECTED_young)



library(config)
library(conflicted)
library(circlize)
library(ggplot2)

library(clusterProfiler)
library(org.Mm.eg.db)



#check_representation_inv(diffcom_results$calico)

#Cutoff exploration
#df_detec_in_young_and_old <- lapply(
#  diffcom_results,
#  remove_undetected_interactions
#)
#calico
#explore_cutoffs_dir_calico <- "test_and_comparison/explore_cutoffs/calico/"
#create_dir(explore_cutoffs_dir_calico)
#explore_calico <- explore_filter_cutoffs(df_detec_in_young_and_old$calico, explore_cutoffs_dir_calico) 
#tms FACS
#explore_cutoffs_dir_tms_facs <- "test_and_comparison/explore_cutoffs/tms_facs/"
#create_dir(explore_cutoffs_dir_tms_facs)
#explore_tms_facs <- explore_filter_cutoffs(df_detec_in_young_and_old$tms_facs, explore_cutoffs_dir_tms_facs) 
#tms Droplet
#explore_cutoffs_dir_tms_droplet <- "test_and_comparison/explore_cutoffs/tms_droplet/"
#create_dir(explore_cutoffs_dir_tms_droplet)
#explore_tms_droplet <- explore_filter_cutoffs(df_detec_in_young_and_old$tms_droplet, explore_cutoffs_dir_tms_droplet)

#explore_calico
#explore_tms_facs
#explore_tms_droplet



#####


######
#Some sub-datatables
diffcom_results_detected_strong <- lapply(
  diffcom_results_detected,
  function(x) {
    dt <- x[!(SIG_TYPE == "FFF"),]
  }
)

diffcom_results_significant <- lapply(
  diffcom_results_detected_strong,
  function(x){
    dt <- x[SIG_TYPE %in% c("FTT", "TFT", "TTTD", "TTTU"),]
  }
)






##
SCAT_facs <- diffcom_results_detected$tms_facs[TISSUE == "SCAT",]

CCI_chord(SCAT_facs, "flat", TRUE, 10)
CCI_chord(SCAT_facs, "Up", TRUE, 5)
CCI_chord(SCAT_facs, "Down", TRUE, 5)

CCI_chord(SCAT_facs, "flat", FALSE, 10)
CCI_chord(SCAT_facs, "Up", FALSE, 5)
CCI_chord(SCAT_facs, "Down", FALSE, 5)

CCI_chord(SCAT_facs, "Down", TRUE, 10)
CCI_chord(SCAT_facs, "Down", FALSE, 10, "../data_scAgeCom/test.png")

gtest <- cowplot::as_grob(CCI_chord(SCAT_facs, "Down", TRUE, 10))

CCI_chord <- function(
  dt,
  direction,
  is_LR,
  filter,
  dir = NULL
) {
  dt_clean <- dt[SIG_TYPE != "FFF", ]
  if(direction == "Up") {
    dt_clean <- dt_clean[SIG_TYPE %in% c("FTT", "TTTU"),  ]
  } else if(direction == "Down") {
    dt_clean <- dt_clean[SIG_TYPE %in% c("TFT", "TTTD"),  ]
  }
  if(is_LR) {
    dt_clean <- merge.data.table(
      x = unique(dt_clean[, c("L_GENE", "R_GENE", "LR_GENES")]),
      y = dt_clean[,.N, by = "LR_GENES"][N >= filter, ],
      by = "LR_GENES",
      all = FALSE
    )[, LR_GENES := NULL]
    dt_clean[, L_GENE := paste0("L-", L_GENE)]
    dt_clean[, R_GENE := paste0("R-", R_GENE)]
  } else{
    dt_clean <- merge.data.table(
      x = unique(dt_clean[, c("L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPES")]),
      y = dt_clean[,.N, by = "LR_CELLTYPES"][N >= filter, ],
      by = "LR_CELLTYPES",
      all = FALSE
    )[, LR_CELLTYPES := NULL]
    dt_clean[, L_CELLTYPE := paste0("L-", L_CELLTYPE)]
    dt_clean[, R_CELLTYPE := paste0("R-", R_CELLTYPE)]
  }
  if(!is.null(dir)) {
    png(dir, width = 2000, height = 2000, res = 200)
  }
  circos.clear()
  chordDiagram(dt_clean, annotationTrack = "grid", 
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(dt_clean))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1.2)
  }, bg.border = NA)
  
  if(!is.null(dir)) {
    dev.off()
  }
}


###############
#############
#########

test_lung_detected <- test_lung[SIG_TYPE %in% c("FTT", "TFT", "TTF", "TTTD", "TTTU"), ]
test_lung_sig <- test_lung[SIG_TYPE %in% c("FTT", "TFT", "TTTD", "TTTU"),]
colnames(test_lung_detected)

df_d <- setDF(test_lung_detected)
df_sig <- setDF(test_lung_sig)


ora_test <- analyze_overrepresentation(df_d, df_sig,
                                       exclude_tissue_overrepresentation=TRUE,
                                       adjust_pvals=FALSE)

orat_test_LR_GENES <- ora_test$L_GENE


#########

as.data.frame.matrix(table(diffcom_results_detected$tms_facs$SIG_TYPE))

#detected interaction by LR_CELLTYPES
table(diffcom_results_detected$tms_facs[TISSUE == "Liver",]$SIG_TYPE)

diffcom_results_detected$tms_facs[, LR_T_CELLTYPES := paste(TISSUE, L_CELLTYPE, TISSUE, R_CELLTYPE, sep = "_" ) ]
distr_CCI_tms_facs <- dcast(diffcom_results_detected$tms_facs[!(SIG_TYPE %in% "FFF"), c("LR_T_CELLTYPES", "SIG_TYPE", "TISSUE") ],
                            TISSUE + LR_T_CELLTYPES ~ SIG_TYPE, value.var = "SIG_TYPE")
distr_CCI_tms_facs[, N := rowSums(.SD), .SDcols = colnames(distr_CCI_tms_facs)[-c(1,2)]]
distr_CCI_tms_facs[, young := rowSums(.SD), .SDcols = c("TFT", "TTF", "TTTD", "TTTU")]
distr_CCI_tms_facs[, old := rowSums(.SD), .SDcols = c("FTT", "TTF", "TTTD", "TTTU")]
distr_CCI_tms_facs[, up := rowSums(.SD), .SDcols = c("FTT", "TTTU")]
distr_CCI_tms_facs[, down := rowSums(.SD), .SDcols = c("TFT", "TTTD")]
distr_CCI_tms_facs[, appear := old - young]
distr_CCI_tms_facs[, net_reg := up - down]
hist(distr_CCI_tms_facs$appear, breaks = 50)
hist(distr_CCI_tms_facs$net_reg, breaks = 50)
hist(distr_CCI_tms_facs[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$appear, breaks = 50)
hist(distr_CCI_tms_facs[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$net_reg, breaks = 50)
ggplot(distr_CCI_tms_facs, aes(x = appear, y = net_reg)) + geom_point()

het_tms_facs <- merge.data.table(
  diffcom_results_detected$tms_facs[SIG_TYPE %in% c("TFT", "TTF", "TTTD", "TTTU"),
                                    .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                    by = "LR_T_CELLTYPES" ],
  diffcom_results_detected$tms_facs[SIG_TYPE %in% c("FTT", "TTF", "TTTD", "TTTU"),
                                    .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                    by = "LR_T_CELLTYPES" ],
  by = "LR_T_CELLTYPES",
  suffixes = c("_young", "_old")
)
het_tms_facs[, c("diff_LR", "diff_L", "diff_R") := .(V1_old - V1_young, V2_old - V2_young, V3_old - V3_young)]
hist(het_tms_facs$diff_LR, breaks = 100)
hist(het_tms_facs$diff_L, breaks = 50)
hist(het_tms_facs$diff_R, breaks = 50)



diffcom_results_detected$tms_droplet[, LR_T_CELLTYPES := paste(TISSUE, L_CELLTYPE, TISSUE, R_CELLTYPE, sep = "_" ) ]
distr_CCI_tms_droplet <- dcast(diffcom_results_detected$tms_droplet[!(SIG_TYPE %in% "FFF"), c("LR_T_CELLTYPES", "SIG_TYPE", "TISSUE") ],
                               TISSUE + LR_T_CELLTYPES ~ SIG_TYPE, value.var = "SIG_TYPE")
distr_CCI_tms_droplet[, N := rowSums(.SD), .SDcols = colnames(distr_CCI_tms_droplet)[-c(1,2)]]
distr_CCI_tms_droplet[, young := rowSums(.SD), .SDcols = c("TFT", "TTF", "TTTD", "TTTU")]
distr_CCI_tms_droplet[, old := rowSums(.SD), .SDcols = c("FTT", "TTF", "TTTD", "TTTU")]
distr_CCI_tms_droplet[, up := rowSums(.SD), .SDcols = c("FTT", "TTTU")]
distr_CCI_tms_droplet[, down := rowSums(.SD), .SDcols = c("TFT", "TTTD")]
distr_CCI_tms_droplet[, appear := old - young]
distr_CCI_tms_droplet[, net_reg := up - down]
hist(distr_CCI_tms_droplet[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$appear, breaks = 50)
hist(distr_CCI_tms_droplet[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$net_reg, breaks = 50)
hist(distr_CCI_tms_droplet$appear, breaks = 50)
hist(distr_CCI_tms_droplet$net_reg, breaks = 50)
ggplot(distr_CCI_tms_droplet, aes(x = appear, y = net_reg)) + geom_point()

het_tms_droplet <- merge.data.table(
  diffcom_results_detected$tms_droplet[SIG_TYPE %in% c("TFT", "TTF", "TTTD", "TTTU"),
                                       .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                       by = "LR_T_CELLTYPES" ],
  diffcom_results_detected$tms_droplet[SIG_TYPE %in% c("FTT", "TTF", "TTTD", "TTTU"),
                                       .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                       by = "LR_T_CELLTYPES" ],
  by = "LR_T_CELLTYPES",
  suffixes = c("_young", "_old")
)
het_tms_droplet[, c("diff_LR", "diff_L", "diff_R") := .(V1_old - V1_young, V2_old - V2_young, V3_old - V3_young)]
hist(het_tms_droplet$diff_LR, breaks = 100)
hist(het_tms_droplet$diff_L, breaks = 50)
hist(het_tms_droplet$diff_R, breaks = 50)

diffcom_results_detected$calico[, LR_T_CELLTYPES := paste(TISSUE, L_CELLTYPE, TISSUE, R_CELLTYPE, sep = "_" ) ]
distr_CCI_calico <- dcast(diffcom_results_detected$calico[!(SIG_TYPE %in% "FFF"), c("LR_T_CELLTYPES", "SIG_TYPE", "TISSUE") ],
                          TISSUE + LR_T_CELLTYPES ~ SIG_TYPE, value.var = "SIG_TYPE")
distr_CCI_calico[, N := rowSums(.SD), .SDcols = colnames(distr_CCI_calico)[-c(1,2)]]
distr_CCI_calico[, young := rowSums(.SD), .SDcols = c("TFT", "TTF", "TTTD", "TTTU")]
distr_CCI_calico[, old := rowSums(.SD), .SDcols = c("FTT", "TTF", "TTTD", "TTTU")]
distr_CCI_calico[, up := rowSums(.SD), .SDcols = c("FTT", "TTTU")]
distr_CCI_calico[, down := rowSums(.SD), .SDcols = c("TFT", "TTTD")]
distr_CCI_calico[, appear := old - young]
distr_CCI_calico[, net_reg := up - down]
hist(distr_CCI_calico$appear, breaks = 50)
hist(distr_CCI_calico$net_reg, breaks = 50)
ggplot(distr_CCI_calico, aes(x = appear, y = net_reg)) + geom_point()

het_calico <- merge.data.table(
  diffcom_results_detected$calico[SIG_TYPE %in% c("TFT", "TTF", "TTTD", "TTTU"),
                                  .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                  by = "LR_T_CELLTYPES" ],
  diffcom_results_detected$calico[SIG_TYPE %in% c("FTT", "TTF", "TTTD", "TTTU"),
                                  .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                                  by = "LR_T_CELLTYPES" ],
  by = "LR_T_CELLTYPES",
  suffixes = c("_young", "_old")
)
het_calico[, c("diff_LR", "diff_L", "diff_R") := .(V1_old - V1_young, V2_old - V2_young, V3_old - V3_young)]
hist(het_calico$diff_LR, breaks = 100)
hist(het_calico$diff_L, breaks = 50)
hist(het_calico$diff_R, breaks = 50)


hist(distr_CCI_tms_facs$N, breaks = 50)
hist(distr_CCI_tms_droplet$N, breaks = 50)
hist(distr_CCI_calico$N, breaks = 30)

#maybe not all genes in calico!!!
hist(distr_CCI_tms_facs[TISSUE == "Kidney", ]$N, breaks = 50)
hist(distr_CCI_tms_droplet[TISSUE == "Kidney", ]$N, breaks = 50)
hist(distr_CCI_calico[TISSUE == "kidney", ]$N, breaks = 50)

#difference in detection
test_lung_sig[SIG_TYPE %in% c("TTTU", "TTTD", "TTF", "TFT"),]
test_lung_sig[SIG_TYPE %in% c("TTTU", "TTTD", "TTF", "FTT"),]


distr_LR_tms_facs <- dcast(diffcom_results_detected$tms_facs[!(SIG_TYPE %in% "FFF"), c("LR_GENES", "SIG_TYPE", "TISSUE") ],
                           LR_GENES ~ SIG_TYPE, value.var = "SIG_TYPE")

distr_LR_tiss_tms_facs <- merge.data.table(
  dcast(diffcom_results_detected$tms_facs[!(SIG_TYPE %in% "FFF"), c("LR_GENES", "SIG_TYPE", "TISSUE") ],
        TISSUE + LR_GENES ~ SIG_TYPE, value.var = "SIG_TYPE"),
  diffcom_results$tms_facs[ , length(unique(LR_CELLTYPES)), by = "TISSUE"],
  by = "TISSUE",
  all.x =  TRUE
)
cols <- c("FTT", "TFT", "TTF", "TTTD", "TTTU")
distr_LR_tiss_tms_facs[, N := rowSums(.SD) , .SDcols = cols]

distr_LR_tiss_tms_facs[, paste0("pct_", c(cols, "N")) := .SD/V1 , .SDcols = c(cols, "N")]

####################



ggplot(diffcom_results$tms_facs, aes(x = LR_SCORE_young)) + geom_histogram(bins = 100) + scale_y_log10()
ggplot(diffcom_results$tms_facs, aes(x = LR_SCORE_old)) + geom_histogram(bins = 150) + scale_y_log10()
ggplot(diffcom_results$tms_facs, aes(x = LR_LOGFC)) + geom_histogram(bins = 150) + scale_y_log10()

ggplot(diffcom_results$tms_facs[LR_DETECTED_young == TRUE,], aes(x = LR_SCORE_young)) + 
  geom_histogram(bins = 100) + scale_y_log10()
ggplot(diffcom_results$tms_facs[LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05,], aes(x = LR_SCORE_old)) + 
  geom_histogram(bins = 100) + scale_y_log10()
ggplot(diffcom_results$tms_facs[(LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05) |
                                  (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05),], aes(x = LR_LOGFC_ABS)) + 
  geom_histogram(bins = 100) + scale_y_log10()

ggplot(diffcom_results$tms_facs[(LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05 & LR_SCORE_old >= 0.85) |
                                  (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05 & LR_SCORE_young >= 0.85),],
       aes(x = LR_LOGFC_ABS)) + 
  geom_histogram(bins = 100) + scale_y_log10()

ggplot(diffcom_results$tms_facs[((LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05 & LR_SCORE_old >= 0.85) |
                                   (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05 & LR_SCORE_young >= 0.85)) &
                                  BH_PVAL_DIFF <= 0.05,],
       aes(x = LR_LOGFC_ABS)) + 
  geom_histogram(bins = 100) + scale_y_log10()

ggplot(diffcom_results$calico[((LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05 & LR_SCORE_old >= 1) |
                                 (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05 & LR_SCORE_young >= 1)) &
                                BH_PVAL_DIFF <= 0.05,],
       aes(x = LR_LOGFC)) + 
  geom_histogram(bins = 100) + scale_y_log10()

ggplot(diffcom_results$tms_droplet[((LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05 & LR_SCORE_old >= 1) |
                                      (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05 & LR_SCORE_young >= 1)) &
                                     BH_PVAL_DIFF <= 0.05,],
       aes(x = LR_LOGFC_ABS)) + 
  geom_histogram(bins = 100) + scale_y_log10()

ggplot(diffcom_results$tms_facs[LR_DETECTED_old == TRUE & BH_PVAL_old | LR_DETECTED_young == TRUE ,], aes(x = LR_LOGFC_ABS)) + 
  geom_histogram(bins = 100) + scale_y_log10()



quantile(diffcom_results$tms_facs[((LR_DETECTED_old == TRUE & BH_PVAL_old <= 0.05 & LR_SCORE_old >= 0.85) |
                                     (LR_DETECTED_young == TRUE & BH_PVAL_young <= 0.05 & LR_SCORE_young >= 0.85)) &
                                    BH_PVAL_DIFF <= 0.05,]$LR_LOGFC_ABS, probs = seq(0,1, 0.1))


########
##non immune cells
diffcom_results_detected$tms_facs[, c("L_T_CELLTYPE", "R_T_CELLTYPE") := .(paste(TISSUE, L_CELLTYPE, sep = "_"), paste(TISSUE, R_CELLTYPE, sep = "_"))]

write.csv(unique(diffcom_results_detected$tms_facs[, .(L_T_CELLTYPE)]),
          file = "../../../../../facs_t_cell_type.csv",
          row.names = FALSE
)

identical(
  sort(unique(diffcom_results_detected$tms_facs[, L_T_CELLTYPE])),
  sort(unique(diffcom_results_detected$tms_facs[, R_T_CELLTYPE]))
)

cell_types_immune_tms_facs <- read.csv("../../../../../facs_t_cell_type_sorted.csv",
                                       header = TRUE, 
                                       stringsAsFactors = FALSE)
table(cell_types_immune_tms_facs$Immune)
anyNA(cell_types_immune_tms_facs)

setDT(cell_types_immune_tms_facs)
ct_nonIm_tms_facs <- cell_types_immune_tms_facs[Immune == FALSE, L_T_CELLTYPE]

tms_facs_nonIm <- diffcom_results_detected$tms_facs[L_T_CELLTYPE %in% ct_nonIm_tms_facs & R_T_CELLTYPE %in% ct_nonIm_tms_facs, ]

distr_nonIm_CCI_tms_facs <- dcast(tms_facs_nonIm[!(SIG_TYPE %in% "FFF"), c("LR_T_CELLTYPES", "SIG_TYPE", "TISSUE") ],
                                  TISSUE + LR_T_CELLTYPES ~ SIG_TYPE, value.var = "SIG_TYPE")
distr_nonIm_CCI_tms_facs[, N := rowSums(.SD), .SDcols = colnames(distr_nonIm_CCI_tms_facs)[-c(1,2)]]
distr_nonIm_CCI_tms_facs[, young := rowSums(.SD), .SDcols = c("TFT", "TTF", "TTTD", "TTTU")]
distr_nonIm_CCI_tms_facs[, old := rowSums(.SD), .SDcols = c("FTT", "TTF", "TTTD", "TTTU")]
distr_nonIm_CCI_tms_facs[, up := rowSums(.SD), .SDcols = c("FTT", "TTTU")]
distr_nonIm_CCI_tms_facs[, down := rowSums(.SD), .SDcols = c("TFT", "TTTD")]
distr_nonIm_CCI_tms_facs[, appear := old - young]
distr_nonIm_CCI_tms_facs[, net_reg := up - down]
hist(distr_nonIm_CCI_tms_facs$appear, breaks = 50)
hist(distr_nonIm_CCI_tms_facs$net_reg, breaks = 50)
hist(distr_nonIm_CCI_tms_facs[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$appear, breaks = 50)
hist(distr_nonIm_CCI_tms_facs[TISSUE %in% c("Kidney", "Spleen", "Lung"),]$net_reg, breaks = 50)
ggplot(distr_nonIm_CCI_tms_facs, aes(x = appear, y = net_reg)) + geom_point()

het_nonIm_tms_facs <- merge.data.table(
  tms_facs_nonIm[SIG_TYPE %in% c("TFT", "TTF", "TTTD", "TTTU"),
                 .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                 by = "LR_T_CELLTYPES" ],
  tms_facs_nonIm[SIG_TYPE %in% c("FTT", "TTF", "TTTD", "TTTU"),
                 .(length(unique(LR_GENES)), length(unique(L_GENE)), length(unique(R_GENE)) ),
                 by = "LR_T_CELLTYPES" ],
  by = "LR_T_CELLTYPES",
  suffixes = c("_young", "_old")
)
het_nonIm_tms_facs[, c("diff_LR", "diff_L", "diff_R") := .(V1_old - V1_young, V2_old - V2_young, V3_old - V3_young)]
hist(het_nonIm_tms_facs$diff_LR, breaks = 100)
hist(het_nonIm_tms_facs$diff_L, breaks = 50)
hist(het_nonIm_tms_facs$diff_R, breaks = 50)

################

test_lung <- test_d[TISSUE == "Lung",]
test_lung <- diffcom_results_detected$tms_facs[TISSUE == "Liver",]
test_lung <- tms_facs_nonIm[TISSUE == "Lung",]

test_lung[SIG_TYPE %in% c("TFT", "FTT", "TTTU", "TTTD") & L_GENE == "Ccl5" & R_GENE == "Cxcr6",]
test_lung[SIG_TYPE %in% c("TFT", "FTT", "TTTU", "TTTD") & L_GENE == "Mrc1" & R_GENE == "Ptprc" & L_CELLTYPE == "hepatocyte",]

diffcom_results_detected$tms_facs[LR_GENES == "Ltb_Ltbr" & LR_CELLTYPES ==  "B cell_Kupffer cell"]
test_lung[LR_GENES == "Mrc1_Ptprc" & LR_CELLTYPES == "endothelial cell of hepatic sinusoid_B cell"]

test_lung <- test_lung[, c("L_GENE", "R_GENE", "L_CELLTYPE", "R_CELLTYPE",
                           "LR_KEEP_young", "LR_KEEP_old", "LR_LOGFC", "LR_CELLTYPES", "LR_GENES", "SIG_TYPE")]
test_lung[, L := paste("L", L_CELLTYPE, L_GENE, sep = "_")]
test_lung[, R := paste("R", R_CELLTYPE, R_GENE, sep = "_")]

test_lung_d <- test_lung[ SIG_TYPE != "FFF", ]
test_lung_up <- test_lung[SIG_TYPE %in% c("FTT", "TTTU"),  ]
test_lung_down <- test_lung[SIG_TYPE %in% c("TFT", "TTTD"),  ]

test_lung_d_up <- merge.data.table(test_lung_d[, .N, by = LR_CELLTYPES],
                                   test_lung_up[, .N, by = LR_CELLTYPES],
                                   by = "LR_CELLTYPES", 
                                   all.x = TRUE
)
test_lung_d_up[is.na(test_lung_d_up)] <- 0
test_lung_d_up[, ratio := N.y/N.x]

test_lung_d_down <- merge.data.table(test_lung_d[, .N, by = LR_CELLTYPES],
                                     test_lung_down[, .N, by = LR_CELLTYPES],
                                     by = "LR_CELLTYPES", 
                                     all.x = TRUE
)
test_lung_d_down[is.na(test_lung_d_down)] <- 0
test_lung_d_down[, ratio := N.y/N.x]



test_lung_chord <- data.table::merge.data.table(y = test_lung_d[, .N, by = LR_CELLTYPES],
                                                x = unique(test_lung_d[, c("L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPES")]),
                                                by = "LR_CELLTYPES",
                                                all = FALSE)[, LR_CELLTYPES := NULL]
test_lung_chord[, L_CELLTYPE := paste0("L_", L_CELLTYPE)]
test_lung_chord[, R_CELLTYPE := paste0("R_", R_CELLTYPE)]
chordDiagram(test_lung_chord, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(test_lung_chord))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()




test_lung_up_chord <- data.table::merge.data.table(y = test_lung_d_up[, c("LR_CELLTYPES", "ratio")],
                                                   x = unique(test_lung_up[, c("L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPES")]),
                                                   by = "LR_CELLTYPES",
                                                   all = FALSE)[, LR_CELLTYPES := NULL]
test_lung_up_chord[, L_CELLTYPE := paste0("L_", L_CELLTYPE)]
test_lung_up_chord[, R_CELLTYPE := paste0("R_", R_CELLTYPE)]
chordDiagram(test_lung_up_chord, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(test_lung_up_chord))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()

test_lung_down_chord <- data.table::merge.data.table(y = test_lung_d_down[, c("LR_CELLTYPES", "ratio")],
                                                     x = unique(test_lung_down[, c("L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPES")]),
                                                     by = "LR_CELLTYPES",
                                                     all = FALSE)[, LR_CELLTYPES := NULL]
test_lung_down_chord[, L_CELLTYPE := paste0("L_", L_CELLTYPE)]
test_lung_down_chord[, R_CELLTYPE := paste0("R_", R_CELLTYPE)]
chordDiagram(test_lung_down_chord, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(test_lung_down_chord))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()

test_lung_up_chord_g <- data.table::merge.data.table(y = test_lung_up[, .N, by = LR_GENES],
                                                     x = unique(test_lung_up[, c("L_GENE", "R_GENE", "LR_GENES")]),
                                                     by = "LR_GENES",
                                                     all = FALSE)[, LR_GENES := NULL]
test_lung_up_chord_g[, L_GENE := paste0("L_", L_GENE)]
test_lung_up_chord_g[, R_GENE := paste0("R_", R_GENE)]
chordDiagram(test_lung_up_chord_g, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(test_lung_up_chord_g))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()

test_lung_down_chord_g <- data.table::merge.data.table(y = test_lung_down[, .N, by = LR_GENES][N > 7, ],
                                                       x = unique(test_lung_down[, c("L_GENE", "R_GENE", "LR_GENES")]),
                                                       by = "LR_GENES",
                                                       all = FALSE)[, LR_GENES := NULL]
test_lung_down_chord_g[, L_GENE := paste0("L_", L_GENE)]
test_lung_down_chord_g[, R_GENE := paste0("R_", R_GENE)]
chordDiagram(test_lung_down_chord_g, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(test_lung_down_chord_g))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()






#try with circos
#test_lung_d_up <- test_lung_d[LR_UP == TRUE, c("L", "R") ]
#test_lung_d_down <- test_lung_d[LR_DOWN == TRUE, c("L", "R") ]

factor_total <- factor(c(sort(unique(test_lung_up$L)), sort(unique(test_lung_up$R))))
factor_LR <- factor(c(rep("L", length(unique(test_lung_up$L))), rep("R", length(unique(test_lung_up$R))) ))
factor_col_LR <- c(rep(rand_color(1), length(unique(test_lung_up$L))), rep(rand_color(1), length(unique(test_lung_up$R))) )

circos.par(
  cell.padding = c(0.00, 0, 0.00, 0),
  gap.degree = 0.,
  start.degree = 270,
  track.height = 0.01)
circos.initialize(factors = factor_total, xlim = c(0,1) )
#circos.track(ylim = c(0,0.1))
#circos.track(ylim = c(0,0.1))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ycenter,
              #CELL_META$sector.index,
              factor_LR[CELL_META$sector.numeric.index],
              col = factor_col_LR[CELL_META$sector.numeric.index]
              #facing = "inside", niceFacing = TRUE
  )
})
for(i in 1:nrow(test_lung_up)) {
  circos.link(test_lung_up$L[[i]], 0, test_lung_up$R[[i]], 0,
              col = "red",
              directional = 0)
}
for(i in 1:nrow(test_lung_down)) {
  circos.link(test_lung_down$L[[i]], 0, test_lung_down$R[[i]], 0,
              col = "blue",
              directional = 1)
}
circos.clear()





########enrichment####
setDT(test_lung_sig)
up_genes <- c(unique(test_lung_sig[SIG_TYPE %in% c("FTT", "TTTU"),]$L_GENE),
              unique(test_lung_sig[SIG_TYPE %in% c("FTT", "TTTU"),]$R_GENE))
down_genes <- c(unique(test_lung_sig[SIG_TYPE %in% c("TFT", "TTTD"),]$L_GENE),
                unique(test_lung_sig[SIG_TYPE %in% c("TFT", "TTTD"),]$R_GENE))
intersect(up_genes, down_genes)

library(clusterProfiler)
library(org.Mm.eg.db)




test_lung_diffnet <- test_lung[, lapply(.SD, sum), by = LR_CELLTYPES, .SDcols = c("LR_UP", "LR_DOWN")]

test_lung_net <- test_lung[, lapply(.SD, sum), by = LR_CELLTYPES, .SDcols = c("LR_KEEP_young", "LR_KEEP_old", "LR_KEEP_DIFF")]


test_lung[LR_CELLTYPES == "B cell_dendritic cell" & (LR_KEEP_old == TRUE | LR_KEEP_young == TRUE),]


table(test$LR_KEEP_DIFF & test$LR_LOGFC > 0)
table(test$LR_KEEP_DIFF & test$LR_LOGFC < 0)
table(test$LR_KEEP_young)
table(test$LR_KEEP_old)





sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC > 0, ]$L_GENE), decreasing = TRUE)
sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC > 0, ]$R_GENE), decreasing = TRUE)
sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC > 0, ]$LR_GENES), decreasing = TRUE)

sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC < 0, ]$L_GENE), decreasing = TRUE)
sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC < 0, ]$R_GENE), decreasing = TRUE)
sort(table(test[LR_KEEP_DIFF == TRUE & LR_LOGFC < 0, ]$LR_GENES), decreasing = TRUE)

#Apoe
test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC > 0,]
test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC < 0,]

sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC > 0,]$TISSUE), decreasing = TRUE)
sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC < 0,]$TISSUE), decreasing = TRUE)

sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC > 0,]$L_CELLTYPE), decreasing = TRUE)
sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC < 0,]$L_CELLTYPE), decreasing = TRUE)

sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC > 0,]$R_GENE), decreasing = TRUE)
sort(table(test[L_GENE == "Apoe" & LR_KEEP_DIFF == TRUE & LR_LOGFC < 0,]$R_GENE), decreasing = TRUE)


#############



ora_facs <- ora_results$tms_facs
ora_LR_facs_all <- ora_facs[Tissue == "All" & Category == "LR_GENES"]
ora_LR_facs_all_up <- ora_LR_facs_all[pval_adjusted_UP <= 0.05 & OR_UP >= 1, c("Value", "pval_adjusted_UP", "OR_UP")]
ora_LR_facs_all_down <- ora_LR_facs_all[pval_adjusted_DOWN <= 0.05 & OR_DOWN >= 1, c("Value", "pval_adjusted_DOWN", "OR_DOWN")]


ora_LR_facs_tiss <- ora_facs[Tissue != "All" & Category == "LR_GENES"]
ora_LR_facs_tiss_up <- ora_LR_facs_tiss[pval_adjusted_UP <= 0.5 & OR_UP >= 1, c("Tissue", "Value", "pval_adjusted_UP", "OR_UP")]
ora_LR_facs_tiss_down <- ora_LR_facs_tiss[pval_adjusted_DOWN <= 0.5 & OR_DOWN >= 1, c("Tissue", "Value", "pval_adjusted_DOWN", "OR_DOWN")]


test2 <- test[Tissue == "All" & Category == "LR_GENES" & pval_adjusted_UP <= 0.05 & OR_UP >= 1]
test3 <- test[Tissue == "All" & Category == "LR_GENES" & pval_adjusted_UP <= 0.05 & OR_UP <= 1]

test4 <- test[Category == "LR_GENES" & pval_adjusted_UP <= 0.05 & OR_UP >= 1 & Tissue != "All"]

test5 <- test[Category == "LR_CELLTYPES"  & Tissue == "Liver"]


fpm_test <- fpm_results$tms_facs$sub1
fpm_test2 <- fpm_results$tms_facs$sub2


table(datasets_filtered$tms_facs$CASE_TYPE)
table(datasets_filtered$tms_droplet$CASE_TYPE)
table(datasets_filtered$calico$CASE_TYPE)

cluster_test <- datasets_filtered$tms_droplet[ TISSUE == "Thymus"]
cluster_test[, SIG_VAL := ifelse(CASE_TYPE %in% c("TTFU", "TTFD"),
                                 0,
                                 ifelse(CASE_TYPE %in% c("TTTD", "TFTD"),
                                        -1,
                                        1))]
cluster_test <- dcast(cluster_test[, c("LR_CELLTYPES", "LR_GENES", "SIG_VAL")],
                      formula = LR_GENES ~ LR_CELLTYPES,
                      value.var = "SIG_VAL")
cluster_test[is.na(cluster_test)] <- 0

mat_cluster <- as.matrix(cluster_test, rownames = 1)
mat_cluster <- mat_cluster[rowSums(abs(mat_cluster)) > 4,]
mat_cluster <- mat_cluster[, colSums(abs(mat_cluster)) > 3]

pheatmap(
  mat_cluster,
  show_rownames = TRUE,
  show_colnames = FALSE
)

pheatmap(mat_cluster, #cellwidth=10, cellheight=10,
         width=11, height=11,
         #color = colorRampPalette(brewer.pal(8, "RdYlBu"))(length(breaks)-1),
         #breaks = breaks,
         na_col = "green",
         scale = "none",
         #cluster_rows = nrows >= 2,
         #cluster_cols = ncols >= 2,
         #cutree_rows = ifelse(nrows >= 5, min(nrows, 3), 1),
         #cutree_cols = ifelse(ncols >= 5, min(ncols, 3), 1),
         display_numbers=FALSE, 
         fontsize=10, 
         #main = title,
         #ylab = "Transmitter cell", 
         #xlab = "Receiver cell",
         legend=TRUE,
         #filename = filename
)

LRall[SYMB_LR == "B2m_Cd3g"]


b2m_facs <- datasets_filtered$tms_facs[L_GENE == "B2m"]

table(datasets_filtered$tms_facs$CASE_TYPE)


ggplot(datasets_filtered$tms_facs, aes(x = LR_LOGFC*log2(exp(1)), y = -log10(BH_PVAL_DIFF + 1E-4), color = TISSUE)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = log2(1.1)) +
  geom_vline(xintercept = -log2(1.1)) +
  geom_vline(xintercept = log2(1.5)) +
  geom_vline(xintercept = -log2(1.5)) +
  xlab(expression(paste(Log[2], "FC"))) +
  ylab(expression(paste(-Log[10], " ", p[BH])))
#geom_density_2d() +

volcano_facs <- cowplot::plot_grid(
  plotlist = lapply(
    unique(datasets_filtered$tms_facs$TISSUE),
    function(tiss) {
      ggplot(datasets_filtered$tms_facs[TISSUE == tiss],
             aes(x = LR_LOGFC*log2(exp(1)), y = -log10(BH_PVAL_DIFF + 1E-4)), color = CASE_TYPE) +
        geom_point() +
        geom_hline(yintercept = -log10(0.05)) +
        geom_vline(xintercept = log2(1.1)) +
        geom_vline(xintercept = -log2(1.1)) +
        geom_vline(xintercept = log2(1.5)) +
        geom_vline(xintercept = -log2(1.5)) +
        xlab(expression(paste(Log[2], "FC"))) +
        ylab(expression(paste(-Log[10], " ", p[BH]))) +
        theme(text=element_text(size=20))
    }
  ),
  ncol = 5,
  align = "v",
  labels = unique(datasets_filtered$tms_facs$TISSUE)
)

volcano_droplet <- cowplot::plot_grid(
  plotlist = lapply(
    unique(datasets_filtered$tms_droplet$TISSUE),
    function(tiss) {
      ggplot(datasets_filtered$tms_droplet[TISSUE == tiss],
             aes(x = LR_LOGFC*log2(exp(1)), y = -log10(BH_PVAL_DIFF + 1E-4)), color = CASE_TYPE) +
        geom_point() +
        geom_hline(yintercept = -log10(0.05)) +
        geom_vline(xintercept = log2(1.1)) +
        geom_vline(xintercept = -log2(1.1)) +
        geom_vline(xintercept = log2(1.5)) +
        geom_vline(xintercept = -log2(1.5)) +
        xlab(expression(paste(Log[2], "FC"))) +
        ylab(expression(paste(-Log[10], " ", p[BH]))) +
        theme(text=element_text(size=20))
    }
  ),
  ncol = 4,
  align = "v",
  labels = unique(datasets_filtered$tms_droplet$TISSUE)
)

ggsave("../../../../../volcano_test_facs.png", plot = volcano_facs, scale =  2.5)

ggplot(datasets_filtered$tms_droplet, aes(x = LR_LOGFC*log10(exp(1)), y = -log10(BH_PVAL_DIFF + 1E-4))) + geom_point() +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = log(1.1)) +
  geom_vline(xintercept = -log(1.1)) +
  geom_vline(xintercept = log(1.5)) +
  geom_vline(xintercept = -log(1.5)) +
  geom_density_2d()



