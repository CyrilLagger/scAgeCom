library(UpSetR)
library(ComplexUpset)
library(scDiffCom)
library(ggplot2)

#load db  
LRdb <- LR6db$LR6db_all
db_sources <- LR6db$LR6db_source
db_names <- names(db_sources)

##### look at which source to remove or to keep ####
#from sctensor: remove all if unique
LR_rm_sctensor <- db_sources$SCTENSOR
#from nichenet: remove ppi (keep kegg_*, pharmacology, ramilowski_known)
LR_rm_nichenet <- db_sources$NICHENET[grepl("ppi", db_sources$NICHENET)]
#from icellnet: keep all (comes only as PMID)
#from cellchat: keep all (comes from PMID, KEGG and PMC)
sum(!grepl("PMID|KEGG|PMC", db_sources$CELLCHAT))
#from cellphonedb: keep all (comes as cpdb)
#from scsr: keep all for now
LR_scsr_source <- sort(unique(unlist(strsplit(db_sources$SCSR, ","))))
LR_scsr_source <- unique(gsub('[[:digit:]]+', '', LR_scsr_source))
LR_scsr_source_manual <- c("cellsignal.com", "fantom5", "HPMR", "HPRD", "IUPHAR", "literature", "PMID", "reactome", "uniprot")
LR_scsr_source_manual2 <- c("cellsignal\\.com", "fantom5", "HPMR", "HPRD", "IUPHAR", "literature", "PMID", "reactome", "uniprot")
sum(!grepl(paste0(LR_scsr_source_manual2, collapse = "|"), db_sources$SCSR))
LR_rm_scsr <- c("uniprot")

#final list to remove
LR_rm <- c(
  LR_rm_sctensor,
  LR_rm_nichenet,
  LR_rm_scsr,
  sapply(LR_rm_sctensor, function(i) {
    sapply(LR_rm_nichenet, function(j) {
      c(paste(i,j, sep = ","), paste(j, i, sep = ","))
    })
  }),
  sapply(LR_rm_sctensor, function(i) {
    sapply(LR_rm_scsr, function(j) {
      c(paste(i,j, sep = ","), paste(j, i, sep = ","))
    })
  }),
  sapply(LR_rm_scsr, function(i) {
    sapply(LR_rm_nichenet, function(j) {
      c(paste(i,j, sep = ","), paste(j, i, sep = ","))
    })
  }),
  sapply(LR_rm_sctensor, function(i) {
    sapply(LR_rm_nichenet, function(j) {
      sapply(LR_rm_scsr, function(k) {
        c(paste(i,j,k, sep = ","), paste(i,k,j, sep = ","), paste(j, i,k, sep = ","), paste(j, k, i, sep = ","),
          paste(k, i, j, sep = ","), paste(k, j, i, sep = ","))
      })
    })
  })
)
LRdb_filtered <- LRdb[!(SOURCE %in% LR_rm)]

#rename source to look at the categories
LRdb_filtered[, SOURCE_CAT := gsub(paste0(c(LR_rm_nichenet, LR_rm_sctensor, "uniprot"), collapse = "|"), "PPI", SOURCE)]
LRdb_filtered[, SOURCE_CAT := gsub("pharmacology", "IUPHAR", SOURCE_CAT)]
LRdb_filtered[, SOURCE_CAT := gsub("kegg", "KEGG", SOURCE_CAT)]
LRdb_filtered[, SOURCE_noDig := gsub(" ", "", gsub('[[:digit:]]+', '', SOURCE))]
temp_char <- LRdb_filtered[3315]$SOURCE_noDig
LRdb_filtered[SOURCE_noDig %in% c("", "; ", " ;", ";", temp_char), SOURCE_CAT := paste0("PMID:", SOURCE_CAT)]
LRdb_filtered[, SOURCE_CAT := gsub("fantom5", "ramilowski", SOURCE_CAT)]

category_sources <- c("PPI", "IUPHAR", "KEGG", "PMID", "CPDB", "ramilowski",
                      "cellsignal.com", "HPMR", "HPRD", "reactome")
sum(!grepl(paste0(category_sources, collapse = "|"), LRdb_filtered$SOURCE_CAT))


source_db_list <- sapply(category_sources, function(i) {
  LRdb_filtered[grepl(i, SOURCE_CAT)]$LR_SORTED
})

UpSetR::upset(fromList(source_db_list), nsets = 9, order.by = "freq", nintersects = 35)

db_list <- sapply(db_names, function(i) {
  LRdb_filtered[grepl(i, DATABASE)]$LR_SORTED
}, 
USE.NAMES = TRUE,
simplify = FALSE)

UpSetR::upset(fromList(db_list), nsets = 6, order.by = "freq", nintersects = 26, scale.intersections = "identity")

#### with ComplexUpset

LR_temp <- LRdb_filtered
LR_temp[, c(category_sources) := lapply(category_sources, function(i) {
  grepl(i, SOURCE_CAT)
})]
LR_temp[, c(db_names) := lapply(db_names, function(i) {
  grepl(i, DATABASE)
})]
LR_temp[, complex := !is.na(LIGAND_2) | !is.na(RECEPTOR_2)]


upset(
  LR_temp,
  category_sources,
  min_size = 20
)

upset(
  LR_temp,
  category_sources,
  min_size = 20,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE,
      aes=aes(fill=complex)
    )
  )
)

upset(
  LR_temp,
  db_names,
  min_size = 20
)


upset(
  LR_temp,
  db_names,
  min_size = 20,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE,
      aes=aes(fill=complex)
    )
  )
)

#look at orthology
#OK:
table(LRdb_filtered$LIGAND_2_CONF)
table(LRdb_filtered$RECEPTOR_3_CONF)
# a few 0
table(LRdb_filtered$LIGAND_1_CONF)
table(LRdb_filtered$RECEPTOR_1_CONF)
table(LRdb_filtered$RECEPTOR_2_CONF)
ftable(LRdb_filtered$LIGAND_1_CONF, LRdb_filtered$RECEPTOR_1_CONF)
ftable(LRdb_filtered$LIGAND_1_CONF, LRdb_filtered$RECEPTOR_1_CONF, LRdb_filtered$RECEPTOR_2_CONF)



library(data.table)
fwrite(LRdb_filtered, file = "../../../../../test.csv")

#annotation with GO
library(clusterProfiler)
library(org.Mm.eg.db)

LR_genes <- unique(unlist(LR6db$LR6db_curated[, c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
LR_genes <- LR_genes[!is.na(LR_genes)]


LR_ggo_bc_l3 <- groupGO(
  gene = LR_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  level = 3,
  readable = FALSE
)@result

LR_ggo_bc_l2 <- groupGO(
  gene = LR_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  level = 2,
  readable = FALSE
)@result

LR_ggo_bc_l4 <- groupGO(
  gene = LR_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  level = 4,
  readable = FALSE
)@result

####

library(biomaRt)

mart <- biomaRt::useMart(
  "ensembl",
  dataset = "mmusculus_gene_ensembl"
)
LR_genes_go <- biomaRt::getBM(
  attributes = c(
    "mgi_symbol",
    "go_id",
    "name_1006", 
    "namespace_1003"
  ),
  filters = "mgi_symbol",
  mart = mart,
  values = LR_genes
)
setDT(LR_genes_go)

LR_genes_go_comb <- LR_genes_go[, list(text = paste(name_1006, collapse=",")), by = mgi_symbol]

LR_genes_go_comb2 <- LR_genes_go[, list(text = paste(mgi_symbol, collapse=",")), by = name_1006]

library(GO.db)
Term("GO:0016021")
t <- c(GOBPOFFSPRING[["GO:0042110"]], "GO:0042110")
Term("GO:0042110")
Term(t[2])

library(data.table)
setDT(LR_genes_go)

LR_genes_go[, count := .N, by = go_id]


LR_genes_go[ go_id == "GO:0050789"]

test <- LRall
test2 <- test[ cpdb == TRUE | scsr == TRUE]
