####################################################
##
## Project: scAgeCom
## 
## October 2020
##
## cyril.lagger@liverpool.ac.uk
## anais.equey@etu.univ-amu.fr
##
## Script to annotate the (curated) ligand-receptor
## interactions.
##
####################################################
##
## SUMMARY
## 1 to 3 - we start by listing all the genes we have in our database
## 4 - then we find the go terms associated with each of our genes
##     (based on ensembl)
## 4 - For each go term, we also add its parents (based on Gene Ontology)
## 5 - Each interaction involves several genes. We take both the union and the
##     intersection of their go terms and add the result to separate columns
##
####################################################

#--------- 1. Load libraries --------#

library(scDiffCom)
library(data.table)
library(ontologyIndex)
library(ontoProc)
library(biomaRt)

#--------- 2. Load the curated ligand-receptor interactions --------#

LR_interactions <- scDiffCom::LR6db$LR6db_curated

#--------- 3. Get all unique genes in the data.table --------#

# Ligand genes
L_genes <- unique(unlist(LR_interactions[, c("LIGAND_1", "LIGAND_2")]))
L_genes <- L_genes[!is.na(L_genes)] # remove NA (empty) values
# Receptor genes
R_genes <- unique(unlist(LR_interactions[, c("RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
R_genes <- R_genes[!is.na(R_genes)] # remove NA (empty) values
# All genes (some genes may appear as both ligand and receptors)
LR_genes <- unique(c(L_genes, R_genes))

#--------- 4. Get the Ensembl terms (GO, family or else) --------#

# create an access to ensembl database via biomart (might take a minute or so)
mart <- biomaRt::useMart(
  "ensembl",
  dataset = "mmusculus_gene_ensembl"
)
# retrieve the information we want for our genes of interest
LR_genes_info <- biomaRt::getBM(
  attributes = c(
    "mgi_symbol",
    "go_id",
    "name_1006",
    #"definition_1006",
    "namespace_1003",
    "family",
    "family_description"
    #"goslim_goa_accession",
    #"goslim_goa_description"
  ),
  filters = "mgi_symbol",
  mart = mart,
  values = LR_genes
)
setDT(LR_genes_info)

#--------- 5. adding ancestors from Gene Ontology --------#

# Load go terms from gene ontology
onto_go_terms = getGeneOnto()
go_names <- onto_go_terms$name
LR_genes_go <- sapply(
  LR_genes,
  function(gene) {
    temp_go <- unique(LR_genes_info[mgi_symbol == gene]$name_1006)
    get_ancestors(onto_go_terms, names(go_names[go_names %in% temp_go]))
  },
  USE.NAMES = TRUE,
  simplify = FALSE
)

#--------- 6. Create table of go terms associated to each LR-interaction --------#

LR_interactions_go_union <- rbindlist(
  apply(
    LR_interactions,
    MARGIN = 1,
    function(row) {
      LIGAND_GO <- unique(c(
        LR_genes_go[[row[["LIGAND_1"]]]],
        LR_genes_go[[row[["LIGAND_2"]]]]
      ))
      RECEPTOR_GO <- unique(c(
        LR_genes_go[[row[["RECEPTOR_1"]]]],
        LR_genes_go[[row[["RECEPTOR_2"]]]],
        LR_genes_go[[row[["RECEPTOR_3"]]]]
      ))
      res_union <- unique(c(LIGAND_GO, RECEPTOR_GO))
      res_union <- data.table(
        LR_SORTED = rep(row[["LR_SORTED"]], length(res_union)),
        GO_union = res_union
      )
    }
  )
)

LR_interactions_go_intersection <- rbindlist(
  apply(
    LR_interactions,
    MARGIN = 1,
    function(row) {
      LIGAND_GO <- unique(c(
        LR_genes_go[[row[["LIGAND_1"]]]],
        LR_genes_go[[row[["LIGAND_2"]]]]
      ))
      RECEPTOR_GO <- unique(c(
        LR_genes_go[[row[["RECEPTOR_1"]]]],
        LR_genes_go[[row[["RECEPTOR_2"]]]],
        LR_genes_go[[row[["RECEPTOR_3"]]]]
      ))
      res_inter <- intersect(LIGAND_GO, RECEPTOR_GO)
      if(length(res_inter) > 0) {
        res_inter <- data.table(
          LR_SORTED = rep(row[["LR_SORTED"]], length(res_inter)),
          GO_intersection = res_inter
        )
      } else {
        res_inter <- NULL
      }
      return(res_inter)
    }
  )
)
