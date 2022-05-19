library(data.table)
source("utils_random_pairs.R")

get_proteomics_dt = function() {
  dt = data.table(
    read.csv("data/proteins/Mouse_lung_proteome_study.csv")
  )
  
  setnames(dt, "Student.s.T.test.Difference_old.young..log2.", "log2FC")
  setnames(dt, "Student.s.T.test.Significant.24m", "is_significant")
  setnames(dt, "X.Log.Student.s.T.test.p.value.24m", "neg_log_p_val")
  setnames(dt, "Student.s.T.test.q.value.24m", "q_val")
  setnames(dt, "Student.s.T.test.Test.statistic.24m", "t_statistic")
  setnames(dt, "Gene.names", "gene")
  
  dt = dt[, list(
    gene, log2FC, t_statistic, neg_log_p_val, q_val, is_significant
  )]
  dt[, is_significant := is_significant == "+"]
  
  # Validater
  all(
    (dt$q_val <= 0.1) == (dt$is_significant)
  )
  return(dt)
}


lri = scDiffCom::LRI_mouse$LRI_curated
lri = subset_simple_lri(lri)
lri_genes = get_genes_lri(lri)

AGE_GROUP=list(
  YOUNG = c("3m"),
  OLD = c("18m", "21m", "24m")
)
s = get_seurat("droplet", "Lung", AGE_GROUP, "female")
res_scDiffCom = run_internal_analysis(
  seurat_obj = s,
  lri_table = lri
)
res_scDiffCom = scDiffCom::FilterCCI(res_scDiffCom, skip_ora = TRUE)
det = res_scDiffCom@cci_table_detected

# Proteomics
dt = get_proteomics_dt()
dt$in_lri = dt$gene %in% lri_genes
sum(dt$gene %in% lri_genes)
contingency = table(dt[, list(is_significant, in_lri)])
fisher = fisher.test(contingency)
dt_sel = dt[is_significant == TRUE & in_lri == TRUE, list(gene, log2FC, q_val)]
# Can recover 17 genes by splitting ";", not worth it.
# recoverable = purrr::map_lgl(dt$gene, function(gene_str) {
#   split = strsplit(gene_str, ";")[[1]]
#   if (length(split) == 0) {
#     return(FALSE)
#   }
#   for (gene in split) {
#     if (gene %in% lri_genes) {
#       return(TRUE)
#     }
#   return(FALSE)
#   }
# })
# sum(recoverable)

# Integrate CCI with Proteomics
det[, IS_LIGAND_1_IN_DT := LIGAND_1 %in% dt$gene]
det[, IS_RECEPTOR_1_IN_DT := RECEPTOR_1 %in% dt$gene]
det[, IS_CCI_IN_DT_SOFT := IS_LIGAND_1_IN_DT | IS_RECEPTOR_1_IN_DT]
det[, IS_CCI_IN_DT_HARD := IS_LIGAND_1_IN_DT & IS_RECEPTOR_1_IN_DT]

get_pvals = function(genes) {
  purrr::map_dbl(genes, function(gene_name) {
    if (gene_name %in% dt$gene) {
      qvals = dt[gene == gene_name, q_val]
      if (length(qvals) > 1) {
        return(max(qvals))
      } else {
        return(qvals)
      }
    } else {
      return(NA)
    }
  })
}

get_log2fc = function(genes) {
  purrr::map_dbl(genes, function(gene_name) {
    if (gene_name %in% dt$gene) {
      log2fcs = dt[gene == gene_name, log2FC]
      if (length(log2fcs) > 1) {
        return(min(log2fcs))
      } else {
        return(log2fcs)
      }
    } else {
      return(NA)
    }
  })
}

ligand_1_pvals = get_pvals(det[, LIGAND_1])
receptor_1_pvals = get_pvals(det[, RECEPTOR_1])
det[, LIGAND_1_DT_P_VAL := ligand_1_pvals]
det[, RECEPTOR_1_DT_P_VAL := receptor_1_pvals]

ligand_1_log2fcs = get_log2fc(det[, LIGAND_1])
receptor_1_log2fcs = get_log2fc(det[, RECEPTOR_1])
det[, LIGAND_1_DT_LOG2FC := ligand_1_log2fcs]
det[, RECEPTOR_1_DT_LOG2FC := receptor_1_log2fcs]


## Analyze
dt_regulation = purrr::pmap_chr(
  det[, list(IS_LIGAND_1_IN_DT,
             LIGAND_1_DT_P_VAL,
             LIGAND_1_DT_LOG2FC,
             IS_RECEPTOR_1_IN_DT,
             RECEPTOR_1_DT_P_VAL,
             RECEPTOR_1_DT_LOG2FC)],
  function(
    IS_LIGAND_1_IN_DT,
    LIGAND_1_DT_P_VAL,
    LIGAND_1_DT_LOG2FC,
    IS_RECEPTOR_1_IN_DT,
    RECEPTOR_1_DT_P_VAL,
    RECEPTOR_1_DT_LOG2FC
    ) {
    if (IS_LIGAND_1_IN_DT & IS_RECEPTOR_1_IN_DT) {
      if ( (LIGAND_1_DT_P_VAL >= 0.1) & (RECEPTOR_1_DT_P_VAL >= 0.1) ) {
        return("FLAT")
      } else if (LIGAND_1_DT_P_VAL < 0.01) {
        if (LIGAND_1_DT_LOG2FC > 0) {
          return("UP_LIGAND")
        } else {
          return("DOWN_LIGAND")
        }
      } else if (RECEPTOR_1_DT_P_VAL < 0.01) {
        if (RECEPTOR_1_DT_LOG2FC > 0) {
          return("UP_RECEPTOR")
        } else {
          return("DOWN_RECEPTOR")
        }
      } else {
        if ( (LIGAND_1_DT_LOG2FC > 0) & (RECEPTOR_1_DT_LOG2FC > 0) ) {
          return("UP_BOTH")
        } else if ( (LIGAND_1_DT_LOG2FC < 0) & (RECEPTOR_1_DT_LOG2FC < 0) ) {
          return("DOWN_BOTH")
        } else {
          "UNCLEAR"
        }
      }
    } else {
      return(NA)
    }
  }
)
det[, DT_REGULATION := dt_regulation]

table(det[, list(REGULATION, DT_REGULATION)])

