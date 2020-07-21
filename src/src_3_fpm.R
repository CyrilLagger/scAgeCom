# Frequent itemsets analysis ---------------------------------------------------

#' analyze_FreqItemSets
#'
#' @param data 
#' @param cols 
#' @param items_to_include 
#' @param target 
#' @param support 
#' @param confidence 
#'
#' @return
#' @export
#'
#' @examples
analyze_FreqItemSets <- function(data, cols, 
                                 items_to_include = NULL,
                                 target = "closed frequent itemsets",
                                 support = 0.05,
                                 confidence = 0.1) {
  
  transactions = convert_data_to_transactions(data, cols)
  
  res = compute_freqitemsets_and_rules(transactions, 
                                       target,
                                       support,
                                       confidence)
  
  freqsets = res[["freqsets"]]
  rules = res[["rules"]]
  
  interesting_rules = get_interesting_rules(rules, transactions, cols)
  
  return(interesting_rules)
}

#' get_interesting_rules
#'
#' @param rules 
#' @param transactions 
#'
#' @return
#' @export
#'
#' @examples
get_interesting_rules <- function(rules, transactions, cols) {
  
  STRING_sign_diff_expression = glue("{s}=Sign", s=cols$DIFFERENTIAL_EXPRESSED)
  STRING_tissue = glue("{s}=", s=cols$TISSUES)
  STRING_ligand_gene = glue("{s}=", s=cols$L_GENE)
  STRING_receptor_gene = glue("{s}=", s=cols$R_GENE)
  STRING_ligand_celltype = glue("{s}=", s=cols$LIGAND_CELLTYPE)
  STRING_receptor_celltype = glue("{s}=", s=cols$RECEPTOR_CELLTYPE)
  
  sub1 = subset(
    rules, 
    subset = rhs %in% STRING_sign_diff_expression
    & lhs %pin% STRING_tissue
    & lhs %pin% STRING_ligand_gene
    & lhs %pin% STRING_receptor_celltype
    # & !(lhs %pin% "Direction=")
  )
  sub1_measure = interestMeasure(sub1,
                                 measure=c("fishersExactTest", "oddsRatio"),
                                 transactions=transactions,
                                 reuse=TRUE)
  sub1_dt = data.table(
    lhs=labels(lhs(sub1)),
    rhs=labels(rhs(sub1)), 
    sub1@quality,
    sub1_measure
  )
  
  
  sub2 = subset(
    rules, 
    subset = rhs %in% STRING_sign_diff_expression
    & lhs %pin% STRING_tissue
    & lhs %pin% STRING_receptor_gene
    & lhs %pin% STRING_ligand_celltype
    # & !(lhs %pin% "Direction=")
  )
  sub2_measure = interestMeasure(sub2,
                                 measure=c("fishersExactTest", "oddsRatio"),
                                 transactions=transactions,
                                 reuse=TRUE)
  sub2_dt = data.table(
    lhs=labels(lhs(sub2)),
    rhs=labels(rhs(sub2)), 
    sub2@quality,
    sub2_measure
  )
  
  subsets = list()
  subsets[["sub1"]] = sub1_dt
  subsets[["sub2"]] = sub2_dt
  
  return(subsets)
}

#' compute_freqitemsets_and_rules
#'
#' @param transactions 
#' @param support 
#' @param confidence 
#'
#' @return
#' @export
#'
#' @examples
compute_freqitemsets_and_rules <- function(transactions,
                                           target,
                                           support,
                                           confidence) {
  
  freqsets = apriori(
    transactions, 
    parameter = list(support = support, 
                     minlen = 1,
                     maxlen = 5,
                     target = target, 
                     ext = TRUE,
                     smax = 1,
                     minval = 0, 
                     maxtime = 0),
    control = list(verbose=FALSE)
  )
  
  rules = ruleInduction(
    freqsets, 
    transactions = transactions,
    confidence = confidence,
    control = list(method = "ptree", verbose = FALSE)
  )
  
  return(list(freqsets = freqsets, rules = rules))
}

#' convert_data_to_transactions
#'
#' @param data 
#' @param cols 
#'
#' @return
#' @export
#'
#' @examples
convert_data_to_transactions <- function(data, cols) {
  
  COL_DIFFERENTIAL_EXPRESSED = cols$DIFFERENTIAL_EXPRESSED
  
  cols_for_items = c(
    cols$TISSUES, 
    cols$LIGAND_CELLTYPE, cols$RECEPTOR_CELLTYPE, 
    cols$L_GENE, cols$R_GENE, 
    cols$DIFFERENTIAL_EXPRESSED, cols$DIFFERENTIAL_DIRECTION
  )
  
  dt_p = data[, ..cols_for_items]
  
  dt_p[, (COL_DIFFERENTIAL_EXPRESSED) := 
         fifelse(get(COL_DIFFERENTIAL_EXPRESSED), "Sign", "Notsign")]
  
  
  transactions = as(dt_p, "transactions")
  return(transactions)
}
