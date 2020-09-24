

analyze_FreqItemSets <- function(
  data,
  items_to_include = NULL,
  target = "closed frequent itemsets",
  support = 0.05,
  confidence = 0.1
) {
  transactions = convert_data_to_transactions(data)
  res = compute_freqitemsets_and_rules(
    transactions, 
    target,
    support,
    confidence
  )
  
  freqsets = res[["freqsets"]]
  rules = res[["rules"]]
  
  interesting_rules = get_interesting_rules(rules, transactions)
  
  return(interesting_rules)
}


get_interesting_rules <- function(rules, transactions, cols) {
  
  STRING_sign_diff_expression = glue("{s}=Sign", s= "DIFFERENTIAL_EXPRESSED")
  STRING_tissue = glue("{s}=", s= "TISSUE")
  STRING_ligand_gene = glue("{s}=", s= "LIGAND_1")
  STRING_receptor_gene = glue("{s}=", s= "RECEPTOR_1")
  STRING_ligand_celltype = glue("{s}=", s= "L_CELLTYPE")
  STRING_receptor_celltype = glue("{s}=", s= "R_CELLTYPE")
  
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


compute_freqitemsets_and_rules <- function(
  transactions,
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


convert_data_to_transactions <- function(data) {
    cols_for_items = c(
    "TISSUE", 
    "L_CELLTYPE", "R_CELLTYPE",
    "LIGAND_1", "RECEPTOR_1", 
    "DIFFERENTIAL_EXPRESSED",
    "DIFFERENTIAL_DIRECTION"
  )
  dt_p = data[, ..cols_for_items]
  dt_p[, DIFFERENTIAL_EXPRESSED := 
         fifelse(DIFFERENTIAL_EXPRESSED, "Sign", "Notsign")]
  transactions = as(dt_p, "transactions")
  return(transactions)
}
