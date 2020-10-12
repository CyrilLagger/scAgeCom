
analyze_ORA <- function(
  data
) {
  dt_ora = ora(data)
  dt_ora_up = ora(data[DIFFERENTIAL_DIRECTION == "UP"])
  dt_ora_down = ora(data[DIFFERENTIAL_DIRECTION == "DOWN"])
  dt_complete = merge(dt_ora, dt_ora_up,
                      by = c("Tissue", "Category", "Value"), 
                      all = TRUE, 
                      suffixes = c("", "_UP"))
  dt_complete = merge(dt_complete, dt_ora_down,
                      by = c("Tissue", "Category", "Value"),
                      all = TRUE,
                      suffixes = c("", "_DOWN"))
  return(dt_complete)
}

ora <- function(data) {
  categories = c(
    "TISSUE",
    "L_CELLTYPE", "R_CELLTYPE", "LR_CELLTYPE",
    "LR_NAME"
  )
  dt_counts = fast_counts(data, categories)
  dt_ora = perform_ora_from_counts(dt_counts)
  return(dt_ora)
}

fast_counts <- function(data, categories) {
  COL_TISSUE <- "TISSUE"
  dt_categories = count_on_all_columns(data, categories)
  dt_tissues = count_significant_on_tissues(data)
  dt_counts = merge(dt_categories, dt_tissues, by=COL_TISSUE, all=TRUE)
    dt_counts = dt_counts[, .(
    Tissue = get(COL_TISSUE),
    Category = Category,
    Value = Value,
    Counts_value_significant = Counts_value_significant,
    Counts_value_notsignificant = Counts_value_notsignificant,
    Counts_notvalue_significant = Counts_significant - Counts_value_significant,
    Counts_notvalue_notsignificant = Counts_notsignificant - Counts_value_notsignificant
  )
  ]
  return(dt_counts)
}

count_on_all_columns <- function(data, categories) {
  dt_categories = list()
  for (category in categories) {
    dt_cat = count_on_one_column(data, category)
    dt_categories[[category]] = dt_cat
  }
  dt_counts = rbindlist(dt_categories, use.names = TRUE)
  return(dt_counts)
}

count_on_one_column <- function(data, column) {
  COL_TISSUE <- "TISSUE"
  dt_val_sig = data[DIFFERENTIAL_EXPRESSED == TRUE,
                    .(Category=column, Counts_value_significant=.N),
                    by=.(get(COL_TISSUE), Value=get(column))]
  dt_val_notsig = data[DIFFERENTIAL_EXPRESSED == FALSE,
                       .(Category=column, Counts_value_notsignificant=.N),
                       by=.(get(COL_TISSUE), Value=get(column))]
  dt_list = list(dt_val_sig, dt_val_notsig)
  lapply(dt_list, function(dt) {
    setnames(dt, new = COL_TISSUE, old = "get")
  })
  dt_cat_final = merge(dt_list[[1]], dt_list[[2]], 
                       by=c(COL_TISSUE, "Value", "Category"),
                       all=TRUE)
  setnafill(dt_cat_final, 
            "const",
            fill=0,
            cols=4:5)
  # Add entries for all organs
  dt_all = dt_cat_final[, lapply(.SD, sum), by=.(Value, Category), .SDcols = !c(COL_TISSUE)]
  dt_all[, (COL_TISSUE) := "All"]
  dt_cat_final = rbindlist(list(dt_cat_final, dt_all), use.names=TRUE)
  return(dt_cat_final)
}

count_significant_on_tissues <- function(data) {
  COL_TISSUE <- "TISSUE"
  # To be merged at last on tissue level + to create "All" entry
  dt_sig = data[DIFFERENTIAL_EXPRESSED == TRUE,
                .(Counts_significant = .N),
                by = .(get(COL_TISSUE))]
  dt_notsig = data[DIFFERENTIAL_EXPRESSED == FALSE,
                   .(Counts_notsignificant = .N),
                   by = .(get(COL_TISSUE))]
  dt_list = list(dt_sig, dt_notsig)
  lapply(dt_list, function(dt) {
    setnames(dt, new = COL_TISSUE, old = "get")  })
  
  dt_tissues_counts = merge(dt_list[[1]], dt_list[[2]], 
                            by=c(COL_TISSUE),
                            all=TRUE)
  setnafill(dt_tissues_counts, "const", fill=0, cols=2:3)
  dt_all = dt_tissues_counts[, lapply(.SD, sum), .SDcols = !c(COL_TISSUE)]
  dt_all[, (COL_TISSUE) := "All"]
  dt_tissues_counts = rbindlist(list(dt_tissues_counts, dt_all), use.names=TRUE)
  return(dt_tissues_counts)
}

perform_ora_from_counts <- function(data) {
  
  # https://stackoverflow.com/questions/11680579/assign-multiple-columns-using-in-data-table-by-group
  
  data[, c("OR", "pval") := 
         vfisher_2sided(Counts_value_significant,
                        Counts_value_notsignificant,
                        Counts_notvalue_significant, 
                        Counts_notvalue_notsignificant)
       ]
  
  data[, c("Kulc_distance", "Imbalance_ratio") := list(
    kulc(Counts_value_significant, 
         Counts_value_notsignificant, 
         Counts_notvalue_significant, 
         Counts_notvalue_notsignificant),
    imbalance_ratio(Counts_value_significant, 
                    Counts_value_notsignificant, 
                    Counts_notvalue_significant,
                    Counts_notvalue_notsignificant)
  )]
  
  # Add adjusted pval
  data[, "pval_adjusted" := .(p.adjust(pval, method="BH"))]
  
  return(data)
}

odds_ratio <- function(Counts_value_significant,
                       Counts_value_notsignificant,
                       Counts_notvalue_significant,
                       Counts_notvalue_notsignificant) {
  
  # There are functions that have better methods for estimating
  #  odds ratio rather than the simple formula below.
  # uses only basic arithmetic operations which are already vectorized
  # or = ((Counts_value_significant * Counts_notvalue_notsignificant)
  #       / (Counts_value_notsignificant * Counts_notvalue_significant))
  
  or = 
    
    return(or)
}


fisher_2sided <- function(Counts_value_significant,
                          Counts_value_notsignificant,
                          Counts_notvalue_significant,
                          Counts_notvalue_notsignificant) {
  
  m = matrix(nrow = 2, ncol = 2, byrow=TRUE,
             data = c(Counts_value_significant, Counts_value_notsignificant,
                      Counts_notvalue_significant, Counts_notvalue_notsignificant)
  )
  
  test = fisher.test(m, alternative="two.sided")
  res = c(test$estimate, test$p.value)
  return(res)
}

vfisher_2sided <- function(Counts_value_significant,
                           Counts_value_notsignificant,
                           Counts_notvalue_significant,
                           Counts_notvalue_notsignificant) {
  
  v = mapply(fisher_2sided, 
             Counts_value_significant,
             Counts_value_notsignificant,
             Counts_notvalue_significant,
             Counts_notvalue_notsignificant)
  
  l = list(OR = v[1, ], pval = v[2, ])
  return(l)
}

kulc <- function(Counts_value_significant,
                 Counts_value_notsignificant,
                 Counts_notvalue_significant,
                 Counts_notvalue_notsignificant) {
  
  # vectorized by using only arithmetic operations
  P_AB = Counts_value_significant / (Counts_value_significant + Counts_notvalue_significant)
  P_BA = Counts_value_significant / (Counts_value_significant + Counts_value_notsignificant)
  avg = (P_AB + P_BA) / 2
  return(avg)
}


imbalance_ratio <- function(Counts_value_significant,
                            Counts_value_notsignificant,
                            Counts_notvalue_significant,
                            Counts_notvalue_notsignificant) {
  # Vectorized
  numerator = abs(Counts_value_notsignificant - Counts_notvalue_significant)
  denominator = (Counts_value_significant 
                 + Counts_value_notsignificant
                 + Counts_notvalue_significant)
  return(numerator / denominator)
  
}



# Plots ------------------------------------------------------------------------

plot_heatmaps <- function(data, cols) {
  
  
  
}

build_heatmaps <- function(dt, cols, dir_path) {
  
  # add "/" if not last char of dir_path
  if (!substr(dir_path, length(dir_path), length(dir_path)) == "/") {
    dir = paste0(dir_path, "/")
  } else {
    dir = dir_path
  }
  
  dir.create(dir, showWarnings = FALSE) 
  
  OR_max = 10
  OR_min = 0.1
  tissues = unique(dt[["Tissue"]])
  ncols = dim(dt)[2]
  
  for (tissue in tissues) {
    
    dt_edge = dt[Tissue == tissue, 2:ncols]
    nrows = dim(dt_edge)[1]
    if (nrows == 1) {message(paste0("One edge for: ", tissue))}
    
    G = graph_from_data_frame(d = dt_edge,
                              directed = TRUE,
                              vertices = NULL)
    
    matrix_all = as_adjacency_matrix(G, attr = "OR", sparse = FALSE)
    matrix_all = clip_matrix(matrix_all, 1, OR_min, OR_max)
    
    matrix_up = as_adjacency_matrix(G, attr = "OR_UP", sparse = FALSE)
    matrix_up = clip_matrix(matrix_up, 1, OR_min, OR_max)
    
    matrix_down = as_adjacency_matrix(G, attr = "OR_DOWN", sparse = FALSE)
    matrix_down = clip_matrix(matrix_down, 1, OR_min, OR_max)
    
    dir.create(paste0(dir, tissue, '/'), showWarnings = FALSE)
    
    explanatory_string = paste0("Row-transmitter, col-receiver; ",
                                "Values are clipped log(OR)")
    eps = 10**(-5)
    generate_heatmap(log(matrix_all + eps),
                     #matrix_all,
                     title = paste0(tissue, ": ", "ALL ", explanatory_string),
                     filename = paste0(dir, tissue, '/', 'ALL.png'))
    generate_heatmap(log(matrix_up + eps), 
                     #matrix_up,
                     title = paste0(tissue, ": ", "UP ", explanatory_string),
                     filename = paste0(dir, tissue, '/', 'UP.png'))
    generate_heatmap(log(matrix_down + eps),
                     #matrix_down,
                     title = paste0(tissue, ": ", "DOWN ", explanatory_string),
                     filename = paste0(dir, tissue, '/', 'DOWN.png'))
  }
}



clip_matrix <- function(matrix, zero_val, min_val, max_val) {
  matrix[matrix == 0] = zero_val
  matrix[matrix < min_val] = min_val
  matrix[matrix > max_val] = max_val
  return(matrix)
}


generate_heatmap <- function(matrix, title, filename) {
  
  message(paste0("Preparing figure: ", title, " ", filename))
  
  nrows = dim(matrix)[1]
  ncols = dim(matrix)[2]
  
  if (nrows < 2 | ncols < 2) {
    message(paste0("Matrix size is < 2x2: ", title))
    return()
  }
  
  breaks = c(seq(-3, -0.5, 0.5), 0, seq(0.5, 3, 0.5))
  pheatmap(matrix, #cellwidth=10, cellheight=10,
           width=11, height=11,
           color = colorRampPalette(brewer.pal(8, "RdYlBu"))(length(breaks)-1),
           breaks = breaks,
           na_col = "green",
           scale = "none",
           cluster_rows = nrows >= 2,
           cluster_cols = ncols >= 2,
           cutree_rows = ifelse(nrows >= 5, min(nrows, 3), 1),
           cutree_cols = ifelse(ncols >= 5, min(ncols, 3), 1),
           display_numbers=FALSE, 
           fontsize=10, 
           main = title,
           ylab = "Transmitter cell", 
           xlab = "Receiver cell",
           legend=TRUE,
           filename = filename)
}


get_celltypes_enrichment <- function(dt, cols) {
  
  COL_LIGAND_RECEPTOR_CELLTYPES = cols$LIGAND_RECEPTOR_CELLTYPES
  
  dt_ctypes = dt[Category == COL_LIGAND_RECEPTOR_CELLTYPES
                 & pval < 0.05]
  dt_ctypes = dt_ctypes[, c("Ligand_cell", "Receptor_cell") := list(
    sub("(.*) : ", "", Value),
    sub(" : (.*)", "", Value)
  )
  ]
  
  cols_to_select = c(
    "Tissue", "Ligand_cell", "Receptor_cell", "OR", "OR_UP", "OR_DOWN"
  )
  
  dt_ctypes = dt_ctypes[, ..cols_to_select]
  
  return(dt_ctypes)
}






