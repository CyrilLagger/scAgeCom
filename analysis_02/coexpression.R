## Failed installation and contains human co-expression networks
# BiocManager::install("annotate")
# BiocManager::install("sva")
# BiocManager::install("GOSim")
# BiocManager::install("WGCNA") 
# devtools::install_github('juanbot/CoExpNets')
# options(download.file.method = "libcurl")
# devtools::install_github('juanbot/CoExpGTEx')

## Try use CoCoCoNet data
# BiocManager::install("rhdf5")

library(igraph)
library(ggplot2)

read_network = function(
  fpath = "data/mouse_HC_AggNet.hdf5"
) {
  # Fields in hdf5 file
  # https://labshare.cshl.edu/shares/gillislab/resource/CoCoCoNet/networks/readme.txt
  # agg contains the co-expression matrix.
  # row/col contains the row and column names of the matrix (genes).
  # chunkSize contains the chunk size used in writing the hdf5 file.
  # 
  # Fields: 
  #   agg 		- FLOAT 	N x N
  # col			- STR 		N
  # row 		- STR		N
  # chunkSize 	- FLOAT 	1
  
  agg = rhdf5::h5read(fpath, "agg")
  # print(paste0("agg dim: ", dim(agg)))
  row = rhdf5::h5read(fpath, "row")
  col = rhdf5::h5read(fpath, "col")
  if (!all(row == col)) {
    stop("rows and columns are not well-ordered.")
  }
  colnames(agg) = col
  rownames(agg) = row
  # chunkSize = rhdf5::h5read(fpath, "chunkSize")

  return(agg)
}

get_network = function(
  fpath = "data/mouse_HC_AggNet.hdf5",
  coexp_threshold = NULL,
  ngenes = NULL
) {
  
  agg = read_network(fpath)
  
  if (!is.null(coexp_threshold)) {
    agg[agg < coexp_threshold] = 0
  }
  
  if (!is.null(ngenes)) {
    agg = agg[1:ngenes, 1:ngenes]
  }
  
  net = igraph::graph_from_adjacency_matrix(
    agg,
    weighted = TRUE,
    diag = FALSE,
    mode = "undirected"
  )
  return(net)
}

compute_degree_df = function(
  graph
) {
  degrees = igraph::degree(graph)
  df_degrees = data.frame(
    gene = names(degrees),
    degree = degrees
  )
  return(df_degrees)
}

compute_closeness_df = function(
  graph
) {
  closeness = igraph::closeness(
    graph,
    normalized = FALSE
  )
  closeness[is.na(closeness)] = 0
  df_closeness = data.frame(
    gene = names(closeness),
    closeness = closeness
  )
  return(closeness)
}

#' add_lri_coexpression
#'
#' Adds NA for genes or edges implied by the `lri` and not found in `net`.
#'
#' @param lri with LIGAND_1_ID, RECEPTOR_1_ID
#' @param net
#'
#' @return lri with `coexpression` column.
add_lri_coexpression = function(
  lri,
  net
) {
  net_vertex_names = V(net)$name
  lri[, coexpression := mapply(function(L, R) {
    if (!(L %in% net_vertex_names) | !(R %in% net_vertex_names)) {
      return(NA)
    }
    
    if (!are.connected(net, L, R)) {
      return(NA)
    } else {
      return(E(net, c(L, R))$weight)
    }
  }, LIGAND_1_ID, RECEPTOR_1_ID)]
  return(lri)
}


#' add_lri_coexpression_adj
#'
#' Different version for add_lri_coexpression working from adj matrix.
#'
#' @param lri with LIGAND_1_ID, RECEPTOR_1_ID
#' @param adj
#'
#' @return lri with `coexpression` column.
add_lri_coexpression_adj = function(
  lri,
  adj
) {
  adj_ids = colnames(adj)
  lri[, coexpression := mapply(function(L, R) {
    if (!(L %in% adj_ids) | !(R %in% adj_ids)) {
      return(NA)
    }
    return(adj[L, R])
  }, LIGAND_1_ID, RECEPTOR_1_ID)]
  return(lri)
}


get_non_lri_coexpression = function(
  lri,
  adj
) {
  adj_ids = colnames(adj)
  mapply()
}


#' add_scdiffcom_results_coexp
#'
#' @param table_detected 
#' @param lri 
#'
#' @return table_detected with additional `coexpression` column
add_scdiffcom_results_coexp = function(
  table_detected,
  lri
) {
  res_coexp = merge(table_detected, lri[, c("LIGAND_1", "RECEPTOR_1", "coexpression")],
                    by = c("LIGAND_1", "RECEPTOR_1"),
                    all.x = TRUE)
  return(res_coexp)
}

#' plot_histogram_weights_gene_subset_vs_rest
#'
#' Builds 2 subgraphs from net defined by the gene subset, then returns
#' histogram with distribution of edge weights in the 2 subgraphs.
#'
#' @param net 
#' @param gene_subset - ensembl ids
#' @param lri - for the `IN_SUBSET` category keep only the edges represented in LRI
#' @param df_weights_savefp
#'
#' @return histogram plot
plot_histogram_weights_gene_subset_vs_rest = function(
  net,
  lri_genes,
  lri_coexpressions,
  df_weights_savefp=NULL
) {
  net_nonsubset = igraph::subgraph(
    net, 
    V(net)$name[!(V(net)$name %in% lri_genes)]
  )
  df_weights = rbind(
    data.frame(weights = lri_coexpressions, type = "IN_LRI"),
    data.frame(weights = E(net_nonsubset)$weight, type = "OUT_LRI"),
  )
  if (!is.null(df_weights_savefp)) {
    write.csv(df_weights, df_weights_savefp)
  }
  print_(
    "Num edges subset: ", length(E(net_subset)),
    " | Num edges non-subset: ", length(E(net_nonsubset))
  )
  g_lri = ggplot(data=df_weights[df_weights$type=="IN_LRI",], aes(x=weights)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity')
  g_nonlri = ggplot(data=df_weights[df_weights$type=="OUT_LRI",], aes(x=weights)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity')
  
  return(list(hist_lri=g_lri, hist_nonlri=g_nonlri))
}

#' plot_scdiffcom_results_by_coexp
#'
#' @param table_detected 
#' @param lri 
#' @param convert_na if TRUE, NA -> 0
#'
#' @return boxplot with coexpression by `REGULATION`
plot_scdiffcom_results_by_coexp = function(
  table_detected,
  lri,  # with coexp
  convert_na = TRUE
) {
  res_coexp = add_scdiffcom_results_coexp(table_detected, lri)
  if (convert_na) {
    res_coexp$coexpression[is.na(res_coexp$coexpression)] = 0
  }
  g = ggplot(res_coexp, aes(x=REGULATION, y=coexpression)) + 
    geom_boxplot()
  return(g)
}

# df_degrees = compute_degree_df(net)
# df_degrees[["subset"]] = df_degrees$gene %in% subset
# 
# ggplot(data=df_degrees, aes(x=get("degree"), fill=get("subset"))) +
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#   scale_fill_manual(values=c("#69b3a2", "#404080")) +
#   # theme_ipsum() +
#   labs(fill="")
# # 
# ggplot(data=df_degrees[df_degrees$subset == TRUE,], aes(x=degree, fill=subset)) +
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#   scale_fill_manual(values=c("#69b3a2", "#404080")) +
#   # theme_ipsum() +
#   labs(fill="")


# net_s = igraph::subgraph(net, subset)

## Weighted sampling
# TODO
# LRI sampling based on co-expression:
#   - sample Ligand | Receptor ~ centrality
#   - sample the pair randomly | based on closeness to first gene
# sample(x, n, prob=())



