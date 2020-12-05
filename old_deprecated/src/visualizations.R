# Heatmap plots ------------------------------------------------------------------------

build_heatmaps <- function(dt, dir_path) {
  
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
    
    print(colnames(dt_edge))
    
    nrows = dim(dt_edge)[1]
    if (nrows == 1) {message(paste0("One edge for: ", tissue))}
    
    G = graph_from_data_frame(d = dt_edge,
                              directed = TRUE,
                              vertices = NULL)
    
    matrix_all = as_adjacency_matrix(G, attr = "OR_DIFF", sparse = FALSE)
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


get_celltypes_enrichment <- function(dt, use_adjpval=NULL) {
  
  message('get_celltype_enrichment: upgrade warning: OR -> OR_DIFF')
  
  COL_LIGAND_RECEPTOR_CELLTYPES = "LR_CELLTYPE"
  # COL_pval = 'pval'
  COL_pval = 'pval_DIFF'
  COL_pval_UP = 'pval_UP'
  COL_pval_DOWN = 'pval_DOWN'
  
  # These indicate p-values that have been readjusted
  #  at celltype level for keeping the signal.
  COL_pval_readj = 'pval_readj_on_celltypes'
  COL_pval_readj_UP = 'pval_readj_on_celltypes_UP'
  COL_pval_readj_DOWN = 'pval_readj_on_celltypes_DOWN'
  
  dt_ctypes = dt[
    Category == COL_LIGAND_RECEPTOR_CELLTYPES
    ][,
      c(COL_pval_readj, COL_pval_readj_UP, COL_pval_readj_DOWN) := .(
        p.adjust(get(COL_pval), method='BH'),
        p.adjust(get(COL_pval_UP), method='BH'),
        p.adjust(get(COL_pval_DOWN), method='BH')
      )
      ]
  
  dt_ctypes = dt_ctypes[, c("Ligand_cell", "Receptor_cell") := list(
    sub("_.*", "", Value),  # substitute w/ regex
    sub(".*_", "", Value)
  )
  ]
  
  cols_to_select = c(
    "Tissue", "Ligand_cell", "Receptor_cell", "OR_DIFF", "OR_UP", "OR_DOWN",
    COL_pval, COL_pval_UP, COL_pval_DOWN,
    COL_pval_readj, COL_pval_readj_UP, COL_pval_readj_DOWN
  )
  
  dt_ctypes = dt_ctypes[, ..cols_to_select]
  
  return(dt_ctypes)
}

# Bipartite graph ----

analyze_Graph <- function(
  dt_ora, 
  dt_filtered,
  tissue, 
  config=GRAPH_CONFIG, 
  use_adjpval=TRUE,
  disperse=TRUE,
  dir=NULL,
  analysis_name=NULL) {
  
  ora_has_tissue = 'Tissue' %in% names(dt_ora)
  filt_has_tissue = 'Tissue' %in% names(dt_filtered)
  
  if( !(ora_has_tissue & filt_has_tissue) ) {
    stop(paste0('analyze_Graph: Filtered and ORA must have a tissue specified.',
                ' Insert a dummy tissue for compatibility with current code.'))
  }
  
  message('Solve the low statistical power issue by controlling num interactions
          in the BH adjustment.')
  
  if( is.null(dt_filtered) ) {stop('analyze_Graph: dt_filtered is NULL.')}
  
  # Process ORA results and construct graph
  G = construct_graph(dt_ora, tissue, dt_filtered=dt_filtered)
  
  # Setup
  G = setup_graph(G, use_adjpval=use_adjpval, disperse=disperse)
  
  # Save as edge table and plot
  if ( !is.null(dir) ) {
    
    if ( is.null(analysis_name) ) {
      stop('analyze_Graph: Analysis name not provided.')
    }
    
    # save_run_metadata()
    
    subdirs = c('edge_tables', 'plots')
    create_analysis_dirs(dir, analysis_name, subdirs)
    
    write_as_edge_table(
      G, 
      path = file.path(dir, analysis_name, 'edge_tables', paste0(tissue, '.tsv'))
    )
    
    plot_graph(
      G,
      config, 
      path = file.path(dir, analysis_name, 'plots', paste0(tissue, '.pdf'))
    )
  } else {
    
    if ( !is.null(analysis_name) ) {
      stop('analyze_Graph: Analysis name not null.')
    }
    
    plot_graph(
      G, 
      config, 
      path=NULL)
  }
  
}

create_analysis_dirs <- function(dir,
                                 analysis_name,
                                 subdirs) {
  
  analysis_dir = file.path(dir, analysis_name)
  dir.create(analysis_dir)
  
  map(
    subdirs,
    ~ dir.create( file.path(analysis_dir, .x) )
  )
  
}

color_edge <- function(edge, 
                       config=GRAPH_CONFIG,
                       use_adjpval=TRUE) {
  
  COLOR_UP_AND_DOWN = config$EDGE_COLORING$COLOR_BOTH
  COLOR_UP = config$EDGE_COLORING$COLOR_UP
  COLOR_DOWN = config$EDGE_COLORING$COLOR_DOWN
  COLOR_ROBUST = config$EDGE_COLORING$COLOR_ROBUST
  COLOR_NONE = config$EDGE_COLORING$COLOR_NONE
  
  CUTOFF_CHANGE_UP = config$EDGE_COLORING$CUTOFF_UP
  CUTOFF_CHANGE_DOWN = config$EDGE_COLORING$CUTOFF_DOWN
  CUTOFF_ROBUST_UP = 0.99
  CUTOFF_ROBUST_DOWN = 0.99
  message('color_edge: some hard-coded cutoff.')
  
  message(paste0(
    'color_edge: refactor to set edge state in the ',
    'edge data table to be reused within further analysis.')
  )
  
  if (use_adjpval) {
    is_pval_UP = edge$pval_adj_UP < 0.05
    is_pval_DOWN = edge$pval_adj_DOWN < 0.05
  } else {
    is_pval_UP = edge$pval_UP < 0.05
    is_pval_DOWN = edge$pval_DOWN < 0.05
  }
  is_OR_UP_CHANGE = edge$OR_UP >= CUTOFF_CHANGE_UP
  is_OR_DOWN_CHANGE = edge$OR_DOWN >= CUTOFF_CHANGE_DOWN
  
  is_OR_UP_ROBUST = edge$OR_UP <= CUTOFF_ROBUST_UP
  is_OR_DOWN_ROBUST = edge$OR_DOWN <= CUTOFF_ROBUST_DOWN
  
  if (is_OR_UP_CHANGE & is_OR_DOWN_CHANGE & is_pval_UP & is_pval_DOWN) {
    return(COLOR_UP_AND_DOWN)
  } else if (is_OR_UP_CHANGE & is_pval_UP) {
    return(COLOR_UP)
  } else if (is_OR_DOWN_CHANGE & is_pval_DOWN) {
    return(COLOR_DOWN)
  } else if (is_OR_UP_ROBUST & is_OR_DOWN_ROBUST & is_pval_UP & is_pval_DOWN) {
    return(COLOR_ROBUST)
  } else {
    return(COLOR_NONE)
  }
  
  # if ((edge$OR_UP >= CUTOFF_UP) & (edge$OR_DOWN >= CUTOFF_DOWN)) {
  #   return(COLOR_UP_AND_DOWN)
  # } else if (edge$OR_UP >= CUTOFF_UP) {
  #   return(COLOR_UP)
  # } else if (edge$OR_DOWN >= CUTOFF_DOWN) {
  #   return(COLOR_DOWN)
  # } else {
  #   return(COLOR_NONE)
  # }
}

weight_edge <- function(edge, 
                        min_width, 
                        max_width, 
                        min_value,
                        max_value, 
                        config=GRAPH_CONFIG) {
  
  # TODO: weight_edge and color edge can be combined in process_edge
  # TODO: scale edge width by OR magnitude?
  # (c + min_width)*(max_width/min_width)
  
  CUTOFF_UP = config$EDGE_COLORING$CUTOFF_UP
  CUTOFF_DOWN = config$EDGE_COLORING$CUTOFF_DOWN
  
  # CUTOFF_CHANGE_UP = config$EDGE_COLORING$CUTOFF_UP
  # CUTOFF_CHANGE_DOWN = config$EDGE_COLORING$CUTOFF_DOWN
  # CUTOFF_ROBUST_UP = 0.99
  # CUTOFF_ROBUST_DOWN = 0.99
  # message('weight_edge: some hard-coded cutoff.')
  # 
  # is_pval_UP = edge$pval_adj_UP < 0.05
  # is_pval_DOWN = edge$pval_adj_DOWN < 0.05
  # 
  # is_OR_UP_CHANGE = edge$OR_UP >= CUTOFF_CHANGE_UP
  # is_OR_DOWN_CHANGE = edge$OR_DOWN >= CUTOFF_CHANGE_DOWN
  # 
  # is_OR_UP_ROBUST = edge$OR_UP <= CUTOFF_ROBUST_UP
  # is_OR_DOWN_ROBUST = edge$OR_DOWN <= CUTOFF_ROBUST_DOWN
  
  
  if ((edge$OR_UP >= CUTOFF_UP) & (edge$OR_DOWN >= CUTOFF_DOWN)) {
    return(edge$OR)
  } else if (edge$OR_UP >= CUTOFF_UP) {
    return(edge$OR_UP)
  } else if (edge$OR_DOWN >= CUTOFF_DOWN) {
    return(edge$OR_DOWN)
  } else {
    return(1)
  }
}

construct_graph <- function(dt_ora, tissue, dt_filtered=NULL) {
  # Extract celltypes overrepresentation
  dt_ctypes = get_celltypes_enrichment(dt_ora, cols)
  
  # Subset the tissue
  if( tissue == 'All' ) {
    stop('construct_graph: Not implemented for All. Gets in conflict with some
         downstream merging from filtered datasets that don\'t represent
         explicitely a tissue named "All"')
  }
  dt_edge = dt_ctypes[Tissue == tissue, .SD, .SDcols = !c('Tissue')]
  
  message('construct_graph: colnames are hard-coded in get_celltypes_enrichment')
  # Label the transmitter (L) and receiver (R) cells
  dt_edge = dt_edge[, .(
    'Ligand_cell' = paste0(Ligand_cell, ' (L)'),
    'Receptor_cell' = paste0(Receptor_cell, ' (R)'),
    # 'OR' = OR,
    'OR_DIFF' = OR_DIFF,
    'OR_UP' = OR_UP,
    'OR_DOWN' = OR_DOWN,
    # 'pval' = pval,
    'pval_DIFF' = pval_DIFF,
    'pval_UP' = pval_UP,
    'pval_DOWN' = pval_DOWN,
    'pval_adj' = pval_readj_on_celltypes,
    'pval_adj_UP' =  pval_readj_on_celltypes_UP,
    'pval_adj_DOWN' = pval_readj_on_celltypes_DOWN
  )]
  
  # Handle NAs
  message(paste0('process_edge_dt: #NAs in ORA:celltypes for ',
                 tissue, ':', sum(is.na(dt_edge)),'; NA->1'))
  dt_edge[is.na(dt_edge)] = 1
  
  # Handle Inf
  # Shouldn't have Inf with pseudocounts
  # CLIP_VAL = 50
  if (sum(dt_edge==Inf) > 0) {
    stop('construct_graph: Inf values in dt_edge')
    # message(paste0('process_edge_dt: Inf value present in edge data.table for',
    #                tissue, '; Clipping to ', CLIP_VAL))
    # dt_edge[dt_edge == Inf] = min(dt_edge[, list(OR, OR_UP, OR_DOWN)], CLIP_VAL)
  }
  
  G = graph_from_data_frame(d = dt_edge,
                            directed = TRUE,
                            vertices = NULL)
  
  # Set graph tissue attribute
  G$tissue = tissue
  G$ora = dt_ora
  G$dt = dt_filtered
  G$interacts_counts = count_interactions_cellpairs_tissue(dt_filtered, tissue)
  G$celltype_counts = count_celltypes_tissue(dt_filtered, tissue)
  
  message('construct_graph: igraph G holds as graph attributes dt_ora and dt_filtered.')
  return(G)
}

infer_vertex_types <- function(vertex_names) {
  print(vertex_names)
  vertex_types = map_lgl(
    strsplit(vertex_names, '[()]'),
    ~ .x[2] == 'L'
  )
  print(vertex_types)
  return(vertex_types)
}

setup_graph <- function(G, config=GRAPH_CONFIG, use_adjpval, disperse) {
  G = setup_vertices(G, config)
  G = setup_edges(G, config, use_adjpval)
  G = setup_layout(G, config, disperse)
  return(G)
}

# To be moved
count_celltypes_tissue <- function(dt, tissue) {
  # @dt - dt_filtered
  
  # dt_t = dt[TISSUE == tissue]
  dt_t = copy(dt)
  
  counts_dt = funion(
    dt_t[, .(
      name = L_CELLTYPE,
      counts = (L_NCELLS_OLD + L_NCELLS_YOUNG)/2  # average counts in young-old
    )],
    dt_t[, .(
      name = R_CELLTYPE,
      counts = (R_NCELLS_OLD + R_NCELLS_YOUNG)/2  # average counts in young-old
    )],
  )
  counts_dt = counts_dt[order(counts_dt$name)]
  return(counts_dt)
}

add_vertex_size <- function(G) {
  
  V(G)$vertex.size = 30
  MAXSIZE = 20
  MINSIZE = 5
  
  vertex_names = V(G)$name
  
  # counts_dt = count_celltypes_tissue(
  #   G$dt,
  #   G$tissue
  # )
  
  celltypes_to_vertexnames <- function(counts_dt) {
    
    vertex_dt = funion(
      counts_dt[, .(
        name = paste0(name, ' (L)'),
        counts = counts
      )],
      counts_dt[, .(
        name = paste0(name, ' (R)'),
        counts = counts
      )]
    )
    
    return(vertex_dt)
  }
  
  # vertex_counts_dt = celltypes_to_vertexnames(counts_dt)
  vertex_counts_dt = celltypes_to_vertexnames(G$celltype_counts)
  vertex_counts_dt = vertex_counts_dt[vertex_counts_dt$name %in% vertex_names]
  
  # Order data.table by vertex_names order and extract counts in correct order.
  ord = match(vertex_names, vertex_counts_dt$name)
  vertex_counts = vertex_counts_dt[ord, counts]
  
  V(G)$vertex.counts = vertex_counts
  
  library('scales')
  vertex_sizes = rescale(
    V(G)$vertex.counts, 
    to = c(MINSIZE, MAXSIZE)
  )
  
  V(G)$vertex.size = vertex_sizes
  
  return(G)
}


count_interactions_cellpairs_tissue <- function(dt_filtered, tissue) {
  return(
    dt_filtered[
      # TISSUE == tissue,
      ,
      .(count = .N),
      by = .(L_CELLTYPE, R_CELLTYPE)
      ][, .(
        'Ligand_celltype' = L_CELLTYPE,
        'Receptor_celltype' = R_CELLTYPE,
        'num_interacts' = count
      )]
  )
}


add_edge_width <- function(G) {
  
  WIDTH_MIN = 0
  WIDTH_MAX = 10
  
  # or_max = mapply(max, E(G)$OR_UP, E(G)$OR_DOWN)
  
  edge_dt = data.table(
    as_data_frame(G)[, c('from', 'to')]
  )
  names(edge_dt) = c('Ligand_celltype', 'Receptor_celltype')
  
  edge_dt[, Ligand_celltype := 
            sub('(*) [(]L[)]', '', Ligand_celltype)
          ][,
            Receptor_celltype := 
              sub('(*) [(]R[)]', '', Receptor_celltype)
            ]
  
  edge_dt = merge(
    edge_dt, G$interacts_counts, 
    by = c('Ligand_celltype', 'Receptor_celltype'),
    all.x = TRUE,
    sort = FALSE)
  
  edge_dt[is.na(edge_dt)] = 0
  message('add_edge_width: Some NAs occur, unclear reason. For now NA->0.')
  
  # Correct order is accomplished by the merge.
  num_interacts = edge_dt[, num_interacts]
  
  E(G)$width = rescale(
    sqrt(num_interacts),
    to = c(WIDTH_MIN, WIDTH_MAX)
  )
  message('add_edge_width: widths are scaled non-linearly.')
  
  return(G)
}

setup_vertices <- function(G, config) {
  V(G)$vertex_types = infer_vertex_types(V(G)$name)
  G = add_vertex_size(G)
  return(G)
}

setup_edges <- function(G, config, use_adjpval) {
  E(G)$arrow.size = config$EDGE_STYLE$ARROW_SIZE
  E(G)$edge.color <- map_chr(E(G), ~ color_edge(.x, use_adjpval=use_adjpval))
  G = add_edge_width(G)
  return(G)
}

setup_layout <- function(G, config, disperse) {
  layout = layout_as_bipartite(
    G, 
    types = V(G)$vertex_types, 
    hgap = config$LAYOUT$HGAP,
    vgap = config$LAYOUT$VGAP
  )
  layout = layout[, 2:1]  # horizontal to vertical
  layout = sort_layout_vertices(layout, G, config, disperse)
  G$layout = layout
  return(G)
}

#### Functions for layout ####
get_edges_from_vertex = function(v_name, G) {
  # return(E(G)[from(v_name)])
  return(incident(G, v_name, mode='out'))
}

get_edges_to_vertex = function(v_name, G) {
  return(incident(G, v_name, mode='in'))
}

count_up_and_down_edges = function(edges, config) {
  color_up = config$EDGE_COLORING$COLOR_UP
  color_down = config$EDGE_COLORING$COLOR_DOWN
  
  counts = table(edges$edge.color)
  num_up = as.integer(counts[color_up])
  num_down = as.integer(counts[color_down])
  
  num_up = ifelse(is.na(num_up), 0, num_up)
  num_down = ifelse(is.na(num_down), 0, num_down)
  
  return(list(N_UP=num_up, N_DOWN=num_down))
}

vertex_sort_key_atomic = function(v_name, from, G, config) {
  
  if(from) {
    edges = get_edges_from_vertex(v_name, G)
  }
  else {
    edges = get_edges_to_vertex(v_name, G)
  }
  counts = count_up_and_down_edges(edges, config)
  sort_key = as.integer(counts$N_UP - counts$N_DOWN)
  return(sort_key)
}

vertex_sort_keys = function(v_names, from, G, config) {
  sort_keys = map_dbl(
    v_names,
    ~ vertex_sort_key_atomic(.x, from, G, config)
  )
  return(sort_keys)
}


# layout modifiers
sort_layout_vertices = function(layout, G, config, disperse) {
  
  vertex_names = V(G)$name
  num_v = length(vertex_names)
  midpoint = num_v/2
  ligand_vertices = vertex_names[1:midpoint]
  receptor_vertices = vertex_names[(midpoint+1):num_v]
  
  ligand_keys = vertex_sort_keys(ligand_vertices, from=TRUE, G, config)
  receptor_keys = vertex_sort_keys(receptor_vertices, from=FALSE, G, config)
  
  # TODO: refactor into clear code
  # The reorders layout rows to achieve sorting
  layout_copy = copy(layout)
  new_layout = copy(layout)
  new_layout[order(ligand_keys), ] = layout_copy[1:midpoint, ]
  new_layout[midpoint + order(receptor_keys), ] = layout_copy[(midpoint+1):num_v, ]
  
  if(disperse) {
    ## Make as function
    # Get hgap
    vgap = config$LAYOUT$HGAP
    scale_factor = 3
    
    new_layout[1:midpoint, 2][ligand_keys==0] = new_layout[1:9, 2][ligand_keys==0] + scale_factor*vgap
    new_layout[1:midpoint, 2][ligand_keys>0] = new_layout[1:9, 2][ligand_keys>0] + 2*scale_factor*vgap
    
    new_layout[(midpoint+1):num_v, 2][receptor_keys==0] = new_layout[(midpoint+1):num_v, 2][receptor_keys==0] + scale_factor*vgap
    new_layout[(midpoint+1):num_v, 2][receptor_keys>0] = new_layout[(midpoint+1):num_v, 2][receptor_keys>0] + 2*scale_factor*vgap
  }
  # new_layout = disperse_layout_vertices(new_layout, G, config)
  return(new_layout)
}
#### Write and plot functions ####
write_as_edge_table <- function(G, path) {
  
  if (!is.null(path)) {
    df_edge = as_data_frame(G, what='edges')
    write.table(
      df_edge, 
      file = path,
      sep='\t',
      quote = FALSE,
      row.names = FALSE
    )
  } else {
    stop("write_as_edge_table: path argument not provided.")
  }
  
}

plot_graph <- function(G,
                       config=GRAPH_CONFIG,
                       path=NULL,
                       show_legend=FALSE) {
  
  # Plot
  MAIN = paste0('Communication profile with age: ', G$tissue)
  MARGIN = 0.5
  SUBTITLE = ''
  # SUBTITLE = paste0('Difference graph of overrepresented celltype ',
  #                   'communication that get\n altered with ageing')

  if( !is.null(path) ) {pdf(path)}
  
  plot.igraph(
    G, 
    layout = G$layout,
    vertex.color = config$VERTEX_STYLE$COLOR,
    vertex.size = V(G)$vertex.size,
    vertex.label.dist = config$VERTEX_STYLE$LABEL_DIST,
    vertex.label.cex = config$VERTEX_STYLE$LABEL_CEX,
    edge.color = E(G)$edge.color,
    edge.width = E(G)$edge.width,
    main = MAIN,
    margin = MARGIN,
    sub = SUBTITLE
  )
  
  if (show_legend) {
    LEGEND_TITLE = 'Edge color legend'
    LEGEND_COLORS = c(
      config$EDGE_COLORING$COLOR_NONE, 
      config$EDGE_COLORING$COLOR_UP,
      config$EDGE_COLORING$COLOR_DOWN,
      config$EDGE_COLORING$COLOR_BOTH
    )
    
    legend(x=-1.5, y=-1.1, 
           
           legend = config$LEGEND$LEGEND_LABELS,
           col = LEGEND_COLORS,
           title = LEGEND_TITLE,
           
           bty = 'o',
           pch = config$LEGEND$PCH,
           cex = config$LEGEND$CEX,
           pt.cex = config$LEGEND$PT.CEX,
           bg = config$LEGEND$BG,
           ncol = config$LEGEND$NCOL
    )
  }
  if( !is.null(path) ) {dev.off()}
}
