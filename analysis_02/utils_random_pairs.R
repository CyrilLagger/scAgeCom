#' run_internal_analysis
#'
#' Runs an internal scDiffCom analysis on the `seruat_obj` using `lri_table` for pairing ligands and receptors.
#'
#' @param seurat_obj 
#' @param lri_table 
#'
#' @return 
run_internal_analysis <- function(
  seurat_obj,
  lri_table
) {
  params = list(
    LRI_species = "mouse",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = list(
      column_name = "age_group",
      cond1_name = "YOUNG",
      cond2_name = "OLD"
    ),  
    object_name = "scDiffCom_object",
    seurat_assay = "RNA",
    seurat_slot = "data",
    log_scale = FALSE,
    score_type = "geometric_mean",
    threshold_min_cells = 5,
    threshold_pct = 0.1,
    iterations = 1000,
    threshold_quantile_score = 0.2,
    threshold_p_value_specificity = 0.05,
    threshold_p_value_de = 0.05,
    threshold_logfc = log(1.5),
    return_distributions = FALSE,
    seed = 42,
    verbose = TRUE
  )
  
  res = scDiffCom:::run_internal_raw_analysis(
    seurat_object = seurat_obj,
    LRI_table = lri_table,
    LRI_species = "mouse",
    params = params
  )
  
  return(res)
}

#' get_genes_from_seurat
#'
#' @param seurat_obj 
#'
#' @return list of genes in `seurat_obj`
get_genes_from_seurat = function(seurat_obj) {
  return(seurat_obj@assays$RNA@counts@Dimnames[[1]])
}

#' seurat_shuffle
#'
#' @param seurat_obj 
#'
#' @return a seurat object with shuffled genes
seurat_shuffle <- function(
  seurat_obj
) {
  
  # From https://github.com/satijalab/seurat/issues/1049
  RenameGenesSeurat <- function(obj, newnames) { 
    # Replace gene names in different slots of a Seurat object. 
    # Run this before integration. 
    # It only changes obj@assays$RNA@counts, @data and @scale.data.
    print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
    RNA <- obj@assays$RNA
    
    if (nrow(RNA) == length(newnames)) {
      if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
      if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
      if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
    } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
    obj@assays$RNA <- RNA
    return(obj)
  }
  
  shuffled_names = sample(get_genes_from_seurat(seurat_obj))
  shuffled_obj = copy(seurat_obj)
  shuffled_obj = RenameGenesSeurat(shuffled_obj, newnames = shuffled_names)
  return(shuffled_obj)
}

#' subset_simple_lri
#'
#' @param lri 
#'
#' @return lri subset with simple LRI
subset_simple_lri <- function(lri) {
  return(
    lri[(is.na(lri$LIGAND_2)) & (is.na(lri$RECEPTOR_2)) & (is.na(lri$RECEPTOR_3))]
  )
}

#' get_genes_lri
#'
#' Extracts lists with different categories of genes (defined in `gene_type`) from an `lri` dataset.
#'
#' @param lri 
#' @param gene_type in {LR, L, R, LR_simple, L_simple, R_simple}
#'
#' @return list of genes with specified `gene_type`
get_genes_lri <- function(
  lri,
  gene_type = "LR"
) {
  
  if (gene_type == "LR") {
    return(
      unique(c(
        lri$LIGAND_1, lri$LIGAND_2, lri$RECEPTOR_1, lri$RECEPTOR_2, lri$RECEPTOR_3
      ))
    )    
  } else if (gene_type == "L") {
    return(
      unique(c(
        lri$LIGAND_1, lri$LIGAND_2
      ))
    )
  } else if (gene_type == "R") {
    return(
      unique(c(
        lri$RECEPTOR_1, lri$RECEPTOR_2, lri$RECEPTOR_3
      ))
    )
  } else {
    lri_simple = subset_simple_lri(lri)
    if (gene_type == "LR_simple") {
      return(unique(c(
        lri_simple$LIGAND_1, lri_simple$LIGAND_2, lri_simple$RECEPTOR_1, lri_simple$RECEPTOR_2, lri_simple$RECEPTOR_3
      )))
    } else if (gene_type == "L_simple") {
      return(unique(c(
        lri_simple$LIGAND_1, lri_simple$LIGAND_2
      )))
    } else if (gene_type == "R_simple") {
      return(unique(c(
        lri_simple$RECEPTOR_1, lri_simple$RECEPTOR_2, lri_simple$RECEPTOR_3
      )))
    }
  }
}

#' count_lri_gene_types
#'
#' @param lri 
#'
#' @return list with counts of various `gene_type` (see `get_genes_lri`)
count_lri_gene_types <- function(lri) {
  return(list(
    LR = length(get_genes_lri(lri, "LR")),
    L = length(get_genes_lri(lri, "L")),
    R = length(get_genes_lri(lri, "R")),
    LR_simple = length(get_genes_lri(lri, "LR_simple")),
    L_simple = length(get_genes_lri(lri, "L_simple")),
    R_simple = length(get_genes_lri(lri, "R_simple"))
  ))
}

#' create_random_LRIs
#'
#' A general function for generating different kinds of LRI datasets.
#'
#' @param seurat_obj 
#' @param lri 
#' @param random_lri_nrow 
#' @param build_pairs_from_same_genes 
#' @param build_pairs_from_same_genes_type 
#'
#' @return
create_random_LRIs <- function(
  seurat_obj,
  lri,
  random_lri_nrow = NULL,
  build_pairs_from_same_genes = FALSE,
  build_pairs_from_same_genes_type = NULL
) {
  
  seurat_g = get_genes_from_seurat(seurat_obj)
  lri_g = get_genes_lri(lri, "LR")
  
  nonLRI_genes_in_seurat_mask = !(seurat_g %in% lri_g)
  nonLRI_genes_in_seurat = seurat_g[nonLRI_genes_in_seurat_mask]
  # num_nonLRI_genes_in_seurat = sum(nonLRI_genes_in_seurat_mask)
  
  lri_simple_with_seurat_genes = (
    subset_simple_lri(lri)[ (LIGAND_1 %in% seurat_g) & (RECEPTOR_1 %in% seurat_g) ]
  )
  num_lri_simple_with_seurat_genes = nrow(lri_simple_with_seurat_genes)
  
  if (!is.null(random_lri_nrow)) {
    N = random_lri_nrow
  } else {
    N = num_lri_simple_with_seurat_genes
  }
  
  if (build_pairs_from_same_genes) {
    if (is.null(build_pairs_from_same_genes_type)) {
      Ligand_1 = Receptor_1 = sample(seurat_g, size = N, replace = TRUE)
    } else if (build_pairs_from_same_genes_type == "Ligands") {
      L_sample = sample(
        seurat_g[seurat_g %in% get_genes_lri(lri, "L_simple")],
        size = N,
        replace = FALSE
      )
      Ligand_1 = Receptor_1 = L_sample
    } else if (build_pairs_from_same_genes_type == "Receptors") {
      R_sample = sample(
        seurat_g[seurat_g %in% get_genes_lri(lri, "R_simple")],
        size = N,
        replace = FALSE
      )
      Ligand_1 = Receptor_1 = R_sample
    } else if (build_pairs_from_same_genes_type == "Non-LR") {
      Ligand_1 = Receptor_1 = sample(nonLRI_genes_in_seurat, size = N, replace = FALSE)
    } else {stop()}
  } else {
    Ligand_1 = sample(nonLRI_genes_in_seurat, size = N, replace = TRUE)
    Receptor_1 = sample(nonLRI_genes_in_seurat, size = N, replace = TRUE)    
  }
  
  return(data.table(
    LRI = paste(Ligand_1, Receptor_1, sep = ":"),
    LIGAND_1 = Ligand_1,
    LIGAND_2 = NA,
    RECEPTOR_1 = Receptor_1,
    RECEPTOR_2 = NA,
    RECEPTOR_3 = NA
  ))
}

print_ = function(...) {
  print(paste0(...))
}

plot_volcano = function(res) {
  plot(res@cci_table_detected$LOGFC, -log10(res@cci_table_detected$P_VALUE_DE))
}

#' dge_by_celltypes
#'
#' @param seurat_obj 
#'
#' @return dataframe with DGE results by celltype
dge_by_celltypes = function(seurat_obj) {
  
  cell_types = unique(seurat_obj@meta.data$cell_type)
  Idents(object = seurat_obj) = "cell_type"
  
  markers = data.frame()
  for (ct in cell_types) {
    print(ct)
    
    metadata = seurat_obj[[]]
    age_groups = metadata[metadata$cell_type == ct, "age_group"]
    
    if ((sum(age_groups == "YOUNG") < 3)
        | sum(age_groups == "OLD") < 3) {
      next()
    }
    
    markers_ct = Seurat::FindMarkers(
      seurat_obj,
      ident.1 = "OLD",
      group.by = "age_group",
      subset.ident = ct,
      # logfc.threshold = 0.25, # default
      min.pct = 0.1,
      min.cells.feature = 5,
      # base = exp(1),
      # logfc.threshold = 1.5  # too stringent for some tissues
      # test.use = "MAST"
    )
    
    if (nrow(markers_ct) == 0) {
      warning(glue("No markers found for {ct}."))
      next()
    }
    
    markers_ct['gene'] = rownames(markers_ct)
    markers_ct['cell_type'] = ct
    rownames(markers_ct) = NULL
    
    markers = rbind(markers, markers_ct)
  }
  return(data.table::data.table(markers))
}

drop_genes_from_seurat = function(
  seurat_obj,
  genes_to_drop
){
  counts = GetAssayData(seurat_obj, assay = "RNA")
  s = subset(
    seurat_obj, 
    features = rownames(counts)[
      which(!(rownames(counts) %in% genes_to_drop))
    ]
  )
}

#' simulate_random_gene_drop
#'
#' @param seurat_obj 
#' @param genes_to_drop - list of genes from which to drop in every simulation run
#' @param num_sim - number of runs in simulation
#' @param lri - lri dataset to be used
#'
#' @return dataframe with simulation results (1 row per run)
simulate_random_gene_drop = function(
  seurat_obj,
  genes_to_drop,
  num_sim,
  lri = lri_simple_in_seurat,
  with_shuffling = TRUE
) {
  len_genes_to_drop = length(genes_to_drop)
  nums_to_drop = sample(1:len_genes_to_drop, num_sim, replace = TRUE)
  t = c()
  for (i in 1:num_sim) {
    
    print_("i = ", i)
    num_to_drop = nums_to_drop[i]
    
    sample_genes_to_drop = sample(genes_to_drop, num_to_drop)
    seurat_obj_with_dropped_genes = drop_genes_from_seurat(seurat_sample_tms_liver, sample_genes_to_drop)
    
    if (with_shuffling) {
      seurat_obj_with_dropped_genes = seurat_shuffle(
        seurat_obj_with_dropped_genes
      )
    }
    
    res_shuffled = run_internal_analysis(
      seurat_obj = seurat_obj_with_dropped_genes,
      lri_table = lri
    )
    res_shuffled = FilterCCI(res_shuffled, skip_ora = TRUE)
    
    t_i = table(res_shuffled@cci_table_detected$REGULATION)
    t_i = c(t_i, num_to_drop)
    
    t = rbind(t, t_i)
  }
  return(t)
}

#' read_drop_simulation_results
#'
#' @param which in {dge_drop, stable_drop, dge_drop_no_shuffle, stable_drop_no_shuffle}
#'
#' @return df with saved simulation results for each type `which`
read_drop_simulation_results = function(
  which = "dge_drop"
) {
  base_path = "data/"
  dge_drop_filepath = "regulation_vs_dge_drop.csv"
  stable_drop_filepath = "regulation_vs_dge_stable_drop.csv"
  
  if (which == "dge_drop") {
    filepath = dge_drop_filepath
  } else if (which == "stable_drop") {
    filepath = stable_drop_filepath
  } else if (which == "dge_drop_no_shuffle") {
    filepath = "regulation_vs_dge_drop_no_shuffle.csv"
  } else if (which == "stable_drop_no_shuffle") {
    filepath = "regulation_vs_dge_stable_drop_no_shuffle.csv"
  }
  
  df = read.csv(paste0(base_path, filepath), row.names = NULL)
  names(df) = c(names(df)[1:length(names(df))-1], "DROP")
  if ( !("UP" %in% names(df)) ) {
    df$UP = 0
  }
  df$SC = df$UP + df$DOWN
  return(df)
}

#' map_ensembl_ids
#'
#' Uses biomart to map gene names to ensembl IDs.
#'
#' @param gene_names 
#' @param ensembl 
#'
#' @return df: names - ensembl ids
map_ensembl_ids = function(
  gene_names,
  ensembl = NULL
) {
  if (is.null(ensembl)) {
    ensembl = biomaRt::useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version=95)
  }
  gene_names_to_ids = biomaRt::getBM(attributes=c("external_gene_name", 'ensembl_gene_id'), 
                                     filters = 'external_gene_name',
                                     values = gene_names,
                                     mart = ensembl)
  return(gene_names_to_ids)
}


#' add_lri_ensembl_ids
#'
#' Adds ensembl ids to lri datasets consisting of simple LRIs.
#' Not implemented for complex LRIs - raises Exception.
#'
#' @param lri 
#'
#' @return lri with new columns: LIGAND_1_ID, RECEPTOR_1_ID
add_lri_ensembl_ids = function(
  lri
) {
  if (any(!is.na(lri$LIGAND_2)) | (any(!is.na(lri$RECEPTOR_2))) | (any(!is.na(lri$RECEPTOR_2)))) {
    stop("Not implemented for complex interactions")
  }
  
  ligand_1_ids = map_ensembl_ids(lri$LIGAND_1)
  ligand_1_ids["LIGAND_1_ID"] = ligand_1_ids["ensembl_gene_id"]
  ligand_1_num_dupl = sum(duplicated(ligand_1_ids$external_gene_name))
  ligand_1_ids = ligand_1_ids[!duplicated(ligand_1_ids$external_gene_name), ]
  
  receptor_1_ids = map_ensembl_ids(lri$RECEPTOR_1)
  receptor_1_ids["RECEPTOR_1_ID"] = receptor_1_ids["ensembl_gene_id"]
  receptor_1_num_dupl = sum(duplicated(receptor_1_ids$external_gene_name))
  receptor_1_ids = receptor_1_ids[!duplicated(receptor_1_ids$external_gene_name), ]
  
  warning(paste0("Removed duplicates: Ligand 1: ", ligand_1_num_dupl, " | Receptor 1: ", receptor_1_num_dupl))
  
  lri = merge(lri, ligand_1_ids, by.x = "LIGAND_1", by.y = "external_gene_name", all.X = TRUE)
  lri = merge(lri, receptor_1_ids, by.x = "RECEPTOR_1", by.y = "external_gene_name", all.X = TRUE)
  return(lri)
}


simulate_random_seurat_regulation = function(
  seurat_sample_obj,
  lri,
  save_path,
  num_samples,
  read_only=FALSE
) {
  if (file.exists(save_path)) {
    t = read.csv(save_path)
    t$X = NULL
    t = data.table(t)
  } else {
    t = c()
  }

  if (!read_only) {
    for (i in 1:num_samples) {
      print(paste0("simulate_random_seurat_regulation: i = ", i))
      res_shuffled = run_internal_analysis(
        seurat_obj = seurat_shuffle(seurat_sample_obj),
        lri_table = lri
      )
      res_shuffled = FilterCCI(res_shuffled, skip_ora = TRUE)
      t = rbind(t, rbind(table(res_shuffled@cci_table_detected$REGULATION)), fill=TRUE)
      t = data.table(t)
      write.csv(t, save_path)
    }
  }
  return(t)
}


get_seurat = function(
  INPUT,
  TISSUE,
  AGE_GROUP,
  SEX
) {
  if (INPUT == "droplet") {
    seurat_sample_obj = readRDS("../data/seurat_shared_tms_droplet.rds")
    seurat_sample_obj = subset(x = seurat_sample_obj, subset = tissue == TISSUE)
    seurat_sample_obj[["age_group"]] = ifelse(
      seurat_sample_obj[[]]$age %in% AGE_GROUP$YOUNG, 
      "YOUNG", 
      ifelse(
        seurat_sample_obj[[]]$age %in% AGE_GROUP$OLD, 
        "OLD", 
        "IGNORE"
      )
    )
    seurat_sample_obj = subset(x = seurat_sample_obj, subset = age_group != "IGNORE")
    seurat_sample_obj = subset(x = seurat_sample_obj, subset = sex %in% SEX)
    seurat_sample_obj[["cell_type"]] = seurat_sample_obj$cell_type_scagecom
    seurat_sample_obj[["cell_abbreviation"]] = seurat_sample_obj$cell_abbreviation_scagecom
  } else if (INPUT == "facs") {
    seurat_sample_obj = readRDS("../data/seurat_shared_tms_facs.rds")
    seurat_sample_obj = subset(x = seurat_sample_obj, subset = tissue == TISSUE)
    seurat_sample_obj[["age_group"]] = ifelse(
      seurat_sample_obj[[]]$age %in% AGE_GROUP$YOUNG, 
      "YOUNG", 
      ifelse(
        seurat_sample_obj[[]]$age %in% AGE_GROUP$OLD, 
        "OLD", 
        "IGNORE"
      )
    )
    seurat_sample_obj = subset(x = seurat_sample_obj, subset = age_group != "IGNORE")
    seurat_sample_obj = subset(x = seurat_sample_obj, subset = sex %in% SEX)
    seurat_sample_obj[["cell_type"]] = seurat_sample_obj$cell_type_scagecom
    seurat_sample_obj[["cell_abbreviation"]] = seurat_sample_obj$cell_abbreviation_scagecom
  } else if (INPUT == "sample") {
    load("../data/seurat_sample_tms_liver.rda")
    seurat_sample_obj = seurat_sample_tms_liver
  } else {
    stop("INPUT not recognized.")
  }
  return(seurat_sample_obj)
}
