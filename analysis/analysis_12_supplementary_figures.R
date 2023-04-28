####################################################
##
## Project: scAgeCom
##
## lagger.cyril@gmail.com
## ursu_eugen@hotmail.com
## anais.equey@gmail.com
##
## Supplementary Figure preparation
##
####################################################
##

## secretomics association Supplementary Figure ####

## Prepare Figure Upset Plot Validation dataset

setnames(
  val_dt_selection,
  old = c(
    "pancreas", "huvec", "macro", "mscat",
    "neurons", "glia", "cardio"
  ),
  new = c(
    "hPDE (Li 2022)",
    "hUVEC (Zhao 2020)",
    "mBMM (Meissner 2013)",
    "mMSC-AT (Acar 2020)",
    "mNeuron (Tushaus 2020)",
    "mGlial (Tushaus 2020)",
    "rCM (Kuhn 2020)"
  )
)

figp_val_upset <- ComplexUpset::upset(
  as.data.frame(val_dt_selection),
  colnames(val_dt_selection)[-c(1, 7, 9, 10)],
  name = "",
  set_sizes = ComplexUpset::upset_set_size()
  + ylab("# of LR genes"),
  base_annotations = list(
    "Intersection Size" = ComplexUpset::intersection_size(
      counts = TRUE,
      bar_number_threshold = 100,
      text = list(size = 8),
      colour = "coral",
      fill = "coral"
    ) + theme(
      axis.title.y = element_text(
        margin = margin(t = 0, r = -200, b = 0, l = 0)
      )
    )
  ),
  themes = ComplexUpset::upset_default_themes(
    text = element_text(size = 40),
    plot.title = element_text(size = 32)
  ),
  min_size = 8
) + ggtitle(
  "Ligand/Receptor genes matched to secreted proteins"
)
ggsave(
  paste0(
    path_scagecom_output,
    "fig_val_upset2.png"
  ),
  plot = figp_val_upset,
  width = 2100,
  height = 1200,
  units = "px",
  scale = 3
)

## mMSC-AT

sfig_dt_val_mscat_at <- dt_val_mscat_at[, .(mean(OR), exp(mean(log(BH)))), by = "category"]

fig_ora_mmscat_at <- ggplot(
  sfig_dt_val_mscat_at,
  aes(
    x = V1,
    y = reorder(category, V1),
    size = -log10(
      V2 + min(sfig_dt_val_mscat_at[V2 != 0]$V2)
    )
  )
) + geom_point(
) + geom_vline(
  xintercept = 1,
  linetype = "dashed",
  size = 1.5
) + scale_size(
  name = "-log10(Adj. P-Value)",
  breaks = round(fivenum(-log10(
    sfig_dt_val_mscat_at$V2 + min(sfig_dt_val_mscat_at[V2 != 0]$V2)
  ))),
  labels = round(fivenum(-log10(
    sfig_dt_val_mscat_at$V2 + min(sfig_dt_val_mscat_at[V2 != 0]$V2)
  ))),
  limit = c(0, 50)
) + xlab(
  "Odds Ratio"
) + ylab(
  "Emitter Cell Type from TMS Adipose Tissues"
  #) + theme_minimal(
) + ggtitle(
  "Association with the mMSC-AT secretome"
) + theme(
  text = element_text(size = 18),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  axis.text = element_text(size = 19)
)
ggsave(
  paste0(
    path_scagecom_output,
    "fig_ora_mmscat_at2.png"
  ),
  plot = fig_ora_mmscat_at,
  width = 2100,
  height = 1200,
  units = "px",
  scale = 1.6
)

## mNeuron

sfig_dt_val_neuron_brain <- dt_val_neuron_brain[, .(mean(OR), exp(mean(log(BH)))), by = "category"]

fig_ora_neuron_brain <- ggplot(
  sfig_dt_val_neuron_brain,
  aes(
    x = V1,
    y = reorder(category, V1),
    size = -log10(
      V2 + min(sfig_dt_val_neuron_brain[V2 != 0]$V2)
    )
  )
) + geom_point(
) + geom_vline(
  xintercept = 1,
  linetype = "dashed",
  size = 1.5
) + scale_size(
  name = "-log10(Adj. P-Value)",
  breaks = round(fivenum(-log10(
    sfig_dt_val_neuron_brain$V2 + min(sfig_dt_val_neuron_brain[V2 != 0]$V2)
  ))),
  labels = round(fivenum(-log10(
    sfig_dt_val_neuron_brain$V2 + min(sfig_dt_val_neuron_brain[V2 != 0]$V2)
  ))),
  limit = c(0, 20)
) + xlab(
  "Odds Ratio"
) + ylab(
  "Emitter Cell Type from TMS Brain Samples"
  #) + theme_minimal(
) + ggtitle(
  "Association with the mNeuron secretome"
) + theme(
  text = element_text(size = 18),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  axis.text = element_text(size = 19)
)

ggsave(
  paste0(
    path_scagecom_output,
    "fig_ora_neuron_brain2.png"
  ),
  plot = fig_ora_neuron_brain,
  width = 2100,
  height = 1200,
  units = "px",
  scale = 1.6
)

## rCM

sfig_dt_val_cardio_heart <- dt_val_cardio_heart[, .(mean(OR), exp(mean(log(BH)))), by = "category"]

fig_ora_cardio_heart <- ggplot(
  sfig_dt_val_cardio_heart,
  aes(
    x = V1,
    y = reorder(category, V1),
    size = -log10(
      V2 + min(sfig_dt_val_cardio_heart[V2 != 0]$V2)
    )
  )
) + geom_point(
) + geom_vline(
  xintercept = 1,
  linetype = "dashed",
  size = 1.5
) + scale_size(
  name = "-log10(Adj. P-Value)",
  breaks = round(fivenum(-log10(
    sfig_dt_val_cardio_heart$V2 + min(sfig_dt_val_cardio_heart[V2 != 0]$V2)
  ))),
  labels = round(fivenum(-log10(
    sfig_dt_val_cardio_heart$V2 + min(sfig_dt_val_cardio_heart[V2 != 0]$V2)
  ))),
  limit = c(0, 80)
) + xlab(
  "Odds Ratio"
) + ylab(
  "Emitter Cell Type from TMS Heart Samples"
  #) + theme_minimal(
) + ggtitle(
  "Association with the rCM secretome"
) + theme(
  text = element_text(size = 18),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  axis.text = element_text(size = 19)
)

ggsave(
  paste0(
    path_scagecom_output,
    "fig_ora_cardio_heart2.png"
  ),
  plot = fig_ora_cardio_heart,
  width = 2100,
  height = 1200,
  units = "px",
  scale = 1.6
)

## hPDE

sfig_dt_val_pancreas_pct <- dt_val_pancreas_pct[, .(mean(OR), exp(mean(log(BH)))), by = "category"]

fig_ora_hpde_pancreas <- ggplot(
  sfig_dt_val_pancreas_pct,
  aes(
    x = V1,
    y = reorder(category, V1),
    size = -log10(
      V2 + min(sfig_dt_val_pancreas_pct[V2 != 0]$V2)
    )
  )
) + geom_point(
) + geom_vline(
  xintercept = 1,
  linetype = "dashed",
  size = 1.5
) + scale_size(
  name = "-log10(Adj. P-Value)",
  breaks = round(fivenum(-log10(
    sfig_dt_val_pancreas_pct$V2 + min(sfig_dt_val_pancreas_pct[V2 != 0]$V2)
  ))),
  labels = round(fivenum(-log10(
    sfig_dt_val_pancreas_pct$V2 + min(sfig_dt_val_pancreas_pct[V2 != 0]$V2)
  ))),
  limit = c(0, 20)
) + xlab(
  "Odds Ratio"
) + ylab(
  "Emitter Cell Type from TMS Pancreas Samples"
  #) + theme_minimal(
) + ggtitle(
  "Association with the hPDE secretome"
) + theme(
  text = element_text(size = 18),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  axis.text = element_text(size = 19)
)

ggsave(
  paste0(
    path_scagecom_output,
    "fig_ora_hpde_pancreas2.png"
  ),
  plot = fig_ora_hpde_pancreas,
  width = 2100,
  height = 1200,
  units = "px",
  scale = 1.6
)


## scDiffCom vs SDGA Supplementary Figure ####

dt_sdea_fig <- dt_cci_sdea_comp[
  ,
  c("dataset", "scDiffCom_regulation", "ligand_deg", "receptor_deg")
]
dt_sdea_fig[, "CCI Not DE" := scDiffCom_regulation == "NO"]
dt_sdea_fig[, "CCI Up" := scDiffCom_regulation == "UP"]
dt_sdea_fig[, "CCI Down" := scDiffCom_regulation == "DOWN"]
dt_sdea_fig[, "Ligand Not DE" := ligand_deg == "NOT DE"]
dt_sdea_fig[, "Ligand Up" := ligand_deg == "UP"]
dt_sdea_fig[, "Ligand Down" := ligand_deg == "DOWN"]
dt_sdea_fig[, "Receptor Not DE" := receptor_deg == "NOT DE"]
dt_sdea_fig[, "Receptor Up" := receptor_deg == "UP"]
dt_sdea_fig[, "Receptor Down" := receptor_deg == "DOWN"]

data.table::setcolorder(
  dt_sdea_fig,
  neworder = c(
    "dataset", "scDiffCom_regulation", "ligand_deg", "receptor_deg",
    "Receptor Down", "Receptor Up", "Receptor Not DE",
    "Ligand Down", "Ligand Up", "Ligand Not DE",
    "CCI Down", "CCI Up", "CCI Not DE"
  )
)

sfig_sdea_comp <- upset(
  as.data.frame(dt_sdea_fig),
  colnames(dt_sdea_fig)[5:13],
  name = "",
  base_annotations = list(
    "Intersection Size" = intersection_size(
      counts = TRUE,
      text = list(
        vjust = -0.8,
        hjust = 0.2,
        angle = 45
        #size = 4,
      ),
      text_colors = c(on_background = "black")
    )
  ),
  set_sizes = FALSE,
  themes = upset_default_themes(
    text = element_text(size = 30),
    plot.title = element_text(size = 30)
  ),
  min_size = 0,
  max_size = 50000,
  sort_sets = FALSE,
  stripes = upset_stripes(
    colors = c(
      rep("grey95", 3),
      rep("cornsilk1", 3),
      rep("grey95", 3)
    )
  ),
  queries = list(
    upset_query(
        intersect = c(
          "CCI Not DE",
          "Ligand Not DE",
          "Receptor Up"
        ),
        color = "#76b5c5", fill = "#76b5c5"
    ),
    upset_query(
        intersect = c(
          "CCI Not DE",
          "Ligand Not DE",
          "Receptor Down"
        ),
        color = "#76b5c5", fill = "#76b5c5"
    ),
    upset_query(
        intersect = c(
          "CCI Not DE",
          "Ligand Up",
          "Receptor Not DE"
        ),
        color = "#76b5c5", fill = "#76b5c5"
    ),
    upset_query(
        intersect = c(
          "CCI Not DE",
          "Ligand Up",
          "Receptor Up"
        ),
        color = "#e28743", fill = "#e28743"
    ),
    upset_query(
        intersect = c(
          "CCI Not DE",
          "Ligand Up",
          "Receptor Down"
        ),
        color = "#76b5c5", fill = "#76b5c5"
    ),
    upset_query(
        intersect = c(
          "CCI Not DE",
          "Ligand Down",
          "Receptor Not DE"
        ),
        color = "#76b5c5", fill = "#76b5c5"
    ),
    upset_query(
        intersect = c(
          "CCI Not DE",
          "Ligand Down",
          "Receptor Up"
        ),
        color = "#76b5c5", fill = "#76b5c5"
    ),
    upset_query(
        intersect = c(
          "CCI Not DE",
          "Ligand Down",
          "Receptor Down"
        ),
        color = "#e28743", fill = "#e28743"
    ),
    upset_query(
        intersect = c(
          "CCI Up",
          "Ligand Not DE",
          "Receptor Not DE"
        ),
        color = "#eab676", fill = "#eab676"
    ),
    upset_query(
        intersect = c(
          "CCI Up",
          "Ligand Not DE",
          "Receptor Up"
        ),
        color = "#76b5c5", fill = "#76b5c5"
    ),
    upset_query(
        intersect = c(
          "CCI Up",
          "Ligand Not DE",
          "Receptor Down"
        ),
        color = "#873e23", fill = "#873e23"
    ),
    upset_query(
        intersect = c(
          "CCI Up",
          "Ligand Up",
          "Receptor Not DE"
        ),
        color = "#76b5c5", fill = "#76b5c5"
    ),
    upset_query(
        intersect = c(
          "CCI Up",
          "Ligand Up",
          "Receptor Up"
        ),
        color = "#1e81b0", fill = "#1e81b0"
    ),
    upset_query(
        intersect = c(
          "CCI Up",
          "Ligand Up",
          "Receptor Down"
        ),
        color = "#76b5c5", fill = "#76b5c5"
    ),
    upset_query(
        intersect = c(
          "CCI Up",
          "Ligand Down",
          "Receptor Not DE"
        ),
        color = "#873e23", fill = "#873e23"
    ),
    upset_query(
        intersect = c(
          "CCI Up",
          "Ligand Down",
          "Receptor Up"
        ),
        color = "#76b5c5", fill = "#76b5c5"
    ),
    # upset_query(
    #     intersect = c(
    #       "CCI Up",
    #       "Ligand Down",
    #       "Receptor Down"
    #     ),
    #     color = "#873e23", fill = "#873e23"
    # ),
    upset_query(
        intersect = c(
          "CCI Down",
          "Ligand Not DE",
          "Receptor Not DE"
        ),
        color = "#eab676", fill = "#eab676"
    ),
    upset_query(
        intersect = c(
          "CCI Down",
          "Ligand Not DE",
          "Receptor Up"
        ),
        color = "#873e23", fill = "#873e23"
    ),
    upset_query(
        intersect = c(
          "CCI Down",
          "Ligand Not DE",
          "Receptor Down"
        ),
        color = "#76b5c5", fill = "#76b5c5"
    ),
    upset_query(
        intersect = c(
          "CCI Down",
          "Ligand Up",
          "Receptor Not DE"
        ),
        color = "#873e23", fill = "#873e23"
    ),
    # upset_query(
    #     intersect = c(
    #       "CCI Down",
    #       "Ligand Up",
    #       "Receptor Up"
    #     ),
    #     color = "#873e23", fill = "#873e23"
    # ),
    upset_query(
        intersect = c(
          "CCI Down",
          "Ligand Up",
          "Receptor Down"
        ),
        color = "#76b5c5", fill = "#76b5c5"
    ),
    upset_query(
        intersect = c(
          "CCI Down",
          "Ligand Down",
          "Receptor Not DE"
        ),
        color = "#76b5c5", fill = "#76b5c5"
    ),
    upset_query(
        intersect = c(
          "CCI Down",
          "Ligand Down",
          "Receptor Up"
        ),
        color = "#76b5c5", fill = "#76b5c5"
    ),
    upset_query(
        intersect = c(
          "CCI Down",
          "Ligand Down",
          "Receptor Down"
        ),
        color = "#1e81b0", fill = "#1e81b0"
    )
  )
) + ggtitle(
  paste0(
    "Comparison of differential expression results: ",
    "CCI level vs gene level"
  )
) + scale_color_manual(
  name = "Category",
  breaks = c(
    "Consistent", "Resolve ambiguity",
    "Additive effect", "Subtractive effect",
    "Inconsistent"
  ),
  values = c(
    "Consistent" = "#1e81b0",
    "Resolve ambiguity" = "#76b5c5",
    "Additive effect" = "#eab676",
    "Subtractive effect" = "#e28743",
    "Inconsistent" = "#873e23"
  )
) + theme(
  legend.position = c(0.8, 2.8)
)
sfig_sdea_comp

ggsave(
  paste0(
    path_scagecom_output,
    "sfig_sdea_comp2.png"
  ),
  plot = sfig_sdea_comp,
  width = 2100,
  height = 1200,
  units = "px",
  scale = 2
)
#manual save: 2100x1200

## Immune processes Supplementary Figure ####

# prepare GO/KEGG table
sf_immune_gokegg <- unique(
  aging_dt_immune_sex[
    ORA_CATEGORY %in% c("GO_TERMS", "KEGG_PWS", "ER_CELLFAMILIES"),
    c("ORA_CATEGORY", "VALUE", "ORA_REGULATION_age", "Overall (Union)_age")
  ]
)[order(-`Overall (Union)_age`)]
sf_immune_gokegg[
  ,
  ORA_CATEGORY := ifelse(
    ORA_CATEGORY == "GO_TERMS",
    "GO",
    ifelse(
       ORA_CATEGORY == "KEGG_PWS",
       "KEGG",
       "Cell Family"
    )
  )
]
setnames(
  sf_immune_gokegg,
  old = colnames(sf_immune_gokegg),
  new = c(
    "Category",
    "Term",
    "Age Regulation",
    "Over-represented in # tissues"
  )
)

sf_immune_gokegg_sex <- unique(
  aging_dt_immune_sex[
    ORA_CATEGORY %in% c("GO_TERMS", "KEGG_PWS"),
    c(
      "ORA_CATEGORY",
      "VALUE",
      "ORA_REGULATION_age",
      "Overall Male (Union)",
      "Overall Female (Union)"
    )
  ]
)[order(-`Overall Male (Union)`)]
sf_immune_gokegg_sex[
  ,
  ORA_CATEGORY := ifelse(
    ORA_CATEGORY == "GO_TERMS",
    "GO",
    "KEGG"
  )
]
setnames(
  sf_immune_gokegg_sex,
  old = colnames(sf_immune_gokegg_sex),
  new = c(
    "Category",
    "Term",
    "Age Regulation",
    "Over-represented in # male tissues",
    "Over-represented in # female tissues"
  )
)

# remove redundant terms
sf_immune_remove <- c(
  "lymphocyte differentiation",
  "regulation of viral life cycle",
  "regulation of immune response",
  "regulation of immune system process",
  "positive regulation of immune response	",
  "modulation by symbiont of entry into host",
  "regulation of mature B cell apoptotic process",
  "regulation of lymphocyte activation",
  "regulation of leukocyte activation",
  "regulation of cytokine production involved in immune response",
  "positive regulation of adaptive immune response",
  "regulation of leukocyte apoptotic process",
  "regulation of lymphocyte apoptotic process",
  "negative regulation of mature B cell apoptotic process",
  "regulation of B cell apoptotic process",
  "regulation of adaptive immune response",
  "negative regulation of B cell apoptotic process",
  "leukocyte differentiation",
  "leukocyte migration"
)

sf_immune_gokegg <- sf_immune_gokegg[
  !Term %in% sf_immune_remove
]

sf_immune_gokegg %>% kbl(
  align = rep("c", 4)
) %>% kable_styling(
  "striped",
  full_width = FALSE,
  html_font = "arial"
) %>% kable_styling(
  font_size = 18
) %>% save_kable(
  paste0(path_scagecom_output, "sfig_immune_a.png"),
  vwidth = 750,
  zoom = 2
)

sf_immune_gokegg_sex <- sf_immune_gokegg_sex[
  !Term %in% sf_immune_remove
]

sf_immune_gokegg_sex %>% kbl(
  align = rep("c", 5)
) %>% kable_styling(
  "striped",
  full_width = FALSE,
  html_font = "arial"
) %>% kable_styling(
  font_size = 18
) %>% save_kable(
  paste0(path_scagecom_output, "sfig_immune_a_sex_2.png"),
  vwidth = 950,
  zoom = 2
)

# prepare LRI table

sf_immune_lri <- unique(
  aging_dt_immune_sex[
    ORA_CATEGORY %in% c("LRI"),
    c(
      "VALUE", "ORA_REGULATION_age",
      "Overall (Union)_age", "pubmed", "HAGR",
      "summary_val"
      )
  ]
)[order(-`Overall (Union)_age`)]

sf_immune_lri[
  grepl("H2", VALUE),
  pubmed := paste0(pubmed, "*")
]
setnames(
  sf_immune_lri,
  old = colnames(sf_immune_lri),
  new = c(
    "LRI",
    "Age Regulation",
    "Over-represented in # tissues",
    "Associated to # aging PubMed articles",
    "In HAGR",
    "Secretomics detection (either Ligand or Receptor)"
  )
)

sf_immune_lri %>% kbl(
  align = rep("c", 6)
) %>% kable_styling(
  "striped",
  full_width = FALSE,
  html_font = "arial"
) %>% kable_styling(
  font_size = 18
) %>% save_kable(
  paste0(path_scagecom_output, "sfig_immune_b2.png"),
  vwidth = 1350,
  zoom = 2
)

## Lipid Metabolism Supplementary Figure ####

# prepare GO/KEGG table
sf_lipmed_gokegg <- unique(
  aging_dt_lipmed_sex[
    ORA_CATEGORY %in% c("GO_TERMS", "KEGG_PWS"),
    c(
      "ORA_CATEGORY",
      "VALUE",
      "ORA_REGULATION_age",
      "Overall Male (Union)",
      "Overall Female (Union)"
    )
  ]
)[order(-`Overall Male (Union)`)]
sf_lipmed_gokegg[
  ,
  ORA_CATEGORY := ifelse(
    ORA_CATEGORY == "GO_TERMS",
    "GO",
    "KEGG"
  )
]
setnames(
  sf_lipmed_gokegg,
  old = colnames(sf_lipmed_gokegg),
  new = c(
    "Category",
    "Term",
    "Age Regulation",
    "Over-represented in # male tissues",
    "Over-represented in # female tissues"
  )
)
# remove redundant terms

sf_lipmed_gokegg <- sf_lipmed_gokegg[
  !Term %in% c(
    "cellular lipid metabolic process",
    "positive regulation of lipid metabolic process",
    "regulation of lipid localization",
    "unsaturated fatty acid biosynthetic process",
    "unsaturated fatty acid metabolic process"
  )
]

sf_lipmed_gokegg %>% kbl(
  align = rep("c", 5)
) %>% kable_styling(
  "striped",
  full_width = FALSE,
  html_font = "arial"
) %>% kable_styling(
  font_size = 18
) %>% save_kable(
  paste0(path_scagecom_output, "sfig_lipmed_a2.png"),
  vwidth = 750,
  zoom = 2
)

# prepare LRI table

sf_lipmed_lri <- unique(
  aging_dt_lipmed_sex[order(VALUE)][
    ORA_CATEGORY %in% c("LRI"),
    c(
      "VALUE", "ORA_REGULATION_age",
      "Overall Male (Union)",
      "Overall Female (Union)",
      "pubmed", "HAGR",
      "summary_val"
      )
  ]
)

setnames(
  sf_lipmed_lri,
  old = colnames(sf_lipmed_lri),
  new = c(
    "LRI",
    "Age Regulation",
    "Over-represented in # male tissues",
    "Over-represented in # female tissues",
    "Associated to # aging PubMed articles",
    "In HAGR",
    "Secretomics detection (either Ligand or Receptor)"
  )
)

sf_lipmed_lri[
  ,
  LRI_reg := paste(
    LRI,
    `Age Regulation`,
    sep = "_"
  )
]

sf_lipmed_lri <- sf_lipmed_lri[
  !LRI_reg %in% c(
    "Apoa1:Ldlr_UP",
    "Apoa1:Ldlr_DOWN",
    "Apoe:Ldlr_DOWN",
    "App:Cav1_UP",
    "App:Notch1_UP",
    "App:Notch2_UP",
    "App:Slc45a3_UP",
    "App:Slc45a3_DOWN",
    "App:Sorl1_UP",
    "F13a1:Itgb1_UP",
    "Mmp2:Itgb1_UP",
    "Psen1:Ncstn_UP",
    "Rps27a:Ldlr_DOWN",
    "Rps27a:Ldlr_DOWN",
    "Tgm2:Adgrg1_UP",
    "Tgm2:Adgrg1_DOWN"
  )
]

sf_lipmed_lri[is.na(sf_lipmed_lri)] <- 0

sf_lipmed_lri[, -c(8)] %>% kbl(
  align = rep("c", 7)
) %>% kable_styling(
  "striped",
  full_width = FALSE,
  html_font = "arial"
) %>% kable_styling(
  font_size = 18
) %>% save_kable(
  paste0(path_scagecom_output, "sfig_lipmed_b2.png"),
  vwidth = 1450,
  zoom = 2
)

## Apoe Brain

sfig_apoe_brain <- fun_plot_lrfc_cci(
  dt_cci_full[
    dataset == "TMS FACS (male)" &
      tissue == "Brain" &
      LIGAND_1 == "Apoe"
  ],
  "Brain - TMS FACS (male) - Apoe Ligand",
  "Apoe Log2(FC)"
)

save_image(
  sfig_apoe_brain, 
  paste0(
    path_scagecom_output,
    "sfig_apoe_ligand2.svg"
  ),
  scale = 1,
  width = 900,
  height = 700
)


## ECM Supplementary figures ####

# prepare GO/KEGG table
sf_ecm_gokegg <- unique(
  aging_dt_ecm_sex[
    ORA_CATEGORY %in% c("GO_TERMS", "KEGG_PWS"),
    c(
      "ORA_CATEGORY",
      "VALUE",
      "ORA_REGULATION_age",
      "Overall Male (Union)",
      "Overall Female (Union)"
    )
  ]
)[order(-`Overall Male (Union)`)]
sf_ecm_gokegg[
  ,
  ORA_CATEGORY := ifelse(
    ORA_CATEGORY == "GO_TERMS",
    "GO",
    "KEGG"
  )
]
setnames(
  sf_ecm_gokegg,
  old = colnames(sf_ecm_gokegg),
  new = c(
    "Category",
    "Term",
    "Age Regulation",
    "Over-represented in # male tissues",
    "Over-represented in # female tissues"
  )
)
# remove redundant terms

sf_ecm_gokegg <- sf_ecm_gokegg[
  !Term %in% c(
    "extracellular region",
    "focal adhesion",
    "cell junction organization",
    "cell-substrate junction",
    "cell junction assembly",
    "basement membrane organization",
    "extracellular structure organization",
    "extracellular matrix organization",
    "biological adhesion"
  )
]

sf_ecm_gokegg %>% kbl(
  align = rep("c", 5)
) %>% kable_styling(
  "striped",
  full_width = FALSE,
  html_font = "arial"
) %>% kable_styling(
  font_size = 18
) %>% save_kable(
  paste0(path_scagecom_output, "sfig_ecm_a2.png"),
  vwidth = 750,
  zoom = 2
)

# prepare LRI table

sf_ecm_lri <- unique(
  aging_dt_ecm_sex[order(VALUE)][
    ORA_CATEGORY %in% c("LRI"),
    c(
      "VALUE", "ORA_REGULATION_age",
      "Overall Male (Union)",
      "Overall Female (Union)",
      "pubmed", "HAGR",
      "summary_val"
    )
  ]
)

sf_ecm_lri[, summary_val := ifelse(
  is.na(summary_val),
  "",
  summary_val
)]

setnames(
  sf_ecm_lri,
  old = colnames(sf_ecm_lri),
  new = c(
    "LRI",
    "Age Regulation",
    "Over-represented in # male tissues",
    "Over-represented in # female tissues",
    "Associated to # aging PubMed articles",
    "In HAGR",
    "Secretomics detection (either Ligand or Receptor)"
  )
)

sf_ecm_lri[
  ,
  LRI_reg := paste(
    LRI,
    `Age Regulation`,
    sep = "_"
  )
]

sf_ecm_lri <- sf_ecm_lri[
  !LRI_reg %in% c(
    "Adam15:Itga9_DOWN",
    "Adam15:Itgb1_DOWN",
    "Adam9:Itgav_DOWN",
    "Adam9:Itgb5_DOWN",
    "Col4a1:Itga9_Itgb1_DOWN",
    "Col4a1:Itgb1_Itga1_DOWN",
    "Col4a2:Itga9_Itgb1_DOWN",
    "Tgm2:Itgb1_DOWN",
    "Col1a1:Sdc4_DOWN",
    "Col1a2:Cd36_DOWN",
    "Col1a2:Cd44_DOWN",
    "Ccn4:Itgb1_DOWN",
    "Col5a2:Itgb1_Itga1_DOWN",
    "Col6a1:Cd44_DOWN"
  )
]

sf_ecm_lri[is.na(sf_ecm_lri)] <- 0

sf_ecm_lri[, -c(8)] %>% kbl(
  align = rep("c", 7)
) %>% kable_styling(
  "striped",
  full_width = FALSE,
  html_font = "arial"
) %>% kable_styling(
  font_size = 18
) %>% save_kable(
  paste0(path_scagecom_output, "sfig_ecm_b2.png"),
  vwidth = 1450,
  zoom = 2
)

# prepare ER table

sf_ecm_er <- unique(
  aging_dt_ecm_sex[order(VALUE)][
    ORA_CATEGORY %in% c("ER_CELLFAMILIES"),
    c(
      "VALUE", "ORA_REGULATION_age",
      "Overall Male (Union)",
      "Overall Female (Union)"
    )
  ]
)

setnames(
  sf_ecm_er,
  old = colnames(sf_ecm_er),
  new = c(
    "Emmitter-Receiver Cell Families",
    "Age Regulation",
    "Over-represented in # male tissues",
    "Over-represented in # female tissues"
  )
)


sf_ecm_er[is.na(sf_ecm_er)] <- 0

sf_ecm_er %>% kbl(
  align = rep("c", 7)
) %>% kable_styling(
  "striped",
  full_width = FALSE,
  html_font = "arial"
) %>% kable_styling(
  font_size = 18
) %>% save_kable(
  paste0(path_scagecom_output, "sfig_ecm_c2.png"),
  vwidth = 1450,
  zoom = 2
)

## Growth/Dev/Angio Supplementary figures ####

# prepare GO/KEGG table
sf_dev_gokegg <- unique(
  aging_dt_dev_sex[
    ORA_CATEGORY %in% c("GO_TERMS", "KEGG_PWS"),
    c(
      "ORA_CATEGORY",
      "VALUE",
      "ORA_REGULATION_age",
      "Overall Male (Union)",
      "Overall Female (Union)"
    )
  ]
)[order(-`Overall Male (Union)`)]
sf_dev_gokegg[
  ,
  ORA_CATEGORY := ifelse(
    ORA_CATEGORY == "GO_TERMS",
    "GO",
    "KEGG"
  )
]
setnames(
  sf_dev_gokegg,
  old = colnames(sf_dev_gokegg),
  new = c(
    "Category",
    "Term",
    "Age Regulation",
    "Over-represented in # male tissues",
    "Over-represented in # female tissues"
  )
)
# remove redundant terms

sf_dev_gokegg <- sf_dev_gokegg[
  !Term %in% c(
    "developmental growth involved in morphogenesis",
    "tissue morphogenesis"
  )
]

sf_dev_gokegg %>% kbl(
  align = rep("c", 5)
) %>% kable_styling(
  "striped",
  full_width = FALSE,
  html_font = "arial"
) %>% kable_styling(
  font_size = 18
) %>% save_kable(
  paste0(path_scagecom_output, "sfig_dev_a2.png"),
  vwidth = 750,
  zoom = 2
)

# prepare LRI table

sf_dev_lri <- unique(
  aging_dt_dev_sex[order(VALUE)][
    ORA_CATEGORY %in% c("LRI"),
    c(
      "VALUE", "ORA_REGULATION_age",
      "Overall Male (Union)",
      "Overall Female (Union)",
      "pubmed", "HAGR",
      "summary_val"
    )
  ]
)

sf_dev_lri[, summary_val := ifelse(
  is.na(summary_val),
  "",
  summary_val
)]

setnames(
  sf_dev_lri,
  old = colnames(sf_dev_lri),
  new = c(
    "LRI",
    "Age Regulation",
    "Over-represented in # male tissues",
    "Over-represented in # female tissues",
    "Associated to # aging PubMed articles",
    "In HAGR",
    "Secretomics detection (either Ligand or Receptor)"
  )
)

sf_dev_lri[
  ,
  LRI_reg := paste(
    LRI,
    `Age Regulation`,
    sep = "_"
  )
]

sf_dev_lri <- sf_dev_lri[
  !LRI_reg %in% c(
    ""
  )
]

sf_dev_lri[is.na(sf_dev_lri)] <- 0

sf_dev_lri[, -c(8)] %>% kbl(
  align = rep("c", 7)
) %>% kable_styling(
  "striped",
  full_width = FALSE,
  html_font = "arial"
) %>% kable_styling(
  font_size = 18
) %>% save_kable(
  paste0(path_scagecom_output, "sfig_dev_b2.png"),
  vwidth = 1450,
  zoom = 2
)

## Slpi Supplementary Figure ####

sfig_slpi <- plot_KEYWORD_summary(
  ora_keyword_summary = shiny_list_full$ORA_KEYWORD_SUMMARY,
  ora_keyword_template = shiny_list_full$ORA_KEYWORD_TEMPLATE,
  category = "Ligand-Receptor Interaction",
  keyword = "Slpi:Plscr1"
)

ggsave(
  paste0(
    path_scagecom_output,
    "sfig_slpi_a.png"
  ),
  plot = sfig_slpi,
  width = 2000,
  height = 1400,
  units = "px",
  scale = 2.7
)

