####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## additional figure preparation
##
####################################################
##

## scDiffCom vs SDGA Supplementary Figure ####

dt_sdea_fig <- dt_cci_sdea_comp[
  ,
  c("dataset", "scDiffCom_regulation", "ligand_deg", "receptor_deg")
]
dt_sdea_fig[, "CCI Not DE" := scDiffCom_regulation == "NOT DE"]
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
    text = element_text(size = 22),
    plot.title = element_text(size = 22)
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
    "Additive effect", "Substractive effect",
    "Inconsistent"
  ),
  values = c(
    "Consistent" = "#1e81b0",
    "Resolve ambiguity" = "#76b5c5",
    "Additive effect" = "#eab676",
    "Substractive effect" = "#e28743",
    "Inconsistent" = "#873e23"
  )
) + theme(
  legend.position = c(0.8, 2.8)
)
sfig_sdea_comp

ggsave(
  paste0(
    path_scagecom_output,
    "sfig_sdea_comp.png"
  ),
  plot = sfig_sdea_comp,
  width = 2100,
  height = 1200,
  units = "px",
  scale = 3
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
  paste0(path_scagecom_output, "sfig_immune_a_sex.png"),
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
  paste0(path_scagecom_output, "sfig_immune_b.png"),
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
  paste0(path_scagecom_output, "sfig_lipmed_a.png"),
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
  paste0(path_scagecom_output, "sfig_lipmed_b.png"),
  vwidth = 1450,
  zoom = 2
)

## Extended Data LRI distribution (not included) ####

ExtFig1 <- ggplot(
  data = NLRI_template,
  aes(
    y = tissue,
    x = N
  )
) + geom_boxplot(
  outlier.shape = NA
) + facet_wrap(
  ~ dataset,
  ncol = 5
) + ggplot2::scale_y_discrete(
    limits = sort(
      unique(NLRI_table$tissue),
      decreasing = TRUE
    )
) + xlab(
  "Number of detected LRIs per cell-type pair"
) + geom_jitter(
  size = 0.2
) + theme(
  #text = element_text(size = 40, face = "bold"),
  #axis.text.x = element_text(size = 30),
  #axis.text.y = element_text(size = 36, face = "bold"),
  #axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
  panel.spacing = unit(2, "lines")
)
ExtFig1
