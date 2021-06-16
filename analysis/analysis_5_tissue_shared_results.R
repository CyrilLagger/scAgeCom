####################################################
##
## Project: scAgeCom
##
## Last update - June 2021
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## preparation of results shared between tissues
##
####################################################
##

## Libraries ####

library(scDiffCom)
library(data.table)

## load tissue specific results ####

data_4_tissue_specific_results <- readRDS(
  "../data_scAgeCom/analysis/outputs_data/data_4_tissue_specific_results.rds"
)

## create a tissue-counts table for all OR keywords ####

ORA_KEYWORD_COUNTS <- rbindlist(
  list(
    dcast.data.table(
      copy(data_4_tissue_specific_results$ORA_table[
        ORA_REGULATION %in% c("UP", "UP:DOWN"),
        c("Dataset", "Tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ])[, ORA_REGULATION := "UP"][
        ,
        .N,
        by = c("Dataset", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      ORA_CATEGORY + ORA_REGULATION + VALUE ~ Dataset,
      value.var = "N",
      fill = 0
    )[
      unique(
        copy(data_4_tissue_specific_results$ORA_table[
          ORA_REGULATION %in% c("UP", "UP:DOWN"),
          c("Dataset", "Tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
        ])[, ORA_REGULATION := "UP"][
          ,
          -c("Dataset")
        ]
      )[
        ,
        .N,
        by = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      on = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE"),
      `Overall (Union)` := i.N
    ],
    dcast.data.table(
      copy(data_4_tissue_specific_results$ORA_table[
        ORA_REGULATION %in% c("DOWN", "UP:DOWN"),
        c("Dataset", "Tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ])[, ORA_REGULATION := "DOWN"][
        ,
        .N,
        by = c("Dataset", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      ORA_CATEGORY + ORA_REGULATION + VALUE ~ Dataset,
      value.var = "N",
      fill = 0
    )[
      unique(
        copy(data_4_tissue_specific_results$ORA_table[
          ORA_REGULATION %in% c("DOWN", "UP:DOWN"),
          c("Dataset", "Tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
        ])[, ORA_REGULATION := "DOWN"][
          ,
          -c("Dataset")
        ]
      )[
        ,
        .N,
        by = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      on = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE"),
      `Overall (Union)` := i.N
    ],
    dcast.data.table(
      copy(data_4_tissue_specific_results$ORA_table[
        ORA_REGULATION %in% c("FLAT"),
        c("Dataset", "Tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ])[, ORA_REGULATION := "FLAT"][
        ,
        .N,
        by = c("Dataset", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      ORA_CATEGORY + ORA_REGULATION + VALUE ~ Dataset,
      value.var = "N",
      fill = 0
    )[
      unique(
        copy(data_4_tissue_specific_results$ORA_table[
          ORA_REGULATION %in% c("FLAT"),
          c("Dataset", "Tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
        ])[, ORA_REGULATION := "FLAT"][
          ,
          -c("Dataset")
        ]
      )[
        ,
        .N,
        by = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      on = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE"),
      `Overall (Union)` := i.N
    ]
  )
)

ORA_KEYWORD_COUNTS[
  data.table(
    old_category = c(
      "LRI",
      "LIGAND_COMPLEX",
      "RECEPTOR_COMPLEX",
      "ER_CELLTYPES",
      "EMITTER_CELLTYPE",
      "RECEIVER_CELLTYPE",
      "GO_TERMS",
      "KEGG_PWS",
      "ER_CELLFAMILIES",
      "EMITTER_CELLFAMILY",
      "RECEIVER_CELLFAMILY"
    ),
    new_category = c(
      "Ligand-Receptor Interaction",
      "Ligand",
      "Receptor",
      "Emitter-Receiver Cell Type",
      "Emitter Cell Type",
      "Receiver Cell Type",
      "GO Term",
      "KEGG Pathway",
      "Emitter-Receiver Cell Type Family",
      "Emitter Cell Type Family",
      "Receiver Cell Type Family"
    )
  ),
  on = c("ORA_CATEGORY==old_category"),
  ORA_CATEGORY := i.new_category
]

setcolorder(
  ORA_KEYWORD_COUNTS,
  c(
    "ORA_CATEGORY",
    "ORA_REGULATION",
    "VALUE",
    "Overall (Union)",
    "TMS FACS (male)",
    "TMS FACS (female)",
    "TMS Droplet (male)",
    "TMS Droplet (female)",
    "Calico Droplet (male)"
  )
)

ORA_KEYWORD_COUNTS[
  unique(
    data_4_tissue_specific_results$ORA_table[
      ORA_CATEGORY == "GO_TERMS",
      c("VALUE", "ASPECT", "GO Level")
    ]
  ),
  on = "VALUE",
  c("ASPECT", "GO Level") := mget(paste0("i.", c("ASPECT", "GO Level")))
]

ORA_KEYWORD_COUNTS[, `GO Level` := factor(`GO Level`)]

setkey(ORA_KEYWORD_COUNTS)
setorder(
  ORA_KEYWORD_COUNTS,
  -`Overall (Union)`
)

ORA_KEYWORD_COUNTS[
  ,
  `Overall (Union)`:= factor(
    ifelse(
      
      `Overall (Union)` < 10 & `Overall (Union)` > 0 ,
      paste0("0", `Overall (Union)`, "/23"),
      paste0(`Overall (Union)`, "/23")
    )
  ),
]
ORA_KEYWORD_COUNTS[
  ,
  `TMS FACS (male)`:= factor(
    ifelse(
      `TMS FACS (male)` < 10 & `TMS FACS (male)` > 0 ,
      paste0("0", `TMS FACS (male)`, "/21"),
      paste0(`TMS FACS (male)`, "/21")
    )
  ),
]
ORA_KEYWORD_COUNTS[
  ,
  `TMS FACS (female)`:= factor(
    ifelse(
      `TMS FACS (female)` < 10 & `TMS FACS (female)` > 0 ,
      paste0("0", `TMS FACS (female)`, "/19"),
      paste0(`TMS FACS (female)`, "/19")
    )
  ),
]
ORA_KEYWORD_COUNTS[
  ,
  `TMS Droplet (male)`:= factor(
    paste0(`TMS Droplet (male)`, "/6")
  ),
]
ORA_KEYWORD_COUNTS[
  ,
  `TMS Droplet (female)` :=  factor(
    paste0(`TMS Droplet (female)`, "/9")
  ),
]
ORA_KEYWORD_COUNTS[
  ,
  `Calico Droplet (male)` := factor(
    paste0(`Calico Droplet (male)`, "/3")
  ),
]

## Create a "tissue vs dataset" summary table for each OR keyword ####

ORA_KEYWORD_SUMMARY <- melt.data.table(
  dcast.data.table(
    data_4_tissue_specific_results$ORA_table[
      ,
      c("ORA_CATEGORY", "VALUE", "Tissue", "Dataset" , "ORA_REGULATION")
    ],
    ORA_CATEGORY + VALUE + Tissue ~ Dataset,
    value.var = "ORA_REGULATION",
    fill = "Not Detected"
  ),
  id.vars = c("ORA_CATEGORY", "VALUE", "Tissue"),
  variable.name = "Dataset",
  value.name = "ORA_REGULATION"
)[
  ,
  Dataset_Tissue := paste(
    Dataset,
    Tissue,
    sep = "_"
  )
][
  ,
  ORA_REGULATION := ifelse(
    Dataset_Tissue %in% unique(
      data_4_tissue_specific_results$CCI_table[
        ,
        c("Dataset", "Tissue")
      ]
    )[
      , 
      Dataset_Tissue := paste(
        Dataset, Tissue, sep = "_"
      )
    ]$Dataset_Tissue,
    ORA_REGULATION,
    "No Data"
  )
]

ORA_KEYWORD_SUMMARY[
  data.table(
    old_category = c(
      "LRI",
      "LIGAND_COMPLEX",
      "RECEPTOR_COMPLEX",
      "ER_CELLTYPES",
      "EMITTER_CELLTYPE",
      "RECEIVER_CELLTYPE",
      "GO_TERMS",
      "KEGG_PWS",
      "ER_CELLFAMILIES",
      "EMITTER_CELLFAMILY",
      "RECEIVER_CELLFAMILY"
    ),
    new_category = c(
      "Ligand-Receptor Interaction",
      "Ligand",
      "Receptor",
      "Emitter-Receiver Cell Type",
      "Emitter Cell Type",
      "Receiver Cell Type",
      "GO Term",
      "KEGG Pathway",
      "Emitter-Receiver Cell Type Family",
      "Emitter Cell Type Family",
      "Receiver Cell Type Family"
    )
  ),
  on = c("ORA_CATEGORY==old_category"),
  ORA_CATEGORY := i.new_category
]

ORA_KEYWORD_SUMMARY[
  unique(
    data_4_tissue_specific_results$ORA_table[
      ORA_CATEGORY == "GO_TERMS",
      c("VALUE", "ASPECT", "GO Level")
    ]
  ),
  on = "VALUE",
  c("ASPECT", "GO Level") := mget(paste0("i.", c("ASPECT", "GO Level")))
]
ORA_KEYWORD_SUMMARY[, `GO Level` := factor(`GO Level`)]

ORA_KEYWORD_TEMPLATE <- unique(
  ORA_KEYWORD_SUMMARY[
    ORA_REGULATION != "No Data"
    ,
    c("Tissue", "Dataset")
  ]
)

## create vectors to access categories in shiny ####

ALL_ORA_CATEGORIES_GLOBAL <- c(
  "By GO/KEGG",
  "By Genes",
  "By Cell Type Families"
)

ALL_ORA_CATEGORIES_KEYWORD <- c(
  "Ligand-Receptor Interaction",
  "Ligand",
  "Receptor",
  "GO Term",
  "KEGG Pathway",
  "Emitter-Receiver Cell Type Family",
  "Emitter Cell Type Family",
  "Receiver Cell Type Family"
)

## save all results ####

data_5_tissue_shared_results <- list(
  ORA_KEYWORD_COUNTS = ORA_KEYWORD_COUNTS,
  ORA_KEYWORD_SUMMARY = ORA_KEYWORD_SUMMARY,
  ORA_KEYWORD_TEMPLATE = ORA_KEYWORD_TEMPLATE,
  ALL_ORA_CATEGORIES_GLOBAL = ALL_ORA_CATEGORIES_GLOBAL,
  ALL_ORA_CATEGORIES_KEYWORD = ALL_ORA_CATEGORIES_KEYWORD
)

saveRDS(
  data_5_tissue_shared_results,
  "../data_scAgeCom/analysis/outputs_data/data_5_tissue_shared_results.rds"
)
