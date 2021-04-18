####################################################
##
## Project: scAgeCom
##
## Last update - April 2021
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
library(ggplot2)
library(htmltools)
library(plotly)

## load tissue specific results

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
      "LRIs",
      "Ligand Gene(s)",
      "Receptor Gene(s)",
      "ER Cell Types",
      "Emitter Cell Types",
      "Receiver Cell Types",
      "GO Terms",
      "KEGG Pathways",
      "ER Cell Families",
      "Emitter Cell Families",
      "Receiver Cell Families"
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
    "Calico2019"
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
  Calico2019 := factor(
    paste0(Calico2019, "/3")
  ),
]



# ORA_KEYWORD_COUNTS[
#   ,
#   `TMS Droplet (male)`:= ifelse(
#     `TMS Droplet (male)` < 10 & `TMS Droplet (male)` > 0 ,
#     paste0("0", `TMS Droplet (male)`, "/6"),
#     paste0(`TMS Droplet (male)`, "/6")
#   ),
# ]
# ORA_KEYWORD_COUNTS[
#   ,
#   `TMS Droplet (female)` := ifelse(
#     `TMS Droplet (female)` < 10 & `TMS Droplet (female)` > 0 ,
#     paste0("0", `TMS Droplet (female)`, "/9"),
#     paste0(`TMS Droplet (female)`, "/9")
#   ),
# ]
# ORA_KEYWORD_COUNTS[
#   ,
#   Calico2019 := ifelse(
#     Calico2019 < 10 & Calico2019 > 0 ,
#     paste0("0", Calico2019, "/3"),
#     paste0(Calico2019, "/3")
#   ),
# ]

## create utility functions for ORA_KEYWORD_COUNTS table on shiny ####

build_KEYWORD_COUNTS_display <- function(
  ora_keyword_counts,
  category,
  regulation,
  go_aspect = NULL
) {
  dt <- ora_keyword_counts[
    ORA_CATEGORY == category &
      ORA_REGULATION == regulation
  ]
  setnames(dt, old = "VALUE", new = category)
  if (category == "GO Terms") {
    temp_aspect <- ifelse(
      go_aspect == "Biological Process",
      "biological_process",
      ifelse(
        go_aspect == "Molecular Function",
        "molecular_function",
        "cellular_component"
      )
    )
    dt <- dt[ASPECT == temp_aspect]
    dt <- dt[, -c(1,2,10)]
  } else {
    dt <- dt[, -c(1,2,10,11)]
  }
  DT <- DT::datatable(
    data = dt,
    filter = list(
      position ="top",
      clear = FALSE,
      plain = FALSE
    ),
    class = "display compact",
    options =list(
      pageLength = 10,
      columnDefs = list(
        list(width = '300px', targets = c(1))
      )
    ),
    caption = tags$caption(
      style = 'caption-side: top; text-align: center; color:black; font-size:100% ;',
      paste0(
        "Number of tissues in which ",
        category,
        " are over-represented among ",
        regulation,
        " CCIs."
      )
    )
  )  %>%
    DT::formatStyle(
      colnames(dt)[-1],
      `text-align` = 'center'
    )
  if (category == "GO Terms") {
    DT <- DT %>% DT::formatStyle(c(7), `border-right` = "solid 2px")
  }
  DT
}

# build_KEYWORD_COUNTS_display(
#   ora_keyword_counts = ORA_KEYWORD_COUNTS,
#   category = "LRIs",
#   regulation = "DOWN",
#   go_aspect = "Biological Process"
# )


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
      "LRIs",
      "Ligand Gene(s)",
      "Receptor Gene(s)",
      "ER Cell Types",
      "Emitter Cell Types",
      "Receiver Cell Types",
      "GO Terms",
      "KEGG Pathways",
      "ER Cell Families",
      "Emitter Cell Families",
      "Receiver Cell Families"
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

## create utility functions for ORA_KEYWORD_SUMMARY table on shiny ####

plot_KEYWORD_SUMMARY <- function(
  ora_keyword_summary,
  ora_keyword_template,
  category,
  keyword
) {
  dt <- copy(ora_keyword_summary)[
    ORA_CATEGORY == category
  ]
  if (!(keyword %in% dt$VALUE)) {
    stop("`keyword` not found")
  }
  dt <- dt[
    VALUE == keyword
  ]
  dt <- copy(ora_keyword_template)[
    dt,
    on = c("Tissue", "Dataset"),
    Regulation := i.ORA_REGULATION
  ]
  dt$Dataset <- gsub(" ", "\n", dt$Dataset)
  dt[is.na(dt)] <- "Not Detected"
  p <- ggplot(dt) +
    geom_tile(
      aes(
        Dataset,
        Tissue,
        fill = Regulation,
        width = 0.9,
        height = 0.9
      ),
      colour = "black"
    ) +
    scale_fill_manual(
      values = c(
        "No Data" = "transparent",
        "Not Over-represented" = "white",
        "Not Detected" = "gray",
        "UP" = "red",
        "DOWN" = "blue",
        "FLAT" = "green",
        "UP:DOWN" = "yellow"
      )
    ) +
    ggtitle(
      substr(
        paste0(
          "Over-representation of ",
          keyword
        ),
        1,
        50
      )
    ) +
    scale_x_discrete(
      limits = c(
        "TMS\nFACS\n(male)",
        "TMS\nFACS\n(female)" ,
        "TMS\nDroplet\n(male)",
        "TMS\nDroplet\n(female)",
        "Calico2019"
      ),
      guide = guide_axis(n.dodge = 2)
    ) +
    scale_y_discrete(
      limits = sort(
        unique(dt$Tissue),
        decreasing = TRUE
      )
    ) +
    xlab("") +
    ylab("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text=element_text(size = 12)) +
    theme(axis.text=element_text(size = 12))
  plotly::ggplotly(
    p,
    source = "TCA_PLOT_KEYWORD_SUMMARY",
    tooltip = c("Dataset", "Tissue", "Regulation")
  )
}

plot_KEYWORD_SUMMARY(
  ORA_KEYWORD_SUMMARY,
  ORA_KEYWORD_TEMPLATE,
  "GO Terms",
  "A band"
)

## create vectors to access categories in shiny ####

ALL_ORA_CATEGORIES_GLOBAL <- c(
  "By GO/KEGG",
  "By Genes",
  "By Cell Type Families"
)


ALL_ORA_CATEGORIES_KEYWORD <- c(
  "LRIs",
  "Ligand Gene(s)",
  "Receptor Gene(s)",
  #"ER Cell Types",
  #"Emitter Cell Types",
  #"Receiver Cell Types",
  "GO Terms",
  "KEGG Pathways",
  "ER Cell Families",
  "Emitter Cell Families",
  "Receiver Cell Families"
)

## save all results ####

data_5_tissue_shared_results <- list(
  ORA_KEYWORD_COUNTS = ORA_KEYWORD_COUNTS,
  ORA_KEYWORD_SUMMARY = ORA_KEYWORD_SUMMARY,
  ORA_KEYWORD_TEMPLATE = ORA_KEYWORD_TEMPLATE,
  ALL_ORA_CATEGORIES_GLOBAL = ALL_ORA_CATEGORIES_GLOBAL,
  ALL_ORA_CATEGORIES_KEYWORD = ALL_ORA_CATEGORIES_KEYWORD,
  build_KEYWORD_COUNTS_display = build_KEYWORD_COUNTS_display,
  plot_KEYWORD_SUMMARY = plot_KEYWORD_SUMMARY
)

saveRDS(
  data_5_tissue_shared_results,
  "../data_scAgeCom/analysis/outputs_data/data_5_tissue_shared_results.rds"
)
