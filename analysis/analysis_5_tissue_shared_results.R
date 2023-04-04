####################################################
##
## Project: scAgeCom
##
## lagger.cyril@gmail.com
## ursu_eugen@hotmail.com
## anais.equey@gmail.com
##
## Process cross-tissue results
##
####################################################
##

## Create a tissue-counts table for all ORA keywords ####

dt_ora_key_counts <- rbindlist(
  list(
    dcast.data.table(
      copy(dt_ora_full[!grepl("mixed", dataset)][
        ORA_REGULATION %in% c("UP", "UP:DOWN"),
        c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ])[, ORA_REGULATION := "UP"][
        ,
        .N,
        by = c("dataset", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      ORA_CATEGORY + ORA_REGULATION + VALUE ~ dataset,
      value.var = "N",
      fill = 0
    )[
      unique(
        copy(dt_ora_full[!grepl("mixed", dataset)][
          ORA_REGULATION %in% c("UP", "UP:DOWN"),
          c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
        ])[, ORA_REGULATION := "UP"][
          ,
          -c("dataset")
        ]
      )[
        ,
        .N,
        by = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      on = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE"),
      `Overall (Union)` := i.N
    ][
      unique(
        copy(dt_ora_full[
          dataset %in% c(
            "Calico Droplet (male)", "TMS Droplet (male)", "TMS FACS (male)"
          )
        ][
          ORA_REGULATION %in% c("UP", "UP:DOWN"),
          c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
        ])[, ORA_REGULATION := "UP"][
          ,
          -c("dataset")
        ]
      )[
        ,
        .N,
        by = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      on = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE"),
      `Overall Male (Union)` := i.N
    ][
      unique(
        copy(dt_ora_full[
          dataset %in% c(
            "TMS Droplet (female)", "TMS FACS (female)"
          )
        ][
          ORA_REGULATION %in% c("UP", "UP:DOWN"),
          c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
        ])[, ORA_REGULATION := "UP"][
          ,
          -c("dataset")
        ]
      )[
        ,
        .N,
        by = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      on = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE"),
      `Overall Female (Union)` := i.N
    ],
    dcast.data.table(
      copy(dt_ora_full[!grepl("mixed", dataset)][
        ORA_REGULATION %in% c("DOWN", "UP:DOWN"),
        c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ])[, ORA_REGULATION := "DOWN"][
        ,
        .N,
        by = c("dataset", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      ORA_CATEGORY + ORA_REGULATION + VALUE ~ dataset,
      value.var = "N",
      fill = 0
    )[
      unique(
        copy(dt_ora_full[!grepl("mixed", dataset)][
          ORA_REGULATION %in% c("DOWN", "UP:DOWN"),
          c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
        ])[, ORA_REGULATION := "DOWN"][
          ,
          -c("dataset")
        ]
      )[
        ,
        .N,
        by = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      on = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE"),
      `Overall (Union)` := i.N
    ][
      unique(
        copy(dt_ora_full[
          dataset %in% c(
            "Calico Droplet (male)", "TMS Droplet (male)", "TMS FACS (male)"
          )
        ][
          ORA_REGULATION %in% c("DOWN", "UP:DOWN"),
          c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
        ])[, ORA_REGULATION := "DOWN"][
          ,
          -c("dataset")
        ]
      )[
        ,
        .N,
        by = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      on = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE"),
      `Overall Male (Union)` := i.N
    ][
      unique(
        copy(dt_ora_full[
          dataset %in% c(
            "TMS Droplet (female)", "TMS FACS (female)"
          )
        ][
          ORA_REGULATION %in% c("DOWN", "UP:DOWN"),
          c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
        ])[, ORA_REGULATION := "DOWN"][
          ,
          -c("dataset")
        ]
      )[
        ,
        .N,
        by = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      on = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE"),
      `Overall Female (Union)` := i.N
    ],
    dcast.data.table(
      copy(dt_ora_full[!grepl("mixed", dataset)][
        ORA_REGULATION %in% c("FLAT"),
        c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ])[, ORA_REGULATION := "FLAT"][
        ,
        .N,
        by = c("dataset", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      ORA_CATEGORY + ORA_REGULATION + VALUE ~ dataset,
      value.var = "N",
      fill = 0
    )[
      unique(
        copy(dt_ora_full[!grepl("mixed", dataset)][
          ORA_REGULATION %in% c("FLAT"),
          c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
        ])[, ORA_REGULATION := "FLAT"][
          ,
          -c("dataset")
        ]
      )[
        ,
        .N,
        by = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      on = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE"),
      `Overall (Union)` := i.N
    ][
      unique(
        copy(dt_ora_full[
          dataset %in% c(
            "Calico Droplet (male)", "TMS Droplet (male)", "TMS FACS (male)"
          )
        ][
          ORA_REGULATION %in% c("FLAT"),
          c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
        ])[, ORA_REGULATION := "FLAT"][
          ,
          -c("dataset")
        ]
      )[
        ,
        .N,
        by = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      on = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE"),
      `Overall Male (Union)` := i.N
    ][
      unique(
        copy(dt_ora_full[
          dataset %in% c(
            "TMS Droplet (female)", "TMS FACS (female)"
          )
        ][
          ORA_REGULATION %in% c("FLAT"),
          c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
        ])[, ORA_REGULATION := "FLAT"][
          ,
          -c("dataset")
        ]
      )[
        ,
        .N,
        by = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      on = c("ORA_CATEGORY", "ORA_REGULATION", "VALUE"),
      `Overall Female (Union)` := i.N
    ]
  )
)

## Create a "tissue vs dataset" summary table for all ORA keywords ####

dt_ora_key_summary <- melt.data.table(
  dcast.data.table(
    dt_ora_full[!grepl("mixed", dataset)][
      ,
      c("ORA_CATEGORY", "VALUE", "tissue", "dataset", "ORA_REGULATION")
    ],
    ORA_CATEGORY + VALUE + tissue ~ dataset,
    value.var = "ORA_REGULATION",
    fill = "Not Detected"
  ),
  id.vars = c("ORA_CATEGORY", "VALUE", "tissue"),
  variable.name = "dataset",
  value.name = "ORA_REGULATION"
)[
  ,
  dataset_tissue := paste(
    dataset,
    tissue,
    sep = "_"
  )
][
  ,
  ORA_REGULATION := ifelse(
    dataset_tissue %in% unique(
      dt_ora_full[!grepl("mixed", dataset)][
        ,
        c("dataset", "tissue")
      ]
    )[
      ,
      dataset_tissue := paste(
        dataset, tissue, sep = "_"
      )
    ]$dataset_tissue,
    ORA_REGULATION,
    "No Data"
  )
]

## Create a diffsex tissue-counts table for all ORA keywords ####

dt_ora_key_counts_diffsex <- rbindlist(
  list(
    dcast.data.table(
      copy(dt_ora_full_diffsex[!grepl("combined", dataset)][
        ORA_REGULATION %in% c("UP", "UP:DOWN"),
        c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ])[, ORA_REGULATION := "UP"][
        ,
        .N,
        by = c("dataset", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      ORA_CATEGORY + ORA_REGULATION + VALUE ~ dataset,
      value.var = "N",
      fill = 0
    )[
      unique(
        copy(dt_ora_full_diffsex[!grepl("combined", dataset)][
          ORA_REGULATION %in% c("UP", "UP:DOWN"),
          c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
        ])[, ORA_REGULATION := "UP"][
          ,
          -c("dataset")
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
      copy(dt_ora_full_diffsex[!grepl("combined", dataset)][
        ORA_REGULATION %in% c("DOWN", "UP:DOWN"),
        c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ])[, ORA_REGULATION := "DOWN"][
        ,
        .N,
        by = c("dataset", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      ORA_CATEGORY + ORA_REGULATION + VALUE ~ dataset,
      value.var = "N",
      fill = 0
    )[
      unique(
        copy(dt_ora_full_diffsex[!grepl("combined", dataset)][
          ORA_REGULATION %in% c("DOWN", "UP:DOWN"),
          c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
        ])[, ORA_REGULATION := "DOWN"][
          ,
          -c("dataset")
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
      copy(dt_ora_full_diffsex[!grepl("combined", dataset)][
        ORA_REGULATION %in% c("FLAT"),
        c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ])[, ORA_REGULATION := "FLAT"][
        ,
        .N,
        by = c("dataset", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
      ],
      ORA_CATEGORY + ORA_REGULATION + VALUE ~ dataset,
      value.var = "N",
      fill = 0
    )[
      unique(
        copy(dt_ora_full_diffsex[!grepl("combined", dataset)][
          ORA_REGULATION %in% c("FLAT"),
          c("dataset", "tissue", "ORA_CATEGORY", "ORA_REGULATION", "VALUE")
        ])[, ORA_REGULATION := "FLAT"][
          ,
          -c("dataset")
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

## Create a diffsex "tissue vs dataset" summary table for all ORA keywords ####

dt_ora_key_summary_diffsex <- melt.data.table(
  dcast.data.table(
    dt_ora_full_diffsex[
      ,
      c("ORA_CATEGORY", "VALUE", "tissue", "dataset", "ORA_REGULATION")
    ],
    ORA_CATEGORY + VALUE + tissue ~ dataset,
    value.var = "ORA_REGULATION",
    fill = "Not Detected"
  ),
  id.vars = c("ORA_CATEGORY", "VALUE", "tissue"),
  variable.name = "dataset",
  value.name = "ORA_REGULATION"
)[
  ,
  dataset_tissue := paste(
    dataset,
    tissue,
    sep = "_"
  )
][
  ,
  ORA_REGULATION := ifelse(
    dataset_tissue %in% unique(
      dt_ora_full_diffsex[
        ,
        c("dataset", "tissue")
      ]
    )[
      ,
      dataset_tissue := paste(
        dataset, tissue, sep = "_"
      )
    ]$dataset_tissue,
    ORA_REGULATION,
    "No Data"
  )
]
