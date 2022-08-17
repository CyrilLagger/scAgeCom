####################################################
##
## Project: scAgeCom
##
## lagger.cyril@gmail.com
## ursu_eugen@hotmail.com
## anais.equey@gmail.com
##
## Extract aging results to discuss
##
####################################################
##


## Diagram ####

table(dt_cci_age$REGULATION)
table(dt_cci_sex$REGULATION)

dt_cci_age[
    ,
    tissue_cci := paste(tissue, CCI, sep = "_")
]
dt_cci_sex[
    ,
    tissue_cci := paste(tissue, CCI, sep = "_")
]

test1 <- dt_cci_age[grepl("TMS FACS", dataset)]$tissue_cci
test2 <- dt_cci_sex[grepl("TMS FACS", dataset)]$tissue_cci

any(duplicated(test1))

length(test1)
length(test2)
length(intersect(test1, test2))

test <- merge.data.table(
    dcast.data.table(
        dt_cci_age[
            grepl("TMS FACS", dataset) &
            tissue_cci %in% dt_cci_sex[grepl("TMS FACS", dataset)]$tissue_cci,
            c("dataset", "tissue_cci", "REGULATION")
        ],
        formula = tissue_cci ~ dataset,
        value.var = "REGULATION",
        fill = "NOT_DETECTED"
    ),
    dcast.data.table(
        dt_cci_sex[
            grepl("TMS FACS", dataset) &
            tissue_cci %in% dt_cci_age[grepl("TMS FACS", dataset)]$tissue_cci,
            c("dataset", "tissue_cci", "REGULATION")
        ],
        formula = tissue_cci ~ dataset,
        value.var = "REGULATION",
        fill = "NOT_DETECTED"
    ),
    by = "tissue_cci",
    all = TRUE
)

test[
   `TMS FACS (young)` == "FLAT" &
   `TMS FACS (male)` == "UP" &
   `TMS FACS (female)` == "DOWN"
][, 1]

test[
   `TMS FACS (young)` == "FLAT" &
   `TMS FACS (male)` == "DOWN" &
   `TMS FACS (female)` == "UP"
][, 1]

test[
   `TMS FACS (young)` == "FLAT" &
   `TMS FACS (male)` == "UP" &
   `TMS FACS (female)` == "FLAT"
][, 1]

test[
   `TMS FACS (young)` == "FLAT" &
   `TMS FACS (male)` == "DOWN" &
   `TMS FACS (female)` == "FLAT"
][, 1]

test[
   `TMS FACS (young)` == "FLAT" &
   `TMS FACS (male)` == "FLAT" &
   `TMS FACS (female)` == "UP"
][, 1]

table(test[
   `TMS FACS (young)` == "FLAT" &
   `TMS FACS (male)` == "FLAT" &
   `TMS FACS (female)` == "DOWN"
]$`TMS FACS (old)`)


test2 <- test[
    ,
    .N,
    by = c("TMS FACS (female)", "TMS FACS (male)", "TMS FACS (young)", "TMS FACS (old)")
][order(-N)]

test3 <- test[grepl("App:Lrp10", tissue_cci, fixed = TRUE)]
test3 <- test[grepl("Apoe:Sdc4", tissue_cci, fixed = TRUE)]
test4 <- test3[
    ,
    .N,
    by = c("TMS FACS (female)", "TMS FACS (male)", "TMS FACS (young)", "TMS FACS (old)")
][order(-N)]


test <- dt_ora_key_summary_diffsex[
    VALUE == "App" &
    !grepl("combined", dataset)]

test2 <- dt_ora_key_summary[
    VALUE == "App" &
    tissue == "Brain"]


## Distribution of cells across age vs sex for each tissue ####

dt_md_tms_age_sex_dc <- dcast.data.table(
    dt_md_tms_age_sex,
    dataset + tissue + cell_ontology_final ~ sex + age_group,
    value.var = "N",
    fill = 0
)
dt_md_tms_age_sex_dc[
    ,
    OR := (female_OLD / male_OLD) / (female_YOUNG / male_YOUNG)
]
dt_md_tms_age_sex_dc[
  ,
  tissue := ifelse(
    tissue == "BAT", "Adipose_Brown",
    ifelse(
      tissue == "GAT", "Adipose_Gonadal",
      ifelse(
        tissue == "MAT", "Adipose_Mesenteric",
        ifelse(
          tissue == "SCAT", "Adipose_Subcutaneous",
          ifelse(
            tissue == "Brain_Myeloid", "Brain",
            ifelse(
              tissue == "Brain_Non-Myeloid", "Brain",
              tissue
            )
          )
        )
      )
    )
  )
]
dt_md_tms_age_sex_dc[
  ,
  dataset_male := ifelse(
    dataset == "tms_facs",
    paste("TMS FACS (male)", tissue, cell_ontology_final, sep = "_"),
    paste("TMS Droplet (male)", tissue, cell_ontology_final, sep = "_")
  )
]
dt_md_tms_age_sex_dc[
  ,
  dataset_female := ifelse(
    dataset == "tms_facs",
    paste("TMS FACS (female)", tissue, cell_ontology_final, sep = "_"),
    paste("TMS Droplet (female)", tissue, cell_ontology_final, sep = "_")
  )
]
dt_md_tms_age_sex_dc[
    ,
    male_scd := dataset_male %in% unique(
        paste(
            dt_cci_rel[grepl("(male)", dataset, fixed = TRUE)]$dataset,
            dt_cci_rel[grepl("(male)", dataset, fixed = TRUE)]$tissue,
            dt_cci_rel[
                grepl("(male)", dataset, fixed = TRUE)
            ]$EMITTER_CELLTYPE,
            sep = "_"
        )
    )
]
dt_md_tms_age_sex_dc[
    ,
    female_scd := dataset_female %in% unique(
        paste(
            dt_cci_rel[grepl("(female)", dataset, fixed = TRUE)]$dataset,
            dt_cci_rel[grepl("(female)", dataset, fixed = TRUE)]$tissue,
            dt_cci_rel[
                grepl("(female)", dataset, fixed = TRUE)
            ]$EMITTER_CELLTYPE,
            sep = "_"
        )
    )
]

fwrite(
    dt_md_tms_age_sex_dc,
    paste0(
        path_scagecom_output,
        "Supplementary_DATA_sex_distribution.csv"
    )
)

## Define tissues and cell types of interest for sex-related analyses ####

dt_md_tms_age_sex_dc[
    OR <= 5 & OR >= 0.2 & male_scd & female_scd
][, .N, by = c("dataset", "tissue")][order(-N)][1:10]

#tms facs lung

dt_md_tms_age_sex_dc[
    dataset == "tms_facs" & tissue == "Lung" &
    OR <= 5 & OR >= 0.2 & male_scd & female_scd
]$cell_ontology_final


test <- table(
    dt_cci_full[
        grepl("TMS FACS", dataset) & tissue == "Lung" &
        EMITTER_CELLTYPE %in% dt_md_tms_age_sex_dc[
            dataset == "tms_facs" & tissue == "Lung" &
            OR <= 5 & OR >= 0.2 & male_scd & female_scd
        ]$cell_ontology_final &
        RECEIVER_CELLTYPE %in% dt_md_tms_age_sex_dc[
            dataset == "tms_facs" & tissue == "Lung" &
            OR <= 5 & OR >= 0.2 & male_scd & female_scd
        ]$cell_ontology_final
    ]$dataset,
        dt_cci_full[
        grepl("TMS FACS", dataset) & tissue == "Lung" &
        EMITTER_CELLTYPE %in% dt_md_tms_age_sex_dc[
            dataset == "tms_facs" & tissue == "Lung" &
            OR <= 5 & OR >= 0.2 & male_scd & female_scd
        ]$cell_ontology_final &
        RECEIVER_CELLTYPE %in% dt_md_tms_age_sex_dc[
            dataset == "tms_facs" & tissue == "Lung" &
            OR <= 5 & OR >= 0.2 & male_scd & female_scd
        ]$cell_ontology_final
    ]$REGULATION
)


table(
    dt_cci_full[
        grepl("TMS FACS", dataset) & tissue == "Lung" &
        EMITTER_CELLTYPE %in% dt_md_tms_age_sex_dc[
            dataset == "tms_facs" & tissue == "Lung" &
            OR <= 5 & OR >= 0.2 & male_scd & female_scd
        ]$cell_ontology_final &
        RECEIVER_CELLTYPE %in% dt_md_tms_age_sex_dc[
            dataset == "tms_facs" & tissue == "Lung" &
            OR <= 5 & OR >= 0.2 & male_scd & female_scd
        ]$cell_ontology_final
    ]$dataset,
        dt_cci_full[
        grepl("TMS FACS", dataset) & tissue == "Lung" &
        EMITTER_CELLTYPE %in% dt_md_tms_age_sex_dc[
            dataset == "tms_facs" & tissue == "Lung" &
            OR <= 5 & OR >= 0.2 & male_scd & female_scd
        ]$cell_ontology_final &
        RECEIVER_CELLTYPE %in% dt_md_tms_age_sex_dc[
            dataset == "tms_facs" & tissue == "Lung" &
            OR <= 5 & OR >= 0.2 & male_scd & female_scd
        ]$cell_ontology_final
    ]$REGULATION
)


table(dt_cci_full_diffsex$dataset)

test2 <- table(
    dt_cci_full_diffsex[
        grepl("TMS FACS", dataset) & tissue == "Lung" &
        EMITTER_CELLTYPE %in% dt_md_tms_age_sex_dc[
            dataset == "tms_facs" & tissue == "Lung" &
            OR <= 5 & OR >= 0.2 & male_scd & female_scd
        ]$cell_ontology_final &
        RECEIVER_CELLTYPE %in% dt_md_tms_age_sex_dc[
            dataset == "tms_facs" & tissue == "Lung" &
            OR <= 5 & OR >= 0.2 & male_scd & female_scd
        ]$cell_ontology_final
    ]$dataset,
        dt_cci_full_diffsex[
        grepl("TMS FACS", dataset) & tissue == "Lung" &
        EMITTER_CELLTYPE %in% dt_md_tms_age_sex_dc[
            dataset == "tms_facs" & tissue == "Lung" &
            OR <= 5 & OR >= 0.2 & male_scd & female_scd
        ]$cell_ontology_final &
        RECEIVER_CELLTYPE %in% dt_md_tms_age_sex_dc[
            dataset == "tms_facs" & tissue == "Lung" &
            OR <= 5 & OR >= 0.2 & male_scd & female_scd
        ]$cell_ontology_final
    ]$REGULATION
)

unique(
    dt_cci_full_diffsex[grepl("TMS FACS", dataset) & tissue == "Lung"][, c("EMITTER_CELLTYPE", "dataset")]
)[, .N, by = c("dataset")]

unique(
    dt_cci_full[grepl("TMS FACS", dataset) & tissue == "Lung"][, c("EMITTER_CELLTYPE", "dataset")]
)[, .N, by = c("dataset")]

test <- table(
    dt_cci_full[
        grepl("TMS FACS", dataset) & tissue == "Lung"
    ]$dataset,
        dt_cci_full[
        grepl("TMS FACS", dataset) & tissue == "Lung"
    ]$REGULATION
)

test <- table(
    dt_cci_full[
        grepl("TMS FACS", dataset) & tissue == "Lung"
    ]$dataset,
        dt_cci_full[
        grepl("TMS FACS", dataset) & tissue == "Lung"
    ]$REGULATION
)

intersect(
    dt_cci_full[
        dataset == "TMS FACS (male)" & tissue == "Lung" &
        REGULATION == "UP"
    ]$CCI,
    dt_cci_full[
        dataset == "TMS FACS (female)" & tissue == "Lung" &
        REGULATION == "UP"
    ]$CCI
)

intersect(
    dt_cci_full[
        dataset == "TMS FACS (male)" & tissue == "Lung" &
        REGULATION == "DOWN"
    ]$CCI,
    dt_cci_full[
        dataset == "TMS FACS (female)" & tissue == "Lung" &
        REGULATION == "DOWN"
    ]$CCI
)

intersect(
    dt_cci_full[
        dataset == "TMS FACS (male)" & tissue == "Marrow" &
        REGULATION == "UP"
    ]$CCI,
    dt_cci_full[
        dataset == "TMS FACS (female)" & tissue == "Marrow" &
        REGULATION == "UP"
    ]$CCI
)

intersect(
    dt_cci_full[
        dataset == "TMS FACS (male)" & tissue == "Marrow" &
        REGULATION == "DOWN"
    ]$CCI,
    dt_cci_full[
        dataset == "TMS FACS (female)" & tissue == "Marrow" &
        REGULATION == "DOWN"
    ]$CCI
)

test <- lapply(
    unique(dt_cci_rel$tissue),
    function(tiss) {
        cci_male_down <- dt_cci_rel[
            dataset == "TMS FACS (male)" & tissue == tiss &
            REGULATION == "DOWN"
        ]$CCI
        cci_female_down <- dt_cci_rel[
            dataset == "TMS FACS (female)" & tissue == tiss &
            REGULATION == "DOWN"
        ]$CCI
        cci_male_up <- dt_cci_rel[
            dataset == "TMS FACS (male)" & tissue == tiss &
            REGULATION == "UP"
        ]$CCI
        cci_female_up <- dt_cci_rel[
            dataset == "TMS FACS (female)" & tissue == tiss &
            REGULATION == "UP"
        ]$CCI
        data.table(
            tissue = tiss,
            n_male_down = length(cci_male_down),
            n_male_up = length(cci_male_up),
            n_female_down = length(cci_female_down),
            n_female_up = length(cci_female_up),
            n_inter_down = length(intersect(cci_male_down, cci_female_down)),
            n_inter_up = length(intersect(cci_male_up, cci_female_up))
        )
    }
)
test <- rbindlist(test)
test
