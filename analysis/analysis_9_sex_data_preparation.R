####################################################
##
## Project: scAgeCom
##
## lagger.cyril@gmail.com
## ursu_eugen@hotmail.com
## anais.equey@gmail.com
##
## Combine sex and age analysis
##
####################################################
##

## prepare Calico cell number distribution ####

dt_metadata_calico <- dt_metadata[
  grepl("calico", dataset)
]
dt_metadata_calico[
  ,
  age_group := ifelse(
    age %in% c("young"), "YOUNG", "OLD"
  )
]
dt_metadata_calico[
  ,
  tissue := ifelse(
    grepl("kidney", dataset),
    "Kidney",
    ifelse(
      grepl("lung", dataset),
      "Lung",
      "Spleen"
    )
  )
]
dt_metadata_calico[, dataset := "calico_droplet"]
dt_metadata_calico[, sex := "male"]
dt_md_tms_age_sex_calico <- dt_metadata_calico[
  ,
  .N,
  by = c("dataset", "tissue", "cell_ontology_final", "sex", "age_group")
][order(dataset, tissue, sex, age_group)]

## Distribution of cells across age vs sex for each tissue ####

dt_md_tms_age_sex_full <- rbindlist(
  l = list(
    dt_md_tms_age_sex,
    dt_md_tms_age_sex_calico
  )
)

dt_md_tms_age_sex_dc <- dcast.data.table(
  dt_md_tms_age_sex_full,
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
    ifelse(
      dataset == "tms_droplet",
      paste("TMS Droplet (male)", tissue, cell_ontology_final, sep = "_"),
      paste("Calico Droplet (male)", tissue, cell_ontology_final, sep = "_")
    )
  )
]
dt_md_tms_age_sex_dc[
  ,
  dataset_female := ifelse(
    dataset == "tms_facs",
    paste("TMS FACS (female)", tissue, cell_ontology_final, sep = "_"),
    ifelse(
      dataset == "tms_droplet",
      paste("TMS Droplet (female)", tissue, cell_ontology_final, sep = "_"),
      paste("Calico Droplet (female)", tissue, cell_ontology_final, sep = "_")
    )
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
        "Supplementary_DATA_sex_distribution2.csv"
    )
)

dt_md_tms_age_sex_dc <- fread(
  paste0(
    path_scagecom_output,
    "Supplementary_DATA_sex_distribution2.csv"
  )
)

## Extract CCIs detected in both sex and age analysis ####

dt_cci_age[, tissue_cci := paste(tissue, CCI, sep = "_")]
dt_cci_sex[, tissue_cci := paste(tissue, CCI, sep = "_")]

dt_cci_sexage_facs <- merge.data.table(
    dcast.data.table(
        dt_cci_age[
            grepl("TMS FACS", dataset) &
            tissue_cci %in% dt_cci_sex[grepl("TMS FACS", dataset)]$tissue_cci,
            c("dataset", "LRI", "tissue", "tissue_cci", "REGULATION")
        ],
        formula = tissue_cci + tissue + LRI ~ dataset,
        value.var = "REGULATION",
        fill = "NOT_DETECTED"
    ),
    dcast.data.table(
        dt_cci_sex[
            grepl("TMS FACS", dataset) &
            tissue_cci %in% dt_cci_age[grepl("TMS FACS", dataset)]$tissue_cci,
            c("dataset", "LRI", "tissue", "tissue_cci", "REGULATION")
        ],
        formula = tissue_cci + tissue + LRI ~ dataset,
        value.var = "REGULATION",
        fill = "NOT_DETECTED"
    ),
    by = c("tissue_cci", "LRI", "tissue"),
    all = TRUE
)

dt_cci_sexage_droplet <- merge.data.table(
    dcast.data.table(
        dt_cci_age[
            grepl("TMS Droplet", dataset) &
            tissue_cci %in%
            dt_cci_sex[grepl("TMS Droplet", dataset)]$tissue_cci,
            c("dataset", "LRI", "tissue", "tissue_cci", "REGULATION")
        ],
        formula = tissue_cci + tissue + LRI ~ dataset,
        value.var = "REGULATION",
        fill = "NOT_DETECTED"
    ),
    dcast.data.table(
        dt_cci_sex[
            grepl("TMS Droplet", dataset) &
            tissue_cci %in%
            dt_cci_age[grepl("TMS Droplet", dataset)]$tissue_cci,
            c("dataset", "LRI", "tissue", "tissue_cci", "REGULATION")
        ],
        formula = tissue_cci + tissue + LRI ~ dataset,
        value.var = "REGULATION",
        fill = "NOT_DETECTED"
    ),
    by = c("tissue_cci", "LRI", "tissue"),
    all = TRUE
)

dt_cci_sexage_facs[
    ,
    .N,
    by = c(
        "TMS FACS (female)",
        "TMS FACS (male)",
        "TMS FACS (young)",
        "TMS FACS (old)"
    )
][order(-N)][1:20]

dt_cci_sexage_facs[
   `TMS FACS (young)` %in% c("NSC", "FLAT") &
   `TMS FACS (male)` == "UP" &
   `TMS FACS (female)` %in% c("DOWN", "FLAT", "NCS")
]$LRI

sort(table(dt_cci_sexage_facs[
   `TMS FACS (young)` %in% c("NSC", "FLAT") &
   `TMS FACS (male)` == "UP" &
   `TMS FACS (female)` %in% c("DOWN", "FLAT")
]$LRI))

## Consider for example lung results ####

dt_cci_sexage_facs_lung <- dt_cci_sexage_facs[tissue == "Lung"]
