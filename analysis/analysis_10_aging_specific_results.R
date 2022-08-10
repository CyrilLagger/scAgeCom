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

## Annotate ORA with HAGR and PMID ####

dt_ora_lri <- dt_ora_rel[
    ORA_CATEGORY == "LRI" & ORA_REGULATION != "Not Over-represented",
    c("dataset", "tissue", "VALUE", "ORA_REGULATION")
]
dt_ora_lri[
    dt_lri_mouse,
    on = "VALUE==LRI",
    (pmid_colnames) := mget(paste0("i.", pmid_colnames))
]
dt_ora_lri[
    dt_lri_mouse,
    on = "VALUE==LRI",
    (hagr_colnames) := mget(paste0("i.", hagr_colnames))
]

dt_ora_ligand <- dt_ora_rel[
    ORA_CATEGORY == "LIGAND_COMPLEX" &
     ORA_REGULATION != "Not Over-represented",
    c("dataset", "tissue", "VALUE", "ORA_REGULATION")
]
dt_ora_ligand[
    ,
    L1 := unlist(lapply(
        strsplit(VALUE, "_"),
        function(i) {
            return(i[[1]])
        }
    ))
]
dt_ora_ligand[
    ,
    L2 := unlist(lapply(
        strsplit(VALUE, "_"),
        function(i) {
            if (length(i) >= 2) {
                return(i[[2]])
            } else {
                return(NA)
            }
        }
    ))
]
dt_ora_ligand[
    dt_lri_mouse,
    on = "L1==LIGAND_1",
    L1_N_agepmid := i.L1_N_agepmid
]
dt_ora_ligand[
    dt_lri_mouse,
    on = "L2==LIGAND_2",
    L2_N_agepmid := i.L2_N_agepmid
]
dt_ora_ligand[
    dt_lri_mouse,
    on = "L1==LIGAND_1",
    HAGR_L1 := i.HAGR_L1
]
dt_ora_ligand[
    dt_lri_mouse,
    on = "L2==LIGAND_2",
    HAGR_L2 := i.HAGR_L2
]

dt_ora_receptor <- dt_ora_rel[
    ORA_CATEGORY == "RECEPTOR_COMPLEX" &
     ORA_REGULATION != "Not Over-represented",
    c("dataset", "tissue", "VALUE", "ORA_REGULATION")
]
dt_ora_receptor[
    ,
    R1 := unlist(lapply(
        strsplit(VALUE, "_"),
        function(i) {
            return(i[[1]])
        }
    ))
]
dt_ora_receptor[
    ,
    R2 := unlist(lapply(
        strsplit(VALUE, "_"),
        function(i) {
            if (length(i) >= 2) {
                return(i[[2]])
            } else {
                return(NA)
            }
        }
    ))
]
dt_ora_receptor[
    ,
    R3 := unlist(lapply(
        strsplit(VALUE, "_"),
        function(i) {
            if (length(i) >= 3) {
                return(i[[3]])
            } else {
                return(NA)
            }
        }
    ))
]
dt_ora_receptor[
    dt_lri_mouse,
    on = "R1==RECEPTOR_1",
    R1_N_agepmid := i.R1_N_agepmid
]
dt_ora_receptor[
    dt_lri_mouse,
    on = "R2==RECEPTOR_2",
    R2_N_agepmid := i.R2_N_agepmid
]
dt_ora_receptor[
    dt_lri_mouse,
    on = "R3==RECEPTOR_3",
    R3_N_agepmid := i.R3_N_agepmid
]
dt_ora_receptor[
    dt_lri_mouse,
    on = "R1==RECEPTOR_1",
    HAGR_R1 := i.HAGR_R1
]
dt_ora_receptor[
    dt_lri_mouse,
    on = "R2==RECEPTOR_2",
    HAGR_R2 := i.HAGR_R2
]
dt_ora_receptor[
    dt_lri_mouse,
    on = "R3==RECEPTOR_3",
    HAGR_R3 := i.HAGR_R3
]

## Top Ligand and receptors vs aging resources ####

dt_ora_ligand_counts_res <- copy(
    shiny_dt_ora_key_counts[
        ORA_CATEGORY == "Ligand"
    ]
)
dt_ora_ligand_counts_res[
    dt_ora_ligand,
    on = "VALUE",
    (pmid_colnames[1:2]) := mget(paste0("i.", pmid_colnames[1:2]))
]
dt_ora_ligand_counts_res[
    ,
    top_LRI := sapply(
        VALUE,
        function(i) {
            temp_lri <- dt_lri_mouse[LIGAND_1 == i | LIGAND_2 == i]$LRI
            paste(head(unique(shiny_dt_ora_key_counts[
                ORA_CATEGORY == "Ligand-Receptor Interaction"
            ][VALUE %in% temp_lri]$VALUE)), collapse = ",")
        }
    )
]

dt_ora_receptor_counts_res <- copy(
    shiny_dt_ora_key_counts[
        ORA_CATEGORY == "Receptor"
    ]
)
dt_ora_receptor_counts_res[
    dt_ora_receptor,
    on = "VALUE",
    (pmid_colnames[3:5]) := mget(paste0("i.", pmid_colnames[3:5]))
]
dt_ora_receptor_counts_res[
    ,
    top_LRI := sapply(
        VALUE,
        function(i) {
            temp_lri <- dt_lri_mouse[
                RECEPTOR_1 == i | RECEPTOR_2 == i | RECEPTOR_3 == i
            ]$LRI
            paste(head(unique(shiny_dt_ora_key_counts[
                ORA_CATEGORY == "Ligand-Receptor Interaction"
            ][VALUE %in% temp_lri]$VALUE)), collapse = ",")
        }
    )
]

dt_ora_lr_counts_res <- rbindlist(
    list(
        dt_ora_ligand_counts_res,
        dt_ora_receptor_counts_res
    ),
    fill = TRUE
)

fwrite(
    dt_ora_lr_counts_res,
    paste0(
        path_scagecom_output,
        "temp_dt_ora_lr_counts_res.csv"
    )
)

## Inflammation/immune ####

# select GO terms of interest in more than 10 tissues
# from shiny_dt_ora_key_counts

aging_go_immune_topm10 <- c(
    "immune response", "immune system process",
    "lymphocyte differentiation", "T cell differentiation",
    "regulation of immune response",
    "regulation of viral life cycle", "regulation of viral process",
    "leukocyte differentiation",
    "positive regulation of adaptive immune response",
    "regulation of immune system process",
    "cytokine-mediated signaling pathway",
    "T cell activation",
    "chemokine-mediated signaling pathway",
    "defense response",
    "leukocyte chemotaxis",
    "lymphocyte migration",
    "modulation by symbiont of entry into host",
    "positive regulation of immune effector process",
    "positive regulation of immune response",
    "regulation of lymphocyte apoptotic process",
    "regulation of viral entry into host cell",
    "regulation of adaptive immune response",
    "regulation of leukocyte activation",
    "T cell migration",
    "leukocyte migration",
    "negative regulation of B cell apoptotic process",
    "negative regulation of lymphocyte apoptotic process",
    "negative regulation of mature B cell apoptotic process",
    "regulation of B cell activation",
    "regulation of B cell apoptotic process",
    "regulation of chemokine (C-X-C motif) ligand 2 production",
    "regulation of cytokine production involved in immune response",
    "regulation of leukocyte apoptotic process",
    "regulation of leukocyte migration",
    "regulation of lymphocyte activation",
    "regulation of mature B cell apoptotic process"
)
table(
    aging_go_immune_topm10 %in% dt_ora_key_counts$VALUE
)
aging_go_immune_topm10 <- dt_ora_key_counts[
    ORA_CATEGORY == "GO_TERMS" &
    VALUE %in% aging_go_immune_topm10 &
    ORA_REGULATION == "UP" &
    `Overall (Union)` >= 10
]

# select KEGG pws of interest in more than 10 tissues

aging_kegg_immune_topm10 <- c(
    "Human immunodeficiency virus 1 infection",
    "Epstein-Barr virus infection",
    "Human T-cell leukemia virus 1 infection",
    "Antigen processing and presentation",
    "Chemokine signaling pathway",
    "Viral protein interaction with cytokine and cytokine receptor",
    "Staphylococcus aureus infection",
    "Cytokine-cytokine receptor interaction",
    "NF-kappa B signaling pathway",
    "Natural killer cell mediated cytotoxicity"
)
table(
    aging_kegg_immune_topm10 %in% dt_ora_key_counts$VALUE
)
aging_kegg_immune_topm10 <- dt_ora_key_counts[
    ORA_CATEGORY == "KEGG_PWS" &
    VALUE %in% aging_kegg_immune_topm10 &
    ORA_REGULATION == "UP" &
    `Overall (Union)` >= 10
]

# ER_celltype overrepresentation

aging_er_immune <- dt_ora_key_counts[
    ORA_CATEGORY == "ER_CELLFAMILIES" &
    ORA_REGULATION == "UP" &
    grepl("leukocyte", VALUE) &
    `Overall (Union)` >= 10
]

#select LRI based on go/kegg terms 

aging_dt_lri_fom_go_immune <- dt_ora_key_counts[
    ORA_CATEGORY == "LRI" &
    VALUE %in% scDiffCom::LRI_mouse$LRI_curated_GO[
        GO_ID %in% scDiffCom::gene_ontology_level[
            NAME %in% aging_go_immune_topm10$VALUE
        ]$ID
    ]$LRI &
    ORA_REGULATION == "UP"
]
aging_dt_lri_fom_kegg_immune <- dt_ora_key_counts[
    ORA_CATEGORY == "LRI" &
    VALUE %in% scDiffCom::LRI_mouse$LRI_curated_KEGG[
        KEGG_NAME %in% aging_kegg_immune_topm10
    ]$LRI
]

aging_lri_from_gokegg_immune <- c(
    "B2m:Cd3g", "B2m:Cd3d", "H2-D1:Cd8b1",
    "H2-K1:Cd8b1", "B2m:Cd247", "Mif:Cd74",
    "H2-Q6:Cd8b1", "Ccl5:Ccrl2", "H2-K1:Cd8a",
    "H2-T23:Cd8b1", "B2m:Klrd1", "Ccl11:Cxcr3",
    "H2-D1:Cd8a", "Ccl5:Ccr1", "Ccl5:Ccr5", "Hmgb1:Thbd",
    "Tnfsf12:Tnfrsf12a", "Ccl11:Cxcr3",
    #"H2-M3:Cd8b1",
    #"H2-T22:Cd8b1",
    "Slpi:Plscr1",
    "Hmgb1:Ager"
)
table(
    aging_lri_from_gokegg_immune %in% dt_ora_key_counts$VALUE
)
aging_dt_lri_immune <- dt_ora_key_counts[
    ORA_CATEGORY == "LRI" &
    ORA_REGULATION == "UP" &
    VALUE %in% aging_lri_from_gokegg_immune
]

# add aging information
aging_dt_lri_immune <- merge.data.table(
  aging_dt_lri_immune,
  dt_lri_mouse_val_clean[, c("LRI", "pubmed", "HAGR"), with = FALSE],
  by.x = "VALUE",
  by.y = "LRI",
  all.x = TRUE,
  all.y = FALSE,
  sort = FALSE
)
#careful H2-XX (HLA in human) should be more represented

# add validation results
aging_dt_lri_immune <- merge.data.table(
  aging_dt_lri_immune,
  dt_lri_mouse_val_clean[, c(
    "LRI",
    "summary_val"),
    with = FALSE],
  by.x = "VALUE",
  by.y = "LRI",
  all.x = TRUE,
  all.y = FALSE,
  sort = FALSE
)

# full table
aging_dt_immune <- rbindlist(
    list(
        aging_dt_lri_immune,
        aging_go_immune_topm10,
        aging_kegg_immune_topm10,
        aging_er_immune
    ),
    fill = TRUE,
    use.names = TRUE
)

# add sex results

aging_dt_immune_sex <- merge.data.table(
    aging_dt_immune,
    dt_ora_key_counts_diffsex[
        ORA_REGULATION %in% c("UP", "DOWN")
    ],
    by = c("ORA_CATEGORY", "VALUE"),
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE,
    suffixes = c("_age", "_sex")
)

#most interesting new finding is Slpi:Plscr1

pmid_aging_lri[["Slpi"]]
pmid_aging_lri[["Plscr1"]]

dt_ora_key_summary[VALUE == "Slpi:Plscr1" & ORA_REGULATION == "UP"]
table(dt_cci_rel[LRI == "Slpi:Plscr1"]$REGULATION)

table(dt_cci_age[
    LRI == "Slpi:Plscr1" & REGULATION == "UP"
]$dataset_tissue)

dt_cci_age[
    LRI == "Slpi:Plscr1" & REGULATION == "UP" &
    grepl("arrow", tissue)
]$CCI

sort(table(dt_cci_age[
    LRI == "Slpi:Plscr1" & REGULATION == "UP"
]$EMITTER_CELLTYPE))

sort(table(dt_cci_age[
    LRI == "Slpi:Plscr1" & REGULATION == "UP"
]$RECEIVER_CELLTYPE))

sort(table(dt_cci_age[
    LRI == "Slpi:Plscr1" & REGULATION == "UP"
]$ER_CELLTYPE))

## Lipid metabolism ####

# select GO terms of interest in more than 10 tissues

aging_go_lipmed_candidates <- unique(
        dt_ora_key_counts[
        ORA_CATEGORY == "GO_TERMS" &
        (grepl("fat", VALUE) | grepl("lip", VALUE) |
         grepl("chole", VALUE))
    ]$VALUE
)
aging_go_lipmet_topm10 <- dt_ora_key_counts[
    ORA_CATEGORY == "GO_TERMS" &
    VALUE %in% aging_go_lipmed_candidates &
    ORA_REGULATION %in% c("UP", "DOWN") &
    `Overall (Union)` >= 10
]

# select KEGG pws of interest in more than 10 tissues #

aging_kegg_lipmed_candidates <- unique(
        dt_ora_key_counts[
        ORA_CATEGORY == "KEGG_PWS" &
        (grepl("fat", VALUE) | grepl("lip", VALUE) |
         grepl("Chole", VALUE))
    ]$VALUE
)
aging_kegg_lipmet_topm10 <- dt_ora_key_counts[
    ORA_CATEGORY == "KEGG_PWS" &
    VALUE %in% aging_kegg_lipmed_candidates &
    ORA_REGULATION %in% c("UP", "DOWN") &
    `Overall (Union)` >= 10
]

#select LRI based on go/kegg terms

aging_dt_lri_fom_go_lipmet <- dt_ora_key_counts[
    ORA_CATEGORY == "LRI" &
    VALUE %in% scDiffCom::LRI_mouse$LRI_curated_GO[
        GO_ID %in% scDiffCom::gene_ontology_level[
            NAME %in% aging_go_lipmet_topm10$VALUE
        ]$ID
    ]$LRI &
    ORA_REGULATION %in% c("UP", "DOWN")
]
aging_dt_lri_fom_kegg_lipmet <- dt_ora_key_counts[
    ORA_CATEGORY == "LRI" &
    VALUE %in% scDiffCom::LRI_mouse$LRI_curated_KEGG[
        KEGG_NAME %in% aging_kegg_lipmet_topm10$VALUE
    ]$LRI
]
aging_lri_from_gokegg_lipmed <- c(
    "Apoe:Sdc4", "Apoe:Lrp1", "App:Lrp10", "Apoe:Ldlr", "Apoa1:Ldlr",
    "F13a1:Itgb1", "App:Cav1", "App:Ncstn", "App:Notch2",
    "App:Sorl1", "Mmp2:Itgb1", "Psen1:Ncstn",
    "C3:Lrp1", "Psap:Lrp1", "App:Vldlr", "Anxa1:Fpr2",
    "Rps27a:Ldlr", "Tgm2:Adgrg1", "App:Notch1",
    "App:Slc45a3"
)
table(
    aging_lri_from_gokegg_lipmed %in% dt_ora_key_counts$VALUE
)
aging_dt_lri_lipmed <- dt_ora_key_counts[
    ORA_CATEGORY == "LRI" &
    ORA_REGULATION %in% c("UP", "DOWN") &
    VALUE %in% aging_lri_from_gokegg_lipmed
]

#select Ligands based on go/kegg terms
aging_dt_ligand_fom_go_lipmet <- dt_ora_key_counts[
    ORA_CATEGORY == "LIGAND_COMPLEX" &
    VALUE %in% unlist(strsplit(scDiffCom::LRI_mouse$LRI_curated_GO[
        GO_ID %in% scDiffCom::gene_ontology_level[
            NAME %in% aging_go_lipmet_topm10$VALUE
        ]$ID
    ]$LRI, ":")) &
    ORA_REGULATION %in% c("UP", "DOWN")
]
aging_dt_ligand_fom_kegg_lipmet <- dt_ora_key_counts[
    ORA_CATEGORY == "LIGAND_COMPLEX" &
    VALUE %in% unlist(strsplit(scDiffCom::LRI_mouse$LRI_curated_KEGG[
        KEGG_NAME %in% aging_kegg_lipmet_topm10$VALUE
    ]$LRI, ":"))
]

# add aging information
aging_dt_lri_lipmed <- merge.data.table(
  aging_dt_lri_lipmed,
  dt_lri_mouse_val_clean[, c("LRI", "pubmed", "HAGR"), with = FALSE],
  by.x = "VALUE",
  by.y = "LRI",
  all.x = TRUE,
  all.y = FALSE,
  sort = FALSE
)

# add validation results
aging_dt_lri_lipmed <- merge.data.table(
  aging_dt_lri_lipmed,
  dt_lri_mouse_val_clean[, c(
    "LRI",
    "summary_val"),
    with = FALSE],
  by.x = "VALUE",
  by.y = "LRI",
  all.x = TRUE,
  all.y = FALSE,
  sort = FALSE
)

# full table
aging_dt_lipmed <- rbindlist(
    list(
        aging_dt_lri_lipmed,
        aging_go_lipmet_topm10,
        aging_kegg_lipmet_topm10
    ),
    fill = TRUE,
    use.names = TRUE
)

# add sex results

aging_dt_lipmed_sex <- merge.data.table(
    aging_dt_lipmed,
    dt_ora_key_counts_diffsex[
        ORA_REGULATION %in% c("UP", "DOWN")
    ],
    by = c("ORA_CATEGORY", "VALUE"),
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE,
    suffixes = c("_age", "_sex")
)

# App:Lrp10 sex-dependent

dt_ora_key_summary[
    VALUE == "App:Lrp10" &
    grepl("TMS", dataset)
]

test <- merge.data.table(
    dcast.data.table(
        dt_ora_key_summary[
            VALUE == "App:Lrp10" &
            grepl("TMS", dataset)
        ],
        tissue ~ dataset,
        value.var = "ORA_REGULATION"
    ),
    dcast.data.table(
        dt_ora_key_summary_diffsex[
            VALUE == "App:Lrp10" &
            grepl("TMS", dataset) &
            !grepl("combined", dataset)
        ],
        tissue ~ dataset,
        value.var = "ORA_REGULATION"
    ),
    by = "tissue",
    all = TRUE
)

dt_ora_key_summary_diffsex[
    VALUE == "App:Lrp10" &
    grepl("TMS", dataset) &
    !grepl("combined", dataset)
]

dt_ora_key_summary_diffsex[
    VALUE == "App:Lrp10" &
    ORA_REGULATION %in% c("UP", "DOWN")
]

dt_ora_key_summary[
    VALUE == "App" &
    dataset_tissue %in% c(
        "TMS FACS (male)_Brain",
        "TMS FACS (female)_Brain"
    )
]

dt_ora_key_summary_diffsex[
    VALUE == "App" &
    dataset_tissue %in% c(
        "TMS FACS (young)_Brain",
        "TMS FACS (old)_Brain"
    )
]

dt_ora_key_summary[
    VALUE == "Apoe" &
    dataset_tissue %in% c(
        "TMS FACS (male)_Brain",
        "TMS FACS (female)_Brain"
    )
]

dt_ora_key_summary_diffsex[
    VALUE == "Apoe" &
    dataset_tissue %in% c(
        "TMS FACS (young)_Brain",
        "TMS FACS (old)_Brain"
    )
]

#Apoe:Sdc4 is interesting based on longevity and lipid


# interest in the following up LRI for (not really lipid, move below)
# cancer, angiongenesis, etc
# L linked to vascualar aging, R not linked to aging much
dt_lri_mouse[LRI == "Tgm2:Adgrg1"]
# actually Tgm2:Itgb1 is down in 5 tissues, same or not?
# other low Agrn:Itgb1, Ptdss1:Scarb1, Ptprm:Ptprm, Adam10:Tspan14
# App:Slc45a3
dt_lri_mouse[LRI == "Ptprm:Ptprm"]
test <- dt_cci_rel[
    LRI == "Ptprm:Ptprm" & REGULATION == "DOWN"
]

# compare to aging ref

aging_lri_from_gokegg_lipmed <- dt_lri_mouse[
    LRI %in% aging_lri_from_gokegg_lipmed
]

## ECM ####

# select GO terms of interest in more than 10 tissues
# from shiny_dt_ora_key_counts

aging_go_ecm_topm10 <- c(
    "basement membrane",
    "collagen-containing extracellular matrix",
    "extracellular matrix",
    "extracellular matrix organization",
    "extracellular structure organization",
    "biological adhesion",
    "cell adhesion",
    "cell-substrate adhesion",
    "collagen binding",
    "cell junction",
    "cell junction organization",
    "cell-matrix adhesion",
    "cell-substrate junction",
    "extracellular region",
    "basement membrane organization",
    "cell junction assembly",
    "focal adhesion"
)
table(
    aging_go_ecm_topm10 %in% dt_ora_key_counts$VALUE
)
aging_go_ecm_topm10 <- dt_ora_key_counts[
    ORA_CATEGORY == "GO_TERMS" &
    VALUE %in% aging_go_ecm_topm10 &
    ORA_REGULATION == "DOWN" &
    `Overall (Union)` >= 10
]

# select KEGG pws of interest in more than 10 tissues

aging_kegg_ecm_topm10 <- c(
    "ECM-receptor interaction",
    "Focal adhesion"
)
table(
    aging_kegg_ecm_topm10 %in% dt_ora_key_counts$VALUE
)
aging_kegg_ecm_topm10 <- dt_ora_key_counts[
    ORA_CATEGORY == "KEGG_PWS" &
    VALUE %in% aging_kegg_ecm_topm10 &
    ORA_REGULATION == "DOWN" &
    `Overall (Union)` >= 10
]

# ER_celltype overrepresentation

aging_er_ecm <- dt_ora_key_counts[
    ORA_CATEGORY == "ER_CELLFAMILIES" &
    ORA_REGULATION == "DOWN" &
    grepl("connective", VALUE) &
    `Overall (Union)` >= 7
]

#select LRI based on go/kegg terms

aging_dt_lri_fom_go_ecm <- dt_ora_key_counts[
    ORA_CATEGORY == "LRI" &
    VALUE %in% scDiffCom::LRI_mouse$LRI_curated_GO[
        GO_ID %in% scDiffCom::gene_ontology_level[
            NAME %in% aging_go_ecm_topm10$VALUE
        ]$ID
    ]$LRI &
    ORA_REGULATION == "DOWN"
]
aging_dt_lri_fom_kegg_ecm <- dt_ora_key_counts[
    ORA_CATEGORY == "LRI" &
    VALUE %in% scDiffCom::LRI_mouse$LRI_curated_KEGG[
        KEGG_NAME %in% aging_kegg_ecm_topm10$VALUE
    ]$LRI &
    ORA_REGULATION == "DOWN"
]

aging_lri_from_gokegg_ecm <- c(
    "Sell:Cd34", "Adam17:Itgb1", "Adam9:Itgb1",
    "Nid1:Itgb1", "Adam15:Itga9", "Adam15:Itgb1",
    "Hspg2:Itgb1", "Mmp2:Itgb1", "Postn:Itgb1",
    "Col4a2:Itgb1_Itga1", "Adam9:Itgb5", "Alcam:Nrp1",
    "Col15a1:Itgb1_Itga1",
    #"Fam3c:Lifr",
    "Fbln1:Itgb1",
    "Gcg:Dpp4", "Tgm2:Itgb1", "Adam9:Itgav", "Ccn4:Itgb1",
    "Cdh1:Cdh1", "Cdh1:Itga1_Itgb1", "Ceacam1:Ceacam1",
    "Col1a1:Sdc4", "Col1a2:Sdc4", "Col3a1:Ddr2",
    "Col4a1:Itga9_Itgb1", "Col4a1:Itgb1_Itga1",
    "Col4a2:Itga9_Itgb1", "Col4a2:Sdc4",
    "Col5a2:Itgb1_Itga1", "Col6a1:Cd44",
    "Col6a2:Itgb1_Itga1", "Col4a2:Itgb5",
    "Col1a1:Cd44", "Col4a1:Cd44", "Col1a2:Cd36",
    "Col1a2:Cd44"
)
table(
    aging_lri_from_gokegg_ecm %in% dt_ora_key_counts$VALUE
)
aging_dt_lri_ecm <- dt_ora_key_counts[
    ORA_CATEGORY == "LRI" &
    ORA_REGULATION == "DOWN" &
    VALUE %in% aging_lri_from_gokegg_ecm
]

# add aging information
aging_dt_lri_ecm <- merge.data.table(
  aging_dt_lri_ecm,
  dt_lri_mouse_val_clean[, c("LRI", "pubmed", "HAGR"), with = FALSE],
  by.x = "VALUE",
  by.y = "LRI",
  all.x = TRUE,
  all.y = FALSE,
  sort = FALSE
)

# add validation results
aging_dt_lri_ecm <- merge.data.table(
  aging_dt_lri_ecm,
  dt_lri_mouse_val_clean[, c(
    "LRI",
    "summary_val"),
    with = FALSE],
  by.x = "VALUE",
  by.y = "LRI",
  all.x = TRUE,
  all.y = FALSE,
  sort = FALSE
)

# full table
aging_dt_ecm <- rbindlist(
    list(
        aging_dt_lri_ecm,
        aging_go_ecm_topm10,
        aging_kegg_ecm_topm10,
        aging_er_ecm
    ),
    fill = TRUE,
    use.names = TRUE
)

# add sex results

aging_dt_ecm_sex <- merge.data.table(
    aging_dt_ecm,
    dt_ora_key_counts_diffsex[
        ORA_REGULATION %in% c("UP", "DOWN")
    ],
    by = c("ORA_CATEGORY", "VALUE"),
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE,
    suffixes = c("_age", "_sex")
)

#most interesting new finding is

dt_lri_mouse[LRI == "Alcam:Nrp1"]

###############################################
## New LRI of interest ####

ORA_LIGAND_N[
    ORA_REGULATION != "FLAT" &
    LIGAND_1_agepmid <= 5 &
    HAGR_L1 == FALSE &
    N >= 5
][1:20]
ORA_RECEPTOR_N[
    ORA_REGULATION != "FLAT" &
    RECEPTOR_1_agepmid <= 5 &
    HAGR_R1 == FALSE &
    N >= 5
][1:20]

## Egfl7 #####

ORA_LRI_N[
    grepl("Egfl", VALUE)
]

# Egfl7:Notch1/2/3 upregulated.
# Egfl7 not aging characterized
# Notch is linked to pubmed aging (thymocyte development, angiogenesis
# down implies less regeneration in muscle, senescence, etc, LRI_mouse_pmid_aging$Notch1). Notch
# is related to many processes.
# Here: angiogenesis (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3286203/) and
# (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3031397/)
# neural stem cell renewal (https://www.nature.com/articles/ncb1896)
# Egfl7 acts as a Notch antagonist (reduce target gene expression)/opposite as well
#/depent on the micorenviornment
# Overexpression of Egfl7 can lead to abnormal vessel patterning
# In NSC, Egfl7 reduces proliferation and trigger differentiation,
# here maybe a similar effect on MSC in adipose tissues

ccis_egfl7_up <- data_4_tissue_specific_results$CCI_table[
    LIGAND_1 == "Egfl7" &
    `Age Regulation` == "UP"
]

#appears up-regulaed in adipose tissues, brain, heart and mammary gland a lot
ccis_egfl7_up[, .N, by = c("Dataset", "Tissue")][order(-N)]
#secreted by endothelial cells
ccis_egfl7_up[, .N, by = c("Tissue", "Emitter Cell Type")][order(-N)]
#towards MSC in adipose, endocardial cells and fibroblast of cardiac tissues in heart, T cells
ccis_egfl7_up[, .N, by = c("Tissue", "Receiver Cell Type")][order(-N)][1:10]
#towards pericyte and oligodendrocyte precursor cells in Brain
ccis_egfl7_up[, .N, by = c("Tissue", "Receiver Cell Type")][order(-N)][11:20]

ccis_egfl7_up[, .N,
 by = c("Dataset", "Tissue", "Emitter Cell Type", "Receiver Cell Type")
][order(-N)][1:10]

## F13a1 ####

ORA_LRI_N[
    grepl("F13a1", VALUE)
]

#downregulated with integrins (https://pubmed.ncbi.nlm.nih.gov/10816592/)
#The interaction of α4β1 and α9β1 with FXIII could be biologically significant.
# Both α4β1 and α9β1 are highly expressed on specific populations of leukocytes,
# where they play a prominent role in transendothelial leukocyte migration.
# As a member of the coagulation cascade, FXIII is likely to be enriched at sites
# of vascular injury and inflammation, where its interaction with α4β1 and α9β1 
#could promote leukocyte extravasation.
LRI_mouse_dt[LIGAND_1 == "F13a1"]

#not well documented with aging (AMD papers)
LRI_combined_pmid_aging$F13a1

ccis_f13a1_down <- data_4_tissue_specific_results$CCI_table[
    LIGAND_1 == "F13a1" &
    `Age Regulation` == "DOWN"
]

#appears down-regulaed in Marrow, adipose tissue and other
ccis_f13a1_down[, .N, by = c("Dataset", "Tissue")][order(-N)]
#secreted by (pro-)monocyte/macrophages/myeloid principally
ccis_f13a1_down[, .N, by = c("Tissue", "Emitter Cell Type")][order(-N)]
#towards immune cell types
ccis_f13a1_down[, .N, by = c("Tissue", "Receiver Cell Type")][order(-N)]

#not well documented interaction with integrins

## Gpi1 ####

ORA_LRI_N[
    grepl("Gpi1", VALUE)
]

LRI_combined_pmid_aging$Gpi1
LRI_mouse_dt[LIGAND_1 == "Gpi1"]

## Tln1 (probably false positive as intracellular with Itgb1) ####

## Col6a2, H2-Q6 # to do ####

## Slpi Plscr1 ####
ORA_LRI_N[
    grepl("Slpi", VALUE)
]

#Slpi with Plscr1 not linked to aging

LRI_mouse_dt[LRI == "Slpi:Pls"]

ccis_slpi_up <- data_4_tissue_specific_results$CCI_table[
    LIGAND_1 == "Slpi" &
    `Age Regulation` == "UP"
]

## Rpsa
ORA_LRI_N[
    grepl("Rpsa", VALUE)
]

ccis_rpsa_up <- data_4_tissue_specific_results$CCI_table[
    RECEPTOR_1 == "Rpsa" &
    `Age Regulation` == "UP"
]

## Gpc1

## Tnfrsf21