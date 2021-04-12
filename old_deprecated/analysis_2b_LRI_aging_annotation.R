
library(scDiffCom)
library(data.table)


LRI_mouse_curated <- copy(LRI_mouse$LRI_curated)

## Add cellAge annotation ####

hagr_cellage <- setDT(
  read.csv(
    "../data_scAgeCom/analysis/inputs_data/hagr_cellAge.csv",
    sep = ";",
    header = TRUE
  )
)
hagr_cellage[
  scDiffCom:::get_orthologs(
    hagr_cellage_human$gene_name,
    input_species = "human",
    one2one = FALSE
  ),
  on = c("gene_name==human_symbol"),
  mouse_gene_name := i.mouse_symbol
]
hagr_cellage$mouse_gene_name[is.na(hagr_cellage$mouse_gene_name)] <- "no_ortholog"

LRI_mouse_curated[
  hagr_cellage,
  on = c("LIGAND_1==mouse_gene_name"),
  LIGAND_1_CELLAGE := i.senescence_effect
]
LRI_mouse_curated[
  hagr_cellage,
  on = c("LIGAND_2==mouse_gene_name"),
  LIGAND_2_CELLAGE := i.senescence_effect
]
LRI_mouse_curated[
  hagr_cellage,
  on = c("RECEPTOR_1==mouse_gene_name"),
  RECEPTOR_1_CELLAGE := i.senescence_effect
]
LRI_mouse_curated[
  hagr_cellage,
  on = c("RECEPTOR_2==mouse_gene_name"),
  RECEPTOR_2_CELLAGE := i.senescence_effect
]
LRI_mouse_curated[
  hagr_cellage,
  on = c("RECEPTOR_3==mouse_gene_name"),
  RECEPTOR_3_CELLAGE := i.senescence_effect
]

LRI_mouse_curated[,
  IS_IN_CELLAGE := ifelse(
    LIGAND_1 %in% hagr_cellage$mouse_gene_name |
      LIGAND_2 %in% hagr_cellage$mouse_gene_name |
      RECEPTOR_1 %in% hagr_cellage$mouse_gene_name |
      RECEPTOR_2 %in% hagr_cellage$mouse_gene_name |
      RECEPTOR_3 %in% hagr_cellage$mouse_gene_name,
    TRUE,
    FALSE
  )
]

# LRI_mouse_curated[, .N, by = IS_IN_CELLAGE]

## Add GenAge annotation ####

hagr_genage <- setDT(
  read.csv(
    "../data_scAgeCom/analysis/inputs_data/hagr_genage_models.csv",
    sep = ",",
    header = TRUE
  )
)
hagr_genage <- hagr_genage[organism == "Mus musculus"]
hagr_genage[, EFFECT := paste(
  lifespan.effect,
  longevity.influence,
  sep = ":"
)]

LRI_mouse_curated[
  hagr_genage,
  on = c("LIGAND_1==symbol"),
  LIGAND_1_GENAGE := i.EFFECT
]
LRI_mouse_curated[
  hagr_genage,
  on = c("LIGAND_2==symbol"),
  LIGAND_2_GENAGE := i.EFFECT
]
LRI_mouse_curated[
  hagr_genage,
  on = c("RECEPTOR_1==symbol"),
  RECEPTOR_1_GENAGE := i.EFFECT
]
LRI_mouse_curated[
  hagr_genage,
  on = c("RECEPTOR_2==symbol"),
  RECEPTOR_2_GENAGE := i.EFFECT
]
LRI_mouse_curated[
  hagr_genage,
  on = c("RECEPTOR_3==symbol"),
  RECEPTOR_3_GENAGE := i.EFFECT
]

LRI_mouse_curated[,
                  IS_IN_GENAGE := ifelse(
                    LIGAND_1 %in% hagr_genage$symbol |
                      LIGAND_2 %in% hagr_genage$symbol |
                      RECEPTOR_1 %in% hagr_genage$symbol |
                      RECEPTOR_2 %in% hagr_genage$symbol |
                      RECEPTOR_3 %in% hagr_genage$symbol,
                    TRUE,
                    FALSE
                  )
]

LRI_mouse_curated[, .N, by = c("IS_IN_GENAGE", "IS_IN_CELLAGE")]


## Add SASP/SEN annotation ####

senesence_markers <- setDT(
  read.csv(
    "../data_scAgeCom/analysis/inputs_data/senescence_markers.csv",
    sep = ",",
    header = TRUE
  )
)

LRI_mouse_curated[
  senesence_markers,
  on = c("LIGAND_1==gene"),
  LIGAND_1_SENESCENCE := i.category
]
LRI_mouse_curated[
  senesence_markers,
  on = c("LIGAND_2==gene"),
  LIGAND_2_SENESCENCE := i.category
]
LRI_mouse_curated[
  senesence_markers,
  on = c("RECEPTOR_1==gene"),
  RECEPTOR_1_SENESCENCE := i.category
]
LRI_mouse_curated[
  senesence_markers,
  on = c("RECEPTOR_2==gene"),
  RECEPTOR_2_SENESCENCE := i.category
]
LRI_mouse_curated[
  senesence_markers,
  on = c("RECEPTOR_3==gene"),
  RECEPTOR_3_SENESCENCE := i.category
]
