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
## collect and prepare all results for shiny
##
####################################################
##

## Libraries ####

library(scDiffCom)
library(data.table)
library(htmltools)
library(bsplus)

## glossary (bsplus accordion) ####

shiny_bsplus_glossary <- bs_accordion(
  id = "HELP_GLOSSARY"
) %>%
  bs_append(
    title = "LRI: Ligand-Receptor Interaction",
    content = paste(
      "TODO"
    )
  ) %>%
  bs_append(
    title = "CCI: Cell-Cell Interaction",
    content = paste(
      "TODO"
    )
  ) %>%
  bs_append(
    title = "ORA Score",
    content = paste(
      "TODO"
    )
  ) %>%
  bs_append(
    title = "TODO",
    content = paste(
      "TODO"
    )
  )

## Layout config ####

shiny_layout <- list(
  color_theme = "color: rgb(20 120 206)",
  style_main_title = paste(
    "font-size: 26px;"
  ),
  style_intro_title = paste(
    "margin: 20px auto;",
    "font-size: 26px;",
    "text-align: center;"
  ),
  style_accordion_titles = paste(
    "font-size: 20px;",
    "margin-top: 0;",
    "margin-bottom: 0;"
  ),
  style_intro_text = paste(
    "margin:auto;",
    "font-size: 15px;",
    "text-align: justify;"
  ),
  style_intro_text_accordion = paste(
    "margin:auto;",
    "font-size: 15px;",
    "text-align: justify;"
  )
)

## Shiny html ####

shiny_text <- list(
  main_title = paste(
    "A Murine ATLAS of Age-related Changes in Intercellular Communication",
    " (Development Website!)."
  ),
  intro_title = paste(
    "Welcome to scAgeCom!"
  )
)

shiny_html_content <- list(
  main_title = tags$span(
    div(
      style = shiny_layout$style_main_title,
      shiny_text$main_title
    )
  ),
  intro_title = tags$p(
    div(
      style = shiny_layout$style_intro_title,
      shiny_text$intro_title
    )
  ),
  intro_overview = tags$div(
    style = shiny_layout$style_intro_text,
    tags$p(
      tags$b(
        paste(
          "Explore how intercellular communication changes with age",
          "in 23 mouse tissues."
        )
      )
    ),
    tags$p(
      paste(
        "This website includes both tissue specific results and",
        "a global comparison of the changes shared accross",
        "the different organs.",
        "The full methodology behind this analysis is described in our",
        "(manuscript in preparation) but the most important steps",
        "are highligthed below."
      )
    )
  ),
  intro_code_title = tags$h3(
    "Methodology and Code",
    style = shiny_layout$style_accordion_titles
  ),
  intro_code_text = tags$div(
    style = shiny_layout$style_intro_text_accordion,
    tags$h4("The package"),
    tags$p(
      "We built the R package",
      tags$a(
        href = "https://github.com/CyrilLagger/scDiffCom",
        "scDiffCom",
        target = "_blank"
      ),
      paste(
        "to inverstigate how intercellular communication varies between",
        "two biological conditions of interest (young/old, sick/healthy, etc.).",
        "Here, we have used scDiffCom to study differences in intercellular",
        "comunication between young and old murine cells. Please visit"
      ),
      tags$a(
        href = "https://cyrillagger.github.io/scDiffCom/articles/scDiffCom-vignette.html",
        "this tutorial",
        target = "_blank"
      ),
      "to apply this package to your own data."
    ),
    tags$h4("How does it work?"),
    tags$p(
      paste(
        "scDiffCom can be applied to any scRNA-seq dataset",
        "(defined as a "
      ),
      tags$a(
        href = "https://satijalab.org/seurat/index.html",
        "Seurat object",
        target = "_blank",
        .noWS = "outside"
      ),
      paste(
        "). It uses its own integrated database of ligand-receptor interactions",
        "to list all the potential signals between the cell types of the dataset."
      )
    ),
    tags$p(
      paste(
        "Statistical tests are then performed to only retain biologically",
        "significant interactions and to measure their change between",
        "the two conditions. Finally, an over-representation test allows the",
        "package to infer the dominant",
        "changing patterns on a gene-centric, cell type-centric or",
        "annotation-centric (GO/KEGG) level."
      ),
      .noWS = c("after-begin", "before-end")
    ),
    tags$h4("Find our scripts on GitHub"),
    tags$ol(
      tags$li(
        "Package for differential analysis on any dataset:",
        tags$a(
          href = "https://github.com/CyrilLagger/scDiffCom",
          "scDiffCom",
          target = "_blank"
        )
      ),
      tags$li(
        "Code for this aging analysis:",
        tags$a(
          href = "https://github.com/CyrilLagger/scAgeCom",
          "scAgeCom",
          target = "_blank"
        )
      ),
      tags$li(
        "Code for this website:",
        tags$a(
          href = "https://github.com/CyrilLagger/scAgeComShiny",
          "scAgeComShiny",
          target = "_blank"
        )
      )
    ),
  ),
  intro_scrna_title =  tags$h3(
    "Aging Single Cell Datasets",
    style = shiny_layout$style_accordion_titles
  ),
  intro_scrna_text = tags$div(
    style = shiny_layout$style_intro_text_accordion,
    tags$p(
      paste(
        "We based our analysis on several murine scRNA-seq datasets",
        "provided by the two studies below:"
      )
    ),
    tags$ul(
      tags$li(
        "Tabula Muris Senis (TMS):",
        tags$ul(
          tags$li(
            tags$a(
              href = "https://tabula-muris-senis.ds.czbiohub.org/",
              "webpage",
              target = "_blank"
            )
          ),
          tags$li(
            tags$a(
              href = "https://www.nature.com/articles/s41586-020-2496-1",
              "Nature article",
              target = "_blank",
              .noWS = "outside"
            )
          )
        )
      ),
      tags$li(
        "Calico2019:",
        tags$ul(
          tags$li(
            tags$a(
              href = "https://mca.research.calicolabs.com/",
              "webpage",
              target = "_blank"
            )
          ),
          tags$li(
            tags$a(
              href = "https://genome.cshlp.org/content/29/12/2088",
              "Genome Research article",
              target = "_blank",
              .noWS = "outside"
            )
          )
        )
      )
    ),
    tags$p(
      "Overall, this represents data for 23 organs (give here a summary picture)"
    )
  ),
  intro_lri_title = tags$h3(
    "Collection of Ligand-Receptor Interactions",
    style = shiny_layout$style_accordion_titles
  ),
  intro_lri_text = tags$div(
    style = shiny_layout$style_intro_text_accordion,
    tags$p(
      "scDiffCom internally relies on a list of ligand-receptor interactions",
      "that we have processed and combined from the eight studies in the list",
      "opposite. We carefully took into account the interactions involving",
      "heteromeric complexes (see our manuscript in preparation)."
    ),
    tags$p(
      "You will find below the mouse interactions for the scAgeCom analysis.",
      "This table and its human equivalent can also be directly accessed",
      "in the scDiffCom package."
    )
  ),
  intro_lri_db_list = tags$ul(
    tags$li(
      "CellChat:",
      tags$a(
        href = "http://www.cellchat.org/", "webpage",
        target = "_blank"
      ),
      " and ",
      tags$a(
        href = "https://www.nature.com/articles/s41467-021-21246-9",
        "Nature Communications article",
        target = "_blank",
        .noWS = "outside"
      ),
      "."
    ),
    tags$li(
      "CellPhoneDB:",
      tags$a(
        href = "https://www.cellphonedb.org/",
        "webpage",
        target = "_blank"
      ),
      " and ",
      tags$a(
        href = "https://www.nature.com/articles/s41596-020-0292-x",
        "Nature Protocol article",
        target = "_blank",
        .noWS = "outside"
      ),
      "."
    ),
    tags$li(
      "CellTalkDB:",
      tags$a(
        href = "http://tcm.zju.edu.cn/celltalkdb/",
        "webpage",
        target = "_blank"
      ),
      " and ",
      tags$a(
        href = "https://academic.oup.com/bib/advance-article-abstract/doi/10.1093/bib/bbaa269/5955941",
        "Briefings in Bioinformatics article",
        target = "_blank",
        .noWS = "outside"
      )
    ),
    tags$li(
      "connectomeDB2020:",
      tags$a(
        href = "https://github.com/forrest-lab/NATMI",
        "webpage",
        target = "_blank"
      ),
      " and ",
      tags$a(
        href = "https://www.nature.com/articles/s41467-020-18873-z",
        "Nature Communications article",
        target = "_blank",
        .noWS = "outside"
      ),
      "."
    ),
    tags$li(
      "ICELLNET:",
      tags$a(
        href = "https://github.com/soumelis-lab/ICELLNET",
        "webpage",
        target = "_blank"
      ),
      " and ",
      tags$a(
        href = "https://www.nature.com/articles/s41467-021-21244-x",
        "Nature Communications article",
        target = "_blank",
        .noWS = "outside"
      ),
      "."
    ),
    tags$li(
      "NicheNet:",
      tags$a(
        href = "https://github.com/saeyslab/nichenetr",
        "webpage",
        target = "_blank"
      ),
      " and ",
      tags$a(
        href = "https://www.nature.com/articles/s41592-019-0667-5",
        "Nature Methods article",
        target = "_blank",
        .noWS = "outside"
      ),
      "."
    ),
    tags$li(
      "SingleCellSignalR:",
      tags$a(
        href = "http://www.bioconductor.org/packages/release/bioc/html/SingleCellSignalR.html",
        "webpage",
        target = "_blank"
      ),
      " and ",
      tags$a(
        href = "https://academic.oup.com/nar/article/48/10/e55/5810485",
        "Nucleic Acids Research article",
        target = "_blank",
        .noWS = "outside"
      ),
      "."
    ),
    tags$li(
      "scTensor:",
      tags$a(
        href = "https://github.com/rikenbit/scTensor",
        "webpage",
        target = "_blank"
      ),
      " and ",
      tags$a(
        href = "https://www.biorxiv.org/content/10.1101/566182v1",
        "bioRxiv article",
        target = "_blank",
        .noWS = "outside"
      ),
      "."
    )
  ),
  intro_contact_title = tags$h3(
    "Contact",
    style = shiny_layout$style_accordion_titles
  ),
  intro_contact_text = tags$div(
    style = shiny_layout$style_intro_text_accordion,
    tags$p(
      "Please contact ..."
    )
  )
)

## save all results ####

scAgeCom_shiny_data <- c(
  readRDS(
    "../data_scAgeCom/analysis/outputs_data/data_2_LRI_data_preparation.rds"
  ),
  readRDS(
    "../data_scAgeCom/analysis/outputs_data/data_4_tissue_specific_results.rds"
  ),
  readRDS(
    "../data_scAgeCom/analysis/outputs_data/data_5_tissue_shared_results.rds"
  ),
  list(
    shiny_bsplus_glossary = shiny_bsplus_glossary,
    shiny_html_content = shiny_html_content
  )
)

saveRDS(
  scAgeCom_shiny_data,
  "../data_scAgeCom/analysis/outputs_data/scAgeCom_shiny_data.rds"
)

