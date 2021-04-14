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

## Layout config ####

shiny_layout <- list(
  color_theme = "color: rgb(20 120 206)",
  style_intro_title = paste(
    "width: 60%;",
    "margin:auto;",
    "font-size: 26px;",
    "text-align: center;"
  ),
  style_intro_text = paste(
    "width: 60%;",
    "margin:auto;",
    "font-size: 18px;",
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
  ),
  intro_method_title = paste(
    "Overview of the analysis"
  )
)

shiny_html_content <- list(
  main_title = shiny_text$main_title,
  intro_title = tags$p(
    div(
      style = shiny_layout$style_intro_title,
      shiny_text$intro_title
    )
  ),
  intro_overview = tags$div(
    style = shiny_layout$style_intro_text,
    tags$p(
      paste(
        "This website offers the opportunity to explore how "
      ),
      tags$b(
        "intercellular communication changes with age in 23 mouse tissues",
        .noWS = "outside"
      ),
      paste(
        ". It includes both tissue specific results and patterns shared accross",
        "several organs."
      )
    ),
    tags$p(
      "The full methodology behind this analysis is described in our",
      "(manuscript in preparation). The most important steps are highlighted",
      "below and we also recommend to look at the Glossary section",
      "of the site to make the best use of it."
    )
  ),
  intro_method_title = tags$p(
    div(
      style = shiny_layout$style_intro_title,
      shiny_text$intro_method_title
    )
  ),
  intro_code_title = tags$h3(
    "Methodology and Code",
    style = shiny_layout$color_theme
  ),
  intro_code_text = tags$div(
    style = shiny_layout$style_intro_text,
    tags$p(
      paste(
        "In order to investigate how intercellular communication varies",
        "between two biological conditions of interest, we built scDiffCom.",
        "This R package is available on"
      ),
      tags$a(
        href = "https://github.com/CyrilLagger/scDiffCom",
        "GitHub",
        target = "_blank"
      ),
      paste(
        "(and soon on CRAN) and a tutorial is provided"
      ),
      tags$a(
        href = "https://cyrillagger.github.io/scDiffCom/articles/scDiffCom-vignette.html",
        "here",
        target = "_blank"
      ),
      "to explain how to apply it to your own data.",
      tags$p(
        paste(
          "In summary, scDiffCom uses as inputs a user-defined scRNA-seq",
          "dataset (in the form of a "
        ),
        tags$a(
          href = "https://satijalab.org/seurat/index.html",
          "Seurat object",
          target = "_blank",
          .noWS = "outside"
        ),
        paste(
          ") and a internal collection of curated ligand-receptor interactions.",
          "It then considers all potential signalling patterns that may be",
          "mediated by these ligands and receptors between the cell types",
          "of the dataset. A series of statistical tests is performed",
          "to first retain only those interactions that are likely to",
          "correspond to biological signals and second to determine how they",
          "change between two conditions given on the cells. Finally, an",
          "over-representation test allows the package to infer the dominant",
          "changing patterns at either a gene-centric, cell type-centric or",
          "annotation-centric (GO/KEGG) level."
        ),
        .noWS = c("after-begin", "before-end")
      )
    ),
    tags$p(
      paste(
        "Here, we have specifically used scDiffCom to study differences in",
        "intercellular communication between young and old murine cells (see",
        "the data description below). The scripts used for both the analyis",
        "and this website are available on GitHub:"
      ),
      tags$a(
        href = "https://github.com/CyrilLagger/scAgeCom",
        "scAgeCom",
        target = "_blank"
      ),
      "and ",
      tags$a(
        href = "https://github.com/CyrilLagger/scAgeComShiny",
        "scAgeComShiny.",
        target = "_blank",
        .noWS = "outside"
      ),
      "."
    ),
  ),
  intro_scrna_title =  tags$h3(
    "Aging Single Cell Datasets",
    style = shiny_layout$color_theme
  ),
  intro_scrna_text = tags$div(
    style = shiny_layout$style_intro_text,
    tags$p(
      paste(
        "This analyis was made possible thanks to aging scRNA-seq datasets",
        "from two previously published studies:"
      )
    ),
    tags$ol(
      tags$li(
        "Tabula Muris Senis (TMS):",
        tags$a(
          href = "https://tabula-muris-senis.ds.czbiohub.org/",
          "webpage",
          target = "_blank"
        ),
        " and ",
        tags$a(
          href = "https://www.nature.com/articles/s41586-020-2496-1",
          "Nature article",
          target = "_blank",
          .noWS = "outside"
        ),
        "."
      ),
      tags$li(
        "Calico2019:",
        tags$a(
          href = "https://mca.research.calicolabs.com/",
          "webpage",
          target = "_blank"
        ),
        " and ",
        tags$a(
          href = "https://genome.cshlp.org/content/29/12/2088",
          "Genome Research article",
          target = "_blank",
          .noWS = "outside"
        ),
        "."
      )
    ),
    tags$p(
      "Overall, this represents data for 23 organs (give here a summary picture)"
    )
  ),
  intro_lri_title = tags$h3(
    "Collection of Ligand-Receptor Interactions",
    style = shiny_layout$color_theme
  ),
  intro_lri_text =   tags$div(
    style = shiny_layout$style_intro_text,
    tags$p(
      "scDiffCom internally relies on a list of ligand-receptor interactions",
      "that we have retrieved and processed from eight previous studies:"
    ),
    tags$ol(
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
          "bioRxiv article.",
          target = "_blank",
          .noWS = "outside"
        ),
        "."
      )
    ),
    tags$p(
      paste(
        "Note that when combining these databases, we carefully took into",
        "account the interactions involving heteromeric complexes (see our",
        "manuscript in preparation)."
      )
    ),
    tags$p(
      paste(
        "We give below the mouse interactions used as input for the",
        "scAgeCom analysis. This table and its human equivalent can also",
        "directly be accessed from scDiffcom in R."
      )
    )
  ),
  intro_contact_title = tags$h3(
    "Contact",
    style = shiny_layout$color_theme
  ),
  intro_contact_text = tags$div(
    style = shiny_layout$style_intro_text,
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
    shiny_html_content = shiny_html_content
  )
)

saveRDS(
  scAgeCom_shiny_data,
  "../data_scAgeCom/analysis/outputs_data/scAgeCom_shiny_data.rds"
)

