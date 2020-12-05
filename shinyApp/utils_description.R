get_description_html <- function(
  input
) {
  renderUI(
    {
      color_theme <- "color: rgb(20 120 206)"
      tags$div(
        style = "width: 60%; margin:auto; font-size: 18px; text-align: justify;",
        tags$h2(
          "Overview of the analysis", 
          style = color_theme
        ),
        tags$p(
          "This project provides a comprehensive investigation of age-related changes in mouse intercellular communication.
          It combines scRNA-seq data and curated ligand-receptor interactions with a novel analysis technique that
          allows to statistically infere differentially expressed cell-cell interaction patterns."
        ),
        tags$p(
          "Link to the future paper:..."
        ),
        tags$h3(
          "Single-cell Datasets",
          style = color_theme
        ),
        tags$p(
          "We have leveraged transcritomics single-data from two previous studies:"
        ),
        tags$ol(
          tags$li(
            "Tabula Muris Senis (TMS):",
            tags$a(href= "https://tabula-muris-senis.ds.czbiohub.org/", "see there webpage", target="_blank"),
            " and ",
            tags$a(href= "https://www.nature.com/articles/s41586-020-2496-1", "Nature article.", target="_blank")
            ),
          tags$li(
            "Calico 2019:",
            tags$a(href= "https://mca.research.calicolabs.com/", "see there webpage", target="_blank"),
            " and ",
            tags$a(href= "https://genome.cshlp.org/content/29/12/2088", "Genome Research article.", target="_blank")
            )
        ),
        tags$h3(
          "Ligand-receptor Databases",
          style = "color: rgb(20 120 206);"
        ),
        tags$p(
          "We have collected and compared curated ligand-receptor interactions from 6 previous analysis:"
        ),
        tags$ol(
          tags$li(
            "CellPhoneDB (CPDB):",
            tags$a(href= "https://www.cellphonedb.org/", "see there webpage", target="_blank"),
            " and ",
            tags$a(href= "https://www.nature.com/articles/s41596-020-0292-x", "Nature Protocol article.", target="_blank")
          ),
          tags$li(
            "SingleCellSignalR (SCSR):",
            tags$a(href= "http://www.bioconductor.org/packages/release/bioc/html/SingleCellSignalR.html", "see there webpage", target="_blank"),
            " and ",
            tags$a(href= "https://academic.oup.com/nar/article/48/10/e55/5810485", "Nucleic Acids Research article.", target="_blank")
          ),
          tags$li(
            "NicheNet:",
            tags$a(href= "https://github.com/saeyslab/nichenetr", "see there webpage", target="_blank"),
            " and ",
            tags$a(href= "https://www.nature.com/articles/s41592-019-0667-5", "Nature Methods article.", target="_blank")
          ),
          tags$li(
            "ICELLNET:",
            tags$a(href= "https://github.com/soumelis-lab/ICELLNET", "see there webpage", target="_blank"),
            " and ",
            tags$a(href= "https://www.biorxiv.org/content/10.1101/2020.03.05.976878v1", "bioRxiv article.", target="_blank")
          ),
          tags$li(
            "CellChat:",
            tags$a(href= "http://www.cellchat.org/", "see there webpage", target="_blank"),
            " and ",
            tags$a(href= "https://www.biorxiv.org/content/10.1101/2020.07.21.214387v1", "bioRxiv article.", target="_blank")
          ),
          tags$li(
            "scTensor:",
            tags$a(href= "https://github.com/rikenbit/scTensor", "see there webpage", target="_blank"),
            " and ",
            tags$a(href= "https://www.biorxiv.org/content/10.1101/566182v1", "bioRxiv article.", target="_blank")
          )
        ),
        tags$h3(
          "Code - scDiffCom", 
          style = color_theme
        ),
        tags$p(
          "Inspired by several analysis techniques from the aferomentionned tools, we have build an R package,
          scDiffCom, that allows to statistically assess the differential expression of cell-cell interactions between two conditions of interest."
        ),
        tags$p(
          "The package is currently available on GitHub",
          tags$a(href= "https://github.com/CyrilLagger/scDiffCom", "here.", target="_blank"),
          "It is not intended to only be used for the aging analysis presented here, but can in theory be applied to any scRNA-seq datasets."
        ),
        tags$p(
          "A vignette will follow soon..."
        ),
        tags$p(
          "The code specifically related to this aging analysis and the ShinyApp are also available on GitHub",
          tags$a(href= "https://github.com/CyrilLagger/scAgeCom", "here.", target="_blank")
        ),
        tags$h3(
          "Team and Acknowledgement", 
          style = color_theme
        )
      )
    }
  )
}