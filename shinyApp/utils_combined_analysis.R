get_TCA_title <- function(
  input
) {
  renderUI(
    {
      req(input$TCA_DATASET_CHOICE)
      color_theme <- "color: rgb(20 120 206)"
      tags$p(
        "Global Analysis of the ",
        span(input$TCA_DATASET_CHOICE, style = color_theme)
      )
    }
  )
}

get_TCA_summary_type <- function(
  input
) {
  renderUI({
    if (input$TCA_SUMMARY_TYPE_CHOICE == "Table") {
      DT::dataTableOutput("TCA_SUMMARY_TABLE")
    } else if (input$TCA_SUMMARY_TYPE_CHOICE == "ORA Plot") {
      plotOutput("TCA_ORA_PLOT", height = "600px")
    }
  })
}

get_TCA_summary_intro <- function(
  input
) {
  renderUI(
    {
      req(input$TCA_DATASET_CHOICE)
      color_theme <- "color: rgb(20 120 206)"
      tags$div(
        tags$p(
          "Analysis of signals that are shared between several tissues."
        ),
        tags$p(
          "You can also see a plot of the most over-represented signals. Note that a term with a large ORA Score ",
          "does not necessarily take part in a lot of tissues (this information is given in the table)."
        ),
        tags$p(
          "We will clarify the names of the column in the tables. For now, here is what they mean:",
          "ID, global ORA Score Down, global ORA Score Up, Number of tissues in which the id is ",
          "down over-represented, same but for up over-represented, number of cell-cell interactions in which the ID ",
          "takes part, same but only in young cells, same but only in old cells."
        ),
        style = "font-size: 1.5em; text-align: justify;"
      )
    }
  )
}

get_TCA_summary_table <- function(
  input
) {
  DT::renderDT({
    req(input$TCA_DATASET_CHOICE, input$TCA_SUMMARY_TABLE_CHOICE, input$TCA_ORA_TYPE_CHOICE)
    categ <- input$TCA_SUMMARY_TABLE_CHOICE
    if (categ == "Tissues") {
      categ <- "ID"
      old_cols <- c("ID", "ORA_SCORE_DOWN", "ORA_SCORE_UP", "N_CELLTYPES", "N_CCI")
      new_cols <- c("Tissue", "ORA Score Down", "ORA Score Up", "Cell-Type N.", "CCI N.")
    } else if (categ == "Cell-Type Families") {
      categ <- "ER_CELL_FAMILY_GLOBAL"
      old_cols <- c("ER_CELL_FAMILY","ORA_SCORE_DOWN", "ORA_SCORE_UP", "N_LR")#, "N_LR_YOUNG", "N_LR_OLD")
      new_cols <- c("Cell-Type Family", "ORA Score Down", "ORA Score Up", "LR N.")#, "LR N. Young", "LR N. Old" )
    } else if (categ ==  "LR Genes") {
      categ <- "LR_GENES_GLOBAL" 
      old_cols <- c("LR_GENES", "ORA_SCORE_DOWN", "ORA_SCORE_UP", "N_ID", "N_ID_DOWN", "N_ID_UP", "N_ORA_ID_DOWN", "N_ORA_ID_UP",
                    "N_ER")#, "N_ER_YOUNG", "N_ER_OLD")
      new_cols <- c("LR Genes", "ORA Score Down", "ORA Score Up", "Tissue N.", "Tissue N. Down", "Tissue N. Up",
                    "ORA Tissue N. Down", "ORA Tissue N. Up", "ER N.")#, "ER N. Young", "ER N. Old")
    } else if (categ == "GO Terms") {
      categ <- "GO_TERMS_GLOBAL"
      old_cols <- c("GO_NAME", "ORA_SCORE_DOWN", "ORA_SCORE_UP", "N_ORA_ID_DOWN", "N_ORA_ID_UP", "N_ER")#, "N_ER_YOUNG", "N_ER_OLD")
      new_cols <- c("GO Term", "ORA Score Down", "ORA Score Up", "ORA Tissue N. Down",
                    "ORA Tissue N. Up", "CCI N.")#, "CCI N. Young", "CCI N. Old")
    } else if (categ == "KEGG Pathways") {
      categ <- "KEGG_PWS_GLOBAL"
      old_cols <- c("KEGG_NAME", "ORA_SCORE_DOWN", "ORA_SCORE_UP", "N_ORA_ID_DOWN", "N_ORA_ID_UP", "N_ER")#, "N_ER_YOUNG", "N_ER_OLD")
      new_cols <- c("GO Term", "ORA Score Down", "ORA Score Up", "ORA Tissue N. Down",
                    "ORA Tissue N. Up", "CCI N.")#, "CCI N. Young", "CCI N. Old")
    }
    dt <-  CCI_SUMMARY[[input$TCA_DATASET_CHOICE]][[categ]]
    dt <- dt[, old_cols, with = FALSE]
    setnames(dt, old_cols, new_cols)
    if (input$TCA_ORA_TYPE_CHOICE == "Up") {
      setorder(dt, -`ORA Score Up`, na.last = TRUE)
    } else if (input$TCA_ORA_TYPE_CHOICE == "Down") {
      setorder(dt, -`ORA Score Down`, na.last = TRUE)
    } 
    show_DT(
      dt,
      new_cols,
      new_cols[-1]
    )
  })
}

plot_TCA_ora <- function(
  input
) {
  renderPlot({
    req(input$TCA_DATASET_CHOICE, input$TCA_SUMMARY_TABLE_CHOICE, input$TCA_ORA_TYPE_CHOICE)
    categ <- input$TCA_SUMMARY_TABLE_CHOICE
    if (categ == "Tissues") {
      categ <- "ID"
    } else if (categ == "Cell-Type Families") {
      categ <- "ER_CELL_FAMILY"
    } else if (categ ==  "LR Genes") {
      categ <- "LR_GENES" 
    } else if (categ == "GO Terms") {
      categ <- "GO_TERMS"
    } else if (categ == "KEGG Pathways") {
      categ <- "KEGG_PWS"
    }
    obj <- DATASETS_COMBINED[[input$TCA_DATASET_CHOICE]]
    req(obj)
    if (input$TCA_ORA_TYPE_CHOICE == "Up") {
      reg <- "UP"
    } else if (input$TCA_ORA_TYPE_CHOICE == "Down") {
      reg <- "DOWN"
    } else if ( input$TCA_ORA_TYPE_CHOICE == "Stable") {
      reg <- "FLAT"
    }
    p <- PlotORA(
      object = obj,
      subID = NULL,
      category = categ,
      regulation = reg,
      max_terms_show = 30,
      global = TRUE,
      OR_threshold = 1,
      p_value_threshold = 0.05,
      stringent = FALSE
    )
    return(p)
  })
}


get_TCA_cci_intro <- function(
  input
) {
  renderUI(
    {
      req(input$TCA_DATASET_CHOICE)
      color_theme <- "color: rgb(20 120 206)"
      tags$div(
        tags$p(
          "Table with all interactions from the selected tissues, or Volcan Plot."
        ),
        style = "font-size: 1.5em; text-align: justify;"
      )
    }
  )
}

get_TCA_cci_details <- function(
  input
) {
  renderUI({
    if (input$TCA_CCI_DETAILS_CHOICE == "CCI Table") {
      DT::dataTableOutput("TCA_INTERACTION_TABLE")
    } else if (input$TCA_CCI_DETAILS_CHOICE == "Volcano Plot") {
      plotOutput("TCA_VOLCANO_PLOT", brush = "TCA_VOLCANO_brush", height = "600px")
    } else if (input$TCA_CCI_DETAILS_CHOICE == "Score Plot") {
      plotOutput("TCA_SCORES_PLOT", brush = "TCA_SCORES_brush", height = "600px")
    }
  })
}

choose_TCA_tissue <- function(
  input
) {
  renderUI({
    choices <- sort(unique(DATASETS_COMBINED[[input$TCA_DATASET_CHOICE]]@cci_detected$ID))
    pickerInput(
      inputId = "TCA_TISSUE_CHOICE",
      label = "Tissue",
      choices = choices,
      selected = choices,
      options = list(`actions-box` = TRUE),
      multiple = TRUE
    )
  })
}

choose_TCA_regulation <- function(
  input
) {
  renderUI({
    choices <- c("UP", "DOWN", "FLAT")
    pickerInput(
      inputId = "TCA_REGULATION_CHOICE",
      label = "Regulation",
      choices = choices,
      selected = choices,
      options = list(`actions-box` = TRUE),
      multiple = TRUE
    )
  })
}

get_TCA_cci_text <- function(
  input
) {
  renderUI({
    if (input$TCA_CCI_DETAILS_CHOICE == "CCI Table") {
      NULL
    } else if (input$TCA_CCI_DETAILS_CHOICE == "Volcano Plot") {
      verbatimTextOutput("TCA_VOLCANO_TEXTOUTPUT")
    } else if (input$TCA_CCI_DETAILS_CHOICE == "Score Plot") {
      verbatimTextOutput("TCA_SCORES_TEXTOUTPUT")
    }
  })
}

get_TCA_interaction_table <- function(
  input
) {
  DT::renderDT({
    req(input$TCA_DATASET_CHOICE, input$TCA_TISSUE_CHOICE, input$TCA_REGULATION_CHOICE)
    dt <-  DATASETS_COMBINED[[input$TCA_DATASET_CHOICE]]@cci_detected[
      ID %in% input$TCA_TISSUE_CHOICE & REGULATION %in% input$TCA_REGULATION_CHOICE
      ]
    req(dt)
    setorder(
      dt,
      -LOG2FC,
      `BH_P_VALUE_DE`
    )
    show_DT(
      dt,
      c("ID", "LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "LOG2FC", "BH_P_VALUE_DE",
        "REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD"),
      c("LOG2FC", "BH_P_VALUE_DE", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD"),
      "Table of all detected CCIs"
    )
  })
}

plot_TCA_VOLCANO <- function(
  input
) {
  renderPlot({
    req(input$TCA_DATASET_CHOICE, input$TCA_TISSUE_CHOICE, input$TCA_REGULATION_CHOICE)
    dt <-  DATASETS_COMBINED[[input$TCA_DATASET_CHOICE]]@cci_detected[
      ID %in% input$TCA_TISSUE_CHOICE & REGULATION %in% input$TCA_REGULATION_CHOICE
      ]
    req(dt)
    dt[, mlog10_pval := -log10(`BH_P_VALUE_DE` + 1E-4)]
    dt <- dt[, c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC", "mlog10_pval", "REGULATION")]
    show_volcano(dt)
  })
}

get_TCA_VOLCANO_text <- function(
  input
) {
  renderPrint({
    req(input$TCA_DATASET_CHOICE, input$TCA_TISSUE_CHOICE, input$TCA_REGULATION_CHOICE)
    dt <-  DATASETS_COMBINED[[input$TCA_DATASET_CHOICE]]@cci_detected[
      ID %in% input$TCA_TISSUE_CHOICE & REGULATION %in% input$TCA_REGULATION_CHOICE
      ]
    req(dt)
    dt[, mlog10_pval := -log10(`BH_P_VALUE_DE` + 1E-4)]
    dt <- dt[, c("ID", "LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC", "mlog10_pval", "REGULATION")]
    brushedPoints(dt, input$TCA_VOLCANO_brush, xvar = "LOG2FC", yvar = "mlog10_pval")
  })
}

plot_TCA_SCORES <- function(
  input
) {
  renderPlot({
    req(input$TCA_DATASET_CHOICE, input$TCA_TISSUE_CHOICE, input$TCA_REGULATION_CHOICE)
    dt <-  DATASETS_COMBINED[[input$TCA_DATASET_CHOICE]]@cci_detected[
      ID %in% input$TCA_TISSUE_CHOICE & REGULATION %in% input$TCA_REGULATION_CHOICE
      ]
    dt <- dt[, c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC", "REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD")]
    show_scores(dt)
  })
}

get_TCA_SCORES_text <- function(
  input
) {
  renderPrint({
    req(input$TCA_DATASET_CHOICE, input$TCA_TISSUE_CHOICE, input$TCA_REGULATION_CHOICE)
    dt <-  DATASETS_COMBINED[[input$TCA_DATASET_CHOICE]]@cci_detected[
      ID %in% input$TCA_TISSUE_CHOICE & REGULATION %in% input$TCA_REGULATION_CHOICE
      ]
    req(dt)
    dt <- dt[, c("ID", "LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC", "REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD")]
    brushedPoints(dt, input$TCA_SCORES_brush)
  })
}



