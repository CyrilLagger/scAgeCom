
## TCA overall ####

choose_TCA_dataset <- function(
  input
) {
  renderUI({
    pickerInput(
      inputId = "TCA_DATASET_CHOICE",
      choices = names(DATASETS_COMBINED),
      options = list(
        `actions-box` = TRUE
      ),
      multiple = FALSE,
      inline = FALSE
    )
  })
}

get_TCA_title <- function(
  input
) {
  renderUI(
    {
      tags$p(
        div(style="display: inline-block;", "Please choose a dataset: "),
        div(style="display: inline-block;margin-top: 25px;", uiOutput("TCA_DATASET_CHOICE", inline = TRUE))
      )
    }
  )
}

get_TCA_overview_table <- function(
  input
) {
  DT::renderDT({
    req(input$TCA_DATASET_CHOICE)
    dt <-  copy(CCI_SUMMARY[[input$TCA_DATASET_CHOICE]][["ID"]])
    req(dt)
    dt[, N_TISSUES := 1]
    cols_to_keep <- c(
      "N_TISSUES", "N_CELLTYPES", "N_CCI", "N_CCI_FLAT", "N_CCI_DOWN", "N_CCI_UP", 
      "N_CCI_NON_SIGNIFICANT_CHANGE"
    )
    dt <- dt[, lapply(.SD, sum), .SDcols = cols_to_keep]
    new_col_names <- c(
      "Total tissues", "Total cell-types", "Total CCI", "Flat CCI",
      "Down CCI", "Up CCI", "Non-significant CCI"
    )
    setnames(dt, cols_to_keep, new_col_names)
    show_DT(
      data = dt,
      cols_to_show = new_col_names,
      cols_numeric = NULL,
      table_title = NULL,
      options = list(dom = "t"),
      callback = htmlwidgets::JS(
        "var tips = [
      'Number of tissues in the dataset of interest',
      'Number of cell-types in all tissues',
      'Number of cell-cell interactions detected in all tissues',
      'Number of stable CCIs with age',
      'Number of down-regulated CCIs with age',
      'Number of up-regulated CCIs with age',
      'Number of CCIs changing non-significantly with age'
      ],
      header = table.columns().header();
      for (var i = 0; i < tips.length; i++) {
      $(header[i]).attr('title', tips[i]);
      }"),
      rownames = FALSE
    )
  })
}

## TCA CCIs sidebar ####

download_TCA_table <- function(
  input
) {
  dt <- reactive(
    {
      dt <-  DATASETS_COMBINED[[input$TCA_DATASET_CHOICE]]@cci_detected
      setorder(
        dt,
        -LOG2FC_ALL,
        `BH_P_VALUE_DE`
      )
      cols_to_keep <- c("ID", "LRI", "Emitter Cell Type", "Receiver Cell Type", "LOG2FC", "Adj. p-value",
                        "Age-Regulation", "Score Young", "Score Old")
      setnames(
        dt,
        old = c("ID", "LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "LOG2FC_ALL", "BH_P_VALUE_DE",
                "REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD"),
        new = cols_to_keep
      )
      dt[, cols_to_keep, with = FALSE]
    }
  )
  downloadHandler(
    filename = function() {
      paste0(
        "cci_table_",
        tolower(gsub(" ", "_", input$TCA_DATASET_CHOICE, fixed = TRUE)),
        ".csv"
      )
    },
    content = function(file) {
      fwrite(dt(), file)
    }
  )
}

choose_TCA_tissue <- function(
  input
) {
  renderUI({
    req(input$TCA_DATASET_CHOICE)
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

choose_TCA_lri <- function(
  input 
) {
  renderUI({
    req(input$TCA_DATASET_CHOICE)
    ALL_LRI_LABEL = 'All LRI'
    choices <-
      c(ALL_LRI_LABEL,
        sort(unique(DATASETS_COMBINED[[input$TCA_DATASET_CHOICE]]@cci_detected[["LR_GENES"]]))
      )
    selectizeInput(
      inputId = 'TCA_LRI_CHOICE',
      label = 'Filter by Ligand-Receptor Interactions',
      choices = choices,
      selected = ALL_LRI_LABEL,
      multiple = TRUE,
      options = list(allowEmptyOption = TRUE,
                     placeholder = 'Type LRIs')
    )
  })
}

## TCA CCIs mainpanel ####

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
    } else if (input$TCA_CCI_DETAILS_CHOICE == "LRI-FC Plot") {
      plotOutput("TCA_LRIFC_PLOT", brush = "TCA_LRIFC_brush", height = "600px")
    }
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
    } else if (input$TCA_CCI_DETAILS_CHOICE == "LRI-FC Plot") {
      verbatimTextOutput("TCA_LRIFC_TEXTOUTPUT")
    }
  })
}

get_TCA_interaction_table <- function(
  input
) {
  DT::renderDT({
    req(input$TCA_DATASET_CHOICE, input$TCA_TISSUE_CHOICE)#, input$TCA_REGULATION_CHOICE)
    dt <-  copy(DATASETS_COMBINED[[input$TCA_DATASET_CHOICE]]@cci_detected[
      ID %in% input$TCA_TISSUE_CHOICE])# & REGULATION %in% input$TCA_REGULATION_CHOICE]
    req(dt)
    if (!('All LRI' %in% input$TCA_LRI_CHOICE)) {
      dt <- dt[`LR_GENES` %in% input$TCA_LRI_CHOICE]
    }
    setorder(
      dt,
      -LOG2FC_ALL,
      `BH_P_VALUE_DE`
    )
    setnames(
      dt,
      old = c("ID", "LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "LOG2FC_ALL", "BH_P_VALUE_DE",
              "REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD"),
      new = c("Tissue", "LRI", "Emitter Cell Type", "Receiver Cell Type", "LOG2FC", "Adj. p-value",
              "Age-Regulation", "Score Young", "Score Old")
    )
    show_DT(
      dt,
      c("Tissue", "LRI", "Emitter Cell Type", "Receiver Cell Type", "LOG2FC", "Adj. p-value",
        "Age-Regulation", "Score Young", "Score Old"),
      c("LOG2FC", "Adj. p-value", "Score Young", "Score Old"),
      "Cell-Cell Interation Table"
    )
  })
}

plot_TCA_VOLCANO <- function(
  input
) {
  renderPlot({
    req(input$TCA_DATASET_CHOICE, input$TCA_TISSUE_CHOICE)#, input$TCA_REGULATION_CHOICE)
    dt <-  DATASETS_COMBINED[[input$TCA_DATASET_CHOICE]]@cci_detected[
      ID %in% input$TCA_TISSUE_CHOICE]# & REGULATION %in% input$TCA_REGULATION_CHOICE]
    req(dt)
    xlims <- c(min(dt[["LOG2FC_ALL"]]), max(dt[["LOG2FC_ALL"]]))
    ylims <- c(0, 4)
    if (!('All LRI' %in% input$TCA_LRI_CHOICE)) {
      dt <- dt[`LR_GENES` %in% input$TCA_LRI_CHOICE]
    }
    dt[, minus_log10_pval := -log10(`BH_P_VALUE_DE` + 1E-4)]
    setnames(
      dt,
      old = c("LOG2FC_ALL", "REGULATION"),
      new = c("LOG2FC", "Age Regulation")
    )
    dt <- dt[, c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC", "minus_log10_pval", "Age Regulation")]
    show_volcano(
      data = dt,
      xlims = xlims,
      ylims = ylims
    )
  })
}

get_TCA_VOLCANO_text <- function(
  input
) {
  renderPrint({
    req(input$TCA_DATASET_CHOICE, input$TCA_TISSUE_CHOICE)#, input$TCA_REGULATION_CHOICE)
    dt <-  DATASETS_COMBINED[[input$TCA_DATASET_CHOICE]]@cci_detected[
      ID %in% input$TCA_TISSUE_CHOICE]# & REGULATION %in% input$TCA_REGULATION_CHOICE]
    req(dt)
    if (!('All LRI' %in% input$TCA_LRI_CHOICE)) {
      dt <- dt[`LR_GENES` %in% input$TCA_LRI_CHOICE]
    }
    dt[, minus_log10_pval := -log10(`BH_P_VALUE_DE` + 1E-4)]
    dt <- dt[, c("ID", "LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "REGULATION", "LOG2FC_ALL", "minus_log10_pval")]
    setnames(
      dt, 
      old = c("ID", "LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "REGULATION", "LOG2FC_ALL"),
      new = c("Tissue", "LRI", "Emitter Cell Type", "Receiver Cell Type", "Age Regulation" ,"LOG2FC")
    )
    brushedPoints(dt, input$TCA_VOLCANO_brush, xvar = "LOG2FC", yvar = "minus_log10_pval")
  })
}

plot_TCA_SCORES <- function(
  input
) {
  renderPlot({
    req(input$TCA_DATASET_CHOICE, input$TCA_TISSUE_CHOICE)#, input$TCA_REGULATION_CHOICE)
    dt <-  DATASETS_COMBINED[[input$TCA_DATASET_CHOICE]]@cci_detected[
      ID %in% input$TCA_TISSUE_CHOICE]# & REGULATION %in% input$TCA_REGULATION_CHOICE]
    
    min_young <- min(dt[CCI_SCORE_YOUNG > 0][["CCI_SCORE_YOUNG"]])
    min_old <- min(dt[CCI_SCORE_OLD > 0][["CCI_SCORE_OLD"]])
    dt[, CCI_SCORE_YOUNG := ifelse(CCI_SCORE_YOUNG == 0, min_young, CCI_SCORE_YOUNG)]
    dt[, CCI_SCORE_OLD := ifelse(CCI_SCORE_OLD == 0, min_old, CCI_SCORE_OLD)]
    xlims <- c(min(dt[["CCI_SCORE_YOUNG"]]), max(dt[["CCI_SCORE_YOUNG"]]))
    ylims <- c(min(dt[["CCI_SCORE_OLD"]]), max(dt[["CCI_SCORE_OLD"]]))
    req(dt)
    if (!('All LRI' %in% input$TCA_LRI_CHOICE)) {
      dt <- dt[`LR_GENES` %in% input$TCA_LRI_CHOICE]
    }
    dt <- dt[, c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC_ALL", "REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD")]
    setnames(
      dt,
      old = c("REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD", "LOG2FC_ALL"),
      new = c("Age Regulation","Score Young", "Score Old", "LOG2FC")
    )
    show_scores(data = dt, xlims = xlims, ylims = ylims)
  })
}

get_TCA_SCORES_text <- function(
  input
) {
  renderPrint({
    req(input$TCA_DATASET_CHOICE, input$TCA_TISSUE_CHOICE)#, input$TCA_REGULATION_CHOICE)
    dt <-  DATASETS_COMBINED[[input$TCA_DATASET_CHOICE]]@cci_detected[
      ID %in% input$TCA_TISSUE_CHOICE]# & REGULATION %in% input$TCA_REGULATION_CHOICE]
    min_young <- min(dt[CCI_SCORE_YOUNG > 0][["CCI_SCORE_YOUNG"]])
    min_old <- min(dt[CCI_SCORE_OLD > 0][["CCI_SCORE_OLD"]])
    dt[, CCI_SCORE_YOUNG := ifelse(CCI_SCORE_YOUNG == 0, min_young, CCI_SCORE_YOUNG)]
    dt[, CCI_SCORE_OLD := ifelse(CCI_SCORE_OLD == 0, min_old, CCI_SCORE_OLD)]
    req(dt)
    if (!('All LRI' %in% input$TCA_LRI_CHOICE)) {
      dt <- dt[`LR_GENES` %in% input$TCA_LRI_CHOICE]
    }
    dt <- dt[, c("ID", "LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC_ALL", "REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD")]
    setnames(
      dt, 
      old = c("ID", "LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD", "LOG2FC_ALL"),
      new = c("Tissue", "LRI", "Emitter Cell Type", "Receiver Cell Type", "Age Regulation", "Score Young", "Score Old", "LOG2FC")
    )
    brushedPoints(dt, input$TCA_SCORES_brush)
  })
}

plot_TCA_LRIFC <- function(
  input
) {
  renderPlot({
    req(input$TCA_DATASET_CHOICE, input$TCA_TISSUE_CHOICE)#, input$TCA_REGULATION_CHOICE)
    dt <-  DATASETS_COMBINED[[input$TCA_DATASET_CHOICE]]@cci_detected[
      ID %in% input$TCA_TISSUE_CHOICE]# & REGULATION %in% input$TCA_REGULATION_CHOICE]
    dt[, LOG2FC_L := log2(
      pmin(L1_EXPRESSION_OLD, L2_EXPRESSION_OLD, na.rm = TRUE)
      /
        pmin(L1_EXPRESSION_YOUNG, L2_EXPRESSION_YOUNG, na.rm = TRUE)
    ) ]
    dt[, LOG2FC_R := log2(
      pmin(R1_EXPRESSION_OLD, R2_EXPRESSION_OLD, R3_EXPRESSION_OLD, na.rm = TRUE)
      /
        pmin(R1_EXPRESSION_YOUNG, R2_EXPRESSION_YOUNG, R3_EXPRESSION_YOUNG, na.rm = TRUE)
    ) ]
    max_L <- max(dt[is.finite(LOG2FC_L)][["LOG2FC_L"]])
    min_L <- min(dt[is.finite(LOG2FC_L)][["LOG2FC_L"]])
    max_L <- max(max_L, -min_L)
    min_L <- min(-max_L, min_L)
    dt[, LOG2FC_L := ifelse(is.infinite(LOG2FC_L) & LOG2FC_L > 0, max_L,
                            ifelse(is.infinite(LOG2FC_L) & LOG2FC_L < 0, min_L, LOG2FC_L))]
    max_R <- ceiling(max(dt[is.finite(LOG2FC_R)][["LOG2FC_R"]]))
    min_R <- floor(min(dt[is.finite(LOG2FC_R)][["LOG2FC_R"]]))
    max_R <- max(max_R, -min_R)
    min_R <- min(-max_R, min_R)
    dt[, LOG2FC_R := ifelse(is.infinite(LOG2FC_R) & LOG2FC_R > 0, max_R,
                            ifelse(is.infinite(LOG2FC_R) & LOG2FC_R < 0, min_R, LOG2FC_R))]
    xlims <- c(min_L, max_L)
    ylims <- c(min_R, max_R)
    req(dt)
    if (!('All LRI' %in% input$TCA_LRI_CHOICE)) {
      dt <- dt[`LR_GENES` %in% input$TCA_LRI_CHOICE]
    }
    dt <- dt[, c("ID", "LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC_ALL", "REGULATION", "LOG2FC_L", "LOG2FC_R")]
    setnames(
      dt,
      old = c("ID", "REGULATION", "LOG2FC_ALL"),
      new = c("Tissue", "Age Regulation", "LOG2FC")
    )
    show_LRIFC(data = dt, xlims = xlims, ylims = ylims)
  })
}

get_TCA_LRIFC_text <- function(
  input
) {
  renderPrint({
    req(input$TCA_DATASET_CHOICE, input$TCA_TISSUE_CHOICE)#, input$TCA_REGULATION_CHOICE)
    dt <-  DATASETS_COMBINED[[input$TCA_DATASET_CHOICE]]@cci_detected[
      ID %in% input$TCA_TISSUE_CHOICE]# & REGULATION %in% input$TCA_REGULATION_CHOICE]
    dt[, LOG2FC_L := log2(
      pmin(L1_EXPRESSION_OLD, L2_EXPRESSION_OLD, na.rm = TRUE)
      /
        pmin(L1_EXPRESSION_YOUNG, L2_EXPRESSION_YOUNG, na.rm = TRUE)
    ) ]
    dt[, LOG2FC_R := log2(
      pmin(R1_EXPRESSION_OLD, R2_EXPRESSION_OLD, R3_EXPRESSION_OLD, na.rm = TRUE)
      /
        pmin(R1_EXPRESSION_YOUNG, R2_EXPRESSION_YOUNG, R3_EXPRESSION_YOUNG, na.rm = TRUE)
    ) ]
    max_L <- max(dt[is.finite(LOG2FC_L)][["LOG2FC_L"]])
    min_L <- min(dt[is.finite(LOG2FC_L)][["LOG2FC_L"]])
    max_L <- max(max_L, -min_L)
    min_L <- min(-max_L, min_L)
    dt[, LOG2FC_L := ifelse(is.infinite(LOG2FC_L) & LOG2FC_L > 0, max_L,
                            ifelse(is.infinite(LOG2FC_L) & LOG2FC_L < 0, min_L, LOG2FC_L))]
    max_R <- ceiling(max(dt[is.finite(LOG2FC_R)][["LOG2FC_R"]]))
    min_R <- floor(min(dt[is.finite(LOG2FC_R)][["LOG2FC_R"]]))
    max_R <- max(max_R, -min_R)
    min_R <- min(-max_R, min_R)
    dt[, LOG2FC_R := ifelse(is.infinite(LOG2FC_R) & LOG2FC_R > 0, max_R,
                            ifelse(is.infinite(LOG2FC_R) & LOG2FC_R < 0, min_R, LOG2FC_R))]
    xlims <- c(min_L, max_L)
    ylims <- c(min_R, max_R)
    if (!('All LRI' %in% input$TCA_LRI_CHOICE)) {
      dt <- dt[`LR_GENES` %in% input$TCA_LRI_CHOICE]
    }
    dt <- dt[, c("ID", "LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC_ALL", "REGULATION", "LOG2FC_L", "LOG2FC_R")]
    setnames(
      dt, 
      old = c("ID", "LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "REGULATION", "LOG2FC_ALL"),
      new = c("Tissue", "LRI", "Emitter Cell Type", "Receiver Cell Type", "Age Regulation", "LOG2FC")
    )
    brushedPoints(dt, input$TCA_LRIFC_brush)
  })
}

## TCA GLOBAL sidebar ####

## TCA GlOBAL mainpanel ####

get_TCA_global_details <- function(
  input
) {
  renderUI({
    if (input$TCA_GLOBAL_DETAILS_CHOICE == "Table") {
      DT::dataTableOutput("TCA_GLOBAL_TABLE")
    } else if (input$TCA_GLOBAL_DETAILS_CHOICE == "ORA Plot") {
      plotOutput("TCA_GLOBAL_ORA_PLOT", height = "600px")
    }
  })
}

get_TCA_global_table <- function(
  input
) {
  DT::renderDT({
    req(input$TCA_DATASET_CHOICE, input$TCA_GLOBAL_TABLE_CHOICE, input$TCA_GLOBAL_ORA_TYPE_CHOICE)
    categ <- input$TCA_GLOBAL_TABLE_CHOICE
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
    if (input$TCA_GLOBAL_ORA_TYPE_CHOICE == "Up") {
      setorder(dt, -`ORA Score Up`, na.last = TRUE)
    } else if (input$TCA_GLOBAL_ORA_TYPE_CHOICE == "Down") {
      setorder(dt, -`ORA Score Down`, na.last = TRUE)
    } 
    show_DT(
      dt,
      new_cols,
      new_cols[-1]
    )
  })
}

plot_TCA_global_ora <- function(
  input
) {
  renderPlot({
    req(input$TCA_DATASET_CHOICE, input$TCA_GLOBAL_TABLE_CHOICE, input$TCA_GLOBAL_ORA_TYPE_CHOICE)
    categ <- input$TCA_GLOBAL_TABLE_CHOICE
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
    if (input$TCA_GLOBAL_ORA_TYPE_CHOICE == "Up") {
      reg <- "UP"
    } else if (input$TCA_GLOBAL_ORA_TYPE_CHOICE == "Down") {
      reg <- "DOWN"
    } else if ( input$TCA_GLOBAL_ORA_TYPE_CHOICE == "Stable") {
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

## deprecated ####
# 
# get_TCA_cci_intro <- function(
#   input
# ) {
#   renderUI(
#     {
#       req(input$TCA_DATASET_CHOICE)
#       color_theme <- "color: rgb(20 120 206)"
#       tags$div(
#         tags$p(
#           "Table with all interactions from the selected tissues, or Volcan Plot."
#         ),
#         style = "font-size: 1.5em; text-align: justify;"
#       )
#     }
#   )
# }
# 
# choose_TCA_regulation <- function(
#   input
# ) {
#   renderUI({
#     choices <- c("UP", "DOWN", "FLAT")
#     pickerInput(
#       inputId = "TCA_REGULATION_CHOICE",
#       label = "Regulation",
#       choices = choices,
#       selected = choices,
#       options = list(`actions-box` = TRUE),
#       multiple = TRUE
#     )
#   })
# }
# 
# get_TCA_global_intro <- function(
#   input
# ) {
#   renderUI(
#     {
#       req(input$TCA_DATASET_CHOICE)
#       color_theme <- "color: rgb(20 120 206)"
#       tags$div(
#         tags$p(
#           "Analysis of signals that are shared between several tissues."
#         ),
#         tags$p(
#           "You can also see a plot of the most over-represented signals. Note that a term with a large ORA Score ",
#           "does not necessarily take part in a lot of tissues (this information is given in the table)."
#         ),
#         tags$p(
#           "We will clarify the names of the column in the tables. For now, here is what they mean:",
#           "ID, global ORA Score Down, global ORA Score Up, Number of tissues in which the id is ",
#           "down over-represented, same but for up over-represented, number of cell-cell interactions in which the ID ",
#           "takes part, same but only in young cells, same but only in old cells."
#         ),
#         style = "font-size: 1.5em; text-align: justify;"
#       )
#     }
#   )
# }


