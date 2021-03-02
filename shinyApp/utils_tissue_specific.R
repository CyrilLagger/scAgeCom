choose_TSA_dataset <- function(
  input
) {
  renderUI({
    pickerInput(
      inputId = "TSA_DATASET_CHOICE",
      #label = "Dataset",
      choices = names(DATASETS_COMBINED),
      options = list(
        `actions-box` = TRUE
      ),
      multiple = FALSE,
      inline = FALSE
    )
  })
}

choose_TSA_tissue <- function(
  input
) {
  renderUI({
    req(input$TSA_DATASET_CHOICE)
    choices <- sort(unique(DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected$ID))
    pickerInput(
      inputId = "TSA_TISSUE_CHOICE",
      #label = "Tissue",
      choices = choices,
      options = list(`actions-box` = TRUE),
      multiple = FALSE,
      inline = FALSE
    )
  })
}

get_TSA_title <- function(
  input
) {
  renderUI(
    {
      #req(input$TSA_TISSUE_CHOICE, input$TSA_DATASET_CHOICE)
      #req(input$TSA_DATASET_CHOICE)
      #color_theme <- "color: rgb(20 120 206)"
      # tags$p(
      #   "Analysis of the ",
      #   span(input$TSA_TISSUE_CHOICE, style = color_theme),
      #   " from the ",
      #   span(input$TSA_DATASET_CHOICE, style = color_theme)
      #   )
      # tags$p(
      #   div(style="display: inline-block;", "Please choose a dataset ("),
      #   div(style="display: inline-block;margin-top: 25px;", uiOutput("TSA_DATASET_CHOICE", inline = TRUE)),
      #   div(style="display: inline-block;", ") and a tissue (" ),
      #   div(style="display: inline-block;margin-top: 25px;", uiOutput("TSA_TISSUE_CHOICE", inline = TRUE)),
      #   div(style="display: inline-block;", ")" )
      # )
      # 
      tags$p(
        div(style="display: inline-block;", "Please choose a dataset and a tissue: "),
        div(style="display: inline-block;margin-top: 25px;", uiOutput("TSA_DATASET_CHOICE", inline = TRUE)),
        div(style="display: inline-block;margin-top: 25px;", uiOutput("TSA_TISSUE_CHOICE", inline = TRUE)),
      )
      
    }
  )
}

get_TSA_overview_intro <- function(
  input
) {
  renderUI(
    {
      req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE)
      color_theme <- "color: rgb(20 120 206)"
      tags$div(
        tags$p(
          "We give the number of detected cell-cell interactions (CCI) in the ",
          span(input$TSA_TISSUE_CHOICE, style = color_theme), 
          ", and how many of them change significantly with age.",
          "Note that some CCI can appear or disappear with age whereas others can be down- or up-regulated",
          " but still detected in both young and old samples."
        ),
        tags$p("We then present several networks that provide an overview of the intercellular communication in this tissue."),
        style = "font-size: 1.5em; text-align: justify;"
      )
    }
  )
}

get_TSA_overview_table <- function(
  input
) {
  DT::renderDT({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE)
    dt <-  CCI_SUMMARY[[input$TSA_DATASET_CHOICE]][["ID"]][
      ID == input$TSA_TISSUE_CHOICE
    ]
    req(dt)
    old_col_names <- c(
      "N_CELLTYPES", "N_CCI", "N_CCI_FLAT", "N_CCI_DOWN", "N_CCI_UP", 
      "N_CCI_NON_SIGNIFICANT_CHANGE"
    )
    new_col_names <- c(
      "Total cell-types", "Total CCI", "Flat CCI",
      "Down CCI", "Up CCI", "Non-significant CCI"
    )
    dt <- dt[, old_col_names, with = FALSE]
    setnames(dt, old_col_names, new_col_names)
    show_DT(
      data = dt,
      cols_to_show = new_col_names,
      cols_numeric = NULL,
      table_title = NULL, # "Cell-Cell Interaction Counts",
      options = list(dom = "t"),#list(pageLength = 10, lengthChange = FALSE),
      callback = JS(
        "var tips = [
      'Number of cell-types in the tissue of interest',
      'Number of cell-cell interactions detected in the tissue of interest',
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

get_TSA_network_intro <- function(
  input
) {
  renderUI(
    {
      req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE)
      color_theme <- "color: rgb(20 120 206)"
      tags$div(
        tags$p(
          "Please choose between different network types (more info to come about each network).",
          "You can interact with the networks. When the network type is 'ORA' you can find more information ",
          "about each egde by hovering it."
        ),
        style = "font-size: 1.5em; text-align: justify;"
      )
    }
  )
}

plot_TSA_network <- function(
  input
) {
  renderVisNetwork({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, 
        input$TSA_NETWORK_TYPE_CHOICE, input$TSA_NETWORK_LAYOUT_CHOICE)
    obj <- DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]
    net_type <- input$TSA_NETWORK_TYPE_CHOICE
    net_lay <- input$TSA_NETWORK_LAYOUT_CHOICE
    if (net_type == "Young") {
      net_type2 <- "COUNTS_COND1"
    } else if (net_type == "Old") {
      net_type2 <- "COUNTS_COND2"
    } else if (net_type == "Change") {
      net_type2 <- "COUNTS_DIFF"
    } else if (net_type == "ORA") {
      net_type2 <- "ORA"
    }
    if (net_lay == "Circular") {
      net_lay2 <- "celltypes"
    } else if (net_lay == "Bipartite") {
      net_lay2 <- "bipartite"
    }
    req(obj)
    BuildNetwork(
      object = obj,
      network_type = net_type2,
      network_layout = net_lay2,
      ID = input$TSA_TISSUE_CHOICE
    )
    # ora_tables <- scDiffCom::get_ora_tables(obj)
    # names(ora_tables) <- c("GO Terms", "LR_GENES", "LR_CELLTYPE", "Cell Families")
    # ora_tables <- lapply(
    #   ora_tables,
    #   function(ORA_dt) {
    #     dt <- copy(ORA_dt)
    #     setnames(
    #       dt,
    #       new = c("OR_UP", "pval_adjusted_UP", "OR_DOWN", "pval_adjusted_DOWN",
    #               "OR_FLAT", "pval_adjusted_FLAT"),
    #       old = c("OR_UP", "BH_P_VALUE_UP", "OR_DOWN", "BH_P_VALUE_DOWN",
    #               "OR_FLAT", "BH_P_VALUE_FLAT")
    #     )
    #     return(dt)
    #   }
    # )
    # obj <- scDiffCom::set_ora_tables(obj, ora_tables)
    # cci_table_filtered <- copy(scDiffCom:::get_cci_table_filtered(obj))
    # setnames(
    #   cci_table_filtered,
    #   new = c("L_CELLTYPE", "R_CELLTYPE", "BH_PVAL_DIFF", "LR_NAME"),
    #   old = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "Adj. P-Value", "LR_GENES")
    # )
    # obj <- scDiffCom:::set_cci_table_filtered(obj, cci_table_filtered)
    # req(obj)
    # scDiffCom::build_celltype_bipartite_graph(
    #   object = obj,
    #   disperse = FALSE,
    #   dir = NULL
    # )
  })
}

get_TSA_cci_intro <- function(
  input
) {
  renderUI(
    {
      req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE)
      color_theme <- "color: rgb(20 120 206)"
      tags$div(
        tags$p(
          "We provide the full table of all detected CCI in this tissue, including relevant scores. ",
          "In addition, a Volcano Plot highlights how they are distributed with respect to their significance and ",
          "strength of change with age, whereas a Score Plot indicates their prevalence in both young and old samples. ",
          "The latter graph can also be useful to observe strong signals that do not change with age."
        ),
        style = "font-size: 1.5em; text-align: justify;"
      )
    }
  )
}

choose_TSA_emitter <- function(
  input
) {
  renderUI({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE)
    choices <- sort(unique(DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected[
      ID == input$TSA_TISSUE_CHOICE
    ][["EMITTER_CELLTYPE"]]))
    pickerInput(
      inputId = "TSA_EMITTER_CHOICE",
      label = "Filter by Emitter Cell Types",
      choices = choices,
      selected = choices,
      options = list(`actions-box` = TRUE),
      multiple = TRUE
    )
  })
}

choose_TSA_receiver <- function(
  input
) {
  renderUI({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE)
    choices <- sort(unique(DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected[
      ID == input$TSA_TISSUE_CHOICE
    ][["RECEIVER_CELLTYPE"]]))
    pickerInput(
      inputId = "TSA_RECEIVER_CHOICE",
      label = "Filter by Receiver Cell Types",
      choices = choices,
      selected = choices,
      options = list(`actions-box` = TRUE),
      multiple = TRUE
    )
  })
}

choose_TSA_lri <- function(
  input 
) {
  renderUI({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE)
    ALL_LRI_LABEL = 'All LRI'
    choices <-
      c(ALL_LRI_LABEL,
        sort(unique(DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected[ID == input$TSA_TISSUE_CHOICE][["LR_GENES"]]))
      )
    selectizeInput(inputId = 'TSA_LRI_CHOICE',
                   label = 'Filter by Ligand-Receptor Interactions',
                   choices = choices,
                   selected = ALL_LRI_LABEL,
                   multiple = TRUE,
                   options = list(allowEmptyOption = TRUE,
                                  placeholder = 'Type LRIs')
    )
  })
}

get_TSA_slider_log2fc <- function(
  input
) {
  renderUI({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE)
    max_val <- ceiling(max(DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected[
      ID == input$TSA_TISSUE_CHOICE
    ][["LOG2FC"]]))
    sliderInput(
      inputId = "TSA_SLIDER_LOG2FC",
      label = "Filter by LOG2FC",
      min = 0,
      max = max_val,
      value = 0,
      step = 0.01
    )
  })
}

get_TSA_cci_details <- function(
  input
) {
  renderUI({
    if (input$TSA_CCI_DETAILS_CHOICE == "CCI Table") {
      DT::dataTableOutput("TSA_INTERACTION_TABLE")
    } else if (input$TSA_CCI_DETAILS_CHOICE == "Volcano Plot") {
      plotOutput("TSA_VOLCANO_PLOT", brush = "TSA_VOLCANO_brush", height = "600px")
    } else if (input$TSA_CCI_DETAILS_CHOICE == "Score Plot") {
      plotOutput("TSA_SCORES_PLOT", brush = "TSA_SCORES_brush", height = "600px")
    } else if (input$TSA_CCI_DETAILS_CHOICE == "LRI-FC Plot") {
      plotOutput("TSA_LRIFC_PLOT", brush = "TSA_LRIFC_brush", height = "600px")
    }
  })
}

get_TSA_cci_text <- function(
  input
) {
  renderUI({
    if (input$TSA_CCI_DETAILS_CHOICE == "CCI Table") {
      NULL
    } else if (input$TSA_CCI_DETAILS_CHOICE == "Volcano Plot") {
      verbatimTextOutput("TSA_VOLCANO_TEXTOUTPUT")
    } else if (input$TSA_CCI_DETAILS_CHOICE == "Score Plot") {
      verbatimTextOutput("TSA_SCORES_TEXTOUTPUT")
    } else if (input$TSA_CCI_DETAILS_CHOICE == "LRI-FC Plot") {
      verbatimTextOutput("TSA_LRIFC_TEXTOUTPUT")
    }
  })
}

get_TSA_interaction_table <- function(
  input
) {
  DT::renderDT({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_EMITTER_CHOICE, input$TSA_RECEIVER_CHOICE,
        input$TSA_LRI_CHOICE, input$TSA_SLIDER_PVALUE, input$TSA_SLIDER_LOG2FC)
    dt <-  DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected[
      ID == input$TSA_TISSUE_CHOICE
    ]
    req(dt)
    if ('All LRI' %in% input$TSA_LRI_CHOICE) {
      dt <- dt[
        `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
          `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
          `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
          abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    } else {
      dt <- dt[
        `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
          `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
          `LR_GENES` %in% input$TSA_LRI_CHOICE &
          `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
          abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    }
    setorder(
      dt,
      -LOG2FC,
      `BH_P_VALUE_DE`
    )
    setnames(
      dt,
      old = c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "LOG2FC", "BH_P_VALUE_DE",
              "REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD"),
      new = c("LRI", "Emitter Cell Type", "Receiver Cell Type", "LOG2FC", "Adj. p-value",
              "Age-Regulation", "Score Young", "Score Old")
    )
    show_DT(
      dt,
      c("LRI", "Emitter Cell Type", "Receiver Cell Type", "LOG2FC", "Adj. p-value",
        "Age-Regulation", "Score Young", "Score Old"),
      c("LOG2FC", "Adj. p-value", "Score Young", "Score Old"),
      "Cell-Cell Interation Table"
    )
  })
}

plot_TSA_VOLCANO <- function(
  input
) {
  renderPlot({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_EMITTER_CHOICE, input$TSA_RECEIVER_CHOICE,
        input$TSA_LRI_CHOICE, input$TSA_SLIDER_PVALUE, input$TSA_SLIDER_LOG2FC)
    dt <-  DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected[
      ID == input$TSA_TISSUE_CHOICE
    ]
    req(dt)
    xlims <- c(min(dt[["LOG2FC"]]), max(dt[["LOG2FC"]]))
    ylims <- c(0, 4)
    if ('All LRI' %in% input$TSA_LRI_CHOICE) {
      dt <- dt[
        `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
          `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
          `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
          abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    } else {
      dt <- dt[
        `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
          `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
          `LR_GENES` %in% input$TSA_LRI_CHOICE &
          `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
          abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    }
    dt[, minus_log10_pval := -log10(`BH_P_VALUE_DE` + 1E-4)]
    dt <- dt[, c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC", "minus_log10_pval", "REGULATION")]
    setnames(
      dt,
      old = "REGULATION",
      new = "Age Regulation"
    )
    show_volcano(
      data = dt,
      xlims = xlims,
      ylims = ylims
    )
  })
}

get_TSA_VOLCANO_text <- function(
  input
) {
  renderPrint({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_EMITTER_CHOICE, input$TSA_RECEIVER_CHOICE,
        input$TSA_LRI_CHOICE, input$TSA_SLIDER_PVALUE, input$TSA_SLIDER_LOG2FC)
    dt <-  DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected[
      ID == input$TSA_TISSUE_CHOICE
    ]
    req(dt)
    if ('All LRI' %in% input$TSA_LRI_CHOICE) {
      dt <- dt[
        `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
          `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
          `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
          abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    } else {
      dt <- dt[
        `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
          `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
          `LR_GENES` %in% input$TSA_LRI_CHOICE &
          `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
          abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    }
    dt[, minus_log10_pval := -log10(`BH_P_VALUE_DE` + 1E-4)]
    dt <- dt[, c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "REGULATION", "LOG2FC", "minus_log10_pval")]
    setnames(
      dt, 
      old = c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "REGULATION"),
      new = c("LRI", "Emitter Cell Type", "Receiver Cell Type", "Age Regulation")
    )
    brushedPoints(dt, input$TSA_VOLCANO_brush, xvar = "LOG2FC", yvar = "minus_log10_pval")
  })
}

plot_TSA_SCORES <- function(
  input
) {
  renderPlot({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_EMITTER_CHOICE, input$TSA_RECEIVER_CHOICE,
        input$TSA_LRI_CHOICE, input$TSA_SLIDER_PVALUE, input$TSA_SLIDER_LOG2FC)
    dt <-  DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected[
      ID == input$TSA_TISSUE_CHOICE
    ]
    min_young <- min(dt[CCI_SCORE_YOUNG > 0][["CCI_SCORE_YOUNG"]])
    min_old <- min(dt[CCI_SCORE_OLD > 0][["CCI_SCORE_OLD"]])
    dt[, CCI_SCORE_YOUNG := ifelse(CCI_SCORE_YOUNG == 0, min_young, CCI_SCORE_YOUNG)]
    dt[, CCI_SCORE_OLD := ifelse(CCI_SCORE_OLD == 0, min_old, CCI_SCORE_OLD)]
    xlims <- c(min(dt[["CCI_SCORE_YOUNG"]]), max(dt[["CCI_SCORE_YOUNG"]]))
    ylims <- c(min(dt[["CCI_SCORE_OLD"]]), max(dt[["CCI_SCORE_OLD"]]))
    req(dt)
    if ('All LRI' %in% input$TSA_LRI_CHOICE) {
      dt <- dt[
        `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
          `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
          `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
          abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    } else {
      dt <- dt[
        `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
          `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
          `LR_GENES` %in% input$TSA_LRI_CHOICE &
          `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
          abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    }
    dt <- dt[, c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC", "REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD")]
    setnames(
      dt,
      old = c("REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD"),
      new = c("Age Regulation","Score Young", "Score Old")
    )
    show_scores(data = dt, xlims = xlims, ylims = ylims)
  })
}

get_TSA_SCORES_text <- function(
  input
) {
  renderPrint({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_EMITTER_CHOICE, input$TSA_RECEIVER_CHOICE,
        input$TSA_LRI_CHOICE, input$TSA_SLIDER_PVALUE, input$TSA_SLIDER_LOG2FC)
    dt <-  DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected[
      ID == input$TSA_TISSUE_CHOICE
    ]
    min_young <- min(dt[CCI_SCORE_YOUNG > 0][["CCI_SCORE_YOUNG"]])
    min_old <- min(dt[CCI_SCORE_OLD > 0][["CCI_SCORE_OLD"]])
    dt[, CCI_SCORE_YOUNG := ifelse(CCI_SCORE_YOUNG == 0, min_young, CCI_SCORE_YOUNG)]
    dt[, CCI_SCORE_OLD := ifelse(CCI_SCORE_OLD == 0, min_old, CCI_SCORE_OLD)]
    req(dt)
    if ('All LRI' %in% input$TSA_LRI_CHOICE) {
      dt <- dt[
        `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
          `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
          `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
          abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    } else {
      dt <- dt[
        `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
          `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
          `LR_GENES` %in% input$TSA_LRI_CHOICE &
          `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
          abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    }
    dt <- dt[, c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC", "REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD")]
    setnames(
      dt, 
      old = c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD"),
      new = c("LRI", "Emitter Cell Type", "Receiver Cell Type", "Age Regulation", "Score Young", "Score Old")
    )
    brushedPoints(dt, input$TSA_SCORES_brush)
  })
}

plot_TSA_LRIFC <- function(
  input
) {
  renderPlot({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_EMITTER_CHOICE, input$TSA_RECEIVER_CHOICE,
        input$TSA_LRI_CHOICE, input$TSA_SLIDER_PVALUE, input$TSA_SLIDER_LOG2FC)
    dt <-  DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected[
      ID == input$TSA_TISSUE_CHOICE
    ]
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
    if ('All LRI' %in% input$TSA_LRI_CHOICE) {
      dt <- dt[
        `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
          `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
          `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
          abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    } else {
      dt <- dt[
        `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
          `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
          `LR_GENES` %in% input$TSA_LRI_CHOICE &
          `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
          abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    }
    dt <- dt[, c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC", "REGULATION", "LOG2FC_L", "LOG2FC_R")]
    setnames(
      dt,
      old = c("REGULATION"),
      new = c("Age Regulation")
    )
    show_LRIFC(data = dt, xlims = xlims, ylims = ylims)
  })
}

get_TSA_LRIFC_text <- function(
  input
) {
  renderPrint({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_EMITTER_CHOICE, input$TSA_RECEIVER_CHOICE,
        input$TSA_LRI_CHOICE, input$TSA_SLIDER_PVALUE, input$TSA_SLIDER_LOG2FC)
    dt <-  DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected[
      ID == input$TSA_TISSUE_CHOICE
    ]
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
    if ('All LRI' %in% input$TSA_LRI_CHOICE) {
      dt <- dt[
        `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
          `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
          `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
          abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    } else {
      dt <- dt[
        `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
          `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
          `LR_GENES` %in% input$TSA_LRI_CHOICE &
          `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
          abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    }
    dt <- dt[, c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC", "REGULATION", "LOG2FC_L", "LOG2FC_R")]
    setnames(
      dt, 
      old = c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "REGULATION"),
      new = c("LRI", "Emitter Cell Type", "Receiver Cell Type", "Age Regulation")
    )
    brushedPoints(dt, input$TSA_LRIFC_brush)
  })
}

get_TSA_ora_intro <- function(
  input
) {
  renderUI(
    {
      req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE)
      color_theme <- "color: rgb(20 120 206)"
      tags$div(
        tags$p(
          "We present the result of our over-representation analysis. Terms from various categories can be associated to each CCI, such as  ",
          "the genes and the cell-types involved in the interactions but also GO terms or cell-type families. For each ",
          "category, we ask (via a Fisher's exact test) if specific terms are more associated with the up, down or stable CCI ",
          "rather than the other CCI. The ORA Score is given as log(Odds Ratio) times -log(Adj. P-Value). ",
          "For instance, if the interacting pair of genes A:B has a large (say) Up ORA Score, it means that it takes more part in ",
          "up-regulated CCI than what it is expected by chance and is thus a potentially interesting aging-related signal."
        ),
        tags$p(
          "Note that the absence of over-represented terms in a given category does not imply the absence of interesting biological",
          " process related this category."
        ),
        style = "font-size: 1.5em; text-align: justify;"
      )
    }
  )
}

choose_TSA_ORA_category <- function(
  input
) {
  renderUI({
    req(input$TSA_DATASET_CHOICE)
    #choices <- names(DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@ora_default)
    choices <- c("KEGG Pathways", "GO Terms", "LRIs", "Cell Types", "Cell Families")
    pickerInput(
      inputId = "TSA_ORA_CATEGORY_CHOICE",
      label = "Category",
      choices = choices,
      options = list(`actions-box` = TRUE),
      multiple = FALSE
    )
  })
}

get_TSA_ORA_slider_or <- function(
  input
) {
  renderUI({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_ORA_CATEGORY_CHOICE, input$TSA_ORA_TYPE_CHOICE)
    ora_table <- DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@ora_default[[input$TSA_ORA_CATEGORY_CHOICE]][
      ID == input$TSA_TISSUE_CHOICE
    ]
    req(ora_table)
    if(input$TSA_ORA_TYPE_CHOICE == "Up") {
      max_val <- ceiling(max(ora_table[["OR_UP"]]))
    } else if(input$TSA_ORA_TYPE_CHOICE == "Down") {
      max_val <- ceiling(max(ora_table[["OR_DOWN"]]))
    } else if(input$TSA_ORA_TYPE_CHOICE == "Stable") {
      max_val <- ceiling(max(ora_table[["OR_FLAT"]]))
    }
    sliderInput(
      inputId = "TSA_ORA_SLIDER_OR",
      label = "Filter by Odds Ratio",
      min = 1,
      max = max_val,
      value = 1,
      step = 0.01
    )
  })
}

get_TSA_ora_details <- function(
  input
) {
  renderUI({
    if (input$TSA_ORA_DETAILS_CHOICE == "ORA Table") {
      DT::dataTableOutput("TSA_ORA_TABLE")
    } else if (input$TSA_ORA_DETAILS_CHOICE == "ORA Plot") {
      plotOutput("TSA_ORA_PLOT", height = "800px")
    }
  })
}

get_TSA_ORA_table <- function(
  input
) {
  DT::renderDataTable({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_ORA_CATEGORY_CHOICE, input$TSA_ORA_TYPE_CHOICE)
    replacement <- data.table(
      VALUE_OLD = c("KEGG_PWS", "GO_TERMS", "LR_GENES", "ER_CELLTYPES", "ER_CELL_FAMILY"),
      VALUE_NEW = c("KEGG Pathways", "GO Terms", "LRIs", "Cell Types", "Cell Families")
    )
    replacement <- replacement[VALUE_NEW == input$TSA_ORA_CATEGORY_CHOICE ][["VALUE_OLD"]]
    dt <- DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@ora_default[[replacement]][
      ID == input$TSA_TISSUE_CHOICE
    ]
    req(dt)
    if(input$TSA_ORA_TYPE_CHOICE == "Up") {
      dt <- dt[`OR_UP` >= 1 & BH_P_VALUE_UP <= 0.05, c("VALUE", "ORA_SCORE_UP", "OR_UP", "BH_P_VALUE_UP")]
      # dt <- dt[`BH_P_VALUE_UP` <= input$TSA_ORA_SLIDER_PVALUE &
      #            `OR_UP` >= input$TSA_ORA_SLIDER_OR]
    } else if(input$TSA_ORA_TYPE_CHOICE == "Down") {
      dt <- dt[`OR_DOWN` >= 1 & BH_P_VALUE_DOWN <= 0.05, c("VALUE","ORA_SCORE_DOWN", "OR_DOWN", "BH_P_VALUE_DOWN")]
      # dt <- dt[`BH_P_VALUE_DOWN` <= input$TSA_ORA_SLIDER_PVALUE &
      #            `OR_DOWN` >= input$TSA_ORA_SLIDER_OR]
    } else if(input$TSA_ORA_TYPE_CHOICE == "Flat") {
      dt <- dt[`OR_FLAT` >= 1 & BH_P_VALUE_FLAT <= 0.05, c("VALUE","ORA_SCORE_FLAT", "OR_FLAT", "BH_P_VALUE_FLAT")]
      # dt <- dt[`BH_P_VALUE_FLAT` <= input$TSA_ORA_SLIDER_PVALUE &
      #            `OR_FLAT` >= input$TSA_ORA_SLIDER_OR]
    }
    setnames(
      dt,
      old = colnames(dt),
      new = c(input$TSA_ORA_CATEGORY_CHOICE, "ORA Score", "Odds Ratio", "Adj. p-value")
    )
    setorder(dt, -`ORA Score`)
    cols_numeric <- c("ORA Score", "Odds Ratio", "Adj. p-value")
    show_DT(
      dt,
      cols_to_show = colnames(dt),
      cols_numeric = cols_numeric,
      table_title = paste0(
        input$TSA_ORA_CATEGORY_CHOICE,
        " over-represented among ",
        input$TSA_ORA_TYPE_CHOICE,
        "-regulated CCIs"
      )
    )
  })
}

plot_TSA_ORA <- function(
  input
) {
  renderPlot({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_ORA_CATEGORY_CHOICE, input$TSA_ORA_TYPE_CHOICE)
    replacement <- data.table(
      VALUE_OLD = c("KEGG_PWS", "GO_TERMS", "LR_GENES", "ER_CELLTYPES", "ER_CELL_FAMILY"),
      VALUE_NEW = c("KEGG Pathways", "GO Terms", "LRIs", "Cell Types", "Cell Families")
    )
    replacement <- replacement[VALUE_NEW == input$TSA_ORA_CATEGORY_CHOICE ][["VALUE_OLD"]]
    if (input$TSA_ORA_TYPE_CHOICE == "Up") {
      reg <- "UP"
    } else if (input$TSA_ORA_TYPE_CHOICE == "Down") {
      reg <- "DOWN"
    } else if ( input$TSA_ORA_TYPE_CHOICE == "Flat") {
      reg <- "FLAT"
    }
    obj <- DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]
    req(obj)
    p <- PlotORA(
      object = obj,
      subID = input$TSA_TISSUE_CHOICE,
      category = replacement,
      regulation = reg,
      max_terms_show = 20,
      OR_threshold = 1,
      p_value_threshold = 0.05,
      stringent = FALSE
    )
    if (is.character(p)) {
      return(p)
    }
    p <- p + ggtitle(
      paste0(
        "Top-20 ",
        input$TSA_ORA_CATEGORY_CHOICE,
        " over-represented among ",
        input$TSA_ORA_TYPE_CHOICE,
        "-regulated CCIs"
      )
    ) + 
      ylab("") +
      theme(text=element_text(size=20)) +
      theme(plot.title = element_text(hjust = 0.5, size = 16)) 
    return(p)
  })
}

