choose_TSA_tissue <- function(
  input
) {
  renderUI({
    choices <- sort(unique(DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected$ID))
    pickerInput(
      inputId = "TSA_TISSUE_CHOICE",
      label = "Tissue",
      choices = choices,
      options = list(`actions-box` = TRUE),
      multiple = FALSE
    )
  })
}

get_TSA_title <- function(
  input
) {
  renderUI(
    {
      req(input$TSA_TISSUE_CHOICE, input$TSA_DATASET_CHOICE)
      color_theme <- "color: rgb(20 120 206)"
      tags$p(
        "Analysis of the ",
        span(input$TSA_TISSUE_CHOICE, style = color_theme),
        " from the ",
        span(input$TSA_DATASET_CHOICE, style = color_theme)
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
      "Number of cell-types", "Total CCI", "Flat CCI",
      "Down CCI", "Up CCI", "Non-significant CCI"
    )
    dt <- dt[, old_col_names, with = FALSE]
    setnames(dt, old_col_names, new_col_names)
    show_DT(
      data = dt,
      cols_to_show = new_col_names,
      cols_numeric = NULL,
      table_title = "Cell-Cell Interaction Counts",
      options = list(dom = "t"),#list(pageLength = 10, lengthChange = FALSE),
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
      label = "EMITTER_CELLTYPE",
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
      label = "RECEIVER_CELLTYPE",
      choices = choices,
      selected = choices,
      options = list(`actions-box` = TRUE),
      multiple = TRUE
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
      label = "LOG2FC Threshold",
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
    }
  })
}

get_TSA_interaction_table <- function(
  input
) {
  DT::renderDT({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_EMITTER_CHOICE, input$TSA_RECEIVER_CHOICE,
        input$TSA_SLIDER_PVALUE, input$TSA_SLIDER_LOG2FC)
    dt <-  DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected[
      ID == input$TSA_TISSUE_CHOICE
      ]
    req(dt)
    dt <- dt[
      `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
        `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
        `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
        abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    setorder(
      dt,
      -LOG2FC,
      `BH_P_VALUE_DE`
    )
    show_DT(
      dt,
      c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "LOG2FC", "BH_P_VALUE_DE",
        "REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD"),
      c("LOG2FC", "BH_P_VALUE_DE", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD"),
      "Table of all detected CCIs"
    )
  })
}

plot_TSA_VOLCANO <- function(
  input
) {
  renderPlot({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_EMITTER_CHOICE, input$TSA_RECEIVER_CHOICE,
        input$TSA_SLIDER_PVALUE, input$TSA_SLIDER_LOG2FC)
    dt <-  DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected[
      ID == input$TSA_TISSUE_CHOICE
      ]
    req(dt)
    dt <- dt[
      `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
        `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
        `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
        abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    dt[, mlog10_pval := -log10(`BH_P_VALUE_DE` + 1E-4)]
    dt <- dt[, c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC", "mlog10_pval", "REGULATION")]
    show_volcano(dt)
  })
}

get_TSA_VOLCANO_text <- function(
  input
) {
  renderPrint({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_EMITTER_CHOICE, input$TSA_RECEIVER_CHOICE,
        input$TSA_SLIDER_PVALUE, input$TSA_SLIDER_LOG2FC)
    dt <-  DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected[
      ID == input$TSA_TISSUE_CHOICE
      ]
    req(dt)
    dt <- dt[
      `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
        `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
        `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
        abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    dt[, mlog10_pval := -log10(`BH_P_VALUE_DE` + 1E-4)]
    dt <- dt[, c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC", "mlog10_pval", "REGULATION")]
    brushedPoints(dt, input$TSA_VOLCANO_brush, xvar = "LOG2FC", yvar = "mlog10_pval")
  })
}

plot_TSA_SCORES <- function(
  input
) {
  renderPlot({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_EMITTER_CHOICE, input$TSA_RECEIVER_CHOICE,
        input$TSA_SLIDER_PVALUE, input$TSA_SLIDER_LOG2FC)
    dt <-  DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected[
      ID == input$TSA_TISSUE_CHOICE
      ]
    dt <- dt[
      `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
        `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
        `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
        abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    dt <- dt[, c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC", "REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD")]
    show_scores(dt)
  })
}

get_TSA_SCORES_text <- function(
  input
) {
  renderPrint({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_EMITTER_CHOICE, input$TSA_RECEIVER_CHOICE,
        input$TSA_SLIDER_PVALUE, input$TSA_SLIDER_LOG2FC)
    dt <-  DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@cci_detected[
      ID == input$TSA_TISSUE_CHOICE
      ]
    req(dt)
    dt <- dt[
      `EMITTER_CELLTYPE` %in% input$TSA_EMITTER_CHOICE &
        `RECEIVER_CELLTYPE` %in% input$TSA_RECEIVER_CHOICE &
        `BH_P_VALUE_DE` <= input$TSA_SLIDER_PVALUE &
        abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    dt <- dt[, c("LR_GENES", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",  "LOG2FC", "REGULATION", "CCI_SCORE_YOUNG", "CCI_SCORE_OLD")]
    brushedPoints(dt, input$TSA_SCORES_brush)
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
    choices <- c("KEGG_PWS", "GO_TERMS", "LR_GENES", "ER_CELLTYPES", "ER_CELL_FAMILY")
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
      label = "Odds Ratio Threshold",
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
    dt <- DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]@ora_default[[input$TSA_ORA_CATEGORY_CHOICE]][
      ID == input$TSA_TISSUE_CHOICE
      ]
    req(dt)
    if(input$TSA_ORA_TYPE_CHOICE == "Up") {
      dt <- dt[`OR_UP` >= 1, c("VALUE", "ORA_SCORE_UP", "OR_UP", "BH_P_VALUE_UP")]
      dt <- dt[`BH_P_VALUE_UP` <= input$TSA_ORA_SLIDER_PVALUE &
                 `OR_UP` >= input$TSA_ORA_SLIDER_OR]
      setorder(dt, `BH_P_VALUE_UP`)
      cols_numeric <- c("ORA_SCORE_UP", "OR_UP", "BH_P_VALUE_UP")
    } else if(input$TSA_ORA_TYPE_CHOICE == "Down") {
      dt <- dt[`OR_DOWN` >= 1, c("VALUE","ORA_SCORE_DOWN", "OR_DOWN", "BH_P_VALUE_DOWN")]
      dt <- dt[`BH_P_VALUE_DOWN` <= input$TSA_ORA_SLIDER_PVALUE &
                 `OR_DOWN` >= input$TSA_ORA_SLIDER_OR]
      setorder(dt, `BH_P_VALUE_DOWN`)
      cols_numeric <- c("ORA_SCORE_DOWN", "OR_DOWN", "BH_P_VALUE_DOWN")
    } else if(input$TSA_ORA_TYPE_CHOICE == "Stable") {
      dt <- dt[`OR_FLAT` >= 1, c("VALUE","ORA_SCORE_FLAT", "OR_FLAT", "BH_P_VALUE_FLAT")]
      dt <- dt[`BH_P_VALUE_FLAT` <= input$TSA_ORA_SLIDER_PVALUE &
                 `OR_FLAT` >= input$TSA_ORA_SLIDER_OR]
      setorder(dt, `BH_P_VALUE_FLAT`)
      cols_numeric <- c("ORA_SCORE_FLAT", "OR_FLAT", "BH_P_VALUE_FLAT")
    }
    show_DT(
      dt,
      cols_to_show = colnames(dt),
      cols_numeric = cols_numeric,
      table_title = "Over-representation Table"
    )
  })
}

plot_TSA_ORA <- function(
  input
) {
  renderPlot({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_ORA_CATEGORY_CHOICE, input$TSA_ORA_TYPE_CHOICE)
    obj <- DATASETS_COMBINED[[input$TSA_DATASET_CHOICE]]
    req(obj)
    if (input$TSA_ORA_TYPE_CHOICE == "Up") {
      reg <- "UP"
    } else if (input$TSA_ORA_TYPE_CHOICE == "Down") {
      reg <- "DOWN"
    } else if ( input$TSA_ORA_TYPE_CHOICE == "Stable") {
     reg <- "FLAT"
    }
    p <- PlotORA(
      object = obj,
      subID = input$TSA_TISSUE_CHOICE,
      category = input$TSA_ORA_CATEGORY_CHOICE,
      regulation = reg,
      max_terms_show = 20,
      OR_threshold = input$TSA_ORA_SLIDER_OR,
      p_value_threshold = min(0.05, input$TSA_ORA_SLIDER_PVALUE),
      stringent = FALSE
    )
    return(p)
    #if (is.null(p)) {
    #  
    #}
    #p +
    #  ggtitle("Over-representation Plot")
  })
}

