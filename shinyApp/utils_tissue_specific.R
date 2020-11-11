choose_TSA_tissue <- function(
  input
) {
  renderUI({
    choices <- sort(names(DATASETS_light[[input$TSA_DATASET_CHOICE]]))
    pickerInput(
      inputId = "TSA_TISSUE_CHOICE",
      label = "Tissue",
      choices = choices,
      options = list(`actions-box` = TRUE),
      multiple = FALSE
    )
  })
}

choose_TSA_emitter <- function(
  input
) {
  renderUI({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE)
    obj <- DATASETS_light[[input$TSA_DATASET_CHOICE]][[input$TSA_TISSUE_CHOICE]]
    req(obj)
    choices <- sort(unique(scDiffCom:::get_cci_table_filtered(obj)[["Emitter Cell Type"]]))
    pickerInput(
      inputId = "TSA_EMITTER_CHOICE",
      label = "Emitter Cell Type",
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
    obj <- DATASETS_light[[input$TSA_DATASET_CHOICE]][[input$TSA_TISSUE_CHOICE]]
    req(obj)
    choices <- sort(unique(scDiffCom:::get_cci_table_filtered(obj)[["Receiver Cell Type"]]))
    pickerInput(
      inputId = "TSA_RECEIVER_CHOICE",
      label = "Receiver Cell Type",
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
    obj <- DATASETS_light[[input$TSA_DATASET_CHOICE]][[input$TSA_TISSUE_CHOICE]]
    req(obj)
    max_val <- ceiling(max(scDiffCom:::get_cci_table_filtered(obj)[["LOG2FC"]]))
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

choose_TSA_ORA_category <- function(
  input
) {
  renderUI({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE)
    obj <- DATASETS_light[[input$TSA_DATASET_CHOICE]][[input$TSA_TISSUE_CHOICE]]
    req(obj)
    choices <- names(scDiffCom:::get_ora_tables(obj))
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
    obj <- DATASETS_light[[input$TSA_DATASET_CHOICE]][[input$TSA_TISSUE_CHOICE]]
    req(obj)
    ora_table <- scDiffCom:::get_ora_tables(obj)[[input$TSA_ORA_CATEGORY_CHOICE]]
    req(ora_table)
    if(input$TSA_ORA_TYPE_CHOICE == "Up") {
      max_val <- ceiling(max(ora_table[["Odds Ratio Up"]]))
    } else if(input$TSA_ORA_TYPE_CHOICE == "Down") {
      max_val <- ceiling(max(ora_table[["Odds Ratio Down"]]))
    } else if(input$TSA_ORA_TYPE_CHOICE == "Stable") {
      max_val <- ceiling(max(ora_table[["Odds Ratio Stable"]]))
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

get_TSA_interaction_table <- function(
  input
) {
  DT::renderDataTable({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_EMITTER_CHOICE, input$TSA_RECEIVER_CHOICE,
        input$TSA_SLIDER_PVALUE, input$TSA_SLIDER_LOG2FC)
    obj <- DATASETS_light[[input$TSA_DATASET_CHOICE]][[input$TSA_TISSUE_CHOICE]]
    req(obj)
    dt <- scDiffCom:::get_cci_table_filtered(obj)
    req(dt)
    dt <- dt[
      `Emitter Cell Type` %in% input$TSA_EMITTER_CHOICE &
        `Receiver Cell Type` %in% input$TSA_RECEIVER_CHOICE &
        `Adj. P-Value` <= input$TSA_SLIDER_PVALUE &
        abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
        ]
    setorder(
      dt,
      -LOG2FC,
      `Adj. P-Value`
    )
    show_DT(
      dt,
      cols_to_show_DATA,
      cols_numeric_DATA
      )
  })
}

plot_TSA_VOLCANO <- function(
  input
) {
  renderPlot({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_EMITTER_CHOICE, input$TSA_RECEIVER_CHOICE,
        input$TSA_SLIDER_PVALUE, input$TSA_SLIDER_LOG2FC)
    obj <- DATASETS_light[[input$TSA_DATASET_CHOICE]][[input$TSA_TISSUE_CHOICE]]
    req(obj)
    dt <- scDiffCom:::get_cci_table_filtered(obj)
    req(dt)
    dt <- dt[
      `Emitter Cell Type` %in% input$TSA_EMITTER_CHOICE &
        `Receiver Cell Type` %in% input$TSA_RECEIVER_CHOICE &
        `Adj. P-Value` <= input$TSA_SLIDER_PVALUE &
        abs(LOG2FC) >= input$TSA_SLIDER_LOG2FC
      ]
    show_volcano(dt)
  })
}

get_TSA_ORA_table <- function(
  input
) {
  DT::renderDataTable({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_ORA_CATEGORY_CHOICE, input$TSA_ORA_TYPE_CHOICE)
    obj <- DATASETS_light[[input$TSA_DATASET_CHOICE]][[input$TSA_TISSUE_CHOICE]]
    req(obj)
    dt <- scDiffCom:::get_ora_tables(obj)[[input$TSA_ORA_CATEGORY_CHOICE]]
    req(dt)
    if(input$TSA_ORA_TYPE_CHOICE == "Up") {
      dt <- dt[`Odds Ratio Up` >= 1, c("Value", "Odds Ratio Up", "Adj. P-Value Up")]
      dt <- dt[`Adj. P-Value Up` <= input$TSA_ORA_SLIDER_PVALUE &
                 `Odds Ratio Up` >= input$TSA_ORA_SLIDER_OR]
      setorder(dt, `Adj. P-Value Up`)
      cols_numeric <- c("Odds Ratio Up", "Adj. P-Value Up")
    } else if(input$TSA_ORA_TYPE_CHOICE == "Down") {
      dt <- dt[`Odds Ratio Down` >= 1, c("Value", "Odds Ratio Down", "Adj. P-Value Down")]
      dt <- dt[`Adj. P-Value Down` <= input$TSA_ORA_SLIDER_PVALUE &
                 `Odds Ratio Down` >= input$TSA_ORA_SLIDER_OR]
      setorder(dt, `Adj. P-Value Down`)
      cols_numeric <- c("Odds Ratio Down", "Adj. P-Value Down")
    } else if(input$TSA_ORA_TYPE_CHOICE == "Stable") {
      dt <- dt[`Odds Ratio Stable` >= 1, c("Value", "Odds Ratio Stable", "Adj. P-Value Stable")]
      dt <- dt[`Adj. P-Value Stable` <= input$TSA_ORA_SLIDER_PVALUE &
                 `Odds Ratio Stable` >= input$TSA_ORA_SLIDER_OR]
      setorder(dt, `Adj. P-Value Stable`)
      cols_numeric <- c("Odds Ratio Stable", "Adj. P-Value Stable")
    }
    show_DT(
      dt,
      cols_to_show = colnames(dt),
      cols_numeric = cols_numeric
    )
  })
}

plot_TSA_ORA <- function(
  input
) {
  renderPlot({
    req(input$TSA_DATASET_CHOICE, input$TSA_TISSUE_CHOICE, input$TSA_ORA_CATEGORY_CHOICE, input$TSA_ORA_TYPE_CHOICE)
    obj <- DATASETS_light[[input$TSA_DATASET_CHOICE]][[input$TSA_TISSUE_CHOICE]]
    req(obj)
    if (input$TSA_ORA_TYPE_CHOICE == "Up") {
      OR_val <- "Odds Ratio Up"
      pval_val <- "Adj. P-Value Up"
      ORA_score_val <- "ORA_score_UP"
    } else if (input$TSA_ORA_TYPE_CHOICE == "Down") {
      OR_val <- "Odds Ratio Down"
      pval_val <- "Adj. P-Value Down"
      ORA_score_val <- "ORA_score_DOWN"
    } else if ( input$TSA_ORA_TYPE_CHOICE == "Stable") {
      OR_val <- "Odds Ratio Stable"
      pval_val <- "Adj. P-Value Stable"
      ORA_score_val <- "ORA_score_FLAT"
    }
    scDiffCom:::plot_ORA(
      object = obj,
      category = input$TSA_ORA_CATEGORY_CHOICE,
      OR_val,
      pval_val,
      ORA_score_val,
      max_value = 20,
      OR_cutoff = input$TSA_ORA_SLIDER_OR,
      pval_cutoff = input$TSA_ORA_SLIDER_PVALUE
    )
  })
}


