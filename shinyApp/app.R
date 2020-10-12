

## Load libraries ####

library(shiny)
library(data.table)
library(DT)
library(ggplot2)
library(shinyWidgets)

## Global options ####

#to display NA in DT
options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))

## Load data ####

#LR database
LR6db <- readRDS("data/LR6db_curated.rds")
cols_to_show_LR6db <- c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3",
                        "DATABASE", "SOURCE", "ANNOTATION", "FAMILY", "SUBFAMILY")

#scDiffCom results
scdiffcom_data <- readRDS("data/a4_data_diffcom_all_filtered.rds")
lapply(scdiffcom_data, function(data) {
  data$results[, L_TCT := paste(TISSUE, L_CELLTYPE, sep = "_")]
  data$results[, R_TCT := paste(TISSUE, R_CELLTYPE, sep = "_")]
  data$results[, LOG2FC := LOGFC*log2(exp(1))]
  data$results[, Direction := ifelse(
    CASE_TYPE == "TFTD", "Down (disappears)",
    ifelse(
      CASE_TYPE == "TTTD", "Down",
      ifelse(
        CASE_TYPE == "TTTU", "Up",
        ifelse(
          CASE_TYPE == "FTTU", "Up (appears)",
          "Flat"
        )
      )
    )
  )]
  setnames(
    data$results,
    old = c("L_CELLTYPE", "R_CELLTYPE", "BH_PVAL_DIFF"),
    new = c("Emitter cell type", "Receiver cell type", "Adj. P-value")
  )
  return(data)
})
#cols_to_show_scdiffcom <- c("TISSUE", "LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3", "Emitter cell type", "Receiver cell type",
#                            "LOG2FC", "Adj. P-value", "Direction")

cols_to_show_scdiffcom <- c("TISSUE", "LR_NAME", "Emitter cell type", "Receiver cell type",
                          "LOG2FC", "Adj. P-value", "Direction")
cols_numeric_scdiffcom <- c("LOG2FC", "Adj. P-value")
#names(scdiffcom_data) <- lapply(scdiffcom_data, function(i) {i$id})

#ORA data
ora_data <- readRDS("data/a4_data_ora.rds")
ora_data <- lapply(ora_data, function(data) {
  setnames(
    data,
    old = c("OR", "pval_adjusted", "OR_UP", "pval_adjusted_UP", "OR_DOWN", "pval_adjusted_DOWN"),
    new = c("Odds Ratio", "Adj. P-value", "Odds Ratio Up", "Adj. P-value Up", "Odds Ratio Down", "Adj. P-value Down")
  )
  return(data)
})
cols_to_show_ora <- c("Tissue", "Category", "Value", "Odds Ratio",
                      "Adj. P-value", "Odds Ratio Up", "Adj. P-value Up",
                      "Odds Ratio Down", "Adj. P-value Down")
cols_numeric_ora<- c("Odds Ratio", "Adj. P-value", "Odds Ratio Up", "Adj. P-value Up",
                     "Odds Ratio Down", "Adj. P-value Down")

## Utility functions ####

show_DT <- function(data, cols_to_show, cols_numeric) {
  res <- DT::datatable(data[, ..cols_to_show],
                options = list(
                  pageLength = 10
                )
  ) 
  if(!is.null(cols_numeric)) {
    res <- DT::formatSignif(res, columns = cols_numeric, digits = 3 )
  }
  return(res)
}

show_volcano <- function(data) {
  ggplot(data, aes(x = LOG2FC, y = -log10(`Adj. P-value` + 1E-4))) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05)) +
    geom_vline(xintercept = log2(1.1)) +
    geom_vline(xintercept = -log2(1.1)) +
    geom_vline(xintercept = log2(1.5)) +
    geom_vline(xintercept = -log2(1.5)) +
    xlab(expression(paste(Log[2], "FC"))) +
    ylab(expression(paste(-Log[10], " ", p[BH]))) +
    theme(text=element_text(size=20))
}

ui <- fluidPage(
  navbarPage(
    "Age-related changes in mouse intercellular communication. (For testing only!)",
    tabPanel(
      "Main Results",
      sidebarLayout(
        sidebarPanel(
          selectInput(
            inputId = "dataset_selection",
            label = "Dataset",
            choices = names(scdiffcom_data)
          ),
          uiOutput("tissue_selection"),
          uiOutput("ligand_cell_selection"),
          uiOutput("receptor_cell_selection"),
          sliderInput(
            inputId = "slider_pvalue",
            label = "P-value filter",
            min = 0, 
            max = 1,
            value = 1
          ),
          uiOutput("log2fc_slider"),
          width = 3
        ),
        mainPanel(
          tabsetPanel(
            type = "tabs",
            tabPanel("Table", DT::dataTableOutput("scdiffcom_table")),
            tabPanel("Volcano Plot", plotOutput("scdiffcom_volcano"))
          )
        )
      )),
    tabPanel(
      "ORA",
      sidebarLayout(
        sidebarPanel(
          selectInput(
            inputId = "ora_dataset_selection",
            label = "Dataset",
            choices = names(scdiffcom_data)
          ),
          uiOutput("ora_tissue_selection"),
          uiOutput("ora_category_selection"),
          selectInput(
            inputId = "ora_direction",
            label = "Direction",
            choices = list("Up", "Down", "Both")
          ),
          selectInput(
            inputId = "ora_correlation",
            label = "Correlation",
            choices = list("Change", "Stable")
          ),
          sliderInput(
            inputId = "ora_slider_pvalue",
            label = "P-value filter",
            min = 0, 
            max = 1,
            value = 1
          ),
          width = 3
        ),
        mainPanel(
          tabsetPanel(
            type = "tabs",
            tabPanel("ORA", DT::dataTableOutput("ora_table")),
            tabPanel("HeatMaps  (ToDo)")
          )
        )
      )
    ),
    tabPanel(
      "Ligand-Receptor DB",
      sidebarLayout(
        sidebarPanel(
          pickerInput(
            inputId = "LR6db_select",
            "Database",
            choices = c("SCSR", "CELLPHONEDB", "CELLCHAT", "NICHENET", "ICELLNET", "SCTENSOR"),
            selected = c("SCSR", "CELLPHONEDB", "CELLCHAT", "NICHENET", "ICELLNET", "SCTENSOR"),
            options = list(`actions-box` = TRUE),
            multiple = T
          ),
          width = 2
        ),
        mainPanel(
          DT::dataTableOutput("LR6db_table")
        )
      )
    )
  )
)

server <- function(input, output) {
  output$tissue_selection <- renderUI({
    pickerInput(
      inputId = "tissue_select",
      "Tissue",
      choices = scdiffcom_data[[input$dataset_selection]]$tissues,
      selected = scdiffcom_data[[input$dataset_selection]]$tissues,
      options = list(`actions-box` = TRUE),
      multiple = T
    )
  })
  output$ligand_cell_selection <- renderUI({
    pickerInput(
      inputId = "L_cell_select",
      "Emitter cell Type",
      choices = unique(scdiffcom_data[[input$dataset_selection]]$results[TISSUE %in% input$tissue_select, L_TCT]),
      selected = unique(scdiffcom_data[[input$dataset_selection]]$results[TISSUE %in% input$tissue_select, L_TCT]),
      options = list(`actions-box` = TRUE),
      multiple = T
    )
  })
  output$receptor_cell_selection <- renderUI({
    pickerInput(
      inputId = "R_cell_select",
      "Receiver cell Type",
      choices = unique(scdiffcom_data[[input$dataset_selection]]$results[TISSUE %in% input$tissue_select, R_TCT]),
      selected = unique(scdiffcom_data[[input$dataset_selection]]$results[TISSUE %in% input$tissue_select, R_TCT]),
      options = list(`actions-box` = TRUE),
      multiple = T
    )
  })
  output$log2fc_slider <- renderUI({
    sliderInput(
      inputId = "slider_log2fc",
      label = "LOG2FC filter",
      min = 0, 
      max = max(ceiling(abs(scdiffcom_data[[input$dataset_selection]]$results$LOG2FC))),
      value = 0,
      step = 0.01
    )
  })
  output$scdiffcom_table <- DT::renderDataTable({
    req(input$tissue_select, input$L_cell_select, input$R_cell_select)
    show_DT(
      setorder(
        scdiffcom_data[[input$dataset_selection]]$results[TISSUE %in% input$tissue_select &
                                                  L_TCT %in% input$L_cell_select &
                                                  R_TCT %in% input$R_cell_select &
                                                  `Adj. P-value` <= input$slider_pvalue &
                                                  abs(LOG2FC) >= input$slider_log2fc],
        -LOG2FC,
        `Adj. P-value`
      ),
            cols_to_show_scdiffcom, cols_numeric_scdiffcom)
  })
  output$scdiffcom_volcano <- renderPlot({
    req(input$tissue_select, input$L_cell_select, input$R_cell_select)
    show_volcano(scdiffcom_data[[input$dataset_selection]]$results[TISSUE %in% input$tissue_select &
                                                           L_TCT %in% input$L_cell_select &
                                                           R_TCT %in% input$R_cell_select &
                                                           `Adj. P-value` <= input$slider_pvalue &
                                                           abs(LOG2FC) >= input$slider_log2fc])
  })
  output$ora_tissue_selection <- renderUI({
    pickerInput(
      inputId = "ora_tiss_select",
      "Tissue",
      choices = unique(ora_data[[input$ora_dataset_selection]][, Tissue]),
      options = list(`actions-box` = TRUE),
      multiple = FALSE
    )
  })
  output$ora_category_selection <- renderUI({
    pickerInput(
      inputId = "ora_cat_select",
      "Category",
      choices = unique(ora_data[[input$ora_dataset_selection]][Tissue %in% input$ora_tiss_select, Category]),
      options = list(`actions-box` = TRUE),
      multiple = FALSE
    )
  })
  output$ora_table <- DT::renderDataTable({
    req(input$ora_tiss_select, input$ora_cat_select)
    data <- ora_data[[input$ora_dataset_selection]][Tissue %in% input$ora_tiss_select &
                                                  Category %in% input$ora_cat_select]
    if(input$ora_direction == "Up") {
      if(input$ora_correlation == "Change") {
        data <- data[`Odds Ratio Up` >= 1, c("Tissue", "Category", "Value", "Odds Ratio Up", "Adj. P-value Up")]
      } else {
        data <- data[`Odds Ratio Up` <= 1, c("Tissue", "Category", "Value", "Odds Ratio Up", "Adj. P-value Up")]
      }
      data <- data[`Adj. P-value Up` <= input$ora_slider_pvalue]
      setorder(data, `Adj. P-value Up`)
      cols_numeric <- c("Odds Ratio Up", "Adj. P-value Up")
    } else if(input$ora_direction == "Down") {
      if(input$ora_correlation == "Change") {
        data <- data[`Odds Ratio Down` >= 1, c("Tissue", "Category", "Value", "Odds Ratio Down", "Adj. P-value Down")]
      } else {
        data <- data[`Odds Ratio Down` <= 1, c("Tissue", "Category", "Value", "Odds Ratio Down", "Adj. P-value Down")]
      }
      data <- data[`Adj. P-value Down` <= input$ora_slider_pvalue]
      setorder(data, `Adj. P-value Down`)
      cols_numeric <- c("Odds Ratio Down", "Adj. P-value Down")
    } else if(input$ora_direction == "Both") {
      if(input$ora_correlation == "Change") {
        data <- data[`Odds Ratio` >= 1, c("Tissue", "Category", "Value", "Odds Ratio", "Adj. P-value")]
      } else {
        data <- data[`Odds Ratio` <= 1, c("Tissue", "Category", "Value", "Odds Ratio", "Adj. P-value")]
      }
      data <- data[`Adj. P-value` <= input$ora_slider_pvalue]
      setorder(data, `Adj. P-value`)
      cols_numeric <- c("Odds Ratio", "Adj. P-value")
    }
    show_DT(
      data,
      cols_to_show = colnames(data),
      cols_numeric = cols_numeric
      )
  })
  output$LR6db_table <- DT::renderDataTable({
   req(input$LR6db_select)
    LR6db <- LR6db[apply(sapply(input$LR6db_select, function(i) {
      grepl(i, LR6db$DATABASE)
    }),
    MARGIN = 1,
    any
    )]
    
    
    # 
    # if(!("SCSR" %in% input$LR6db_select)) {
    #   LR6db <- LR6db[!grepl("SCSR", DATABASE)]
    # }
    # if(!("CELLPHONEDB" %in% input$LR6db_select)) {
    #   LR6db <- LR6db[!grepl("CELLPHONEDB", DATABASE)]
    # }
    # if(!("NICHENET" %in% input$LR6db_select)) {
    #   LR6db <- LR6db[!grepl("NICHENET", DATABASE)]
    # }
    # if(!("ICELLNET" %in% input$LR6db_select)) {
    #   LR6db <- LR6db[!grepl("ICELLNET", DATABASE)]
    # }
    # if(!("CELLCHAT" %in% input$LR6db_select)) {
    #   LR6db <- LR6db[!grepl("CELLCHAT", DATABASE)]
    # }
    # if(!("SCTENSOR" %in% input$LR6db_select)) {
    #   LR6db <- LR6db[!grepl("SCTENSOR", DATABASE)]
    # }
    show_DT(
      setorder(
        LR6db,
        LIGAND_1, LIGAND_2, RECEPTOR_1, RECEPTOR_2, RECEPTOR_3
      ),
      #LR6db,
      cols_to_show_LR6db, NULL)
  })
}

shinyApp(ui = ui, server = server)
