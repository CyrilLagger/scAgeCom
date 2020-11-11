####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk - October 2020
##
####################################################
##

## Load libraries ####

library(shiny)
library(data.table)
library(DT)
library(ggplot2)
library(shinyWidgets)

## Global options ####

#to display NA in DT
options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))

## Load full data ####

LR6db <- readRDS("data/LR6db_curated.rds")
cols_to_show_LR6db <- c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3",
                        "DATABASE", "SOURCE", "ANNOTATION", "FAMILY", "SUBFAMILY")

DATASETS_light <- readRDS("data/a4_data_results_all_cases.rds")
names(DATASETS_light) <- c("Calico Data", "TMS Droplet Data", "TMS FACS Data")

## Process data ####

# Add columns and change colname

scdiffcom_data_filtered <- rbindlist(
  lapply(
    scdiffcom_data,
    function(dataset) {
      rbindlist(
        lapply(
          dataset,
          function(tiss) {
            scDiffCom:::get_cci_table_filtered(tiss)
          }
        ),
        use.names = TRUE,
        idcol = "TISSUE",
        fill = TRUE
      )
    }
  ),
  use.names = TRUE,
  idcol = "DATASET"
)



scdiffcom_data_filtered[, c("L_TCT", "R_TCT", "LOG2FC") := list(
  paste(TISSUE, L_CELLTYPE, sep = ": "),
  paste(TISSUE, R_CELLTYPE, sep = ": "),
  LOGFC*log2(exp(1))
)]
dataset_name_conversion <- data.table(
  old_names = unique(scdiffcom_data_filtered$DATASET),
  Dataset = c("Calico Data", "TMS Droplet Data", "TMS FACS Data")
)
scdiffcom_data_filtered[dataset_name_conversion, on = "DATASET==old_names", Dataset := i.Dataset]
setnames(
  scdiffcom_data_filtered,
  old = c("L_CELLTYPE", "R_CELLTYPE", "BH_PVAL_DIFF", "LR_NAME"),
  new = c("Emitter cell type", "Receiver cell type", "Adj. P-value", "LR_GENES")
)
cols_to_show_scdiffcom <- c("TISSUE", "LR_GENES", "Emitter cell type", "Receiver cell type",
                          "LOG2FC", "Adj. P-value", "REGULATION")
cols_numeric_scdiffcom <- c("LOG2FC", "Adj. P-value")

## Process ORA data as a unique list for each dataset and each tissue ####

scdiffcom_data_ora <-  lapply(
  scdiffcom_data,
  function(dataset) {
      lapply(
        dataset,
        function(tiss) {
          scDiffCom:::get_ora_tables(tiss)
        }
      )
  }
)



scdiffcom_data_ora <-  lapply(
  scdiffcom_data_ora,
  function(dataset) {
    lapply(
      dataset,
      function(tiss) {
        temp <- tiss[c("GO", "LR_NAME", "LR_CELLTYPE", "LR_CELL_FAMILY")]
        names(temp) <- c("GO Terms", "Ligand-receptor Genes", "Cell Types", "Cell Families")
        lapply(
          temp,
          function(ORA_dt) {
            setnames(
              ORA_dt,
              old = c("OR_UP", "pval_adjusted_UP", "OR_DOWN", "pval_adjusted_DOWN",
                      "OR_FLAT", "pval_adjusted_FLAT"),
              new = c("Odds Ratio Up", "Adj. P-value Up", "Odds Ratio Down", "Adj. P-value Down",
                      "Odds Ratio Stable", "Adj. P-value Stable")
            )
            return(ORA_dt)
          }
        )
      }
    )
  }
)

#cols_to_show_ora <- c("Value",
#                      "Odds Ratio Up", "Adj. P-value Up", "Odds Ratio Down", "Adj. P-value Down",
#                      "Odds Ratio", "Adj. P-value", "Odds Ratio Stable", "Adj. P-value Stable")
#cols_numeric_ora<- c("Odds Ratio Up", "Adj. P-value Up", "Odds Ratio Down", "Adj. P-value Down",
#                     "Odds Ratio", "Adj. P-value", "Odds Ratio Stable", "Adj. P-value Stable")

## Utility functions ####

show_DT <- function(
  data,
  cols_to_show,
  cols_numeric
) {
  res <- DT::datatable(
    data[, ..cols_to_show],
    options = list(
      pageLength = 10
    )
  ) 
  if(!is.null(cols_numeric)) {
    res <- DT::formatSignif(res, columns = cols_numeric, digits = 3 )
  }
  return(res)
}

show_volcano <- function(
  data
  ) {
  ggplot(data, aes(x = LOG2FC, y = -log10(`Adj. P-value` + 1E-3))) +
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

## Shiny  UI code ####

ui <- fluidPage(
  navbarPage(
    "Age-related changes in mouse intercellular communication. (Alpha version)",
    tabPanel(
      "Main Results",
      sidebarLayout(
        sidebarPanel(
          selectInput(
            inputId = "dataset_selection",
            label = "Dataset",
            choices = c("TMS FACS Data", "TMS Droplet Data", "Calico Data")
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
            choices = c("TMS FACS Data", "TMS Droplet Data", "Calico Data")
          ),
          uiOutput("ora_tissue_selection"),
          uiOutput("ora_category_selection"),
          selectInput(
            inputId = "ora_direction",
            label = "Direction",
            choices = list("Up", "Down", "Stable")
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
            tabPanel("ORA Table", DT::dataTableOutput("ora_table")),
            tabPanel("ORA Plots", plotOutput("ora_plots")),
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

## Shiny server code ####

server <- function(input, output) {
  output$tissue_selection <- renderUI({
    pickerInput(
      inputId = "tissue_select",
      "Tissue",
      choices = sort(unique(scdiffcom_data_filtered[Dataset %in% input$dataset_selection][["TISSUE"]])),
      selected = sort(unique(scdiffcom_data_filtered[Dataset %in% input$dataset_selection][["TISSUE"]])),
      options = list(`actions-box` = TRUE),
      multiple = T
    )
  })
  output$ligand_cell_selection <- renderUI({
    pickerInput(
      inputId = "L_cell_select",
      "Emitter cell Type",
      choices = sort(unique(scdiffcom_data_filtered[Dataset %in% input$dataset_selection & TISSUE %in% input$tissue_select, L_TCT])),
      selected = sort(unique(scdiffcom_data_filtered[Dataset %in% input$dataset_selection & TISSUE %in% input$tissue_select, L_TCT])),
      options = list(`actions-box` = TRUE),
      multiple = T
    )
  })
  output$receptor_cell_selection <- renderUI({
    pickerInput(
      inputId = "R_cell_select",
      "Receiver cell Type",
      choices = sort(unique(scdiffcom_data_filtered[Dataset %in% input$dataset_selection & TISSUE %in% input$tissue_select, R_TCT])),
      selected = sort(unique(scdiffcom_data_filtered[Dataset %in% input$dataset_selection & TISSUE %in% input$tissue_select, R_TCT])),
      options = list(`actions-box` = TRUE),
      multiple = T
    )
  })
  output$log2fc_slider <- renderUI({
    sliderInput(
      inputId = "slider_log2fc",
      label = "LOG2FC filter",
      min = 0, 
      max = min(max(ceiling(abs(scdiffcom_data_filtered[Dataset %in% input$dataset_selection][["LOG2FC"]]))), 12),
      value = 0,
      step = 0.01
    )
  })
  output$scdiffcom_table <- DT::renderDataTable({
    req(input$tissue_select, input$L_cell_select, input$R_cell_select)
    show_DT(
      setorder(
        scdiffcom_data_filtered[
          Dataset %in% input$dataset_selection &
            TISSUE %in% input$tissue_select &
            L_TCT %in% input$L_cell_select &
            R_TCT %in% input$R_cell_select &
            `Adj. P-value` <= input$slider_pvalue &
            abs(LOG2FC) >= input$slider_log2fc
          ],
        -LOG2FC,
        `Adj. P-value`
      ),
      cols_to_show_scdiffcom, cols_numeric_scdiffcom)
  })
  output$scdiffcom_volcano <- renderPlot({
    req(input$tissue_select, input$L_cell_select, input$R_cell_select)
    show_volcano(scdiffcom_data_filtered[
      Dataset %in% input$dataset_selection &
        TISSUE %in% input$tissue_select &
        L_TCT %in% input$L_cell_select &
        R_TCT %in% input$R_cell_select &
        `Adj. P-value` <= input$slider_pvalue &
        abs(LOG2FC) >= input$slider_log2fc
      ])
  })
  output$ora_tissue_selection <- renderUI({
    pickerInput(
      inputId = "ora_tiss_select",
      "Tissue",
      choices = sort(names(scdiffcom_data_ora[[input$ora_dataset_selection]])),
      options = list(`actions-box` = TRUE),
      multiple = FALSE
    )
  })
  output$ora_category_selection <- renderUI({
    pickerInput(
      inputId = "ora_cat_select",
      "Category",
      choices = names(scdiffcom_data_ora[[input$ora_dataset_selection]][[input$ora_tiss_select]]),
      options = list(`actions-box` = TRUE),
      multiple = FALSE
    )
  })
  output$ora_table <- DT::renderDataTable({
    req(input$ora_tiss_select, input$ora_cat_select, input$ora_direction)
    data <- scdiffcom_data_ora[[input$ora_dataset_selection]][[input$ora_tiss_select]][[input$ora_cat_select]]
    if(input$ora_direction == "Up") {
      data <- data[`Odds Ratio Up` >= 1, c("Value", "Odds Ratio Up", "Adj. P-value Up")]
      data <- data[`Adj. P-value Up` <= input$ora_slider_pvalue]
      setorder(data, `Adj. P-value Up`)
      cols_numeric <- c("Odds Ratio Up", "Adj. P-value Up")
    } else if(input$ora_direction == "Down") {
      data <- data[`Odds Ratio Down` >= 1, c("Value", "Odds Ratio Down", "Adj. P-value Down")]
      data <- data[`Adj. P-value Down` <= input$ora_slider_pvalue]
      setorder(data, `Adj. P-value Down`)
      cols_numeric <- c("Odds Ratio Down", "Adj. P-value Down")
    } else if(input$ora_direction == "Either direction") {
      data <- data[`Odds Ratio` >= 1, c("Value", "Odds Ratio", "Adj. P-value")]
      data <- data[`Adj. P-value` <= input$ora_slider_pvalue]
      setorder(data, `Adj. P-value`)
      cols_numeric <- c("Odds Ratio", "Adj. P-value")
    } else if(input$ora_direction == "Stable") {
      data <- data[`Odds Ratio Stable` >= 1, c("Value", "Odds Ratio Stable", "Adj. P-value Stable")]
      data <- data[`Adj. P-value Stable` <= input$ora_slider_pvalue]
      setorder(data, `Adj. P-value Stable`)
      cols_numeric <- c("Odds Ratio Stable", "Adj. P-value Stable")
    }
    show_DT(
      data,
      cols_to_show = colnames(data),
      cols_numeric = cols_numeric
      )
  })
  output$ora_plots <- renderPlot({
    req(input$ora_tiss_select, input$ora_cat_select, input$ora_direction)
    scDiffCom:::plot_ORA(
      object = obj,category = , ORA_type = , max_value = , OR_cutoff = , pval_cutoff = a
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
      show_DT(
      setorder(
        LR6db,
        LIGAND_1, LIGAND_2, RECEPTOR_1, RECEPTOR_2, RECEPTOR_3
      ),
      #LR6db,
      cols_to_show_LR6db, NULL)
  })
}

## Shiny final call ####

shinyApp(
  ui = ui,
  server = server
)
