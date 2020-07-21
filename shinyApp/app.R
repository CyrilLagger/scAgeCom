library(shiny)
library(data.table)
library(DT)
library(ggplot2)
library(shinyWidgets)

options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))

diffcom_data <- readRDS("data/analysis_4_data_diffcom_filter_new.rds")
diffcom_data <- lapply(diffcom_data, function(data) {
  data[, L_TCT := paste(TISSUE, L_CELLTYPE, sep = "_")]
  data[, R_TCT := paste(TISSUE, R_CELLTYPE, sep = "_")]
  data[, LOG2FC := LR_LOGFC*log2(exp(1))]
  return(data)
})
cols_to_show_diffcom <- c("TISSUE", "L_GENE", "R_GENE", "L_CELLTYPE", "R_CELLTYPE",
                          "LOG2FC", "BH_PVAL_DIFF", "CASE_TYPE")
cols_numeric_diffcom <- c("LOG2FC", "BH_PVAL_DIFF")

ora_data <- readRDS("data/analysis_4_data_ora.rds")
cols_to_show_ora <- c("Tissue", "Category", "Value", "OR", "pval_adjusted", "OR_UP", "pval_adjusted_UP",
                      "OR_DOWN", "pval_adjusted_DOWN")
cols_numeric_ora<- c("OR", "pval_adjusted", "OR_UP", "pval_adjusted_UP",
                     "OR_DOWN", "pval_adjusted_DOWN")



TISSUE_DATASET <- list(
  calico = c("Kidney", "Lung", "Spleen"),
  calico_sub = c("Kidney", "Lung", "Spleen"),
  tms_facs = c("Aorta", "BAT", "Bladder", "Brain_Myeloid",
               "Brain_Non-Myeloid", "Diaphragm", "GAT",
               "Heart", "Kidney", "Large_Intestine",
               "Limb_Muscle", "Liver", "Lung", "Mammary_Gland",
               "Marrow", "MAT", "Pancreas", "SCAT", "Skin",
               "Spleen", "Thymus", "Tongue", "Trachea"),
  tms_droplet = c("Bladder", "Heart_and_Aorta", "Kidney",
                  "Limb_Muscle", "Liver", "Lung",
                  "Mammary_Gland", "Marrow", "Spleen",
                  "Thymus", "Tongue")
)

show_DT <- function(data, cols_to_show, cols_numeric) {
  DT::datatable(data[, ..cols_to_show],
                options = list(
                  pageLength = 10
                )
  ) %>%
    DT::formatSignif(columns = cols_numeric, digits = 3 )
}

show_volcano <- function(data) {
  ggplot(data, aes(x = LOG2FC, y = -log10(BH_PVAL_DIFF + 1E-4))) +
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
            choices = list(
              "Tabula Muris Senis FACS" = "tms_facs",
              "Tabula Muris senis Droplet" = "tms_droplet",
              "Calico 2019" = "calico"
            )
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
          width = 3
        ),
        mainPanel(
          tabsetPanel(
            type = "tabs",
            tabPanel("Table", DT::dataTableOutput("diffcom_table")),
            tabPanel("Volcano Plot", plotOutput("diffcom_volcano"))
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
            choices = list(
              "Tabula Muris Senis FACS" = "tms_facs",
              "Tabula Muris senis Droplet" = "tms_droplet",
              "Calico 2019" = "calico"
            )
          ),
          uiOutput("ora_tissue_selection"),
          uiOutput("ora_category_selection"),
          selectInput(
            inputId = "ora_direction",
            label = "Direction",
            choices = list("UP", "DOWN", "BOTH")
          ),
          selectInput(
            inputId = "ora_correlation",
            label = "Correlation",
            choices = list("Change", "Stable")
          ),
          width = 3
        ),
        mainPanel(
          tabsetPanel(
            type = "tabs",
            tabPanel("ORA", DT::dataTableOutput("ora_table")),
            tabPanel("HeatMaps")
          )
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
      choices = TISSUE_DATASET[[input$dataset_selection]],
      selected = TISSUE_DATASET[[input$dataset_selection]],
      options = list(`actions-box` = TRUE),
      multiple = T
    )
  })
  output$ligand_cell_selection <- renderUI({
    pickerInput(
      inputId = "L_cell_select",
      "Ligand Cell Type",
      choices = unique(diffcom_data[[input$dataset_selection]][TISSUE %in% input$tissue_select, L_TCT]),
      selected = unique(diffcom_data[[input$dataset_selection]][TISSUE %in% input$tissue_select, L_TCT]),
      options = list(`actions-box` = TRUE),
      multiple = T
    )
  })
  output$receptor_cell_selection <- renderUI({
    pickerInput(
      inputId = "R_cell_select",
      "Receptor Cell Type",
      choices = unique(diffcom_data[[input$dataset_selection]][TISSUE %in% input$tissue_select, R_TCT]),
      selected = unique(diffcom_data[[input$dataset_selection]][TISSUE %in% input$tissue_select, R_TCT]),
      options = list(`actions-box` = TRUE),
      multiple = T
    )
  })
  output$diffcom_table <- DT::renderDataTable({
    req(input$tissue_select, input$L_cell_select, input$R_cell_select)
    show_DT(
      setorder(
        diffcom_data[[input$dataset_selection]][TISSUE %in% input$tissue_select &
                                                  L_TCT %in% input$L_cell_select &
                                                  R_TCT %in% input$R_cell_select &
                                                  BH_PVAL_DIFF <= input$slider_pvalue],
        -LOG2FC,
        BH_PVAL_DIFF
      ),
            cols_to_show_diffcom, cols_numeric_diffcom)
  })
  output$diffcom_volcano <- renderPlot({
    req(input$tissue_select, input$L_cell_select, input$R_cell_select)
    show_volcano(diffcom_data[[input$dataset_selection]][TISSUE %in% input$tissue_select &
                                                           L_TCT %in% input$L_cell_select &
                                                           R_TCT %in% input$R_cell_select &
                                                           BH_PVAL_DIFF <= input$slider_pvalue])
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
    if(input$ora_direction == "UP") {
      if(input$ora_correlation == "Change") {
        data <- data[OR_UP >= 1, c("Tissue", "Category", "Value", "OR_UP", "pval_adjusted_UP")]
      } else {
        data <- data[OR_UP <= 1, c("Tissue", "Category", "Value", "OR_UP", "pval_adjusted_UP")]
      }
      setorder(data, pval_adjusted_UP)
      cols_numeric <- c("OR_UP", "pval_adjusted_UP")
    } else if(input$ora_direction == "DOWN") {
      if(input$ora_correlation == "Change") {
        data <- data[OR_DOWN >= 1, c("Tissue", "Category", "Value", "OR_DOWN", "pval_adjusted_DOWN")]
      } else {
        data <- data[OR_DOWN <= 1, c("Tissue", "Category", "Value", "OR_DOWN", "pval_adjusted_DOWN")]
      }
      setorder(data, pval_adjusted_DOWN)
      cols_numeric <- c("OR_DOWN", "pval_adjusted_DOWN")
    } else {
      if(input$ora_correlation == "Change") {
        data <- data[OR >= 1, c("Tissue", "Category", "Value", "OR", "pval_adjusted")]
      } else {
        data <- data[OR <= 1, c("Tissue", "Category", "Value", "OR", "pval_adjusted")]
      }
      setorder(data, pval_adjusted)
      cols_numeric <- c("OR", "pval_adjusted")
    }
    show_DT(
      data,
      cols_to_show = colnames(data),
      cols_numeric = cols_numeric
      )
  })
}

shinyApp(ui = ui, server = server)
