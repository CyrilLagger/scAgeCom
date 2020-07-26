library(shiny)
library(data.table)
library(DT)
library(ggplot2)
library(shinyWidgets)

# To Do: the file is becoming longer and longer... I will reorganize it in different subfiles. 

#to display NA in DT
options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))

#Load data

#LR database
LRdb <- readRDS("data/LRdb_data.rds")
setnames(
  LRdb,
  old = c("GENESYMB_L", "GENESYMB_R", "SYMB_LR",
          "CONF_L", "TYPE_L", "CONF_R", "TYPE_R"),
  new = c("Ligand", "Receptor", "LR",
          "Ortho. conf. L", "Ortho. type L", "Ortho. conf. R", "Ortho. type R")
)
cols_to_show_LRdb <- c("Ligand", "Receptor", "LR", "source_scsr", "source_cpdb", "source_nichenet", "source_sctensor",
                       "Ortho. conf. L", "Ortho. type L", "Ortho. conf. R", "Ortho. type R")

#scDiffCom data
#diffcom_data <- readRDS("data/analysis_4_data_diffcom_filter_new.rds")
diffcom_data <- readRDS("data/analysis_4_data_diffcom_filter_less_stringent.rds")
diffcom_data <- lapply(diffcom_data, function(data) {
  data[, L_TCT := paste(TISSUE, L_CELLTYPE, sep = "_")]
  data[, R_TCT := paste(TISSUE, R_CELLTYPE, sep = "_")]
  data[, LOG2FC := LR_LOGFC*log2(exp(1))]
  data[, Direction := ifelse(
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
    data,
    old = c("L_GENE", "R_GENE", "L_CELLTYPE", "R_CELLTYPE", "BH_PVAL_DIFF"),
    new = c("Ligand", "Receptor", "Emitter cell type", "Receiver cell type", "Adj. P-value")
  )
  return(data)
})
cols_to_show_diffcom <- c("TISSUE", "Ligand", "Receptor", "Emitter cell type", "Receiver cell type",
                          "LOG2FC", "Adj. P-value", "Direction")
cols_numeric_diffcom <- c("LOG2FC", "Adj. P-value")

#ORA data
#ora_data <- readRDS("data/analysis_4_data_ora.rds")
ora_data <- readRDS("data/analysis_4_data_ora_less_stringent.rds")
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

#List of tissue for each dataset
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

#utility functions
show_DT <- function(data, cols_to_show, cols_numeric) {
  DT::datatable(data[, ..cols_to_show],
                options = list(
                  pageLength = 10
                )
  ) %>%
    DT::formatSignif(columns = cols_numeric, digits = 3 )
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
          uiOutput("log2fc_slider"),
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
            inputId = "LRdb_select",
            "Database",
            choices = c("scsr", "cpdb", "nichenet", "sctensor"),
            selected = c("scsr", "cpdb", "nichenet", "sctensor"),
            options = list(`actions-box` = TRUE),
            multiple = T
          ),
          width = 2
        ),
        mainPanel(
          DT::dataTableOutput("LRdb_table")
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
      "Emitter cell Type",
      choices = unique(diffcom_data[[input$dataset_selection]][TISSUE %in% input$tissue_select, L_TCT]),
      selected = unique(diffcom_data[[input$dataset_selection]][TISSUE %in% input$tissue_select, L_TCT]),
      options = list(`actions-box` = TRUE),
      multiple = T
    )
  })
  output$receptor_cell_selection <- renderUI({
    pickerInput(
      inputId = "R_cell_select",
      "Receiver cell Type",
      choices = unique(diffcom_data[[input$dataset_selection]][TISSUE %in% input$tissue_select, R_TCT]),
      selected = unique(diffcom_data[[input$dataset_selection]][TISSUE %in% input$tissue_select, R_TCT]),
      options = list(`actions-box` = TRUE),
      multiple = T
    )
  })
  output$log2fc_slider <- renderUI({
    sliderInput(
      inputId = "slider_log2fc",
      label = "LOG2FC filter",
      min = 0, 
      max = max(ceiling(abs(diffcom_data[[input$dataset_selection]]$LOG2FC))),
      value = 0,
      step = 0.01
    )
  })
  output$diffcom_table <- DT::renderDataTable({
    req(input$tissue_select, input$L_cell_select, input$R_cell_select)
    show_DT(
      setorder(
        diffcom_data[[input$dataset_selection]][TISSUE %in% input$tissue_select &
                                                  L_TCT %in% input$L_cell_select &
                                                  R_TCT %in% input$R_cell_select &
                                                  `Adj. P-value` <= input$slider_pvalue &
                                                  abs(LOG2FC) >= input$slider_log2fc],
        -LOG2FC,
        `Adj. P-value`
      ),
            cols_to_show_diffcom, cols_numeric_diffcom)
  })
  output$diffcom_volcano <- renderPlot({
    req(input$tissue_select, input$L_cell_select, input$R_cell_select)
    show_volcano(diffcom_data[[input$dataset_selection]][TISSUE %in% input$tissue_select &
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
  output$LRdb_table <- DT::renderDataTable({
    req(input$LRdb_select)
    if(!("scsr" %in% input$LRdb_select)) {
      LRdb <- LRdb[scsr == FALSE]
    }
    if(!("cpdb" %in% input$LRdb_select)) {
      LRdb <- LRdb[cpdb == FALSE]
    }
    if(!("nichenet" %in% input$LRdb_select)) {
      LRdb <- LRdb[nichenet == FALSE]
    }
    if(!("sctensor" %in% input$LRdb_select)) {
      LRdb <- LRdb[sctensor == FALSE]
    }
    show_DT(
      setorder(
        LRdb,
        LR
      ),
      cols_to_show_LRdb, NULL)
  })
}

shinyApp(ui = ui, server = server)
