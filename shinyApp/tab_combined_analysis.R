tab_combined_analysis <- tabPanel(
  title = "Combined Analysis",
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "TCA_DATASET_CHOICE",
        label = "Dataset",
        choices =  c("TMS FACS Data", "TMS Droplet Data", "Calico Data")#, "TMS FACS Data (Sub)")
      ),
      hr(),
      conditionalPanel(
        condition = "input.active_TCA_panel == 'TCA_SUMMARY'",
        selectInput(
          inputId = "TCA_SUMMARY_TYPE_CHOICE",
          label = "Type of Information",
          choices = c("Table", "ORA Plot")
        ),
        selectInput(
          inputId = "TCA_SUMMARY_TABLE_CHOICE",
          label = "Category",
          choices = c("KEGG Pathways", "GO Terms", "LR Genes", "Tissues", "Cell-Type Families")
        ),
        selectInput(
          inputId = "TCA_ORA_TYPE_CHOICE",
          label = "ORA Regulation",
          choices = c("Up", "Down", "Stable")
        )
      ),
      conditionalPanel(
        condition = "input.active_TCA_panel == 'TCA_INTERACTION_ANALYSIS'",
        selectInput(
          inputId = "TCA_CCI_DETAILS_CHOICE",
          label = "Type of Information",
          choices = c("CCI Table", "Volcano Plot", "Score Plot")
        ),
        hr(),
        uiOutput("TCA_TISSUE_CHOICE"),
        uiOutput("TCA_REGULATION_CHOICE")
        #uiOutput("TCA_EMITTER_CHOICE"),
        #uiOutput("TCA_RECEIVER_CHOICE"),
        #sliderInput(
        #  inputId = "TCA_SLIDER_PVALUE",
        #  label = "P-value Threshold",
        #  min = 0, 
        #  max = 1,
        #  value = 1
        #),
        #uiOutput("TCA_SLIDER_LOG2FC")
      ),
      width = 3
    ),
    mainPanel(
      titlePanel(htmlOutput("TCA_TITLE")),
      tabsetPanel(
        type = "tabs",
        tabPanel(
          title = "Summary Tables",
          fluidRow(
            column(width = 12, htmlOutput("TCA_SUMMARY_INTRO"), style = "padding:50px"),
            column(width = 12, uiOutput("TCA_SUMMARY_TYPE"))
          ),
          value = "TCA_SUMMARY"
        ),
        tabPanel(
          title = "Detailed Interactions",
          fluidRow(
            column(width = 12, htmlOutput("TCA_CCI_INTRO"), style = "padding:50px"),
            column(width = 12, uiOutput("TCA_CCI_DETAILS")),
            column(width = 12, uiOutput("TCA_CCI_TEXTOUTPUT"), style = "padding:10px")
          ),
          value = "TCA_INTERACTION_ANALYSIS"
        ),
        id = "active_TCA_panel"
      )
    )
  ))