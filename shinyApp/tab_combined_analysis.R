tab_combined_analysis <- tabPanel(
  title = "Global Analysis",
  fluidRow(
    column(width = 6, titlePanel(htmlOutput("TCA_TITLE")), offset = 3),
    column(width = 6, DT::dataTableOutput("TCA_OVERVIEW_TABLE"), style = "padding-bottom: 50px", offset = 3)
  ),
  tabsetPanel(
    type = "tabs",
    tabPanel(
      title = "Cell-Cell Interactions",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput(
            inputId = "TCA_CCI_DETAILS_CHOICE",
            label = "Choose between Table or Plots",
            choices = c("CCI Table", "Volcano Plot", "Score Plot", "LRI-FC Plot")
          ),
          hr(),
          downloadButton("TCA_DOWNLOAD_TABLE", "Download Table"),
          hr(),
          uiOutput("TCA_TISSUE_CHOICE"),
          #uiOutput("TCA_EMITTER_CHOICE"),
          #uiOutput("TCA_RECEIVER_CHOICE"),
          uiOutput("TCA_LRI_CHOICE")#,
          # sliderInput(
          #   inputId = "TCA_SLIDER_PVALUE",
          #   label = "Filter by Adj. p-value",
          #   min = 0,
          #   max = 1,
          #   value = 1
          # ),
          #uiOutput("TCA_SLIDER_LOG2FC")
        ),
        mainPanel(
          fluidRow(
            column(width = 12, uiOutput("TCA_CCI_DETAILS"), style = "padding:50px"),
            column(width = 12, uiOutput("TCA_CCI_TEXTOUTPUT"))
            #column(width = 12, htmlOutput("TCA_CCI_INTRO"), style = "padding:50px"),
          )
        )
      ),
      value = "TCA_INTERACTION_ANALYSIS"
    ),
    tabPanel(
      title = "Common and Global Changes",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput(
            inputId = "TCA_GLOBAL_DETAILS_CHOICE",
            label = "Type of Information",
            choices = c("Table", "ORA Plot")
          ),
          hr(),
          selectInput(
            inputId = "TCA_GLOBAL_TABLE_CHOICE",
            label = "Category",
            choices = c("KEGG Pathways", "GO Terms", "LR Genes", "Tissues", "Cell-Type Families")
          ),
          selectInput(
            inputId = "TCA_GLOBAL_ORA_TYPE_CHOICE",
            label = "ORA Regulation",
            choices = c("Up", "Down", "Stable")
          )
          #downloadButton("TCA_DOWNLOAD_TABLE", "Download Table"),
        ),
        mainPanel(
          fluidRow(
            #column(width = 12, htmlOutput("TCA_SUMMARY_INTRO"), style = "padding:50px"),
            column(width = 12, uiOutput("TCA_GLOBAL_DETAILS"), style = "padding:50px")
          )
        )
      ),
      value = "TCA_GLOBAL_CHANGES"
    ),
    id = "active_TCA_panel"
  )
)
