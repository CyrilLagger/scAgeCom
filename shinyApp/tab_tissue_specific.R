tab_tissue_specific <- tabPanel(
  title = "Tissue Specific Analysis",
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "TSA_DATASET_CHOICE",
        label = "Dataset",
        choices = c("TMS FACS Data", "TMS Droplet Data", "Calico Data")
      ),
      uiOutput("TSA_TISSUE_CHOICE"),
      hr(),
      conditionalPanel(
        condition = "input.active_TSA_panel == 'TSA_INTERACTION_ANALYSIS'",
        uiOutput("TSA_EMITTER_CHOICE"),
        uiOutput("TSA_RECEIVER_CHOICE"),
        sliderInput(
          inputId = "TSA_SLIDER_PVALUE",
          label = "P-value Threshold",
          min = 0, 
          max = 1,
          value = 1
        ),
        uiOutput("TSA_SLIDER_LOG2FC")
      ),
      conditionalPanel(
        condition = "input.active_TSA_panel == 'TSA_ORA'",
        uiOutput("TSA_ORA_CATEGORY_CHOICE"),
        selectInput(
          inputId = "TSA_ORA_TYPE_CHOICE",
          label = "Regulation",
          choices = list("Up", "Down", "Stable")
        ),
        sliderInput(
          inputId = "TSA_ORA_SLIDER_PVALUE",
          label = "P-value Threshold",
          min = 0, 
          max = 1,
          value = 1
        ),
        uiOutput("TSA_ORA_SLIDER_OR")
      ),
      width = 3,
      style = "position: fixed; width: 20%"
    ),
    mainPanel(
      titlePanel(htmlOutput("TSA_TITLE")),
      tabsetPanel(
        type = "tabs",
        tabPanel(
          title = "Overview",
          fluidRow(
            br(),
            column(12,htmlOutput("TSA_OVERVIEW")),
            br(),
            br(),
            column(12, visNetworkOutput("TSA_NETWORK_PLOT", height = "800px"))
          ),
          value = "TSA_OVERVIEW"
        ),
        tabPanel(
          title = "Interaction Analysis",
          fluidRow(
            br(),
            column(12, plotOutput("TSA_VOLCANO_PLOT", brush = "TSA_VOLCANO_brush", height = "600px")),
            verbatimTextOutput("TSA_VOLCANO_TEXTOUTPUT"),
            br(),
            br(),
            column(12, plotOutput("TSA_SCORES_PLOT",brush = "TSA_SCORES_brush", height = "600px")),
            verbatimTextOutput("TSA_SCORES_TEXTOUTPUT"),
            br(),
            br(),
            column(12, DT::dataTableOutput("TSA_INTERACTION_TABLE"))
          ),
          value = "TSA_INTERACTION_ANALYSIS"
        ),
        tabPanel(
          title = "Over-representation",
          fluidRow(
            br(),
            column(12, plotOutput("TSA_ORA_PLOT", height = "800px")),
            br(),
            br(),
            br(),
            br(),
            column(12, DT::dataTableOutput("TSA_ORA_TABLE"))
          ),
          value = "TSA_ORA"
        ),
        #tabPanel(
        #  title = "Network Graph",
        #  plotOutput("TSA_NETWORK_PLOT", height = "800px"),
        #  value = "TSA_NETWORK_PLOT"
        #),
        id = "active_TSA_panel"
      )
    )
  )
)
