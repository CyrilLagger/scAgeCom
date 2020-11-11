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
      conditionalPanel(
        condition = "input.active_TSA_panel == 'TSA_INTERACTION_TABLE' || input.active_TSA_panel == 'TSA_VOLCANO_PLOT' ",
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
        condition = "input.active_TSA_panel == 'TSA_ORA_TABLE' || input.active_TSA_panel == 'TSA_ORA_PLOT' ",
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
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        type = "tabs",
        tabPanel(
          title = "Interaction Table",
          DT::dataTableOutput("TSA_INTERACTION_TABLE"),
          value = "TSA_INTERACTION_TABLE"
        ),
        tabPanel(
          title = "Volcano Plot",
          plotOutput("TSA_VOLCANO_PLOT", height = "600px"),
          value = "TSA_VOLCANO_PLOT"
        ),
        tabPanel(
          title = "ORA Table",
          DT::dataTableOutput("TSA_ORA_TABLE"),
          value = "TSA_ORA_TABLE"
        ),
        tabPanel(
          title = "ORA Plot",
          plotOutput("TSA_ORA_PLOT", height = "800px"),
          value = "TSA_ORA_PLOT"
        ),
        tabPanel(
          title = "Network"
        ),
        id = "active_TSA_panel"
      )
    )
  )
)