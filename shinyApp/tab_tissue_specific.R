tab_tissue_specific <- tabPanel(
  title = "Tissue-Specific Analysis",
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "TSA_DATASET_CHOICE",
        label = "Dataset",
        choices = c("Calico Data",
                    "TMS Droplet Data Female", "TMS Droplet Data Male", "TMS Droplet Data Mixed",
                    "TMS FACS Data Female", "TMS FACS Data Male", "TMS FACS Data Mixed") #c("TMS FACS Data", "TMS Droplet Data", "Calico Data")
      ),
      uiOutput("TSA_TISSUE_CHOICE"),
      hr(),
      conditionalPanel(
        condition = "input.active_TSA_panel == 'TSA_OVERVIEW'",
        selectInput(
          inputId = "TSA_NETWORK_TYPE_CHOICE",
          label = "Network Type",
          choices = c("Young", "Old", "Change", "ORA")
        ),
        selectInput(
          inputId = "TSA_NETWORK_LAYOUT_CHOICE",
          label = "Layout Type",
          choices = c("Circular", "Bipartite")
        )
      ),
      conditionalPanel(
        condition = "input.active_TSA_panel == 'TSA_INTERACTION_ANALYSIS'",
        selectInput(
          inputId = "TSA_CCI_DETAILS_CHOICE",
          label = "Type of Information",
          choices = c("CCI Table", "Volcano Plot", "Score Plot")
        ),
        hr(),
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
        selectInput(
          inputId = "TSA_ORA_DETAILS_CHOICE",
          label = "Type of Information",
          choices = c("ORA Table", "ORA Plot")
        ),
        hr(),
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
      width = 3#,
      #style = "position: fixed; width: 20%"
    ),
    mainPanel(
      titlePanel(htmlOutput("TSA_TITLE")),
      tabsetPanel(
        type = "tabs",
        tabPanel(
          title = "Overview",
          fluidRow(
            column(width = 12, htmlOutput("TSA_OVERVIEW_INTRO"), style = "padding:50px"),
            column(width = 12, DT::dataTableOutput("TSA_OVERVIEW_TABLE"), style = "padding-bottom: 50px"),
            column(width = 12, htmlOutput("TSA_NETWORK_INTRO"), style = "padding:50px"),
            column(width = 12, visNetworkOutput("TSA_NETWORK_PLOT", height = "800px"))
          ),
          value = "TSA_OVERVIEW"
        ),
        tabPanel(
          title = "Detailed Interactions",
          fluidRow(
            column(width = 12, htmlOutput("TSA_CCI_INTRO"), style = "padding:50px"),
            column(width = 12, uiOutput("TSA_CCI_DETAILS")),
            column(width = 12, uiOutput("TSA_CCI_TEXTOUTPUT"), style = "padding:10px")
          ),
          value = "TSA_INTERACTION_ANALYSIS"
        ),
        tabPanel(
          title = "Over-representation Analysis",
          fluidRow(
            column(width = 12, htmlOutput("TSA_ORA_INTRO"), style = "padding:50px"),
            column(width = 12, uiOutput("TSA_ORA_DETAILS"))
          ),
          value = "TSA_ORA"
        ),
        id = "active_TSA_panel"
      )
    )
  )
)
