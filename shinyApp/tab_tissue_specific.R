tab_tissue_specific <- tabPanel(
  title = "Tissue-Specific Analysis",
  #titlePanel(htmlOutput("TSA_TITLE")),
  fluidRow(
    column(width = 8, titlePanel(htmlOutput("TSA_TITLE")), offset = 2),
    column(width = 6, DT::dataTableOutput("TSA_OVERVIEW_TABLE"), style = "padding-bottom: 50px", offset = 3)
  ),
  tabsetPanel(
    type = "tabs",
    tabPanel(
      title = "Cell-Cell Interactions",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput(
            inputId = "TSA_CCI_DETAILS_CHOICE",
            label = "Choose between Table or Plots",
            choices = c("CCI Table", "Volcano Plot", "Score Plot", "LRI-FC Plot")
          ),
          uiOutput("TSA_EMITTER_CHOICE"),
          uiOutput("TSA_RECEIVER_CHOICE"),
          uiOutput("TSA_LRI_CHOICE"),
          sliderInput(
            inputId = "TSA_SLIDER_PVALUE",
            label = "Filter by Adj. p-value",
            min = 0,
            max = 1,
            value = 1
          ),
          uiOutput("TSA_SLIDER_LOG2FC")
        ),
        mainPanel(
          fluidRow(
            # column(width = 3,
            #        selectInput(
            #          inputId = "TSA_CCI_DETAILS_CHOICE",
            #          label = "Type of Information",
            #          choices = c("CCI Table", "Volcano Plot", "Score Plot")
            #        )
            # ),
            # column(width = 3, uiOutput("TSA_EMITTER_CHOICE")),
            # column(width = 3, uiOutput("TSA_RECEIVER_CHOICE")),
            # column(width = 3, uiOutput("TSA_LRI_CHOICE")),
            column(width = 12, uiOutput("TSA_CCI_DETAILS"), style = "padding:50px"),
            column(width = 12, uiOutput("TSA_CCI_TEXTOUTPUT"), style = "padding:50px"),
            column(width = 12, htmlOutput("TSA_CCI_INTRO"), style = "padding:50px")
          )
        )
      ),
      value = "TSA_INTERACTION_ANALYSIS"
    ),
    tabPanel(
      title = "Network Representation",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput(
            inputId = "TSA_NETWORK_TYPE_CHOICE",
            label = "Choose the Type of Network",
            choices = c("Young", "Old", "Change", "ORA")
          ),
          selectInput(
            inputId = "TSA_NETWORK_LAYOUT_CHOICE",
            label = "Choose between Conventional or Bipartite Layout",
            choices = c("Circular", "Bipartite")
          )
        ),
        mainPanel(
          fluidRow(
            column(width = 12, visNetworkOutput("TSA_NETWORK_PLOT", height = "800px")),
            column(width = 12, htmlOutput("TSA_NETWORK_INTRO"), style = "padding:50px"),
            column(width = 12, htmlOutput("TSA_OVERVIEW_INTRO"), style = "padding:50px")
          )
          
        )
      ),
      value = "TSA_OVERVIEW"
    ),
    tabPanel(
      title = "Over-Representation Analysis",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput(
            inputId = "TSA_ORA_DETAILS_CHOICE",
            label = "Choose between Table or Plot",
            choices = c("ORA Table", "ORA Plot")
          ),
          uiOutput("TSA_ORA_CATEGORY_CHOICE"),
          selectInput(
            inputId = "TSA_ORA_TYPE_CHOICE",
            label = "Age Regulation",
            choices = list("Up", "Down", "Flat")
          )
        ),
        mainPanel(
          fluidRow(
            column(width = 12, uiOutput("TSA_ORA_DETAILS"), style = "padding:50px"),
            column(width = 12, htmlOutput("TSA_ORA_INTRO"), style = "padding:50px")
          )
        )
      ),
      value = "TSA_ORA"
    ),
    id = "active_TSA_panel"
  )
  #,
  # sidebarLayout(
  #   sidebarPanel(
  #     #uiOutput("TSA_DATASET_CHOICE"),
  #     #uiOutput("TSA_TISSUE_CHOICE"),
  #     hr(),
  #     conditionalPanel(
  #       condition = "input.active_TSA_panel == 'TSA_OVERVIEW'",
  #       selectInput(
  #         inputId = "TSA_NETWORK_TYPE_CHOICE",
  #         label = "Network Type",
  #         choices = c("Young", "Old", "Change", "ORA")
  #       ),
  #       selectInput(
  #         inputId = "TSA_NETWORK_LAYOUT_CHOICE",
  #         label = "Layout Type",
  #         choices = c("Circular", "Bipartite")
  #       )
  #     ),
  #     conditionalPanel(
  #       condition = "input.active_TSA_panel == 'TSA_INTERACTION_ANALYSIS'",
  #       # selectInput(
  #       #   inputId = "TSA_CCI_DETAILS_CHOICE",
  #       #   label = "Type of Information",
  #       #   choices = c("CCI Table", "Volcano Plot", "Score Plot")
  #       # ),
  #       hr(),
  #       # uiOutput("TSA_EMITTER_CHOICE"),
  #       # uiOutput("TSA_RECEIVER_CHOICE"),
  #       # uiOutput("TSA_LRI_CHOICE"),
  #       sliderInput(
  #         inputId = "TSA_SLIDER_PVALUE",
  #         label = "P-value Threshold",
  #         min = 0, 
  #         max = 1,
  #         value = 1
  #       ),
  #       uiOutput("TSA_SLIDER_LOG2FC")
  #     ),
  #     conditionalPanel(
  #       condition = "input.active_TSA_panel == 'TSA_ORA'",
  #       selectInput(
  #         inputId = "TSA_ORA_DETAILS_CHOICE",
  #         label = "Type of Information",
  #         choices = c("ORA Table", "ORA Plot")
  #       ),
  #       hr(),
  #       uiOutput("TSA_ORA_CATEGORY_CHOICE"),
  #       selectInput(
  #         inputId = "TSA_ORA_TYPE_CHOICE",
  #         label = "Regulation",
  #         choices = list("Up", "Down", "Stable")
  #       ),
  #       sliderInput(
  #         inputId = "TSA_ORA_SLIDER_PVALUE",
  #         label = "P-value Threshold",
  #         min = 0, 
  #         max = 1,
  #         value = 1
  #       ),
  #       uiOutput("TSA_ORA_SLIDER_OR")
  #     ),
  #     width = 3#,
  #     #style = "position: fixed; width: 20%"
  #   )
  # )
)
