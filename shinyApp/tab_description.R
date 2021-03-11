tab_description <- tabPanel(
  title = "Introduction",
  fluidRow(
    column(width = 6, titlePanel(htmlOutput("INTRO_TITLE")), offset = 3),
  ),
  tabsetPanel(
    type = "tabs",
    tabPanel(
      title = "Overview",
      fluidRow(
        column(width = 12,htmlOutput("INTRO_OVERVIEW"), style = "padding:50px")
      ),
      value = "INTRO_OVERVIEW"
    ),
    tabPanel(
      title = "Methodology",
      value = "INTRO_METHOD"
    ),
    tabPanel(
      title = "Single-cell Data",
      fluidRow(
        column(width = 12, htmlOutput("INTRO_SCRNA_HTML"), style = "padding:50px")
      ),
      value = "INTRO_SCRNA_DATA"
    ),
    tabPanel(
      title = "Ligand-Receptor Database",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput(
            inputId = "INTRO_LRI_DETAILS_CHOICE",
            label = "Choose between Table or Plot",
            choices = c("Ligand-Receptor Table", "Upset Plot 1", "Upset Plot 2")
          ),
          conditionalPanel(
            condition = "input.INTRO_LRI_DETAILS_CHOICE == 'Ligand-Receptor Table'",
            pickerInput(
              inputId = "LRdb_DATABASE",
              label = "Databases",
              choices = c("CellChat", "CellPhoneDB", "CellTalkDB", "connectomeDB2020",
                          "ICELLNET", "NicheNet", "SingleCellSignalR", "scTensor"),
              selected = c("CellChat", "CellPhoneDB", "CellTalkDB", "connectomeDB2020",
                           "ICELLNET", "NicheNet", "SingleCellSignalR", "scTensor"),
              options = list(`actions-box` = TRUE),
              multiple = TRUE
            )
          )
          # conditionalPanel(
          #   condition = "input.active_LRdb_panel=='LRdb_TABLE'",
          # ),
          # conditionalPanel(
          #   condition = "input.active_LRdb_panel=='LRdb_UPSET_PLOT'",
          #   selectInput(
          #     inputId = "LRdb_PLOT_CHOICE",
          #     label = "Graph Type",
          #     choices = c("By Databases", "By Sources")
          #   )
          # )
        ),
        mainPanel(
          fluidRow(
            column(width = 12, htmlOutput("INTRO_LRI_HTML"), style = "padding:50px"),
            column(width = 12, uiOutput("INTRO_LRI_DETAILS"), style = "padding:50px")
          )
        )
      ),
      value = "INTRO_LRI_DATA"
    ),
    id = "active_INTRO_panel"
  )
)