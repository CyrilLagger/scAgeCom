tab_LR6db <- tabPanel(
  title = "Ligand-Receptor Database",
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        condition = "input.active_LR6db_panel=='LR6db_TABLE'",
        pickerInput(
          inputId = "LR6db_DATABASE",
          label = "Database",
          choices = c("CELLCHAT", "CELLPHONEDB", "ICELLNET", "NICHENET", "SCSR", "SCTENSOR"),
          selected = c("CELLCHAT", "CELLPHONEDB", "ICELLNET", "NICHENET", "SCSR", "SCTENSOR"),
          options = list(`actions-box` = TRUE),
          multiple = TRUE
        )
      ),
      conditionalPanel(
        condition = "input.active_LR6db_panel=='LR6db_UPSET_PLOT'",
        selectInput(
          inputId = "LR6db_PLOT_CHOICE",
          label = "Plot type",
          choices = c("By Database", "By Source")
        )
      ),
      width = 2
    ),
    mainPanel(
      tabsetPanel(
        type = "tabs",
        tabPanel(
          title = "Table",
          DT::dataTableOutput("LR6db_TABLE"),
          value = "LR6db_TABLE"
          ),
        tabPanel(
          title = "Upset Plot",
          plotOutput("LR6db_UPSET_PLOT"),
          value = "LR6db_UPSET_PLOT"
          ),
        id = "active_LR6db_panel"
      )
    )
  )
)