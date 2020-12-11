tab_LRdb <- tabPanel(
  title = "Ligand-Receptor Database",
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        condition = "input.active_LRdb_panel=='LRdb_TABLE'",
        pickerInput(
          inputId = "LRdb_DATABASE",
          label = "Database",
          choices = c("CONNECTOMEDB", "CELLCHAT", "CELLPHONEDB", "CELLTALK", "ICELLNET", "NICHENET", "SCSR", "SCTENSOR"),
          selected = c("CONNECTOMEDB", "CELLCHAT", "CELLPHONEDB", "CELLTALK", "ICELLNET", "NICHENET", "SCSR", "SCTENSOR"),
          options = list(`actions-box` = TRUE),
          multiple = TRUE
        )
      ),
      conditionalPanel(
        condition = "input.active_LRdb_panel=='LRdb_UPSET_PLOT'",
        selectInput(
          inputId = "LRdb_PLOT_CHOICE",
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
          DT::dataTableOutput("LRdb_TABLE"),
          value = "LRdb_TABLE"
          ),
        tabPanel(
          title = "Upset Plot",
          plotOutput("LRdb_UPSET_PLOT"),
          value = "LRdb_UPSET_PLOT"
          ),
        id = "active_LRdb_panel"
      )
    )
  )
)