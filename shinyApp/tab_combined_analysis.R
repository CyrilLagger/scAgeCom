tab_combined_analysis <- tabPanel(
  title = "Combined Analysis",
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "TCA_DATASET_CHOICE",
        label = "Dataset",
        choices =  c("TMS FACS Data", "TMS Droplet Data", "Calico Data")
      ),
      selectInput(
        inputId = "TCA_CATEGORY_CHOICE",
        label = "Category",
        choices =  c("ID", "ER_CELL_FAMILY_GLOBAL", "LR_GENES_GLOBAL", "GO_TERMS_GLOBAL")
      ),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        type = "tabs",
        tabPanel("Table", DT::dataTableOutput("TCA_TABLE"))#,
        #tabPanel("Volcano Plot", plotOutput("scdiffcom_volcano"))
      )
    )
  ))