tab_combined_analysis <- tabPanel(
  title = "Combined Analysis",
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "TCA_DATASET_CHOICE",
        label = "Dataset",
        choices = c("TMS FACS Data", "TMS Droplet Data")
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