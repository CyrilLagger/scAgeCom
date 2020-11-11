tabPanel(
  title = "Global Analysis",
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "dataset_selection",
        label = "Dataset",
        choices = c("TMS FACS Data", "TMS Droplet Data", "Calico Data")
      ),
      uiOutput("tissue_selection"),
      uiOutput("ligand_cell_selection"),
      uiOutput("receptor_cell_selection"),
      sliderInput(
        inputId = "slider_pvalue",
        label = "P-value filter",
        min = 0, 
        max = 1,
        value = 1
      ),
      uiOutput("log2fc_slider"),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        type = "tabs",
        tabPanel("Table", DT::dataTableOutput("scdiffcom_table")),
        tabPanel("Volcano Plot", plotOutput("scdiffcom_volcano"))
      )
    )
  ))