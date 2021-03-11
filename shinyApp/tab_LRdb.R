# tab_LRdb <- tabPanel(
#   title = "Ligand-Receptor Databases",
#   sidebarLayout(
#     sidebarPanel(
#       conditionalPanel(
#         condition = "input.active_LRdb_panel=='LRdb_TABLE'",
#         pickerInput(
#           inputId = "LRdb_DATABASE",
#           label = "Databases",
#           choices = c("CellChat", "CellPhoneDB", "CellTalkDB", "connectomeDB2020",
#                       "ICELLNET", "NicheNet", "SingleCellSignalR", "scTensor"),
#           selected = c("CellChat", "CellPhoneDB", "CellTalkDB", "connectomeDB2020",
#                        "ICELLNET", "NicheNet", "SingleCellSignalR", "scTensor"),
#           options = list(`actions-box` = TRUE),
#           multiple = TRUE
#         )
#       ),
#       conditionalPanel(
#         condition = "input.active_LRdb_panel=='LRdb_UPSET_PLOT'",
#         selectInput(
#           inputId = "LRdb_PLOT_CHOICE",
#           label = "Graph Type",
#           choices = c("By Databases", "By Sources")
#         )
#       ),
#       width = 2
#     ),
#     mainPanel(
#       tabsetPanel(
#         type = "tabs",
#         tabPanel(
#           title = "Table",
#           DT::dataTableOutput("LRdb_TABLE"),
#           value = "LRdb_TABLE"
#           ),
#         tabPanel(
#           title = "Graph Summary",
#           plotOutput("LRdb_UPSET_PLOT", height = "600px"),
#           value = "LRdb_UPSET_PLOT"
#           ),
#         id = "active_LRdb_panel"
#       )
#     )
#   )
# )