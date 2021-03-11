# show_LRdb_table <- function(
#   input
# ) 
# {
#   DT::renderDataTable({
#     req(input$LRdb_DATABASE)
#     dt <- LRdb_curated[
#       apply(
#         sapply(
#           input$LRdb_DATABASE,
#           function(i) {
#             grepl(i, LRdb_curated$DATABASE)
#           }
#         ),
#         MARGIN = 1,
#         any
#       )
#       ]
#     cols_to_show_LRdb <- c(
#       "LIGAND_1", "LIGAND_2",
#       "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3",
#       "DATABASE", "SOURCE"
#     )
#     dt <- dt[, cols_to_show_LRdb, with = FALSE]
#     setcolorder(dt, cols_to_show_LRdb)
#     options_LRdb <- list(
#       pageLength = 10,
#       columnDefs = list(
#         list(
#           targets = c(6,7),
#           render = htmlwidgets::JS(
#             "function(data, type, row, meta) {",
#             "return type === 'display' && data.length > 20 ?",
#             "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
#             "}")
#         )
#       )
#     )
#     show_DT(
#       data = dt,
#       cols_to_show = cols_to_show_LRdb,
#       cols_numeric = NULL,
#       table_title = "Table of Ligand-Receptor Interactions",
#       options = options_LRdb
#     )
#   })
# }
# 
# show_LRdb_upset <- function(
#   input
# ) {
#   renderPlot({
#     req(input$LRdb_PLOT_CHOICE)
#     LRdb_sources <- c(
#       "FANTOM5", "HPMR", "HPRD", "PMID",
#       "CellPhoneDB", "KEGG", "IUPHAR", "reactome",
#       "cellsignal.com", "PPI"
#     )
#     LRdb_DBS <- c(
#       "CellChat", "CellPhoneDB", "CellTalkDB", "connectomeDB2020",
#       "ICELLNET", "NicheNet", "SingleCellSignalR", "scTensor"
#     )
#     dt <- LRdb_curated
#     dt[, COMPLEX := !is.na(LIGAND_2) | !is.na(RECEPTOR_2)]
#     dt[, c(LRdb_DBS) := lapply(LRdb_DBS, function(i) {
#       ifelse(grepl(i, DATABASE), TRUE, FALSE)
#     })]
#     dt[, c(LRdb_sources) := lapply(LRdb_sources, function(i) {
#       ifelse(grepl(i, SOURCE), TRUE, FALSE)
#     })]
#     if(input$LRdb_PLOT_CHOICE == "By Databases") {
#       ComplexUpset::upset(
#         dt,
#         LRdb_DBS,
#         base_annotations=list(
#           'Intersection size'=intersection_size(
#             counts=TRUE,
#             aes=aes(fill=COMPLEX),
#             text = list(position = position_stack(vjust = 0.0)),
#             bar_number_threshold = 100
#           )
#         ),
#         themes=upset_default_themes(text=element_text(size=20)),
#         min_size = 39
#       ) +
#         ggtitle("Overlap by Databases of Origin")
#     } else if(input$LRdb_PLOT_CHOICE == "By Sources") {
#       ComplexUpset::upset(
#         dt,
#         LRdb_sources,
#         base_annotations=list(
#           'Intersection size'=intersection_size(
#             counts=TRUE,
#             aes=aes(fill=COMPLEX),
#             text = list(position = position_stack(vjust = 0.0)),
#             bar_number_threshold = 100
#           )
#         ),
#         themes=upset_default_themes(text=element_text(size=20)),
#         min_size = 30
#       ) +
#         ggtitle("Overlap by Sources of Origin")
#     }
#   })
# }