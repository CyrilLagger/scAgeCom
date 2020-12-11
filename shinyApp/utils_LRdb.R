show_LRdb_table <- function(
  input
) 
{
  DT::renderDataTable({
    req(input$LRdb_DATABASE)
    dt <- LRdb_curated[
      apply(
        sapply(
          input$LRdb_DATABASE,
          function(i) {
            grepl(i, LRdb_curated$DATABASE)
          }
        ),
        MARGIN = 1,
        any
      )
      ]
    show_DT(
      data = setorder(
        dt,
        LIGAND_1, LIGAND_2, RECEPTOR_1, RECEPTOR_2, RECEPTOR_3
      ),
      cols_to_show = cols_to_show_LRdb,
      NULL
    )
  })
}

show_LRdb_upset <- function(
  input
) {
  NULL
  # renderPlot({
  #   req(input$LRdb_PLOT_CHOICE)
  #   category_sources <- c("PPI", "IUPHAR", "KEGG", "PMID", "CPDB", "Ramilowski",
  #                         "cellsignal.com", "HPMR", "HPRD", "reactome")
  #   category_DBs <- c("CELLCHAT", "CELLPHONEDB", "ICELLNET", "NICHENET", "SCSR", "SCTENSOR")
  #   dt <- LRdb_curated
  #   dt[, c(category_sources) := lapply(
  #     category_sources,
  #     function(i) {
  #       ifelse(grepl(i, SOURCE), 1, 0)
  #     }
  #   )]
  #   dt[, c(category_DBs) := lapply(
  #     category_DBs,
  #     function(i) {
  #       ifelse(grepl(i, DATABASE), 1, 0)
  #     }
  #   )]
  #   dt[, complex := !is.na(LIGAND_2) | !is.na(RECEPTOR_2)]
  #   if(input$LRdb_PLOT_CHOICE == "By Database") {
  #     ComplexUpset::upset(
  #       dt,
  #       category_DBs,
  #       base_annotations=list(
  #         'Intersection size'=intersection_size(
  #           counts=TRUE,
  #           aes=aes(fill=complex),
  #           text = list(position = position_stack(vjust = 0.0)),
  #           bar_number_threshold = 100
  #         )
  #       ),
  #       min_size = 20
  #     )
  #   } else if(input$LRdb_PLOT_CHOICE == "By Source") {
  #     ComplexUpset::upset(
  #       dt,
  #       category_sources,
  #       base_annotations=list(
  #         'Intersection size'=intersection_size(
  #           counts=TRUE,
  #           aes=aes(fill=complex),
  #           text = list(position = position_stack(vjust = 0.0)),
  #           bar_number_threshold = 100
  #         )
  #       ),
  #       min_size = 30
  #     )
  #   }
  # })
}