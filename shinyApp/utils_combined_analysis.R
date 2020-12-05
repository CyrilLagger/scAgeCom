get_TCA_table <- function(
  input
) {
  DT::renderDT({
    req(input$TCA_DATASET_CHOICE)
    dt <- LR_GENES_summary[[input$TCA_DATASET_CHOICE]]
    req(dt)
    #setorder(
    #  dt,
    #  -LOG2FC,
    #  `Adj. P-Value`
    #)
    show_DT(
      dt,
      colnames(dt)[c(1,2,3,4,5,6,7,11,14,15)],
      #cols_numeric_DATA,
      colnames(dt)[c(2,3,4,5,6,7,11,14)],
      "Summary Table"
    )
  })
}