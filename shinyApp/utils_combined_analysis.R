get_TCA_table <- function(
  input
) {
  DT::renderDT({
    req(input$TCA_DATASET_CHOICE, input$TCA_CATEGORY_CHOICE)
    dt <-  CCI_SUMMARY[[input$TCA_DATASET_CHOICE]][[input$TCA_CATEGORY_CHOICE]]
    show_DT(
      dt,
      colnames(dt),
      NULL
    )
  })
}