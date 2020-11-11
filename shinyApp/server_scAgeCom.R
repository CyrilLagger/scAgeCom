server_scAgeCom <- function(
  input,
  output
) {
  source("utils_LR6db.R", local = TRUE)
  output$LR6db_TABLE <- show_LR6db_table(input)
  output$LR6db_UPSET_PLOT <- show_LR6db_upset(input)
  source("utils_tissue_specific.R", local = TRUE)
  output$TSA_TISSUE_CHOICE <- choose_TSA_tissue(input)
  output$TSA_EMITTER_CHOICE <- choose_TSA_emitter(input)
  output$TSA_RECEIVER_CHOICE <- choose_TSA_receiver(input)
  output$TSA_SLIDER_LOG2FC <- get_TSA_slider_log2fc(input)
  output$TSA_ORA_CATEGORY_CHOICE <- choose_TSA_ORA_category(input)
  output$TSA_ORA_SLIDER_OR <- get_TSA_ORA_slider_or(input)
  output$TSA_INTERACTION_TABLE <- get_TSA_interaction_table(input)
  output$TSA_VOLCANO_PLOT <- plot_TSA_VOLCANO(input)
  output$TSA_ORA_TABLE <- get_TSA_ORA_table(input)
  output$TSA_ORA_PLOT <- plot_TSA_ORA(input)
}

