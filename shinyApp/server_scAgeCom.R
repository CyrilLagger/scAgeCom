server_scAgeCom <- function(
  input,
  output
) {
  source("utils_description.R")
  output$description_html <- get_description_html(input)
  source("utils_tissue_specific.R", local = TRUE)
  output$TSA_TITLE <- get_TSA_title(input)
  output$TSA_OVERVIEW <- get_TSA_overview(input)
  output$TSA_NETWORK_PLOT <- plot_TSA_network(input)
  output$TSA_TISSUE_CHOICE <- choose_TSA_tissue(input)
  output$TSA_EMITTER_CHOICE <- choose_TSA_emitter(input)
  output$TSA_RECEIVER_CHOICE <- choose_TSA_receiver(input)
  output$TSA_SLIDER_LOG2FC <- get_TSA_slider_log2fc(input)
  output$TSA_ORA_CATEGORY_CHOICE <- choose_TSA_ORA_category(input)
  output$TSA_ORA_SLIDER_OR <- get_TSA_ORA_slider_or(input)
  output$TSA_INTERACTION_TABLE <- get_TSA_interaction_table(input)
  output$TSA_VOLCANO_PLOT <- plot_TSA_VOLCANO(input)
  output$TSA_VOLCANO_TEXTOUTPUT <- get_TSA_VOLCANO_text(input)
  output$TSA_SCORES_PLOT <- plot_TSA_SCORES(input)
  output$TSA_SCORES_TEXTOUTPUT <- get_TSA_SCORES_text(input)
  output$TSA_ORA_TABLE <- get_TSA_ORA_table(input)
  output$TSA_ORA_PLOT <- plot_TSA_ORA(input)
  source("utils_LRdb.R", local = TRUE)
  output$LRdb_TABLE <- show_LRdb_table(input)
  output$LRdb_UPSET_PLOT <- show_LRdb_upset(input)
  source("utils_combined_analysis.R", local = TRUE)
  output$TCA_TABLE <- get_TCA_table(input)
}

