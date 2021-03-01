server_scAgeCom <- function(
  input,
  output
) {
  source("utils_description.R")
  output$description_html <- get_description_html(input)
  source("utils_tissue_specific.R", local = TRUE)
  output$TSA_TISSUE_CHOICE <- choose_TSA_tissue(input)
  output$TSA_TITLE <- get_TSA_title(input)
  # tab overview
  output$TSA_OVERVIEW_INTRO <- get_TSA_overview_intro(input)
  output$TSA_OVERVIEW_TABLE <- get_TSA_overview_table(input)
  output$TSA_NETWORK_INTRO <- get_TSA_network_intro(input)
  output$TSA_NETWORK_PLOT <- plot_TSA_network(input)
  # tab detailed interactions
  output$TSA_CCI_INTRO <- get_TSA_cci_intro(input)
  output$TSA_EMITTER_CHOICE <- choose_TSA_emitter(input)
  output$TSA_RECEIVER_CHOICE <- choose_TSA_receiver(input)
  output$TSA_SLIDER_LOG2FC <- get_TSA_slider_log2fc(input)
  output$TSA_CCI_DETAILS <- get_TSA_cci_details(input)
  output$TSA_CCI_TEXTOUTPUT <- get_TSA_cci_text(input)
  output$TSA_INTERACTION_TABLE <- get_TSA_interaction_table(input)
  output$TSA_VOLCANO_PLOT <- plot_TSA_VOLCANO(input)
  output$TSA_VOLCANO_TEXTOUTPUT <- get_TSA_VOLCANO_text(input)
  output$TSA_SCORES_PLOT <- plot_TSA_SCORES(input)
  output$TSA_SCORES_TEXTOUTPUT <- get_TSA_SCORES_text(input)
  #tab over-representation
  output$TSA_ORA_INTRO <- get_TSA_ora_intro(input)
  output$TSA_ORA_CATEGORY_CHOICE <- choose_TSA_ORA_category(input)
  output$TSA_ORA_SLIDER_OR <- get_TSA_ORA_slider_or(input)
  output$TSA_ORA_DETAILS <- get_TSA_ora_details(input)
  output$TSA_ORA_TABLE <- get_TSA_ORA_table(input)
  output$TSA_ORA_PLOT <- plot_TSA_ORA(input)
  source("utils_combined_analysis.R", local = TRUE)
  output$TCA_TITLE <- get_TCA_title(input)
  output$TCA_SUMMARY_TYPE <- get_TCA_summary_type(input)
  output$TCA_SUMMARY_INTRO <- get_TCA_summary_intro(input)
  output$TCA_SUMMARY_TABLE <- get_TCA_summary_table(input)
  output$TCA_ORA_PLOT <- plot_TCA_ora(input)
  #
  output$TCA_CCI_INTRO <- get_TCA_cci_intro(input)
  output$TCA_CCI_DETAILS <- get_TCA_cci_details(input)
  output$TCA_TISSUE_CHOICE <- choose_TCA_tissue(input)
  output$TCA_REGULATION_CHOICE <- choose_TCA_regulation(input)
  output$TCA_CCI_TEXTOUTPUT <- get_TCA_cci_text(input)
  output$TCA_INTERACTION_TABLE <- get_TCA_interaction_table(input)
  output$TCA_VOLCANO_PLOT <- plot_TCA_VOLCANO(input)
  output$TCA_VOLCANO_TEXTOUTPUT <- get_TCA_VOLCANO_text(input)
  output$TCA_SCORES_PLOT <- plot_TCA_SCORES(input)
  output$TCA_SCORES_TEXTOUTPUT <- get_TCA_SCORES_text(input)

    source("utils_LRdb.R", local = TRUE)
  output$LRdb_TABLE <- show_LRdb_table(input)
  output$LRdb_UPSET_PLOT <- show_LRdb_upset(input)
}

