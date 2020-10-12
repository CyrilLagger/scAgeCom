

# need to be done properly. 

diffcom_results <- readRDS("../../data_scAgeCom/analysis/analysis_4_data_diffcom_filter_new.rds")

LRkeep <- LRall[scsr == TRUE | cpdb == TRUE]

table(diffcom_results$tms_facs$CASE_TYPE)

LR_det_facs <- unique(diffcom_results$tms_facs$LR_GENES)

LR_up_facs <- unique(diffcom_results$tms_facs[CASE_TYPE %in% c("FTTU", "TTTU")]$LR_GENES)
LR_down_facs <- unique(diffcom_results$tms_facs[CASE_TYPE %in% c("TFTD", "TTTD")]$LR_GENES)
intersect(LR_up_facs, LR_down_facs)

LR_strong_up_facs <- unique(diffcom_results$tms_facs[CASE_TYPE %in% c("FTTU", "TTTU") & LR_LOGFC >= log(2)]$LR_GENES)
LR_strong_down_facs <- unique(diffcom_results$tms_facs[CASE_TYPE %in% c("TFTD", "TTTD") & LR_LOGFC <= -log(2)]$LR_GENES)
intersect(LR_strong_up_facs, LR_strong_down_facs)

LR_up_droplet <- unique(diffcom_results$tms_droplet[CASE_TYPE %in% c("FTTU", "TTTU")]$LR_GENES)
LR_down_droplet <- unique(diffcom_results$tms_droplet[CASE_TYPE %in% c("TFTD", "TTTD")]$LR_GENES)
intersect(LR_up_droplet, LR_down_droplet)

intersect(LR_up_droplet, LR_up_facs)
intersect(LR_down_droplet, LR_down_facs)

LRkeep[GENESYMB_L == "Igf1"]

LRkeep[GENESYMB_L == "Kl"]

diffcom_results$tms_facs[L_GENE]