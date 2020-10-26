float_equal = function(a, b, eps = 10**(-5)) {
    return( abs(a-b) < eps )
}

test_src_1_filtering = function() {
    
    DATASET_ID = "TMS Droplet Female"
    # Setup expected values
    NUM_COLUMNS = 52
    PVAL_DIFF_MEAN = 0.2952629
    BHPVAL_DIFF_MEAN = 0.4493505
    LR_SCORE_COND1_MEAN = 1.068921
    LR_SCORE_COND2_MEAN = 1.040462
    LR_DETECT_AND_SIGN_COND1_SUM = 54889
    LR_DETECT_AND_SIGN_COND2_SUM = 53515
    DIFFERENTIAL_EXPRESSED_SUM = 14387
    DIFFERENTIAL_DIRECTION_COUNT_UP = 31037
    CASE_TYPE_COUNT_TTTU = 3030
    CASE_TYPE_COUNT_FFF = 19151
    CASE_TYPE_COUNT_FTTU = 2865
    
    
    # Setup test conditions
    CCI_all <- readRDS("../data_scAgeCom/analysis/a4_data_diffcom_all_filtered.rds")
    data = CCI_all[[DATASET_ID]]$results
    conditions = CCI_all[[DATASET_ID]]$conds
    cutoff_quantile = 0.25
    cutoff_logFC_abs = log(1.1)
    
    # Define colnames
    # LR_SCORE_COND1_COLNAME = paste0("LR_SCORE_", conditions[[1]])
    # LR_SCORE_COND2_COLNAME = paste0("LR_SCORE_", conditions[[2]])
    # LR_DETECT_AND_SIGN_COND1_COLNAME = paste0(
    #     "LR_DETECTED_AND_SIGNIFICANT_IN_", conditions[[1]]
    # )
    # LR_DETECT_AND_SIGN_COND2_COLNAME = paste0(
    #     "LR_DETECTED_AND_SIGNIFICANT_IN_", conditions[[2]]
    # )
    # DIFFERENTIAL_EXPRESSED_COLNAME = "DIFFERENTIAL_EXPRESSED"
    # DIFFERENTIAL_DIRECTION_COLNAME = "DIFFERENTIAL_DIRECTION"
    # CASE_TYPE_COLNAME = "CASE_TYPE"
    
    # Compute and compare results to expected
    res = analyze_CCI_per_tissue(
        data, conditions, 
        cutoff_quantile, cutoff_logFC_abs
    )
    
    TESTS = list(
        NUM_COLS = dim(res)[2] == NUM_COLUMNS,
        PVAL_DIFF = float_equal(
            mean(res$PVAL_DIFF), PVAL_DIFF_MEAN
        ),
        BH_PVAL_DIFF = float_equal(
            mean(res$BH_PVAL_DIFF), BHPVAL_DIFF_MEAN
        ),
        LR_SCORE_COND1 = float_equal(
            mean(res$LR_SCORE_YOUNG), LR_SCORE_COND1_MEAN
        ),
        LR_SCORE_COND2 = float_equal(
            mean(res$LR_SCORE_OLD), LR_SCORE_COND2_MEAN
        ),
        LR_DETECT_AND_SIGN_COND1 = (sum(res$LR_DETECTED_AND_SIGNIFICANT_IN_YOUNG) 
                                    == LR_DETECT_AND_SIGN_COND1_SUM),
        LR_DETECT_AND_SIGN_COND2 = (sum(res$LR_DETECTED_AND_SIGNIFICANT_IN_OLD)
                                    == LR_DETECT_AND_SIGN_COND2_SUM),
        DIFFERENTIAL_EXPRESSED = (sum(res$DIFFERENTIAL_EXPRESSED) 
                                  == DIFFERENTIAL_EXPRESSED_SUM),
        DIFFERENTIAL_DIRECTION = (sum(res$DIFFERENTIAL_DIRECTION == 'UP')
                                  == DIFFERENTIAL_DIRECTION_COUNT_UP),
        CASE_TYPE_TTTU = sum(res$CASE_TYPE == 'TTTU') == CASE_TYPE_COUNT_TTTU,
        CASE_TYPE_FFF = sum(res$CASE_TYPE == 'FFF') == CASE_TYPE_COUNT_FFF,
        CASE_TYPE_FTTU = sum(res$CASE_TYPE == 'FTTU') == CASE_TYPE_COUNT_FTTU
    )
    
    if( all(as.logical(TESTS)) ) {
        message("test_src_1_filtering: PASSED")
    } else {
        print( TESTS )
        stop("test_src_1_filtering: FAILED")
    }
}
