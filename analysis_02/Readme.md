# Analysis_02


- 01: Random gene pairs analysis.
  - 01a_compute_regulation_stats.R
  Shuffle seurat genes. Results in "data/shuffled_regulation".
  - 01b_analyze_regulation_stats.R
  Statistical analysis and plots in code (# detections in random and increased detection when using LRI).
- 02: Coexpression analysis.
  - 02_compute_coexpression.R
  Computationally expensive processing of co-expression and plots at "data/coexpression_plots".
  - 02b_compute_scdiffcom_with_coexpression.R
  Precomputation for 02c. Results at "data/datasets_detections".
  - 02c_analyze_scdiffcom_coexpression.R
  Compare scDiffCom results with coexpression. Plots in code. No clear pattern, not used in the paper.
- 03: Connections between DGE and LRI.
  - 03a_compute_expression.R
  Precomputation for 03b. Results at "data/expression".
  - 03b_analyze_expression.R
  Plots in code. LRI have a tendency to be in DGE across datasets.
- 04: Compare our results with those from other omics studies.
  - 04a_proteins.R: compare with a lung proteomics study
  - 04b_lung_scrnaseq.R
  - 04c_lung_scdiffcom.R
- 05_GO_Terms_levels.R
  Short program to check the median GO level across datasets.
- 06_sex_analysis.R
  Analyse differences across sex and age.
- 07_complex_LRI.R


Support files
- Readme.md
- coexpression.R
- random_gene_pairs_analysis.R
- utils_random_pairs.R