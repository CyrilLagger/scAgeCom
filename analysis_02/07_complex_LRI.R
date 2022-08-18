library(glue)
library(data.table)
library(scDiffCom)


## 1. Ligands and receptors in complex LRIs
lri = scDiffCom::LRI_mouse$LRI_curated
lri_complex_L = lri[!is.na(LIGAND_2)]
lri_complex_R = lri[!is.na(RECEPTOR_2)]
# lri[!is.na(RECEPTOR_3)] - very few

print(glue("# complex ligands LRIs: {nrow(lri_complex_L)}"))  # Integrins, some IL, Inhba
print(glue("# complex receptors LRIs: {nrow(lri_complex_R)}"))

lri_complex_R_uniq = unique(lri_complex_R[, list(RECEPTOR_1, RECEPTOR_2, RECEPTOR_3)])
print(glue("# multi-subunit receptors: {nrow(lri_complex_R_uniq)}"))
R_units = c(lri_complex_R_uniq$RECEPTOR_1, lri_complex_R_uniq$RECEPTOR_2, lri_complex_R_uniq$RECEPTOR_3)
tab = table(R_units)
tab = tab[order(-tab)]

R_units_sub = substr(R_units, 1, 3)
tab_sub = table(R_units_sub)
tab_sub = tab_sub[order(-tab_sub)]

N = c()
for (sub in names(tab_sub)) {
  bool_mask = (
    (!is.na(lri_complex_R_uniq$RECEPTOR_1) & substr(lri_complex_R_uniq$RECEPTOR_1, 1, 3) == sub)
    | (!is.na(lri_complex_R_uniq$RECEPTOR_2) & substr(lri_complex_R_uniq$RECEPTOR_2, 1, 3) == sub)
    | (!is.na(lri_complex_R_uniq$RECEPTOR_3) & substr(lri_complex_R_uniq$RECEPTOR_3, 1, 3) == sub)
  )
  N = c(N, sum(bool_mask))
}
names(N) = names(tab_sub)
dt_sub = data.table(Subunit=names(N), Num_complex_R=N)
dt_sub$Perc_complex_R = dt_sub$Num_complex_R / nrow(lri_complex_R_uniq)
dt_sub = dt_sub[order(-Num_complex_R)]

## 2. Simple LRIs
lri = scDiffCom::LRI_mouse$LRI_curated
lri_simple_L = lri[is.na(LIGAND_2)]
lri_simple_R = lri[is.na(RECEPTOR_2)]

lri_simple_R_uniq = unique(lri_simple_R[, list(RECEPTOR_1)])
tab = table(lri_simple_R_uniq$RECEPTOR_1)
tab = tab[order(-tab)]

R_sub = substr(lri_simple_R_uniq$RECEPTOR_1, 1, 3)
tab_sub = table(R_sub)
tab_sub = tab_sub[order(-tab_sub)]

dt_sub = data.table(Receptor=names(tab_sub), Num=as.integer(tab_sub))
dt_sub$Perc = dt_sub$Num / nrow(lri_simple_R_uniq)
dt_sub = dt_sub[order(-Num)]
