library(data.table)
library(ggplot2)
library(Seurat)
library(fitdistrplus)

set.seed(42)
n_rep <- 10000

## function to simulation permutation #####

get_true_score <- function(mly, mlo, mry, mro, score_type) {
  if (score_type == "ARI") {
    return((mlo + mro - mly - mry)/2)
  } else if (score_type == "GDIFF") {
    return((sqrt(mlo*mro)-sqrt(mly*mry)))
  } else if (score_type == "GLOGFC") {
    return(log(sqrt(mlo*mro)/sqrt(mly*mry)))
  } else {
    stop("error")
  }
}

get_shuffeld_score <- function(ly, lo, ry, ro, score_type, stat_type) {
  if (stat_type == "rep") {
    if (length(ly) > length(lo)) {
      lo <- c(lo, rep(sum(lo)/(length(ly) - length(lo)), length(ly) - length(lo)))
    }
    if (length(ly) < length(lo)) {
      ly <- c(ly, rep(sum(ly)/(length(lo) - length(ly)), length(lo) - length(ly)))
    }
    if (length(ry) > length(ro)) {
      ro <- c(ro, rep(sum(ro)/(length(ry) - length(ro)), length(ry) - length(ro)))
    }
    if (length(ry) < length(ro)) {
      ry <- c(ry, rep(sum(ry)/(length(ro) - length(ry)), length(ro) - length(ry)))
    }
  }
  if (stat_type == "boot") {
    if (length(ly) > length(lo)) {
      lo <- c(lo, sample(lo, length(ly) - length(lo), replace = TRUE))
    }
    if (length(ly) < length(lo)) {
      ly <- c(ly, sample(ly, length(lo) - length(ly), replace = TRUE))
    }
    if (length(ry) > length(ro)) {
      ro <- c(ro, sample(ro, length(ry) - length(ro), replace = TRUE))
    }
    if (length(ry) < length(ro)) {
      ry <- c(ry, sample(ry, length(ro) - length(ry), replace = TRUE))
    }
  }
  if (stat_type == "min") {
    if (length(ly) > length(lo)) {
      ly <- sample(ly, length(lo))
    }
    if (length(ly) < length(lo)) {
      lo <- sample(lo, length(ly))
    }
    if (length(ry) > length(ro)) {
      ry <- sample(ry, length(ro))
    }
    if (length(ry) < length(ro)) {
      ro <- sample(ro, length(ry))
    }
  }
  sample_l <- sample(c(ly, lo))
  sample_r <- sample(c(ry, ro))
  mly <- mean(sample_l[1:length(ly)])
  mlo <- mean(sample_l[(1+length(ly)):(length(ly)+length(lo))])
  mry <- mean(sample_r[1:length(ry)])
  mro <- mean(sample_r[(1+length(ry)):(length(ry)+length(ro))])
  get_true_score(mly, mlo, mry, mro, score_type)
}

## main function that does everything ####

get_everything <- function(n_rep, ly, lo, ry, ro)
  {
  score_types <- c("ARI", "GDIFF", "GLOGFC")
  stat_types <- c("basic", "rep", "boot", "min")
  dt_scores <- data.table(
    score_type = score_types,
    score =  sapply(
      score_types,
      function(score_type) {
            get_true_score(mean(ly), mean(lo), mean(ry), mean(ro), score_type)
          }
    )
  )
  dt_distr <- as.data.table(
    sapply(
      score_types,
      function(score_type) {
        sapply(
          stat_types,
          function(stat_type) {
            replicate(n_rep, get_shuffeld_score(ly, lo, ry, ro, score_type, stat_type))
          },
          simplify = "array"
        )
      },
      simplify = "array"
    )
  )
  setnames(dt_distr, old = c("V1", "V2", "V3", "value"), new = c("rep_id", "stat_type", "score_type", "score"))
  dt_distr[, type := paste(score_type, stat_type, sep = "_")]
  dt_summary <- dt_distr[,  .(
    p_value = sum(abs(score) > abs(dt_scores[score_type == .BY$score_type]$score))/n_rep,
    mean = mean(score),
    sd = sd(score),
    ratio = mean(score)/sd(score)
    ), by = c("stat_type", "score_type")]
  return(list(score = dt_scores, distr = dt_distr, summary = dt_summary))
}

## Real case of interest ####

seurat_droplet_spleen <- readRDS("../data_scAgeCom/test/inputs/seurat_testing_tms_droplet_spleen.rds")
seurat_droplet_spleen <- NormalizeData(seurat_droplet_spleen, assay = "RNA")
seurat_droplet_spleen$age_group <- ifelse(seurat_droplet_spleen$age %in% c('1m', '3m'), 'YOUNG', 'OLD' )
seurat_droplet_spleen$cell_ontology_class <- as.character(seurat_droplet_spleen$cell_ontology_class)

Hmgb1:Tlr2

hist(expm1(as.vector(subset(seurat_droplet_spleen, subset = cell_ontology_class == "plasma cell" & age_group == "YOUNG")$RNA@data["Hmgb1", ])))
hist(expm1(as.vector(subset(seurat_droplet_spleen, subset = cell_ontology_class == "plasma cell" & age_group == "OLD")$RNA@data["Hmgb1", ])))

hist(expm1(as.vector(subset(seurat_droplet_spleen, subset = cell_ontology_class == "granulocyte" & age_group == "YOUNG")$RNA@data["Tlr2", ])))
hist(expm1(as.vector(subset(seurat_droplet_spleen, subset = cell_ontology_class == "granulocyte" & age_group == "OLD")$RNA@data["Tlr2", ])))


ly1 <- expm1(as.vector(subset(seurat_droplet_spleen, subset = cell_ontology_class == "plasma cell" & age_group == "YOUNG")$RNA@data["Hmgb1", ]))
lo1 <- expm1(as.vector(subset(seurat_droplet_spleen, subset = cell_ontology_class == "plasma cell" & age_group == "OLD")$RNA@data["Hmgb1", ]))
ry1 <- expm1(as.vector(subset(seurat_droplet_spleen, subset = cell_ontology_class == "granulocyte" & age_group == "YOUNG")$RNA@data["Tlr2", ]))
ro1 <- expm1(as.vector(subset(seurat_droplet_spleen, subset = cell_ontology_class == "granulocyte" & age_group == "OLD")$RNA@data["Tlr2", ]))

fitdistr(subset(seurat_droplet_spleen, subset = cell_ontology_class == "plasma cell" & age_group == "YOUNG")$RNA@data[4, ], 
         "negative binomial")

## Case 1 ####

ly1 <- rep(1, 10)
lo1 <- rep(5, 100)
ry1 <- rep(0.5, 15)
ro1 <- rep(1.5, 50)

ly1 <- rnbinom(192, mu = 2265.6, size = 0.0533)
lo1 <- rnbinom(192, mu = 273.7, size = 0.0363)
ry1 <- rnbinom(37, mu = 2101.26, size = 0.0247)
ro1 <- rnbinom(37, mu = 84.19, size = 0.0177)

ly1 <- rnbinom(1000, mu = 5, size = 0.2)
lo1 <- rnbinom(10, mu = 5, size = 0.2)
ry1 <- rnbinom(2000, mu = 2, size = 0.5)
ro1 <- rnbinom(20, mu = 2, size = 0.5)

test <- c(ly1, lo1)
test

sum(test == 0)


ly1 <- rnorm(1000, mean = 5, sd = 0.2)
lo1 <- rnorm(50, mean = 5, sd = 0.2)
ry1 <- rnorm(2000, mean = 2, sd = 0.5)
ro1 <- rnorm(50, mean = 2, sd = 0.5)

hist(ry1, breaks = 50)

res_1 <- get_everything(1000, ly1, lo1, ry1, ro1)

res_1$summary
res_1$score

ggplot(res_1$distr[score_type == "ARI" & stat_type %in% c("basic")], aes(x = score, color = type, fill = type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50, fill = "white") +
  geom_density(alpha = 0.6)

ggplot(res_1$distr[score_type == "GDIFF" & stat_type %in% c("basic")], aes(x = score, color = type, fill = type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50, fill = "white") +
  geom_density(alpha = 0.6)

ggplot(res_1$distr[score_type == "ARI" & stat_type %in% c("basic", "rep", "boot", "min")], aes(x = score, color = type, fill = type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50, fill = "white") +
  geom_density(alpha = 0.6)


res_1$score
res_1$summary

ggplot(res_1$distr, aes(x = score, color = type, fill = type)) +
  geom_vline(xintercept = res_1$score$score) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50, fill = "white") +
  geom_density(alpha = 0.6)

ggplot(res_1$distr[score_type == "ARI" & stat_type %in% c("rep", "boot", "basic", "min")], aes(x = score, color = type, fill = type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50, fill = "white") +
  geom_density(alpha = 0.6)

ggplot(res_1$distr[score_type == "ARI" & stat_type %in% c("rep", "boot", "basic", "min")], aes(x = score, color = type, fill = type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50, fill = "white") +
  geom_density(alpha = 0.6)

ggplot(res_1$distr[score_type == "GDIFF" & stat_type %in% c("rep", "boot", "basic", "min")], aes(x = score, color = type, fill = type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50, fill = "white") +
  geom_density(alpha = 0.6)

ggplot(res_1$distr[score_type == "GLOGFC" & stat_type %in% c("rep", "boot", "basic")], aes(x = score, color = type, fill = type)) +
  geom_vline(xintercept = res_1$score[score_type == "GLOGFC"]$score) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50, fill = "white") +
  geom_density(alpha = 0.6) 

## Case 2 ####

ly2 <- rnorm(10, mean = 1, sd = 0.1)
lo2 <- rnorm(100, mean = 5, sd = 0.1)
ry2 <- rnorm(200, 0.5, 0.1)
ro2 <- rnorm(15, 1.5, 0.1)

res_2 <- get_everything(10000, ly2, lo2, ry2, ro2)

res_2$score
res_2$summary

ggplot(res_2$distr, aes(x = score, color = type, fill = type)) +
  geom_vline(xintercept = res_2$score$score) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50, fill = "white") +
  geom_density(alpha = 0.6)

ggplot(res_2$distr[score_type == "ARI" & stat_type %in% c("rep", "boot", "basic")], aes(x = score, color = type, fill = type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50, fill = "white") +
  geom_density(alpha = 0.6)

ggplot(res_2$distr[score_type == "ARI" & stat_type %in% c("rep", "boot", "basic")], aes(x = score, color = type, fill = type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50, fill = "white") +
  geom_density(alpha = 0.6)

ggplot(res_2$distr[score_type == "GDIFF" & stat_type %in% c("rep", "boot", "basic")], aes(x = score, color = type, fill = type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50, fill = "white") +
  geom_density(alpha = 0.6)

ggplot(res_2$distr[score_type == "GLOGFC" & stat_type %in% c("rep", "boot", "basic")], aes(x = score, color = type, fill = type)) +
  geom_vline(xintercept = res_2$score[score_type == "GLOGFC"]$score) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50, fill = "white") +
  geom_density(alpha = 0.6) 











#2




ggplot(res_2, aes(x = score, color = score_type, fill = score_type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50, fill = "white") +
  geom_density(alpha = 0.6)





#3

ly3 <- rnorm(100, mean = 1, sd = 0.1)
lo3 <- rnorm(100, mean = 5, sd = 0.1)
ry3 <- rnorm(50, 0.5, 0.1)
ro3 <- rnorm(50, 1.5, 0.1)

get_true_score(mean(ly3), mean(lo3), mean(ry3), mean(ro3), "ARI")
get_true_score(mean(ly3), mean(lo3), mean(ry3), mean(ro3), "GDIFF")
get_true_score(mean(ly3), mean(lo3), mean(ry3), mean(ro3), "GLOGFC")

res_3 <- data.table(
  score = c(
    replicate(n_rep, get_shuffeld_score(ly3, lo3, ry3, ro3, "ARI")),
    replicate(n_rep, get_shuffeld_score(ly3, lo3, ry3, ro3, "GDIFF")),
    replicate(n_rep, get_shuffeld_score(ly3, lo3, ry3, ro3, "GLOGFC"))
  ),
  score_type = c(rep("ARI", n_rep), rep("GDIFF", n_rep), rep("GLOGFC", n_rep))
)

ggplot(res_3, aes(x = score, color = score_type, fill = score_type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50, fill = "white") +
  geom_density(alpha = 0.6)

sum(abs(res_3[score_type == "ARI"]$score) > abs(get_true_score(mean(ly3), mean(lo3), mean(ry3), mean(ro3), "ARI")))/n_rep
sum(abs(res_3[score_type == "GDIFF"]$score) > abs(get_true_score(mean(ly3), mean(lo3), mean(ry3), mean(ro3), "GDIFF")))/n_rep
sum(abs(res_3[score_type == "GLOGFC"]$score) > abs(get_true_score(mean(ly3), mean(lo3), mean(ry3), mean(ro3), "GLOGFC")))/n_rep

res_3[, mean(score), by = score_type]
res_3[, sd(score), by = score_type]
res_3[, mean(score)/sd(score), by = score_type]

#4

ly4 <- rnorm(100, mean = 5, sd = 0.1)
lo4 <- rnorm(100, mean = 1, sd = 0.1)
ry4 <- rnorm(50, 0.5, 0.1)
ro4 <- rnorm(50, 1.5, 0.1)

get_true_score(mean(ly4), mean(lo4), mean(ry4), mean(ro4), "ARI")
get_true_score(mean(ly4), mean(lo4), mean(ry4), mean(ro4), "GDIFF")
get_true_score(mean(ly4), mean(lo4), mean(ry4), mean(ro4), "GLOGFC")

res_4 <- data.table(
  score = c(
    replicate(n_rep, get_shuffeld_score(ly4, lo4, ry4, ro4, "ARI")),
    replicate(n_rep, get_shuffeld_score(ly4, lo4, ry4, ro4, "GDIFF")),
    replicate(n_rep, get_shuffeld_score(ly4, lo4, ry4, ro4, "GLOGFC"))
  ),
  score_type = c(rep("ARI", n_rep), rep("GDIFF", n_rep), rep("GLOGFC", n_rep))
)

ggplot(res_4, aes(x = score, color = score_type, fill = score_type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50, fill = "white") +
  geom_density(alpha = 0.6)

sum(abs(res_4[score_type == "ARI"]$score) > abs(get_true_score(mean(ly4), mean(lo4), mean(ry4), mean(ro4), "ARI")))/n_rep
sum(abs(res_4[score_type == "GDIFF"]$score) > abs(get_true_score(mean(ly4), mean(lo4), mean(ry4), mean(ro4), "GDIFF")))/n_rep
sum(abs(res_4[score_type == "GLOGFC"]$score) > abs(get_true_score(mean(ly4), mean(lo4), mean(ry4), mean(ro4), "GLOGFC")))/n_rep

res_4[, mean(score), by = score_type]
res_4[, sd(score), by = score_type]
res_4[, mean(score)/sd(score), by = score_type]

#5
ly5 <- rnorm(10, mean = 5, sd = 0.1)
lo5 <- rnorm(100, mean = 1, sd = 0.1)
ry5 <- rnorm(15, 0.5, 0.1)
ro5 <- rnorm(50, 1.5, 0.1)

get_true_score(mean(ly5), mean(lo5), mean(ry5), mean(ro5), "ARI")
get_true_score(mean(ly5), mean(lo5), mean(ry5), mean(ro5), "GDIFF")
get_true_score(mean(ly5), mean(lo5), mean(ry5), mean(ro5), "GLOGFC")

res_5 <- data.table(
  score = c(
    replicate(n_rep, get_shuffeld_score(ly5, lo5, ry5, ro5, "ARI")),
    replicate(n_rep, get_shuffeld_score(ly5, lo5, ry5, ro5, "GDIFF")),
    replicate(n_rep, get_shuffeld_score(ly5, lo5, ry5, ro5, "GLOGFC"))
  ),
  score_type = c(rep("ARI", n_rep), rep("GDIFF", n_rep), rep("GLOGFC", n_rep))
)

ggplot(res_5, aes(x = score, color = score_type, fill = score_type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50, fill = "white") +
  geom_density(alpha = 0.6)

sum(abs(res_5[score_type == "ARI"]$score) > abs(get_true_score(mean(ly5), mean(lo5), mean(ry5), mean(ro5), "ARI")))/n_rep
sum(abs(res_5[score_type == "GDIFF"]$score) > abs(get_true_score(mean(ly5), mean(lo5), mean(ry5), mean(ro5), "GDIFF")))/n_rep
sum(abs(res_5[score_type == "GLOGFC"]$score) > abs(get_true_score(mean(ly5), mean(lo5), mean(ry5), mean(ro5), "GLOGFC")))/n_rep

res_5[, mean(score), by = score_type]
res_5[, sd(score), by = score_type]
res_5[, mean(score)/sd(score), by = score_type]

#6
ly6 <- rnorm(100, mean = 5, sd = 0.1)
lo6 <- rnorm(100, mean = 1, sd = 0.1)
ry6 <- rnorm(50, 0.5, 0.1)
ro6 <- rnorm(50, 1.5, 0.1)

get_true_score(mean(ly6), mean(lo6), mean(ry6), mean(ro6), "ARI")
get_true_score(mean(ly6), mean(lo6), mean(ry6), mean(ro6), "GDIFF")
get_true_score(mean(ly6), mean(lo6), mean(ry6), mean(ro6), "GLOGFC")

res_6 <- data.table(
  score = c(
    replicate(n_rep, get_shuffeld_score(ly6, lo6, ry6, ro6, "ARI")),
    replicate(n_rep, get_shuffeld_score(ly6, lo6, ry6, ro6, "GDIFF")),
    replicate(n_rep, get_shuffeld_score(ly6, lo6, ry6, ro6, "GLOGFC"))
  ),
  score_type = c(rep("ARI", n_rep), rep("GDIFF", n_rep), rep("GLOGFC", n_rep))
)

ggplot(res_6, aes(x = score, color = score_type, fill = score_type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.6, bins = 60, fill = "white") +
  geom_density(alpha = 0.6)

sum(abs(res_6[score_type == "ARI"]$score) > abs(get_true_score(mean(ly6), mean(lo6), mean(ry6), mean(ro6), "ARI")))/n_rep
sum(abs(res_6[score_type == "GDIFF"]$score) > abs(get_true_score(mean(ly6), mean(lo6), mean(ry6), mean(ro6), "GDIFF")))/n_rep
sum(abs(res_6[score_type == "GLOGFC"]$score) > abs(get_true_score(mean(ly6), mean(lo6), mean(ry6), mean(ro6), "GLOGFC")))/n_rep

res_6[, mean(score), by = score_type]
res_6[, sd(score), by = score_type]
res_6[, mean(score)/sd(score), by = score_type]

#7
ly7 <- rnorm(100, mean = 5, sd = 1)
lo7 <- rnorm(100, mean = 1, sd = 0.2)
ry7 <- rnorm(50, 0.5, 0.2)
ro7 <- rnorm(50, 1.5, 0.4)

get_true_score(mean(ly7), mean(lo7), mean(ry7), mean(ro7), "ARI")
get_true_score(mean(ly7), mean(lo7), mean(ry7), mean(ro7), "GDIFF")
get_true_score(mean(ly7), mean(lo7), mean(ry7), mean(ro7), "GLOGFC")

res_7 <- data.table(
  score = c(
    replicate(n_rep, get_shuffeld_score(ly7, lo7, ry7, ro7, "ARI")),
    replicate(n_rep, get_shuffeld_score(ly7, lo7, ry7, ro7, "GDIFF")),
    replicate(n_rep, get_shuffeld_score(ly7, lo7, ry7, ro7, "GLOGFC"))
  ),
  score_type = c(rep("ARI", n_rep), rep("GDIFF", n_rep), rep("GLOGFC", n_rep))
)

ggplot(res_7, aes(x = score, color = score_type, fill = score_type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.7, bins = 70, fill = "white") +
  geom_density(alpha = 0.7)

sum(abs(res_7[score_type == "ARI"]$score) > abs(get_true_score(mean(ly7), mean(lo7), mean(ry7), mean(ro7), "ARI")))/n_rep
sum(abs(res_7[score_type == "GDIFF"]$score) > abs(get_true_score(mean(ly7), mean(lo7), mean(ry7), mean(ro7), "GDIFF")))/n_rep
sum(abs(res_7[score_type == "GLOGFC"]$score) > abs(get_true_score(mean(ly7), mean(lo7), mean(ry7), mean(ro7), "GLOGFC")))/n_rep

res_7[, mean(score), by = score_type]
res_7[, sd(score), by = score_type]
res_7[, mean(score)/sd(score), by = score_type]
