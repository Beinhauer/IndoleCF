################################################################################
# Covariate-adjusted version of statistic.R
#
# All group comparisons (Control vs. CD vs. UC) are adjusted for Age,
# Mesalamine, Steroids and Immunosuppresants use. 
# Not adjusted (see comments): meta analysis of the covariates themselves and
# the paired within-sample tests (Wilcoxon signed-rank, McNemar).
################################################################################

library(data.table)
library(vegan)
library(car)
library(sandwich)

meta <- fread("data/meta.tsv")
meta$disease <- as.factor(meta$disease)
mic <- read.table("data/mic.table")
hrgm <- fread("data/HRGMv2_metadata_clusters.tsv")

dir.create("output/stat_adjusted", recursive = TRUE, showWarnings = FALSE)
out <- function(name) file.path("output/stat_adjusted", name)


################################################################################
##                                                                             #
##                     Covariates and helper functions                         #
##                                                                             #
################################################################################

covars <- c("Age", "Mesalamine", "Steroids", "Immunosuppresants")
meds <- c("Mesalamine", "Steroids", "Immunosuppresants")

# Missing medication is counted as "No": 22 controls (no IBD medication) and
# 3 IBD samples (PRISM_7762, PRISM_8794, PRISM_8802) with single missing values.
meta_cov <- meta[, c("mg.id", covars), with = FALSE]
for (m in meds){
  meta_cov[is.na(get(m)) | get(m) == "", (m) := "No"]
  meta_cov[, (m) := factor(get(m), levels = c("No", "Yes"))]
}
meta_cov <- meta_cov[complete.cases(meta_cov)]

# Add covariates to a table with column mg.id (complete cases only)
add_covars <- function(dt){
  dt <- as.data.table(dt)
  drop <- intersect(covars, names(dt))
  if (length(drop) > 0) dt <- dt[, !drop, with = FALSE]
  merge(dt, meta_cov, by = "mg.id")
}

# Covariates that vary in the given data (e.g. medication is constant in
# controls and can not be used within the control group)
active_covars <- function(dt){
  covars[vapply(covars, function(v) uniqueN(dt[[v]]) > 1, logical(1))]
}

adj_formula <- function(y, dt, group = "disease"){
  reformulate(c(active_covars(dt), group), response = y)
}

# Pairwise group differences from a fitted model, adjusted for the covariates
# in the model (estimate = group1 - group2 on the model scale)
pairwise_contrasts <- function(fit, V = vcov(fit), group = "disease",
                               df = fit$df.residual){
  lv <- fit$xlevels[[group]]
  b <- coef(fit)
  level_vec <- function(l){
    L <- setNames(numeric(length(b)), names(b))
    nm <- paste0(group, l)
    if (nm %in% names(b)) L[nm] <- 1
    L
  }
  res <- rbindlist(lapply(combn(lv, 2, simplify = FALSE), function(p){
    L <- level_vec(p[1]) - level_vec(p[2])
    est <- sum(L * b)
    se <- sqrt(drop(t(L) %*% V %*% L))
    stat <- est / se
    pval <- if (is.finite(df)) 2 * pt(-abs(stat), df) else 2 * pnorm(-abs(stat))
    data.table(group1 = p[1], group2 = p[2], estimate = est, se = se,
               statistic = stat, df = df, pval = pval)
  }))
  res[, padj := p.adjust(pval, method = "BH")]
  res[, signif := ifelse(padj < 0.05, "*", "ns")]
  res
}

# ANCOVA (type II) for the group comparison. A Levene test (median-centred)
# on the model residuals checks variance homogeneity between the groups:
# p >= 0.05 -> classic standard errors, p < 0.05 -> HC3-robust standard
# errors. rank-based ANCOVA used for comparisons on non-normal distributed data.
adj_lm_test <- function(dt, y, group = "disease", rank = FALSE){
  dt <- add_covars(dt)
  dt[[group]] <- droplevels(as.factor(dt[[group]]))
  if (rank){
    dt[["y_rank"]] <- frank(dt[[y]])
    y <- "y_rank"
  }
  fit <- lm(adj_formula(y, dt, group), data = dt)
  levene_p <- leveneTest(resid(fit), dt[[group]])[1, "Pr(>F)"]
  hc3 <- levene_p < 0.05
  V <- if (hc3) vcovHC(fit, type = "HC3") else vcov(fit)
  tab <- if (hc3) suppressMessages(Anova(fit, type = 2, vcov. = V)) else
    Anova(fit, type = 2)
  f_col <- grep("^F", colnames(tab), value = TRUE)[1]
  anova <- as.data.table(tab, keep.rownames = "term")
  anova[, partial_eta2 := get(f_col) * Df / (get(f_col) * Df + fit$df.residual)]
  global <- data.table(statistic = tab[group, f_col], df1 = tab[group, "Df"],
                       df2 = fit$df.residual, pval = tab[group, "Pr(>F)"],
                       partial_eta2 = anova[term == group, partial_eta2],
                       n = nrow(dt), levene_p = levene_p,
                       se = ifelse(hc3, "HC3", "classic"),
                       covariates = paste(active_covars(dt), collapse = ", "))
  pairwise <- pairwise_contrasts(fit, V = V, group = group)
  if (rank){
    pairwise[, r := statistic / sqrt(statistic^2 + df)]
  } else {
    pairwise[, cohens_d := estimate / sigma(fit)]
  }
  pairwise[, se_type := ifelse(hc3, "HC3", "classic")]
  list(global = global, pairwise = pairwise, anova = anova)
}

# Rank-based ANCOVA:
adj_rank_test <- function(dt, y, group = "disease"){
  adj_lm_test(dt, y, group, rank = TRUE)
}

# Logistic regression with likelihood-ratio test: 
# y must be 0/1. Pairwise comparisons are fitted on the two groups only;
# odds_ratio = odds of y in group1 vs. group2. "separation" marks comparisons
# where a group (or a group within a covariate level) has only 0 or only 1:
# odds ratio and CI are then unreliable, the LRT p-value is still usable.
adj_logistic_test <- function(dt, y, group = "disease"){
  dt <- add_covars(dt)
  dt[[group]] <- droplevels(as.factor(dt[[group]]))
  has_sep <- function(d, fit = NULL){
    se <- if (is.null(fit)) 0 else
      sqrt(diag(vcov(fit)))[grep(paste0("^", group), names(coef(fit)))]
    any(tapply(d[[y]], d[[group]], function(v) all(v == v[1]))) || any(se > 5)
  }
  lr_test <- function(d){
    fit <- glm(adj_formula(y, d, group), data = d, family = binomial)
    tab <- drop1(fit, scope = reformulate(group), test = "LRT")
    list(fit = fit, stat = tab[group, "LRT"], df = tab[group, "Df"],
         pval = tab[group, "Pr(>Chi)"])
  }
  g <- lr_test(dt)
  global <- data.table(statistic = g$stat, df = g$df, pval = g$pval,
                       delta_r2_mcfadden = g$stat / g$fit$null.deviance,
                       n = nrow(dt), separation = has_sep(dt, g$fit),
                       covariates = paste(active_covars(dt), collapse = ", "))
  pairwise <- rbindlist(lapply(combn(levels(dt[[group]]), 2, simplify = FALSE),
                               function(p){
    d <- dt[dt[[group]] %in% p]
    d[[group]] <- factor(d[[group]], levels = c(p[2], p[1]))
    r <- lr_test(d)
    nm <- paste0(group, p[1])
    ci <- confint.default(r$fit)[nm, ]
    data.table(group1 = p[1], group2 = p[2],
               odds_ratio = exp(coef(r$fit)[nm]), ci_low = exp(ci[1]),
               ci_high = exp(ci[2]), statistic = r$stat, pval = r$pval,
               separation = has_sep(d, r$fit))
  }))
  pairwise[, padj := p.adjust(pval, method = "BH")]
  pairwise[, signif := ifelse(padj < 0.05, "*", "ns")]
  list(global = global, pairwise = pairwise)
}

# Partial Spearman correlation of x and y, adjusted for the covariates that
# vary in dt (all variables rank-transformed, then partial correlation)
partial_spearman <- function(dt, x, y){
  dt <- add_covars(dt)
  dt <- dt[complete.cases(dt[, c(x, y), with = FALSE])]
  cv <- active_covars(dt)
  rx <- frank(dt[[x]])
  ry <- frank(dt[[y]])
  if (length(cv) > 0){
    Z <- model.matrix(reformulate(cv), data = dt)[, -1, drop = FALSE]
    Z <- apply(Z, 2, frank)
    rx <- resid(lm(rx ~ Z))
    ry <- resid(lm(ry ~ Z))
  }
  r <- cor(rx, ry)
  n <- nrow(dt)
  dfree <- n - 2 - length(cv)
  stat <- r * sqrt(dfree / (1 - r^2))
  data.table(statistic = stat, rho = round(r, digits = 3),
             pval = signif(2 * pt(-abs(stat), dfree), digits = 2), df = dfree,
             n = n, covariates = paste(cv, collapse = ", "))
}

# Effect size of the paired Wilcoxon signed-rank test: matched-pairs
# rank-biserial correlation (W+ - W-) / (W+ + W-) from V = W+ and the number
# of non-zero differences (positive: x > y)
rank_biserial <- function(V, x, y){
  n <- sum(x - y != 0)
  2 * V / (n * (n + 1) / 2) - 1
}

# Change function of pairwise.t.test to include degree of freedom and
# t statistic in output

pairwise.t.test.w.stats <- function (x, g, p.adjust.method = p.adjust.methods,
                                     pool.sd = !paired,
                                     paired = FALSE,
                                     alternative =
                                       c("two.sided", "less", "greater"),
                                     ...)
{
  if (paired && pool.sd)
    stop("pooling of SD is incompatible with paired tests")
  DNAME <- paste(deparse1(substitute(x)), "and", deparse1(substitute(g)))
  g <- factor(g)
  p.adjust.method <- match.arg(p.adjust.method)
  alternative <- match.arg(alternative)
  if (pool.sd) {
    METHOD <- "t tests with pooled SD"
    xbar <- tapply(x, g, mean, na.rm = TRUE)
    s <- tapply(x, g, sd, na.rm = TRUE)
    n <- tapply(!is.na(x), g, sum)
    degf <- n - 1
    total.degf <- sum(degf)
    pooled.sd <- sqrt(sum(ifelse(degf, s^2, 0) * degf)/total.degf)
    compare.levels.full <- function(i, j) {
      dif <- xbar[i] - xbar[j]
      se.dif <- pooled.sd * sqrt(1/n[i] + 1/n[j])
      t.val <- dif/se.dif
      p.val <- if (alternative == "two.sided")
        2 * pt(-abs(t.val), total.degf)
      else pt(t.val, total.degf, lower.tail = (alternative == "less"))
      list(t = t.val, df = total.degf, p = p.val)
    }
  }
  else {
    METHOD <- if (paired)
      "paired t tests"
    else "t tests with non-pooled SD"
    compare.levels.full <- function(i, j) {
      xi <- x[as.integer(g) == i]
      xj <- x[as.integer(g) == j]
      test <- t.test(xi, xj, paired = paired, alternative = alternative, ...)
      list(t = test$statistic, df = test$parameter, p = test$p.value)
    }
  }
  compare.p <- function(i, j) compare.levels.full(i, j)$p
  compare.t <- function(i, j) compare.levels.full(i, j)$t
  compare.df <- function(i, j) compare.levels.full(i, j)$df

  PVAL <- pairwise.table(compare.p, levels(g), p.adjust.method)
  TVAL <- pairwise.table(compare.t, levels(g), "none")
  DFVAL <- pairwise.table(compare.df, levels(g), "none")

  ans <- list(method = METHOD, p.value = PVAL, statistic = TVAL,
              parameter = DFVAL, p.adjust.method = p.adjust.method)
  class(ans) <- "pairwise.htest"
  ans
}


################################################################################
##                                                                             #
##                              Meta analysis                                  #
##                                                                             #
################################################################################

# Not adjusted: this section describes the cohort itself (Age and medications 
# considered as covariates)

# One-way ANOVA per variable and, for significant ANOVAs, pairwise t-tests with
# pooled SD and FDR correction
meta_vars <- c("Age", "Total_Reads", "unmapped")

meta_anova <- rbindlist(lapply(meta_vars, function(v){
  s <- summary(aov(meta[[v]] ~ meta$disease))[[1]]
  data.table(variable = v, statistic = s[1, "F value"], df1 = s[1, "Df"],
             df2 = s[2, "Df"], pval = s[1, "Pr(>F)"],
             eta2 = s[1, "Sum Sq"] / sum(s[, "Sum Sq"]))
}))
fwrite(meta_anova, file = out("meta_anova.tsv"), sep = "\t")

meta_pairwise <- rbindlist(lapply(meta_anova[pval < 0.05, variable], function(v){
  pw <- pairwise.t.test.w.stats(meta[[v]], meta$disease, p.adjust = "fdr",
                                paired = FALSE)
  xbar <- tapply(meta[[v]], meta$disease, mean)
  pooled_sd <- sigma(lm(meta[[v]] ~ meta$disease))
  rbindlist(lapply(combn(levels(meta$disease), 2, simplify = FALSE), function(p){
    # pairwise.table holds group2 - group1 in row p[2], column p[1]
    data.table(variable = v, group1 = p[1], group2 = p[2],
               statistic = -pw$statistic[p[2], p[1]],
               df = pw$parameter[p[2], p[1]], padj = pw$p.value[p[2], p[1]],
               cohens_d = (xbar[[p[1]]] - xbar[[p[2]]]) / pooled_sd)
  }))
}))
meta_pairwise[, signif := ifelse(padj < 0.05, "*", "ns")]
fwrite(meta_pairwise, file = out("meta_pairwise.tsv"), sep = "\t")

# Medication use per group
med_table <- meta_cov[meta[, .(mg.id, disease)], on = "mg.id"][
  , lapply(.SD, function(v) sum(v == "Yes", na.rm = TRUE)), by = disease,
  .SDcols = meds]
fwrite(med_table, file = out("medication_per_group.tsv"), sep = "\t")

# Therapy administration counts CD vs. UC: Pearson's Chi-squared test per
# medication (2x2 table, without Yates' continuity correction)
meta_ibd <- meta_cov[meta[disease %in% c("CD", "UC"), .(mg.id, disease)],
                     on = "mg.id"]
med_chisq <- rbindlist(lapply(meds, function(m){
  tab <- table(factor(meta_ibd$disease, levels = c("CD", "UC")), meta_ibd[[m]])
  ct <- chisq.test(tab, correct = FALSE)
  data.table(medication = m,
             CD_yes = tab["CD", "Yes"], CD_n = sum(tab["CD", ]),
             UC_yes = tab["UC", "Yes"], UC_n = sum(tab["UC", ]),
             statistic = unname(ct$statistic), df = unname(ct$parameter),
             pval = ct$p.value,
             cramers_v = sqrt(unname(ct$statistic) /
                                (sum(tab) * (min(dim(tab)) - 1))),
             min_expected = min(ct$expected))
}))
med_chisq[, padj := p.adjust(pval, method = "BH")]
med_chisq[, signif := ifelse(padj < 0.05, "*", "ns")]
fwrite(med_chisq, file = out("medication_chisq_CD_UC.tsv"), sep = "\t")


################################################################################
##                                                                             #
##                      Microbial diversity analysis                           #
##                                                                             #
################################################################################

#-------------------------------------------------------#
# 1. Beta-diversity as Bray-Curtis distance (with PCoA) #

bc_dist <- vegdist(mic, method = "bray")

meta_bd <- add_covars(meta[match(rownames(mic), meta$mg.id), .(mg.id, disease)])
meta_bd <- meta_bd[match(intersect(rownames(mic), meta_bd$mg.id), mg.id)]
bc_dist_bd <- as.dist(as.matrix(bc_dist)[meta_bd$mg.id, meta_bd$mg.id])

adonis_result <- adonis2(adj_formula("bc_dist_bd", meta_bd), data = meta_bd,
                         permutations = 9999, by = "margin")
fwrite(as.data.table(adonis_result, keep.rownames = "term"),
       file = out("beta_div_permanova.tsv"), sep = "\t")

groups <- c("CD", "UC", "Control")
pairs <- combn(groups, 2, simplify = FALSE)

pairwise_results <- lapply(pairs, function(g){
  idx <- meta_bd$disease %in% g
  sub_dist <- as.dist(as.matrix(bc_dist_bd)[idx, idx])
  sub_meta <- meta_bd[idx]
  sub_meta$disease <- droplevels(as.factor(sub_meta$disease))
  f <- adj_formula("sub_dist", sub_meta)
  environment(f) <- environment()
  res <- adonis2(f, data = sub_meta, permutations = 9999, by = "margin")
  data.table(group1 = g[1], group2 = g[2], pval = res["disease", "Pr(>F)"],
             f = res["disease", "F"], df = res["disease", "Df"],
             df_resid = res["Residual", "Df"], n = nrow(sub_meta),
             r2 = res["disease", "R2"])
})

pairwise_results <- rbindlist(pairwise_results)
pairwise_results[, padj := p.adjust(pval, method = "BH")]
fwrite(pairwise_results, file = out("beta_div_pval_pairwise.tsv"), sep = "\t")

# PCoA calculation (visualisation only, not adjusted).
pcoa_bc <- cmdscale(bc_dist, k=2, eig=TRUE)
pcoa_dt <- data.table(
  PC1 = pcoa_bc$points[, 1],
  PC2 = pcoa_bc$points[, 2],
  disease = meta$disease[match(rownames(mic), meta$mg.id)]
)

#-------------------------------------------------------#
#            2. Alpha-diversity                         #

alpha_div <- data.table(
  mg.id = rownames(mic),
  shannon = diversity(mic, index = "shannon"),
  observed = specnumber(mic)
)
alpha_div <- merge(alpha_div, meta, by = "mg.id")
alpha_div$disease <- as.factor(alpha_div$disease)

ad_test <- adj_rank_test(alpha_div, "shannon")
fwrite(ad_test$global, file = out("alpha_diversity_global.tsv"), sep = "\t")
fwrite(ad_test$pairwise, file = out("alpha_diversity_pval.tsv"), sep = "\t")

#-------------------------------------------------------#
#            3. Firmicutes/Bacteroidetes ratio          #

fb_ratio <- fread("output/tables/FB_ratio.tsv")
fb_ratio$disease <- as.factor(fb_ratio$disease)

fb_test <- adj_rank_test(fb_ratio, "FB_ratio")
fwrite(fb_test$global, file = out("fb_ratio_global.tsv"), sep = "\t")
fwrite(fb_test$pairwise, file = out("fb_ratio_pairwise.tsv"), sep = "\t")

#-------------------------------------------------------#
#            4. Phyla abundance                         #
phyla_per_sample <- fread("output/tables/phyla_abun_per_sample.tsv")
phyla_per_sample$disease <- as.factor(phyla_per_sample$disease)

phyla_test <- list()
global_test <- list()
for (i in unique(phyla_per_sample$Phylum)){
  sub <- phyla_per_sample[Phylum == i]
  sub$disease <- droplevels(as.factor(sub$disease))
  if (length(unique(sub$disease)) < 2) next
  res <- adj_rank_test(sub, "Abundance")
  global_test[[as.character(i)]] <- cbind(Phylum = i, res$global)
  if (res$global$pval < 0.05){
    phyla_test[[as.character(i)]] <- cbind(Phylum = i, res$pairwise)
  }
}

global_test <- rbindlist(global_test)
phyla_test <- rbindlist(phyla_test)

fwrite(phyla_test, file = out("phyla_abun_pairwise_test.tsv"), sep = "\t")
fwrite(global_test, file = out("phyla_abun_global_test.tsv"), sep ="\t")

################################################################################
##                                                                             #
##                        Indole salvager analysis                             #
##                                                                             #
################################################################################

dt_salvager_abun <- fread("output/tables/salvager_abundance.tsv")
dt_salvager_abun$disease <- as.factor(dt_salvager_abun$disease)

#------------------------------------------------------------------------------#
# 1. Total relative abundance                                                  #

salv_res <- adj_rank_test(dt_salvager_abun, "Abundance")
fwrite(salv_res$global, file = out("salvager_test_global.tsv"), sep ="\t")
fwrite(salv_res$pairwise, file = out("salvager_test.tsv"), sep ="\t")

#------------------------------------------------------------------------------#
# 2. Mean relative abundance                                                   #

salv_res_mean <- adj_rank_test(dt_salvager_abun, "Mean")
fwrite(salv_res_mean$global, file = out("salvager_test_mean_global.tsv"),
       sep ="\t")
fwrite(salv_res_mean$pairwise, file = out("salvager_test_mean.tsv"), sep ="\t")

#------------------------------------------------------------------------------#
# 3. Count - 1 % vs. 0.5 % ACO (per group)                                     #

# Not adjusted: paired comparison within the same sample, the covariates are
# identical in both measurements and can not confound the difference.
wilcox_res <- list()
for (i in unique(dt_salvager_abun$disease)){
  sub <- dt_salvager_abun[disease == i]
  wtest <- wilcox.test(sub$n_01, sub$n_005, alternative = "less", paired = TRUE,
                       conf.int = TRUE)
  wilcox_res[[i]] <- data.table(disease = i, statistic = wtest$statistic,
                                ci_low = wtest$conf.int[1],
                                ci_high = wtest$conf.int[2],
                                pval = wtest$p.value,
                                r_rb = rank_biserial(wtest$statistic, sub$n_01,
                                                     sub$n_005),
                                n = sum(sub$n_01 != sub$n_005))
}

wilcox_res <- rbindlist(wilcox_res)
fwrite(wilcox_res, file = out("count_salvager_test.tsv"), sep = "\t")

#------------------------------------------------------------------------------#
# 4. Correlation: Indole salvager ~ alpha-diversity                            #

is_ad <- merge(dt_salvager_abun, alpha_div[, .(mg.id, shannon)], by =  "mg.id",
               all.x = TRUE)
is_ad_cor <- merge(is_ad, meta[, .(mg.id, Age, unmapped)], by = "mg.id",
                   all.x = TRUE)
is_ad_cor$disease <- factor(is_ad_cor$disease,
                            levels = c("Control", "CD", "UC"))

results_all_pcor <- data.table()
for (i in levels(is_ad_cor$disease)){
  dt_spear <- cbind(disease = i,
                    partial_spearman(is_ad_cor[disease == i], "Abundance",
                                     "shannon"))
  results_all_pcor <- rbind(results_all_pcor, dt_spear)
}

fwrite(results_all_pcor, file = out("partial_cor_is_ad.tsv"), sep = "\t")

#------------------------------------------------------------------------------#
# 5. Salvaging species prevalence - Con vs. CD vs. UC                          #

salvager_prev_all <- fread("output/tables/salvager_prevalence_all.tsv")
salvager_prev_all$disease <- as.factor(salvager_prev_all$disease)
salvager_prev_all[, present := as.integer(prevalent > 0)]

prop_test_res <- list()
pairwise_res <- list()
for (i in unique(salvager_prev_all$cluster)){
  sub <- salvager_prev_all[cluster == i]
  if (uniqueN(sub$present) < 2 || uniqueN(sub$disease) < 2){
    prop_test_res[[as.character(i)]] <- data.table(cluster = i, pval = NA_real_)
    next
  }
  res <- adj_logistic_test(sub, "present")
  prop_test_res[[as.character(i)]] <- cbind(cluster = i, res$global)
  pairwise_res[[as.character(i)]] <- cbind(cluster = i, res$pairwise)
}

prop_test_res <- rbindlist(prop_test_res, fill = TRUE)
prop_test_res[, padj := p.adjust(pval, method = "BH")]

# Pairwise comparisons only for clusters with significant global test,
# BH-adjusted across all of these comparisons
sig_cluster <- prop_test_res[!is.na(padj) & padj < 0.05, cluster]
pairwise_res <- rbindlist(pairwise_res)[cluster %in% sig_cluster]
pairwise_res[, padj := p.adjust(pval, method = "BH")]
pairwise_res[, signif := ifelse(padj < 0.05, "*", "ns")]

pairwise_res <- merge(pairwise_res, hrgm[, .(`HRGMv2 Cluster`, Species)],
                      by.x = "cluster", by.y = "HRGMv2 Cluster", all.x = TRUE)

fwrite(prop_test_res, file = out("IS_species_prev_test.tsv"), sep ="\t")
fwrite(pairwise_res, file = out("IS_species_pairwise_prev_test.tsv"), sep ="\t")

################################################################################
##                                                                             #
##                        Indole producer analysis                             #
##                                                                             #
################################################################################

dt_prod_abun <- fread("output/tables/producer_abudance.tsv")
dt_prod_abun$disease <- as.factor(dt_prod_abun$disease)

#------------------------------------------------------------------------------#
# 1. Total relative abundance                                                  #

prod_res <- adj_rank_test(dt_prod_abun, "Abundance")
fwrite(prod_res$global, file = out("producer_test_global.tsv"), sep="\t")
fwrite(prod_res$pairwise, file = out("producer_test.tsv"), sep="\t")

#------------------------------------------------------------------------------#
# 2. Mean relative abundance                                                   #

prod_res_mean <- adj_rank_test(dt_prod_abun, "Mean")
fwrite(prod_res_mean$global, file = out("producer_test_mean_global.tsv"),
       sep="\t")
fwrite(prod_res_mean$pairwise, file = out("producer_test_mean.tsv"), sep="\t")


################################################################################
##                                                                             #
##                      Indole cross-feeding analysis                          #
##                                                                             #
################################################################################

dt_long <- fread("output/tables/crossfeeding_scen0.tsv")

#------------------------------------------------------------------------------#
# 1. Normalized indole cross-feeding flux - Con vs. CD vs. UC (per ACO)        #

global_res <- list()
pair_res <- list()

for (h in unique(dt_long$abun_co)){
  sub <- dt_long[abun_co == h & cf == TRUE]
  sub$disease <- as.factor(sub$disease)
  if (uniqueN(sub$disease) > 1){
    res <- adj_rank_test(sub, "crossfeeding")
    global_res[[as.character(h)]] <- cbind(abun_co = h, res$global)
    if (res$global$pval < 0.05){
      pair_res[[as.character(h)]] <- cbind(abun_co = h, res$pairwise)
    }
  }
}

global_res <- rbindlist(global_res)
pair_res <- rbindlist(pair_res, fill = TRUE)

fwrite(global_res, file = out("crossfeeding_global_test.tsv"), sep = "\t")
fwrite(pair_res, file = out("crossfeeding_pairwise_test.tsv"), sep = "\t")

#------------------------------------------------------------------------------#
# 2. Normalized indole cross-feeding flux - 1% vs. 0.5% ACO (per group)        #

# Not adjusted: paired comparison within the same sample (see above).
wilcox_res <- list()
for (k in unique(dt_long$disease)){
  wide <- dcast(dt_long[disease == k], mg.id ~ abun_co,
                value.var = "crossfeeding")
  cn <- setdiff(colnames(wide), "mg.id")
  if (length(cn) == 2){
    keep <- wide[[cn[1]]] != 0 & wide[[cn[2]]] != 0
    wtest <- wilcox.test(wide[[cn[1]]][keep], wide[[cn[2]]][keep], paired = TRUE)
    wilcox_res[[k]] <- data.table(disease = k, v = wtest$statistic,
                                  pval = wtest$p.value,
                                  r_rb = rank_biserial(wtest$statistic,
                                                       wide[[cn[1]]][keep],
                                                       wide[[cn[2]]][keep]),
                                  n = sum(keep))
  }
}

wilcox_res <- rbindlist(wilcox_res)
wilcox_res$disease <- as.factor(wilcox_res$disease)
wilcox_res$disease <- factor(wilcox_res$disease,
                             levels = c("Control", "CD", "UC"))

setorder(wilcox_res, disease)

cf_test_wilx <- data.table(
  disease = wilcox_res$disease,
  group1 = cn[1],
  group2 = cn[2],
  v = wilcox_res$v,
  n = wilcox_res$n,
  pval = wilcox_res$pval,
  r_rb = wilcox_res$r_rb,
  signif = ifelse(wilcox_res$pval < 0.05, "*", "ns"))

fwrite(cf_test_wilx, file = out("wilcoxon_crossfeeding_abuntest.tsv"),
       sep = "\t" )

#------------------------------------------------------------------------------#
# 3. Indole cross-feeding prevalence - Con vs. CD vs. UC (per ACO)             #

prop_test_res <- list()
prop_pair_res <- list()
for (i in unique(dt_long$abun_co)){
  sub <- dt_long[abun_co == i]
  sub[, cf_present := as.integer(as.logical(cf))]
  if (uniqueN(sub$cf_present) < 2 || uniqueN(sub$disease) < 2){
    prop_test_res[[as.character(i)]] <- data.table(abun_co = i, pval = NA_real_)
    next
  }
  res <- adj_logistic_test(sub, "cf_present")
  prop_test_res[[as.character(i)]] <- cbind(abun_co = i, res$global)
  prop_pair_res[[as.character(i)]] <- cbind(abun_co = i, res$pairwise)
}

prop_test_res <- rbindlist(prop_test_res, fill = TRUE)
prop_pair_res <- rbindlist(prop_pair_res)

fwrite(prop_test_res, file = out("icf_prevalence_logistic.tsv"), sep ="\t")
fwrite(prop_pair_res, file = out("icf_prevalence_logistic_pairwise.tsv"),
       sep ="\t")

# Prevalence per group and ACO for plotting (descriptive, not adjusted) with
# exact (Clopper-Pearson) 95% confidence interval
icf_prev <- dt_long[, .(n = .N, n_cf = sum(as.logical(cf))),
                    by = .(disease, abun_co)]
icf_prev[, c("prevalence", "ci_low", "ci_high") := {
  ci <- mapply(function(x, n) binom.test(x, n)$conf.int, n_cf, n)
  list(n_cf / n * 100, ci[1, ] * 100, ci[2, ] * 100)
}]
icf_prev[, disease := factor(disease, levels = c("Control", "CD", "UC"))]
setorder(icf_prev, abun_co, disease)
fwrite(icf_prev, file = "output/tables/icf_prevalence_ci.tsv", sep = "\t")

#------------------------------------------------------------------------------#
# 4. Indole cross-feeding prevalence- 1% vs. 0.5% ACO (per group)              #

# Not adjusted: paired comparison within the same sample (see above).
mcnemar_res <- list()
for (j in unique(dt_long$disease)){
  wide <- dcast(dt_long[disease == j], mg.id ~ abun_co, value.var = "cf")
  cn <- setdiff(colnames(wide), "mg.id")
  res <- tryCatch({
    tab <- table(factor(wide[[cn[1]]], levels = c(TRUE, FALSE)),
                 factor(wide[[cn[2]]], levels = c(TRUE, FALSE)))
    mt <- mcnemar.test(tab)
    discordant <- c(tab[1, 2], tab[2, 1])
    data.table(disease = j,  statistic = mt$statistic, df = mt$parameter,
               pval = mt$p.value,
               cohens_g = abs(discordant[1] / sum(discordant) - 0.5))
  }, error = function(e) {
    data.table(disease = j, pval = message(e))
  })
  mcnemar_res[[j]] <- res
}

mcnemar <- rbindlist(mcnemar_res)
fwrite(mcnemar, file = out("icf_prevalence_mcnemar.tsv"), sep ="\t")

#------------------------------------------------------------------------------#
# 5. Correlation: Indole cross-feeding ~ indole salvager abundance             #

results_pc <- data.table()
for (i in unique(dt_long$disease)){
  cf <- dt_long[disease == i]
  is <- dt_salvager_abun[disease == i]
  for (j in unique(cf$abun_co)){
    cf_aco <- cf[abun_co == j]
    dt <- merge(cf_aco, is[, .(mg.id, Abundance)], by = "mg.id")
    dt_spear <- cbind(disease = i, abun_co = j,
                      partial_spearman(dt, "crossfeeding", "Abundance"))
    results_pc <- rbind(results_pc, dt_spear)
  }
}

fwrite(results_pc, file = out("partial_cor_IS_ICF.tsv"), sep = "\t")

################################################################################
##                                                                             #
##                     Indole derivative producer analysis                     #
##                                                                             #
################################################################################

ipa_dt_prod_abun <- fread("output/tables/ipa_producer_abundance.tsv")
ipa_dt_prod_abun$disease <- as.factor(ipa_dt_prod_abun$disease)
ila_dt_prod_abun <- fread("output/tables/ila_producer_abundance.tsv")
ila_dt_prod_abun$disease <- as.factor(ila_dt_prod_abun$disease)
iaa_dt_prod_abun <- fread("output/tables/iaa_producer_abundance.tsv")
iaa_dt_prod_abun$disease <- as.factor(iaa_dt_prod_abun$disease)
ratio_dt <- fread("output/tables/ratio_salvager_id_producer.tsv")
ratio_dt$disease <- as.factor(ratio_dt$disease)

#------------------------------------------------------------------------------#
# 1. IPA producer                                                              #

ipa_res <- adj_lm_test(ipa_dt_prod_abun, "Abundance")
fwrite(ipa_res$global, file = out("ipa_producer_test_global.tsv"), sep ="\t")
fwrite(ipa_res$pairwise, file = out("ipa_producer_test.tsv"), sep ="\t")

#------------------------------------------------------------------------------#
# 2. ILA producer                                                              #

ila_res <- adj_lm_test(ila_dt_prod_abun, "Abundance")
fwrite(ila_res$global, file = out("ila_producer_test_global.tsv"), sep ="\t")
fwrite(ila_res$pairwise, file = out("ila_producer_test.tsv"), sep ="\t")

#------------------------------------------------------------------------------#
# 3. IAA producer                                                              #

iaa_res <- adj_rank_test(iaa_dt_prod_abun, "Abundance")
fwrite(iaa_res$global, file = out("iaa_producer_test_global.tsv"), sep ="\t")
fwrite(iaa_res$pairwise, file = out("iaa_producer_test.tsv"), sep ="\t")

#------------------------------------------------------------------------------#
# 4. Correlation: Indole salvager  ~ Indole derivative producer                #

global_all <- data.table()
pair_all <- data.table()

for (i in unique(ratio_dt$type)){
  res <- adj_rank_test(ratio_dt[type == i], "ratio")
  global_all <- rbind(global_all, cbind(type = i, res$global))
  pw <- cbind(type = i, res$pairwise)
  if (res$global$pval > 0.05){
    pw[, c("statistic", "pval", "padj", "r") := NA_real_]
    pw[, signif := "ns"]
  }
  pair_all <- rbind(pair_all, pw)
}

fwrite(global_all, file = out("ratio_salvager_id_prod_global_test.tsv"),
       sep ="\t")
fwrite(pair_all, file = out("ratio_salvager_id_prod_pairwise_test.tsv"),
       sep ="\t")

################################################################################
##                                                                             #
##                             SCFA analysis                                   #
##                                                                             #
################################################################################

scfa <- fread("output/tables/scfa_production.tsv")

#------------------------------------------------------------------------------#
# 1. Acetate                                                                   #

a_res <- adj_lm_test(scfa[react_name == "Acetate"], "norm_flux")

fwrite(a_res$anova, file = out("acetate_test.tsv"), sep = "\t")
fwrite(a_res$global, file = out("acetate_test_global.tsv"), sep = "\t")
fwrite(a_res$pairwise, file = out("acetate_test_pairwise.tsv"), sep = "\t")

#------------------------------------------------------------------------------#
# 2. Butyrate and propionate                                                   #

global_res <- list()
pair_res <- list()

for (h in unique(scfa[react_name != "Acetate"]$react_name)){
  sub <- scfa[react_name == h]
  sub$disease <- as.factor(sub$disease)
  if (uniqueN(sub$disease) > 1){
    res <- adj_rank_test(sub, "norm_flux")
    global_res[[as.character(h)]] <- cbind(react_name = h, res$global)
    if (res$global$pval < 0.05){
      pair_res[[as.character(h)]] <- cbind(react_name = h, res$pairwise)
    }
  }
}

global_res <- rbindlist(global_res)
pair_res <- rbindlist(pair_res, fill = TRUE)

fwrite(global_res, file = out("butyrate_propionate_global_test.tsv"),
       sep = "\t")
fwrite(pair_res, file = out("butyrate_propionate_test.tsv"), sep = "\t")


#------------------------------------------------------------------------------#
# 3. SCFA production prevalence - Con vs. CD vs. UC (per SCFA)                 #

# CoM.R only stores fluxes with abs(flux) > 1e-6, so samples without production
# are missing from scfa. Here they are added with norm_flux = 0.
scfa_zero <- CJ(mg.id = meta$mg.id, react_name = unique(scfa$react_name))
scfa_zero <- merge(scfa_zero, scfa[, .(mg.id, react_name, norm_flux)],
                   by = c("mg.id", "react_name"), all.x = TRUE)
scfa_zero[is.na(norm_flux), norm_flux := 0]
scfa_zero <- merge(scfa_zero, meta[, .(mg.id, disease)], by = "mg.id")
scfa_zero[, producing := as.integer(norm_flux > 0)]

# Prevalence per group for plotting (descriptive, not adjusted) with exact
# (Clopper-Pearson) 95% confidence interval
scfa_prev <- scfa_zero[, .(n = .N, n_producing = sum(producing)),
                       by = .(react_name, disease)]
scfa_prev[, c("prevalence", "ci_low", "ci_high") := {
  ci <- mapply(function(x, n) binom.test(x, n)$conf.int, n_producing, n)
  list(n_producing / n * 100, ci[1, ] * 100, ci[2, ] * 100)
}]
setorder(scfa_prev, react_name, disease)
fwrite(scfa_prev, file = "output/tables/scfa_prevalence.tsv", sep = "\t")

prev_global <- list()
prev_pair <- list()
for (h in unique(scfa_zero$react_name)){
  res <- adj_logistic_test(scfa_zero[react_name == h], "producing")
  prev_global[[h]] <- cbind(react_name = h, res$global)
  if (res$global$pval < 0.05){
    prev_pair[[h]] <- cbind(react_name = h, res$pairwise)
  }
}

fwrite(rbindlist(prev_global), file = out("scfa_prevalence_test_global.tsv"),
       sep = "\t")
fwrite(rbindlist(prev_pair), file = out("scfa_prevalence_test_pairwise.tsv"),
       sep = "\t")
