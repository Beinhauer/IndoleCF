library(data.table)
library(dunn.test)
library(vegan)
library(rstatix)

meta <- fread("data/meta.tsv")
meta$disease <- as.factor(meta$disease)
mic <- read.table("data/mic.table")


dir.create("output/stat")



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

meta$disease <- as.factor(meta$disease)
anova_age <- aov(meta$Age ~ meta$disease)
summary(anova_age)
pw_age <- pairwise.t.test.w.stats(meta$Age, meta$disease, p.adjust="fdr",
                                  paired = FALSE)

anova_tR <- aov(meta$Total_Reads ~ meta$disease)
summary(anova_tR)

anova_unmap <- aov(meta$unmapped ~ meta$disease)
summary(anova_unmap)
pw_unm <- pairwise.t.test.w.stats(meta$unmapped, meta$disease, p.adjust="fdr",
                                  paired = FALSE)


################################################################################
##                                                                             #
##                      Microbial diversity analysis                           #
##                                                                             #
################################################################################

#-------------------------------------------------------#
# 1. Beta-diversity as Bray-Curtis distance (with PCoA) #

bc_dist <- vegdist(mic, method = "bray")

meta_bd <- meta[match(rownames(mic), meta$mg.id)]
adonis_result <- adonis2(bc_dist ~ disease, data=meta_bd, permutations = 9999)

groups <- c("CD", "UC", "Control")
pairs <- combn(groups, 2, simplify = FALSE)

pairwise_results <- lapply(pairs, function(g){
  idx <- meta_bd$disease %in% g
  sub_dist <- as.dist(as.matrix(bc_dist)[idx, idx])
  sub_meta <- meta_bd[idx]
  sub_meta$disease <- as.factor(sub_meta$disease)
  res <- adonis2(sub_dist ~ disease, data = sub_meta, permutations = 9999)
  data.table(group1 = g[1], group2 = g[2], pval= res$`Pr(>F)`[1],
             f = res$`F`[1],df = res$Df)
})

pairwise_results <- rbindlist(pairwise_results)
pairwise_results[, padj := p.adjust(pval, method = "BH")]
fwrite(pairwise_results, file = "output/stat/beta_div_pval_pairwise.tsv",
       sep = "\t")

# PCoA calculation
pcoa_bc <- cmdscale(bc_dist, k=2, eig=TRUE)
pcoa_dt <- data.table(
  PC1 = pcoa_bc$points[, 1],
  PC2 = pcoa_bc$points[, 2],
  disease = meta$disease
)
fwrite(pcoa_dt, file = "output/tables/pcoa.tsv", sep = "\t")

#-------------------------------------------------------#
#            2. Alpha-diversity                         #

alpha_div <- data.table(
  mg.id = rownames(mic),
  shannon = diversity(mic, index = "shannon"),
  observed = specnumber(mic)
)
alpha_div <- merge(alpha_div, meta, by = "mg.id")
alpha_div$disease <- as.factor(alpha_div$disease)

fwrite(alpha_div, file = "output/tables/alpha_diversity.tsv", sep = "\t")

kw_ad <- kruskal.test(shannon ~ disease, data = alpha_div)
ad_dunn <- dunn.test(alpha_div$shannon, alpha_div$disease, method = "bh",
                     altp = TRUE)
p_values <- data.table(group1 = c("CD", "CD", "Control"),
                       group2 = c("Control", "UC", "UC"),
                       p_value = ad_dunn$altP.adjusted)
p_values[, annotations := ifelse(p_value < 0.05, "*", "ns")]
fwrite(p_values, file = "output/stat/alpha_diversity_pval.tsv", sep = "\t")

#-------------------------------------------------------#
#            3. Firmicutes/Bacteroidetes ratio          #

fb_ratio <- fread("output/tables/FB_ratio.tsv")
fb_ratio$disease <- as.factor(fb_ratio$disease)

kw <- kruskal.test(FB_ratio ~ disease, data = fb_ratio)

#-------------------------------------------------------#
#            4. Phyla abundance                         #
phyla_per_sample <- fread("output/tables/phyla_abun_per_sample.tsv")
phyla_per_sample$disease <- as.factor(phyla_per_sample$disease)

phyla_test <- list()
kw_test <- list()
for (i in unique(phyla_per_sample$Phylum)){
  sub <- phyla_per_sample[Phylum == i]
  sub$disease <- droplevels(as.factor(sub$disease))
  if (length(unique(sub$disease)) < 2) next
  kw <- kruskal.test(sub$Abundance ~ sub$disease)
  kw_test[[as.character(i)]] <- data.table(Phylum = i, pval = kw$p.value,
                                           statistic = kw$statistic,
                                           df = kw$parameter)
  if (kw$p.value < 0.05){
    dn <- dunn.test(sub$Abundance, sub$disease, method = "bh", altp = TRUE)
    phyla_test[[as.character(i)]] <- data.table(Phylum = i,
                                                Comparison = dn$comparisons,
                                                pval = dn$altP.adjusted,
                                                zstatistic = dn$Z)
    next
  }
}

kw_test <- rbindlist(kw_test)
phyla_test <- rbindlist(phyla_test)
phyla_test[, signif := ifelse(pval < 0.05, "*", "ns")]

fwrite(phyla_test, file = "output/tables/phyla_abun_dn_test.tsv", sep = "\t")
fwrite(kw_test, file = "output/tables/phyla_abun_kw_test.tsv", sep ="\t")

################################################################################
##                                                                             #
##                        Indole salvager analysis                             #
##                                                                             #
################################################################################

dt_salvager_abun <- fread("output/tables/salvager_abundance.tsv")
dt_salvager_abun$disease <- as.factor(dt_salvager_abun$disease)

#------------------------------------------------------------------------------#
# 1. Total relative abundance                                                  #

kw <- kruskal.test(dt_salvager_abun$Abundance ~ dt_salvager_abun$disease)
dn <- as.data.table(dunn.test(dt_salvager_abun$Abundance,
                              dt_salvager_abun$disease,  method = "bh",
                              altp= TRUE))

salvager_test <- dn[, .(group1 = c("CD", "CD", "Control"),
                        group2 = c("Control", "UC", "UC"),
                        z = dn$Z, pval = dn$altP.adjusted,
                        signif = ifelse(dn$altP.adjusted < 0.05, "*", "ns"))]

write(salvager_test, file = "output/stat/salavager_test.tsv", sep ="\t")

#------------------------------------------------------------------------------#
# 2. Mean relative abundance                                                   #

kw_mean <- kruskal.test(dt_salvager_abun$Mean ~ dt_salvager_abun$disease)
dn_mean <- as.data.table(dunn.test(dt_salvager_abun$Mean,
                                   dt_salvager_abun$disease,method = "bh",
                                   altp = TRUE))

salvager_test_mean <- dn_mean[, .(group1 = c("CD", "CD", "Control"),
                                  group2 = c("Control", "UC", "UC"),
                                  z = dn_mean$Z, pval = dn_mean$altP.adjusted,
                                  signif = ifelse(
                                    dn_mean$altP.adjusted < 0.05, "*", "ns"))]

fwrite(salvager_test_mean, file = "output/stat/salavager_test_mean.tsv",
       sep ="\t")

#------------------------------------------------------------------------------#
# 3. Count - 1 % vs. 0.5 % ACO (per group)                                     #

wilcox_res <- list()
for (i in unique(dt_salvager_abun$disease)){
  sub <- dt_salvager_abun[disease == i]
  wtest <- wilcox.test(sub$n_01, sub$n_005, alternative = "less", paired = TRUE,
                       conf.int = TRUE)
  wilcox_res[[i]] <- data.table(disease = i, statistic = wtest$statistic,
                                ci = wtest$conf.int, pval = wtest$p.value)
}

wilcox_res <- rbindlist(wilcox_res)
fwrite(wilcox_res, file = "output/stat/count_salvager_test.tsv", sep = "\t")

#------------------------------------------------------------------------------#
# 4. Correlation: Indole salvager ~ alpha-diversity                            #

is_ad <- merge(dt_salvager_abun, alpha_div[, .(mg.id, shannon)], by =  "mg.id",
               all.x = TRUE)
is_ad_cor <- merge(is_ad, meta[, .(mg.id, Age, unmapped)], by = "mg.id",
                   all.x = TRUE)
is_ad_cor$disease <- as.factor(is_ad_cor$disease)
is_ad_cor$disease <- factor(is_ad_cor$disease,
                            levels = c("Control", "CD", "UC"))
fwrite(is_ad_cor, file = "output/tables/spearman_salvagager_abun_adiv.tsv",
       sep ="\t")

results_all_pcor <- data.table()
for (i in unique(is_ad_cor$disease)){
  dt <- is_ad_cor[disease == i]
  spear <- cor.test(dt$Abundance, dt$shannon, method = "spearman")
  dt_spear <- data.table(disease = i, statistic = spear$statistic,
                         rho = round(spear$estimate, digits = 3),
                         pval = signif(spear$p.value, digits = 2),
                         df = spear$parameter)
  results_all_pcor <- rbind(results_all_pcor, dt_spear)
}

fwrite(results_all_pcor, file = "output/stat/partial_cor_is_ad.tsv", sep = "\t")

#------------------------------------------------------------------------------#
# 5. Salvaging species prevalence - Con vs. CD vs. UC                          #

salvager_prev_all <- fread("output/tables/salvager_prevalence_all.tsv")
salvager_prev_all$disease <- as.factor(salvager_prev_all$disease)

prop_test_res <- list()
for (i in unique(salvager_prev_all$cluster)){
  tab <- table(factor(salvager_prev_all[cluster == i, disease],
                      levels = unique(salvager_prev_all$disease)),
               factor(salvager_prev_all[cluster == i, prevalent],
                      levels = c(0, 100)))
  if (any(rowSums(tab)==0) || any(colSums(tab)==0)){
    prop_test_res[[as.character(i)]] <- data.table(cluster = i, pval = NA_real_)
    next
  }
  ft <- fisher.test(tab, simulate.p.value = FALSE)
  prop_test_res[[as.character(i)]] <- data.table(cluster = i, pval = ft$p.value)
}

prop_test_res <- rbindlist(prop_test_res)
prop_test_res[, padj := p.adjust(pval, method = "BH")]

sig_cluster <- prop_test_res[!is.na(padj) & padj < 0.05, cluster]
pairs <- combn(unique(salvager_prev_all$disease), 2, simplify = FALSE)

pairwise_res <- list()
k <- 1

for (i in sig_cluster){
  for (p in pairs){
    sub <- salvager_prev_all[cluster == i & disease %in% p]
    tab <- table(factor(sub$disease, levels = p),
                 factor(sub$prevalent, levels = c(0, 100)))
    if (any(rowSums(tab) == 0) || any(colSums(tab) == 0 )){
      pairwise_res[[k]] <- data.table(cluster = i, group1 = p[1], group2 = p[2],
                                      pval = NA_real_)
      k <- k+1
      next
    }

    ft <- fisher.test(tab, conf.int = TRUE)
    pairwise_res[[k]] <- data.table(cluster = i, group1 = p[1], group2 = p[2],
                                    pval = ft$p.value, estimate = ft$estimate,
                                    ci_low = ft$conf.int[1],
                                    ci_high = ft$conf.int[2])
    k <- k+1
  }
}
pairwise_res <- rbindlist(pairwise_res)
pairwise_res[, padj := p.adjust(pval, method = "BH")]

pairwise_res <- merge(pairwise_res, hrgm[, .(`HRGMv2 Cluster`, Species)],
                      by.x = "cluster", by.y = "HRGMv2 Cluster", all.x = TRUE)
pairwise_res[, signif := ifelse(padj < 0.05, "*", "ns")]

fwrite(pairwise_res, file = "output/stat/IS_species_pairwise_prev_test.tsv",
       sep ="\t")

################################################################################
##                                                                             #
##                        Indole producer analysis                             #
##                                                                             #
################################################################################

dt_prod_abun <- fread("output/tables/producer_abudance.tsv")
dt_prod_abun$disease <- as.factor(dt_prod_abun$disease)

#------------------------------------------------------------------------------#
# 1. Total relative abundance                                                  #

kw <- kruskal.test(dt_prod_abun$Abundance ~ dt_prod_abun$disease)
prod_test <- data.table(group1 = c("CD", "CD", "Control"),
                        group2 = c("Control", "UC", "UC"),
                        signif = c("ns", "ns", "ns"))
fwrite(prod_test, file = "output/stat/producer_test.tsv", sep="\t")

#------------------------------------------------------------------------------#
# 2. Mean relative abundance                                                   #

kw_mean <- kruskal.test(dt_prod_abun$Mean ~ dt_prod_abun$disease)
prod_test_mean <- data.table(group1 = c("CD", "CD", "Control"),
                             group2 = c("Control", "UC", "UC"),
                             signif = c("ns", "ns", "ns"))
fwrite(prod_test_mean, file = "output/sta/producer_test_mean.tsv", sep="\t")


################################################################################
##                                                                             #
##                      Indole cross-feeding analysis                          #
##                                                                             #
################################################################################

dt_long <- fread("output/tables/crossfeeding_scen2.tsv")

#------------------------------------------------------------------------------#
# 1. Normalized indole cross-feeding flux - Con vs. CD vs. UC (per ACO)        #

kw_res <- list()
dunn_res <- list()

for (h in unique(dt_long$abun_co)){
  sub <- dt_long[abun_co == h & cf == TRUE]
  sub$disease <- as.factor(sub$disease)
  if (uniqueN(sub$disease) > 1){
    kw <- kruskal.test(crossfeeding ~ disease, data = sub)
    kw_res[[as.character(h)]] <- data.table(abun_co = h, pval = kw$p.value,
                                            statistic = kw$statistic,
                                            df = kw$parameter)
    if (kw$p.value < 0.05){
      dunn <- dunn.test(x = sub$crossfeeding, g = sub$disease, method = "bh",
                        alpha = 0.05, altp = TRUE)
      dt_dunn <- data.table(comparison = dunn$comparisons,
                            Z = dunn$Z,
                            pval = dunn$altP,
                            padj = dunn$altP.adjusted)
      dt_dunn[, abun_co := h]
      dunn_res[[as.character(h)]] <- dt_dunn
    }
  }
}

kw_res <- rbindlist(kw_res)
dunn_res <- rbindlist(dunn_res, fill = TRUE)

cf_test_dn <- data.table(abun_co = c("1", "1", "1"),
                         group1 = c("CD", "CD", "Control"),
                         group2 = c("Control", "UC", "UC"), z = dunn_res$Z,
                         padj = dunn_res$padj,
                         signif = ifelse(dunn_res$padj <= 0.05, "*", "ns"))

fwrite(cf_test_dn, file = "output/stat/dunn_crossfeeding_1abunco.tsv",
       sep = "\t")

#------------------------------------------------------------------------------#
# 2. Normalized indole cross-feeding flux - 1% vs. 0.5% ACO (per group)        #

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
                                  ci = wtest$conf.int, n = sum(keep))
  }
}

wilcox_res <- rbindlist(wilcox_res)
wilcox_res$disease <- as.factor(wilcox_res$disease)
wilcox_res$disease <- factor(wilcox_res$disease,
                             levels = c("Control", "CD", "UC"))

cf_test_wilx <- data.table(
  disease = c("Control", "CD", "UC"),
  group1 = c("1", "1", "1"),
  group2 = c("0.5", "0.5", "0.5"),
  pval = wilcox_res$pval,
  signif = ifelse(wilcox_res$pval < 0.5, "*", "ns"))

fwrite(cf_test_wilx, file = "output/stat/wilcoxon_crossfeeding_abuntest.tsv",
       sep = "\t" )

#------------------------------------------------------------------------------#
# 3. Indole cross-feeding prevalence - Con vs. CD vs. UC (per ACO)             #

prop_test_res <- list()
for (i in unique(dt_long$abun_co)){
  tab <- table(factor(dt_long[abun_co == i, disease],
                      levels = unique(dt_long$disease)),
               factor(dt_long[abun_co == i, cf], levels = c(TRUE, FALSE)))
  if (any(rowSums(tab)==0) || any(colSums(tab)==0)){
    prop_test_res[[as.character(i)]] <- data.table(abun_co = i, pval = NA_real_)
    next
  }
  ft <- fisher.test(tab, simulate.p.value = FALSE, conf.int = TRUE)
  prop_test_res[[as.character(i)]] <- data.table(abun_co = i,
                                                 pval = ft$p.value,
                                                 estimate = ft$estimate,
                                                 ci = ft$conf.int)
}

prop_test_res <- rbindlist(prop_test_res)


fwrite(prop_test_res, file ="output/stat/icf_prevalence_fisher.tsv", sep ="\t")

#------------------------------------------------------------------------------#
# 4. Indole cross-feeding prevalence- 1% vs. 0.5% ACO (per group)              #

mcnemar_res <- list()
for (j in unique(dt_long$disease)){
  wide <- dcast(dt_long[disease == j], mg.id ~ abun_co, value.var = "cf")
  cn <- setdiff(colnames(wide), "mg.id")
  res <- tryCatch({
    tab <- table(factor(wide[[cn[1]]], levels = c(TRUE, FALSE)),
                 factor(wide[[cn[2]]], levels = c(TRUE, FALSE)))
    mt <- mcnemar.test(tab)
    data.table(disease = j,  statistic = mt$statistic, df = mt$parameter,
               pval = mt$p.value)
  }, error = function(e) {
    data.table(disease = j, pval = message(e))
  })
  mcnemar_res[[j]] <- res
}

mcnemar <- rbindlist(mcnemar_res)
fwrite(mcnemar, file ="output/stat/icf_prevalence_mcnemar.tsv", sep ="\t")

#------------------------------------------------------------------------------#
# 5. Correlation: Indole cross-feeding ~ indole salvager abundance             #

abun_salv_pc <- merge(dt_salvager_abun, meta[, .(mg.id, Age, unmapped)],
                      by = "mg.id", all.x = TRUE)

results_pc <- data.table()
for (i in unique(dt_long$disease)){
  cf <- dt_long[disease == i]
  is <- abun_salv_pc[disease == i]
  for (j in unique(cf$abun_co)){
    cf_aco <- cf[abun_co == j]
    dt <- merge(cf_aco, is[, .(mg.id, Abundance, Age, unmapped)], by = "mg.id")
    spear <- cor.test(dt$crossfeeding, dt$Abundance,
                       method = "spearman")
    dt_spear <- data.table(disease = i, abun_co = j,
                           statistic = spear$statistic,
                           df = spear$parameter,
                           rho = round(spear$estimate, digits = 3),
                           pval = signif(spear$p.value, digits = 2))
    results_pc <- rbind(results_pc, dt_spear)
  }
}

fwrite(results_pc, file = "output/stat/partial_cor_IS_ICF.tsv", sep = "\t")

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

ipa_welch_anova <- oneway.test(Abundance ~ disease, data = ipa_dt_prod_abun,
                               var.equal = FALSE)
ipa_ph <- games_howell_test(ipa_dt_prod_abun, Abundance ~ disease,
                            conf.level = 0.95, detailed = TRUE)
fwrite(ipa_ph, file = "output/tables/ipa_producer_test_gh.tsv", sep ="\t")

#------------------------------------------------------------------------------#
# 2. ILA producer                                                              #

ila_welch_anova <- oneway.test(Abundance ~ disease, data = ila_dt_prod_abun,
                               var.equal = FALSE)
ila_ph <- games_howell_test(ila_dt_prod_abun, Abundance ~ disease,
                            conf.level = 0.95, detailed = TRUE)
fwrite(ila_ph, file = "output/tables/ila_producer_test_gh.tsv", sep ="\t")

#------------------------------------------------------------------------------#
# 3. IAA producer                                                              #

iaa_kw <- kruskal.test(iaa_dt_prod_abun$Abundance ~ iaa_dt_prod_abun$disease)
iaa_dn <- as.data.table(dunn.test(iaa_dt_prod_abun$Abundance,
                                  iaa_dt_prod_abun$disease, method = "bh",
                                  altp = TRUE))
iaa_test <- iaa_dn[, .(group1 = c("CD", "CD", "Control"),
                       group2 = c("Control", "UC", "UC"),
                       pval = iaa_dn$altP.adjusted,
                       signif = ifelse(iaa_dn$altP.adjusted < 0.05, "*", "ns"))]
fwrite(iaa_test, file = "output/tables/iaa_producer_test.tsv", sep ="\t")

#------------------------------------------------------------------------------#
# 4. Correlation: Indole salvager  ~ Indole derivative producer                #

kw_all <- data.table()
dn_all <- data.table()

for (i in unique(ratio_dt$type)){
  kw <- kruskal.test(ratio ~ disease, data = ratio_dt[type == i])
  kw_dt <- data.table(type = i, pval = kw$p.value, statistic = kw$statistic,
                      df = kw$parameter)
  kw_all <- rbind(kw_all, kw_dt)
  if (kw$p.value <= 0.05){
    dn <- dunn.test(ratio_dt[type == i]$ratio, ratio_dt[type == i]$disease,
                    method = "bh", altp = TRUE)
    dn_dt <- data.table(type = i, group1 = c("CD", "CD", "Control"),
                        group2 = c("Control", "UC", "UC"), z = dn$Z,
                        pval = dn$altP.adjusted,
                        signif = ifelse(dn$altP.adjusted < 0.05, "*", "ns"))
    dn_all <- rbind(dn_all, dn_dt)
  } else {
    dn_dt <- data.table(type = i, group1 = c("CD", "CD", "Control"),
                        group2 = c("Control", "UC", "UC"),
                        z = NA, pval = NA, signif = "ns")
    dn_all <- rbind(dn_all, dn_dt)
  }
}

fwrite(kw_all, file = "output/tables/ratio_salvager_id_prod_kw_test.tsv",
       sep ="\t")
fwrite(dn_all, file = "output/tables/ratio_salvager_id_prod_dn_test.tsv",
       sep ="\t")

################################################################################
##                                                                             #
##                             SCFA analysis                                   #
##                                                                             #
################################################################################

scfa <- fread("output/tables/scfa_production.tsv")

#------------------------------------------------------------------------------#
# 1. Acetate                                                                   #

linmod <- lm(norm_flux ~ disease, data = scfa[react_name == "Acetate"])
a_anova <- anova(linmod)

fwrite(a_anova, file = "output/stat/acetate_test.tsv", sep = "\t")

#------------------------------------------------------------------------------#
# 2. Butyrate and propionate                                                   #

kw_res <- list()
dunn_res <- list()

for (h in unique(scfa[react_name != "Acetate"]$react_name)){
  sub <- scfa[react_name == h]
  sub$disease <- as.factor(sub$disease)
  if (uniqueN(sub$disease) > 1){
    kw <- kruskal.test(norm_flux ~ disease, data = sub)
    kw_res[[as.character(h)]] <- data.table(react_name = h, pval = kw$p.value,
                                            statistic = kw$statistic,
                                            df = kw$parameter)
    if (kw$p.value < 0.05){
      dunn <- dunn.test(x = sub$norm_flux, g = sub$disease, method = "bh",
                        alpha = 0.05, altp = TRUE)
      dt_dunn <- data.table(comparison = dunn$comparisons,
                            Z = dunn$Z,
                            pval = dunn$altP,
                            padj = dunn$altP.adjusted)
      dt_dunn[, react_name := h]
      dunn_res[[as.character(h)]] <- dt_dunn
    }
  }
}

kw_res <- rbindlist(kw_res)
dunn_res <- rbindlist(dunn_res, fill = TRUE)

fwrite(dunn_res, file = "output/stat/butyrate_propionate_test.tsv", sep = "\t")

