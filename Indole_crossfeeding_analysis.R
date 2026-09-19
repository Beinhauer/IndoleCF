# ANALYSIS OF INDOLE CROSS-FEEDING:
# cross-feeding calculation for all scenarios
# sensitivity towards adjustment of trp uptake costs, pFBA coefficient and
# coupling constraints
# Community simulations using different scenarios performed with CoM.R

library(data.table)
library(dunn.test)

meta <- fread("data/meta.tsv")
all_ids <- unique(meta$mg.id)

##---------- scenario -------------##
scen <- c("scen0", "scen1", "scen2", "scen3", "scen4", "scen5")

kw_all_res <- data.table()
dunn_all_res <- data.table()
wilcox_all_res <- data.table()

for (i in unique(scen)){
  dt_01 <- fread(
    paste0("results/ICF_0.01/", i, "/all_samples_combined_", i, ".tsv"))
  dt_005 <- fread(
    paste0("results/ICF_0.005/all_samples_combined", i, ".tsv"))

  indole_cf_01 <- dt_01[grepl("cpd00359", react),
                        .(mg.id, name, flux_indole = flux, growth)]
  indole_cf_005 <- dt_005[grepl("cpd00359", react),
                          .(mg.id, name, flux_indole = flux, growth)]
  indole_cf_01 <- indole_cf_01[, norm_flux := flux_indole / growth ]
  indole_cf_005 <- indole_cf_005[, norm_flux := flux_indole / growth ]

  cf_per_sample_01 <- indole_cf_01[, .(
    sum_prod = sum(norm_flux[norm_flux > 0], na.rm=TRUE),
    sum_cons = abs(sum(norm_flux[norm_flux < 0], na.rm =TRUE))
  ), by = .(mg.id)]
  cf_per_sample_005 <- indole_cf_005[, .(
    sum_prod = sum(norm_flux[norm_flux > 0], na.rm=TRUE),
    sum_cons = abs(sum(norm_flux[norm_flux < 0], na.rm =TRUE))
  ), by = .(mg.id)]
  cf_per_sample_01[, crossfeeding := pmin(sum_prod, sum_cons)]
  cf_per_sample_005[, crossfeeding := pmin(sum_prod, sum_cons)]
  cf_per_sample_01[, abun_co := 1]
  cf_per_sample_005[, abun_co := 0.5]
  cf_per_sample_comb <- rbind(cf_per_sample_005, cf_per_sample_01)
  cf_per_sample_comb$abun_co <- as.factor(cf_per_sample_comb$abun_co)

  full_grid <- CJ(mg.id = all_ids,
                  abun_co = unique(cf_per_sample_comb$abun_co))

  dt_long <- merge(full_grid, cf_per_sample_comb, by = c("mg.id", "abun_co"),
                   all.x = TRUE)
  dt_long[is.na(crossfeeding), crossfeeding := 0]
  dt_long <- merge(dt_long, meta[, .(mg.id, disease)], by = "mg.id")
  dt_long[, cf := crossfeeding != 0]

  # write table with crossfeeding fluxes per sample for different scenarios
  fwrite(dt_long, file = paste0("output/tables/crossfeeding_", i,".tsv"),
         sep = "\t")

  # kruskal-wallis test with dunn's post hoc test
  kw_res <- list()
  dunn_res <- list()
  for (h in unique(dt_long$abun_co)){
    sub <- dt_long[abun_co == h & cf == TRUE]
    sub$disease <- as.factor(sub$disease)
    if (uniqueN(sub$disease) > 1){
      kw <- kruskal.test(crossfeeding ~ disease, data = sub)
      grp_means <- sub[, .(mean_cf = mean(crossfeeding),
                           sd_cf = sd(crossfeeding)),
                       by = disease]
      kw_res[[as.character(h)]] <- data.table(abun_co = h, pval = kw$p.value,
                                              disease = grp_means$disease,
                                              mean_cf = grp_means$mean_cf,
                                              sd_cf = grp_means$sd_cf)
    if (kw$p.value < 0.05){
      dunn <- dunn.test(x = sub$crossfeeding, g = sub$disease,
                        method = "bh", alpha = 0.05, altp = TRUE)
      grp_means <- sub[, .(mean_cf = mean(crossfeeding),
                           sd_cf = sd(crossfeeding)), by = disease]
      dt_dunn <- data.table(comparison = dunn$comparisons, Z = dunn$Z,
                            pval = dunn$altP, padj = dunn$altP.adjusted)
      dt_dunn[, `:=`(mean_g1 =
                       grp_means$mean_cf[match(tstrsplit(comparison," - ")[[1]],
                                                       grp_means$disease)],
                     mean_g2 =
                       grp_means$mean_cf[match(tstrsplit(comparison, " - ")[[2]],
                                                       grp_means$disease)])]
      dt_dunn[, abun_co := h]
       dunn_res[[as.character(h)]] <- dt_dunn
      }
    }
  }
  kw_res <- rbindlist(kw_res, fill = TRUE)
  kw_res[, scen := i]
  kw_all_res <- rbind(kw_all_res, kw_res, fill = TRUE)
  dunn_res <- rbindlist(dunn_res, fill = TRUE)
  dunn_res[, scen := i]
  dunn_all_res <- rbind(dunn_all_res, dunn_res, fill = TRUE)

  # wilcoxon rank sum test
  wilcox_res <- list()
  for (k in unique(dt_long$disease)){
    wide <- dcast(dt_long[disease == k], mg.id ~ abun_co,
                  value.var = "crossfeeding")
    cn <- setdiff(colnames(wide), "mg.id")
    if (length(cn) == 2){
      keep <- wide[[cn[1]]] != 0 & wide[[cn[2]]] != 0
      wtest <- wilcox.test(wide[[cn[1]]][keep], wide[[cn[2]]][keep],
                           paired = TRUE)
      wilcox_res[[k]] <- data.table(disease = k, pval = wtest$p.value,
                                    n = sum(keep),
                                    mean_co1 = mean(wide[[cn[1]]][keep]),
                                    sd_co1 = sd(wide[[cn[1]]][keep]),
                                    mean_co2 = mean(wide[[cn[2]]][keep]),
                                    sd_co2 = sd(wide[[cn[2]]][keep]))
      setnames(wilcox_res[[k]], c("mean_co1","mean_co2"), paste0("mean_", cn))
      setnames(wilcox_res[[k]], c("sd_co1","sd_co2"), paste0("sd_", cn))
    }
  }
  wilcox_res <- rbindlist(wilcox_res, fill = TRUE)
  wilcox_res[, scen := i]
  wilcox_all_res <- rbind(wilcox_all_res, wilcox_res, fill = TRUE)
}

fwrite(kw_all_res, file = "output/tables/sensitivity_analysis_KW.tsv",
       sep = "\t")
fwrite(dunn_all_res, file = "output/tables/sensitivity_analysis_DN.tsv",
       sep = "\t")
fwrite(wilcox_all_res, file = "output/tables/sensitivity_analysis_WX.tsv",
       sep = "\t")


# calculate indole cross-feeding prevalene per disease and abundance cut-off
# for standard scenario

dt_long <- fread("output/tables/crossfeeding_scen0.tsv")
cf_prev_prop <- dt_long[, .(prop = mean(cf)*100), by =c("disease", "abun_co")]

fwrite(cf_prev_prop, file = "output/tables/prevalence_ICF.tsv", sep = "\t")

