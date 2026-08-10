# Construction and simulation of microbial community models

library(qs2)
library(cobrar)
library(cobrarCPLEX)
library(data.table)
library(stringr)

COBRAR_SETTINGS("SOLVER","cplex")

hrgm <- fread("data/HRGMv2_metadata_clusters.tsv")

models <- qs_read("data/models.qs2")
mic <- read.table("data/mic.table")
meta <- fread("data/meta.tsv")

# all simulations repeated in a sensitivity analysis with the following parameters
# standard/scenario 0: bwcost trp 400, pFBAcoeff 1e-3, cpl_c 50, cpl_u 0.001
# scenario 1: bwcost trp 200, -//-
# scenario 2: bwcost trp 600, -//-
# scenario 3: bwcost trp 400, pFBAcoeff 0.5e-3, -//-
# scenario 4: bwcost trp 400, pFBAcoeff 1.5e-3, -//-
# scenario 5: bwcost trp 400, pFBAcoeff 1e-3, cpl_c 500, cpl_c 0.01

##### ------------------------------ Simulation with abundance cut-off at 0.005 ----------------------------------------#############

abun_co <- 0.005

#########################################################################################
##
##  1. Identification of minimal lower bounds throughout all individual microbial
##  models to create comparable nutritional environment
##
#########################################################################################

# 1.1 Create files with exchange reactions for all communities

for (i in seq_along(meta$mg.id)) {
  mgx_spl <- meta$mg.id[i]
  rel_models_idx <- which(mic[mgx_spl,] >= abun_co)
  rel_models_dt  <- data.table(model = colnames(mic)[rel_models_idx],
                               model_idx = rel_models_idx,
                               rel_abun = unlist(mic[mgx_spl, rel_models_idx]))
  rel_models_dt[, rel_abun := rel_abun / sum(rel_abun)]

  cmod <- joinModels(models[rel_models_dt$model], mergeLB = "min",
                     abun = rel_models_dt$rel_abun)

  cmod <- fixBMRatios(cmod)
  cmod <- fluxBMCoupling(cmod, cpl_c = 50, cpl_u = 0.001)
  nutri <- subset(findExchReact(cmod), exchange_type=="Community" & grepl("^EX_",react_id))
  nutri <- merge(nutri, mgx_spl)
  fwrite (nutri, file = paste0("results/Nutrition_0.005/", mgx_spl,".tsv"), sep = "\t")
}

# 1.2 extract minimum lower bound for each reaction

ExchR <- list.files(path = "results/Nutrition_0.005/", pattern = "\\.tsv", full.names=TRUE)
all_ExchR <- as.data.table(rbindlist(lapply(ExchR, fread), fill=TRUE))
fwrite(all_ExchR, file = paste0("results/Nutrition_0.005/exchangeReactions.tsv"), sep = "\t")

minlb <- all_ExchR[, .SD[which.min(lb)], by = react_id]
minlb <- minlb[, .(react_id = gsub("^EX_|_e0|_c0$", "", react_id), met_name, lb)]
fwrite(minlb, file = paste0("results/Nutrition_0.005/minLowerBounds.tsv"), sep ="\t")

##########################################################################
##
##  2. Simulation
##
##########################################################################

minlb <- fread("results/Nutrition_0.005/minLowerBounds.tsv")
minlb <- as.data.table(minlb)

for (i in seq_along(meta$mg.id)) {
  mgx_spl <- meta$mg.id[i]
  rel_models_idx <- which(mic[mgx_spl,] >= abun_co)
  rel_models_dt  <- data.table(model = colnames(mic)[rel_models_idx],
                               model_idx = rel_models_idx,
                               rel_abun = unlist(mic[mgx_spl, rel_models_idx]))
  rel_models_dt[, rel_abun := rel_abun / sum(rel_abun)]

  cmod <- joinModels(models[rel_models_dt$model], mergeLB = "min",
                     abun = rel_models_dt$rel_abun)
  cmod <- fixBMRatios(cmod)
  cmod <- fluxBMCoupling(cmod, cpl_c = 50, cpl_u = 0.001)
  # set nutritional environment
  for (j in seq_len(nrow(minlb))) {
    cpd_id <- paste0("EX_", minlb$react_id[j], "_e0")
    if (cpd_id %in% cmod@react_id) {
      cmod <- changeBounds(cmod, react = cpd_id, lb = minlb$lb[j])
    } else {
      warning("reaction not found: ", cpd_id)
    }
  }
  # increase cost of trp uptake and weight of flux minimization
  bwcosts <- rep(1, react_num(cmod))
  bwcosts[grep("^M[0-9]+_EX_cpd00065", cmod@react_id)] <- 400
  sol <- pfbaHeuristic(cmod, costcoefbw = bwcosts, pFBAcoeff = 1e-3)
  dtsol <- data.table(growth = sol@obj,
                      flux = sol@fluxes,
                      react = cmod@react_id,
                      react_name = cmod@react_name)
  dtsol[, index := str_extract(string = react, pattern = "^M[0-9]+")]
  dtsol[, mg.id := mgx_spl]
  dtsol <- merge(dtsol, cmod@community[, c("index","name")], all.x = TRUE)
  dtsol <- merge(dtsol, meta[mg.id == mgx_spl, .(mg.id, disease)], by = "mg.id", all.x = TRUE)
  dtsol <- dtsol[grepl("bio1|cpd00359|cpd27289|cpd00065|cpd00054|cpd00211|cpd00029|cpd00141|rxn01964|rxn00474|rxn00473|rxn08751", react)]
  dtsol <- dtsol[abs(flux) > 1e-6]
  fwrite (dtsol, file = paste0("results/Icf_Im_pM_abun0.005/scen/", mgx_spl, ".tsv"), sep = "\t")
}

# merge into one tsv
files <- list.files(path = "results/Icf_Im_pM_abun0.005/scen/", pattern = "\\.tsv", full.names=TRUE)
all_data <- rbindlist(lapply(files, fread), fill=TRUE)
fwrite (all_data, file = paste0("results/Icf_Im_pM_abun0.005/scen/all_samples_combined.tsv"), sep = "\t")


#################################################################################################################################
##### ------------------------------ Simulation with abundance cut-off at 0.01 ----------------------------------------#############

abun_co <- 0.01

#########################################################################################
##
##  1. Identification of minimal lower bounds throughout all individual microbial
##  models to create comparable nutritional environment
##
#########################################################################################

# 1.1 Create files with exchange reactions for all communities

for (i in seq_along(meta$mg.id)) {
  mgx_spl <- meta$mg.id[i]
  rel_models_idx <- which(mic[mgx_spl,] >= abun_co)
  rel_models_dt  <- data.table(model = colnames(mic)[rel_models_idx],
                               model_idx = rel_models_idx,
                               rel_abun = unlist(mic[mgx_spl, rel_models_idx]))
  rel_models_dt[, rel_abun := rel_abun / sum(rel_abun)]

  cmod <- joinModels(models[rel_models_dt$model], mergeLB = "min",
                     abun = rel_models_dt$rel_abun)

  cmod <- fixBMRatios(cmod)
  cmod <- fluxBMCoupling(cmod, cpl_c = 50, cpl_u = 0.001)
  nutri <- subset(findExchReact(cmod), exchange_type=="Community" & grepl("^EX_",react_id))
  nutri <- merge(nutri, mgx_spl)
  fwrite (nutri, file = paste0("results/Nutrition_0.01/", mgx_spl,".tsv"), sep = "\t")
}

# 1.2 extract minimum lower bound for each reaction

ExchR <- list.files(path = "results/Nutrition_0.01/", pattern = "\\.tsv", full.names=TRUE)
all_ExchR <- as.data.table(rbindlist(lapply(ExchR, fread), fill=TRUE))
fwrite(all_ExchR, file = paste0("results/Nutrition_0.01/exchangeReactions.tsv"), sep = "\t")

minlb <- all_ExchR[, .SD[which.min(lb)], by = react_id]
minlb <- minlb[, .(react_id = gsub("^EX_|_e0|_c0$", "", react_id), met_name, lb)]
fwrite(minlb, file = paste0("results/Nutrition_0.01/minLowerBounds.tsv"), sep ="\t")

##########################################################################
##
##  2. Simulation
##
##########################################################################

minlb <- fread("results/Nutrition_0.01/minLowerBounds.tsv")
minlb <- as.data.table(minlb)

for (i in seq_along(meta$mg.id)) {
  mgx_spl <- meta$mg.id[i]
  rel_models_idx <- which(mic[mgx_spl,] >= abun_co)
  rel_models_dt  <- data.table(model = colnames(mic)[rel_models_idx],
                               model_idx = rel_models_idx,
                               rel_abun = unlist(mic[mgx_spl, rel_models_idx]))
  rel_models_dt[, rel_abun := rel_abun / sum(rel_abun)]
  cmod <- joinModels(models[rel_models_dt$model], mergeLB = "min",
                     abun = rel_models_dt$rel_abun)
  cmod <- fixBMRatios(cmod)
  cmod <- fluxBMCoupling(cmod, cpl_c = 50, cpl_u = 0.001)
  # set nutritional environment
  for (j in seq_len(nrow(minlb))) {
    cpd_id <- paste0("EX_", minlb$react_id[j], "_e0")
    if (cpd_id %in% cmod@react_id) {
      cmod <- changeBounds(cmod, react = cpd_id, lb = minlb$lb[j])
    } else {
      warning("reaction not found: ", cpd_id)
    }
  }
  # increase cost of trp uptake and weight of flux minimization
  bwcosts <- rep(1, react_num(cmod))
  bwcosts[grep("^M[0-9]+_EX_cpd00065", cmod@react_id)] <- 400
  sol <- pfbaHeuristic(cmod, costcoefbw = bwcosts, pFBAcoeff = 1e-3)
  dtsol <- data.table(growth = sol@obj,
                      flux = sol@fluxes,
                      react = cmod@react_id,
                      react_name = cmod@react_name)
  dtsol[, index := str_extract(string = react, pattern = "^M[0-9]+")]
  dtsol[, mg.id := mgx_spl]
  dtsol <- merge(dtsol, cmod@community[, c("index","name")], all.x = TRUE)
  dtsol <- merge(dtsol, meta[mg.id == mgx_spl, .(mg.id, disease)], by = "mg.id", all.x = TRUE)
  dtsol <- dtsol[grepl("bio1|cpd00359|cpd27289|cpd00054|cpd00211|cpd00029|cpd00141|rxn01964|rxn00474|rxn00473|rxn08751", react)]
  dtsol <- dtsol[abs(flux) > 1e-6]
  fwrite (dtsol, file = paste0("results/Icf_Im_pM_abun0.01/scen/", mgx_spl, ".tsv"), sep = "\t")
}


# merge into one tsv
files <- list.files(path = "results/Icf_Im_pM_abun0.01/scen/", pattern = "\\.tsv", full.names=TRUE)
all_data <- rbindlist(lapply(files, fread), fill=TRUE)
fwrite (all_data, file = paste0("results/Icf_Im_pM_abun0.01/scen/all_samples_combined.tsv"), sep = "\t")

