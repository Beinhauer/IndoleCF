# ANALYSIS OF INDOLE SALVAGERS:
# Indole salvager identification through growth simulations with adjusted lower
# bounds of tryptophan and indole

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

abun_co <- 0.005

# create string of clusters above detection abundance in at least one sample
cluster <- colnames(mic)[apply(mic, 2, function(x) any(x >= abun_co))]

# -------- 1. Simulate growth with Trp --------#

results_wT <- data.table()

for (i in cluster){
  model <- models[[i]]
  sol <- pfba(model)
  dtflux <- data.table(cluster = i,
                       growth = sol@obj)
  results_wT <- rbind(results_wT, dtflux)
}

# -------- 2. Simulate growth withot Trp --------#

results_nT <- data.table()

for (i in cluster){
  model <- models[[i]]
  if("EX_cpd00065_e0" %in% model@react_id){
    model <- changeBounds(model, react = "EX_cpd00065_e0", lb=0)
  }
  sol <- pfba(model)
  dtflux <- data.table(cluster = i,
                       growth = sol@obj)
  results_nT <- rbind(results_nT, dtflux)
}

# -------- 3. Simulate growth of auxotrophs with indole --------#

# Trp auxotrophs
auxcluster <- results_nT$cluster[which(results_nT$growth == 0)]

results_IS <- data.table()

for (i in auxcluster){
  model <- models[[i]]
  if("EX_cpd00065_e0" %in% model@react_id){
    model <- changeBounds(model, react = "EX_cpd00065_e0", lb=0)
  }
  if("EX_cpd00359_e0" %in% model@react_id){
    model <- changeBounds(model, react = "EX_cpd00359_e0", lb=-10)
  }
  sol <- pfba(model)
  dtflux <- data.table(cluster = i,
                       growth = sol@obj)
  results_IS <- rbind(results_IS, dtflux)
}

results_IS <- merge(results_IS,
                    hrgm[, .(`HRGMv2 Cluster`, Phylum, Genus, Species)],
                    by.x = "cluster", by.y = "HRGMv2 Cluster")

# subset  indole salvagers
results_onlyIS <- subset(results_IS, growth > 0)
results_onlyIS <- cbind(results_onlyIS, X = "x") #was macht das?

rel_salvager <- intersect(results_onlyIS$cluster, colnames(mic))


# calculate relative abundance of salvagers in samples

salvager_abun_long <- melt(
  mic[, c("mg.id", rel_salvager), with = FALSE],
  id.vars = "mg.id",
  variable.name = "cluster",
  value.name = "abundance"
)

salvager_abun_long <- merge(salvager_abun_long,
                            hrgm[, .(`HRGMv2 Cluster`, Phylum, Genus, Species)],
                            by.x = "cluster", by.y = "HRGMv2 Cluster")
salvager_abun_long <- merge(salvager_abun_long, meta[, .(mg.id, disease)],
                            by = "mg.id")
salvager_abun_long <- as.data.table(salvager_abun_long)

dt_salvager_abun <- salvager_abun_long[, .(Abundance = sum(abundance),
                                           Mean = fifelse(
                                             sum(abundance != 0) == 0, 0,
                                             mean(abundance[abundance != 0])),
                                           n_01 = sum(abundance >= 0.01),
                                           n_005 = sum(abundance >= 0.005)),
                                       by = .(mg.id, disease)]
dt_salvager_abun[, type := "Salvager"]

fwrite(dt_salvager_abun, file = "output/tables/salvager_abundance.tsv",
       sep ="\t")


# calculate prevalence of salvagers in samples

salvager_prev_all <- salvager_abun_long[, .(prevalent = mean(abundance > 0)*100),
                                        by = .(cluster, mg.id)]
salvager_prev_all <- merge(salvager_prev_all, meta[, .(mg.id, disease)],
                           by = "mg.id")
fwrite(salvager_prev_all, file = "output/tables/salvager_prevalence_all.tsv",
       sep = "\t")


# calculate prevalence of salvagers per disease

salvager_prev <- salvager_prev_all[, .(prevalence = mean(prevalent > 0)*100),
                                   by = .(cluster, disease)]
salvager_prev <- merge(salvager_prev,
                       hrgm[, .(`HRGMv2 Cluster`, Phylum, Genus, Species)],
                       by.x = "cluster", by.y = "HRGMv2 Cluster")
salvager_prev[, Species := gsub("s__", "", Species)]
salvager_prev[, Phylum := gsub("p__", "", Phylum)]

fwrite(salvager_prev, file = "output/tables/salvager_prevalence.tsv", sep = "\t")
