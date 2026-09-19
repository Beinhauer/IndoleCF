# ANALYSIS OF INDOLE PRODUCER ABUNDANCES:
# Indole producer identified by presence of synthesis pathways (last enzyme of
# pathway)
# Additional calculation of indole salvager to individual indole derivative
# producer abundance

library(qs2)
library(cobrar)
library(cobrarCPLEX)
library(data.table)
library(stringr)
library(vegan)
library(dunn.test)
library(rstatix)

COBRAR_SETTINGS("SOLVER","cplex")

hrgm <- fread("data/HRGMv2_metadata_clusters.tsv")

models <- qs_read("data/models.qs2")
mic <- as.data.table(read.table("data/mic.table"), keep.rownames = "mg.id")
meta <- fread("data/meta.tsv")

is <- fread("output/tables/salvager_abundance.tsv")

#############################################################################
##
## Indole derivative producer abundances
##
#############################################################################

# 1. IPA producer

ipa_pattern <- "rxn45361"
ipa_producer <- names(models)[sapply(models,
                                     function(m) any(grepl(ipa_pattern,
                                                           m@react_id)))]

ipa_prod_abun <- as.data.table(melt(
  mic[, c("mg.id", ipa_producer), with = FALSE],
  id.vars = "mg.id",
  variable.name = "cluster",
  value.name = "abundance"
))

ipa_prod_abun <- merge(ipa_prod_abun, meta[, .(mg.id, disease)],
                       by = "mg.id", all.x = TRUE)

ipa_dt_prod_abun <- ipa_prod_abun[, .(Abundance = sum(abundance),
                                      Mean = mean(abundance)),
                                  by = .(mg.id, disease)]
ipa_dt_prod_abun[, type := "IPA Producer"]

fwrite(ipa_dt_prod_abun, file = "output/tables/ipa_producer_abundance.tsv",
       sep ="\t")

# 1.1 indole salvager/ipa producer ratio
dt_ipa_is <- cbind(is[, .(mg.id, disease, salvager = Abundance)],
                   ipa_dt_prod_abun[, .("ipa producer" = Abundance)])
dt_ipa_is[, ratio := (dt_ipa_is$salvager/dt_ipa_is$`ipa producer`)]
dt_ipa_is[, type := "Indole salvager/ipa producer"]

# 2. ILA producer

ila_pattern <- "rxn01447"

ila_producer <- names(models)[sapply(models,
                                     function(m) any(grepl(ila_pattern,
                                                           m@react_id)))]

ila_prod_abun <- as.data.table(melt(
  mic[, c("mg.id", ila_producer), with = FALSE],
  id.vars = "mg.id",
  variable.name = "cluster",
  value.name = "abundance"
))

ila_prod_abun <- merge(ila_prod_abun, meta[, .(mg.id, disease)],
                       by = "mg.id", all.x = TRUE)

ila_dt_prod_abun <- ila_prod_abun[, .(Abundance = sum(abundance),
                                      Mean = mean(abundance)),
                                  by = .(mg.id, disease)]
ila_dt_prod_abun[, type := "ILA Producer"]

fwrite(ila_dt_prod_abun, file = "output/tables/ila_producer_abundance.tsv",
       sep ="\t")

# 2.1 indole salvager/ila producer ratio

dt_ila_is <- cbind(is[, .(mg.id, disease, salvager = Abundance)],
                   ila_dt_prod_abun[, .("ila producer" = Abundance)])
dt_ila_is[, ratio := (dt_ila_is$salvager/dt_ila_is$`ila producer`)]
dt_ila_is[, type := "Indole salvager/ila producer"]

# 3. IAA producer

iaaH_pattern <- "rxn02222" # only production
iaadH_pattern <- "rxn01937"

iaaH_producer <- names(models)[sapply(models,
                                      function(m) any(grepl(iaaH_pattern,
                                                            m@react_id)))]
iaadH_producer <- names(models)[sapply(models,
                                       function(m) any(grepl(iaadH_pattern,
                                                             m@react_id)))]
iaa_producer <- unique(c(iaadH_producer, iaaH_producer))

iaa_prod_abun <- as.data.table(melt(
  mic[, c("mg.id", iaa_producer), with = FALSE],
  id.vars = "mg.id",
  variable.name = "cluster",
  value.name = "abundance"
))

iaa_prod_abun <- merge(iaa_prod_abun, meta[, .(mg.id, disease)],
                       by = "mg.id", all.x = TRUE)

iaa_dt_prod_abun <- iaa_prod_abun[, .(Abundance = sum(abundance),
                                      Mean = mean(abundance)),
                                  by = .(mg.id, disease)]
iaa_dt_prod_abun[, type := "IAA Producer"]

fwrite(iaa_dt_prod_abun, file = "output/tables/iaa_producer_abundance.tsv",
       sep ="\t")

# 2.1 indole salvager/iaa producer ratio

dt_iaa_is <- cbind(is[, .(mg.id, disease, salvager = Abundance)],
                   iaa_dt_prod_abun[, .("iaa producer" = Abundance)])
dt_iaa_is[, ratio := (dt_iaa_is$salvager/dt_iaa_is$`iaa producer`)]
dt_iaa_is[, type := "Indole salvager/iaa producer"]

dt <- rbind(dt_ipa_is[, .(mg.id, disease, ratio, type)],
            dt_ila_is[, .(mg.id, disease, ratio, type)],
            dt_iaa_is[, .(mg.id, disease, ratio, type)])

fwrite(dt, file = "output/tables/ratio_salvager_id_producer.tsv", sep ="\t")
