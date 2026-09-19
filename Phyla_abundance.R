library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dunn.test)


hrgm <- fread("data/HRGMv2_metadata_clusters.tsv")
mic <- as.data.table(read.table("data/mic.table"), keep.rownames = "mg.id")
meta <- fread("data/meta.tsv")


mic_long <- reshape2::melt(mic, id.vars = "mg.id", variable.name = "HRGMv2 Cluster", value.name = "Abundance")
mic_long <- merge(mic_long, hrgm[`HRGMv2 Cluster`== `HRGMv2 Cluster`, .(`HRGMv2 Cluster`, Phylum)],
                  by = "HRGMv2 Cluster", all.x = TRUE)
mic_long <- as.data.table(mic_long)

phyla_per_sample <- mic_long[, .(Abundance = sum(Abundance)), by = .(mg.id, Phylum)]
phyla_per_sample <- merge(phyla_per_sample, meta[mg.id == mg.id, .(mg.id, disease)], by = "mg.id", all.x = TRUE)

fwrite(phyla_per_sample, file = "output/tables/phyla_abun_per_sample.tsv", sep = "\t")

