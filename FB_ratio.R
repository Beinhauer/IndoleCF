library(data.table)
library(stringr)
library(dplyr)
library(tidyr)

hrgm <- fread("data/HRGMv2_metadata_clusters.tsv")
mic <- read.table("data/mic.table")
meta <- fread("data/meta.tsv")

mic <- cbind(rownames(mic), data.frame(mic, row.names=NULL))
names(mic)[1] <- "mg.id"

mic <- merge(mic, meta[mg.id == mg.id, .(mg.id, disease)], by = "mg.id",
             all.x = TRUE)

rownames(mic) <- NULL
mic_t <- mic |>
  tibble::column_to_rownames(var="mg.id") |>
  t() |>
  as.data.frame() |>
  mutate(across(everything(), as.numeric))

mic_phyla <- cbind(rownames(mic_t), data.frame(mic_t, row.names=NULL))
names(mic_phyla)[1] <- "HRGMv2 Cluster"
mic_phyla <- left_join(mic_phyla, hrgm[, c("HRGMv2 Cluster", "Phylum")],
                       by = "HRGMv2 Cluster")

sample_disease <- mic[, c("mg.id", "disease")]


## Summarize all Firmicutes types

mic_phyla <- mic_phyla %>%
  mutate(Phylum = recode(Phylum, "p__Firmicutes_A" = "p__Firmicutes",
                         "p__Firmicutes_B" = "p__Firmicutes",
                         "p__Firmicutes_C" = "p__Firmicutes"))

phyla_sums <- mic_phyla |>
  select(-`HRGMv2 Cluster`) |>
  group_by(Phylum) |>
  summarise(across(everything(), sum)) |>
  pivot_longer(-Phylum, names_to = "mg.id", values_to = "Abundance")|>
  left_join(sample_disease, by = "mg.id")


##########################################################################
##
##  Calculate Firmicutes/Bacteroides ratio
##
##########################################################################

setDT(phyla_sums)
fb_ratio <- phyla_sums[Phylum %in% c("p__Firmicutes", "p__Bacteroidota"),
                       .(total = sum(Abundance)), by = .(mg.id, disease, Phylum)]
fb_ratio <- dcast(fb_ratio, mg.id + disease ~ Phylum, value.var="total", fill =0)
fb_ratio[, FB_ratio := p__Firmicutes / p__Bacteroidota]
fb_ratio[is.infinite(FB_ratio) | is.nan(FB_ratio), FB_ratio := NA]

fwrite(fb_ratio, file = "output/tables/FB_ratio.tsv", sep = "\t")

