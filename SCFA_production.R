library(data.table)
library(dunn.test)
library(dplyr)

dt <- fread("results/ICF_0.01/all_samples_combined_scen0.tsv")

meta <- fread("data/meta.tsv")

scfa <- subset(dt,
               react == "EX_cpd00029_e0"
               | react == "EX_cpd00211_e0"
               | react == "EX_cpd00141_e0")
scfa <- scfa[, norm_flux := flux / growth ]

scfa <- scfa %>%
  mutate(react_name = recode(react_name, "Acetate-e0 exchange" = "Acetate",
                             "Butyrate-e0 exchange" = "Butyrate",
                             "Propionate-e0 exchange" = "Propionate"))
scfa$disease <- as.factor(scfa$disease)

fwrite(scfa, file = "output/tables/scfa_production.tsv", sep ="\t")
