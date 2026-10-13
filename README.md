# IndoleCF
*In silico* modeling approach to analyze differences in indole cross-feeding between healthy and dysbiotic (IBD-associated) microbiomes. Data for the analysis was acquired from the study by Franzosa et al. [^1]. Individual gapseq models were supplied by Prof. Dr. Silvio Waschina (Christian-Albrechts-Universität zu Kiel, dept. Nutriinformatics) and may be requested there.

Please note that all simulations and statistical analysis require the following R packages: qs2, cobrar, cobrarCPLEX (IBM's ILOG CPLEX solver), data.table, stringr, dplyr, tidyr, dunn.test, vegan and rstatix. 

---
## 1. Microbial community simulations
Microbial communities can be simulated using CoM.R for the defined relative bacterial abundance cut-offs (1% and 0.5%). First step in the simulations is the curation of the underlying nutritional environment based on present and absent metabolic pathways. All simulations can be repeated in six different scenarios with varying simulation parameters to allow a subsequent sensitivity analysis.
## 2. Indole cross-feeding analysis
Indole cross-feeding flux and prevalence can be evaluated via Indole_crossfeeding_analysis.R for all previously simulated scenarios.
## 3. Indole salvager simulations
Indole salvagers can be identified by simulating the individual microbial growth rates in different medias (without tryptophan & indole, without tryptophan but with indole) with the script Indole_salvager_simulations.R. 
## 4. Additional analysis
Indole derivative (IAA, ILA & IPA) producer abundance can be examined by searching for the presence of the respective (last) enzymes in the corresponding biosynthesis pathway via Indole_derivative_producer_analysis.R. To evaluate the underlying microbial compositions Phyla_abundance.R, FB_ratio.R and further analysis via statistic.R can be performed. Additionally, SCFA_production.R allows for an analysis of flux predictions for the production of acetate, butyrate and propionate.


[^1]: Franzosa EA, Sirota-Madi A, Avila-Pacheco J, Fornelos N, Haiser HJ, Reinker S, et al. Gut microbiome structure and metabolic activity in inflammatory bowel disease. Nat Microbiol. 2019 Feb;4(2):293–305. doi:10.1038/s41564-018-0306-4 PubMed PMID: 30531976; PubMed Central PMCID: PMC6342642.
