# ============================================================================
# HiSSE Model Analysis Script
# ============================================================================
# This script tests whether binary traits (0/1) influence species diversification rates.
# Specifically, it estimates:
#   - Speciation rates under different trait states
#   - Extinction rates under different trait states
#   - Transition rates between trait states
#
# Usage:
#   Input a time-calibrated phylogeny (ultrametric tree) and a CSV file containing
#   species traits, then construct BiSSE likelihood functions and estimate parameters
#   via maximum likelihood.
#
# Dependencies:
#   ape           - Phylogenetic tree reading and processing
#   diversitree   - Implementation of BiSSE/MuSSE/QuaSSE models
#   phytools      - Phylogenetic analysis utilities
#   hisse         - Hidden state speciation and extinction models
#
# Core functions:
#   make.bisse()  - Construct BiSSE likelihood function
#   hisse()       - Fit HiSSE/BiSSE models
# ============================================================================

# ----------------------------------------------------------------------------
# 0. Load required packages
# ----------------------------------------------------------------------------

# rm(list = ls())  # Clear all objects from the R environment

library(ape)          # Phylogenetic tree processing
library(diversitree)  # SSE model implementation
library(phytools)     # Phylogenetic analysis utilities
library(AICcmodavg)   # AICc calculations
library(hisse)        # Hidden state SSE models

# ----------------------------------------------------------------------------
# 1. Configure parameters
# ----------------------------------------------------------------------------

# (1) Phylogenetic tree file path
tree_file <- "./06-Trait-dependent_diversification/Spiraea_no_outgroup.dated.ultrametric.tre"

# (2) Trait CSV file path
trait_csv_file <- "./06-Trait-dependent_diversification/Spiraea_traits.csv"

# CSV file format:
#   First column: Species names
#   Subsequent columns: Trait values (binary 0/1)

# Example:
# Species,Inflorescence_type,Follicle_inflation,...
# Spiraea_bella,0,0,...
# Spiraea_fritschiana,0,0,...

# (3) Sampling proportions
sampled_species <- 56
total_species <- 110
sampling_proportion <- sampled_species / total_species
sampling.f <- c(sampling_proportion, sampling_proportion)

# (4) Output tsv file
output_file <- "./06-Trait-dependent_diversification/HiSSE_analysis_output.tsv"

# ----------------------------------------------------------------------------
# 2. Read tree and trait data
# ----------------------------------------------------------------------------

# Read phylogenetic tree
tree <- read.tree(tree_file)

# Read trait data
traits <- read.csv(trait_csv_file, header = TRUE)

# Inspect first few rows to verify correct reading
head(traits)

# ----------------------------------------------------------------------------
# Check if tree is ultrametric
#
# BiSSE models require an ultrametric tree (all tips equidistant from root),
# i.e., a time-calibrated phylogeny.
# ----------------------------------------------------------------------------
cat("Checking ultrametric tree...\n")

if (!is.ultrametric(tree)) {
  cat("Tree is NOT ultrametric — converting using force.ultrametric\n")
  tree <- force.ultrametric(tree, method = "extend")
  # Alternative: use chronos() for time calibration
  # tree <- chronos(tree, quiet = TRUE)
} else {
  cat("Tree is ultrametric ✓\n")
}

# ----------------------------------------------------------------------------
# Get all trait columns (excluding the first column which is Species)
# ----------------------------------------------------------------------------
trait_columns <- colnames(traits)[-1]  # Exclude first column (Species)
cat("Traits to analyze:", paste(trait_columns, collapse = ", "), "\n")

# ----------------------------------------------------------------------------
# 3. Define AIC weight calculation function
# ----------------------------------------------------------------------------
calc_aic_weights <- function(aic_vec) {
  delta_aic <- aic_vec - min(aic_vec)
  weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
  return(list(delta_aic = delta_aic, weights = weights))
}

# ----------------------------------------------------------------------------
# 4. Build transition rate matrices
# ----------------------------------------------------------------------------
trans.rates.bisse <- TransMatMaker.old(hidden.states = FALSE)
trans.rates.bisse <- ParEqual(trans.rates.bisse, c(1, 2))
trans.rates.hisse <- TransMatMaker.old(hidden.states = TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, c(3, 5, 8, 10))
trans.rates.hisse <- ParEqual(trans.rates.hisse, c(1, 2, 1, 3, 1, 4, 1, 5, 1, 6, 1, 7, 1, 8))

# ----------------------------------------------------------------------------
# 5. Process each trait and fit models
# ----------------------------------------------------------------------------

# Initialize results data frame
all_results <- data.frame(
  Trait = character(),
  Model = character(),
  lnL = numeric(),
  AIC = numeric(),
  AICc = numeric(),
  Delta_AIC = numeric(),
  AIC_weight = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each trait
for (trait_name in trait_columns) {
  cat("\n========================================\n")
  cat("Processing trait:", trait_name, "\n")
  cat("========================================\n")
  
  tip_states <- traits[[trait_name]]
  names(tip_states) <- traits$Species
  
  common <- intersect(tree$tip.label, names(tip_states))
  tree_bisse <- keep.tip(tree, common)
  tip_states <- tip_states[common]
  
  unique_values <- unique(tip_states)
  if (!all(unique_values %in% c(0, 1))) {
    cat("WARNING: Trait", trait_name, "is not binary. Skipping.\n")
    next
  }
  
  # -------------------------------------------------------------
  # HiSSE input
  # -------------------------------------------------------------
  hisse_data <- data.frame(
    sp = names(tip_states),
    state = as.numeric(tip_states),
    stringsAsFactors = FALSE
  )
  
  # Model fitting
  cat("\nFitting BiSSE_full model...\n")
  hisse.BiSSE_full.fit <- hisse(phy = tree_bisse, data = hisse_data, f = sampling.f, 
                                hidden.states = FALSE, turnover = c(1, 2), 
                                eps = c(1, 1), trans.rate = trans.rates.bisse)
  
  cat("Fitting BiSSE_null model...\n")
  hisse.BiSSE_null.fit <- hisse(phy = tree_bisse, data = hisse_data, f = sampling.f, 
                                hidden.states = FALSE, turnover = c(1, 1), 
                                eps = c(1, 1), trans.rate = trans.rates.bisse)
  
  cat("Fitting HiSSE_full model...\n")
  hisse.HiSSE_full.fit <- hisse(phy = tree_bisse, data = hisse_data, f = sampling.f, 
                                hidden.states = TRUE, turnover = c(1, 2, 3, 4), 
                                eps = c(1, 1, 1, 1), trans.rate = trans.rates.hisse)
  
  cat("Fitting CID_2 model...\n")
  hisse.CID_2.fit <- hisse(phy = tree_bisse, data = hisse_data, f = sampling.f, 
                           hidden.states = TRUE, turnover = c(1, 1, 2, 2), 
                           eps = c(1, 1, 1, 1), trans.rate = trans.rates.hisse)
  
  cat("Fitting CID_4 model...\n")
  hisse.CID_4.fit <- hisse.null4.old(phy = tree_bisse, data = hisse_data, f = sampling.f, 
                                     eps.anc = rep(1, 8))
  
  model_list <- list(
    BiSSE_full = hisse.BiSSE_full.fit,
    BiSSE_null = hisse.BiSSE_null.fit,
    HiSSE_full = hisse.HiSSE_full.fit,
    CID_2      = hisse.CID_2.fit,
    CID_4      = hisse.CID_4.fit
  )
  
  # Extract AICc and calculate Delta AIC
  aicc_values <- sapply(model_list, function(x) x$AICc)
  aic_results <- calc_aic_weights(aicc_values)
  
  trait_results <- data.frame(
    Trait       = trait_name,
    Model       = names(model_list),
    lnL         = sapply(model_list, function(x) x$loglik),
    AIC         = sapply(model_list, function(x) x$AIC),
    AICc        = aicc_values,
    Delta_AICc  = aic_results$delta_aic,
    AICc_weight = aic_results$weights,
    stringsAsFactors = FALSE
  )
  
  all_results <- rbind(all_results, trait_results)
  print(trait_results)
}

# ----------------------------------------------------------------------------
# 6. Save results to TSV file
# ----------------------------------------------------------------------------
write.table(all_results, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n========================================\n")
cat("All results saved to:", output_file, "\n")
cat("========================================\n")
print(all_results)

