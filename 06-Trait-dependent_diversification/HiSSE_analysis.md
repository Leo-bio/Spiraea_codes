## HiSSE analysis

### 1. Configure parameters

Before running the script [`HiSSE_analysis.R`](https://github.com/Leo-bio/Spiraea_codes/edit/main/06-Trait-dependent_diversification/HiSSE_analysis.R), you need to configure the following parameters which are written in the beginning of the script.
  - **`tree_file`**: Phylogenetic tree file path (ultramatric tree)
  - **`trait_csv_file`**:  Trait CSV file path (First column: Species names; Subsequent columns: Trait values (binary 0/1))
  - **`sampled_species`**: The number of taxon included in your phylogenetic tree.
  - **`total_species`**: The number of accepted species within your target taxa around the world.
  - **`output_file`**: The path to the output model fitting result (tsv format).
```
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
```

### 2. Run the script
```
cd <path to Spiraea_codes>
Rscript HiSSE_analysis.R
```

### 3. Check the output
The output tsv file documents the AIC， AICc, ΔAIC, and wAIC value of five models (full BiSSE, null BiSSE, full HiSSE, CID-2 and CID-4) among all binary traits.




