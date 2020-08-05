# PhageLoader

Scripts for pre-processing and calculating enrichment for phage count data between two groups.

Required Input File :
  - Read count matrix where rows are peptides and columns are samples (CSV format)

Pipeline :

```{bash}

# Preprocesses input csv for loading into enrichment
pipeline/preprocess.sh <input_csv> <output_directory>

# Run Z-Score Enrichment on a processed dataset
bin/enrichment -i <input_arma> -n <sample_names> ... -o <output_prefix>

# Can run many different options on the same preprocessed dataset with minimal overhead
for z_score in {1..15}; do
  for c_min in {1..15}; do
    bin/enrichment -i <input_arma> -n <sample_names> -z $z_score -c $c_min -o <output_prefix>
  done;
done;

```
