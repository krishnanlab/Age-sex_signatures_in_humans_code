# Age-sex_signatures_in_humans_code
This GitHub repository contains the code to reproduce the results from our manuscript "**Human pan-body age- and sex-specific molecular phenomena inferred from public transcriptome data using machine learning**", which can be found [here](https://www.biorxiv.org/content/10.1101/2023.01.12.523796v2).

The accompanying web app, which visualizes the gene signatures and gene enrichment results, can be found [here](http://mlgenesignatures.org/).

# Data
The expression data used in this paper can be downloaded from zenodo [here](https://zenodo.org/records/10056218) or, if the repository is cloned, by running the bash script `data/expression/download_expression_data.sh`. Running this bash script, or downloading the data into this directory, is required to run the age or sex prediction scripts. However, the final results, performance, and weights for the models used in the paper can be found in the results directory. The gene enrichment results are stored in `zscore_gene_enrichment_analysis/results/`.
