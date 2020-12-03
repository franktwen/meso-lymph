# meso-lymph

This repo contains scripts for naming genes from ExomeDepth outputs.

ExomeDepth outputs are placed in `data`. 
The script `scripts/name_genes.R` gets gene names from Ensembl using the the R `biomaRt` package and merges them with the ExomeDepth output files. The named output files are written to `outputs`.
The output files contain the following new columns:

`hgnc_symbol` gene name

`start_gene`  gene starting point from ensembl

`end_gene`    gene ending point from ensembl

`start_intersect` used for merging data based on coordinates

`end_intersect` used for merging data based on coordinates

`filename` the name of the file

`rowId` a unique row indicator for each row in the original ExomeDepth file (e.g. 1-10000 for each file)

Note that unnamed reads are not included in the output.
