### R code for obtaining GAPIT input data

Before running this code, make sure you have obtained the following input files:

* marker information columns from Python (script - create_marker_information.py):
  * marker_information_ID.csv.gz
  * marker_information_CHR.csv.gz
  * marker_information_POS.csv.gz
* chunks of SNP table from Python (script - create_SNP_chunks.py)
  * 100 chunks of SNP_table_100_[...].csv.gz
* csv of common accession names from Python - acc_excerpt.csv (script - extract_common_accs.py)
* csv of target gene names from Python - gene_output_for_gapit.csv (script - extract_GO_genes.py)
