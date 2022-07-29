# Lab rotation project at UZH

### Code Workflow

* [Python] Open GSEA and AraGWAS data
  * Identify and extract common accessions
  * Remove top 5% and bottom 5% variance SNPs (top 5% - most accessions' genome varies greatly from the reference, bottom 5% - most accessions' genome doesn't vary much from the reference)
  * Extract eligible SNPs for common accessions (SNP table) and the genomic positions of each SNP (marker information table)
* [Python] Select genes based on GO categories present in the GSEA gene list
  * Run the gene list through [ShinyGO](http://bioinformatics.sdstate.edu/go/) to obtain GO categories that are covered by the GSEA gene list
  * Download the GO category summary with gene lists for each category
  * Select categories with less than 250 genes, representing specific functions in A. thaliana
  * Focus on categories related to plant immune response and photoperiodism
  * Export gene lists for GAPIT
* [R] Prepare data for GAPIT
  * Trim GSEA table to only include common accessions and genes of interest as phenotype in GAPIT
  * Transform read counts in GSEA table into log2(n + 0.01) to decrease the spread of phenotype data
  * Assemble marker information from 3 columns (ID, Chr, Pos)
  * Assemble SNP table from chunks (0-99 / init + 1-99)
  * Make 1 GAPIT test run to obtain kinship table
* [R] Run GAPIT
  * GAPIT run on 528+123 gene log2 read counts for 681 accessions
  * Runtime around 7-8 days (~15-20 minutes per gene)
* [Python] Analyze coexpression patterns from GSEA and GAPIT prediction data
  * Focus on one GO term at a time (based on a specific gene list)
  * Create coexpression and mutual rank tables from GSEA data or assembled GAPIT prediction data
  * Shortlist top 10/20 mutual ranks or correlation values above 0.5 and turn into binary matrices
  * Investigate top connections in shortlisted matrices as well as in matrices of common connections in GSEA and GAPIT data
