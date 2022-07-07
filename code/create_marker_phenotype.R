#################################################
# Create phenotype info for GAPIT
#################################################

print('Reading GSEA data')

GSEA_data <- read.table('GSE80744_normCounts.tsv', header = TRUE)
gene_ids <- GSEA_data['gene_id']
acc_names <- t(read.csv('acc_excerpt.csv')) # common accessions between GSEA and AraGWAS
GSEA_data <- as.data.frame(t(GSEA_data[,acc_names]))
colnames(GSEA_data) <- t(gene_ids)
GSEA_data <- cbind(rownames(GSEA_data), GSEA_data)
colnames(GSEA_data)[1] <- 'acc_id'

print('Phenotype data prepared')
saveRDS(GSEA_data, 'GSEA_data_for_GAPIT.rds')

#################################################
# Load the necessary marker data from Python output
#################################################
# make marker information table
print('Reading marker information')

marker_ID <- read.table("marker_information/marker_information_ID.csv.gz",
                   header = FALSE)
marker_CHR <- read.table("marker_information/marker_information_CHR.csv.gz", 
                    header = FALSE)
marker_POS <- read.table("marker_information/marker_information_POS.csv.gz", 
                    header = FALSE)

marker_info <- cbind(marker_ID, marker_CHR)
marker_info <- cbind(marker_info, marker_POS)
colnames(marker_info) <- c('ID','CHR','POS')

print('Marker information generated')
saveRDS(marker_info, 'marker_information_for_GAPIT.rds')
