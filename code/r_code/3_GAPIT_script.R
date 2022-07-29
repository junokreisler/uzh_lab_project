marker_info <- readRDS('marker_information_for_GAPIT.rds')
GSEA_data <- readRDS('GSEA_data_for_GAPIT_log2.rds')
SNP_table <- readRDS('SNP_Table_full_AraGWAS.rds')
print('SNP table loaded.')

library(GAPIT3) # running in a miniconda environment for R

print('Reading gene shortlist for analysis')
gene_shortlist <- na.omit(t(read.csv2('gene_output_for_gapit.csv', 
                                      sep = ',', header = FALSE,)))
print('GAPIT loaded, loading first prediction results')

Pred_table <- read.csv('GAPIT_output/GAPIT.MLM.AT1G01010.Pred.result.csv')$Pred

# the initial prediction column was obtained from the GAPIT run on the very first gene in the phenotype table.
# it acts as an example output, as well as a dummy initial file, onto which further results are attached during the loop.
# the column of the dummy file is deleted at line 40 because AT1G01010 is not in the target gene list.

print('Initial prediction column loaded. Starting BLUP loop...')

Kinship_matrix <- read.csv('GAPIT_output/GAPIT.Kin.Zhang.csv', header = FALSE) # prevents kinship recalculation at every run

for (i in gene_shortlist) {
  t_start <- Sys.time()
  BLUP_test <- GAPIT(
    Y=GSEA_data[,c('acc_id',i)], #[,c(acc_id, gene_name)]
    GD=SNP_table,
    GM=marker_info,
    KI=Kinship_matrix,
    # SNP.MAF=0.05,
    model="gBLUP")
  Pred_table <- cbind(Pred_table, BLUP_test$Pred$Prediction)
  print('one BLUP done.')
  print(Sys.time()-t_start)
  print('-----------------------------------------------------')
  print(i)
  print('DONE')
  print('-----------------------------------------------------')
}
Pred_table <- Pred_table[,-1] # removing the dummy / init column
saveRDS(Pred_table, 'Pred_table_shortlist_cut.rds')
write.csv(Pred_table,'data/Pred_table.csv')
print('saved prediciton table for all shortlisted genes')
