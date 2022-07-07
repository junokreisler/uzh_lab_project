marker_info <- readRDS('marker_information_for_GAPIT.rds')
GWAS_data <- readRDS('GWAS_data_for_GAPIT.rds')
SNP_table <- readRDS('SNP_Table_full_AraGWAS.rds')
print('SNP table loaded.')
library(GAPIT3)
print('Reading gene shortlist for analysis')
gene_shortlist <- na.omit(t(read.csv2('gene_output_for_gapit.csv', 
                                      sep = ',', header = FALSE,)))
print('GAPIT loaded, loading first prediction results')
Pred_table <- read.csv('GAPIT_output/GAPIT.MLM.AT1G01010.Pred.result.csv')$Pred
print('Initial prediction column loaded. Starting BLUP loop...')
Kinship_matrix <- read.csv('GAPIT_output/GAPIT.Kin.Zhang.csv', header = FALSE)
for (i in gene_shortlist) {
  t_start <- Sys.time()
  BLUP_test <- GAPIT(
    Y=GWAS_data[,c('acc_id',i)], #[,c(acc_id, gene)]
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
Pred_table <- Pred_table[,-1]
saveRDS(Pred_table, 'Pred_table_shortlist_cut.rds')
write.csv(Pred_table,'data/Pred_table.csv')
print('saved prediciton table for all shortlisted genes')
