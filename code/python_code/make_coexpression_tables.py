import numpy as np
import pandas as pd
import os
import matrix_network_functions as func
from sklearn.metrics import jaccard_score as js


# 1)   from an input file, read the list of genes of interest
# 1.1) input file format - comma-separated word list

# ONLY CHANGE THE FOLLOWING TWO
gene_list_path = 'flowering_list_full.txt'
output_folder = 'flowering_output'

try:
    os.mkdir(output_folder)
except:
    print('Output folder', output_folder, 'already exists')

gene_list = func.import_gene_list(gene_list_path)

acc_list = func.load_accs()

# 2)   based on the genes of interest,
# 2.1) extract read count data from GSEA data

GSEA_filename = 'GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv.gz'
GSEA_data = pd.read_csv(GSEA_filename, delimiter='\t')

target_genes = GSEA_data[GSEA_data['gene_id'].isin(gene_list)]
select_columns = ['gene_id'] + acc_list

target_matrix = target_genes[select_columns].set_index('gene_id')
target_matrix.to_csv(output_folder+'/1_GSEA_readcount_'+gene_list_path+'.csv')

target_log2_matrix = np.log2(target_matrix+0.01)
target_log2_matrix.to_csv(output_folder+'/1_GSEA_readcount_log2_'+gene_list_path+'.csv')

# 2.2) extract the Pred values from gene files and put them into one table
pred_filenames = []
gene_list_sorted = list(target_matrix.index)

for gene in gene_list_sorted:
    pred_filenames.append('GAPIT.MLM.'+gene+'.Pred.result.csv')

pred_table = assemble_prediction_table(gene_list_sorted)
pred_table.to_csv(output_folder+'/2_GAPIT_prediction_log2_'+gene_list_path+'.csv')

# 3)   export
# 3.1) coexpression tables

target_log2_coexpression = target_log2_matrix.T.astype(float).corr(method = 'pearson')
target_log2_coexpression.to_csv(output_folder+'/3_GSEA_correlation_'+gene_list_path+'.csv')

pred_coexpression = pred_table.astype(float).corr(method = 'pearson')
pred_coexpression.to_csv(output_folder+'/3_GAPIT_correlation_'+gene_list_path+'.csv')

# 3.2) mutual rank tables

target_log2_mr = len(gene_list) - func.mutual_rank(target_log2_coexpression.rank(axis=0), target_log2_coexpression.rank(axis=1))
target_log2_mr.to_csv(output_folder+'/4_GSEA_mutual_rank_'+gene_list_path+'.csv')

pred_log2_mr = len(gene_list) - func.mutual_rank(pred_coexpression.rank(axis=0), pred_coexpression.rank(axis=1))
pred_log2_mr.to_csv(output_folder+'/4_GAPIT_mutual_rank_'+gene_list_path+'.csv')

# export binary tables and common connection table

target_corr_bin = func.shortlist_to_binary_matrix(target_log2_coexpression, 0.5, binary = True, abs_corr = False)
target_mr_bin = func.shortlist_to_binary_matrix(target_log2_mr, 10, binary = True)
target_corr_bin.to_csv(output_folder+'/5_GSEA_binary_shortlist_corr_'+gene_list_path+'.csv')
target_mr_bin.to_csv(output_folder+'/5_GSEA_binary_shortlist_mr_'+gene_list_path+'.csv')


pred_corr_bin = func.shortlist_to_binary_matrix(pred_coexpression, 0.5, binary = True, abs_corr = False)
pred_mr_bin = func.shortlist_to_binary_matrix(pred_log2_mr, 10, binary = True)
pred_corr_bin.to_csv(output_folder+'/5_GAPIT_binary_shortlist_corr_'+gene_list_path+'.csv')
pred_mr_bin.to_csv(output_folder+'/5_GAPIT_binary_shortlist_mr_'+gene_list_path+'.csv')

common_connections_corr = target_corr_bin * pred_corr_bin
common_connections_mr = target_mr_bin * pred_mr_bin
common_connections_corr.to_csv(output_folder+'/5_COMMON_binary_shortlist_corr_'+gene_list_path+'.csv')
common_connections_mr.to_csv(output_folder+'/5_COMMON_binary_shortlist_mr_'+gene_list_path+'.csv')


similarity = js(np.array(target_corr_bin), np.array(pred_corr_bin), average = 'samples')
