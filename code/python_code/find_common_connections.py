import pandas as pd
import numpy as np
import matrix_network_functions as func
from sklearn.metrics import jaccard_score as js

# 3.1) coexpression tables

target_log2_coexpression = pd.read_csv(output_folder+'/3_GSEA_correlation_'+gene_list_path+'.csv').set_index('gene_id')
pred_log2_coexpression = pd.read_csv(output_folder+'/3_GAPIT_correlation_'+gene_list_path+'.csv').set_index('Unnamed: 0')

# 3.2) mutual rank tables

target_log2_mr = pd.read_csv(output_folder+'/4_GSEA_mutual_rank_'+gene_list_path+'.csv').set_index('gene_id')
pred_log2_mr = pd.read_csv(output_folder+'/4_GAPIT_mutual_rank_'+gene_list_path+'.csv').set_index('Unnamed: 0')

# target correlation/mr binary matrices

target_corr_bin = func.shortlist_to_binary_matrix(target_log2_coexpression, 0.5, binary = True, abs_corr = False)
target_mr_bin = func.shortlist_to_binary_matrix(target_log2_mr, 10, binary = True)
target_corr_bin.to_csv(output_folder+'/5_GSEA_binary_shortlist_corr_'+gene_list_path+'.csv')
target_mr_bin.to_csv(output_folder+'/5_GSEA_binary_shortlist_mr_'+gene_list_path+'.csv')

# prediction correlation/mr binary matrices

pred_corr_bin = func.shortlist_to_binary_matrix(pred_log2_coexpression, 0.5, binary = True, abs_corr = False)
pred_mr_bin = func.shortlist_to_binary_matrix(pred_log2_mr, 10, binary = True)
pred_corr_bin.to_csv(output_folder+'/5_GAPIT_binary_shortlist_corr_'+gene_list_path+'.csv')
pred_mr_bin.to_csv(output_folder+'/5_GAPIT_binary_shortlist_mr_'+gene_list_path+'.csv')

# investigate how similar the top connections are

target_corr_top_nodes = func.top_binary_nodes_dict(target_corr_bin)
target_mr_top_nodes = func.top_binary_nodes_dict(target_mr_bin)
pred_corr_top_nodes = func.top_binary_nodes_dict(pred_corr_bin)
pred_mr_top_nodes = func.top_binary_nodes_dict(pred_mr_bin)

# common connections obtained by multiplying binary matrices

common_connections_corr = target_corr_bin * pred_corr_bin
common_connections_mr = target_mr_bin * pred_mr_bin
common_connections_corr.to_csv(output_folder+'/5_COMMON_binary_shortlist_corr_'+gene_list_path+'.csv')
common_connections_mr.to_csv(output_folder+'/5_COMMON_binary_shortlist_mr_'+gene_list_path+'.csv')

similarity = js(np.array(target_corr_bin), np.array(pred_corr_bin), average = 'samples')

common_connections_corr_top_nodes = func.top_binary_nodes_dict(common_connections_corr)
common_connections_mr_top_nodes = func.top_binary_nodes_dict(common_connections_mr)

list_of_dicts = [target_corr_top_nodes,target_mr_top_nodes,
                 pred_corr_top_nodes,pred_mr_top_nodes,
                 common_connections_corr_top_nodes,common_connections_mr_top_nodes]

list_of_rownames = ['target_corr_top_nodes', 'target_mr_top_nodes',
                    'pred_corr_top_nodes', 'pred_mr_top_nodes',
                    'common_connections_corr_top_nodes','common_connections_mr_top_nodes']

common_connection_summary = func.create_top_connection_summary(list_of_dicts,list_of_rownames)
common_connection_summary.to_csv(output_folder+'/7_gene_connection_hubs_'+gene_list_path+'.csv')
