import networkx as nx
import pandas as pd
import os
import matrix_network_functions as func

# 1) upload coexpression files created by make_coexpression_table.py

gene_list_path = 'flowering_list_full.txt'
output_folder = 'flowering_output'

try:
    os.mkdir(output_folder)
except:
    print('Output folder', output_folder, 'already exists')

target_log2_coexpression = pd.read_csv(output_folder+'/3_GSEA_correlation_'+gene_list_path+'.csv').set_index('gene_id')
pred_coexpression = pd.read_csv(output_folder+'/3_GAPIT_correlation_'+gene_list_path+'.csv').set_index('Unnamed: 0')

target_log2_mr = pd.read_csv(output_folder+'/4_GSEA_mutual_rank_'+gene_list_path+'.csv').set_index('gene_id')
pred_log2_mr = pd.read_csv(output_folder+'/4_GSEA_mutual_rank_'+gene_list_path+'.csv').set_index('gene_id')

target_corr_bin = pd.read_csv(output_folder+'/5_GSEA_binary_shortlist_corr_'+gene_list_path+'.csv')
target_mr_bin = pd.read_csv(output_folder+'/5_GSEA_binary_shortlist_mr_'+gene_list_path+'.csv')
pred_corr_bin = pd.read_csv(output_folder+'/5_GAPIT_binary_shortlist_corr_'+gene_list_path+'.csv')
pred_mr_bin = pd.read_csv(output_folder+'/5_GAPIT_binary_shortlist_mr_'+gene_list_path+'.csv')

# 2) create graph objects and export them as .gml

# select correlation treshold for the correlation matrix
corr_threshold = 0.5
# select top edges from the mutual rank matrix
top_n_mr = 20

graph_target_corr = nx.Graph()
graph_target_corr.add_weighted_edges_from(func.dataframe_to_edgelist(target_log2_coexpression, corr_threshold, mr = False))

nx.write_gml(graph_target_corr, output_folder+'/6_GSEA_graph_corr_'+gene_list_path+'.gml')

graph_pred_corr = nx.Graph()
graph_pred_corr.add_weighted_edges_from(func.dataframe_to_edgelist(pred_coexpression, corr_threshold, mr = False))

nx.write_gml(graph_pred_corr, output_folder+'/6_GAPIT_graph_corr_'+gene_list_path+'.gml')

graph_target_mr = nx.Graph()
graph_target_mr.add_weighted_edges_from(func.dataframe_to_edgelist(target_log2_mr, top_n_mr, mr = True))

nx.write_gml(graph_target_corr, output_folder+'/6_GSEA_graph_mr_'+gene_list_path+'.gml')

graph_pred_mr = nx.Graph()
graph_pred_mr.add_weighted_edges_from(func.dataframe_to_edgelist(pred_log2_mr, top_n_mr, mr = True))

nx.write_gml(graph_pred_corr, output_folder+'/6_GAPIT_graph_mr_'+gene_list_path+'.gml')
