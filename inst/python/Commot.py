import os
import gc
import ot
import pickle
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt

import seaborn as sns

import commot as ct


def run_Commot(
    st_data_path,
    metadata_path = None,
    label_key = None,
    save_path = '.',
    species = 'mouse',
    signaling_type = 'Secreted Signaling',
    database = 'CellChat',
    min_cell_pct = 0.05,
    dis_thr = 500,
    n_permutations = 1000
):
    adata = sc.read_visium(st_data_path)
    adata.var_names_make_unique()
    adata.raw = adata
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    
    if label_key is None:
        adata_metadata = adata.copy()
        sc.pp.highly_variable_genes(adata_metadata, 
                                    min_mean=0.0125, 
                                    max_mean=3, 
                                    min_disp=0.5)
        adata_metadata = adata_metadata[:, adata_metadata.var.highly_variable]
        sc.tl.pca(adata_metadata, svd_solver='arpack')
        sc.pp.neighbors(adata_metadata, n_neighbors=10, n_pcs=40)
        sc.tl.umap(adata_metadata)
        sc.tl.leiden(adata_metadata, resolution=0.6)
#         sc.pl.spatial(adata_metadata, color='leiden')
        
        label_key = 'leiden'
        adata.obs['leiden'] = adata_metadata.obs['leiden']
        
        del adata_metadata
    else:
        print('Loading metadata from ' + metadata_path)
        metadata = pd.read_csv(metadata_path, index_col=0)
        
        adata = adata[metadata.index, :]
        adata.obs[label_key] = pd.Categorical(metadata[label_key].astype(str))
    
    ct.tl.cluster_position(adata, clustering=label_key)
    
    df = ct.pp.ligand_receptor_database(species=species, 
                                        signaling_type=signaling_type, 
                                        database=database)
    df_filtered = ct.pp.filter_lr_database(df, 
                                           adata, 
                                           min_cell_pct=min_cell_pct)
    ct.tl.spatial_communication(
        adata,
        database_name=database,
        df_ligrec=df_filtered,
        dis_thr=dis_thr,
        heteromeric=True, 
        pathway_sum=True
    )
    
    # Calculate the strength of pathways within each cluster
    print('Calculating the strength of pathways within each cluster')
    pathway_names = np.unique(df_filtered.loc[:, 2])
    
    within_cluster_commot_matrix = {}
    within_cluster_commot_pvalue = {}
    
    for pathway_name in pathway_names:
        ct.tl.cluster_communication(
            adata, 
            database_name=database, 
            pathway_name=pathway_name, 
            clustering=label_key,
            n_permutations=n_permutations
        )
        communication_pvalue = adata.uns['commot_cluster-' + label_key + '-' + 
                                         database + '-' + pathway_name]['communication_pvalue']
        communication_matrix = adata.uns['commot_cluster-' + label_key + '-' + 
                                         database + '-' + pathway_name]['communication_matrix']

        communication_matrix.to_csv(save_path + r'/Pathway_matrix/' + pathway_name + '.csv')
        communication_pvalue.to_csv(save_path + r'/Pathway_pvalue/' + pathway_name + '.csv')
                                     
        within_cluster_pathway_matrix = []
        within_cluster_pathway_pvalue = []
        for i in range(communication_matrix.shape[0]):
            within_cluster_pathway_matrix.append(communication_matrix.iloc[i, i])
            within_cluster_pathway_pvalue.append(communication_pvalue.iloc[i, i])
        within_cluster_commot_matrix[pathway_name] = within_cluster_pathway_matrix
        within_cluster_commot_pvalue[pathway_name] = within_cluster_pathway_pvalue
        
        cluster_index = communication_matrix.index
        
    within_cluster_commot_matrix = pd.DataFrame.from_dict(within_cluster_commot_matrix)
    within_cluster_commot_pvalue = pd.DataFrame.from_dict(within_cluster_commot_pvalue)
    within_cluster_commot_matrix.index = cluster_index
    within_cluster_commot_pvalue.index = cluster_index
    
    within_cluster_commot_matrix.to_csv(save_path + r'/Pathway_within_cluster_matrix.csv')
    within_cluster_commot_pvalue.to_csv(save_path + r'/Pathway_within_cluster_pvalue.csv')
    
    # Calculate the strength of LR pairs within each cluster
    print('Calculating the strength of LR pairs within each cluster')
    n_LR_pairs = df_filtered.shape[0]

    within_cluster_commot_matrix = {}
    within_cluster_commot_pvalue = {}
    
    for i in range(n_LR_pairs):
        ct.tl.cluster_communication(
            adata, 
            database_name=database, 
            lr_pair=[df_filtered.loc[i, 0], 
                     df_filtered.loc[i, 1]], 
            clustering=label_key,
            n_permutations=n_permutations
        )
        
        LR_pair = df_filtered.loc[i, 0] + '-' + df_filtered.loc[i, 1]
        communication_pvalue = adata.uns['commot_cluster-' + label_key + '-' + 
                                         database + '-' + LR_pair]['communication_pvalue']
        communication_matrix = adata.uns['commot_cluster-' + label_key + '-' + 
                                         database + '-' + LR_pair]['communication_matrix']
        
        communication_matrix.to_csv(save_path + r'/LR_matrix/' + LR_pair + '.csv')
        communication_pvalue.to_csv(save_path + r'/LR_pvalue/' + LR_pair + '.csv')
                                         
        within_cluster_pathway_matrix = []
        within_cluster_pathway_pvalue = []
        for i in range(communication_matrix.shape[0]):
            within_cluster_pathway_matrix.append(communication_matrix.iloc[i, i])
            within_cluster_pathway_pvalue.append(communication_pvalue.iloc[i, i])
        within_cluster_commot_matrix[LR_pair] = within_cluster_pathway_matrix
        within_cluster_commot_pvalue[LR_pair] = within_cluster_pathway_pvalue
        
        cluster_index = communication_matrix.index
    
    within_cluster_commot_matrix = pd.DataFrame.from_dict(within_cluster_commot_matrix)
    within_cluster_commot_pvalue = pd.DataFrame.from_dict(within_cluster_commot_pvalue)
    within_cluster_commot_matrix.index = cluster_index
    within_cluster_commot_pvalue.index = cluster_index

    within_cluster_commot_matrix.to_csv(save_path + r'/LR_within_cluster_matrix.csv')
    within_cluster_commot_pvalue.to_csv(save_path + r'/LR_within_cluster_pvalue.csv')
    
    adata.write_h5ad(save_path + '/adata_Commot.h5ad')
    
    # # Draw the strength of each pathway on the image
    # print('Drawing the strength of each pathway on the image')
    # 
    # os.chdir(save_path)
    # if not os.path.exists('./pathway_plot'):
    #     os.mkdir('./pathway_plot')
    # os.chdir('./pathway_plot')
    # 
    # for pathway_name in pathway_names:
    #     if 'commot_cluster-' + label_key + '-' + database + '-' + pathway_name in adata.uns_keys():
    #         ct.pl.plot_cluster_communication_network(
    #             adata, 
    #             uns_names=['commot_cluster-' + label_key + '-' + database + '-' + pathway_name],
    #             clustering=label_key,
    #             nx_node_pos='cluster', 
    #             nx_pos_idx=np.array([0, 1]), 
    #             nx_bg_pos=True, 
    #             nx_bg_ndsize=0.25, 
    #             p_value_cutoff=5e-2,
    #             nx_node_size=0.5,
    #             nx_edge_width_lb_quantile = 0.05,
    #             nx_edge_width_ub_quantile = 0.95,
    #             nx_edge_width_min = 5,
    #             nx_edge_width_max = 10,
    #             filename=pathway_name+'_cluster_spatial.pdf', 
    #             nx_node_cmap='Light24'
    #         )
