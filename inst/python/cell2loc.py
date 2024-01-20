import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
import scvi

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

def run_cell2loc(
    st_data_path,
    sc_h5ad_path,
    save_path = '.',
    labels_key = 'seurat_clusters',
    sc_batch_key = None,
    st_batch_key = None, 
    use_gpu = True, 
    cell_count_cutoff=5, 
    cell_percentage_cutoff2=0.03, 
    nonz_mean_cutoff=1.12, 
    sc_max_epochs=1000,
    st_max_epochs=10000,
    species='human'
):

    ref_run_name = f'{save_path}/reference_signatures'
    run_name = f'{save_path}/cell2location_map'

    if species == 'human':
        mt_genes_head = 'MT-'
    elif species == 'mouse':
        mt_genes_head = 'mt-'

    ## st data
    adata_vis = sc.read_visium(st_data_path)
    adata_vis.var_names_make_unique()
    adata_vis.var['SYMBOL'] = adata_vis.var_names

    adata_vis.var['MT_gene'] = [gene.startswith(mt_genes_head) for gene in adata_vis.var['SYMBOL']]
    adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

    ## sc ref
    adata_ref = sc.read_h5ad(sc_h5ad_path)
    adata_ref.var_names_make_unique()
    adata_ref.var['SYMBOL'] = adata_ref.var_names

    adata_ref.var['MT_gene'] = [gene.startswith(mt_genes_head) for gene in adata_ref.var['SYMBOL']]
    adata_ref.obsm['MT'] = adata_ref[:, adata_ref.var['MT_gene'].values].X.toarray()
    adata_ref = adata_ref[:, ~adata_ref.var['MT_gene'].values]

    from cell2location.utils.filtering import filter_genes
    selected = filter_genes(
        adata_ref, 
        cell_count_cutoff=cell_count_cutoff, 
        cell_percentage_cutoff2=cell_percentage_cutoff2, 
        nonz_mean_cutoff=nonz_mean_cutoff
    )
    adata_ref = adata_ref[:, selected].copy()

    ## training
    cell2location.models.RegressionModel.setup_anndata(
        adata=adata_ref, 
        batch_key=sc_batch_key,
        labels_key=labels_key
    )
    from cell2location.models import RegressionModel
    mod = RegressionModel(adata_ref)

    mod.train(max_epochs=sc_max_epochs, use_gpu=use_gpu)
    adata_ref = mod.export_posterior(
        adata_ref, 
        sample_kwargs={'num_samples': 1000, 
                       'batch_size': 2500, 
                       'use_gpu': use_gpu}
    )

    # Save model
    mod.save(f"{ref_run_name}", overwrite=True)

    adata_file = f"{ref_run_name}/sc.h5ad"
    adata_ref.write(adata_file)

    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_ref.uns['mod']['factor_names']

    intersect = np.intersect1d(adata_vis.var_names, 
                            inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    cell2location.models.Cell2location.setup_anndata(
        adata=adata_vis, 
        batch_key=st_batch_key
    )
    mod = cell2location.models.Cell2location(
        adata_vis, cell_state_df=inf_aver,
        N_cells_per_location=30,
        detection_alpha=20
    )
    mod.train(max_epochs=st_max_epochs,
            batch_size=None,
            train_size=1,
            use_gpu=use_gpu)

    adata_vis = mod.export_posterior(
        adata_vis, 
        sample_kwargs={'num_samples': 1000, 
                       'batch_size': mod.adata.n_obs, 
                       'use_gpu': use_gpu}
    )

    # Save model
    mod.save(f"{run_name}", overwrite=True)

    ## save results
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
    adata_file = f"{save_path}/st.h5ad"
    
    del adata_vis.obsm
    del adata_vis.uns
    del adata_vis.var
    
    adata_vis.write(adata_file)
