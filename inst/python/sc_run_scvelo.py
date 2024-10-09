# numpy version == 1.23.5
# pandas version == 1.3.5
import scvelo as scv
import scanpy as sc
import anndata as ad
import pandas as pd
import os
import pickle
scv.settings.set_figure_params('scvelo')
sdata = ad.read_mtx(os.path.join(outputDir, 'Step12.Construct_trajectories','scVelo','splicedmat.mtx'))
sdata.layers['spliced'] = sdata.X
ndata = ad.read_mtx(os.path.join(outputDir, 'Step12.Construct_trajectories','scVelo','umat.mtx'))
sdata.layers['unspliced'] = ndata.X
mdata = ad.read_mtx(os.path.join(outputDir, 'Step12.Construct_trajectories','scVelo','amat.mtx'))
sdata.layers['ambiguous'] = mdata.X
pd_obs = pd.read_csv(os.path.join(outputDir, 'Step12.Construct_trajectories','scVelo','scvelo_combined.obs.csv'))
pd_obsm = pd.read_csv(os.path.join(outputDir, 'Step12.Construct_trajectories','scVelo','cell.embeddings.csv'))
sdata.obs = pd_obs
sdata.obsm["X_umap"] = pd_obsm.values
gene_info = pd.read_csv(os.path.join(outputDir, 'Step12.Construct_trajectories','scVelo','geneInfo.csv'))
gene_info.set_index(["x"], inplace=True)
sdata.var = gene_info

scv.pp.filter_and_normalize(sdata, min_counts=10, min_counts_u=10, 
                            min_cells=10, min_cells_u=10, 
                            min_shared_counts=None, min_shared_cells=None, 
                            n_top_genes=None, flavor='seurat', log=True, copy=False)

sc.pp.pca(sdata)
sc.pp.neighbors(sdata, n_pcs=30, n_neighbors=50)

scv.pp.moments(sdata,n_pcs=None, n_neighbors=None)

scv.tl.recover_dynamics(sdata, n_jobs=10)

# scv.pl.scatter(sdata, color='latent_time', fontsize=24, size=100,
#               color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1])

scv.tl.velocity_graph(sdata, basis='umap')

scv.tl.velocity_embedding(sdata, basis='umap')

scv.pl.velocity_embedding_grid(sdata, basis='umap',color=['cellTypes'], save=os.path.join(outputDir, 'Step12.Construct_trajectories','scVelo','velocity_embedding_grid.pdf'),
                                                             palette=['#4292c6','#08519c','#c994c7'])

scv.pl.velocity_embedding_stream(sdata, basis='umap',color=['cellTypes'], save=os.path.join(outputDir, 'Step12.Construct_trajectories','scVelo','velocity_embedding_stream.pdf'),
                                palette=['#4292c6','#08519c','#c994c7'])

# sdata.obs['root_cells'] = sdata.obs['cellTypes']

# scv.tl.recover_latent_time(sdata)

# scv.tl.velocity_confidence(sdata)
# keys = 'velocity_length', 'velocity_confidence'
# scv.pl.scatter(sdata, c=keys, cmap='coolwarm', perc=[5, 95], save='velocity_length_confidence.pdf')
