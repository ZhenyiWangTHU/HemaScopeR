import os
import pandas as pd
import scanpy as sc
import stereo as st
import warnings
warnings.filterwarnings('ignore')

def read_stereo(
    data_path,
    save_path = '.',
    data_type = 'gem',
    sep = '\t',
    bin_type = 'bins',
    bin_size = 100,
    is_sparse = True,
    gene_list: list = None, 
    region: list = None
):
    if data_type is 'gem':
        adata = st.io.read_gem(
            file_path=data_path，
            sep=sep,
            bin_type=bin_type,
            bin_size=bin_size,
            is_sparse=is_sparse
        )
    elif data_path is 'gef':
        adata = st.io.read_gef(
            file_path=data_path,
            bin_type=bin_type,
            bin_size=bin_size,
            is_sparse=is_sparse,
            gene_list=gene_list,
            region=region,
        )
    elif data_path is 'h5ad':
        adata = st.io.read_stereo_h5ad(file_path=data_path)
    else:
        print(data_type + ' is not supported.')
        exit(1)
    
    adata = st.io.stereo_to_anndata(adata)
    sc.write(os.path.join(save_path, 'stereo_adata.h5ad'), 
             adata=adata)
    adata.obs.to_csv(os.path.join(save_path, 'stereo_metadata.csv'))
