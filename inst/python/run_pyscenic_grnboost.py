import os
import pandas as pd
from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names
#from dask.distributed import Client, LocalCluster

if __name__ == '__main__':
    # local_cluster = LocalCluster(n_workers=31, 
    #                              threads_per_worker=1)
    # custom_client = Client(local_cluster)
    # print(custom_client)
    ex_matrix = pd.read_csv("./int/1.1_exprMatrix_filtered_t.txt", sep='\t')
    tf_names = load_tf_names("./int/1.1_inputTFs.txt")
    network = grnboost2(expression_data=ex_matrix,tf_names=tf_names)
    network.to_csv('./03.grnboot2_network_links.tsv',sep='\t',header=True,index=False)