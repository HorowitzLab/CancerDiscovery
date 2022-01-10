#!/usr/bin/env python

"""
The initial processing and analysis script for the CYTOF data for the Cancer Discovery BCG Story
It follows the following order
- Import
- Concatenation
- Phenograph
- Differential Expression
TODO: Batch correction? Not sure what batches were run. 
TODO: implement meld
TODO: Ask Amir about any processing steps that were already done in cytobank. Not sure the quality of the data but assuming berengere already filtered and whatnot. 
"""

__author__ = "Daniel Ranti"
__credits__ = "Yuan-Sho Wang"
__license__ = "Open Access"
__version__ = "1.0.1"
__maintainer__ = "Daniel Ranti"
__email__ = "daniel.l.ranti@gmail.com"
__status__ = "Development"

import os
import pandas as pd
import numpy as np
import scipy
import diffxpy.api as de

from anndata import AnnData
import scanpy as sc
import graphtools as gt
import phate
# import magic
import scprep
import meld
import cmocean
import sklearn
from FlowCytometryTools import FCMeasurement

# setting defaults for matplotlib font sizes
import seaborn as sns
import matplotlib.pyplot as plt
plt.rc('font', size=14)

# making sure plots & clusters are reproducible
np.random.seed(42)

# Listing and importing the FCS files
fcs_list = []
directory = os.fsencode('.')
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".fcs") : 
        fcs_list.append(FCMeasurement(ID=file, datafile=filename))

# Too many channels - limiting to those conjugated to metals we want        
channels_of_interest = ['111Cd_LILRB1','112Cd_Granzyme_A','113Cd_CD38','114Cd_CD11b','115In_IFNg','116Cd_CD57','141Pr_cKit','142Nd_KLRG1','143Nd_Granzyme_K','144Nd_CD69','145Nd_NKG2D_','146Nd_DNAM1','147Sm_NKp80','148Nd_KIR3DL1_L2','149Sm_KIR3DL1','150Nd_IL2','151Eu_CD107a','152Sm_TNFa','153Eu_CD14','154Sm_MIP1b','155Gd_NKp46','156Gd_Tim3','158Gd_KIR2DL1','159Tb_CD56','160Gd_NKG2A','161Dy_NKp44','162Dy_CD27','163Dy_Eomes','164Dy_NKG2C','165Ho_KIR2DL3','166Er_KIR2DL1_S1','167Er_KIR2DL2_L3','168Er_TRAIL','169Tm_Tbet','170Er_CD3','171Yb_CD127','172Yb_Perforin','173Yb_Granzyme_B','174Yb_TIGIT','175Lu_XCL1','176Yb_CD200R','209Bi_CD16',]

# Import Berengere's data
clinical = pd.read_excel('id_list_cytof_blca_bcg.xlsx')

# Debugging - need to remove one non-specific AHBCG17
clinical.drop(53, inplace=True)

# Generating a list of anndata structures
adata_list = []
for fcs in fcs_list:
    filter_list = []
    for sample in clinical.sampleID:
        if str(sample) in str(fcs.ID):
            filter_list.append(sample)
    matched_clinical = clinical[clinical.sampleID.isin(filter_list)]
    match_num = matched_clinical.shape[0]
    if match_num == 1:      
        temp = AnnData(fcs.data[channels_of_interest])
        for key,value in matched_clinical.to_dict(orient="records")[0].items():
            temp.obs[key] = value
        temp.var_names_make_unique() 
    adata_list.append(temp)    

# NEED TO ASK ABOUT BATCH CORRECTION!!! 
# import scanpy as sc 
# help(sc.pp.combat)

# Concatenating multiple matrices
cumulative_mtx = adata_list[0]
for mtx in adata_list[1:]:
    cumulative_mtx = AnnData.concatenate(
        cumulative_mtx,
        mtx,
        join='outer'
    )

# Running Phenograph
k = 30
sc.tl.pca(cumulative_mtx, n_comps = 10)
communities, graph, Q = sce.tl.phenograph(cumulative_mtx.obsm['X_pca'], k = k)
cumulative_mtx.obs['PhenoGraph_clusters'] = pd.Categorical(communities)
cumulative_mtx.uns['PhenoGraph_Q'] = Q
cumulative_mtx.uns['PhenoGraph_k'] = k
# Writing for better future use
cumulative_mtx.write_h5ad('cytof_bcg_phenograph.h5ad')
