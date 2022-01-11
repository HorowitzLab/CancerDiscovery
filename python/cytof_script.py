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

# Standard Imports
import os
import pandas as pd
import numpy as np
import scipy
import logging

# Third Party Imports
from anndata import AnnData
import scanpy as sc
import scanpy.external as sce
import graphtools as gt
import phate
# import diffxpy.api as de
# import magic
import scprep
import meld
import cmocean
import sklearn
from FlowCytometryTools import FCMeasurement

# Plotting Imports
import seaborn as sns
import matplotlib.pyplot as plt
plt.rc('font', size=14)

# FOR LOGGING
logger = logging.getLogger('spatial_script')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

# making sure plots & clusters are reproducible
np.random.seed(42)

def _replicate_normalize_densities(sample_densities, replicate):
    replicates = np.unique(replicate)
    sample_likelihoods = sample_densities.copy()
    for rep in replicates:
        curr_cols = sample_densities.columns[[col.endswith(rep) for col in sample_densities.columns]]
        sample_likelihoods[curr_cols] = sklearn.preprocessing.normalize(sample_densities[curr_cols], norm='l1')
    return sample_likelihoods

def _plot_jitter_by_cluster(metadata, sample_cmap, cluster_key, condition_key):
    fig, ax = plt.subplots(1, figsize=(10,10))

    # See example usage: https://scprep.readthedocs.io/en/stable/examples/jitter.html
    scprep.plot.jitter(metadata[cluster_key], metadata['chd_likelihood'], c=metadata[condition_key], 
                    cmap=sample_cmap,legend=False, plot_means=False, xlabel=False, ylabel='Mean chd likelihood',
                    ax=ax)

    ### This code will plot the ratio of pre:post cells per cluster
    means = metadata.groupby('clusterID')['genotype'].mean()
    ax.scatter(means.index, means - np.mean(metadata['genotype']) + 0.5, color='#e0d8f0', edgecolor='k', s=100)

    # Axis tick labels
    ax.set_xticklabels(metadata.set_index('clusterID')['cluster'].drop_duplicates().sort_index(), rotation=90)
    ax.set_ylim(0,1)

    fig.tight_layout()
    return fig


def run_meld_cytof(combined_adata, condition1, condition2, cluster_key, condition_key='sample_type'):
    '''
    The entire meld algorithm start to finish will be contained in here. 
    combined_adata: Anndata structure containing your conditions
    condition1: string 
    condition2: string2
    condition_key
    '''
    logger.info('MELD: Sample Counts')
    sample_cmap = {
        condition1: '#fb6a4a',
        condition2: '#de2d26',
#         'chdC': '#a50f15',
#         'tyrA': '#6baed6',
#         'tyrB': '#3182bd',
#         'tyrC': '#08519c'
    }
    fig, ax = plt.subplots(1)
    metadata = combined_adata.obs
    data = combined_adata.to_df()

    groups, counts = np.unique(metadata[condition_key], return_counts=True)
    for i, c in enumerate(counts):
        ax.bar(i, c, color=sample_cmap[groups[i]])

    ax.set_xticks(np.arange(i+1))
    ax.set_xticklabels(groups)
    ax.set_ylabel('Cells Counts Per Condition')
    
    fig.tight_layout()
    fig.save_fig('meld_condition_counts_{}_{}_{}.png'.format(condition_key,condition1, condition2),dpi=200)


    # Meld Run Phate
    logger.info('MELD: Run Phate')
    data_pca = scprep.reduce.pca(data)
    phate_op = phate.PHATE(n_jobs=-1)
    data_phate = phate_op.fit_transform(data_pca)
    fig = scprep.plot.scatter2d(data_phate, c=metadata[condition_key], cmap=sample_cmap, 
                      legend_anchor=(1,1), figsize=(6,5),dpi=200, s=10, label_prefix='PHATE', ticks=False)
    fig.save_fig('meld_phate_{}_{}_{}.png'.format(condition_key,condition1, condition2),dpi=200)
    
    # Run Meld
    logger.info('MELD: Run MELD')
    metadata['condition'] = [1 if sl.startswith(condition2) else 0 for sl in metadata[condition_key]]
    metadata['condition_name'] = ['Pre-BCG' if g == 0 else 'Post-BCG' for g in metadata['condition']]
    metadata['replicate'] = metadata['sample_file_id']

    # G = gt.Graph(data_pca, knn=int(top_result['knn']), use_pygsp=True)
    # TODO: need to implement the parameter search. Using published paramters here. 
    meld_op = meld.MELD(beta=67, knn=7)
    sample_densities = meld_op.fit_transform(data_pca, sample_labels=metadata[condition_key])
    # Normalizing across samples. 
    logger.info('MELD: Normalizing Samples')
    sample_likelihoods = _replicate_normalize_densities(sample_densities, metadata['replicate'])

    fig, axes = plt.subplots(1,3, figsize=(13,4))
    experimental_samples = [condition1, condition2]

    for i, ax in enumerate(axes):
        curr_sample = experimental_samples[i]
        scprep.plot.scatter2d(data_phate, c=sample_likelihoods[curr_sample], cmap=meld.get_meld_cmap(),
                            vmin=0, vmax=1,
                            title=curr_sample, ticks=False, ax=ax)

    fig.tight_layout()
    fig.save_fig('meld_phate_vfc__{}_{}_{}.png'.format(condition_key,condition1, condition2),dpi=200)
    fig = _plot_jitter_by_cluster(metadata, sample_cmap, cluster_key=cluster_key, condition_key=condition_key)
    fig.save_fig('meld_jitter_vfc__{}_{}_{}.png'.format(condition_key,condition1, condition2),dpi=200)
    
if __name__ == '__main__':
    # Listing and importing the FCS files; extracting their condition as well
    fcs_list = []
    condition_list = []
    directory = os.fsencode('.')
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".fcs"): 
            fcs_list.append(FCMeasurement(ID=file, datafile=filename))
            condition=''
            if 'K562' in filename:
                condition='K562'
            elif 'RAJI' in filename: 
                condition='RAJI'
            elif 'MED' in filename:
                condition='MED'
            condition_list.append(condition)

    # Too many channels - limiting to those conjugated to metals we want        
    channels_of_interest = ['111Cd_LILRB1','112Cd_Granzyme_A','113Cd_CD38','114Cd_CD11b','115In_IFNg','116Cd_CD57','141Pr_cKit','142Nd_KLRG1','143Nd_Granzyme_K','144Nd_CD69','145Nd_NKG2D_','146Nd_DNAM1','147Sm_NKp80','148Nd_KIR3DL1_L2','149Sm_KIR3DL1','150Nd_IL2','151Eu_CD107a','152Sm_TNFa','153Eu_CD14','154Sm_MIP1b','155Gd_NKp46','156Gd_Tim3','158Gd_KIR2DL1','159Tb_CD56','160Gd_NKG2A','161Dy_NKp44','162Dy_CD27','163Dy_Eomes','164Dy_NKG2C','165Ho_KIR2DL3','166Er_KIR2DL1_S1','167Er_KIR2DL2_L3','168Er_TRAIL','169Tm_Tbet','170Er_CD3','171Yb_CD127','172Yb_Perforin','173Yb_Granzyme_B','174Yb_TIGIT','175Lu_XCL1','176Yb_CD200R','209Bi_CD16',]

    # Import Berengere's data
    clinical = pd.read_excel('id_list_cytof_blca_bcg.xlsx')

    # Debugging - need to remove one non-specific AHBCG17
    clinical.drop(53, inplace=True)

    # Generating a list of anndata structures
    adata_list = []
    for fcs, condition in zip(fcs_list, condition_list):
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
            temp.obs['condition'] = condition
            temp.obs['sample_file_id'] = str(fcs.ID)
            temp.var_names_make_unique() 
        adata_list.append(temp)    


    # TODO: ASK ABOUT BATCH CORRECTION!!! 
    # import scanpy as sc 
    # help(sc.pp.combat)
    logger.info('concatenating structures')
    # Concatenating multiple matrices
    pre_bcg = []
    recurrence = []
    for mtx in adata_list:
        sample_type = list(set(mtx.obs.sample_type))[0]
        if sample_type == 'I1D1':
            pre_bcg.append(mtx)
        if sample_type in 'at recurrence':
            recurrence.append(mtx)

    pre_mtx = pre_bcg[0]
    for mtx in pre_bcg[1:]:
        pre_mtx = AnnData.concatenate(
            pre_mtx,
            mtx,
            join='outer'
        )

    recurrent_mtx = recurrence[0]
    for mtx in recurrence[1:]:
        recurrent_mtx = AnnData.concatenate(
            recurrent_mtx,
            mtx,
            join='outer'
        )

    # limiting the mtx to med only.
    pre_mtx_medOnly = pre_mtx[pre_mtx.obs['condition'] == 'MED']
    recurrence_mtx_medOnly = recurrent_mtx[recurrent_mtx.obs['condition'] == 'MED']

    combined_mtx = AnnData.concatenate(
        pre_mtx_medOnly,
        recurrence_mtx_medOnly,
        join='outer'
    )

    # combined_mtx = AnnData.concatenate(
    #         recurrent_mtx,
    #         pre_mtx,
    #         join='outer'
    #     )

    # Running Phenograph
    logger.info('Running Phenograph')
    k = 30
    sc.tl.pca(combined_mtx, n_comps = 10)
    communities, graph, Q = sce.tl.phenograph(combined_mtx.obsm['X_pca'], k = k)
    combined_mtx.obs['PhenoGraph_clusters'] = pd.Categorical(communities)
    combined_mtx.uns['PhenoGraph_Q'] = Q
    combined_mtx.uns['PhenoGraph_k'] = k
    # Writing for future analysis

    # Drawing UMAP
    logger.info('Drawing UMAP')
    sc.pp.neighbors(combined_mtx, n_neighbors=30, n_pcs=10)
    sc.tl.umap(combined_mtx)
    sc.pl.umap(
        combined_mtx, 
        color=['PhenoGraph_clusters', 'sample_type'],
        title="PhenoGraph Assigned Clusters: Cytof Data Pre BCG vs Recurrence",
        save='phenograph test 01102022.png'
    )

    # combined_mtx.write_h5ad('combined_mtx_MED_only.h5ad')

    # MELD
    logger.info('Running Meld')
    run_meld_cytof(
        combined_adata=combined_mtx, 
        condition1='I1D1', 
        condition2='recurrence', 
        condition_key='sample_type', 
        cluster_key='PhenoGraph_clusters'
        )
