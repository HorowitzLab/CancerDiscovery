#!/usr/bin/env python

"""The processing and analysis script for Cancer Discovery BCG NKG2A Story

This script follows the following analysis format
- Initial processing (normalization, log transformation, min and max cell filtering, neighbors analysis, and leiden clustering)
- Scanorama batch correction
- Concatenation
- 
"""

__author__ = "Daniel Ranti"
__credits__ = "Yuan-Sho Wang"
__license__ = "Open Access"
__version__ = "1.0.1"
__maintainer__ = "Daniel Ranti"
__email__ = "daniel.l.ranti@gmail.com"
__status__ = "Development"

# Standard library
import os
import logging
import datetime

# Third party imports
import pandas as pd
import numpy as np
import scipy

import anndata as ad
from anndata import AnnData
import scanorama
import tangram
import scanpy as sc
import scanpy.external as sce
import squidpy as sq
import tangram as tg

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text

# Local Imports
from spatial_tools import *

# FOR LOGGING
logger = logging.getLogger('simple_example')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

# FOR LATER ON WHEN WE WANT LOGFILES
# logging.basicConfig(filename='example.log', encoding='utf-8', level=logging.DEBUG)
# logging.debug('This message should go to the log file')
# logging.info('So should this')
# logging.warning('And this, too')
# logging.error('And non-ASCII stuff, too, like Øresund and Malmö')

# Scanpy Figure Settings
sc.set_figure_params(facecolor="white", figsize=(8, 8), transparent=True, fontsize=14, dpi_save=200, dpi=100)
sc.settings.verbosity = 3

def create_paths(spatial_path='wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01392_JohnSfakianos'):
    # Creating a list of targets along with my path to the visium data
    nmibc_path_list = []
    sample_ids = []
    for item in os.listdir(spatial_path):
        if os.path.isdir(os.path.join(spatial_path, item)) & ('bladder' in str(item).lower()):
            sample_ids.append(str(item))
            nmibc_path_list.append(spatial_path+'/'+item+'/outs')
    nmibc_sections = read_10x_path_list(nmibc_path_list, sample_ids)
    return nmibc_sections, sample_ids

def initial_processing(spatial_object):
    sc.pp.normalize_total(spatial_object, inplace=True)
    sc.pp.log1p(spatial_object)
    sc.pp.highly_variable_genes(spatial_object, flavor="seurat", n_top_genes=2000, inplace=True)
    sc.pp.filter_cells(spatial_object, min_counts=500,inplace=True)
    sc.pp.filter_cells(spatial_object, max_counts=35000, inplace=True)
    if len(spatial_object.obs) > 0: 
        sc.pp.neighbors(spatial_object, n_neighbors=15, n_pcs=50)
        sc.tl.leiden(spatial_object, key_added="clusters")
        return spatial_object

def object_concatenation(corrected_spatial):
    combined_2 = corrected_spatial[0].concatenate(
        corrected_spatial[1],
        uns_merge="unique",
        batch_categories=[
            k
            for d in [
                corrected_spatial[0].uns["spatial"],
                corrected_spatial[1].uns["spatial"],
            ]
            for k, v in d.items()
        ],
    )

    combined_1 = corrected_spatial[2].concatenate(
        corrected_spatial[3],
        uns_merge="unique",
        batch_categories=[
            k
            for d in [
                corrected_spatial[2].uns["spatial"],
                corrected_spatial[3].uns["spatial"],
            ]
            for k, v in d.items()
        ],
    )

    combined = combined_1.concatenate(
        combined_2,
        uns_merge="unique",
        batch_categories=[
            k
            for d in [
                combined_1[0].uns["spatial"],
                combined_2[1].uns["spatial"],
            ]
            for k, v in d.items()
        ],
    )
    return combined


def run_tangram(spatial, singleCell):
    sc.tl.rank_genes_groups(singleCell, groupby="cell_subclass", use_raw=False)

    markers_df = pd.DataFrame(singleCell.uns["rank_genes_groups"]["names"]).iloc[0:100, :]
    genes_sc = np.unique(markers_df.melt().value.values)
    genes_st = spatial.var_names.values
    genes = list(set(genes_sc).intersection(set(genes_st)))
    
    tg.pp_adatas(singleCell, spatial, genes=genes)

    ad_map = tg.map_cells_to_space(
        singleCell,
        spatial,
        mode="constrained",
        target_count=spatial.obs.cell_count.sum(),
        density_prior=np.array(spatial.obs.cell_count) / spatial.obs.cell_count.sum(),
        num_epochs=1000,
        device="cpu",
    )
    tg.project_cell_annotations(ad_map, spatial, annotation="cell_subclass")
    return spatial

    
if __name__ == '__main__': 
    # Quality control for each NMIBC figure
    logger.info('Starting QC')
    nmibc_sections, sample_ids = create_paths()
    processed = []
    for spatial_object in nmibc_sections:
        temp = initial_processing(spatial_object)
        processed.append(temp)
    
    # Picking out samples 8, 9, 10, 11
    logger.info('Starting concatenation and scanorama')
    nmibc_good_quality = [processed[i]for i in [0,4,5,9]]
    sample_ids_good_quality = [sample_ids[i]for i in [0,4,5,9]]
    corrected_spatial = scanorama.correct_scanpy(nmibc_good_quality, return_dimred=True)
    combined = object_concatenation(corrected_spatial)
    
    # neighbors / umap / leiden with all samples combined
    logger.info('Starting neighbors / umap / leiden')
    sc.pp.neighbors(combined, use_rep="X_scanorama")
    sc.tl.umap(combined)
    sc.tl.leiden(combined, key_added="clusters")

    # Tangram / labelling analysis
    logger.info('Starting tangram and deconvolution')
    labelled_sc = ad.read_h5ad('/sc/arion/projects/nmibc_bcg/data/spatial/scSeq_labelled/alice_converted.h5ad')

    combined_labelled = run_tangram(combined, labelled_sc)
    logger.info('Tangram Done. Saving the h5ad.')
    combined_labelled.write_h5ad('labelled_spatial_tangram.h5ad')
