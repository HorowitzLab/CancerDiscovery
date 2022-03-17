"""The processing and analysis script for Cancer Discovery BCG NKG2A Story

This script follows the following analysis format
- Initial processing (normalization, log transformation, min and max cell filtering, neighbors analysis, and leiden clustering)
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

# import tangram
import scanpy as sc
import scanpy.external as sce
import squidpy as sq
import tangram as tg

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# from adjustText import adjust_text

# Local Imports
from spatial_tools import *
import customTGplotting

# FOR LOGGING
logger = logging.getLogger("SPATIAL_SCRIPT")
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)

# FOR LATER ON WHEN WE WANT LOGFILES
# logging.basicConfig(filename='example.log', encoding='utf-8', level=logging.DEBUG)
# logging.debug('This message should go to the log file')
# logging.info('So should this')
# logging.warning('And this, too')
# logging.error('And non-ASCII stuff, too, like Øresund and Malmö')

# Scanpy Figure Settings
sc.set_figure_params(
    facecolor="white",
    figsize=(8, 8),
    transparent=True,
    fontsize=14,
    dpi_save=200,
    dpi=100,
)
sc.settings.verbosity = 3


def cell_seg_count(adata):
    sq.im.process(img=adata, layer="image", method="smooth")
    sq.im.segment(
        img=img, layer="image_smooth", method="watershed", channel=0,
    )


def create_paths(
    spatial_path="wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01392_JohnSfakianos",
):
    """
    A very specific function to walk down the wangy spatial path and return scRNA objects
    Not modular at all lol 
    """
    # Creating a list of targets along with my path to the visium data
    nmibc_path_list = []
    sample_ids = []
    for item in os.listdir(spatial_path):
        if os.path.isdir(os.path.join(spatial_path, item)) & (
            "bladder" in str(item).lower()
        ):
            sample_ids.append(str(item))
            nmibc_path_list.append(spatial_path + "/" + item + "/outs")
    nmibc_sections = read_10x_path_list(nmibc_path_list, sample_ids)
    return nmibc_sections, sample_ids


def initial_processing(spatial_obs):
    """
    Process each spatial object according to standard processing guidelines
        - normalize data 
        - log transform
        - filter cells 250 < genes < 35000
    spatial_obs: a list of adata objects
    returns spatial_obs_filtered (nonzero adata objects filtered)
    """

    # Filtering each spatial object inplace
    # normalize; log transform; filter 1000 < x < 35000
    # Previously tried 250, but seems like a higher threshold is needed for
    # QC purposes
    # keep objects that passed filtering threshold
    spatial_obs_filtered = []
    for adata in spatial_obs:
        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)
        sc.pp.filter_cells(adata, min_counts=1000, inplace=True)
        sc.pp.filter_cells(adata, max_counts=35000, inplace=True)
        if len(adata.obs) > 0:
            spatial_obs_filtered.append(adata)
    return spatial_obs_filtered


def object_concatenation(adata_list):
    """
    Take in a list of spatial adatas and concatenate them
    return a single concatenated obect
    """
    # Concatenate all objects into a single AnnData Object
    concatenated = adata_list[0]
    for item in adata_list[1:]:
        #print(item.shape, item.shape[0])
        try: 
            if item.shape[0] > 0:
                concatenated = concatenated.concatenate(
                    item,
                    uns_merge="unique",
                    batch_categories=[
                        k
                        for d in [concatenated[0].uns["spatial"], item[1].uns["spatial"],]
                        for k, v in d.items()
                    ],
                )
        except (ValueError,IndexError):
            print('INDEX / VALUE ERROR')
            print(item)
            
    return concatenated


def run_tangram(spatial, singleCell, ngenes):
    '''
    spatial: the adata object you want to label
    singleCell: the single cell object with your cells of interest
    ngenes: the number of genes you want to take as your markers of interest
    '''
    sc.tl.rank_genes_groups(singleCell, groupby="cell.type", use_raw=False)

    markers_df = pd.DataFrame(singleCell.uns["rank_genes_groups"]["names"]).iloc[
        0:ngenes, :
    ]
    genes_sc = np.unique(markers_df.melt().value.values)
    genes_st = spatial.var_names.values
    genes = list(set(genes_sc).intersection(set(genes_st)))    
    tg.pp_adatas(singleCell, spatial, genes=genes)

    ad_map = tg.map_cells_to_space(
        singleCell,
        spatial,
        mode="clusters",
        cluster_label="cell.type",
        num_epochs=500,
        device="cpu",
    )
    tg.project_cell_annotations(ad_map, spatial, annotation="cell.type")
    return spatial


if __name__ == "__main__":
    # New Spatial Files
    jan22_spatial_files = [
        "13538",
        "16235",
        "20121",
        "27173",
        "28537",
        "30150",
        "37236",
        "37249",
    ]

    # current data dir
    data_dir = "/sc/arion/projects/nmibc_bcg/CancerDiscovery/data/spatial/"

    wangy_path = (
        "jan2022exp_realigned/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD005713_AmirHorowitz/"
    )
    target = data_dir + wangy_path
    nmibc_path_list = []
    sample_ids = []
    for item in os.listdir(target):
        if os.path.isdir(os.path.join(target, item)) & (
            str(item) in jan22_spatial_files
        ):
            sample_ids.append(str(item))
            nmibc_path_list.append(target + "/" + item + "/outs")
    visium_new = read_10x_path_list(nmibc_path_list, sample_ids)
    print(visium_new)

    # Old spatial_files
    wangy_path = "wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01392_JohnSfakianos"
    target = data_dir + wangy_path
    nmibc_path_list = []
    sample_ids = []
    for item in os.listdir(target):
        if os.path.isdir(os.path.join(target, item)) & ("bladder" in str(item).lower()):
            sample_ids.append(str(item))
            nmibc_path_list.append(target + "/" + item + "/outs")
    visium_old = read_10x_path_list(nmibc_path_list, sample_ids)
    print(visium_old)

    # Combine the extracted files from the experimental dirs
    spatial_obs = visium_new + visium_old
    spatial_obs_filtered = initial_processing(spatial_obs)

    ##########################################################
    ##########################################################

    # Importing the clinical data
    clinical_info = pd.read_excel(data_dir + "spatial_clinical_anon.xlsx")
    clinical_info["sample_id"] = clinical_info["BRP ID#"]
    clinical_info["timepoint"] = clinical_info["Characterization"]

    # Concatenate all objects into a single AnnData Object
    concatenated = object_concatenation(spatial_obs_filtered)
    sc.pp.filter_genes(concatenated, min_counts=concatenated.shape[0] * 0.01)

    # Fixing keys.
    for oldkey in concatenated.uns["spatial"].keys():
        if oldkey.startswith("st_count"):
            for newkey in set(concatenated.obs.sample_id):
                if newkey in oldkey:
                    print(newkey)
                    print(oldkey)
                    concatenated.uns["spatial"][newkey] = concatenated.uns[
                        "spatial"
                    ].pop(oldkey)
    # Joining to clinical data
    concatenated.obs["sample_id"] = concatenated.obs["sample_id"].astype("str")
    temp = pd.merge(
        concatenated.obs,
        clinical_info[["sample_id", "timepoint"]].astype(str),
        left_on="sample_id",
        right_on="sample_id",
    )
    temp.index = concatenated.obs.index
    concatenated.obs.timepoint = temp.timepoint
    
    logger.info("Starting tangram and deconvolution")
    labelled_sc = ad.read_h5ad(data_dir + "scSeq_labelled/Spatial/alice_converted.h5ad")
    # ranging the number of genes we use to label dots with
    for ngenes in [100,700,1000]:
        final = run_tangram(concatenated, labelled_sc, ngenes=ngenes)
        logger.info("Tangram Done.")
        df_plot = final.obsm['tangram_ct_pred']
        labels = (df_plot - df_plot.min()) / (df_plot.max() - df_plot.min())
        final.obs = pd.merge(final.obs, labels, left_index=True, right_index=True)
        final.write_h5ad("labelled_spatial_tangram_{}genes.h5ad".format(ngenes))
