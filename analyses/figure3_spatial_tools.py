"""
Spatial analysis tools
"""

__author__ = "Daniel Ranti"
__credits__ = "Yuan-Sho Wang"
__license__ = "Open Access"
__version__ = "1.0.1"
__maintainer__ = "Daniel Ranti"
__email__ = "daniel.l.ranti@gmail.com"
__status__ = "Development"

import scanpy as sc
import scanpy.external as sce

# for reading / unzipping / writing - eventually will write tools with this
import gzip
import shutil

import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import os
import pandas as pd
import numpy as np
import scanpy as sc
from anndata import AnnData
import scipy
import tangram as tg
from matplotlib.lines import Line2D
from adjustText import adjust_text


def merge_cell_idents(master_mtx, subset_mtx, subset_type):
    """
    master_mtx - the matrix without the identities of the cells, but all other obs
    subset_mtx - the sampletype specific matrix with SingleR identities
    subset_type - the cell subset type
    """
    # unique_t_cells = identities.obs.index
    master_mtx = master_mtx[master_mtx.obs["SampleType"].isin([subset_type,])]
    copy = master_mtx[master_mtx.obs.index.isin(subset_mtx.obs.index)]
    # the copy.obs.index df is the same thing if the set of below is = true
    set(copy.obs.index == subset_mtx.obs.index)
    # RESETTING THE PATHOLOGY OBS OF IDENTITIES MATRIX - WHEN CONVERTING BETWEEN SEURAT AND SCANPY, IT GOT TURNED INTO AN INT
    subset_mtx.obs.Pathology = copy.obs.Pathology
    subset_mtx.obs.CellType = copy.obs.CellType
    subset_mtx.obs.SampleType = copy.obs.SampleType
    return subset_mtx


def read_10x_path_list(path_list, sample_ids):
    """
    reads in a path list and a list of samples IDs for each path
    outputs a list of visum adata objects, each labelled with the 
    associated ID
    """
    visium_adata_list = []
    for path, sample_id in zip(path_list, sample_ids):
        temp = sc.read_visium(path)
        temp.obs["sample_id"] = sample_id
        temp.var_names_make_unique()
        sc.pp.calculate_qc_metrics(temp, inplace=True)
        visium_adata_list.append(temp)
    return visium_adata_list


def gen_mpl_labels(
    adata, groupby, exclude=(), ax=None, adjust_kwargs=None, text_kwargs=None
):
    """
    Automated label adjusting for on-plot annotation in scanpy 
    """
    if adjust_kwargs is None:
        adjust_kwargs = {"text_from_points": False}
    if text_kwargs is None:
        text_kwargs = {}

    medians = {}

    for g, g_idx in adata.obs.groupby(groupby).groups.items():
        if g in exclude:
            continue
        medians[g] = np.median(adata[g_idx].obsm["X_umap"], axis=0)

    if ax is None:
        texts = [
            plt.text(x=x, y=y, s=k, **text_kwargs) for k, (x, y) in medians.items()
        ]
    else:
        texts = [ax.text(x=x, y=y, s=k, **text_kwargs) for k, (x, y) in medians.items()]

    adjust_text(texts, **adjust_kwargs)


# Manually label the cells based on a cutoff of 90th percentile)
def manual_cell_label_nmibc(row, label_dict, percentile=90):
    label = []
    for item in label_dict.keys():
        cutoff = np.percentile(combined.obs[item], percentile)
        row_val = row[item]
        if row_val > cutoff:
            label.append(item)
    if len(label) > 2:
        # need to think about this one. for now just returning "multiple"
        if label_dict == lineage_genes:
            return "multiple"
        else:
            return " | ".join(label)
    elif len(label) == 1:
        return label[0]
    else:
        return "unlabelled"


def manual_cell_label(row, percentile=90):
    high_tumor_cutoff = np.percentile(mibc_spatial.obs["Tumor Cells"], percentile)
    high_nk_cutoff = np.percentile(mibc_spatial.obs["NK cells"], percentile)
    high_t_cutoff = np.percentile(mibc_spatial.obs["CD8 T-cells"], percentile)
    tumor_val = row["Tumor Cells"]
    t_val = row["CD8 T-cells"]
    nk_val = row["NK cells"]
    label = ""
    if tumor_val > high_tumor_cutoff:
        label = label + ("high tumor")
    if nk_val > high_nk_cutoff:
        label = label + (" high NK")
    if t_val > high_t_cutoff:
        label = label + (" high T")
    if label == "":
        label = "Unassigned"

    return label


def label_tertiles(row):
    row_dist = row["nearest_tumor_cell"]
    if row["manual_label_NKtumorT"] == "CD8 T-cells":
        row_dist = int(row_dist)
        if row_dist <= 112:
            label = "Closest (First Tertile)"
        elif (row_dist > 112) & (row_dist <= 231):
            label = "Middle (Second Tertile)"
        elif row_dist > 231:
            label = "Furthest (Third tertile)"
        return label


def label_nk_tertiles(row):
    """
    Based on nk_cells.obs['nearest_tumor_cell'].quantile([0.33, 0.66])
    first tertile cutoff: 170
    third tertile cutoff: 636
    """
    row_dist = row["nearest_tumor_cell"]
    row_dist = int(row_dist)
    if row_dist <= 170:
        label = "Closest (First Tertile)"
    elif (row_dist > 170) & (row_dist <= 636):
        label = "Middle (Second Tertile)"
    elif row_dist > 636:
        label = "Furthest (Third tertile)"
    return label
