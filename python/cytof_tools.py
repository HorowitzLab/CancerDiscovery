#!/usr/bin/env python

"""
Preprocessing functions and various analyses tools for the CYTOF data
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
import logging

# Third Party Imports
from anndata import AnnData
import anndata as ad
import scanpy as sc
from FlowCytometryTools import FCMeasurement

# FOR LOGGING
logger = logging.getLogger("CYTOF_tools")
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)


def _prepare_cytof(fcs_list, clinical_data, condition_list, channels=None):
    """
    fcs_list: a list of FCS structures 
    clinical_data: a pandas dataframe to add to the FCS obs. 
    channels: a list of valid channels to filter the FCS.data object by
    condition_list: a condition to add to the obs manually. 
    """
    logger.info("Preparing cytof data")
    adata_list = []
    for fcs, condition in zip(fcs_list, condition_list):
        filter_list = []
        for sample in clinical_data.sampleID:
            if str(sample) in str(fcs.ID):
                filter_list.append(sample)
        matched_clinical = clinical_data[clinical_data.sampleID.isin(filter_list)]
        match_num = matched_clinical.shape[0]
        if match_num == 1:
            if channels:
                temp = AnnData(fcs.data[channels])
            else:
                temp = AnnData(fcs.data)
            for key, value in matched_clinical.to_dict(orient="records")[0].items():
                temp.obs[key] = value
            temp.obs["condition"] = condition
            temp.obs["sample_file_id"] = str(fcs.ID)
            temp.var_names_make_unique()
        adata_list.append(temp)
    return adata_list


def _concatenate_cytof(
    adata_list, condition_list, condition_column, timepoint_list, timepoint_column
):
    """
    adata_list: list of adata structures
    condition_list: list of conditions to filter by
    condition_column: column in adata.obs to filter conditions by
    timepoint_list: list of timepoints to filter by
    timepoint_column: column in adata.obs to filter timepoints by


    """
    logger.info("Concatenating cytof data")
    # Concatente everything
    total_mtx = adata_list[0]
    for mtx in adata_list[1:]:
        total_mtx = AnnData.concatenate(total_mtx, mtx, join="outer")
    filtered_1 = total_mtx[total_mtx.obs[condition_column].isin(condition_list)]
    final = filtered_1[filtered_1.obs[timepoint_column].isin(timepoint_list)]
    return final


def _normalize_cytof(adata, scaling_factor=5):
    """
    https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1917-7#Sec15
    normalization based on the above
    this article used a scaling factor of 5, hence the default
    """
    logger.info("Normalizing cytof data")
    return np.arcsinh(adata.X / scaling_factor)


def _batch_correct_cytof(adata, batch_key, covariates):
    """
    Using a python implementation of:  
        Johnson WE, Rabinovic A, Li C (2007). Adjusting batch effects in microarray
        expression data using Empirical Bayes methods. Biostatistics 8:118-127.  
    COMBAT chosen based on results of analysis showing efficacy of scanpy batch correction
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7338326/
    
    adata
    batch column: column in obs to batch correct on

    """
    logger.info("Batch correcting cytof data")
    return sc.pp.combat(adata, key=batch_key, covariates=covariates)


def preprocess_cytof(fcs_list, clinical_data, filter_dict, batch_key, channels=None):
    """
    fcs_list: a LIST of FCS structures 
    clinical_data: PD DATAFRAME  
    filter_dict requires 5 keys: 
        condition_list
        condition_column
        timepoint_list
        timepoint_column 
    batch_key = the sample ID column
    channels: a LIST of valid channels to filter the FCS data channels by
    """
    condition_list = filter_dict["condition_list"]
    condition_column = filter_dict["condition_column"]
    timepoint_list = filter_dict["timepoint_list"]
    timepoint_column = filter_dict["timepoint_column"]
    sample_column = filter_dict["sample_column"]

    adata_list = _prepare_cytof(
        fcs_list, clinical_data, condition_list, channels=channels,
    )
    concatenated = _concatenate_cytof(
        adata_list, condition_list, condition_column, timepoint_list, timepoint_column
    )
    normalized = _normalize_cytof(concatenated)
    # TODO: any other covariates we need to take into account?
    processed = _batch_correct_cytof(
        normalized, batch_key, covariates=[condition_column, timepoint_column]
    )

    return processed
