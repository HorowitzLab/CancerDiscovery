#!/usr/bin/env python

"""
The initial processing and analysis script for the CYTOF data for the Cancer Discovery BCG Story
The analysis implements the following:
- Import FCS
- Concatenation FCS
- Filter by analysis type
- Preprocess: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1917-7#Sec15
    - Filtered on cytobank to remove annotation incompleteness, doublets, debris, and dead cells. 
    - Normalization by the inverse hyperbolic sine function with a scale factor of 5
    - Batch correction
- Phenograph
- Differential Expression
TODO: PUT ALL MELD FUNCTIONS INTO OWN PYTHON FILE
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
import anndata as ad
import scanpy as sc
import scanpy.external as sce
import graphtools as gt
import phate
import scprep
import meld
import cmocean
import sklearn
from FlowCytometryTools import FCMeasurement
from joblib import Parallel, delayed

# Plotting Imports
import seaborn as sns
import matplotlib.pyplot as plt

plt.rc("font", size=14)
sc.set_figure_params(
    facecolor="white",
    figsize=(8, 8),
    transparent=True,
    fontsize=14,
    dpi_save=200,
    dpi=100,
)
sc.settings.verbosity = 3

# Custom tools
import cytof_tools

# FOR LOGGING
logger = logging.getLogger("CYTOF_analysis_script")
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)

# making sure plots & clusters are reproducible
np.random.seed(42)

# Figure Directory
FIGDIR = "figures/"

def _replicate_normalize_densities(sample_densities, replicate):
    replicates = np.unique(replicate)
    sample_likelihoods = sample_densities.copy()
    for rep in replicates:
        curr_cols = sample_densities.columns[
            [col.endswith(rep) for col in sample_densities.columns]
        ]
        sample_likelihoods[curr_cols] = sklearn.preprocessing.normalize(
            sample_densities[curr_cols], norm="l1"
        )
    return sample_likelihoods


def _plot_jitter_by_cluster(metadata, sample_cmap, cluster_key, condition_key):
    fig, ax = plt.subplots(1, figsize=(10, 10))

    # See example usage: https://scprep.readthedocs.io/en/stable/examples/jitter.html
    scprep.plot.jitter(
        metadata[cluster_key],
        metadata["BCG_likelihood"],
        c=metadata[condition_key],
        #         cmap=sample_cmap,
        legend=False,
        plot_means=False,
        xlabel=False,
        ylabel="Mean BCG-state Likelihood",
        ax=ax,
    )

    ### This code will plot the ratio of pre:post cells per cluster
    means = metadata.groupby(cluster_key)["BCG_likelihood"].mean()
    ax.scatter(means.index, means, color="#e0d8f0", edgecolor="k", s=100)

    # Axis tick labels
    # TODO: Implement sorting better
    #     ax.set_xticklabels(means.sort_values().index.values, rotation=90)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Phenograph Assigned Clusters")
    ax.set_title("Post-BCG State Likelihood by Phenograph Cluster")

    fig.tight_layout()
    return fig

def _simulate_pdf_calculate_likelihood(benchmarker, seed, beta):
    """
    TODO ADD DOCSTRING
    """
    benchmarker.set_seed(seed)
    benchmarker.generate_ground_truth_pdf()

    benchmarker.generate_sample_labels()
    benchmarker.calculate_MELD_likelihood(beta=beta)
    MELD_mse = benchmarker.calculate_mse(benchmarker.expt_likelihood)
    return MELD_mse, seed, beta, benchmarker.graph.knn

def _parameter_search(adata, benchmarker):
    """
    TODO ADD DOCSTRING
    """
    knn_range = np.arange(1, 25)
    beta_range = np.arange(1, 200)

    results = []

    seed = 42

    with Parallel(n_jobs=36) as p:
        for knn in knn_range:
            # doing this outside the parallel loop because building the graph takes the longest
            benchmarker.fit_graph(adata.X, knn=knn)
            print(knn)
            curr_results = p(
                delayed(_simulate_pdf_calculate_likelihood)(benchmarker, seed, beta)
                for seed in range(25)
                for beta in beta_range
            )
            curr_results = pd.DataFrame(
                curr_results, columns=["MSE", "seed", "beta", "knn"]
            )
            results.append(curr_mse)

    results = pd.concat(results, axis=0)

    ax = scprep.plot.scatter(
        results_wide["beta"],
        results_wide["knn"],
        s=50,
        c=results_wide["MSE"],
        vmax=0.006,
        cmap="inferno_r",
    )

    # Highlight the top performing combination with a large red dot
    top_result = results_wide.sort_values("MSE").iloc[0]
    ax.scatter(
        top_result["beta"], top_result["knn"], c="r", s=100, linewidth=1, edgecolor="k"
    )

    return top_result


def run_meld_cytof(
    combined_adata, condition1, condition2, cluster_key, condition_key="timepoint"
):
    """
    The entire meld algorithm start to finish will be contained in here. 
    combined_adata: Anndata structure containing your conditions
    condition1: string 
    condition2: string2
    condition_key
    """
    logger.info("MELD: Sample Counts")
    sample_cmap = {condition1: "#fb6a4a", condition2: "#08519c"}
    fig, ax = plt.subplots(1)
    metadata = combined_adata.obs
    data = combined_adata.to_df()

    groups, counts = np.unique(metadata[condition_key], return_counts=True)
    for i, c in enumerate(counts):
        ax.bar(i, c, color=sample_cmap[groups[i]])

    ax.set_xticks(np.arange(i + 1))
    ax.set_xticklabels(groups)
    ax.set_ylabel("Cells Counts Per Condition")

    fig.tight_layout()
    temp_figname = "meld_condition_counts_{}_{}_{}.png".format(
        condition_key, condition1, condition2
    )
    fig.savefig(FIGDIR + temp_figname, dpi=200)
    logger.info("saved to: {}".format(temp_figname))

    # Meld Run Phate
    logger.info("MELD: Run Phate")

    # TODO: make this an input to the script
    data_libnorm, libsize = scprep.normalize.library_size_normalize(
        data, return_library_size=True
    )
    metadata["library_size"] = libsize

    # Two notes: we already transformed the data, so don't need to take sqrt here
    # CYTOF only has 42 channels so no dimensionality reduction prior to this

    phate_op = phate.PHATE(knn=10, decay=10, n_jobs=-1)
    data_phate = phate_op.fit_transform(data)

    fig = scprep.plot.scatter2d(
        data_phate,
        c=metadata[condition_key],
        cmap=sample_cmap,
        legend_anchor=(1, 1),
        figsize=(6, 5),
        dpi=200,
        s=10,
        filename="meld_phate_{}_{}_{}.png".format(
            condition_key, condition1, condition2
        ),
        label_prefix="PHATE",
        ticks=False,
    )

    # 3D PHATE components are used to create the ground truth PDF
    # benchmarker = meld.Benchmarker()
    # benchmarker.fit_phate(data)

    # Run Meld
    logger.info("MELD: Run MELD")
    metadata["condition"] = [
        1 if sl.startswith(condition2) else 0 for sl in metadata[condition_key]
    ]
    metadata["condition_name"] = [
        "Cond 1" if g == 0 else "Cond 2" for g in metadata["condition"]
    ]
    metadata["replicate"] = metadata[condition_key]
    # TODO: Implement benchmarking to find best parameter combinations
    # TODO: look into MELD's failing unit tests

    meld_op = meld.MELD(beta=67, knn=7)
    sample_densities = meld_op.fit_transform(
        data, sample_labels=metadata[condition_key]
    )
    # Normalize across conditions.
    logger.info("MELD: Normalizing")
    sample_likelihoods = sklearn.preprocessing.normalize(sample_densities, norm="l1")

    # Add to metadata
    fig, axes = plt.subplots(1, 1, figsize=(13, 4))
    experimental_samples = [condition1, condition2]
    metadata["BCG_likelihood"] = sample_likelihoods[:, 1]

    scprep.plot.scatter2d(
        data_phate,
        c=metadata["BCG_likelihood"],
        cmap=meld.get_meld_cmap(),
        vmin=0.2,
        vmax=0.75,
        title="Likelihood of Existing in BCG Refractory State",
        ticks=False,
        ax=axes,
    )

    fig.tight_layout()
    fig.savefig(
        FIGDIR
        + "meld_phate_vfc__{}_{}_{}.png".format(condition_key, condition1, condition2),
        dpi=200,
    )
    fig = _plot_jitter_by_cluster(
        metadata, sample_cmap, cluster_key=cluster_key, condition_key=condition_key
    )
    fig.savefig(
        FIGDIR
        + "meld_jitter_vfc__{}_{}_{}.png".format(condition_key, condition1, condition2),
        dpi=200,
    )
    return metadata


if __name__ == "__main__":
    os.chdir(
        "/sc/arion/projects/nmibc_bcg/CancerDiscovery/data/CancerDiscovery_BCG_CyTOF"
    )
    # Listing and importing the FCS files; extracting their condition as well
    fcs_list = []
    fcs_conditions = []
    directory = os.fsencode(".")
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".fcs"):
            fcs_list.append(FCMeasurement(ID=file, datafile=filename))
            condition = ""
            if "K562" in filename:
                condition = "K562"
            elif "RAJI" in filename:
                condition = "RAJI"
            elif "MED" in filename:
                condition = "MED"
            fcs_conditions.append(condition)

    # Too many channels - limiting to those conjugated to metals we want
    channels_of_interest = [
        "111Cd_LILRB1",
        "112Cd_Granzyme_A",
        "113Cd_CD38",
        "114Cd_CD11b",
        "115In_IFNg",
        "116Cd_CD57",
        # "141Pr_cKit",
        "142Nd_KLRG1",
        "143Nd_Granzyme_K",
        "144Nd_CD69",
        "145Nd_NKG2D_",
        "146Nd_DNAM1",
        "147Sm_NKp80",
        #"148Nd_KIR3DL1_L2",
        #"149Sm_KIR3DL1",
        "150Nd_IL2",
        "151Eu_CD107a",
        "152Sm_TNFa",
        "153Eu_CD14",
        "154Sm_MIP1b",
        "155Gd_NKp46",
        "156Gd_Tim3",
        #"158Gd_KIR2DL1",
        "159Tb_CD56",
        "160Gd_NKG2A",
        "161Dy_NKp44",
        "162Dy_CD27",
        "163Dy_Eomes",
        "164Dy_NKG2C",
        #"165Ho_KIR2DL3",
        #"166Er_KIR2DL1_S1",
        #"167Er_KIR2DL2_L3",
        "168Er_TRAIL",
        #"169Tm_Tbet",
        "170Er_CD3",
        "171Yb_CD127",
        "172Yb_Perforin",
        "173Yb_Granzyme_B",
        "174Yb_TIGIT",
        "175Lu_XCL1",
        "176Yb_CD200R",
        "209Bi_CD16",
    ]

    # Import Berengere's data
    clinical = pd.read_excel("id_list_cytof_blca_bcg.xlsx")
    clinical["timepoint"] = clinical["sample_type"].replace(
        {"I1D1": "Pre-BCG", "at recurrence": "Recurrence", "recurrence": "Recurrence"}
    )

    # Debugging - need to remove one non-specific AHBCG17 
    clinical.drop(53, inplace=True)

    # DEBUG: subsampling to make the process go faster
    # logger.info("Subsampling!!! Take me out!")
    # subsample_index = np.random.choice(adata.shape[0], size=5000, replace=False)
    # adata = adata[subsample_index].copy()
    # fcs_list = fcs_list[0:30]
    # fcs_conditions = fcs_conditions[0:30]

    # PRE VS POST
    filter_dict = {
        "condition_column": "condition",
        "timepoint_column": "timepoint",
        "conditions_of_interest": ["MED"],
        "timepoints_of_interest": ["Pre-BCG", "Recurrence"],
    }
    mtx_prepost = cytof_tools.preprocess_cytof(
        fcs_list,
        fcs_conditions,
        clinical,
        filter_dict,
        batch_key="sampleID",
        channels=channels_of_interest,
    )
    # RECURRENCE MED VS RAJI
    filter_dict = {
        "condition_column": "condition",
        "timepoint_column": "timepoint",
        "conditions_of_interest": ["MED", "RAJI"],
        "timepoints_of_interest": ["Recurrence"],
    }
    mtx_raji = cytof_tools.preprocess_cytof(
        fcs_list,
        fcs_conditions,
        clinical,
        filter_dict,
        batch_key="sampleID",
        channels=channels_of_interest,
    )
    # RECURRENCE MED VS RAJI
    filter_dict = {
        "condition_column": "condition",
        "timepoint_column": "timepoint",
        "conditions_of_interest": ["MED", "K562"],
        "timepoints_of_interest": ["Recurrence"],
    }
    mtx_k562 = cytof_tools.preprocess_cytof(
        fcs_list,
        fcs_conditions,
        clinical,
        filter_dict,
        batch_key="sampleID",
        channels=channels_of_interest,
    )

    # Running Phenograph

    meld_dict = {
        "Pre vs Post": {
            "condition1": "Pre-BCG",
            "condition2": "Recurrence",
            "condition_key": "timepoint",
        },
        "Recurrence Med vs Raji": {
            "condition1": "MED",
            "condition2": "RAJI",
            "condition_key": "condition",
        },
        "Recurrence Med vs K562": {
            "condition1": "MED",
            "condition2": "K562",
            "condition_key": "condition",
        },
    }

    for iteration, adata in zip(
        ["Recurrence Med vs Raji", "Recurrence Med vs K562", "Pre vs Post",],
        [mtx_raji, mtx_k562, mtx_prepost],
    ):
        logger.info("{} Iteration has begun".format(iteration))
        condition_key = meld_dict[iteration]["condition_key"]
        condition1 = meld_dict[iteration]["condition1"]
        condition2 = meld_dict[iteration]["condition2"]
        logger.info("Phenograph")

        k = 30
        sc.tl.pca(adata, n_comps=10)
        communities, graph, Q = sce.tl.phenograph(adata.obsm["X_pca"], k=k)
        adata.obs["PhenoGraph_clusters"] = pd.Categorical(communities)
        adata.uns["PhenoGraph_Q"] = Q
        adata.uns["PhenoGraph_k"] = k

        # Drawing UMAP
        logger.info("Drawing UMAP")
        sc.pp.neighbors(adata, n_neighbors=30, n_pcs=10)
        sc.tl.umap(adata)
        sc.pl.umap(
            adata,
            color=["PhenoGraph_clusters", "sample_type"],
            title="PhenoGraph Assigned Clusters: {}".format(iteration),
            save="phenograph {} 01102022.png".format(iteration),
        )

        # MELD
        logger.info("Running Meld")
        metadata = run_meld_cytof(
            combined_adata=adata,
            condition1=condition1,
            condition2=condition2,
            condition_key=condition_key,
            cluster_key="PhenoGraph_clusters",
        )
        metadata.to_csv(
            "cytof_anndata/cytof_annotated_metadata {} {}.csv".format(condition_key, condition2)
        )
