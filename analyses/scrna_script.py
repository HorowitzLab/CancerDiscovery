#!/usr/bin/env python

"""
The initial processing and analysis script for the CYTOF data for the Cancer Discovery BCG Story
It follows the following order
- Import
- Concatenation
- Phenograph
- Differential Expression
TODO: Batch correction? Not sure what batches were run. 
TODO: Ask Amir about any processing steps that were already done in cytobank. Not sure the quality of the data but assuming berengere already filtered and whatnot. 
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

# FOR LOGGING
logger = logging.getLogger("scRNA_script")
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)

# making sure plots & clusters are reproducible
np.random.seed(42)


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
        cmap=sample_cmap,
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

    fig.tight_layout()
    return fig


from joblib import Parallel, delayed


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
    combined_adata, condition1, condition2, cluster_key, condition_key="sample_type"
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
    fig.savefig(temp_figname, dpi=200)
    logger.info("saved to: {}".format(temp_figname))
    # Meld Run Phate
    logger.info("MELD: Run Phate")
    # Reducing dimensions here because CYTOF is already limited
    # TODO: make this an input to the script
    # sqrt function has a similar form as the log function with the added benefit of being stable at 0.
    data_libnorm, libsize = scprep.normalize.library_size_normalize(
        data, return_library_size=True
    )
    metadata["library_size"] = libsize
    data_sqrt = np.sqrt(data_libnorm)
    data_pca = scprep.reduce.pca(data_sqrt)
    phate_op = phate.PHATE(knn=10, decay=10, n_jobs=-1)
    data_phate = phate_op.fit_transform(data_pca)

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
        "Pre-BCG" if g == 0 else "Post-BCG" for g in metadata["condition"]
    ]
    metadata["replicate"] = metadata[condition_key]
    # TODO: Implement benchmarking
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
        "meld_phate_vfc__{}_{}_{}.png".format(condition_key, condition1, condition2),
        dpi=200,
    )
    fig = _plot_jitter_by_cluster(
        metadata, sample_cmap, cluster_key=cluster_key, condition_key=condition_key
    )
    fig.savefig(
        "meld_jitter_vfc__{}_{}_{}.png".format(condition_key, condition1, condition2),
        dpi=200,
    )
    return metadata


if __name__ == "__main__":
    nmibc = ad.read_h5ad("nmibc_sc_combined.h5ad")
    sample_info = pd.read_excel("scrna_samples.xlsx")
    sample_info[sample_info["Filename"].isin(nmibc.obs["orig.ident"].unique())]
    temp = pd.merge(nmibc.obs, sample_info, left_on="sample.id", right_on="Filename")
    temp["index"] = nmibc.obs.index
    nmibc.obs = temp.set_index("index")
    # Filtering to just BCG and treatment Naive
    nmibc = nmibc[nmibc.obs.Treatment.isin(["BCG", "none"])]
    # subsample for debugging
    # subsample_index = np.random.choice(nmibc.shape[0], size=2000, replace=False)
    # nmibc = nmibc[subsample_index].copy()

    # Running Phenograph
    logger.info("Running Phenograph")
    k = 30
    sc.tl.pca(nmibc, n_comps=10)
    communities, graph, Q = sce.tl.phenograph(nmibc.obsm["X_pca"], k=k)
    nmibc.obs["PhenoGraph_clusters"] = pd.Categorical(communities)
    nmibc.uns["PhenoGraph_Q"] = Q
    nmibc.uns["PhenoGraph_k"] = k
    # Writing for future analysis
    sc.pp.neighbors(nmibc, n_neighbors=30, n_pcs=10)
    sc.tl.umap(nmibc)
    # Drawing UMAP
    logger.info("Drawing UMAP")

    sc.pl.umap(
        nmibc,
        color=["PhenoGraph_clusters", "Treatment", "cell.type"],
        title="PhenoGraph Assigned Clusters: scRNA Data Pre BCG vs Recurrence",
        save="phenograph nmibc 01102022.png",
    )
    # MELD
    logger.info("Running Meld")
    sample = run_meld_cytof(
        combined_adata=nmibc,
        condition1="none",
        condition2="BCG",
        condition_key="Treatment",
        cluster_key="PhenoGraph_clusters",
    )

    nmibc.obs["BCG_likelihood"] = sample["BCG_likelihood"]
    nmibc.obs = nmibc.obs.drop("Unnamed: 13", axis=1)
    nmibc.write_h5ad("test.h5ad")
