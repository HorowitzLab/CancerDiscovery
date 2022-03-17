# CancerDiscovery

Code for all cancer discovery analyses and figures.

The general rule of thumb for the analyses is if the data are too big to practically work with in a notebook, then processing scripts were written to import, preprocess, concatenate, and save a final output file in a workable format. Tools files were written for code that was repetitively run. 

Statistical analyses and figure creation were then run in either R scripts or Jupyter notebooks. Descriptions of the analytical outlines are below. 

## Figure 1: OLINK, HTG, IFN-g Incubation Analyses
- HTG GSEA / Radar Plots
- OLINK longitudinal experiments, validation cohort
- IFN-g Co-Culture Graphs
## Figure 2: Single Cell RNA Seq / Niche Net Analyses
- SC:
    - UMAP clustering (figure panel A)
- NicheNet:
    - Create weighted network
    - Define sender and receiver niches and genes
    - Define target gene set
    - Identify potential ligands
    - Perform ligand activity analysis (figure panel B)
    - Infer receptors and target genes (figure panel C)
    - Generate IFNG-HLA-E network (figure panel D)
- HTG:
    - Plot HLA-E expression by STAT1 tertiles (figure panel E)
## Figure 3: Spatial sequencing
- 
## Figure 4: CYTOF and MELD Analyses
- CYTOF: Follows the steps below
  - Located in the cytof_tools.py script
    - [Normalization](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1917-7#Sec15)
  - Located in the cytof_script.py script
    - [Phenograph for clustering](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1917-7)
    - [MELD](https://www.nature.com/articles/s41587-020-00803-5) 
  - Analysis



