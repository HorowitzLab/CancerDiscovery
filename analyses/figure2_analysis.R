# The processing and analysis script for Figure 2 Cancer Discovery BCG NKG2A Story

# __author__ = "Yuanshuo Alice Wang"
# __credits__ = "Daniel Ranti"
# __license__ = "Open Access"
# __version__ = "1.0.1"
# __maintainer__ = "Yuanshuo Alice Wang"
# __email__ = "yuanshuo.wang@gmail.com"
# __status__ = "Peer Review"

# This script follows the following analysis format
# - import SC data
# - create integrated SC object
# - run dimensionality reduction
# - identify UMAP clusters (figure panel A)

library(Seurat)

wd = "/Users/wangy27/SfakianosHorowitzLab/BCG/SingleCell/"
setwd(wd)

# empty list to hold raw SC data
bc.list = list()
# read SC data objects
bc.list[[1]] = Read10X(data.dir="/Users/wangy27/SfakianosHorowitzLab/BLCA/SingleCell/Data/394_BC_157_tumor-CD45pos/outs/filtered_feature_bc_matrix/")
bc.list[[2]] = Read10X(data.dir="/Users/wangy27/SfakianosHorowitzLab/BLCA/SingleCell/Data/435_168_CD45_Tumor/outs/filtered_feature_bc_matrix/")
bc.list[[3]] = Read10X(data.dir="/Users/wangy27/SfakianosHorowitzLab/BLCA/SingleCell/Data/453_BC_170_CD45pos/outs/filtered_feature_bc_matrix/")
bc.list[[4]] = Read10X(data.dir="/Users/wangy27/SfakianosHorowitzLab/BLCA/SingleCell/Data/454_BC_170_CD45neg/outs/filtered_feature_bc_matrix/")
bc.list[[5]] = Read10X(data.dir="/Users/wangy27/SfakianosHorowitzLab/BLCA/SingleCell/Data/BC228_Tumor/")
bc.list[[6]] = Read10X(data.dir="/Users/wangy27/SfakianosHorowitzLab/BLCA/SingleCell/Data/Nmibc-010421/outs/filtered_feature_bc_matrix/")

# empty list to hold SC objects
sc.list = list()
# names of SC samples
sc.names = c("BC_157_postBCG_resis", "BC_168_postBCG_resis", "BC_170_preBCG_refrac_CD45pos", "BC_170_preBCG_refrac_CD45neg", "BC_228_postBCG_resis", "BC_010421_preBCG")
# empty vector to hold variable genes
var.genes = c()
for(i in 1:length(bc.list)) {
  # create Seurat object
  sc.list[[i]] = CreateSeuratObject(counts = bc.list[[i]], project = sc.names[i], min.cells = 3, min.features = 200)
  # add mitochondrial gene percentage
  sc.list[[i]][["percent.mt"]] = PercentageFeatureSet(sc.list[[i]], pattern = "^MT-")
  
  # subset high quality cells
  # filter cells that have unique feature counts over 2500 or less than 200
  # filter cells that have > 15% mitochondrial counts
  sc.list[[i]] = subset(sc.list[[i]], subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<15)
  
  # noramlize data
  sc.list[[i]] = NormalizeData(sc.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  # find top 2000 variable features
  sc.list[[i]] = FindVariableFeatures(sc.list[[i]], selection.method = "vst", nfeatures = 2000)
  # remember variable genes
  var.genes = c(var.genes, VariableFeatures(sc.list[[i]]))
  
}
# assign names to SC objects
names(sc.list) = sc.names
# remove duplicates from list of variable genes
var.genes = unique(var.genes)

# find anchoring features across list of samples
anchors = FindIntegrationAnchors(object.list=sc.list, dims=1:100)
# combine datasets
sc.combined = IntegrateData(anchorset=anchors, dims=1:100, features.to.integrate=var.genes)

# add patient metadata
patient.id = as.character(Idents(object=sc.combined))
names(patient.id) = names(Idents(object=sc.combined))
sc.combined = AddMetaData(object=sc.combined, metadata=patient.id, col.name="pat.id")

# proceed with analysis using integrated dataset
# shift gene expression (mean at 0) and scale expression values (variance across cells = 1)
all.genes = unique(c(rownames(sc.combined), var.genes))
sc.combined = ScaleData(object=sc.combined, features=all.genes)

# linear dimensional reduction via PCA
sc.combined = RunPCA(sc.combined, features=all.genes, npcs=100)
sc.combined = JackStraw(sc.combined, num.replicate = 100, dims=50)
sc.combined = ScoreJackStraw(sc.combined, dims=1:50)

# clustering
sc.combined = FindNeighbors(sc.combined, dims = 1:50)
sc.combined = FindClusters(sc.combined, resolution = 0.9)

# run UMAP non-linear dimensionality reduction
sc.combined = RunUMAP(sc.combined, dims = 1:50)

# find markers for every cluster compared to all remaining cells, report only the positive ones
sc.combined.markers = FindAllMarkers(sc.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay="RNA")

# add cluster cell type information
new.cluster.ids = c("Tumor", # 0
                    "CD8+ T Cell", # 1
                    "DC", # 2
                    "Tumor", # 3
                    "Tumor", # 4
                    "T Cell", # 5
                    "Monocyte/Macrophage", # 6
                    "Granulocyte", # 7
                    "B Cell", # 8
                    "T Cell", # 9
                    "T Reg", # 10
                    "Fibroblast", # 11
                    "Tumor", # 12
                    "Tumor", # 13
                    "NK Cell", # 14
                    "Tumor", # 15
                    "Epithelium", # 16
                    "Tumor", # 17
                    "Tumor", # 18
                    "Granulocyte", # 19
                    "Plasma Cell", # 20
                    "NK Cell", # 21
                    "22", # 22
                    "CD8+ TRM", # 23
                    "24", 
                    "25", 
                    "26", 
                    "Tumor") # 27
names(new.cluster.ids) = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27")
sc.combined = RenameIdents(sc.combined, new.cluster.ids)
# add cell type information
sc.combined@meta.data$cell.type = Idents(sc.combined)

#####
# NicheNet analysis
#####
# This script follows the following analysis format
# - create weighted network
# - define sender and receiver niches and genes
# - define target gene set
# - identify potential ligands
# - perform ligand activity analysis (figure panel B)
# - infer receptors and target genes (figure panel C)

library(nichenetr)
library(tidyverse)

# Nichenet default networks
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
sig_network = readRDS(url("https://zenodo.org/record/3260758/files/signaling_network.rds"))
gr_network = readRDS(url("https://zenodo.org/record/3260758/files/gr_network.rds"))
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

# add IgG and CD16 ligand-receptor interacttion to lr_network
new_lr_network = lr_network
new_lr_network = new_lr_network %>% add_row(tibble_row(from="IGHG1", to="FCGR3A", source="manual", database="manual"))
new_lr_network = new_lr_network %>% add_row(tibble_row(from="IGHG3", to="FCGR3A", source="manual", database="manual"))


# add pathway commons PPI network to the NicheNet model
input_file = "https://maayanlab.cloud/static/hdfs/harmonizome/data/pc/gene_attribute_edges.txt.gz"
ppi_network = read_tsv(input_file, col_names = TRUE)

ppi_network = ppi_network %>% transmute(from=target,to=source) %>% 
  filter(from %in% geneinfo_human$symbol & to %in% geneinfo_human$symbol) # keep only interactions between genes with oficial gene symbols: optional step
ppi_network = ppi_network %>% mutate(source = "pathway_commons_ppi", database = "pathway_commons") 

# add pathway commons PPI to signaling network
new_sig_network = sig_network %>% bind_rows(ppi_network)

# assign strong weight to new PPI network
new_network_weights_df = tibble(source = "pathway_commons_ppi", weight = 1)
new_source_weights_df = optimized_source_weights_df %>% bind_rows(new_network_weights_df)
new_source_weights_df = new_source_weights_df %>% add_row(tibble_row(source="manual", weight=1))

# make new model
# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks = construct_weighted_networks(lr_network = new_lr_network, sig_network = new_sig_network, gr_network = gr_network, source_weights_df = new_source_weights_df)
# downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks, lr_sig_hub = hyperparameter_list$lr_sig_hub, gr_hub = hyperparameter_list$gr_hub)

# Infer ligand-target regulatory potential scores based on the weighted integrated networks
ligands = as.list(unique(new_lr_network$from))
new_ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR", damping_factor = hyperparameter_list$damping_factor, ltf_cutoff = hyperparameter_list$ltf_cutoff)

ligand_target_matrix = new_ligand_target_matrix
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(new_lr_network %>% distinct(from,to), by = c("from","to"))

# perform NicheNet analysis
DefaultAssay(sc.combined) = "RNA"

# 1. Define a “sender/niche” cell population and a “receiver/target” cell population 
# and determine which genes are expressed in both populations
# receiver
receiver = "Tumor"
expressed_genes_receiver = get_expressed_genes(receiver, seurat_obj=sc.combined, pct = 0.10, assay_oi="RNA")
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

# sender
sender_celltypes = c("CD8+ TRM", "CD8+ T Cell", "NK Cell", "T Cell")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_obj=sc.combined, pct=0.10, assay_oi="RNA")
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

# 2. Define a gene set of interest: 
# these are the genes in the “receiver/target” cell population that are potentially affected by 
# ligands expressed by interacting cells (e.g. genes differentially expressed upon cell-cell interaction)
# Hallmark IFNG response genes
ifng.genes = read.delim("/Users/wangy27/SfakianosHorowitzLab/BCG/GSEA/Hallmark/HALLMARK_INTERFERON_GAMMA_RESPONSE.txt", sep="\t", header=TRUE)[-1, 1]
st.ifng.genes = c("CISH", "IFNG", "IFNGR1", "IFNGR2", "JAK1", "JAK2", "PLA2G2A", "PTPRU", "REG1A", "STAT1", "STATIP1")
geneset_oi = unique(c(ifng.genes, st.ifng.genes, "HLA-E", "HLA-F"))

# 3. Define a set of potential ligands: 
# ligands are expressed by the “sender/niche” cell population and bind a (putative) receptor expressed by the “receiver/target” population
ligands = new_lr_network %>% pull(from) %>% unique() #new_lr_network
receptors = new_lr_network %>% pull(to) %>% unique() #new_lr_network

expressed_ligands = intersect(ligands, expressed_genes_sender)
expressed_receptors = c(intersect(receptors, expressed_genes_receiver))

potential_ligands = new_lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = as.list(potential_ligands), algorithm = "PPR", damping_factor = hyperparameter_list$damping_factor, ltf_cutoff = hyperparameter_list$ltf_cutoff)

# 4 Perform NicheNet ligand activity analysis: 
# rank the potential ligands based on the presence of their target genes in the gene set of interest 
# (compared to the background set of genes)
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
# rank predicted ligands by pearson score
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
# take top 5 ligands for downstream analysis
best_upstream_ligands = ligand_activities %>% top_n(5, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

# 5. Infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis
# Active target gene inference
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names()
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names()

# ligand-target links
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

# Find receptors of top-ranked ligands 
lr_network_top = new_lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()


#####
# HTG Data
#####
# This script follows the following analysis format
# - import HTG expression data
# - plot HLA-E expression by STAT1 tertiles (figure panel E)


library(stats)

wd = "/Users/wangy27/SfakianosHorowitzLab/BCG/HTG/"
setwd(wd)

#####
# get median normalized expression matrix
norm.expr.file = "/Users/wangy27/SfakianosHorowitzLab/BCG/HTG/htg_med_norm_expr.txt"
norm.expr = read.delim(norm.expr.file, sep="\t", row.names=1)

colnames(norm.expr) = norm.expr[1, ]
norm.expr = norm.expr[-1:-3, ]

genes = rownames(norm.expr)
norm.expr = as.data.frame(sapply(norm.expr, as.numeric))
rownames(norm.expr) = genes

# remove columns with NA values from expression matrix
rm.ind = which(colSums(is.na(norm.expr)) > 0)
rm.pat = colnames(norm.expr)[which(colSums(is.na(norm.expr)) > 0)]
norm.expr = norm.expr[, -rm.ind]

# get clinical information
clin.file = "/Users/wangy27/SfakianosHorowitzLab/BCG/HTG/HTG_clin_info.txt"
clin.info = read.delim(clin.file, sep="\t", header=TRUE)

# divide patients by STAT1 expression tertiles
stat1.tert = as.numeric(cut(as.numeric(norm.expr["STAT1", ]), 
                            quantile(norm.expr["STAT1", ], probs=seq(0, 1, 1/3), na.rm = FALSE, names = TRUE), 
                            include.lowest = TRUE))

# boxplot of stat1 expression with samples split by STAT1 expression
# prepare data
stat.expr.tile = data.frame(STAT1.expr = as.numeric(norm.expr["STAT1", ]), 
                            STAT1.tert = as.factor(stat1.tert))




#####
# FIGURES
#####
library(ggplot2)
library(pheatmap)
library(viridis)
library(RColorBrewer)
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(as.matrix(xs), probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
# define resolution
ppi = 100

###
# Figure 2A
###
png("Figure2A.png", width=7*ppi, height=6*ppi)
DimPlot(sc.combined, reduction="umap", label=TRUE, cols=c("red", "yellow", "grey70", "green", "grey80", "grey90", "grey80", "#00B8E7", "grey80", "purple", "grey80", "grey80", "grey80", "magenta", "grey80", "grey80", "grey80", "grey80", "red"), pt.size = 0.5, label.size=7) + theme(text=element_text(size=25), axis.text=element_text(size=25)) + NoLegend()
dev.off()

###
# Figure 2B, part 1
###

# ligand rank heatmap
# plot ligand expression in heatmap
rownames(ligand_activities) = ligand_activities$test_ligand
ligands.plot = ligand_activities[1:5, 1]
ligand.rank = data.frame(IFNG = ligand_activities[ligands.plot, 4])
rownames(ligand.rank) = ligands.plot

breakslist=c(0, 0.137354793496066, 0.146623232174545, 0.155851562398324, 0.164959567258, 0.200386247484357, 0.262131603077395, 0.317818114276245, 0.371485010677034)
png("Figure2B_1.png", height=4*ppi, width=3*ppi)
pheatmap(
  mat               = ligand.rank,
  color             = brewer.pal(n=length(breakslist), name="Greens"), 
  border_color      = "black", 
  show_colnames     = TRUE, 
  show_rownames     = TRUE,
  breaks            = breakslist, 
  cluster_cols      = FALSE, 
  cluster_rows      = FALSE, 
  treeheight_row    = 0, 
  drop_levels       = TRUE, 
  legend            = TRUE, 
  scale             = "none", 
  fontsize          = 20,
  cellwidth = 20, 
  cellheight = 20,
  main              = "" 
)
dev.off()

###
# Figure 2B, part 2
###

# ligand expression heatmap
nk.ligand.expr = AverageExpression(sc.combined[, WhichCells(sc.combined, idents="NK Cell")], group.by="BCG", features=ligands.plot, assays="RNA", slot="data")[[1]]
t.ligand.expr = AverageExpression(sc.combined[, WhichCells(sc.combined, idents=c("T Cell"))], group.by="BCG", features=ligands.plot, assays="RNA", slot="data")[[1]]
cd8t.ligand.expr = AverageExpression(sc.combined[, WhichCells(sc.combined, idents="CD8+ T Cell")], group.by="BCG", features=ligands.plot, assays="RNA", slot="data")[[1]]
trm.ligand.expr = AverageExpression(sc.combined[, WhichCells(sc.combined, idents="CD8+ TRM")], group.by="BCG", features=ligands.plot, assays="RNA", slot="data")[[1]]

ligand.expr = cbind(nk.ligand.expr, 
                    t.ligand.expr, 
                    cd8t.ligand.expr, 
                    trm.ligand.expr)

breakslist=c(0, 0.60443287450129, 0.79105887821123, 1.35239836953228, 2.16072228633471, 4.69898326989307, 5.2705699781725, 7.73613399422055, 8.23918221155433)
png("Figure2B_2.png", height=4*ppi, width=3*ppi)
pheatmap(
  mat               = ligand.expr[, c(1, 3, 5, 7)],
  color             = brewer.pal(n=length(breakslist), name="RdPu"), 
  border_color      = "black",
  show_colnames     = TRUE, 
  show_rownames     = TRUE,
  breaks            = breakslist, 
  gaps_col          = c(1, 2, 3), 
  cluster_cols      = FALSE, 
  cluster_rows      = FALSE, 
  treeheight_row    = 0, 
  drop_levels       = TRUE, 
  legend            = TRUE, 
  scale             = "none", 
  fontsize          = 20,
  cellwidth = 20, 
  cellheight = 20,
  main              = "" 
)
dev.off()

###
# Figure 2B, part 3
###

# heatmap of receptor-ligand interaction scores
receptor.ifng.targets = vis_ligand_receptor_network
receptor.ifng.targets = receptor.ifng.targets[, ligands.plot[which(ligands.plot %in% colnames(receptor.ifng.targets))]]

receptors.plot = rownames(receptor.ifng.targets)
receptor.ligand.mat = receptor.ifng.targets

breakslist=c(0, 0.0046268519359581, 0.19414835240735, 0.352970992493745, 0.477439876642825, 0.98057594888268)
png("Figure2B_3.png", height=5*ppi, width=6*ppi)
pheatmap(
  mat               = t(receptor.ligand.mat),
  color             = brewer.pal(n=length(breakslist), name="RdPu"),
  border_color      = "black",
  show_colnames     = TRUE, 
  show_rownames     = TRUE,
  breaks            = breakslist, 
  cluster_cols      = FALSE, 
  cluster_rows      = FALSE, 
  treeheight_row    = 0, 
  drop_levels       = TRUE, 
  legend            = TRUE, 
  scale             = "none", 
  fontsize          = 20,
  cellwidth = 20, 
  cellheight = 20,
  main              = "" 
)
dev.off()

###
# Figure 2B, part 4
###

# heatmap of receptor expression levels in post-BCG tumor
receptor.ifng.targets = vis_ligand_receptor_network
receptor.ifng.targets = receptor.ifng.targets[, ligands.plot[which(ligands.plot %in% colnames(receptor.ifng.targets))]]

receptors.plot = rownames(receptor.ifng.targets)
receptors.expr = AverageExpression(sc.combined[, WhichCells(sc.combined, idents="Tumor")], group.by="BCG", features=receptors.plot, assays="RNA", slot="data")[[1]]

breakslist=c(0, 0.0692335519689722, 0.337384908711783, 0.454806263263544, 0.646188465376249, 0.844763442593128, 1.3755628461118, 2.4255489055009, 13.580701177635)
png("Figure2B_4.png", height=3*ppi, width=6*ppi)
pheatmap( 
  mat               = t(receptors.expr[, c(2, 1)]), 
  color             = brewer.pal(n=length(breakslist), name="Purples"), 
  border_color      = "black", 
  show_colnames     = TRUE, 
  show_rownames     = TRUE,
  breaks            = breakslist, 
  cluster_cols      = FALSE, 
  cluster_rows      = FALSE, 
  treeheight_row    = 0, 
  drop_levels       = TRUE, 
  legend            = TRUE, 
  scale             = "none", 
  fontsize          = 20,
  cellwidth = 20, 
  cellheight = 20,
  main              = "" 
)
dev.off()

###
# Figure 2C
###

library(circlize)

# prepare ligand type matrix
ligand_type_df = tibble(
  ligand_type = c("IFNG-specific", 
                  "TNF-specific",
                  "PTPRC-specific",
                  "TGFB1-specific", 
                  "HMGB1-specific"),
  ligand = c("IFNG", 
             "TNF", 
             "PTPRC", 
             "TGFB1", 
             "HMGB1"))

# Infer target genes of top-ranked ligands and visualize in a circos plot
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = unique(c(ifng.genes, st.ifng.genes)), ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "hallmark_st_ifng") %>% inner_join(ligand_type_df) 

# show only links with a weight higher than a predefined cutoff: 
# links belonging to the 66% of lowest scores were removed
cutoff_include_ligands = active_ligand_target_links_df$weight %>% quantile(0.66)
active_ligand_target_links_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_ligands)

ligands_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_circos$ligand %>% unique())
targets_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_remove &!ligand %in% ligands_remove)

# Prepare the circos visualization: give each segment of ligands and targets a specific color and order
grid_col_ligand =c("IFNG-specific" = "lawngreen",
                   "TNF-specific" = "royalblue",
                   "TGFB1-specific" = "gold", 
                   "HMGB1-specific" = "purple")
grid_col_target =c("hallmark_st_ifng" = "tomato")

col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(col_tbl_ligand) %>% inner_join(col_tbl_target)
links_circ = circos_links %>% dplyr::select(ligand, target, weight)

ligand_col = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_col = ligand_col$color_ligand_type %>% set_names(ligand_col$ligand)
target_col = circos_links %>% distinct(target,color_target_type)
grid_target_col = target_col$color_target_type %>% set_names(target_col$target)

grid_colors = c(grid_ligand_col, grid_target_col)

# Prepare the circos visualization: order ligands and targets
target_ordering = circos_links$target %>% unique()
ligand_ordering = best_upstream_ligands %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_ordering,target_ordering)

# Prepare the circos visualization: define the gaps between the different segments
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

gaps = c(
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "IFNG-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "TNF-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "TGFB1-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "HMGB1-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "hallmark_st_ifng") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
)

# define order of genes in circos plot
order = c("IFNG  ", "TNF  ", "TGFB1  ", "HMGB1  ",
           "CASP1", "CASP3", "CASP8", "RIPK1", "RIPK2", "XAF1", "FAS", "CD40", "TNFAIP3", "TNFSF10", # apoptosis
           "CCL2", "CCL5", "CCL7", "CXCL9", "CXCL10", "CXCL11", # chemokines
           "ICAM1", "VCAM1", "SELP", # adhesion
           "IFNGR1", "SOCS1", "STAT1", "STAT3", "NFKBIA", "IRF1", "IRF7", "IRF8", # IFN-gamma
           "NLRC5", "IL15", "XCL1", # IFN-gamma related
           "CD86", "CIITA", "PSMB8", # antigen presentation
           "CD38", "PTPN1", "PLA2G2A", "IDO1") # enzymes

# Render the circos plot
png("Figure2C.png", width=11.25*ppi, height=11.5*ppi)
circos.par(gap.degree = gaps)
chordDiagram(links_circ, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_colors, transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow", link.visible = links_circ$weight >= 0, annotationTrack = "grid", # link.visible = links_circ$weight >= cutoff_include_ligands
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 2.25)
}, bg.border = NA) 
circos.clear()
dev.off()

###
# Figure 2D -- rendered in Cytoscape
###

###
# Figure 2E
###

# HTG HLA-E expression boxplot w/ STAT1 tertiles
png("Figure2E.png", width=6*ppi, height=6*ppi)
ggplot(stat.expr.tile, aes(x=STAT1.tert, y=STAT1.expr, fill=STAT1.tert)) + 
  geom_boxplot(show.legend = FALSE) + 
  labs(title="STAT1 Expression by STAT1 Tertiles",x="STAT1 Tertiles", y = "STAT1 Expression") +
  scale_fill_brewer(palette="Blues") + 
  theme_classic()
dev.off()

