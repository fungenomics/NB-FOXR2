---
title: "05 - TF fingerprints for progenitor domains"
author: "Selin Jessa [[selin.jessa@mail.mcgill.ca](mailto:selin.jessa@mail.mcgill.ca)] and Bhavyaa Chandarana [[bhavyaa.chandarana@mail.mcgill.ca](mailto:bhavyaa.chandarana@mail.mcgill.ca)]"
date: "04 November, 2024"
output:
  html_document:
    self_contained: yes
    keep_md: yes
    code_folding: show
    theme: flatly
    css: ../include/style.css
    toc: yes
    toc_depth: 4
    number_sections: true
    df_print: paged
    includes:
      before_body: ../include/header.html
      after_body:  ../include/footer.html
---

<!-- FRONT MATTER, insert configuration info -->


<!-- Load custom CSS/JS for code folding -->
<link rel="stylesheet" type="text/css" href="../include/hideOutput.css">
<script src="../include/hideOutput.js"></script>

***

# Configuration

Configuration of project directory & analysis outputs:

<details><summary>Show full config</summary>

```r
source(here("rr_helpers.R"))

# Set up outputs
message("Document index: ", doc_id)
```

```
## Document index: 05
```

```r
# Specify where to save outputs
out        <- here("output", doc_id); dir.create(out, recursive = TRUE)
figout     <- here("figures", doc_id, "/")
cache      <- paste0("/project/kleinman/bhavyaa.chandarana/cache/NB-FOXR2/public/", doc_id, "/")
```

</details>

The root directory of this project is:

```
## /project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public
```

Outputs and figures will be saved at these paths, relative to project root:

```
## public/output/05
```

```
## public/figures/05//
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

Here, we will investigate expression of transcription factors (TFs) which combinatorially
define anatomical progenitor domains in the developing telencephalon.

To demonstrate that 1) these TFs are sufficient to identify derivatives of 
anatomical progenitor domains and 2) that they persist into adulthood, we will
train an SVM model to predict the progenitor domain of origin of single-cell clusters,
using per-cluster expression and detection rate of each TF. We will perform this
using reference datasets from the developing and adult cortex in both mouse and human.

Among these TFs, we will identify _fingerprints_ of different regions, which are expressed in 
a progenitor domain during development, and maintained in the progeny lineages
of that domain through differentiation and adulthood.

Once we have established these TF fingerprints, we will examine their expression
in NB-FOXR2 and other brain tumors, in bulk and scRNAseq data.



# Libraries


```r
library(here)
library(magrittr)
library(tidyr)
library(dplyr)
library(readr)
library(readxl)
library(stringr)
library(glue)
library(purrr)
library(ggplot2)
library(cowplot)
library(data.table)
library(Seurat)
library(tidymodels)
tidymodels_prefer()
library(grid)
library(gridExtra)

source(here("include/style.R"))
source(here("code/functions/RNAseq.R"))
source(here("code/functions/scRNAseq.R"))

ggplot2::theme_set(theme_min())

square_theme <- theme(aspect.ratio = 1)
large_text <- theme(text = element_text(size = 15))
```


# Load metadata

<div class="fold o">


```r
meta <- read_tsv(here("output/00/metadata_patients_NGS.tsv"))
```

```
## Rows: 176 Columns: 20
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (20): Sample, ID_Patient, ID_Sample, ID, FOXR2_positive, Group, Source, ...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
meta_bulk <- meta %>% filter(RNAseq == "Y")
```



```r
meta_sc <- read_tsv(here("output/00/metadata_patients_sc.tsv"))
```

```
## Rows: 16 Columns: 23
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (23): Sample, ID_Patient, ID_Sample, ID, FOXR2_positive, Group, Source, ...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
sc_samples_foxr2 <- meta_sc %>% filter(Group == "NB-FOXR2") %>% pull(ID)
```

</div>


# Helper functions

Functions to create bubbleplots encoding both average gene expression and detection 
rate per cluster, given a Seurat object.


```r
#' Calculate mean expression and detection rate in select
#' clusters for genes of interest
#'
#' @param seurat Seurat object
#' @param genes character, vector of genes for which to calculate mean exp/det. rate
#' @param cluster_col character, name of column in \code{seurat@meta.data} to use
#' as the cluster labels
#' @param clusters character, (Optional) vector of clusters to use (must be a subset of the
#' values in the \code{cluster_col} column of \code{seurat@meta.data})
#' @param scale logical, whether to scale expression for each gene to [0, 1],
#' across clusters
#' @param min_cells numeric, minimum number of cells required in a cluster for it
#' to be included. Clusters with fewer cells will be dropped.
#'
#' @return dataframe with 5 columns: Cluster, Gene, Expression, Pct1, N_cells
extract_meanexp_pct <- function(seurat, genes, cluster_col,
                                clusters = NULL, scale = TRUE,
                                min_cells = 20) {
    
    genes <- genes[genes %in% rownames(seurat)]
    
    Idents(seurat) <- cluster_col
    
    # subset to required clusters
    if (!is.null(clusters)) seurat <- subset(seurat, idents = clusters)
    
    # get average expression
    data_meanexp <- Seurat::AverageExpression(seurat, features = genes, assays = "RNA") %>% 
        .$RNA %>% 
        t()
    
    # re-scale to [0, 1]
    # if (scale == "minmax") data_meanexp <- minmax_norm(data_meanexp)
    if (scale) data_meanexp <- apply(data_meanexp, 2, scales::rescale)
    
    data_meanexp <- data_meanexp %>%
        as.data.frame() %>% 
        tibble::rownames_to_column(var = "Cluster") %>% 
        gather(Gene, Expression, 2:ncol(.))
    
    # get detection rate
    # from functions in code/functions/scRNAseq.R
    data_pct1 <- calc_pct1(seurat, genes) %>%
        data.frame() %>%
        gather(Gene, Pct1, 2:ncol(.))
    
    # the calc_pct1 function converts "-" in gene names to ".", so here we undo that
    data_pct1$Gene <- gsub("\\.", "-", data_pct1$Gene)
    
    # calculate number of cells per cluster
    data_ncells <- data.frame("Cluster" = seurat@meta.data[[cluster_col]]) %>%
        group_by(Cluster) %>%
        summarize(N_cells = n())
    
    # remove small clusters
    clusters_drop <- data_ncells %>% filter(N_cells < min_cells) %>% pull(Cluster)
    message("@ removing ", length(clusters_drop), " clusters with < ", min_cells, " cells")
    
    # put all outputs together into one dataframe
    data_all <- left_join(data_meanexp, data_pct1, by = c("Cluster", "Gene")) %>% 
        mutate(Gene = factor(Gene, levels = rev(genes))) %>% 
        replace_na(list(Expression = 0, Pct1 = 0)) %>% 
        left_join(data_ncells, by = "Cluster") %>% 
        filter(!(Cluster %in% clusters_drop))
    
    return(data_all)
    
    
}

plot_bubble <- function(data_meanexp_pct, cluster_order, 
                        genes, max_radius = 3) {
    
    data_plot <- data_meanexp_pct %>%
        filter(Pct1 > 0) %>%
        mutate(Gene = factor(Gene, levels = rev(genes))) %>%
        mutate(Cluster = factor(Cluster, levels = cluster_order)) %>% 
        filter(!is.na(Cluster))
    
    data_plot %>%
        ggplot(aes(x = Cluster, y = Gene)) +
        geom_point(aes(size = Pct1, colour = Expression), alpha = 0.8) +
        scale_radius(range = c(0, max_radius), limits = c(0, 1)) +
        scale_color_gradientn(colours = ylrd) +
        theme_min2() +
        rotate_x() +
        theme(panel.grid.major.x = element_line(colour = "grey90"),
              panel.border = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_text(size = 13))
    
}
```

# Defining TFs of interest

Define TFs from literature.


```r
# define mouse and human TFs
tfs <- tribble(~Gene_hg, ~Gene_mm,
               "FOXR2",  "Foxr2",
               "FOXG1",  "Foxg1",
               "PAX6",   "Pax6",
               "EMX1",   "Emx1",
               "EMX2",   "Emx2",
               "TBR1",   "Tbr1",
               "EOMES",  "Eomes",
               "GSX2",   "Gsx2",
               "NKX2-1", "Nkx2-1",
               "LHX6",   "Lhx6",
               "DLX1",   "Dlx1",
               "DLX2",   "Dlx2",
               "DLX5",   "Dlx5",
               "DLX6",   "Dlx6")
```



# TF fingerprints in normal development

## Methods

In this section, we load each reference dataset (mouse adult, mouse development,
human adult, and human fetal). For each dataset, cells are clustered and the loaded
objects contain the cluster labels.

Using the function `extract_meanexp_pct()` above in this document, we will:

1. Remove any clusters with less than 20 cells.
2. Calculate the mean expression (i.e., mean of log-transformed expression values in the
`seurat@assay$RNA@data` slot) For each gene, the mean expression is then scaled to [0, 1]
across clusters. The detection rate is already in the range of [0, 1]. The scaling is done
_within_ each dataset, since different datasets may not have directly comparable expression.
3. Calculate the detection rate (aka, "pct1") for each gene in each cluster, representing the proportion of cells in each cluster where each gene is detected. 
4. In the `data wrangling` section, we assign cells into five main classes: 
"CGE/LGE-derived", "MGE-derived", "Excitatory", "Progenitors", "Non-neuroectodermal", "Glia".
5. Clusters in the three classes "CGE/LGE-derived", "MGE-derived", and "Excitatory" are
retained for downstream analysis.

Acknowledging a few implications for the methods:

- small clusters (< 20 cells) are entirely excluded from the analysis; this is mainly to
avoid unstable statistics based on only a handful of cells. Thus, small clusters are
_not_ included in any calculation or visualization.
- normalization of expression is done _within_ each dataset, and _prior_ to filtering out cells in other lineages



## Load datasets

### Human fetal forebrain A (Yu et al, 2021)


```r
load(here("data/singlecell/references_normal/Yu_NatNeurosci_2021/processed_data/seurat_joint.Rda"))
seurat_human_dev_A <- seurat_yu_subpallium
rm(seurat_yu_subpallium)
```


```r
df_human_dev_A <- extract_meanexp_pct(
    seurat_human_dev_A,
    tfs$Gene_hg,
    cluster_col = "celltype",
    scale = TRUE)
```

```
## @ removing 0 clusters with < 20 cells
```

### Human fetal brain B (Shi et al, 2021)


```r
load(here("data/singlecell/references_normal/Shi_Science_2021/seurat_shi.Rda"))
seurat_human_dev_B <- seurat_shi
rm(seurat_shi)

unique(seurat_human_dev_B$Cluster)
```

```
##  [1] "MGE_1"               "MGE_0"               "progenitor_1"       
##  [4] "progenitor_0"        "progenitor_3"        "CGE_0"              
##  [7] "MGE_2"               "Thalamic neurons_0"  "LGE_0"              
## [10] "Excitatory neuron_2" "Endothelial"         "MGE_3"              
## [13] "LGE_2"               "CGE_1"               "progenitor_4"       
## [16] "Excitatory IPC"      "Thalamic neurons_3"  "CGE_2"              
## [19] "Thalamic neurons_1"  "Excitatory neuron_1" "progenitor_2"       
## [22] "Thalamic neurons_5"  "MGE_4"               "LGE_1"              
## [25] "Microglia"           "Excitatory neuron_0" "Thalamic neurons_2" 
## [28] "Excitatory neuron_3" "Thalamic neurons_4"  "OPC_0"              
## [31] "CGE_3"               "CGE_4"               "OPC_1"              
## [34] "LGE_3"               "CGE_6"               "CGE_5"              
## [37] "progenitor_7"        "progenitor_6"        "progenitor_5"
```


```r
df_human_dev_B <- extract_meanexp_pct(
    seurat_human_dev_B,
    tfs$Gene_hg,
    cluster_col = "Cluster",
    scale = TRUE)
```

```
## @ removing 7 clusters with < 20 cells
```

### Mouse embryonal brain (Jessa et al, 2019 & 2022)


```r
# mean expression and pct1 are already pre-computed
atlas_path      <- "/project/kleinman/selin.jessa/from_hydra/atlas/data/"
pct1_feather    <- file.path(atlas_path, "joint_mouse_extended/pct_per_cluster.feather")
meanexp_feather <- file.path(atlas_path,
                             "joint_mouse_extended/mean_expression_per_cluster.feather")
cluster_col     <- "ID_20210710"
metadata_mouse  <- data.table::fread(
    here("data/singlecell/references_normal/Jessa_NatGenet_2022/metadata_20210710_corrected.tsv"),
    data.table = FALSE)

clusters_keep   <- metadata_mouse %>%
    filter(Location == "Forebrain") %>%
    filter(!grepl("EXCLUDE", Label) & Level1_type != "Unresolved" & !grepl("^Y", Exclude)) %>% 
    filter(N_cells >= 20) %>% 
    pull(Label)

x_mean <- feather::read_feather(meanexp_feather, c("ID_20210710", tfs$Gene_mm)) %>%
    rename(Cluster = ID_20210710) %>% 
    tibble::column_to_rownames(var = "Cluster") %>% 
    apply(2, scales::rescale) %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "Cluster") %>%
    gather(Gene, Expression, 2:ncol(.)) %>% 
    # correct an issue in the Nkx2-1 name
    mutate(Gene = ifelse(Gene == "Nkx2.1", "Nkx2-1", Gene))

x_pct1 <- feather::read_feather(pct1_feather, c("ID_20210710", tfs$Gene_mm)) %>% 
    rename(Cluster = ID_20210710) %>% 
    gather(Gene, Pct1, 2:ncol(.))

df_mouse_dev <- left_join(x_mean, x_pct1, by = c("Cluster", "Gene")) %>%
    filter(Cluster %in% clusters_keep) %>%
    mutate(Gene = factor(Gene, levels = rev(tfs$Gene_mm))) %>% 
    replace_na(list(Expression = 0, Pct1 = 0)) %>% 
    left_join(metadata_mouse %>% select(Label, N_cells),
              by = c("Cluster" = "Label"))
```

### Human adult cortex (Hodge et al, 2019 - Allen human)


```r
load(here("data/singlecell/references_normal/Hodge_Nature_2019__AllenBrainAtlas_human_cortex_SmartSeq/seurat_hodge.Rda"))
seurat_human_adult <- seurat_hodge_cortex
rm(seurat_hodge_cortex)

meta_human_adult <- seurat_human_adult@meta.data
saveRDS(meta_human_adult, file = glue("{out}/cell_metadata_human_adult.Rds"))
```


```r
unique(seurat_human_adult$class_label)
```

```
## [1] "Exclude"       "GABAergic"     "Glutamatergic" "Non-neuronal"
```

```r
seurat_human_adult_clusters_keep <- seurat_human_adult@meta.data %>% 
    distinct(cluster_label, class_label, subclass_label) %>% 
    filter(class_label != "Exclude") %>% 
    pull(cluster_label)

df_human_adult <- extract_meanexp_pct(
    seurat_human_adult,
    tfs$Gene_hg,
    cluster_col = "cluster_label",
    clusters = seurat_human_adult_clusters_keep,
    scale = TRUE)
```

```
## @ removing 3 clusters with < 20 cells
```

### Mouse adult cortex (Yao et al, 2021 - Allen mouse)


```r
load(here("data/singlecell/references_normal/Yao_Cell_2021__AllenBrainAtlas_mouse_cortex_SmartSeq/seurat_yao.Rda"))
seurat_mouse_adult <- seurat_yao_cortex
rm(seurat_yao_cortex)

meta_mouse_adult <- seurat_mouse_adult@meta.data
saveRDS(meta_mouse_adult, file = glue("{out}/cell_metadata_mouse_adult.Rds"))
```



```r
unique(seurat_mouse_adult$class_label)
```

```
## [1] "Glutamatergic" "GABAergic"     "Non-Neuronal"
```

```r
df_mouse_adult <- extract_meanexp_pct(
    seurat_mouse_adult,
    tfs$Gene_mm,
    cluster_col = "cluster_label",
    scale = TRUE)
```

```
## @ removing 80 clusters with < 20 cells
```


## Data wrangling

Reformat each dataset, adding a harmonized "Class" label.
For datasets where a subclass label is available, we'll also display
that below, to be able to annotate plots.

_**NOTE**_: For some datasets (the adult human & mouse datasets from the Allen Brain Atlas), 
the subclass is defined based on the dendrogram, and imperfect (e.g. with a few populations
from one subtype, clustered in the subtree of another subtype).


```r
# ------------------------------------------------------------------------------
df_human_dev_A_class <- df_human_dev_A %>% 
    mutate(Class = case_when(
        Cluster %in% c("CGE1", "LGE1", "LGE2", "LGE3") ~ "CGE/LGE-derived lineage",
        Cluster %in% c("MGE1", "MGE2") ~ "MGE-derived lineage",
        Cluster %in% c("P1", "P2", "P3", "P4", "P5", "P6") ~ "Progenitors",
        Cluster %in% c("ExN1", "ExN2", "ExN3", "ExN4", "ExN5") ~ "Excitatory lineage",
        Cluster %in% c("EC1", "EC2", "MG") ~ "Non-neuroectodermal",
        Cluster == "OPC" ~ "Glia"))

# ------------------------------------------------------------------------------
df_human_dev_B_class <- df_human_dev_B %>% 
    mutate(Class = case_when(
        grepl("CGE|LGE", Cluster)  ~ "CGE/LGE-derived lineage",
        grepl("MGE", Cluster) ~ "MGE-derived lineage",
        grepl("progenitor", Cluster) ~ "Progenitors",
        grepl("Excitatory", Cluster) ~ "Excitatory lineage",
        grepl("Thalamic neurons", Cluster) ~ "Other neurons",
        grepl("OPC", Cluster) ~ "Glia",
        Cluster %in% c("Endothelial", "Microglia") ~ "Non-neuroectodermal"
    ))

# ------------------------------------------------------------------------------
df_mouse_dev_class <- df_mouse_dev %>% 
    left_join(metadata_mouse %>% select(Cluster = Label, Level3_type, Level1_type), by = "Cluster") %>% 
    mutate(Class = case_when(
        Level3_type == "MGE inhibitory neurons" ~ "MGE-derived lineage",
        Level3_type %in% c("Cortical inhibitory neurons", "Striatal spiny neurons") ~ "CGE/LGE-derived lineage",
        Level3_type == "Cortical excitatory neurons" ~ "Excitatory lineage",
        grepl("RGC|IPC|progenitor", Level3_type) | Level1_type == "Progenitors" ~ "Progenitors",
        grepl("[Nn]eurons", Level3_type) ~ "Other neurons",
        Level1_type %in% c("Leptomeningeal", "Blood", "Endothelial", "Vascular", "Immune") ~ "Non-neuroectodermal",
        Level1_type == "Glia" ~ "Glia")) %>% 
    # replace with human gene names for compatibility
    left_join(tfs, by = c("Gene" = "Gene_mm")) %>% 
    select(-Level3_type, -Level1_type, -Gene) %>% 
    rename(Gene = Gene_hg)

# ------------------------------------------------------------------------------
# check the subclass labels
seurat_human_adult@meta.data %>% 
    distinct(cluster_label, subclass_label, class_label) %>% 
    filter(class_label != "Exclude") %>% 
    arrange(class_label, subclass_label, cluster_label) %>% 
    DT::datatable()
```

preserveae8c011d609aa914

```r
df_human_adult_class <- df_human_adult %>% 
    left_join(seurat_human_adult@meta.data %>% distinct(cluster_label, subclass_label, class_label), by = c("Cluster" = "cluster_label")) %>% 
    mutate(Class = case_when(
        class_label == "Glutamatergic" ~ "Excitatory lineage",
        subclass_label %in% c("PVALB", "SST") ~ "MGE-derived lineage",
        subclass_label %in% c("VIP", "LAMP5", "PAX6") ~ "CGE/LGE-derived lineage",
        subclass_label %in% c("Microglia", "Pericyte", "Endothelial", "VLMC") ~ "Non-neuroectodermal",
        subclass_label %in% c("Astrocyte", "OPC", "Oligodendrocyte") ~ "Glia")) %>% 
    dplyr::rename(Subclass = subclass_label) %>% 
    select(-class_label)

# ------------------------------------------------------------------------------
# check the subclass labels
seurat_mouse_adult@meta.data %>% 
    distinct(cluster_label, subclass_label, class_label) %>% 
    arrange(class_label, subclass_label, cluster_label) %>% 
    DT::datatable()
```

preserve79f13230c272a1d3

```r
df_mouse_adult_class <- df_mouse_adult %>% 
    left_join(seurat_mouse_adult@meta.data %>% distinct(cluster_label, subclass_label, class_label), by = c("Cluster" = "cluster_label")) %>% 
    mutate(Class = case_when(
        class_label == "Glutamatergic" ~ "Excitatory lineage",
        subclass_label %in% c("Pvalb", "Sst", "Sst Chodl") ~ "MGE-derived lineage",
        subclass_label %in% c("Vip", "Lamp5", "Sncg") ~ "CGE/LGE-derived lineage",
        subclass_label %in% c("Meis2") ~ "Other neurons",
        subclass_label %in% c("Endo", "Micro-PVM", "VLMC", "SMC-Peri") ~ "Non-neuroectodermal",
        subclass_label %in% c("Oligo", "Astro") ~ "Glia")) %>% 
    # replace with human gene names for compatibility
    left_join(tfs, by = c("Gene" = "Gene_mm")) %>% 
    select(-Gene) %>% 
    rename(Gene = Gene_hg) %>% 
    dplyr::rename(Subclass = subclass_label) %>% 
    select(-class_label)
```

Save intermediates:


```r
save(df_human_dev_A_class, df_human_dev_B_class, df_mouse_dev_class, df_human_adult_class, df_mouse_adult_class,
     file = glue("{out}/training_inputs_perdataset.Rda"))
```

Combine datasets:


```r
# put all datasets together
df_train_all_prefilt <- bind_rows(df_human_dev_A_class %>% mutate(Dataset = "Human dev A"),
                                  df_human_dev_B_class %>% mutate(Dataset = "Human dev B"),
                                  df_mouse_dev_class   %>% mutate(Dataset = "Mouse dev"),
                                  df_human_adult_class %>% mutate(Dataset = "Human adult"),
                                  df_mouse_adult_class %>% mutate(Dataset = "Mouse adult"))

# save a map from cluster ---> class label to be reused later
cluster_class_map <- df_train_all_prefilt %>% distinct(Cluster, Class, Dataset) %>%
    mutate(Class = as.character(Class)) %T>%
    rr_write_tsv(glue("{out}/cluster_class_map.tsv"), "Map from normal brain clusters to broad classes used in analysis") %>%
    select(Cluster, Class) %>%
    tibble::deframe()
```

```
## ...writing description of cluster_class_map.tsv to public/output/05/cluster_class_map.desc
```

```r
df_train_all <- df_train_all_prefilt %>%
    filter(Class %in% c("MGE-derived lineage", "CGE/LGE-derived lineage", "Excitatory lineage")) %>% 
    mutate(Class = factor(Class, levels = c("Excitatory lineage",
                                            "CGE/LGE-derived lineage",
                                            "MGE-derived lineage"))) %>% 
    relocate(Class, .after = "Cluster")

# convert to wide format, so that each sample (datapoint) is in one row, and each
# feature/variable is in a separate column
df_train_wide <- df_train_all %>%
    # remove the column containing # of cells per cluster, and remove the subclass label
    select(-N_cells, -Subclass) %>% 
    # exclude Foxr2 for the SVM portion of the analysis
    filter(Gene != "FOXR2") %>%
    # create two features per gene, Expression and detection rate (Pct1)
    pivot_wider(names_from = "Gene", values_from = c("Expression", "Pct1")) %>% 
    mutate(Dataset = case_when(
        Dataset %in% c("Human dev A", "Human dev B") ~ "Human dev",
        TRUE ~ Dataset
    ))

# sanity check - confirm there's no NAs in the Class labels
nrow(df_train_wide[is.na(df_train_wide$Class), ]) == 0
```

```
## [1] TRUE
```

```r
# sanity check - confirm cluster/class correspondences make sense
df_train_all %>% distinct(Dataset, Cluster, Class, N_cells) %>% arrange(Class)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Dataset"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Cluster"],"name":[2],"type":["chr"],"align":["left"]},{"label":["Class"],"name":[3],"type":["fct"],"align":["left"]},{"label":["N_cells"],"name":[4],"type":["int"],"align":["right"]}],"data":[{"1":"Human dev A","2":"ExN1","3":"Excitatory lineage","4":"2972"},{"1":"Human dev A","2":"ExN3","3":"Excitatory lineage","4":"1698"},{"1":"Human dev A","2":"ExN2","3":"Excitatory lineage","4":"2125"},{"1":"Human dev A","2":"ExN4","3":"Excitatory lineage","4":"1204"},{"1":"Human dev A","2":"ExN5","3":"Excitatory lineage","4":"1195"},{"1":"Human dev B","2":"Excitatory neuron_2","3":"Excitatory lineage","4":"868"},{"1":"Human dev B","2":"Excitatory IPC","3":"Excitatory lineage","4":"1001"},{"1":"Human dev B","2":"Excitatory neuron_1","3":"Excitatory lineage","4":"923"},{"1":"Human dev B","2":"Excitatory neuron_0","3":"Excitatory lineage","4":"944"},{"1":"Human dev B","2":"Excitatory neuron_3","3":"Excitatory lineage","4":"794"},{"1":"Mouse dev","2":"F-e10_CEXN1","3":"Excitatory lineage","4":"261"},{"1":"Mouse dev","2":"F-e12_CEXN1","3":"Excitatory lineage","4":"1199"},{"1":"Mouse dev","2":"F-e13_CEXN1","3":"Excitatory lineage","4":"642"},{"1":"Mouse dev","2":"F-e13_CEXN2","3":"Excitatory lineage","4":"277"},{"1":"Mouse dev","2":"F-e15_CEXN1","3":"Excitatory lineage","4":"1230"},{"1":"Mouse dev","2":"F-e15_CEXN2","3":"Excitatory lineage","4":"552"},{"1":"Mouse dev","2":"F-e16_CEXN1","3":"Excitatory lineage","4":"1138"},{"1":"Mouse dev","2":"F-e16_CEXN2","3":"Excitatory lineage","4":"842"},{"1":"Mouse dev","2":"F-e18_CEXN1","3":"Excitatory lineage","4":"1374"},{"1":"Mouse dev","2":"F-e18_CEXN2","3":"Excitatory lineage","4":"1334"},{"1":"Mouse dev","2":"F-e18_CEXN3","3":"Excitatory lineage","4":"296"},{"1":"Mouse dev","2":"F-p0_CEXN1","3":"Excitatory lineage","4":"803"},{"1":"Mouse dev","2":"F-p0_CEXN2","3":"Excitatory lineage","4":"265"},{"1":"Mouse dev","2":"F-p3_CEXN1","3":"Excitatory lineage","4":"200"},{"1":"Human adult","2":"Exc L2-3 LINC00507 RPL9P17","3":"Excitatory lineage","4":"6507"},{"1":"Human adult","2":"Exc L5-6 THEMIS GPR21","3":"Excitatory lineage","4":"497"},{"1":"Human adult","2":"Exc L5-6 FEZF2 MYBPHL","3":"Excitatory lineage","4":"585"},{"1":"Human adult","2":"Exc L4-5 RORB RPL31P31","3":"Excitatory lineage","4":"2023"},{"1":"Human adult","2":"Exc L4-5 RORB LCN15","3":"Excitatory lineage","4":"1267"},{"1":"Human adult","2":"Exc L6 THEMIS LINC00343","3":"Excitatory lineage","4":"2860"},{"1":"Human adult","2":"Exc L6 FEZF2 FAM95C","3":"Excitatory lineage","4":"1091"},{"1":"Human adult","2":"Exc L4-5 RORB LINC01474","3":"Excitatory lineage","4":"2262"},{"1":"Human adult","2":"Exc L6 FEZF2 KRT17","3":"Excitatory lineage","4":"685"},{"1":"Human adult","2":"Exc L5 RORB SNHG7","3":"Excitatory lineage","4":"410"},{"1":"Human adult","2":"Exc L3-4 RORB SEMA6D","3":"Excitatory lineage","4":"433"},{"1":"Human adult","2":"Exc L3-4 RORB PRSS12","3":"Excitatory lineage","4":"895"},{"1":"Human adult","2":"Exc L5 RORB LINC01202","3":"Excitatory lineage","4":"310"},{"1":"Human adult","2":"Exc L5-6 FEZF2 CYP26B1","3":"Excitatory lineage","4":"103"},{"1":"Human adult","2":"Exc L3-5 RORB CMAHP","3":"Excitatory lineage","4":"626"},{"1":"Human adult","2":"Exc L6 THEMIS C6orf48","3":"Excitatory lineage","4":"99"},{"1":"Human adult","2":"Exc L5-6 THEMIS TMEM233","3":"Excitatory lineage","4":"349"},{"1":"Human adult","2":"Exc L5-6 RORB LINC00320","3":"Excitatory lineage","4":"269"},{"1":"Human adult","2":"Exc L3-4 RORB FOLH1B","3":"Excitatory lineage","4":"364"},{"1":"Human adult","2":"Exc L6 FEZF2 TBC1D26","3":"Excitatory lineage","4":"82"},{"1":"Human adult","2":"Exc L5 FEZF2 SCN7A","3":"Excitatory lineage","4":"44"},{"1":"Human adult","2":"Exc L6 FEZF2 SLITRK6","3":"Excitatory lineage","4":"129"},{"1":"Human adult","2":"Exc L6 FEZF2 P4HA3","3":"Excitatory lineage","4":"40"},{"1":"Human adult","2":"Exc L5-6 FEZF2 ANKRD20A1","3":"Excitatory lineage","4":"273"},{"1":"Human adult","2":"Exc L4-5 RORB HNRNPA1P46","3":"Excitatory lineage","4":"429"},{"1":"Human adult","2":"Exc L6 FEZF2 CPZ","3":"Excitatory lineage","4":"170"},{"1":"Human adult","2":"Exc L3 RORB CARTPT","3":"Excitatory lineage","4":"361"},{"1":"Human adult","2":"Exc L6 FEZF2 ETV4","3":"Excitatory lineage","4":"260"},{"1":"Human adult","2":"Exc L5-6 FEZF2 CABP7","3":"Excitatory lineage","4":"27"},{"1":"Human adult","2":"Exc L5-6 THEMIS THTPA","3":"Excitatory lineage","4":"69"},{"1":"Human adult","2":"Exc L5-6 THEMIS OR1J1","3":"Excitatory lineage","4":"88"},{"1":"Human adult","2":"Exc L3-5 RORB HSPB3","3":"Excitatory lineage","4":"211"},{"1":"Human adult","2":"Exc L3-5 THEMIS ELOF1","3":"Excitatory lineage","4":"94"},{"1":"Human adult","2":"Exc L2-4 RORB GRIK1","3":"Excitatory lineage","4":"1556"},{"1":"Human adult","2":"Exc L3-5 FEZF2 ONECUT1","3":"Excitatory lineage","4":"20"},{"1":"Human adult","2":"Exc L6 FEZF2 TBCC","3":"Excitatory lineage","4":"147"},{"1":"Human adult","2":"Exc L3-5 RORB CD24","3":"Excitatory lineage","4":"101"},{"1":"Human adult","2":"Exc L5-6 THEMIS IL7R","3":"Excitatory lineage","4":"390"},{"1":"Human adult","2":"Exc L4 RORB BHLHE22","3":"Excitatory lineage","4":"1171"},{"1":"Human adult","2":"Exc L4 RORB CACNG5","3":"Excitatory lineage","4":"428"},{"1":"Human adult","2":"Exc L4-5 RORB AIM2","3":"Excitatory lineage","4":"78"},{"1":"Human adult","2":"Exc L4-6 RORB HPCA","3":"Excitatory lineage","4":"240"},{"1":"Human adult","2":"Exc L6 THEMIS EGR3","3":"Excitatory lineage","4":"885"},{"1":"Human adult","2":"Exc L6 FEZF2 VWA2","3":"Excitatory lineage","4":"768"},{"1":"Human adult","2":"Exc L4-5 RORB ASCL1","3":"Excitatory lineage","4":"138"},{"1":"Human adult","2":"Exc L5-6 FEZF2 RSAD2","3":"Excitatory lineage","4":"68"},{"1":"Human adult","2":"Exc L3 LINC00507 PSRC1","3":"Excitatory lineage","4":"360"},{"1":"Human adult","2":"Exc L3-4 RORB RPS3P6","3":"Excitatory lineage","4":"98"},{"1":"Human adult","2":"Exc L5 FEZF2 MORN2","3":"Excitatory lineage","4":"66"},{"1":"Human adult","2":"Exc L3-5 LINC00507 SLN","3":"Excitatory lineage","4":"58"},{"1":"Human adult","2":"Exc L3-5 THEMIS UBE2F","3":"Excitatory lineage","4":"130"},{"1":"Human adult","2":"Exc L3-5 FEZF2 DCN","3":"Excitatory lineage","4":"28"},{"1":"Human adult","2":"Exc L4 RORB CCDC168","3":"Excitatory lineage","4":"256"},{"1":"Human adult","2":"Exc L3 LINC00507 CTXN3","3":"Excitatory lineage","4":"281"},{"1":"Human adult","2":"Exc L3 THEMIS PLA2G7","3":"Excitatory lineage","4":"71"},{"1":"Human adult","2":"Exc L5 FEZF2 DYRK2","3":"Excitatory lineage","4":"35"},{"1":"Mouse adult","2":"253_L5 PT CTX","3":"Excitatory lineage","4":"438"},{"1":"Mouse adult","2":"206_L5 IT CTX","3":"Excitatory lineage","4":"743"},{"1":"Mouse adult","2":"242_L5 PT CTX","3":"Excitatory lineage","4":"202"},{"1":"Mouse adult","2":"203_L5 IT CTX","3":"Excitatory lineage","4":"558"},{"1":"Mouse adult","2":"190_L4/5 IT CTX","3":"Excitatory lineage","4":"269"},{"1":"Mouse adult","2":"252_L5 PT CTX","3":"Excitatory lineage","4":"384"},{"1":"Mouse adult","2":"241_L5 PT CTX","3":"Excitatory lineage","4":"48"},{"1":"Mouse adult","2":"243_L5 PT CTX","3":"Excitatory lineage","4":"39"},{"1":"Mouse adult","2":"207_L5 IT CTX","3":"Excitatory lineage","4":"1008"},{"1":"Mouse adult","2":"192_L4/5 IT CTX","3":"Excitatory lineage","4":"744"},{"1":"Mouse adult","2":"237_L5 PT CTX","3":"Excitatory lineage","4":"144"},{"1":"Mouse adult","2":"189_L4/5 IT CTX","3":"Excitatory lineage","4":"220"},{"1":"Mouse adult","2":"218_L5/6 IT CTX","3":"Excitatory lineage","4":"464"},{"1":"Mouse adult","2":"326_L6 CT CTX","3":"Excitatory lineage","4":"1403"},{"1":"Mouse adult","2":"306_L5 NP CTX","3":"Excitatory lineage","4":"314"},{"1":"Mouse adult","2":"188_L4/5 IT CTX","3":"Excitatory lineage","4":"4273"},{"1":"Mouse adult","2":"182_L2/3 IT CTX","3":"Excitatory lineage","4":"873"},{"1":"Mouse adult","2":"196_L4/5 IT CTX","3":"Excitatory lineage","4":"1610"},{"1":"Mouse adult","2":"304_L5 NP CTX","3":"Excitatory lineage","4":"1148"},{"1":"Mouse adult","2":"226_L6 IT CTX","3":"Excitatory lineage","4":"3020"},{"1":"Mouse adult","2":"221_L6 IT CTX","3":"Excitatory lineage","4":"423"},{"1":"Mouse adult","2":"183_L2/3 IT CTX","3":"Excitatory lineage","4":"3671"},{"1":"Mouse adult","2":"219_L6 IT CTX","3":"Excitatory lineage","4":"47"},{"1":"Mouse adult","2":"260_L6 Car3","3":"Excitatory lineage","4":"1221"},{"1":"Mouse adult","2":"311_L6 NP CT CTX","3":"Excitatory lineage","4":"127"},{"1":"Mouse adult","2":"178_L2/3 IT CTX","3":"Excitatory lineage","4":"407"},{"1":"Mouse adult","2":"323_L6 CT CTX","3":"Excitatory lineage","4":"649"},{"1":"Mouse adult","2":"327_L6 CT CTX","3":"Excitatory lineage","4":"2884"},{"1":"Mouse adult","2":"224_L6 IT CTX","3":"Excitatory lineage","4":"443"},{"1":"Mouse adult","2":"349_L6b CTX","3":"Excitatory lineage","4":"280"},{"1":"Mouse adult","2":"324_L6 CT CTX","3":"Excitatory lineage","4":"154"},{"1":"Mouse adult","2":"227_L6 IT CTX","3":"Excitatory lineage","4":"284"},{"1":"Mouse adult","2":"329_L6 CT CTX","3":"Excitatory lineage","4":"227"},{"1":"Mouse adult","2":"343_L6b CTX","3":"Excitatory lineage","4":"266"},{"1":"Mouse adult","2":"231_L6 IT CTX","3":"Excitatory lineage","4":"45"},{"1":"Mouse adult","2":"200_L4/5 IT CTX","3":"Excitatory lineage","4":"1512"},{"1":"Mouse adult","2":"205_L5 IT CTX","3":"Excitatory lineage","4":"537"},{"1":"Mouse adult","2":"195_L4/5 IT CTX","3":"Excitatory lineage","4":"1089"},{"1":"Mouse adult","2":"179_L2/3 IT CTX","3":"Excitatory lineage","4":"216"},{"1":"Mouse adult","2":"180_L2/3 IT CTX","3":"Excitatory lineage","4":"163"},{"1":"Mouse adult","2":"339_L6b CTX","3":"Excitatory lineage","4":"521"},{"1":"Mouse adult","2":"305_L5 NP CTX","3":"Excitatory lineage","4":"586"},{"1":"Mouse adult","2":"191_L4/5 IT CTX","3":"Excitatory lineage","4":"630"},{"1":"Mouse adult","2":"258_L6 Car3","3":"Excitatory lineage","4":"201"},{"1":"Mouse adult","2":"351_L6b CTX","3":"Excitatory lineage","4":"73"},{"1":"Mouse adult","2":"209_L5 IT TPE-ENT","3":"Excitatory lineage","4":"37"},{"1":"Mouse adult","2":"328_L6 CT CTX","3":"Excitatory lineage","4":"93"},{"1":"Mouse adult","2":"238_L5 PT CTX","3":"Excitatory lineage","4":"45"},{"1":"Mouse adult","2":"348_L6b CTX","3":"Excitatory lineage","4":"121"},{"1":"Mouse adult","2":"157_L2 IT HATA","3":"Excitatory lineage","4":"108"},{"1":"Mouse adult","2":"250_L5 PT CTX","3":"Excitatory lineage","4":"236"},{"1":"Mouse adult","2":"247_L5 PT CTX","3":"Excitatory lineage","4":"90"},{"1":"Mouse adult","2":"244_L5 PT CTX","3":"Excitatory lineage","4":"44"},{"1":"Mouse adult","2":"307_L5 NP CTX","3":"Excitatory lineage","4":"69"},{"1":"Mouse adult","2":"340_L6b CTX","3":"Excitatory lineage","4":"92"},{"1":"Mouse adult","2":"347_L6b CTX","3":"Excitatory lineage","4":"49"},{"1":"Mouse adult","2":"176_L2/3 IT CTX","3":"Excitatory lineage","4":"181"},{"1":"Mouse adult","2":"246_L5 PT CTX","3":"Excitatory lineage","4":"94"},{"1":"Mouse adult","2":"217_L5/6 IT CTX","3":"Excitatory lineage","4":"44"},{"1":"Mouse adult","2":"341_L6b CTX","3":"Excitatory lineage","4":"339"},{"1":"Mouse adult","2":"177_L2/3 IT CTX","3":"Excitatory lineage","4":"45"},{"1":"Mouse adult","2":"251_L5 PT CTX","3":"Excitatory lineage","4":"24"},{"1":"Mouse adult","2":"212_L5 IT TPE-ENT","3":"Excitatory lineage","4":"174"},{"1":"Mouse adult","2":"239_L5 PT CTX","3":"Excitatory lineage","4":"28"},{"1":"Mouse adult","2":"254_L5 PT CTX","3":"Excitatory lineage","4":"48"},{"1":"Mouse adult","2":"259_L6 Car3","3":"Excitatory lineage","4":"517"},{"1":"Mouse adult","2":"211_L5 IT TPE-ENT","3":"Excitatory lineage","4":"61"},{"1":"Mouse adult","2":"236_L3 RSP-ACA","3":"Excitatory lineage","4":"188"},{"1":"Mouse adult","2":"350_L6b CTX","3":"Excitatory lineage","4":"403"},{"1":"Mouse adult","2":"322_L6 CT CTX","3":"Excitatory lineage","4":"227"},{"1":"Mouse adult","2":"321_L6 CT CTX","3":"Excitatory lineage","4":"65"},{"1":"Mouse adult","2":"222_L6 IT CTX","3":"Excitatory lineage","4":"94"},{"1":"Mouse adult","2":"325_L6 CT CTX","3":"Excitatory lineage","4":"45"},{"1":"Mouse adult","2":"330_L6 CT CTX","3":"Excitatory lineage","4":"385"},{"1":"Mouse adult","2":"220_L6 IT CTX","3":"Excitatory lineage","4":"25"},{"1":"Mouse adult","2":"308_L5 NP CTX","3":"Excitatory lineage","4":"80"},{"1":"Mouse adult","2":"131_L2 IT RSPv","3":"Excitatory lineage","4":"511"},{"1":"Mouse adult","2":"197_L4/5 IT CTX","3":"Excitatory lineage","4":"774"},{"1":"Mouse adult","2":"309_L5 NP CTX","3":"Excitatory lineage","4":"31"},{"1":"Mouse adult","2":"255_L5 PT RSP-ACA","3":"Excitatory lineage","4":"78"},{"1":"Mouse adult","2":"208_L5 IT CTX","3":"Excitatory lineage","4":"83"},{"1":"Mouse adult","2":"198_L4/5 IT CTX","3":"Excitatory lineage","4":"52"},{"1":"Mouse adult","2":"320_L6 CT CTX","3":"Excitatory lineage","4":"51"},{"1":"Mouse adult","2":"201_L4/5 IT TPE-ENT","3":"Excitatory lineage","4":"37"},{"1":"Mouse adult","2":"172_L2/3 IT CTX","3":"Excitatory lineage","4":"82"},{"1":"Mouse adult","2":"215_L5/6 IT CTX","3":"Excitatory lineage","4":"65"},{"1":"Mouse adult","2":"184_L2/3 IT TPE","3":"Excitatory lineage","4":"29"},{"1":"Mouse adult","2":"262_L6 Car3","3":"Excitatory lineage","4":"39"},{"1":"Mouse adult","2":"173_L2/3 IT ENTl","3":"Excitatory lineage","4":"73"},{"1":"Mouse adult","2":"171_L2/3 IT CTX","3":"Excitatory lineage","4":"93"},{"1":"Mouse adult","2":"193_L4/5 IT CTX","3":"Excitatory lineage","4":"288"},{"1":"Mouse adult","2":"175_L2/3 IT CTX","3":"Excitatory lineage","4":"22"},{"1":"Mouse adult","2":"334_L6b/CT ENT","3":"Excitatory lineage","4":"228"},{"1":"Mouse adult","2":"166_L2/3 IT ENTl","3":"Excitatory lineage","4":"30"},{"1":"Mouse adult","2":"142_L3 IT ENTm","3":"Excitatory lineage","4":"77"},{"1":"Mouse adult","2":"170_L2/3 IT CTX","3":"Excitatory lineage","4":"90"},{"1":"Mouse adult","2":"337_L6b/CT ENT","3":"Excitatory lineage","4":"43"},{"1":"Mouse adult","2":"146_L3 IT ENTl","3":"Excitatory lineage","4":"65"},{"1":"Mouse adult","2":"332_L6 CT ENTm","3":"Excitatory lineage","4":"108"},{"1":"Mouse adult","2":"331_L6 CT ENTm","3":"Excitatory lineage","4":"275"},{"1":"Mouse adult","2":"336_L6b/CT ENT","3":"Excitatory lineage","4":"22"},{"1":"Mouse adult","2":"139_L2 IT ENTl","3":"Excitatory lineage","4":"166"},{"1":"Mouse adult","2":"344_L6b RHP","3":"Excitatory lineage","4":"24"},{"1":"Mouse adult","2":"143_L3 IT ENTm","3":"Excitatory lineage","4":"228"},{"1":"Mouse adult","2":"148_L2 IT PAR","3":"Excitatory lineage","4":"67"},{"1":"Mouse adult","2":"144_L3 IT ENTl","3":"Excitatory lineage","4":"25"},{"1":"Mouse adult","2":"141_L3 IT ENTm","3":"Excitatory lineage","4":"107"},{"1":"Mouse adult","2":"140_L3 IT ENTm","3":"Excitatory lineage","4":"74"},{"1":"Mouse adult","2":"159_L2/3 IT ENTl","3":"Excitatory lineage","4":"73"},{"1":"Mouse adult","2":"160_L2/3 IT ENTl","3":"Excitatory lineage","4":"28"},{"1":"Mouse adult","2":"164_L2/3 IT ENTl","3":"Excitatory lineage","4":"32"},{"1":"Mouse adult","2":"162_L2/3 IT ENTl","3":"Excitatory lineage","4":"38"},{"1":"Mouse adult","2":"213_L5 IT TPE-ENT","3":"Excitatory lineage","4":"37"},{"1":"Human dev A","2":"CGE1","3":"CGE/LGE-derived lineage","4":"3513"},{"1":"Human dev A","2":"LGE1","3":"CGE/LGE-derived lineage","4":"4348"},{"1":"Human dev A","2":"LGE2","3":"CGE/LGE-derived lineage","4":"1844"},{"1":"Human dev A","2":"LGE3","3":"CGE/LGE-derived lineage","4":"1669"},{"1":"Human dev B","2":"CGE_0","3":"CGE/LGE-derived lineage","4":"3250"},{"1":"Human dev B","2":"LGE_0","3":"CGE/LGE-derived lineage","4":"4117"},{"1":"Human dev B","2":"LGE_2","3":"CGE/LGE-derived lineage","4":"3049"},{"1":"Human dev B","2":"CGE_1","3":"CGE/LGE-derived lineage","4":"2685"},{"1":"Human dev B","2":"CGE_2","3":"CGE/LGE-derived lineage","4":"523"},{"1":"Human dev B","2":"LGE_1","3":"CGE/LGE-derived lineage","4":"3893"},{"1":"Human dev B","2":"CGE_3","3":"CGE/LGE-derived lineage","4":"224"},{"1":"Mouse dev","2":"F-e10_CINHN","3":"CGE/LGE-derived lineage","4":"168"},{"1":"Mouse dev","2":"F-e12_CINHN","3":"CGE/LGE-derived lineage","4":"808"},{"1":"Mouse dev","2":"F-e13_CINHN","3":"CGE/LGE-derived lineage","4":"263"},{"1":"Mouse dev","2":"F-e15_CINHN","3":"CGE/LGE-derived lineage","4":"595"},{"1":"Mouse dev","2":"F-e16_SMSN","3":"CGE/LGE-derived lineage","4":"493"},{"1":"Mouse dev","2":"F-e16_CINHN","3":"CGE/LGE-derived lineage","4":"319"},{"1":"Mouse dev","2":"F-e18_SMSN","3":"CGE/LGE-derived lineage","4":"631"},{"1":"Mouse dev","2":"F-e18_CINHN","3":"CGE/LGE-derived lineage","4":"458"},{"1":"Mouse dev","2":"F-p0_CINHN","3":"CGE/LGE-derived lineage","4":"254"},{"1":"Mouse dev","2":"F-p3_SMSN","3":"CGE/LGE-derived lineage","4":"138"},{"1":"Human adult","2":"Inh L2-5 VIP TOX2","3":"CGE/LGE-derived lineage","4":"258"},{"1":"Human adult","2":"Inh L1 LAMP5 GGT8P","3":"CGE/LGE-derived lineage","4":"54"},{"1":"Human adult","2":"Inh L1 LAMP5 NDNF","3":"CGE/LGE-derived lineage","4":"704"},{"1":"Human adult","2":"Inh L1-3 VIP ZNF322P1","3":"CGE/LGE-derived lineage","4":"443"},{"1":"Human adult","2":"Inh L3 VIP CBLN1","3":"CGE/LGE-derived lineage","4":"148"},{"1":"Human adult","2":"Inh L1-4 LAMP5 DUSP4","3":"CGE/LGE-derived lineage","4":"860"},{"1":"Human adult","2":"Inh L1 SST CXCL14","3":"CGE/LGE-derived lineage","4":"329"},{"1":"Human adult","2":"Inh L1 PAX6 GRIP2","3":"CGE/LGE-derived lineage","4":"43"},{"1":"Human adult","2":"Inh L1-2 VIP PPAPDC1A","3":"CGE/LGE-derived lineage","4":"34"},{"1":"Human adult","2":"Inh L1 PAX6 CA4","3":"CGE/LGE-derived lineage","4":"197"},{"1":"Human adult","2":"Inh L1 ADARB2 ADAM33","3":"CGE/LGE-derived lineage","4":"428"},{"1":"Human adult","2":"Inh L1-4 VIP CHRNA2","3":"CGE/LGE-derived lineage","4":"113"},{"1":"Human adult","2":"Inh L2-6 VIP VIP","3":"CGE/LGE-derived lineage","4":"252"},{"1":"Human adult","2":"Inh L1-6 LAMP5 CA13","3":"CGE/LGE-derived lineage","4":"315"},{"1":"Human adult","2":"Inh L5-6 LAMP5 SFTA3","3":"CGE/LGE-derived lineage","4":"426"},{"1":"Human adult","2":"Inh L1-5 VIP KCNJ2","3":"CGE/LGE-derived lineage","4":"171"},{"1":"Human adult","2":"Inh L1-3 VIP SSTR1","3":"CGE/LGE-derived lineage","4":"120"},{"1":"Human adult","2":"Inh L1 VIP PRSS8","3":"CGE/LGE-derived lineage","4":"130"},{"1":"Human adult","2":"Inh L1-2 PAX6 SCGN","3":"CGE/LGE-derived lineage","4":"33"},{"1":"Human adult","2":"Inh L6 LAMP5 C1QL2","3":"CGE/LGE-derived lineage","4":"69"},{"1":"Human adult","2":"Inh L1-3 VIP GGH","3":"CGE/LGE-derived lineage","4":"189"},{"1":"Human adult","2":"Inh L2-4 VIP DSEL","3":"CGE/LGE-derived lineage","4":"100"},{"1":"Human adult","2":"Inh L1 VIP SOX11","3":"CGE/LGE-derived lineage","4":"25"},{"1":"Human adult","2":"Inh L1 VIP PCDH20","3":"CGE/LGE-derived lineage","4":"89"},{"1":"Human adult","2":"Inh L3-6 VIP KCTD13","3":"CGE/LGE-derived lineage","4":"88"},{"1":"Human adult","2":"Inh L2-4 VIP LGI2","3":"CGE/LGE-derived lineage","4":"85"},{"1":"Human adult","2":"Inh L1-2 VIP RPL41P3","3":"CGE/LGE-derived lineage","4":"89"},{"1":"Human adult","2":"Inh L1-6 VIP RGS16","3":"CGE/LGE-derived lineage","4":"163"},{"1":"Human adult","2":"Inh L1-3 VIP ACHE","3":"CGE/LGE-derived lineage","4":"97"},{"1":"Human adult","2":"Inh L1 VIP TNFAIP8L3","3":"CGE/LGE-derived lineage","4":"37"},{"1":"Human adult","2":"Inh L1-6 VIP PENK","3":"CGE/LGE-derived lineage","4":"57"},{"1":"Human adult","2":"Inh L1-3 PAX6 NABP1","3":"CGE/LGE-derived lineage","4":"52"},{"1":"Human adult","2":"Inh L1-3 VIP CCDC184","3":"CGE/LGE-derived lineage","4":"25"},{"1":"Human adult","2":"Inh L1-6 VIP RCN1","3":"CGE/LGE-derived lineage","4":"49"},{"1":"Mouse adult","2":"44_Vip","3":"CGE/LGE-derived lineage","4":"605"},{"1":"Mouse adult","2":"9_Lamp5","3":"CGE/LGE-derived lineage","4":"1494"},{"1":"Mouse adult","2":"45_Vip","3":"CGE/LGE-derived lineage","4":"463"},{"1":"Mouse adult","2":"52_Vip","3":"CGE/LGE-derived lineage","4":"1017"},{"1":"Mouse adult","2":"8_Lamp5","3":"CGE/LGE-derived lineage","4":"20"},{"1":"Mouse adult","2":"13_Lamp5","3":"CGE/LGE-derived lineage","4":"176"},{"1":"Mouse adult","2":"56_Vip","3":"CGE/LGE-derived lineage","4":"436"},{"1":"Mouse adult","2":"43_Vip","3":"CGE/LGE-derived lineage","4":"796"},{"1":"Mouse adult","2":"16_Lamp5","3":"CGE/LGE-derived lineage","4":"122"},{"1":"Mouse adult","2":"49_Vip","3":"CGE/LGE-derived lineage","4":"205"},{"1":"Mouse adult","2":"46_Vip","3":"CGE/LGE-derived lineage","4":"496"},{"1":"Mouse adult","2":"54_Vip","3":"CGE/LGE-derived lineage","4":"74"},{"1":"Mouse adult","2":"42_Vip","3":"CGE/LGE-derived lineage","4":"248"},{"1":"Mouse adult","2":"7_Lamp5","3":"CGE/LGE-derived lineage","4":"925"},{"1":"Mouse adult","2":"58_Vip","3":"CGE/LGE-derived lineage","4":"87"},{"1":"Mouse adult","2":"11_Lamp5","3":"CGE/LGE-derived lineage","4":"185"},{"1":"Mouse adult","2":"12_Lamp5","3":"CGE/LGE-derived lineage","4":"145"},{"1":"Mouse adult","2":"50_Vip","3":"CGE/LGE-derived lineage","4":"176"},{"1":"Mouse adult","2":"35_Sncg","3":"CGE/LGE-derived lineage","4":"182"},{"1":"Mouse adult","2":"51_Vip","3":"CGE/LGE-derived lineage","4":"148"},{"1":"Mouse adult","2":"3_Lamp5 Lhx6","3":"CGE/LGE-derived lineage","4":"199"},{"1":"Mouse adult","2":"25_Sncg","3":"CGE/LGE-derived lineage","4":"229"},{"1":"Mouse adult","2":"15_Lamp5","3":"CGE/LGE-derived lineage","4":"198"},{"1":"Mouse adult","2":"59_Vip","3":"CGE/LGE-derived lineage","4":"77"},{"1":"Mouse adult","2":"40_Vip","3":"CGE/LGE-derived lineage","4":"48"},{"1":"Mouse adult","2":"33_Sncg","3":"CGE/LGE-derived lineage","4":"87"},{"1":"Mouse adult","2":"32_Sncg","3":"CGE/LGE-derived lineage","4":"48"},{"1":"Mouse adult","2":"47_Vip","3":"CGE/LGE-derived lineage","4":"329"},{"1":"Mouse adult","2":"48_Vip","3":"CGE/LGE-derived lineage","4":"417"},{"1":"Mouse adult","2":"18_Pax6","3":"CGE/LGE-derived lineage","4":"149"},{"1":"Mouse adult","2":"55_Vip","3":"CGE/LGE-derived lineage","4":"42"},{"1":"Mouse adult","2":"10_Lamp5","3":"CGE/LGE-derived lineage","4":"49"},{"1":"Mouse adult","2":"14_Lamp5","3":"CGE/LGE-derived lineage","4":"48"},{"1":"Mouse adult","2":"41_Vip","3":"CGE/LGE-derived lineage","4":"63"},{"1":"Mouse adult","2":"30_Sncg","3":"CGE/LGE-derived lineage","4":"33"},{"1":"Mouse adult","2":"26_Sncg","3":"CGE/LGE-derived lineage","4":"36"},{"1":"Mouse adult","2":"34_Sncg","3":"CGE/LGE-derived lineage","4":"25"},{"1":"Mouse adult","2":"39_Lamp5","3":"CGE/LGE-derived lineage","4":"25"},{"1":"Mouse adult","2":"31_Sncg","3":"CGE/LGE-derived lineage","4":"90"},{"1":"Mouse adult","2":"27_Sncg","3":"CGE/LGE-derived lineage","4":"21"},{"1":"Human dev A","2":"MGE1","3":"MGE-derived lineage","4":"5420"},{"1":"Human dev A","2":"MGE2","3":"MGE-derived lineage","4":"1838"},{"1":"Human dev B","2":"MGE_1","3":"MGE-derived lineage","4":"5623"},{"1":"Human dev B","2":"MGE_0","3":"MGE-derived lineage","4":"6073"},{"1":"Human dev B","2":"MGE_2","3":"MGE-derived lineage","4":"1931"},{"1":"Human dev B","2":"MGE_3","3":"MGE-derived lineage","4":"611"},{"1":"Human dev B","2":"MGE_4","3":"MGE-derived lineage","4":"199"},{"1":"Mouse dev","2":"F-e10_MGINH","3":"MGE-derived lineage","4":"320"},{"1":"Mouse dev","2":"F-e12_MGINH","3":"MGE-derived lineage","4":"1313"},{"1":"Mouse dev","2":"F-e13_MGINH","3":"MGE-derived lineage","4":"498"},{"1":"Mouse dev","2":"F-e15_MGINH","3":"MGE-derived lineage","4":"1098"},{"1":"Mouse dev","2":"F-e16_MGINH","3":"MGE-derived lineage","4":"388"},{"1":"Mouse dev","2":"F-e18_MGINH","3":"MGE-derived lineage","4":"590"},{"1":"Human adult","2":"Inh L4-6 SST MTHFD2P6","3":"MGE-derived lineage","4":"736"},{"1":"Human adult","2":"Inh L5-6 PVALB FAM150B","3":"MGE-derived lineage","4":"177"},{"1":"Human adult","2":"Inh L5 PVALB CNTNAP3P2","3":"MGE-derived lineage","4":"573"},{"1":"Human adult","2":"Inh L1-3 PVALB WFDC2","3":"MGE-derived lineage","4":"202"},{"1":"Human adult","2":"Inh L5-6 SST TH","3":"MGE-derived lineage","4":"82"},{"1":"Human adult","2":"Inh L3-5 SST MAFB","3":"MGE-derived lineage","4":"1016"},{"1":"Human adult","2":"Inh L5-6 SST ISOC1","3":"MGE-derived lineage","4":"191"},{"1":"Human adult","2":"Inh L5-6 SST KLHL14","3":"MGE-derived lineage","4":"81"},{"1":"Human adult","2":"Inh L5-6 PVALB STON2","3":"MGE-derived lineage","4":"325"},{"1":"Human adult","2":"Inh L6 SST NPY","3":"MGE-derived lineage","4":"35"},{"1":"Human adult","2":"Inh L2-4 PVALB C8orf4","3":"MGE-derived lineage","4":"882"},{"1":"Human adult","2":"Inh L6 LHX6 GLP1R","3":"MGE-derived lineage","4":"82"},{"1":"Human adult","2":"Inh L1-6 PVALB SCUBE3","3":"MGE-derived lineage","4":"93"},{"1":"Human adult","2":"Inh L2-4 SST AHR","3":"MGE-derived lineage","4":"137"},{"1":"Human adult","2":"Inh L3-6 PVALB MFI2","3":"MGE-derived lineage","4":"66"},{"1":"Human adult","2":"Inh L3-4 PVALB HOMER3","3":"MGE-derived lineage","4":"379"},{"1":"Human adult","2":"Inh L1-2 PVALB TAC1","3":"MGE-derived lineage","4":"26"},{"1":"Human adult","2":"Inh L4-5 PVALB TRIM67","3":"MGE-derived lineage","4":"83"},{"1":"Mouse adult","2":"113_Pvalb","3":"MGE-derived lineage","4":"997"},{"1":"Mouse adult","2":"91_Sst","3":"MGE-derived lineage","4":"173"},{"1":"Mouse adult","2":"61_Sst Chodl","3":"MGE-derived lineage","4":"212"},{"1":"Mouse adult","2":"114_Pvalb","3":"MGE-derived lineage","4":"109"},{"1":"Mouse adult","2":"107_Pvalb","3":"MGE-derived lineage","4":"201"},{"1":"Mouse adult","2":"66_Sst","3":"MGE-derived lineage","4":"177"},{"1":"Mouse adult","2":"115_Pvalb","3":"MGE-derived lineage","4":"1516"},{"1":"Mouse adult","2":"64_Sst","3":"MGE-derived lineage","4":"301"},{"1":"Mouse adult","2":"73_Sst","3":"MGE-derived lineage","4":"338"},{"1":"Mouse adult","2":"105_Pvalb","3":"MGE-derived lineage","4":"141"},{"1":"Mouse adult","2":"112_Pvalb","3":"MGE-derived lineage","4":"224"},{"1":"Mouse adult","2":"82_Sst","3":"MGE-derived lineage","4":"50"},{"1":"Mouse adult","2":"79_Sst","3":"MGE-derived lineage","4":"171"},{"1":"Mouse adult","2":"67_Sst","3":"MGE-derived lineage","4":"1058"},{"1":"Mouse adult","2":"84_Sst","3":"MGE-derived lineage","4":"301"},{"1":"Mouse adult","2":"81_Sst","3":"MGE-derived lineage","4":"38"},{"1":"Mouse adult","2":"87_Sst","3":"MGE-derived lineage","4":"29"},{"1":"Mouse adult","2":"119_Pvalb Vipr2","3":"MGE-derived lineage","4":"180"},{"1":"Mouse adult","2":"90_Sst","3":"MGE-derived lineage","4":"136"},{"1":"Mouse adult","2":"94_Sst","3":"MGE-derived lineage","4":"498"},{"1":"Mouse adult","2":"102_Sst","3":"MGE-derived lineage","4":"130"},{"1":"Mouse adult","2":"95_Sst","3":"MGE-derived lineage","4":"337"},{"1":"Mouse adult","2":"74_Sst","3":"MGE-derived lineage","4":"26"},{"1":"Mouse adult","2":"89_Sst","3":"MGE-derived lineage","4":"153"},{"1":"Mouse adult","2":"116_Pvalb","3":"MGE-derived lineage","4":"217"},{"1":"Mouse adult","2":"85_Sst","3":"MGE-derived lineage","4":"241"},{"1":"Mouse adult","2":"88_Sst","3":"MGE-derived lineage","4":"81"},{"1":"Mouse adult","2":"111_Pvalb","3":"MGE-derived lineage","4":"111"},{"1":"Mouse adult","2":"86_Sst","3":"MGE-derived lineage","4":"62"},{"1":"Mouse adult","2":"80_Sst","3":"MGE-derived lineage","4":"43"},{"1":"Mouse adult","2":"70_Sst","3":"MGE-derived lineage","4":"77"},{"1":"Mouse adult","2":"65_Sst","3":"MGE-derived lineage","4":"128"},{"1":"Mouse adult","2":"106_Pvalb","3":"MGE-derived lineage","4":"42"},{"1":"Mouse adult","2":"109_Pvalb","3":"MGE-derived lineage","4":"90"},{"1":"Mouse adult","2":"110_Pvalb","3":"MGE-derived lineage","4":"47"},{"1":"Mouse adult","2":"93_Sst","3":"MGE-derived lineage","4":"25"},{"1":"Mouse adult","2":"92_Sst","3":"MGE-derived lineage","4":"28"},{"1":"Mouse adult","2":"68_Sst","3":"MGE-derived lineage","4":"64"},{"1":"Mouse adult","2":"62_Sst Chodl","3":"MGE-derived lineage","4":"38"},{"1":"Mouse adult","2":"108_Pvalb","3":"MGE-derived lineage","4":"21"},{"1":"Mouse adult","2":"72_Sst","3":"MGE-derived lineage","4":"29"},{"1":"Mouse adult","2":"101_Sst","3":"MGE-derived lineage","4":"38"},{"1":"Mouse adult","2":"78_Sst","3":"MGE-derived lineage","4":"29"},{"1":"Mouse adult","2":"117_Pvalb","3":"MGE-derived lineage","4":"78"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
save(df_train_all, df_train_wide, cluster_class_map,
     file = glue("{out}/training_data.Rda"))
```

### TABLE: TF quantification per cluster

Export a supplementary table for the expression and detection rate of TF fingerprint
per cluster in each dataset:


```r
TABLE_df_train <- df_train_all %>% 
    filter(Gene != "FOXR2") %>%
    # replace dataset nicknames with refs
    mutate(Dataset = dplyr::recode(Dataset,
                                   "Human dev A" = "Yu et al 2021",
                                   "Human dev B" = "Shi et al 2021",
                                   "Mouse dev"   = "Jessa et al 2019, 2022",
                                   "Human adult" = "Hodge et al 2019",
                                   "Mouse adult" = "Yao et al 2021")) %>% 
    select(-Subclass)

rr_write_tsv(TABLE_df_train, 
             glue("{out}/TABLE_TF_fingerprints.tsv"),
             "Expression and detection of TF fingerprint genes in datasets, per cluster")
```

```
## ...writing description of TABLE_TF_fingerprints.tsv to public/output/05/TABLE_TF_fingerprints.desc
```

## Evaluate class balance

Here, let's check the composition of each class in terms of species/age, 
and let's check the overall balance of the training dataset in terms of class.


```r
# check the class composition
dim(df_train_wide)
```

```
## [1] 363  29
```

```r
p1 <- df_train_wide %>% 
    group_by(Dataset, Class) %>% 
    count() %>% 
    ggplot(aes(x = Class, y = n)) +
    geom_col(aes(fill = Dataset), position = position_dodge()) +
    rotate_x()

p2 <- df_train_wide %>% 
    group_by(Dataset, Class) %>% 
    count() %>% 
    ggplot(aes(x = Dataset, y = n)) +
    geom_col(aes(fill = Class), position = position_dodge()) +
    rotate_x()

plot_grid(p1, p2, nrow = 1, align = "h", axis = "tblr")
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/05//class_balance-1.png)<!-- -->

## Transcription factor-based classifier

Here, we train a linear SVM on the combined input matrix of TF mean expression and detection
rate. In the multi-class scenario, a linear SVM is basically evaluating several one-vs-all
binary classifiers, one per class, i.e. which classifies the samples (in this case, clusters) as belonging
to that class, or not.

### Train models with $k$-fold CV


```r
# set up a linear SVM for classification using tidymodels/parsnip
# https://parsnip.tidymodels.org/articles/Examples.html#svm_linear-models
svm_spec <- svm_linear() %>% set_mode("classification") %>% set_engine("LiblineaR")
svm_spec
```

```
## Linear Support Vector Machine Model Specification (classification)
## 
## Computational engine: LiblineaR
```

```r
df_train_wide_features <- df_train_wide %>% select(-Cluster, -Dataset)

# set up 4-fold validation,
# stratified so that classes are balanced across folds
# https://www.tidymodels.org/start/resampling/#fit-resamples
table(df_train_wide$Class)
```

```
## 
##      Excitatory lineage CGE/LGE-derived lineage     MGE-derived lineage 
##                     193                      95                      75
```

```r
folds <- vfold_cv(df_train_wide, v = 4, strata = "Class")
folds
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["splits"],"name":[1],"type":["list"],"align":["right"]},{"label":["id"],"name":[2],"type":["chr"],"align":["left"]}],"data":[{"1":"<S3: vfold_split>","2":"Fold1"},{"1":"<S3: vfold_split>","2":"Fold2"},{"1":"<S3: vfold_split>","2":"Fold3"},{"1":"<S3: vfold_split>","2":"Fold4"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# loop over folds, holding out one each time for evaluation
fit_vfold <- imap_dfr(folds$splits, function(fold, i) {
    
    # get idx of held-out samples
    out_id <- setdiff(1:nrow(df_train_wide), fold$in_id)
    
    # fit the SVM on the training idx
    fit_fold_i <- fit(svm_spec,
                      formula = Class ~ .,
                      data = df_train_wide_features[fold$in_id, ])
    
    # predict on the held-out fold
    bind_cols(
        df_train_wide[out_id, ] %>% select(Cluster, Class, Dataset),
        predict(fit_fold_i, df_train_wide_features[out_id, ], type = "class")
    ) %>% 
        mutate(Fold = i)
    
})
```

### Evaluate models

Here, we calculate and plot per-fold accuracy & precision evaluation metrics for the SVM.


```r
multi_metric <- yardstick::metric_set(precision, recall)

# calculate metrics per fold, per class
fit_metrics <- fit_vfold %>%
    group_by(Fold, Class) %>%
    multi_metric(truth = Class, estimate = .pred_class)

# calculate metrics per fold, per dataset
fit_metrics2 <- fit_vfold %>%
    group_by(Fold, Class, Dataset) %>%
    multi_metric(truth = Class, estimate = .pred_class)

save(fit_metrics, fit_metrics2, file = glue("{out}/fit_metrics.Rda"))
```

In the below plots, the horizontal bar denotes median value across folds.


```r
palette_metrics <- c("precision" = "dodgerblue4",
                     "recall" = "darkred")

fit_metrics2 %>% 
  filter(Dataset %in% c("Human adult",
                        "Mouse adult")) %>% 
  mutate(Dataset = factor(Dataset,
                            levels = c("Human adult",
                                        "Mouse adult"))) %>% 
  mutate(Class = factor(Class,
                          levels = c("Excitatory lineage",
                                     "CGE/LGE-derived lineage",
                                     "MGE-derived lineage"))) %>% 
  mutate(.metric = factor(.metric,
                          levels = c("precision",
                                     "recall"))) %>% 
  mutate(Metric_class = interaction(.metric,Class)) %>% 
  #unite("Class_Metric", c(Class, .metric),remove = F) %>% 
  ggplot(aes(x = Metric_class, y = .estimate), plot_num = 1) +
    geom_dotplot(binaxis = "y", stackdir = "center", 
                 aes(fill = .metric), dotsize = 1.2, alpha = 0.8, stroke = 0) +
    facet_grid(~ Dataset) +
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5) +
    scale_color_manual(values = palette_metrics) +
    ggtitle("") +
    scale_y_continuous(breaks=c(0, 0.5, 1), 
                       limits=c(0, 1)) +
    geom_vline(xintercept = c(2.5, 4.5), linetype="dashed", color="grey80")+
    rotate_x()
```

```
## Bin width defaults to 1/30 of the range of the data. Pick better value with
## `binwidth`.
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/05//SVM-dot-plot-median-1.png)<!-- -->


## TF fingerprints

To display the TF "fingerprints", we'll first use the previously computed mean expression/detection
rate per cluster and display these values as bubbleplots. Secondly, we'll re-compute
mean expression/detection rate for each gene, but using the _classes_ within each dataset
instead of individual _clusters_, for simpler visualization.

We'll also generate versions with the minimal set of TFs that are sufficient for
classification and persist into adulthood.


```r
# define mouse and human TFs
tfs_min <- tribble(~Gene_hg, ~Gene_mm,
                   "FOXG1",  "Foxg1",
                   "EMX1",   "Emx1",
                   "TBR1",   "Tbr1",
                   "LHX6",   "Lhx6",
                   "DLX1",   "Dlx1",
                   "DLX2",   "Dlx2",
                   "DLX5",   "Dlx5",
                   "DLX6",   "Dlx6")
```

### Per cluster, by lineage

Let's separate out each class, and then plot clusters in order of species/age. (Used in supplementary figures.)

MGE clusters:


```r
bubble_genes <- tfs_min$Gene_hg

# Remove extra TFs not in minimal set from dataframe used for plotting
# Otherwise an extra unnecessary row called "NA" is present.
df_train_all_fig <- df_train_all %>% 
  filter(Gene %in% bubble_genes)

mge_clusters <- df_train_all %>%
    filter(Class == "MGE-derived lineage") %>%
    mutate(Dataset = factor(Dataset, levels = c("Mouse dev", "Human dev A",
                                                "Human dev B", "Mouse adult", "Human adult"))) %>%
    arrange(Dataset, Subclass, Cluster) %>%
    pull(Cluster) %>%
    unique()

df_train_all_fig %>%
    filter(Class == "MGE-derived lineage") %>%
    plot_bubble(genes = bubble_genes,
                cluster_order = mge_clusters) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(size = 6)) +
    ggtitle("MGE-derived lineages")
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/05//fingerprint_mge_bubble-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/05//fingerprint_mge_bubble...*]~</span>

CGE & LGE clusters:


```r
cge_clusters <- df_train_all %>%
    filter(Class == "CGE/LGE-derived lineage") %>%
    mutate(Dataset = factor(Dataset, levels = c("Mouse dev", "Human dev A",
                                                "Human dev B", "Mouse adult", "Human adult"))) %>%
    arrange(Dataset, Subclass, Cluster) %>%
    pull(Cluster) %>%
    unique()

df_train_all_fig %>%
    filter(Class == "CGE/LGE-derived lineage") %>%
    plot_bubble(genes = bubble_genes,
                cluster_order = cge_clusters) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(size = 6)) +
    ggtitle("CGE/LGE-derived lineages")
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/05//fingerprint_cge_bubble-1.png)<!-- -->

Excitatory neuron clusters:


```r
exn_clusters <- df_train_all %>%
    filter(Class == "Excitatory lineage") %>%
    mutate(Dataset = factor(Dataset, levels = c("Mouse dev", "Human dev A",
                                                "Human dev B", "Mouse adult", "Human adult"))) %>%
    arrange(Dataset, Subclass, Cluster) %>%
    pull(Cluster) %>%
    unique()

df_train_all_fig %>%
    filter(Class == "Excitatory lineage") %>%
    plot_bubble(genes = bubble_genes,
                cluster_order = exn_clusters) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(size = 6)) +
    ggtitle("Excitatory lineages")
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/05//fingerprint_exn_bubble-1.png)<!-- -->

### Per aggregated class

Here, we calculate mean expression/detection rate per _class_ instead of per
_cluster_. (Used in main figures.)


```r
# ------------------------------------------------------------------------------
# helper function to make sure each Seurat object has a cluster label
set_class_labels <- function(seurat, cluster_col) {
    
    seurat@meta.data$Class <- plyr::mapvalues(
        as.character(seurat@meta.data[[cluster_col]]),
        from = names(cluster_class_map),
        to = unname(cluster_class_map),
        warn_missing = FALSE)
    
    return(seurat)
    
}

# ------------------------------------------------------------------------------
# do the aggregation for the seurat objects
seurat_human_dev_A <- set_class_labels(seurat_human_dev_A, "celltype")
df_human_dev_A_agg <- extract_meanexp_pct(
    seurat_human_dev_A,
    tfs$Gene_hg,
    cluster_col = "Class",
    scale = TRUE)
```

```
## @ removing 0 clusters with < 20 cells
```

```r
seurat_human_dev_B <- set_class_labels(seurat_human_dev_B, "Cluster")
df_human_dev_B_agg <- extract_meanexp_pct(
    seurat_human_dev_B,
    tfs$Gene_hg,
    cluster_col = "Class",
    scale = TRUE)
```

```
## @ removing 7 clusters with < 20 cells
```

```r
seurat_mouse_adult <- set_class_labels(seurat_mouse_adult, "cluster_label")
df_mouse_adult_agg <- extract_meanexp_pct(
    seurat_mouse_adult,
    tfs$Gene_mm,
    cluster_col = "Class",
    scale = TRUE) %>% 
    left_join(tfs, by = c("Gene" = "Gene_mm")) %>% 
    select(-Gene) %>% 
    rename(Gene = Gene_hg)
```

```
## @ removing 80 clusters with < 20 cells
```

```r
seurat_human_adult <- set_class_labels(seurat_human_adult, "cluster_label")
df_human_adult_agg <- extract_meanexp_pct(
    seurat_human_adult,
    tfs$Gene_hg,
    cluster_col = "Class",
    scale = TRUE)
```

```
## @ removing 10 clusters with < 20 cells
```

The mouse dev brain dataset is not in the form of a Seurat object, so let's 
calculate meanexp/pct1 separately:


```r
# do the aggregation for the mouse dev reference, where we have an indexed
# feather file containing cluster labels and per-cell gene expression values
expr_mouse_dev <- feather::read_feather(
    file.path(atlas_path, "joint_cortex_extended/joint_cortex_extended.embedding_and_genes.feather"),
    columns = c("ID_20210710", tfs$Gene_mm))

expr_mouse_dev_agg <- expr_mouse_dev %>% 
    filter(ID_20210710 %in% clusters_keep) %>% 
    left_join(tibble::enframe(cluster_class_map, "Cluster", "Class"),
              by = c("ID_20210710" = "Cluster")) %>% 
    select(-ID_20210710) %>% 
    relocate(Class, .before = 1)

# calculate mean expression and pct1
meanexpr_mouse_dev_agg <- aggregate(expr_mouse_dev_agg[, 2:ncol(expr_mouse_dev_agg)],
                                    by = list(expr_mouse_dev_agg$Class),
                                    FUN = mean)
names(meanexpr_mouse_dev_agg)[1] <- "Cluster"

# scale to [0, 1] per gene
meanexpr_mouse_dev_agg <- meanexpr_mouse_dev_agg %>% 
    tibble::column_to_rownames("Cluster") %>%  
    apply(2, scales::rescale) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "Cluster")

pct1_mouse_dev_agg <- expr_mouse_dev_agg %>% select(-Class)
# binarize
pct1_mouse_dev_agg[pct1_mouse_dev_agg > 0] <- 1
pct1_mouse_dev_agg <- as.data.table(pct1_mouse_dev_agg)
# Add cluster info for cells
pct1_mouse_dev_agg[, Cluster := as.character(expr_mouse_dev_agg$Class)]
# Get prop of cells expressing a gene, within each cluster
pct1_mouse_dev_agg <- pct1_mouse_dev_agg[, lapply(.SD, prop), by = Cluster] %>% 
    as.data.frame()

# calculate number of cells per cluster
ncells_mouse_dev_agg <- expr_mouse_dev_agg %>% 
    select(Cluster = Class) %>% 
    group_by(Cluster) %>%
    summarize(N_cells = n())

df_mouse_dev_agg <- left_join(
    meanexpr_mouse_dev_agg %>% gather("Gene", "Expression", 2:ncol(.)),
    pct1_mouse_dev_agg %>% gather("Gene", "Pct1", 2:ncol(.)),
    by = c("Cluster", "Gene")) %>% 
    left_join(ncells_mouse_dev_agg) %>% 
    left_join(tfs, by = c("Gene" = "Gene_mm")) %>% 
    select(-Gene) %>% 
    rename(Gene = Gene_hg) 
```

```
## Joining with `by = join_by(Cluster)`
```

```r
save(df_mouse_dev_agg, df_human_dev_A_agg, df_human_dev_B_agg,
     df_mouse_adult_agg, df_human_adult_agg,
     file = glue("{out}/meanexp_pct_aggregated_per_class.Rda"))
```

Combine the aggregated data from all datasets:


```r
df_aggregated <- bind_rows(df_human_dev_A_agg %>% mutate(Dataset = "Human dev A"),
                           df_human_dev_B_agg %>% mutate(Dataset = "Human dev B"),
                           df_mouse_dev_agg   %>% mutate(Dataset = "Mouse dev"),
                           df_human_adult_agg %>% mutate(Dataset = "Human adult"),
                           df_mouse_adult_agg %>% mutate(Dataset = "Mouse adult")) %>%
    filter(Cluster %in% c("MGE-derived lineage", "CGE/LGE-derived lineage", "Excitatory lineage")) %>% 
    mutate(Cluster = factor(Cluster, levels = c("Excitatory lineage",
                                                "CGE/LGE-derived lineage",
                                                "MGE-derived lineage"))) %>% 
    arrange(Cluster) %>% 
    mutate(Cluster = paste0(as.character(Cluster), "\n", Dataset))

agg_order <- df_aggregated %>%  
    pull(Cluster) %>% 
    unique()
```

Minimal set of TFs (used in main figures):


```r
df_aggregated %>%
    filter(Gene %in% tfs_min$Gene_hg) %>% 
    plot_bubble(genes = tfs_min$Gene_hg,
                cluster_order = agg_order,
                max_radius = 8)
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/05//bubbplot_aggregated_by_class_min-1.png)<!-- -->

### TABLE: TF quantification per lineage

Export a supplementary table for the expression and detection rate of TF fingerprint
per lineage in each dataset:


```r
TABLE_df_aggregated <- df_aggregated %>% 
    filter(Gene != "FOXR2") %>%
    # tidy Cluster id
    separate(Cluster, into = c("Cluster", "drop"), sep = "\n") %>% 
    select(-drop) %>% 
    # replace dataset nicknames with refs
    mutate(Dataset = dplyr::recode(Dataset,
                                   "Human dev A" = "Yu et al 2021",
                                   "Human dev B" = "Shi et al 2021",
                                   "Mouse dev"   = "Jessa et al 2019, 2022",
                                   "Human adult" = "Hodge et al 2019",
                                   "Mouse adult" = "Yao et al 2021"))

rr_write_tsv(TABLE_df_aggregated, 
             glue("{out}/TABLE_TF_fingerprints_aggregated.tsv"),
             "Expression and detection of TF fingerprint genes in dev datasets, aggregated by lineage")
```

```
## ...writing description of TABLE_TF_fingerprints_aggregated.tsv to public/output/05/TABLE_TF_fingerprints_aggregated.desc
```



# TF expression in tumor bulk RNAseq

Here, we will inspect expression of the MGE TF fingerprint (developed in above sections) 
as well as key transcription factors and gene markers for extracranial neuroblastoma 
(EC-NB), in bulk RNAseq profiles of human tumors (intracranial pediatric brain tumors as well as EC-NB).


```r
# define mouse and human TFs
tfs <- tribble(~Gene_hg, ~Gene_mm,
               "FOXR2",  "Foxr2",
               "FOXG1",  "Foxg1",
               "PAX6",   "Pax6",
               "EMX1",   "Emx1",
               "EMX2",   "Emx2",
               "TBR1",   "Tbr1",
               "EOMES",  "Eomes",
               "GSX2",   "Gsx2",
               "NKX2-1", "Nkx2-1",
               "LHX6",   "Lhx6",
               "DLX1",   "Dlx1",
               "DLX2",   "Dlx2",
               "DLX5",   "Dlx5",
               "DLX6",   "Dlx6")

# define EC-NB TFs and markers
ecnb_tfs_hg <- c("PHOX2B",
              "HAND1",
              "HAND2",
              "ISL1",
              "GATA3",
              "LMO1",
              "ASCL1",
              "DBH")

# Convert genes from human to mouse using cached biomaRt reference
load(here("data/singlecell/references_genome/biomaRt_mm_to_hg_lds.Rda"))
(ecnb_tfs <- genes_lds %>% filter(HGNC.symbol %in% ecnb_tfs_hg) %>% 
    rename(Gene_mm = MGI.symbol) %>% 
    rename(Gene_hg = HGNC.symbol) %>% 
    arrange(factor(Gene_hg, ecnb_tfs_hg)))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Gene_mm"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Gene_hg"],"name":[2],"type":["chr"],"align":["left"]}],"data":[{"1":"Phox2b","2":"PHOX2B"},{"1":"Hand1","2":"HAND1"},{"1":"Hand2","2":"HAND2"},{"1":"Isl1","2":"ISL1"},{"1":"Gata3","2":"GATA3"},{"1":"Lmo1","2":"LMO1"},{"1":"Ascl1","2":"ASCL1"},{"1":"Dbh","2":"DBH"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# Check that all genes were converted
length(ecnb_tfs$Gene_mm) == length(ecnb_tfs_hg)
```

```
## [1] TRUE
```

```r
# Append ECNB genes to tfs list
tfs_no_ecnb <- tfs
tfs <- rbind(tfs, ecnb_tfs)
```

Define a function to create rows of boxplots:

<details>


```r
# Input:
# - y_maxes: vector of y-axis maximum values, named with the gene symbols to plot
#   --> e.g. c("FOXR2" = 2000, "LHX6" = 4000)
# - scale_to_maxes: T/F, whether to use the y-max values to scale or allow each plot
#   to have a free scale
# - counts: tidy counts matrix which includes columns: gene_symbol, gene_expression, sample, Group
# - palette: palette for Group variable in counts matrix  
# - label_y: OPTIONAL list of indices for box plots that require y axis value labels
#   --> e.g. c(1,2)
#   --> If not provided, only the first box plot will have a label (= c(1))
# - strip_label: OPTIONAL boolean (T/F) whether to add right side strip label for "group"
#   --> default is False (no label)

tf_exp_boxplots <- function(y_maxes,
                            scale_to_maxes = T,
                            counts,
                            palette,
                            label_y = c(1),
                            strip_label = F) {
    
    
    # loop over genes & the maximum y values for each gene,
    # making boxplots for each gene as 1 column to be combined with plot_grid
    # .x = ymaxes (vector values), .y = genes (vector names)
    plots <- imap(y_maxes,
                  function(max, gene_name){
                      
                      p <- counts %>%
                              # subset counts
                              filter(gene_symbol == gene_name) %>%
                                  ggplot(aes(x = factor(1), y = gene_expression)) +
                                  # make boxplot with jittered points
                                  geom_boxplot(aes(fill = Group), outlier.shape = NA, width = 0.3) +
                                  geom_jitter(aes(fill = Group), size = 0.5, width = 0.2, shape = 21) +
                                  scale_fill_manual(values = palette)  +
                                  # make a one-column plot
                                  facet_wrap(~ Group, ncol = 1, strip.position = "right") 
                      
                      if(scale_to_maxes){
                          # subset the plotting area (y-axis) to the specified y-max value
                          # and only show the min and max of the y-axis
                          p <- p + coord_cartesian(ylim = c(0, max)) +
                              scale_y_continuous(breaks = c(0, max))
                      }
                      
                      
                      p <- p + theme_min2() +
                                  no_legend() + ggtitle(gene_name) + xlab(NULL) + ylab(NULL) +
                                  # put the group label on the right instead of the top
                                  theme(strip.text.x = element_text(angle = 0, size = 2)) +
                                  theme(aspect.ratio = 2)
                      
                      return(p)
                      
                  })
    
    if(strip_label){
        # remove the y-axis strip labels in all but the last plot
        plots[1:(length(y_maxes) - 1)] <- map(plots[1:(length(y_maxes) - 1)], ~ .x + theme(strip.text.y = element_blank()))
    } else {
        # remove the y-axis strip labels in ALL plots
        plots[1:length(y_maxes)] <- map(plots[1:length(y_maxes)], ~ .x + theme(strip.text.y = element_blank()))
    }

    # remove the y-axis value labels for all except specified plots
    y_axis_rem <- setdiff(1:length(y_maxes), label_y)
    plots[y_axis_rem] <- map(plots[y_axis_rem], ~ .x + theme(axis.text.y = element_blank()))
    
    plot_grid(plotlist = plots, nrow = 1, align = "h", axis = "tb")
}
```

</details>



## All PBT: bubble plot visualization

Here, we will produce a bubble plot to match those produced above for normal brain, in each pediatric brain tumor type in our cohort.


```r
info_samples <- read_tsv(here("data/RNAseq/pipeline_l3/2023-05-test_pbt/info.samples.tsv"))
```

```
## Rows: 121 Columns: 3
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (3): ID, Nickname, Group
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
counts_lineage <- extract_pipeline_counts(
    path = here("data/RNAseq/pipeline_l3/2023-05-test_pbt/counts/Ensembl.ensGene.exon.norm.tsv.gz"),
    goi = tfs$Gene_hg) %>%
    mutate(gene_symbol = factor(gene_symbol, levels = tfs$Gene_hg)) %>%
    left_join(info_samples, by = c("sample" = "Nickname")) %>%
    mutate(Group = factor(Group, levels = names(palette_groups)))
```

```
## Joining with `by = join_by(gene_symbol)`
```

For this, we will take the median of the DESeq2-normalized expression in each group.


```r
counts_foxr2_medians <- counts_lineage %>%
    group_by(Group, gene_symbol, gene_ensg) %>%
    summarize(median_expr = median(gene_expression)) %>%
    ungroup() %>%
    mutate(log10_median_expr = log10(median_expr))
```

```
## `summarise()` has grouped output by 'Group', 'gene_symbol'. You can override
## using the `.groups` argument.
```

Create a bubbleplot, encoding median expression in each tumor type. Gene expression values where
the log10 expression is <1 are not shown.

Minimal set of TFs: (used in figures)


```r
counts_foxr2_medians %>%
    filter(log10_median_expr >= 1) %>%
    filter(gene_symbol %in% tfs_min$Gene_hg) %>% 
    mutate(gene_symbol = factor(gene_symbol, levels = rev(tfs_min$Gene_hg))) %>%
    ggplot(aes(x = Group, y = gene_symbol), plot_num = 1) +
    geom_point(aes(size = log10_median_expr, colour = log10_median_expr), alpha = 0.8) +
    scale_radius(range = c(0, 10)) +
    scale_color_gradientn(colours = ylrd) +
    theme_min2() +
    rotate_x() +
    theme(panel.grid.major.x = element_line(colour = "grey90"),
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 13)) +
    ggtitle("Normalized expression in \nbulk tumor types")
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/05//bulk_RNAseq_bubble_min-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/05//bulk_RNAseq_bubble_min...*]~</span>


## NB-FOXR2 and DIPG-FOXR2

DIPG-H3K27M tumors with FOXR2 alteration are a useful control to help understand the effect of FOXR2 alteration. For example, if FOXR2 itself activates the MGE TFs, they would also be expressed in DIPG-H3K27M-FOXR2.

Plotting bulk RNA expression of the MGE TF fingerprint in NB-FOXR2 vs. DIPG with FOXR2 alteration.

Since the genes differ in their range of expression, we set a y axis scale of 7000 per gene,
except FOXR2.


```r
y_maxes <- rep(x = 7000,
               times = length(tfs_no_ecnb$Gene_hg))
names(y_maxes) <- tfs_no_ecnb$Gene_hg
y_maxes[which(names(y_maxes) == "FOXR2")] <- 2000

tf_exp_boxplots(y_maxes = y_maxes,
                counts = counts_lineage %>% 
                    filter(Group %in% c("NB-FOXR2",
                                        "DIPG-H3K27M-FOXR2")),
                palette = palette_groups,
                label_y = c(1,2))
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/05//NB-FOXR2_vs_FOXR2_DIPG_7k-1.png)<!-- -->

Print the list of tumors in order, to label the rows of this plot:


```r
names(palette_groups)[c(1,2)]
```

```
## [1] "NB-FOXR2"          "DIPG-H3K27M-FOXR2"
```

## EC-NB 

Here, we will inspect bulk expression of the MGE TF fingerprint in extracranial neuroblastoma (EC-NB), to ensure that they are not activated universally across neuroblastoma, and to ensure the fingerprint is not specifically activated in FOXR2+ EC-NB.

We are using two different datasets of EC-NB to confirm our findings are robust.

1. Gartlgruber et al. *Nature Cancer* 2021
2. TARGET neuroblastoma cohort

### Gartlgruber et al. 2021

Load cohort data:


```r
load(here("output/03/Gartlgruber_et_al_counts.Rda"))
rm(ecnb_counts_norm)
rm(ecnb_counts_vst)

# DESeq norm expression threshold used for FOXR2 +/-
threshold <- 2

foxr2_pos_samples <- ecnb_counts_tidy %>% 
    filter(gene_symbol == "FOXR2") %>% 
    filter(gene_expression > threshold) %>% 
    .$sample

ecnb_counts_tidy <- ecnb_counts_tidy %>% mutate(FOXR2_positive = case_when(sample %in% foxr2_pos_samples ~ "Y",
                                      T ~ "N"))
ecnb_counts_tidy %>% filter(gene_symbol == "FOXR2") %>% select(MYCN, FOXR2_positive) %>% table
```

```
##         FOXR2_positive
## MYCN       N   Y
##   Amp     94  13
##   NonAmp 379  93
```

Set palette & add FOXR2/MYCN status.
Group MYCN-amp together, and split MYCN non-amp by FOXR2+/-


```r
ecnb_counts_tidy %>% filter(gene_symbol == "FOXR2") %>% select(Risk, Stage) %>% table
```

```
##        Stage
## Risk    1-3;4S   4
##   HR        28 192
##   LR/IR    304  37
```

```r
ecnb_counts_tidy <- ecnb_counts_tidy %>%
    mutate(Group = case_when(
        MYCN == "Amp" ~ "EC-NB: MYCN Amp",
        MYCN == "NonAmp" & FOXR2_positive == "Y" ~ "EC-NB: FOXR2+, MYCN NonAmp",
        MYCN == "NonAmp" & FOXR2_positive == "N" ~ "EC-NB: FOXR2-, MYCN NonAmp"))

palette_groups_mut <- palette_groups
palette_groups_mut <- c(palette_groups_mut, "EC-NB: MYCN Amp" = "#A575D3")
palette_groups_mut <- c(palette_groups_mut, "EC-NB: FOXR2+, MYCN NonAmp" = "#f382c6")
palette_groups_mut <- c(palette_groups_mut, "EC-NB: FOXR2-, MYCN NonAmp" = "#215ba3")
```

Subset to stage 4 and high risk samples only:


```r
ecnb_counts_tidy_S4_HR <- ecnb_counts_tidy %>%
                  filter(Risk == "HR") %>%
                  filter(Stage == "4")
```

TF fingerprint, Y-axis set to 2.5k except for FOXR2


```r
y_maxes <- rep(2500, times = length(tfs_no_ecnb$Gene_hg))
names(y_maxes) <- tfs_no_ecnb$Gene_hg
y_maxes[which(names(y_maxes) %in% c("FOXR2"))] <- 450

tf_exp_boxplots(y_maxes = y_maxes,
                counts = ecnb_counts_tidy_S4_HR,
                palette = palette_groups_mut,
                label_y = c(1,2))
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/05//ecnb_boxplot-stage4-highrisk-1.png)<!-- -->

Print the list of tumors in order, to label the rows of this plot:


```r
names(palette_groups_mut)[c(13:15)] %>% rev
```

```
## [1] "EC-NB: FOXR2-, MYCN NonAmp" "EC-NB: FOXR2+, MYCN NonAmp"
## [3] "EC-NB: MYCN Amp"
```


### TARGET cohort

For extra-cranial neuroblastomas from TARGET dataset, data was processed in-house
by Steven Hebert in April 2024. Data was subsetted to only high risk patients which were
Stage 3 or 4. (4S not included.) This gave a total of 128 patient samples.

Original gene annotation in this data differs from our pipeline gene IDs. Therefore, Steven
intersected their Ensembl IDs with Ensembl IDs in our pipeline gene annotation, and 
removed all other gene entries. These modified counts matrices are stored at:

```
/project/kleinman/steven.hebert/from_narval/2024/2023-05-16-FOXR2paper/2024-01-29-revision/EC-NB-TARGET_counts/counts_ensembl_common
```

Load DESeq-normalized counts (e.g. for plotting gene expression box plots)


```r
# load NORM counts
target_norm <- here("data/RNAseq/external_data/TARGET_ECNB/counts_ensembl_common/GENCODEv36.norm.tsv.gz") %>%
    read.table(header = T, row.names = 1, check.names = F, sep = "\t")

target_norm[1:5, 1:5]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["TARGET-30-PAIFXV-01A"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["TARGET-30-PAIPGU-01A"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["TARGET-30-PAISNS-01A"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["TARGET-30-PAITCI-01A"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["TARGET-30-PAITEG-01A"],"name":[5],"type":["dbl"],"align":["right"]}],"data":[{"1":"1612.983582","2":"1337.289410","3":"3232.9537","4":"1649.4543","5":"1697.629906","_rn_":"ENSG00000000003:TSPAN6"},{"1":"9.170469","2":"8.985147","3":"0.0000","4":"0.0000","5":"0.607381","_rn_":"ENSG00000000005:TNMD"},{"1":"1610.945700","2":"2105.519496","3":"2272.0982","4":"4908.1899","5":"2858.335005","_rn_":"ENSG00000000419:DPM1"},{"1":"543.095546","2":"489.690523","3":"816.3660","4":"655.7644","5":"368.680269","_rn_":"ENSG00000000457:SCYL3"},{"1":"177.295732","2":"371.386085","3":"686.3254","4":"309.5681","5":"126.335249","_rn_":"ENSG00000000460:C1orf112"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Create metadata for FOXR2 expression status:


```r
# DESeq norm expression threshold used for FOXR2 +/-
threshold <- 2

target_counts_tidy <- target_norm %>% 
    as.data.frame() %>% 
    rownames_to_column("ID") %>% 
    separate(col = "ID", sep = ":", into = c("ENS", "gene_symbol")) %>% 
    pivot_longer(cols = -c("ENS", "gene_symbol"), names_to = "sample", values_to = "gene_expression") 
    
foxr2_pos_samples_target <- target_counts_tidy %>% 
    filter(gene_symbol == "FOXR2") %>% 
    filter(gene_expression > threshold) %>% 
    .$sample 

target_counts_tidy <- target_counts_tidy %>% 
    mutate(FOXR2_positive = case_when(sample %in% foxr2_pos_samples_target ~ "Y",
                                      T ~ "N")) %>% 
    distinct %>%     
    mutate(Group = case_when(
        FOXR2_positive == "Y" ~ "EC-NB: FOXR2+",
        FOXR2_positive == "N" ~ "EC-NB: FOXR2-")) %>% 
    mutate(Group = factor(Group, levels = c("EC-NB: FOXR2+",
                                            "EC-NB: FOXR2-")))
```

TF fingerprint. 

Since the genes differ in their range of expression, we set a y axis scale of 7000 per gene,
except FOXR2.


```r
palette_groups_mut_2 <- palette_groups
palette_groups_mut_2 <- c(palette_groups_mut_2, "EC-NB: FOXR2+" = "#f382c6")
palette_groups_mut_2 <- c(palette_groups_mut_2, "EC-NB: FOXR2-" = "#215ba3")

y_maxes <- c(rep(7000, times = length(tfs_no_ecnb$Gene_hg)))
names(y_maxes) <- tfs_no_ecnb$Gene_hg
y_maxes[which(names(y_maxes) %in% c("FOXR2"))] <- 1000
# y_maxes[which(names(y_maxes) %in% c("DBH"))] <- 40000

tf_exp_boxplots(y_maxes = y_maxes,
                counts = target_counts_tidy,
                palette = palette_groups_mut_2,
                label_y = c(1,2))
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/05//ecnb_target_boxplot_foxr2_7k-1.png)<!-- -->

Print the list of tumors in order, to label the rows of this plot:


```r
names(palette_groups_mut_2)[c(13:14)]
```

```
## [1] "EC-NB: FOXR2+" "EC-NB: FOXR2-"
```

## NB-FOXR2 only

FOXR2 positivity is also found in some extracranial neuroblastomas, and it is possible that FOXR2 activates key markers and transcription factors of EC-NB even within NB-FOXR2.

To confirm this is not the case, we plot the expression of both MGE TFs & EC-NB TFs in NB-FOXR2.


```r
tfs_mge_ecnb <- tfs$Gene_hg[which(!(tfs$Gene_hg %in% c("PAX6",
                                               "EMX1",
                                               "EMX2",
                                               "TBR1",
                                               "EOMES",
                                               "GSX2",
                                               "DLX1",
                                               "DLX2",
                                               "FOXR2")))]

df <- counts_lineage %>% 
    filter(Group %in% c("NB-FOXR2")) %>% 
    filter(gene_symbol %in% tfs_mge_ecnb)

df$gene_symbol <- factor(df$gene_symbol, 
                         rev(c("FOXG1",
                               "NKX2-1",
                               "LHX6",
                               "DLX5",
                               "DLX6",
                               "PHOX2B",
                               "ISL1",
                               "GATA3",
                               "LMO1",
                               "DBH",
                               "HAND1",
                               "HAND2",
                               "ASCL1")))


df %>% 
    ggplot(aes(x = gene_symbol, y = gene_expression, fill = "red")) +
    geom_boxplot(outlier.shape = NA) +
    #geom_jitter(size = 2, alpha = 0.8) +
    ylab("Normalized expression") + xlab("Gene") +
    scale_fill_manual(values = c("red")) +
    ylim(0,7500)+
    no_legend() + coord_flip() +
    ggtitle("MGE & ECNB TF expression in NB-FOXR2")
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/05//nb-foxr2-boxplot-mge-ecnb-tfs-1.png)<!-- -->

### ASCL1 & HAND2 expression across intracranial tumors

We see in the previous plot that NB-FOXR2 do express ASCL1 and HAND2, characteristic factors of extracranial neuroblastoma. 

Here, we investigate whether these two markers are also expressed across the other bulk intracranial tumors of our cohort, meaning they are not specific to neuroblastoma or FOXR2+ tumors. Outlier points were removed from violin plots.


```r
# load counts 
counts_ASCL1_HAND2 <- extract_pipeline_counts(path = here("data/RNAseq/pipeline_l3/2023-05-test_pbt/counts/Ensembl.ensGene.exon.norm.tsv.gz"),
                                        goi = c("ASCL1", "HAND2")) %>% 
    left_join(info_samples, by = c("sample" = "Nickname")) %>% 
    mutate(Group = factor(Group, levels = names(palette_groups)))
```

```
## Joining with `by = join_by(gene_symbol)`
```


```r
counts_ASCL1_HAND2 %>%
    filter(gene_symbol == "ASCL1") %>% 
    mutate(Group = factor(Group, c("NB-FOXR2",
                                 "HGG-H3.3G34R/V",
                                 "DIPG-H3K27M-FOXR2",
                                 "DIPG-H3K27M",
                                 "HGG-IDH",
                                 "EP-PFA",
                                 "ETMR",
                                 "MB-SHH",
                                 "MB-WNT"))) %>% 
    ggplot(aes(x = Group, y = gene_expression, fill = Group)) +
    geom_boxplot(aes(fill = Group), outlier.shape = NA) +
    coord_cartesian(ylim = c(0,10000)) +
    #geom_jitter(size = 2, alpha = 0.8) +
    scale_fill_manual(values = palette_groups)  +
    ylab("Normalized expression") + xlab("Group") +
    no_legend() + rotate_x() +
    ggtitle("ASCL1 expression")
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/05//ascl1_expression_all-1.png)<!-- -->



```r
counts_ASCL1_HAND2 %>%
    filter(gene_symbol == "HAND2") %>% 
    mutate(Group = factor(Group, c("NB-FOXR2",
                             "HGG-H3.3G34R/V",
                             "DIPG-H3K27M-FOXR2",
                             "DIPG-H3K27M",
                             "HGG-IDH",
                             "EP-PFA",
                             "ETMR",
                             "MB-SHH",
                             "MB-WNT"))) %>% 
    ggplot(aes(x = Group, y = gene_expression, fill = Group)) +
    geom_boxplot(aes(fill = Group), outlier.shape = NA) +
    coord_cartesian(ylim = c(0,5000)) +
    #geom_jitter(size = 2, alpha = 0.8) +
    scale_fill_manual(values = palette_groups)  +
    ylab("Normalized expression") + xlab("Group") +
    no_legend() + rotate_x() +
    ggtitle("HAND2 expression")
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/05//hand2_expression_all-1.png)<!-- -->

## All PBT

Since the genes differ in their range of expression, we set a y axis scale of 7000 per gene,
except FOXR2.

Plot MGE TF Fingerprint in all tumors except DIPG-FOXR2 (which was plotted above):


```r
y_maxes <- rep(x = 7000,
               times = length(tfs_no_ecnb$Gene_hg))
names(y_maxes) <- tfs_no_ecnb$Gene_hg
y_maxes[which(names(y_maxes) == "FOXR2")] <- 2000

tf_exp_boxplots(y_maxes = y_maxes,
                counts = counts_lineage %>% filter(Group != "DIPG-H3K27M-FOXR2"),
                palette = palette_groups,
                label_y = c(1,2))
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/05//boxplots_bulk_TFs_7k-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/05//boxplots_bulk_TFs_7k...*]~</span>

Print the list of tumors in order, to label the rows of this plot:


```r
names(palette_groups)[c(1,3:9)]
```

```
## [1] "NB-FOXR2"       "DIPG-H3K27M"    "HGG-IDH"        "EP-PFA"        
## [5] "HGG-H3.3G34R/V" "ETMR"           "MB-SHH"         "MB-WNT"
```

## OL marker expression in all PBT

We will also inspect expression of key oligodendrocyte markers across intracranial tumor types for comparison with NB-FOXR2.

Load data:


```r
info_samples <- read_tsv(file.path(here("data/RNAseq/pipeline_l3/2023-05-test_pbt/info.samples.tsv")))
```

```
## Rows: 121 Columns: 3
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (3): ID, Nickname, Group
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
counts_bulk <- extract_pipeline_counts(path = file.path(here("data/RNAseq/pipeline_l3/2023-05-test_pbt/counts/Ensembl.ensGene.exon.norm.tsv.gz")),
                                        goi = c("OLIG2","SOX10","PDGFRA")) %>% 
    left_join(info_samples, by = c("sample" = "Nickname")) %>% 
    mutate(Group = factor(Group, levels = names(palette_groups)))
```

```
## Joining with `by = join_by(gene_symbol)`
```

Produce box plots:


```r
counts_bulk %>%
    ggplot(aes(x = Group, y = gene_expression)) +
    geom_boxplot(aes(fill = Group), outlier.shape = NA) +
    geom_jitter(size = 2, alpha = 0.8) +
    scale_fill_manual(values = palette_groups)  +
    facet_wrap(. ~ gene_symbol, scales = "free_y") +
    ylab("Normalized expression") + xlab("Group") +
    no_legend() + rotate_x() 
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/05//box-plot-bulk-rna-OL-markers-1.png)<!-- -->

# Tumor scRNAseq

Here, we will inspect expression of the MGE TF fingerprint (developed in above sections) 
as well as key transcription factors and gene markers for extracranial neuroblastoma 
(EC-NB), in bulk RNAseq profiles of human tumors (intracranial pediatric brain tumors as well as EC-NB).

## Load data

Load per-sample seurat objects:


```r
samples_sc <- map(sc_samples_foxr2, ~ get(load(meta_sc[meta_sc$ID == .x, ]$Path)))
names(samples_sc) <- sc_samples_foxr2
```

Get normal/malignant calls:


```r
load(here("output/04/merged_sc_meta.Rda"))

# correct a discrepancy in how barcodes are named for the multiome sample
samples_sc$`P-6778_S-10155`@meta.data$cell.barcode <- samples_sc$`P-6778_S-10155`@meta.data$gex_barcode
merged_sc_meta <- merged_sc_meta %>%
    mutate(cell.barcode = ifelse(is.na(cell.barcode), gex_barcode, cell.barcode))

# load the malignant normal calls and add them to each seurat object
samples_sc <- map(samples_sc, function(seurat) {

    # sanity check that cell order is preserved in each metadata slot
    merged_barcodes <- merged_sc_meta %>% filter(cell.barcode %in% seurat@meta.data$cell.barcode) %>% pull(cell.barcode)
    all(seurat@meta.data$cell.barcode == merged_barcodes)

    tumor_normal <- merged_sc_meta %>% filter(cell.barcode %in% seurat@meta.data$cell.barcode) %>% pull(Malignant_normal)
    seurat <- AddMetaData(seurat, tumor_normal, "Malignant_normal")

    return(seurat)

})
```

```
## Loading required package: Signac
```

## Generate bubbleplots

Define function to plot bubble plot in each sample, split by normal and 
malignant cells, and annotated with the number of cells in each group:

<details>


```r
plot_bubble_tumor_norm_mal_pseudobulk <- function(seurat, id, tf_list,
                                                  collapse_likely = FALSE,
                                                  remove_gene_labels = FALSE) {

    tfs <- tf_list

    # if specified in arguments, collapse "Likely normal" cells into "Normal" category
    if(collapse_likely){

      categories <- c("Normal", "Malignant")
      new_meta <- seurat@meta.data %>%
        select("Malignant_normal") %>%
        mutate(Malignant_normal = case_when(Malignant_normal %in% c("Normal", "Likely normal") ~ "Normal",
                                            Malignant_normal == "Malignant"                    ~ "Malignant"))
      seurat <- AddMetaData(seurat, new_meta)
      rm(new_meta) # clean up

    } else {

      categories <- c("Normal", "Likely normal", "Malignant")

    }

    # get number of cells in each category
    n_cells_per_category <- seurat@meta.data %>%
      select("Malignant_normal") %>%
      group_by(Malignant_normal) %>%
      summarize("N" = n()) %>%
      mutate(N = glue("n = {N}")) %>%
      dplyr::rename(Cluster = Malignant_normal)

    # extract average expression and detection rate per cluster for the TF
    # fingerprint genes
    df <- extract_meanexp_pct(
        seurat,
        genes = tfs$Gene_hg,
        cluster_col = "Malignant_normal",
        scale = TRUE)

    # add empty rows for the missing genes/clusters
    df_complete <- df %>%
        mutate(Gene = factor(Gene, levels = tfs$Gene_hg),
               Cluster = factor(Cluster, levels = c("dummy", unique(.$Cluster)))) %>%
        complete(Gene, Cluster, fill = list(Expression = 0, Pct1 = 0))

    # add a dummy cluster that expresses all genes, to force showing all genes
    df_complete[df_complete$Cluster == "dummy", ]$Expression <- 0.1
    df_complete[df_complete$Cluster == "dummy", ]$Pct1 <- 0.1

    # make the bubbleplot
    p <- plot_bubble(df_complete %>% filter(Pct1 >= 0.01),
                genes = tfs$Gene_hg,
                cluster_order = c(categories, "dummy"),
                max_radius = 10) +
        ggtitle(id) +
        theme(legend.position = "bottom",
              legend.box = "vertical") +
        geom_text(data = n_cells_per_category,
                  aes(y = (length(tfs$Gene_hg) + 1), label = N),
                  angle = 45, vjust = -1, hjust = 0) +
        expand_limits(y = c(-1, (length(tfs$Gene_hg) + 6 )))

    if(remove_gene_labels){
      p <- p + theme(axis.text.y = element_blank())
    }

    return(p)

}
```

</details>

Minimal TF set:

(Note that only the Malignant cells column for each sample was used in the final figure.)


```r
imap(samples_sc, ~ plot_bubble_tumor_norm_mal_pseudobulk(.x, .y,
                                                         tf_list = tfs_min,
                                                         collapse_likely = TRUE)) %>%
    {plot_grid(plotlist = ., nrow = 1,
               align = "h", axis = "tblr")}
```

```
## @ removing 0 clusters with < 20 cells
## @ removing 0 clusters with < 20 cells
## @ removing 0 clusters with < 20 cells
## @ removing 0 clusters with < 20 cells
## @ removing 0 clusters with < 20 cells
## @ removing 0 clusters with < 20 cells
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/05//bubble-plot-sc-tumors-min-tf-mal-only-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/05//bubble-plot-sc-tumors-min-tf-mal-only...*]~</span>


# FOXR2-transduced human NSCs

Here, we load bulk RNAseq from [Tsai et al. *Cancer Research* 2022](https://pubmed.ncbi.nlm.nih.gov/35802025/)
where the authors transduced H9 human neural stem cells with HA-FOXR2 expression.

We will compare the expression of our TF fingerprint in FOXR2-transduced vs.
control (transduced with Hc-Red control vector) hNSCs, to assess whether 
FOXR2 itself can lead to expression of the TF fingerprint.

## Load data 



```r
# define mouse and human TFs
tfs <- tribble(~Gene_hg, ~Gene_mm,
               "FOXR2",  "Foxr2",
               "FOXG1",  "Foxg1",
               "PAX6",   "Pax6",
               "EMX1",   "Emx1",
               "EMX2",   "Emx2",
               "TBR1",   "Tbr1",
               "EOMES",  "Eomes",
               "GSX2",   "Gsx2",
               "NKX2-1", "Nkx2-1",
               "LHX6",   "Lhx6",
               "DLX1",   "Dlx1",
               "DLX2",   "Dlx2",
               "DLX5",   "Dlx5",
               "DLX6",   "Dlx6")

info_samples_hNSC <- read_tsv(here("data/RNAseq/pipeline_l3/2023-04-FOXR2_hNSCs/info.samples.tsv"))
```

```
## Rows: 10 Columns: 3
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (3): ID, Nickname, Group
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
info_samples_hNSC %>% select(Nickname, Group)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Nickname"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Group"],"name":[2],"type":["chr"],"align":["left"]}],"data":[{"1":"hNSC_FOXR2_1","2":"FOXR2"},{"1":"hNSC_FOXR2_2","2":"FOXR2"},{"1":"hNSC_FOXR2_3","2":"FOXR2"},{"1":"hNSC_FOXR2_4","2":"FOXR2"},{"1":"hNSC_FOXR2_5","2":"FOXR2"},{"1":"hNSC_Ctrl_1","2":"Control"},{"1":"hNSC_Ctrl_2","2":"Control"},{"1":"hNSC_Ctrl_3","2":"Control"},{"1":"hNSC_Ctrl_4","2":"Control"},{"1":"hNSC_Ctrl_5","2":"Control"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# load counts for FOXR2
counts_nsc <- extract_pipeline_counts(path = here("data/RNAseq/pipeline_l3/2023-04-FOXR2_hNSCs/counts/Ensembl.ensGene.exon.norm.tsv.gz"),
                                      goi = tfs$Gene_hg) %>%
    left_join(info_samples_hNSC, by = c("sample" = "Nickname"))
```

```
## Joining with `by = join_by(gene_symbol)`
```

```r
palette_nsc <- c("FOXR2" = "red", "Control" = "gray90")
```

## Generate box plots

Define a function to plot hNSC FOXR2 vs. control TF expression as a box plot:

<details>


```r
# Input:
# - y_maxes: vector of y-axis maximum values, named with the gene symbols to plot
#   --> e.g. c("FOXR2" = 2000, "LHX6" = 4000)
# - counts: tidy counts matrix which includes columns: gene_symbol, gene_expression, sample, Group
# - palette: palette for Group variable in counts matrix  
# - label_y: OPTIONAL list of indices for box plots that require y axis value labels
#   --> e.g. c(1,2)
#   --> If not provided, only the first box plot will have a label (= c(1))
# - strip_label: OPTIONAL boolean (T/F) whether to add right side strip label for "group"
#   --> default is False (no label)

tf_exp_boxplots_hnsc <- function(y_maxes,
                            counts,
                            palette,
                            label_y = c(1),
                            strip_label = F) {
    
    
    # loop over genes & the maximum y values for each gene,
    # making boxplots for each gene as 1 column to be combined with plot_grid
    # .x = ymaxes (vector values), .y = genes (vector names)
    plots <- imap(y_maxes,
                  ~ counts %>%
                      # subset counts
                      filter(gene_symbol == .y) %>%
                      mutate(Group = factor(Group, levels = c("FOXR2", "Control"))) %>% 
                      ggplot(aes(x = factor(1), y = gene_expression)) +
                      # make boxplot with jittered points
                      geom_boxplot(aes(fill = Group), outlier.shape = NA, width = 0.3) +
                      geom_jitter(aes(fill = Group), size = 0.5, width = 0.2, shape = 21) +
                      scale_fill_manual(values = palette)  +
                      facet_wrap(~ gene_symbol, nrow = 1) +
                      # subset the plotting area (y-axis) to the specified y-max value
                      coord_cartesian(ylim = c(0, .x)) +
                      # only show the min and max of the y-axis
                      scale_y_continuous(breaks = c(0, .x)) +
                      theme_min2() +
                      no_legend() + ggtitle(.y) + xlab(NULL) + ylab(NULL) +
                      # put the group label on the right instead of the top
                      theme(strip.text.x = element_text(angle = 0, size = 2)) +
                      theme(aspect.ratio = 2))
    
    if(strip_label){
        # remove the y-axis strip labels in all but the last plot
        plots[1:(length(y_maxes) - 1)] <- map(plots[1:(length(y_maxes) - 1)], ~ .x + theme(strip.text.y = element_blank()))
    } else {
        # remove the y-axis strip labels in ALL plots
        plots[1:length(y_maxes)] <- map(plots[1:length(y_maxes)], ~ .x + theme(strip.text.y = element_blank()))
    }

    # remove the y-axis value labels for all except specified plots
    y_axis_rem <- setdiff(1:length(y_maxes), label_y)
    plots[y_axis_rem] <- map(plots[y_axis_rem], ~ .x + theme(axis.text.y = element_blank()))
    
    plot_grid(plotlist = plots, nrow = 1, align = "h", axis = "tb")
}
```

</details>

Generate plot:

Since the genes differ in their range of expression, we set a y axis scale of 7000 per gene,
except FOXR2.


```r
# Print max value of gene expression (except FOXR2)
counts_nsc %>% 
    filter(gene_symbol != "FOXR2") %>%
    filter(gene_symbol %in% tfs_no_ecnb$Gene_hg) %>%
    arrange(-gene_expression) %>%
    .$gene_expression %>%
    .[1]
```

```
## [1] 8356.808
```

```r
y_maxes_nsc <- rep(7000, times = length(tfs_no_ecnb$Gene_hg))
names(y_maxes_nsc) <- tfs_no_ecnb$Gene_hg
y_maxes_nsc[which(names(y_maxes_nsc) %in% c("FOXR2"))] <- 200000

tf_exp_boxplots_hnsc(y_maxes = y_maxes_nsc,
                counts = counts_nsc %>%
                            mutate(Group = factor(Group, levels = c("FOXR2", "Control"))),
                palette = palette_nsc,
                label_y = c(1,2))
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/05//hnsc_boxplot-7k-no-ecnb-tfs-1.png)<!-- -->


<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2024-11-04 11:33:24
```

The git repository and last commit:

```
## Local:    main /project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public
## Remote:   main @ origin (https://github.com/fungenomics/NB-FOXR2.git)
## Head:     [0e89693] 2024-09-12: Initial commit
```

The random seed was set with `set.seed(100)`

The R session info:
<details>

```
## ─ Session info ───────────────────────────────────────────────────────────────
##  setting  value
##  version  R version 4.1.2 (2021-11-01)
##  os       Rocky Linux 8.10 (Green Obsidian)
##  system   x86_64, linux-gnu
##  ui       X11
##  language (EN)
##  collate  en_US.UTF-8
##  ctype    en_US.UTF-8
##  tz       America/Toronto
##  date     2024-11-04
##  pandoc   1.19.2.1 @ /cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/ (via rmarkdown)
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  ! package          * version    date (UTC) lib source
##  P abind              1.4-5      2016-07-21 [?] CRAN (R 4.1.2)
##  P backports          1.4.1      2021-12-13 [?] CRAN (R 4.1.2)
##  P BiocGenerics       0.40.0     2021-10-26 [?] Bioconductor
##  P BiocManager        1.30.15    2021-05-11 [?] CRAN (R 4.1.2)
##  P BiocParallel       1.28.3     2021-12-09 [?] Bioconductor
##  P Biostrings         2.62.0     2021-10-26 [?] Bioconductor
##  P bit                4.0.4      2020-08-04 [?] CRAN (R 4.1.2)
##  P bit64              4.0.5      2020-08-30 [?] CRAN (R 4.1.2)
##  P bitops             1.0-7      2021-04-24 [?] CRAN (R 4.1.2)
##  P broom            * 1.0.1      2022-08-29 [?] CRAN (R 4.1.2)
##  P bslib              0.3.1      2021-10-06 [?] CRAN (R 4.1.2)
##  P cachem             1.0.6      2021-08-19 [?] CRAN (R 4.1.2)
##  P callr              3.7.6      2024-03-25 [?] RSPM
##  P cellranger         1.1.0      2016-07-27 [?] CRAN (R 4.1.2)
##  P class              7.3-22     2023-05-03 [?] CRAN (R 4.1.2)
##  P cli                3.6.1      2023-03-23 [?] RSPM (R 4.1.2)
##  P cluster            2.1.2      2021-04-17 [?] CRAN (R 4.1.2)
##  P codetools          0.2-18     2020-11-04 [?] CRAN (R 4.1.2)
##  P colorspace         2.0-2      2021-06-24 [?] CRAN (R 4.1.2)
##  P conflicted         1.2.0      2023-02-01 [?] CRAN (R 4.1.2)
##  P cowplot          * 1.1.1      2020-12-30 [?] CRAN (R 4.1.2)
##  P crayon             1.4.2      2021-10-29 [?] CRAN (R 4.1.2)
##  P crosstalk          1.2.0      2021-11-04 [?] CRAN (R 4.1.2)
##  P data.table       * 1.14.2     2021-09-27 [?] CRAN (R 4.1.2)
##  P deldir             1.0-6      2021-10-23 [?] CRAN (R 4.1.2)
##  P devtools           2.4.5      2022-10-11 [?] CRAN (R 4.1.2)
##  P dials            * 1.2.0      2023-04-03 [?] CRAN (R 4.1.2)
##  P DiceDesign         1.9        2021-02-13 [?] CRAN (R 4.1.2)
##  P digest             0.6.35     2024-03-11 [?] CRAN (R 4.1.2)
##  P docopt             0.7.1      2020-06-24 [?] CRAN (R 4.1.2)
##  P dplyr            * 1.1.1      2023-03-22 [?] CRAN (R 4.1.2)
##  P DT                 0.30       2023-10-05 [?] CRAN (R 4.1.2)
##  P ellipsis           0.3.2      2021-04-29 [?] CRAN (R 4.1.2)
##  P evaluate           0.23       2023-11-01 [?] CRAN (R 4.1.2)
##  P fansi              1.0.2      2022-01-14 [?] CRAN (R 4.1.2)
##  P farver             2.1.0      2021-02-28 [?] CRAN (R 4.1.2)
##  P fastmap            1.1.0      2021-01-25 [?] CRAN (R 4.1.2)
##  P fastmatch          1.1-3      2021-07-23 [?] CRAN (R 4.1.2)
##  P feather            0.3.5      2019-09-15 [?] RSPM (R 4.1.2)
##  P fitdistrplus       1.1-6      2021-09-28 [?] CRAN (R 4.1.2)
##  P foreach            1.5.1      2020-10-15 [?] CRAN (R 4.1.2)
##  P fs                 1.5.2      2021-12-08 [?] CRAN (R 4.1.2)
##  P furrr              0.3.1      2022-08-15 [?] CRAN (R 4.1.2)
##  P future             1.25.0     2022-04-24 [?] CRAN (R 4.1.2)
##  P future.apply       1.8.1      2021-08-10 [?] CRAN (R 4.1.2)
##  P generics           0.1.3      2022-07-05 [?] CRAN (R 4.1.2)
##  P GenomeInfoDb       1.30.1     2022-01-30 [?] Bioconductor
##  P GenomeInfoDbData   1.2.4      2023-11-28 [?] Bioconductor
##  P GenomicRanges      1.46.1     2021-11-18 [?] Bioconductor
##  P ggforce            0.3.3      2021-03-05 [?] CRAN (R 4.1.2)
##  P ggplot2          * 3.4.2      2023-04-03 [?] CRAN (R 4.1.2)
##  P ggrepel            0.9.1      2021-01-15 [?] CRAN (R 4.1.2)
##  P ggridges           0.5.3      2021-01-08 [?] CRAN (R 4.1.2)
##  P ggseqlogo          0.1        2017-07-25 [?] CRAN (R 4.1.2)
##  P git2r              0.29.0     2021-11-22 [?] CRAN (R 4.1.2)
##  P globals            0.14.0     2020-11-22 [?] CRAN (R 4.1.2)
##  P glue             * 1.6.2      2022-02-24 [?] CRAN (R 4.1.2)
##  P goftest            1.2-3      2021-10-07 [?] CRAN (R 4.1.2)
##  P gower              1.0.1      2022-12-22 [?] CRAN (R 4.1.2)
##  P GPfit              1.0-8      2019-02-08 [?] CRAN (R 4.1.2)
##  P gridExtra        * 2.3        2017-09-09 [?] CRAN (R 4.1.2)
##  P gtable             0.3.0      2019-03-25 [?] CRAN (R 4.1.2)
##  P hardhat            1.3.0      2023-03-30 [?] CRAN (R 4.1.2)
##  P here             * 1.0.1      2020-12-13 [?] CRAN (R 4.1.2)
##  P highr              0.9        2021-04-16 [?] CRAN (R 4.1.2)
##  P hms                1.1.1      2021-09-26 [?] CRAN (R 4.1.2)
##  P htmltools          0.5.2      2021-08-25 [?] CRAN (R 4.1.2)
##  P htmlwidgets        1.5.4      2021-09-08 [?] CRAN (R 4.1.2)
##  P httpuv             1.6.5      2022-01-05 [?] CRAN (R 4.1.2)
##  P httr               1.4.2      2020-07-20 [?] CRAN (R 4.1.2)
##  P ica                1.0-2      2018-05-24 [?] CRAN (R 4.1.2)
##  P igraph             2.0.3      2024-03-13 [?] CRAN (R 4.1.2)
##  P infer            * 1.0.4      2022-12-02 [?] CRAN (R 4.1.2)
##  P ipred              0.9-14     2023-03-09 [?] CRAN (R 4.1.2)
##  P IRanges            2.28.0     2021-10-26 [?] Bioconductor
##  P irlba              2.3.5      2021-12-06 [?] CRAN (R 4.1.2)
##  P iterators          1.0.13     2020-10-15 [?] CRAN (R 4.1.2)
##  P jquerylib          0.1.4      2021-04-26 [?] CRAN (R 4.1.2)
##  P jsonlite           1.8.8      2023-12-04 [?] CRAN (R 4.1.2)
##  P KernSmooth         2.23-20    2021-05-03 [?] CRAN (R 4.1.2)
##  P knitr              1.37       2021-12-16 [?] CRAN (R 4.1.2)
##  P labeling           0.4.2      2020-10-20 [?] CRAN (R 4.1.2)
##  P later              1.3.0      2021-08-18 [?] CRAN (R 4.1.2)
##  P lattice            0.20-45    2021-09-22 [?] CRAN (R 4.1.2)
##  P lava               1.7.2.1    2023-02-27 [?] CRAN (R 4.1.2)
##  P lazyeval           0.2.2      2019-03-15 [?] CRAN (R 4.1.2)
##  P leiden             0.3.9      2021-07-27 [?] CRAN (R 4.1.2)
##  P lhs                1.1.6      2022-12-17 [?] CRAN (R 4.1.2)
##  P LiblineaR          2.10-24    2024-09-13 [?] CRAN (R 4.1.2)
##  P lifecycle          1.0.3      2022-10-07 [?] CRAN (R 4.1.2)
##  P listenv            0.8.0      2019-12-05 [?] CRAN (R 4.1.2)
##  P lmtest             0.9-39     2021-11-07 [?] CRAN (R 4.1.2)
##  P lsa                0.73.2     2020-05-04 [?] CRAN (R 4.1.2)
##  P lubridate          1.9.2      2023-02-10 [?] CRAN (R 4.1.2)
##  P magrittr         * 2.0.3      2022-03-30 [?] CRAN (R 4.1.2)
##  P MASS               7.3-54     2021-05-03 [?] CRAN (R 4.1.2)
##  P Matrix             1.3-4      2021-06-01 [?] CRAN (R 4.1.2)
##  P matrixStats        0.61.0     2021-09-17 [?] CRAN (R 4.1.2)
##  P memoise            2.0.1      2021-11-26 [?] CRAN (R 4.1.2)
##  P mgcv               1.8-38     2021-10-06 [?] CRAN (R 4.1.2)
##  P mime               0.12       2021-09-28 [?] CRAN (R 4.1.2)
##  P miniUI             0.1.1.1    2018-05-18 [?] CRAN (R 4.1.2)
##  P modeldata        * 1.1.0      2023-01-25 [?] CRAN (R 4.1.2)
##  P munsell            0.5.0      2018-06-12 [?] CRAN (R 4.1.2)
##  P nlme               3.1-153    2021-09-07 [?] CRAN (R 4.1.2)
##  P nnet               7.3-19     2023-05-03 [?] CRAN (R 4.1.2)
##  P parallelly         1.30.0     2021-12-17 [?] CRAN (R 4.1.2)
##  P parsnip          * 1.1.0      2023-04-12 [?] CRAN (R 4.1.2)
##  P patchwork          1.1.1      2020-12-17 [?] CRAN (R 4.1.2)
##  P pbapply            1.5-0      2021-09-16 [?] CRAN (R 4.1.2)
##  P pillar             1.9.0      2023-03-22 [?] RSPM (R 4.1.2)
##  P pkgbuild           1.4.2      2023-06-26 [?] CRAN (R 4.1.2)
##  P pkgconfig          2.0.3      2019-09-22 [?] CRAN (R 4.1.2)
##  P pkgload            1.3.3      2023-09-22 [?] CRAN (R 4.1.2)
##  P plotly             4.10.0     2021-10-09 [?] CRAN (R 4.1.2)
##  P plyr               1.8.6      2020-03-03 [?] CRAN (R 4.1.2)
##  P png                0.1-7      2013-12-03 [?] CRAN (R 4.1.2)
##  P polyclip           1.10-0     2019-03-14 [?] CRAN (R 4.1.2)
##  P prettyunits        1.1.1      2020-01-24 [?] CRAN (R 4.1.2)
##  P processx           3.8.4      2024-03-16 [?] RSPM
##  P prodlim            2023.03.31 2023-04-02 [?] CRAN (R 4.1.2)
##  P profvis            0.3.8      2023-05-02 [?] CRAN (R 4.1.2)
##  P promises           1.2.0.1    2021-02-11 [?] CRAN (R 4.1.2)
##  P ps                 1.7.6      2024-01-18 [?] RSPM
##  P purrr            * 1.0.1      2023-01-10 [?] CRAN (R 4.1.2)
##  P qlcMatrix          0.9.7      2018-04-20 [?] CRAN (R 4.1.2)
##  P R6                 2.5.1      2021-08-19 [?] CRAN (R 4.1.2)
##  P RANN               2.6.1      2019-01-08 [?] CRAN (R 4.1.2)
##  P RColorBrewer     * 1.1-2      2014-12-07 [?] CRAN (R 4.1.2)
##  P Rcpp               1.0.8      2022-01-13 [?] CRAN (R 4.1.2)
##  P RcppAnnoy          0.0.19     2021-07-30 [?] CRAN (R 4.1.2)
##  P RcppRoll           0.3.0      2018-06-05 [?] CRAN (R 4.1.2)
##  P RCurl              1.98-1.5   2021-09-17 [?] CRAN (R 4.1.2)
##  P readr            * 2.1.1      2021-11-30 [?] CRAN (R 4.1.2)
##  P readxl           * 1.3.1      2019-03-13 [?] CRAN (R 4.1.2)
##  P recipes          * 1.0.5      2023-02-20 [?] CRAN (R 4.1.2)
##  P remotes            2.4.2.1    2023-07-18 [?] CRAN (R 4.1.2)
##  P renv               1.0.3      2023-09-19 [?] CRAN (R 4.1.2)
##  P reshape2           1.4.4      2020-04-09 [?] CRAN (R 4.1.2)
##  P reticulate         1.23       2022-01-14 [?] CRAN (R 4.1.2)
##  P rlang              1.1.3      2024-01-10 [?] CRAN (R 4.1.2)
##  P rmarkdown          2.11       2021-09-14 [?] CRAN (R 4.1.2)
##  P ROCR               1.0-11     2020-05-02 [?] CRAN (R 4.1.2)
##  P rpart              4.1-15     2019-04-12 [?] CRAN (R 4.1.2)
##  P rprojroot          2.0.2      2020-11-15 [?] CRAN (R 4.1.2)
##  P rsample          * 1.1.1      2022-12-07 [?] CRAN (R 4.1.2)
##  P Rsamtools          2.10.0     2021-10-26 [?] Bioconductor
##  P rstudioapi         0.13       2020-11-12 [?] CRAN (R 4.1.2)
##  P Rtsne              0.15       2018-11-10 [?] CRAN (R 4.1.2)
##  P S4Vectors          0.32.4     2022-03-24 [?] Bioconductor
##  P sass               0.4.0      2021-05-12 [?] CRAN (R 4.1.2)
##  P scales           * 1.2.1      2022-08-20 [?] CRAN (R 4.1.2)
##  P scattermore        0.7        2020-11-24 [?] CRAN (R 4.1.2)
##  P sctransform        0.3.3      2022-01-13 [?] CRAN (R 4.1.2)
##  P sessioninfo        1.2.2      2021-12-06 [?] CRAN (R 4.1.2)
##  P Seurat           * 4.0.0      2021-01-30 [?] CRAN (R 4.1.2)
##  P SeuratObject     * 4.0.4      2021-11-23 [?] CRAN (R 4.1.2)
##  P shiny              1.7.1      2021-10-02 [?] CRAN (R 4.1.2)
##  P Signac           * 1.3.0      2021-07-12 [?] CRAN (R 4.1.2)
##  P slam               0.1-50     2022-01-08 [?] CRAN (R 4.1.2)
##  P SnowballC          0.7.0      2020-04-01 [?] CRAN (R 4.1.2)
##  P sparsesvd          0.2        2019-07-15 [?] CRAN (R 4.1.2)
##  P spatstat           1.64-1     2020-05-12 [?] CRAN (R 4.1.2)
##  P spatstat.data      2.1-2      2021-12-17 [?] CRAN (R 4.1.2)
##  P spatstat.utils     2.3-0      2021-12-12 [?] CRAN (R 4.1.2)
##  P stringi            1.7.6      2021-11-29 [?] CRAN (R 4.1.2)
##  P stringr          * 1.5.0      2022-12-02 [?] CRAN (R 4.1.2)
##  P survival           3.2-13     2021-08-24 [?] CRAN (R 4.1.2)
##  P tensor             1.5        2012-05-05 [?] CRAN (R 4.1.2)
##  P tibble           * 3.2.1      2023-03-20 [?] RSPM (R 4.1.2)
##  P tidymodels       * 1.0.0      2022-07-13 [?] CRAN (R 4.1.2)
##  P tidyr            * 1.3.0      2023-01-24 [?] CRAN (R 4.1.2)
##  P tidyselect         1.2.0      2022-10-10 [?] CRAN (R 4.1.2)
##  P timechange         0.2.0      2023-01-11 [?] CRAN (R 4.1.2)
##  P timeDate           4022.108   2023-01-07 [?] CRAN (R 4.1.2)
##  P tune             * 1.1.1      2023-04-11 [?] CRAN (R 4.1.2)
##  P tweenr             1.0.2      2021-03-23 [?] CRAN (R 4.1.2)
##  P tzdb               0.3.0      2022-03-28 [?] CRAN (R 4.1.2)
##  P urlchecker         1.0.1      2021-11-30 [?] CRAN (R 4.1.2)
##  P usethis            2.2.2      2023-07-06 [?] CRAN (R 4.1.2)
##  P utf8               1.2.2      2021-07-24 [?] CRAN (R 4.1.2)
##  P uwot               0.1.11     2021-12-02 [?] CRAN (R 4.1.2)
##  P vctrs              0.6.5      2023-12-01 [?] CRAN (R 4.1.2)
##  P viridis          * 0.5.1      2018-03-29 [?] RSPM (R 4.1.2)
##  P viridisLite      * 0.3.0      2018-02-01 [?] CRAN (R 4.1.2)
##  P vroom              1.5.7      2021-11-30 [?] CRAN (R 4.1.2)
##  P withr              2.5.0      2022-03-03 [?] CRAN (R 4.1.2)
##  P workflows        * 1.1.3      2023-02-22 [?] CRAN (R 4.1.2)
##  P workflowsets     * 1.0.1      2023-04-06 [?] CRAN (R 4.1.2)
##  P xfun               0.29       2021-12-14 [?] CRAN (R 4.1.2)
##  P xtable             1.8-4      2019-04-21 [?] CRAN (R 4.1.2)
##  P XVector            0.34.0     2021-10-26 [?] Bioconductor
##  P yaml               2.2.1      2020-02-01 [?] CRAN (R 4.1.2)
##  P yardstick        * 1.1.0      2022-09-07 [?] CRAN (R 4.1.2)
##  P zlibbioc           1.40.0     2021-10-26 [?] Bioconductor
##  P zoo                1.8-9      2021-03-09 [?] CRAN (R 4.1.2)
## 
##  [1] /project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/renv/library/R-4.1/x86_64-pc-linux-gnu
##  [2] /home/kleinman/bhavyaa.chandarana/.cache/R/renv/sandbox/R-4.1/x86_64-pc-linux-gnu/145cef2c
## 
##  P ── Loaded and on-disk path mismatch.
## 
## ──────────────────────────────────────────────────────────────────────────────
```

</details>


***

<!-- END OF END MATTER -->
