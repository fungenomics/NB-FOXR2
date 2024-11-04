---
title: "04 - QC and data preparation for single-cell tumor data"
author: "Selin Jessa [[selin.jessa@mail.mcgill.ca](mailto:selin.jessa@mail.mcgill.ca)] and Bhavyaa Chandarana [[bhavyaa.chandarana@mail.mcgill.ca](mailto:bhavyaa.chandarana@mail.mcgill.ca)]"
date: "01 November, 2024"
output:
  html_document:
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
## Document index: 04
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
## public/output/04
```

```
## public/figures/04//
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

Data for single-cell tumor samples (sc/snRNAseq, scMultiome) have been preprocessed
using our in-house workflow based on Seurat/Signac, which loads Cellranger output,
performs QC to filter cells, and applies normalization, dimensionality reduction,
and clustering. These are run under `data/singlecell/pipeline_sc*` (except 
`pipeline_scMultiome_mm` which is for processing of mouse samples).

For each sample, inferCNV was also applied, to infer normal and malignant cells
based on the single-cell data. This was run inside each sample's pipeline 
directory under `data/singlecell/pipeline_sc*/<sample>/inferCNV`. Malignant and normal
cells were called individually for each sample in the file `data/singlecell/pipeline_sc*/<sample>/inferCNV/window101_exp0.1_refG34normcleaned_HMMnone/infer_cnv_call.html`.

Following the preprocessing, all samples were merged into one object and
normalization/dimensionality reduction/clustering re-computed. This was performed
at `data/singlecell/integrations`. No correction for batch/sample/technology differences
were used.

In this document, we load outputs from these preprocessing steps in the form
of Seurat objects. We will perform post-cluster QC (removing any low quality clusters),
explore the merged dataset, and assign normal and malignant cells based on both
the inferCNV outputs, and the co-clustering of cells from different samples
in the merged object.


# Libraries


```r
library(here)
library(magrittr)
library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(readxl)
library(glue)
library(purrr)
library(ggplot2)
library(cowplot)
library(Seurat)
library(Signac)

source(here("include/style.R"))
source(here("code/functions/scRNAseq.R"))
source(here("code/functions/ssGSEA.R"))

ggplot2::theme_set(theme_min())

square_theme <- theme(aspect.ratio = 1)
large_text <- theme(text = element_text(size = 15))
```


# Load metadata

<div class="fold o">


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

</div>


```r
sc_samples_foxr2 <- meta_sc %>% filter(Group == "NB-FOXR2") %>% pull(ID)
palette_sample <- read_tsv(here("output/00/sc_info.samples.tsv")) %>% 
    select(Sample, Color) %>% tibble::deframe()
```

```
## Rows: 6 Columns: 5
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (5): Path, Sample, Technology, Group, Color
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```



# Post-cluster QC {.tabset}

Load the seurat objects, as produced by the preprocessing workflows.


```r
samples_sc <- map(sc_samples_foxr2, ~ get(load(meta_sc[meta_sc$ID == .x, ]$Path)))
names(samples_sc) <- sc_samples_foxr2
```



```r
# create a function to plot the violin plots for each QC stat
plot_vln_stat_per_cluster <- function(stat, lims = c(0, NA), labels = FALSE, seurat) {
    
    # get the requested statistic
    df <- seurat@meta.data
    df$Stat <- df[[stat]]
    
    # ggplot violin plot
    p1 <- df %>% ggplot(aes(x = seurat_clusters, y = Stat)) +
        geom_violin(aes(fill = seurat_clusters), scale = "width") +
        scale_fill_manual(values = seurat@misc$colours) +
        rotate_x() + no_legend() +
        # zoom in on y-axis if needed,
        # without changing distribution
        # and flip x/y axes
        coord_flip(ylim = lims) +
        xlab(NULL) +
        ggtitle(stat)
    
    # remove sample labels if needed
    # (after coord flip, samples are on the y-axis)
    if (!labels) p1 <- p1 + theme(axis.text.y = element_blank()) + ylab(NULL)
    
    p1
    
}

# wrapper function to plot all QC violin plots for one sample, given the seurat object
plot_qc_per_sample <- function(seurat, nCount_max = 20000, nFeature_max = 7000) {
    
    stats <- replicate(5, c(0, NA), simplify = FALSE)
    
    names(stats) <- c("nCount_RNA",
                      "nFeature_RNA",
                      "percent.mito",
                      "percent.ribo",
                      "G2M.Score")
    
    stats$nCount_RNA <- c(0, nCount_max)
    stats$nFeature_RNA <- c(0, nFeature_max)
    stats$G2M.Score <- c(-0.5, 1.5)
    
    # map over:
    # 1. stat
    # 2. y limits
    # 3. whether to keep labels or not
    pmap(list(names(stats), stats, c(TRUE, rep(FALSE, 4))),
         ~ plot_vln_stat_per_cluster(..1, ..2, ..3, seurat)) %>% 
        {plot_grid(plotlist = ., nrow = 1, align = "h", axis = "tb",
                   rel_widths = c(0.25, rep(0.18, 4)))}
    
}
```

Plot per-cluster QC in each sample (separated into tabs).

## P-2236_S-2236


```r
plot_qc_per_sample(samples_sc$`P-2236_S-2236`)
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/04//qc_P-2236_S-2236-1.png)<!-- -->

In this sample, we flag cluster 7 as low-QC for downstream analysis.

## P-2273_S-2273


```r
plot_qc_per_sample(samples_sc$`P-2273_S-2273`, nCount_max = 8000, nFeature_max = 5000)
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/04//qc_P-2273_S-2273-1.png)<!-- -->

## P-6776_S-10153


```r
plot_qc_per_sample(samples_sc$`P-6776_S-10153`)
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/04//qc_P-6776_S-10153-1.png)<!-- -->

## P-6777_S-10154


```r
plot_qc_per_sample(samples_sc$`P-6777_S-10154`)
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/04//qc_P-6777_S-10154-1.png)<!-- -->

## P-6778_S-10155


```r
plot_qc_per_sample(samples_sc$`P-6778_S-10155`, nCount_max = 25000, nFeature_max = 7000)
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/04//qc_P-6778_S-10155-1.png)<!-- -->

## NBFOXR2_6


```r
plot_qc_per_sample(samples_sc$`NBFOXR2_6`, nCount_max = 10000, nFeature_max = 4000)
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/04//qc_NBFOXR2_6-1.png)<!-- -->

# Naive-merged dataset

## Load data

Load integrated data (without batch correction or harmonization):


```r
load(here("data/singlecell/integrations/NB-FOXR2_naive_join/output/seurat_joint.Rda"))
seurat_joint$Joint_cluster <- Idents(seurat_joint)
```

## Overview

Plot an overview of the merged dataset:


```r
seurat_joint@meta.data <- seurat_joint@meta.data %>% 
    mutate(Technology = case_when(Sample == "NBFOXR2_6" ~ "10X Single Cell 3'",
                                  TRUE ~ as.character(Technology)))

plot_grid(
    umap(seurat_joint, color_by = "Technology", alpha = 0.5, label = FALSE, point_size = 0.5) +
        mod_legend_col() + square_theme + large_text,
    umap(seurat_joint, color_by = "Sample", colors = palette_sample, alpha = 0.5, label = FALSE, point_size = 0.5) +
        mod_legend_col() + square_theme + large_text,
    umap(seurat_joint, legend = TRUE, alpha = 0.5, point_size = 0.5) +
        mod_legend_col() + square_theme + large_text,
    align = "h", axis = "tb", nrow = 1, rel_widths = c(0.35, 0.33, 0.30))
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/04//merged_umap_overview-1.png)<!-- -->

Produce separate rasterized plots for publication.

<details>


```r
umap(seurat_joint, color_by = "Technology", alpha = 0.5, label = FALSE, point_size = 1,
     rasterize = T) +
    mod_legend_col() + 
    theme(text=element_text(size=30),
          legend.position = "bottom") +
    square_theme
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/04//technology_joint_umap-1.png)<!-- -->


```r
umap(seurat_joint, color_by = "Sample", colors = palette_sample, alpha = 0.95, label = FALSE, point_size = 1,
     rasterize = T) +
    mod_legend_col() + 
    theme(text=element_text(size=31),
          legend.position = "bottom") +
    square_theme
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/04//sample_joint_umap-1.png)<!-- -->


```r
umap(seurat_joint, legend = TRUE, alpha = 0.5, point_size = 1,
     rasterize = T) +
    no_legend() +
    theme(text=element_text(size=31)) +
    square_theme
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/04//clusters_joint_umap-1.png)<!-- -->

</details>

## QC of merged data

Plot QC metrics in the merged dataset.


```r
sc_meta <- seurat_joint@meta.data

# harmonize column names for scRNA/Multiome preprocessing
qc_rna <- sc_meta %>%
    filter(Technology == "10X Single Nuclei 3'") %>%
    select(Sample, Technology, nCount_RNA, nFeature_RNA, percent.mito, percent.ribo, G2M.Score, Harmony_cluster = seurat_clusters)

qc_multi <- sc_meta %>%
    filter(Technology == "10X Multiome") %>%
    select(Sample, Technology, nCount_RNA, nFeature_RNA, percent.mito, percent.ribo, G2M.Score,
           nCount_peaks, nFeature_peaks, atac_peak_region_fragments, nucleosome_signal, TSS.enrichment, Harmony_cluster = seurat_clusters)

all_multi <- bind_rows(qc_rna, qc_multi) %>%
    arrange(Technology) %>%
    mutate(Sample = factor(Sample, levels = unique(.$Sample))) %>%
    mutate(Harmony_cluster = factor(Harmony_cluster, levels = as.character(sort(as.numeric(levels(qc_multi$Harmony_cluster))))))


# create a function to plot the violin plots for each QC stat
plot_vln_stat_per_cluster <- function(stat, lims = c(0, NA), labels = FALSE) {

    # get the requested statistic
    df <- all_multi
    df$Stat <- df[[stat]]

    # ggplot violin plot
    p1 <- df %>% ggplot(aes(x = Harmony_cluster, y = Stat)) +
        geom_violin(aes(fill = Harmony_cluster), scale = "width") +
        # scale_fill_manual(values = palette_sample) +
        rotate_x() + no_legend() +
        # zoom in on y-axis if needed,
        # without changing distribution
        # and flip x/y axes
        coord_flip(ylim = lims) +
        xlab(NULL) +
        ggtitle(stat)

    # remove sample labels if needed
    # (after coord flip, samples are on the y-axis)
    if (!labels) p1 <- p1 + theme(axis.text.y = element_blank()) + ylab(NULL)

    p1 + large_text

}

stats <- replicate(5, c(0, NA), simplify = FALSE)

names(stats) <- c("nCount_RNA",
                  "nFeature_RNA",
                  "percent.mito",
                  "percent.ribo",
                  "G2M.Score")

stats$nCount_RNA <- c(0, 25000)
stats$nFeature_RNA <- c(0, 7500)
stats$G2M.Score <- c(-0.5, 1.5)

# number of cells in each cluster
p1 <- sc_meta %>%
    ggplot(aes(x = seurat_clusters)) +
    geom_bar(aes(fill = seurat_clusters)) +
    coord_flip() +
    no_legend() +
    ggtitle("# cells") +
    large_text

p2 <- sc_meta %>%
    ggplot(aes(x = seurat_clusters)) +
    geom_bar(aes(fill = Sample), position = "fill") +
    scale_fill_manual(values = palette_sample) +
    theme(panel.border = element_blank()) +
    rotate_x() +
    ggtitle("Cluster breakdown by sample") +
    coord_flip() +
    no_legend() + 
    large_text


# map over:
# 1. stat
# 2. y limits
# 3. whether to keep labels or not
vln_plots <- pmap(list(names(stats[c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo", "G2M.Score")]),
                       stats[c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo", "G2M.Score")],
                       rep(FALSE, 5)), ~ plot_vln_stat_per_cluster(..1, ..2, ..3))

plot_grid(plotlist = c(list(p1, p2), vln_plots),
          nrow = 1, align = "h", axis = "tb", rel_widths = c(0.2, 0.2, rep(0.15, 5)))
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/04//merged_qc_vln_per_cluster-1.png)<!-- -->



```r
plot_grid(
    umap(seurat_joint, "nCount_RNA",   color_by_type = "continuous", colors = rdbu, alpha = 0.5, label = FALSE, point_size = 0.5) + square_theme + large_text,
    umap(seurat_joint, "nFeature_RNA", color_by_type = "continuous", colors = rdbu, alpha = 0.5, label = FALSE, point_size = 0.5) + square_theme + large_text,
    umap(seurat_joint, "percent.mito", color_by_type = "continuous", colors = rdbu, alpha = 0.5, label = FALSE, point_size = 0.5) + square_theme + large_text,
    umap(seurat_joint, "percent.ribo", color_by_type = "continuous", colors = rdbu, alpha = 0.5, label = FALSE, point_size = 0.5) + square_theme + large_text,
    align = "hv", axis = "rltb")
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/04//merged_umap_qc-1.png)<!-- -->

## Normal vs. malignant calling

If a merged cluster contains >5% Normal cells called by inferCNV, we will label
the whole cluster as Normal.


```r
# get the proportion in each cluster assigned as normal or malignant
prop_category <- prop.table(table(Idents(seurat_joint), seurat_joint$inferCNV), margin = 1)
prop_category
```

```
##     
##         Malignant       Normal
##   0  0.9993148338 0.0006851662
##   1  0.9995025494 0.0004974506
##   2  1.0000000000 0.0000000000
##   3  0.9992619926 0.0007380074
##   4  1.0000000000 0.0000000000
##   5  1.0000000000 0.0000000000
##   6  0.9993442623 0.0006557377
##   7  0.9975649351 0.0024350649
##   8  1.0000000000 0.0000000000
##   9  0.5885608856 0.4114391144
##   10 1.0000000000 0.0000000000
##   11 0.6567567568 0.3432432432
##   12 0.9964664311 0.0035335689
##   13 0.9215686275 0.0784313725
##   14 0.9932432432 0.0067567568
##   15 0.6521739130 0.3478260870
```

```r
# if more than 5% of cells in one cluster are called Normal, call the cells
# in the cluster which weren't as Normal
seurat_joint$Malignant_normal <- map2_chr(
    Idents(seurat_joint), seurat_joint$inferCNV,
    function(cluster, infercnv) {

        if      (infercnv == "Malignant" & prop_category[cluster, "Normal"] > 0.05) "Normal"
        else infercnv

    })

table(seurat_joint$Malignant_normal)
```

```
## 
## Malignant    Normal 
##     33516      1204
```

```r
merged_sc_meta <- seurat_joint@meta.data
save(merged_sc_meta, file = glue("{out}/merged_sc_meta.Rda"))
```

Plot the malignant vs. normal cell labels.


```r
umap(seurat_joint,
     color_by = "Malignant_normal",
     colors = c("Normal" = "blue",
                "Malignant" = "pink"),
     rasterize = TRUE,
     point_size = 1,
     label = FALSE,
     legend = T,
     alpha = 0.5,
     title = glue("Malignant/normal"),
     hide_axes = T) +
    mod_legend_col() +
    theme(text=element_text(size=30),
          legend.position = "bottom") + square_theme
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/04//mal_norm_joint_umap-1.png)<!-- -->



```r
map(c("Malignant", "Normal"),
    function(i) {
        
        # https://github.com/satijalab/seurat/issues/3666#issuecomment-721307360
        cells <- names(which(seurat_joint$Malignant_normal == i))
       
        umap(seurat_joint,
             color_by = "Malignant_normal",
             colors = "#000000",
             cells = cells,
             rasterize = TRUE,
             point_size = 0.5,
             # put cells not selected in grey
             na_color = "grey80",
             label = FALSE,
             legend = FALSE,
             title = glue("{i} - {length(cells)} cells"),
             hide_axes = T) +
            theme(text=element_text(size=30))

    }) %>%
    {cowplot::plot_grid(plotlist = ., nrow = 1)}
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/04//mal_norm_facet_joint_umap-1.png)<!-- -->

# Infer NBFOXR2_6 patient sex

Metadata for patient sample `NBFOXR2_6` did not include patient sex.
Here, we will infer the biological sex of this patient based on expression of 
sex-specific genetic markers compared with the other five single-cell samples 
which include patient sex in metadata:

* `2236` - M
* `2273` - F
* `10153` - F
* `10154` - M
* `10155` - F

We will plot the following sex specific genes across single-cell samples:

* XIST (X chromosome)
* DDX3Y, ZFY, UTY, USP9Y (Y chromosome)

(Publication: [Olney et al. *Biology of Sex Differences* 2020](https://pubmed.ncbi.nlm.nih.gov/32693839/))


```r
VlnPlot(object = seurat_joint,
        features = c("XIST", "DDX3Y", "ZFY", "UTY", "USP9Y"), cols = palette_sample,
        group.by = "Sample",pt.size = 0)
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/04//sex-specific-gene-exp-vln-1.png)<!-- -->

Patient sample `NBFOXR2_6` more closely matches gene expression of patient samples 
assigned male in the metadata. Therefore we also assign this patient male in our metadata.

Based on this information, our single-cell cohort for NB-FOXR2 is balanced
for patient sex.

# TABLE: Human single-cell QC

Export a supplementary table for human single-cell RNA-seq QC.


```r
meta_sc_foxr2 <- meta_sc %>% filter(Group == "NB-FOXR2") %>%
    mutate(scRNAseq_path = gsub("data/", "", scRNAseq_path),
           Note = case_when(
               scRNAseq_path != "-" ~ "GEX from scRNAseq",
               scMultiome_path != "-" & ID != "P-6778_S-10155" ~ "Only GEX used from scMultiome",
               ID == "P-6778_S-10155" ~ "GEX and ATAC from scMultiome"
           ))

meta_sc_foxr2_scRNA_only <- meta_sc_foxr2 %>% filter(ID != "P-6778_S-10155")

# get sample prep info ---------------------------------------------------------
# get paths and IDs:
path_id_map <- meta_sc_foxr2 %>%
    select(ID, scRNAseq_path, scMultiome_path) %>%
    mutate(Path = ifelse(scRNAseq_path != "-", scRNAseq_path, scMultiome_path)) %>%
    select(ID, Path)

scRNA_omega <- read_xlsx(here("data/metadata/20240404-Omegatable-singlecell.xlsx")) %>%
    select(Sample,
           Path,
           Protocol,
           Publication,
           Kit = `kit version`,
           Method_of_dissociation = `Method of dissociation`,
           Starting_material = `Starting material`,
           Starting_material_amount = `tissue (mg)`,
           Targeted_N_cells = `cell/nuclei target`,
           Sequencing_platform = `Instrument`) %>%
    filter(Path %in% meta_sc_foxr2_scRNA_only$scRNAseq_path | Path %in% meta_sc_foxr2_scRNA_only$scMultiome_path) %>%
    left_join(path_id_map, by = "Path")
```

```
## New names:
## • `library yield (nM)` -> `library yield (nM)...36`
## • `copies/ul` -> `copies/ul...37`
## • `copies/ul` -> `copies/ul...38`
## • `library yield (nM)` -> `library yield (nM)...39`
```

```r
# get bioinformatics QC --------------------------------------------------------
scRNAseq_qc <- map2_dfr(meta_sc_foxr2_scRNA_only$ID, meta_sc_foxr2_scRNA_only$sc_preprocessing_path_hydra,
                        ~ data.table::fread(file.path(.y, "seurat_metrics.tsv"), data.table = FALSE) %>%
                            tibble::add_column(ID = .x, ID_omega = .y, .before = 1)) %>%
    left_join(scRNA_omega, by = "ID")

# get Cellranger QC ------------------------------------------------------------
scRNA_cellranger_qc <- map2_dfr(meta_sc_foxr2_scRNA_only$ID, meta_sc_foxr2_scRNA_only$Cellranger_path,
                                ~ data.table::fread(.y, data.table = FALSE) %>%
                                    mutate_all(as.character) %>%
                                    tibble::add_column(ID = .x, .before = 1) %>%
                                    tibble::add_column(Cellranger_path = .y, .after = 1))

# # get processing params --------------------------------------------------------
scRNA_config_params <- imap_dfr(samples_sc[meta_sc_foxr2_scRNA_only$ID], ~ .x@misc$params$config %>%
                                    # drop the parameter which indicates the genes plot during preprocessing
                                    discard_at("genes") %>%
                                    data.frame() %>%
                                    mutate(ID = .y))

scRNA_config_params2 <- scRNA_config_params %>%
    select(ID,
           Min_cells = min_cells,
           Min_features = min_features,
           Vars_to_regress = var_regress,
           N_principal_components = pcs_keep,
           Clustering_resolution = clustering_resolution,
           Seed = seed)

TABLE_scRNAseq_qc <- scRNAseq_qc %>%
    left_join(scRNA_cellranger_qc, by = c("ID")) %>%
    left_join(scRNA_config_params2, by = "ID") %>%
    left_join(meta_sc_foxr2_scRNA_only %>% select(ID, Note)) %>%
    select(
        # Sample info
        ID,
        # QC
        Note,
        Protocol,
        Sequencing_platform,
        Starting_material_amount,
        Method_of_dissociation,
        Targeted_N_cells,
        N_reads = `Number of Reads`,
        N_cells_estimated = `Estimated Number of Cells`,
        N_cells_after_filtering = N_cells_after,
        Reads_mapped_to_genome = `Reads Mapped to Genome`,
        Reads_mapped_to_transcriptome = `Reads Mapped Confidently to Transcriptome`,
        Min_cells, Min_features, Vars_to_regress, N_principal_components, Clustering_resolution, Seed,
        Mean_mitochondrial_content_after_filtering = percent.mito_mean.postQC,
        Mean_UMIs_after_filtering = nCount_RNA_mean.postQC,
        Mean_N_genes_after_filtering = nFeature_RNA_mean.postQC,
        Max_mito_threshold = percent.mito_max.threshold,
        Min_N_genes_threshold = nFeature_RNA_min.threshold,
        Max_N_genes_threshold = nFeature_RNA_max.threshold,
        Max_UMIs_threshold = nCount_RNA_max.threshold) %>%
    arrange(ID)
```

```
## Joining with `by = join_by(ID)`
```

```r
rr_write_tsv(TABLE_scRNAseq_qc,
             glue("{out}/TABLE_scRNAseq_QC.tsv"),
             "Summary of sample info and QC for scRNAseq of human tumor samples")
```

```
## ...writing description of TABLE_scRNAseq_QC.tsv to public/output/04/TABLE_scRNAseq_QC.desc
```

Create a similar table for the multiome sample:


```r
meta_sc_foxr2_Multiome_only <- meta_sc_foxr2 %>% filter(ID == "P-6778_S-10155") %>%
    mutate(Cellranger_path = file.path("/project/kleinman/singlecell_pipeline/data/",
                                       scMultiome_path,
                                       "cellranger-arc_count/summary.csv"))

# get sample prep info ---------------------------------------------------------
multiome_omega <- read_xlsx(here("data/metadata/20240404-Omegatable-singlecell.xlsx")) %>%
    select(Sample,
           Path,
           Protocol,
           Publication,
           Kit = `kit version`,
           Method_of_dissociation = `Method of dissociation`,
           Starting_material = `Starting material`,
           Starting_material_amount = `tissue (mg)`,
           Targeted_N_cells = `cell/nuclei target`,
           Sequencing_platform = `Instrument`) %>%
    filter(Path %in% meta_sc_foxr2_Multiome_only$scMultiome_path) %>%
    left_join(path_id_map, by = "Path")
```

```
## New names:
## • `library yield (nM)` -> `library yield (nM)...36`
## • `copies/ul` -> `copies/ul...37`
## • `copies/ul` -> `copies/ul...38`
## • `library yield (nM)` -> `library yield (nM)...39`
```

```r
# get bioinformatics QC --------------------------------------------------------
multiome_qc <- data.table::fread(file.path(meta_sc_foxr2_Multiome_only$sc_preprocessing_path_hydra, "seurat_metrics.tsv"), data.table = FALSE) %>%
    tibble::add_column(ID = meta_sc_foxr2_Multiome_only$ID,
                       ID_omega = meta_sc_foxr2_Multiome_only$sc_preprocessing_path_hydra,
                       .before = 1) %>%
    left_join(multiome_omega, by = "ID")

# get Cellranger QC ------------------------------------------------------------
multiome_cellranger_qc <- data.table::fread(meta_sc_foxr2_Multiome_only$Cellranger_path, data.table = FALSE) %>%
    mutate_all(as.character) %>%
    tibble::add_column(ID = meta_sc_foxr2_Multiome_only$ID, .before = 1) %>%
    tibble::add_column(Cellranger_path = meta_sc_foxr2_Multiome_only$Cellranger_path, .after = 1)

# # get processing params --------------------------------------------------------
multiome_config_params <- samples_sc$`P-6778_S-10155`@misc$params$config %>%
    discard_at("genes") %>%
    data.frame() %>%
    mutate(ID = "P-6778_S-10155")

multiome_config_params2 <- multiome_config_params %>%
    select(ID,
           Min_cells = min_cells,
           N_principal_components_RNA = rna_pcs_keep,
           N_principal_components_ATAC = atac_pcs_keep,
           Clustering_resolution = clustering_resolution,
           Seed = seed)

TABLE_multiome_qc <- multiome_qc %>%
    left_join(multiome_cellranger_qc, by = c("ID")) %>%
    left_join(multiome_config_params2, by = "ID") %>%
    select(
        # Sample info
        ID,
        # QC
        Protocol,
        Sequencing_platform,
        Starting_material_amount,
        Method_of_dissociation,
        Targeted_N_cells,
        Genome,
        Cellranger_version = `Pipeline version`,
        N_cells_estimated = `Estimated number of cells`,
        N_cells_after_filtering = N_cells_after,
        GEX_Reads_mapped_to_genome = `GEX Reads mapped confidently to genome`,
        GEX_Reads_mapped_to_transcriptome = `GEX Reads mapped confidently to transcriptome`,
        GEX_Mean_mitochondrial_content_after_filtering = percent.mito_mean.postQC,
        GEX_Mean_UMIs_after_filtering = nCount_RNA_mean.postQC,
        GEX_Mean_N_genes_after_filtering = nFeature_RNA_mean.postQC,
        GEX_Max_mito_threshold = percent.mito_max.threshold,
        GEX_Min_N_genes_threshold = nFeature_RNA_min.threshold,
        GEX_Max_N_genes_threshold = nFeature_RNA_max.threshold,
        GEX_Min_UMIs_threshold = nCount_RNA_min.threshold,
        GEX_Max_UMIs_threshold = nCount_RNA_max.threshold,
        ATAC_Fraction_transposition_events_in_peaks_in_cells = `ATAC Fraction of transposition events in peaks in cells`,
        ATAC_Fraction_mapped_confidently = `ATAC Confidently mapped read pairs`,
        ATAC_Fraction_fragments_in_peaks = `ATAC Fraction of high-quality fragments overlapping peaks`,
        Min_cells, Clustering_resolution, N_principal_components_RNA, N_principal_components_ATAC, Seed,
        ATAC_Median_fragments_per_cell = `ATAC Median high-quality fragments per cell`,
        ATAC_Mean_nucleosome_signal_after_filtering = nucleosome_signal_mean.postQC,
        ATAC_Mean_peak_region_fragments_after_filtering = atac_peak_region_fragments_mean.postQC,
        ATAC_Mean_TSS_enrichment_after_filtering = TSS.enrichment_mean.postQC,
        ATAC_Min_peak_region_fragments_threshold = atac_peak_region_fragments_min.threshold,
        ATAC_Max_peak_region_fragments_threshold = atac_peak_region_fragments_max.threshold,
        ATAC_Max_nucleosome_signal_threshold = nucleosome_signal_max.threshold,
        ATAC_Min_TSS_enrichment_threshold = TSS.enrichment_min.threshold,
        ATAC_Max_TSS_enrichment_threshold = TSS.enrichment_max.threshold)

rr_write_tsv(TABLE_multiome_qc,
             glue("{out}/TABLE_scMultiome_QC.tsv"),
             "Summary of sample info and QC for scMultiome of human tumor sample (n=1)")
```

```
## ...writing description of TABLE_scMultiome_QC.tsv to public/output/04/TABLE_scMultiome_QC.desc
```


<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2024-11-01 11:51:38
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
##  date     2024-11-01
##  pandoc   1.19.2.1 @ /cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/ (via rmarkdown)
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  ! package          * version  date (UTC) lib source
##  P abind              1.4-5    2016-07-21 [?] CRAN (R 4.1.2)
##  P beeswarm           0.4.0    2021-06-01 [?] CRAN (R 4.1.2)
##  P BiocGenerics       0.40.0   2021-10-26 [?] Bioconductor
##  P BiocManager        1.30.15  2021-05-11 [?] CRAN (R 4.1.2)
##  P BiocParallel       1.28.3   2021-12-09 [?] Bioconductor
##  P Biostrings         2.62.0   2021-10-26 [?] Bioconductor
##  P bit                4.0.4    2020-08-04 [?] CRAN (R 4.1.2)
##  P bit64              4.0.5    2020-08-30 [?] CRAN (R 4.1.2)
##  P bitops             1.0-7    2021-04-24 [?] CRAN (R 4.1.2)
##  P bslib              0.3.1    2021-10-06 [?] CRAN (R 4.1.2)
##  P cachem             1.0.6    2021-08-19 [?] CRAN (R 4.1.2)
##  P Cairo              1.6-0    2022-07-05 [?] CRAN (R 4.1.2)
##  P callr              3.7.6    2024-03-25 [?] RSPM
##  P cellranger         1.1.0    2016-07-27 [?] CRAN (R 4.1.2)
##  P cli                3.6.1    2023-03-23 [?] RSPM (R 4.1.2)
##  P cluster            2.1.2    2021-04-17 [?] CRAN (R 4.1.2)
##  P codetools          0.2-18   2020-11-04 [?] CRAN (R 4.1.2)
##  P colorspace         2.0-2    2021-06-24 [?] CRAN (R 4.1.2)
##  P cowplot          * 1.1.1    2020-12-30 [?] CRAN (R 4.1.2)
##  P crayon             1.4.2    2021-10-29 [?] CRAN (R 4.1.2)
##  P data.table         1.14.2   2021-09-27 [?] CRAN (R 4.1.2)
##  P deldir             1.0-6    2021-10-23 [?] CRAN (R 4.1.2)
##  P devtools           2.4.5    2022-10-11 [?] CRAN (R 4.1.2)
##  P digest             0.6.35   2024-03-11 [?] CRAN (R 4.1.2)
##  P docopt             0.7.1    2020-06-24 [?] CRAN (R 4.1.2)
##  P dplyr            * 1.1.1    2023-03-22 [?] CRAN (R 4.1.2)
##  P ellipsis           0.3.2    2021-04-29 [?] CRAN (R 4.1.2)
##  P evaluate           0.23     2023-11-01 [?] CRAN (R 4.1.2)
##  P fansi              1.0.2    2022-01-14 [?] CRAN (R 4.1.2)
##  P farver             2.1.0    2021-02-28 [?] CRAN (R 4.1.2)
##  P fastmap            1.1.0    2021-01-25 [?] CRAN (R 4.1.2)
##  P fastmatch          1.1-3    2021-07-23 [?] CRAN (R 4.1.2)
##  P fitdistrplus       1.1-6    2021-09-28 [?] CRAN (R 4.1.2)
##  P fs                 1.5.2    2021-12-08 [?] CRAN (R 4.1.2)
##  P future             1.25.0   2022-04-24 [?] CRAN (R 4.1.2)
##  P future.apply       1.8.1    2021-08-10 [?] CRAN (R 4.1.2)
##  P generics           0.1.3    2022-07-05 [?] CRAN (R 4.1.2)
##  P GenomeInfoDb       1.30.1   2022-01-30 [?] Bioconductor
##  P GenomeInfoDbData   1.2.4    2023-11-28 [?] Bioconductor
##  P GenomicRanges      1.46.1   2021-11-18 [?] Bioconductor
##  P ggbeeswarm         0.7.1    2022-12-16 [?] CRAN (R 4.1.2)
##  P ggforce            0.3.3    2021-03-05 [?] CRAN (R 4.1.2)
##  P ggplot2          * 3.4.2    2023-04-03 [?] CRAN (R 4.1.2)
##  P ggrastr            1.0.1    2021-12-08 [?] CRAN (R 4.1.2)
##  P ggrepel            0.9.1    2021-01-15 [?] CRAN (R 4.1.2)
##  P ggridges           0.5.3    2021-01-08 [?] CRAN (R 4.1.2)
##  P ggseqlogo          0.1      2017-07-25 [?] CRAN (R 4.1.2)
##  P git2r              0.29.0   2021-11-22 [?] CRAN (R 4.1.2)
##  P globals            0.14.0   2020-11-22 [?] CRAN (R 4.1.2)
##  P glue             * 1.6.2    2022-02-24 [?] CRAN (R 4.1.2)
##  P goftest            1.2-3    2021-10-07 [?] CRAN (R 4.1.2)
##  P gridExtra          2.3      2017-09-09 [?] CRAN (R 4.1.2)
##  P gtable             0.3.0    2019-03-25 [?] CRAN (R 4.1.2)
##  P here             * 1.0.1    2020-12-13 [?] CRAN (R 4.1.2)
##  P highr              0.9      2021-04-16 [?] CRAN (R 4.1.2)
##  P hms                1.1.1    2021-09-26 [?] CRAN (R 4.1.2)
##  P htmltools          0.5.2    2021-08-25 [?] CRAN (R 4.1.2)
##  P htmlwidgets        1.5.4    2021-09-08 [?] CRAN (R 4.1.2)
##  P httpuv             1.6.5    2022-01-05 [?] CRAN (R 4.1.2)
##  P httr               1.4.2    2020-07-20 [?] CRAN (R 4.1.2)
##  P ica                1.0-2    2018-05-24 [?] CRAN (R 4.1.2)
##  P igraph             2.0.3    2024-03-13 [?] CRAN (R 4.1.2)
##  P IRanges            2.28.0   2021-10-26 [?] Bioconductor
##  P irlba              2.3.5    2021-12-06 [?] CRAN (R 4.1.2)
##  P jquerylib          0.1.4    2021-04-26 [?] CRAN (R 4.1.2)
##  P jsonlite           1.8.8    2023-12-04 [?] CRAN (R 4.1.2)
##  P KernSmooth         2.23-20  2021-05-03 [?] CRAN (R 4.1.2)
##  P knitr              1.37     2021-12-16 [?] CRAN (R 4.1.2)
##  P labeling           0.4.2    2020-10-20 [?] CRAN (R 4.1.2)
##  P later              1.3.0    2021-08-18 [?] CRAN (R 4.1.2)
##  P lattice            0.20-45  2021-09-22 [?] CRAN (R 4.1.2)
##  P lazyeval           0.2.2    2019-03-15 [?] CRAN (R 4.1.2)
##  P leiden             0.3.9    2021-07-27 [?] CRAN (R 4.1.2)
##  P lifecycle          1.0.3    2022-10-07 [?] CRAN (R 4.1.2)
##  P listenv            0.8.0    2019-12-05 [?] CRAN (R 4.1.2)
##  P lmtest             0.9-39   2021-11-07 [?] CRAN (R 4.1.2)
##  P lsa                0.73.2   2020-05-04 [?] CRAN (R 4.1.2)
##  P magrittr         * 2.0.3    2022-03-30 [?] CRAN (R 4.1.2)
##  P MASS               7.3-54   2021-05-03 [?] CRAN (R 4.1.2)
##  P Matrix             1.3-4    2021-06-01 [?] CRAN (R 4.1.2)
##  P matrixStats        0.61.0   2021-09-17 [?] CRAN (R 4.1.2)
##  P memoise            2.0.1    2021-11-26 [?] CRAN (R 4.1.2)
##  P mgcv               1.8-38   2021-10-06 [?] CRAN (R 4.1.2)
##  P mime               0.12     2021-09-28 [?] CRAN (R 4.1.2)
##  P miniUI             0.1.1.1  2018-05-18 [?] CRAN (R 4.1.2)
##  P munsell            0.5.0    2018-06-12 [?] CRAN (R 4.1.2)
##  P nlme               3.1-153  2021-09-07 [?] CRAN (R 4.1.2)
##  P parallelly         1.30.0   2021-12-17 [?] CRAN (R 4.1.2)
##  P patchwork          1.1.1    2020-12-17 [?] CRAN (R 4.1.2)
##  P pbapply          * 1.5-0    2021-09-16 [?] CRAN (R 4.1.2)
##  P pillar             1.9.0    2023-03-22 [?] RSPM (R 4.1.2)
##  P pkgbuild           1.4.2    2023-06-26 [?] CRAN (R 4.1.2)
##  P pkgconfig          2.0.3    2019-09-22 [?] CRAN (R 4.1.2)
##  P pkgload            1.3.3    2023-09-22 [?] CRAN (R 4.1.2)
##  P plotly             4.10.0   2021-10-09 [?] CRAN (R 4.1.2)
##  P plyr               1.8.6    2020-03-03 [?] CRAN (R 4.1.2)
##  P png                0.1-7    2013-12-03 [?] CRAN (R 4.1.2)
##  P polyclip           1.10-0   2019-03-14 [?] CRAN (R 4.1.2)
##  P prettyunits        1.1.1    2020-01-24 [?] CRAN (R 4.1.2)
##  P processx           3.8.4    2024-03-16 [?] RSPM
##  P profvis            0.3.8    2023-05-02 [?] CRAN (R 4.1.2)
##  P promises           1.2.0.1  2021-02-11 [?] CRAN (R 4.1.2)
##  P ps                 1.7.6    2024-01-18 [?] RSPM
##  P purrr            * 1.0.1    2023-01-10 [?] CRAN (R 4.1.2)
##  P qlcMatrix          0.9.7    2018-04-20 [?] CRAN (R 4.1.2)
##  P R6                 2.5.1    2021-08-19 [?] CRAN (R 4.1.2)
##  P RANN               2.6.1    2019-01-08 [?] CRAN (R 4.1.2)
##  P RColorBrewer     * 1.1-2    2014-12-07 [?] CRAN (R 4.1.2)
##  P Rcpp               1.0.8    2022-01-13 [?] CRAN (R 4.1.2)
##  P RcppAnnoy          0.0.19   2021-07-30 [?] CRAN (R 4.1.2)
##  P RcppRoll           0.3.0    2018-06-05 [?] CRAN (R 4.1.2)
##  P RCurl              1.98-1.5 2021-09-17 [?] CRAN (R 4.1.2)
##  P readr            * 2.1.1    2021-11-30 [?] CRAN (R 4.1.2)
##  P readxl           * 1.3.1    2019-03-13 [?] CRAN (R 4.1.2)
##  P remotes            2.4.2.1  2023-07-18 [?] CRAN (R 4.1.2)
##  P renv               1.0.3    2023-09-19 [?] CRAN (R 4.1.2)
##  P reshape2           1.4.4    2020-04-09 [?] CRAN (R 4.1.2)
##  P reticulate         1.23     2022-01-14 [?] CRAN (R 4.1.2)
##  P rlang              1.1.3    2024-01-10 [?] CRAN (R 4.1.2)
##  P rmarkdown          2.11     2021-09-14 [?] CRAN (R 4.1.2)
##  P ROCR               1.0-11   2020-05-02 [?] CRAN (R 4.1.2)
##  P rpart              4.1-15   2019-04-12 [?] CRAN (R 4.1.2)
##  P rprojroot          2.0.2    2020-11-15 [?] CRAN (R 4.1.2)
##  P Rsamtools          2.10.0   2021-10-26 [?] Bioconductor
##  P Rtsne              0.15     2018-11-10 [?] CRAN (R 4.1.2)
##  P S4Vectors          0.32.4   2022-03-24 [?] Bioconductor
##  P sass               0.4.0    2021-05-12 [?] CRAN (R 4.1.2)
##  P scales             1.2.1    2022-08-20 [?] CRAN (R 4.1.2)
##  P scattermore        0.7      2020-11-24 [?] CRAN (R 4.1.2)
##  P sctransform        0.3.3    2022-01-13 [?] CRAN (R 4.1.2)
##  P sessioninfo        1.2.2    2021-12-06 [?] CRAN (R 4.1.2)
##  P Seurat           * 4.0.0    2021-01-30 [?] CRAN (R 4.1.2)
##  P SeuratObject     * 4.0.4    2021-11-23 [?] CRAN (R 4.1.2)
##  P shiny              1.7.1    2021-10-02 [?] CRAN (R 4.1.2)
##  P Signac           * 1.3.0    2021-07-12 [?] CRAN (R 4.1.2)
##  P slam               0.1-50   2022-01-08 [?] CRAN (R 4.1.2)
##  P SnowballC          0.7.0    2020-04-01 [?] CRAN (R 4.1.2)
##  P sparsesvd          0.2      2019-07-15 [?] CRAN (R 4.1.2)
##  P spatstat           1.64-1   2020-05-12 [?] CRAN (R 4.1.2)
##  P spatstat.data      2.1-2    2021-12-17 [?] CRAN (R 4.1.2)
##  P spatstat.utils     2.3-0    2021-12-12 [?] CRAN (R 4.1.2)
##  P stringi            1.7.6    2021-11-29 [?] CRAN (R 4.1.2)
##  P stringr          * 1.5.0    2022-12-02 [?] CRAN (R 4.1.2)
##  P survival           3.2-13   2021-08-24 [?] CRAN (R 4.1.2)
##  P tensor             1.5      2012-05-05 [?] CRAN (R 4.1.2)
##  P tibble             3.2.1    2023-03-20 [?] RSPM (R 4.1.2)
##  P tidyr            * 1.3.0    2023-01-24 [?] CRAN (R 4.1.2)
##  P tidyselect         1.2.0    2022-10-10 [?] CRAN (R 4.1.2)
##  P tweenr             1.0.2    2021-03-23 [?] CRAN (R 4.1.2)
##  P tzdb               0.3.0    2022-03-28 [?] CRAN (R 4.1.2)
##  P urlchecker         1.0.1    2021-11-30 [?] CRAN (R 4.1.2)
##  P usethis            2.2.2    2023-07-06 [?] CRAN (R 4.1.2)
##  P utf8               1.2.2    2021-07-24 [?] CRAN (R 4.1.2)
##  P uwot               0.1.11   2021-12-02 [?] CRAN (R 4.1.2)
##  P vctrs              0.6.5    2023-12-01 [?] CRAN (R 4.1.2)
##  P vipor              0.4.5    2017-03-22 [?] CRAN (R 4.1.2)
##  P viridis          * 0.5.1    2018-03-29 [?] RSPM (R 4.1.2)
##  P viridisLite      * 0.3.0    2018-02-01 [?] CRAN (R 4.1.2)
##  P vroom              1.5.7    2021-11-30 [?] CRAN (R 4.1.2)
##  P withr              2.5.0    2022-03-03 [?] CRAN (R 4.1.2)
##  P xfun               0.29     2021-12-14 [?] CRAN (R 4.1.2)
##  P xtable             1.8-4    2019-04-21 [?] CRAN (R 4.1.2)
##  P XVector            0.34.0   2021-10-26 [?] Bioconductor
##  P yaml               2.2.1    2020-02-01 [?] CRAN (R 4.1.2)
##  P zlibbioc           1.40.0   2021-10-26 [?] Bioconductor
##  P zoo                1.8-9    2021-03-09 [?] CRAN (R 4.1.2)
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
