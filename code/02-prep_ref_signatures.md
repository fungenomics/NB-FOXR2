---
title: "02 - Prep normal reference signatures"
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
## Document index: 02
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
## public/output/02
```

```
## public/figures/02//
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

Here we'll assemble a set of cell type signatures from the developing brain / adult brain,
building off the collection assembled for our
[Chen et al, Cell, 2020](https://www.cell.com/cell/fulltext/S0092-8674(20)31529-4)
study on G34R/V gliomas.

I will also add additional signatures for other brain datasets.

_NOTE_: This document re-uses a few functions for processing markers into lists
of gene signatures, `filter_markers()` and `signatures_to_list()` which are located
in `code/functions/scRNAseq.R`.


# Libraries


```r
library(here)
library(magrittr)
library(tidyr)
library(dplyr)
library(readr)
library(readxl)
library(glue)
library(purrr)
library(ggplot2)
library(pheatmap)
library(Seurat)

source(here("include/style.R"))
source(here("code/functions/scRNAseq.R"))

ggplot2::theme_set(theme_min())

conflicted::conflicts_prefer(dplyr::filter)
```



# Gather signatures

Set function for conversion ensembl <-> symbol.


```r
signatures_to_list <- function(signatures, type = "ens") {
    
    if (type == "ens") {
        
        signatures %>% 
            select(Cluster, gene) %>%
            split(.$Cluster) %>%
            map(~ pull(.x, gene) %>% symbols2ensembl())
        
    } else if (type == "sym") {
        
        signatures %>% 
            select(Cluster, gene) %>%
            split(.$Cluster) %>%
            map(~ pull(.x, gene))
    }
}

symbols2ensembl <- function(genes) {
    
    load(here("data/singlecell/references_genome/biomaRt_hg_symbol_to_ens_lds.Rda"))
    genes_ens_lds %>% filter(HGNC.symbol %in% genes) %>% pull(Gene.stable.ID)
    
}
```

## Signatures compiled for Chen et al, 2020

In our analysis of H3G34R/V gliomas, we assembled signatures
from several studies on pre and postnatal forebrain, which we will re-use here ([GitHub link](https://fungenomics.github.io/G34-gliomas/bulk_transcriptome_epigenome/analysis/02-GSEA.html)). Copying over the comments from the Chen et al analysis:

> We obtained reference datasets for the forebrain across species and the lifespan:
Jessa et al, 2019; Nowakowski et al, 2017; Velmeshev et al, 2019
We also obtained two reference datasets capturing the adult SVZ:
Anderson et al, 2020; Mizrak et al, 2019



```r
# laod gene signatures
signatures_ens <- readRDS(here("data/singlecell/references_normal/Chen_Cell_2020/signatures_ens.Rds"))
signatures_sym <- readRDS(here("data/singlecell/references_normal/Chen_Cell_2020//signatures_sym.Rds"))

# load signature annotation
cell_type_anno <- read_tsv(here("data/singlecell/references_normal/Chen_Cell_2020/signature_annotations.tsv"))
```

```
## Rows: 340 Columns: 5
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (5): Sample, Age, Cell_type, Cluster, Dataset
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

This set contains signatures from the our mouse forebrain atlas published in Jessa et al,
Nature Genetics, 2019. Since we have recently extended the timepoints covered by this atlas (Jessa et al, Nature Genetics, 2022), we will remove all the 2019 signatures from this collection. In the next section, we load in both 2019 & 2022 signatures.


```r
cell_type_anno_filt <- cell_type_anno %>% 
    filter(Dataset != "Jessa et al 2019")

signatures_ens <- signatures_ens[cell_type_anno_filt$Cluster]
signatures_sym <- signatures_sym[cell_type_anno_filt$Cluster]
```

## Jessa et al, Nature Genetics, 2022

Load the signatures for the extended atlas, combining timepoints published in 2019 and 2022:


```r
# load signatures
jessa_signatures <- readRDS(here("data/singlecell/references_normal/Jessa_NatGenet_2022/joint_mouse_extended.signatures_ID_20210710.Rds"))

# subset to forebrain
jessa_signatures_ens <- jessa_signatures$hg_ens[grepl("^F-", names(jessa_signatures$hg_ens))]
jessa_signatures_sym <- jessa_signatures$hg_sym[grepl("^F-", names(jessa_signatures$hg_sym))]

# load annotation
jessa_cell_type_anno <- read_tsv(here("data/singlecell/references_normal/Jessa_NatGenet_2022/metadata_20210710_with_qc_v2.tsv")) %>%
    filter(Location == "Forebrain") %>% 
    filter(!grepl("EXCLUDE", Label)) %>% 
    filter(Label %in% names(jessa_signatures_ens)) %>% 
    select(Cluster = Label,
           Age = Timepoint,
           Cell_type = Level4_type) %>% 
    mutate(Dataset = "Jessa et al 2022",
           Sample = "Mouse forebrain")
```

```
## Rows: 346 Columns: 22
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (15): Alias, Sample, Location, Timepoint, Colour, Cell_ontological_class...
## dbl  (7): Cluster, N_cells, median_logFC, max_logFC, n_gene_FC_above_1.5, me...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
# sanity check
nrow(jessa_cell_type_anno) == length(jessa_signatures_ens)
```

```
## [1] TRUE
```

```r
nrow(jessa_cell_type_anno) == length(jessa_signatures_sym)
```

```
## [1] TRUE
```


## Aldinger et al, Nat Neurosci, 2021


```r
# load cell annotaiton
aldinger_anno <- read_xlsx(here("data/singlecell/references_normal/Aldinger_NatNeurosci_2021/data/Aldinger_SupplementaryTables.xlsx"), sheet = 8, skip = 2) %>% 
    dplyr::rename(Cluster = `Cluster-Cell Type Abbreviation`, Cell_type = `Cell Type`)

# load cluster markers
load(here("data/singlecell/references_normal/Aldinger_NatNeurosci_2021/aldinger_markers.Rda"))

# convert to signatures
aldinger_signatures <- aldinger_markers %>% 
    mutate(Cluster = paste0("Human fetal cerebellum ", cluster)) %>% 
    dplyr::rename(avg_logFC = avg_log2FC) %>% 
    filter_markers(n_top = 100)

sort(table(aldinger_signatures$Cluster))
```

```
## 
##              Human fetal cerebellum 20-Choroid 
##                                             31 
##              Human fetal cerebellum 05-eCN/UBC 
##                                             98 
##                   Human fetal cerebellum 01-PC 
##                                            100 
##                   Human fetal cerebellum 02-RL 
##                                            100 
##                  Human fetal cerebellum 03-GCP 
##                                            100 
##                   Human fetal cerebellum 04-GN 
##                                            100 
##                  Human fetal cerebellum 06-iCN 
##                                            100 
##                  Human fetal cerebellum 07-PIP 
##                                            100 
##                   Human fetal cerebellum 08-BG 
##                                            100 
##                  Human fetal cerebellum 09-Ast 
##                                            100 
##                 Human fetal cerebellum 10-Glia 
##                                            100 
##                  Human fetal cerebellum 11-OPC 
##                                            100 
##        Human fetal cerebellum 12-Committed OPC 
##                                            100 
##          Human fetal cerebellum 13-Endothelial 
##                                            100 
##            Human fetal cerebellum 14-Microglia 
##                                            100 
##             Human fetal cerebellum 15-Meninges 
##                                            100 
##            Human fetal cerebellum 16-Pericytes 
##                                            100 
##            Human fetal cerebellum 17-Brainstem 
##                                            100 
##                  Human fetal cerebellum 18-MLI 
##                                            100 
##        Human fetal cerebellum 19-Ast/Ependymal 
##                                            100 
## Human fetal cerebellum 21-BS Choroid/Ependymal 
##                                            100
```

```r
aldinger_signatures_ens <- aldinger_signatures %>%
    signatures_to_list("ens")

aldinger_signatures_sym <- aldinger_signatures %>%
    signatures_to_list("sym")

aldinger_cell_type_anno  <- aldinger_anno %>%
    select(-3, -4, -5) %>% 
    mutate(Cluster = paste0("Human fetal cerebellum ", Cluster)) %>% 
    mutate(Dataset = "Aldinger et al 2021",
           Sample = "Human fetal cerebellum",
           Age = "9-21 PCW")
```


## Dong et al, Cancer Cell

This study profiled human fetal embryos and adrenal glands as a reference for 
studies of extra-cranial neuroblastoma.


```r
dong_embryo_markers <- read_xlsx(here("data/singlecell/references_normal/Dong_CancerCell_2020/Dong_TableS3.xlsx"), skip = 2) %>%
    dplyr::rename(Cluster = `Cell  type`, gene = Gene)

dong_sympatho_markers <- read_xlsx(here("data/singlecell/references_normal/Dong_CancerCell_2020/Dong_TableS4.xlsx"), skip = 2) %>%
    dplyr::rename(Cluster = Cell_type, gene = Gene)

dong_signatures <- bind_rows(dong_embryo_markers %>% mutate(Cluster = paste0("Human 4PCW ", Cluster)),
                             dong_sympatho_markers %>% mutate(Cluster = paste0("Human fetal adrenal ", Cluster))) %>%
    # take care of a common cluster name between datasets
    mutate(Cluster = ifelse(Cluster == "Human fetal adrenal SCPs", "Human fetal adrenal SCP", Cluster)) %>% 
    filter_markers(sp = "hg", n_top = 100)

table(dong_signatures$Cluster)
```

```
## 
##                              Human 4PCW EMT NCCs 
##                                              100 
##                           Human 4PCW Floor plate 
##                                              100 
##                  Human 4PCW Forebrain Primordium 
##                                              100 
##                               Human 4PCW INs Pro 
##                                              100 
##                           Human 4PCW Melanoblast 
##                                              100 
##                                 Human 4PCW MNs I 
##                                              100 
##                                Human 4PCW MNs II 
##                                              100 
##                               Human 4PCW MNs III 
##                                              100 
##                               Human 4PCW MNs Pro 
##                                              100 
##                                  Human 4PCW NCCs 
##                                              100 
##                                    Human 4PCW NT 
##                                              100 
##                        Human 4PCW Pre-EMT NCCs I 
##                                              100 
##                       Human 4PCW Pre-EMT NCCs II 
##                                              100 
##                                   Human 4PCW SNs 
##                                              100 
##     Human fetal adrenal Cycling chromaffin cells 
##                                              100 
##       Human fetal adrenal Cycling sympathoblasts 
##                                              100 
## Human fetal adrenal Non-cycling chromaffin cells 
##                                              100 
##   Human fetal adrenal Non-cycling sympathoblasts 
##                                              100 
##                          Human fetal adrenal SCP 
##                                              100
```

```r
dong_signatures_ens <- dong_signatures %>% signatures_to_list("ens")
dong_signatures_sym <- dong_signatures %>% signatures_to_list("sym")

dong_cell_type_anno <- dong_signatures %>%
    distinct(Cluster) %>%
    mutate(Dataset = "Dong et al 2020",
           Age = "Human fetal",
           Cell_type = Cluster,
           Sample = case_when(
               grepl("4PCW", Cluster) ~ "Human 4PCW embryo neural lineages",
               grepl("adrenal", Cluster) ~ "Human 8-14PCW fetal adrenal glands"
           ))
```



## Jansky et al, Nature Genetics, 2021

This study also profiled fetal human adrenal glands and adrenal medulla for 
comparisons/study of neuroblastoma.


```r
jansky_4 <- read_xlsx(here("data/singlecell/references_normal/Jansky_NatGenet_2021/Jansky_et_al_SupplementaryTables1-11.xlsx"), sheet = "Suppl. Table 4", skip = 1)
jansky_5 <- read_xlsx(here("data/singlecell/references_normal/Jansky_NatGenet_2021/Jansky_et_al_SupplementaryTables1-11.xlsx"), sheet = "Suppl. Table 5", skip = 1)

jansky_signatures <- bind_rows(jansky_4, jansky_5) %>% 
    dplyr::rename(Cluster = cluster,
                  avg_logFC = `average logFC`) %>% 
    filter_markers() %>%
    mutate(Cluster = paste0("Human fetal adrenal2 ", Cluster))

jansky_signatures_ens <- jansky_signatures %>% signatures_to_list("ens")
jansky_signatures_sym <- jansky_signatures %>% signatures_to_list("sym")

jansky_cell_type_anno <- jansky_signatures %>%
    distinct(Cluster) %>%
    mutate(Dataset = "Jansky et al 2021",
           Age = "Human fetal 7-17PCW",
           Cell_type = Cluster,
           Sample = "Human fetal adrenal gland")
```


## Kildisiute et al, Science Advances, 2021

This study also profiled human fetal adrenal glands as a reference for 
studies of extra-cranial neuroblastoma.

These markers were called with a custom algorithm by the authors, as described here
by the authors in their Methods:

> Using well-established marker genes of different cell types curated  from the literature (table S2), we assigned a cell type to each cluster.  Where two clusters were annotated as the same cell type, we merged  them together. As further confirmation of our annotation, we next  identified marker genes for each population algorithmically (table S6).  To do this, we used a method that uses the tf-idf metric to identify  genes specific to each population, as implemented in the “quickMarkers”  function in the SoupX R package (29). We further filtered genes to  include only those genes with a P value less than 0.01 after multiple  hypothesis correction (hypergeometric test).

Let's plot the markers produced in this study, comparing for each gene, the TF-IDF metric with the difference in detection rate of the gene in cells inside the cluster compared to all other cells. Each point is one gene.


```r
kildisiute_markers <- read_xlsx(here("data/singlecell/references_normal/Kildisiute_SciAdv_2021/Kildisiute_Tables_S1_to_S12.xlsx"),
                                sheet = "TableS6.algorithmic_markers") %>%
    # filter out tumor datasets
    filter(dataset %in% c("adrenal_gland", "stringent_adrenal")) %>%
    mutate(pct_diff = geneFrequency - geneFrequencyOutsideCluster) %>%
    dplyr::rename(Cluster = cluster) %>%
    filter(Cluster != "Other") %>%
    mutate(Cluster = paste0("Human fetal adrenal ", Cluster))

# quite a large number of markers
table(kildisiute_markers$Cluster, kildisiute_markers$dataset)
```

```
##                                           
##                                            adrenal_gland stringent_adrenal
##   Human fetal adrenal Bridge                        6877                25
##   Human fetal adrenal Chromaffin                    9324                31
##   Human fetal adrenal Cortex                       13593                 0
##   Human fetal adrenal Erythroblasts                  441                 0
##   Human fetal adrenal Leukocytes                    4143                 0
##   Human fetal adrenal Mesenchyme                    5441                 0
##   Human fetal adrenal Podocytes                        0               176
##   Human fetal adrenal SCPs                          1887                57
##   Human fetal adrenal Sympathoblastic               8146                59
##   Human fetal adrenal Vascular Endothelium          2734                 0
```

```r
# explore the dataset
kildisiute_markers %>%
    ggplot(aes(x = pct_diff, y = tfidf)) +
    geom_point(aes(colour = Cluster), alpha = 0.3) +
    geom_smooth(aes(colour = Cluster), method = "lm", se = FALSE) +
    xlab("Diff. in detection rate inside and outside cluster") +
    ylab("TF-IDF metric") +
    ggtitle("Before filtering")
```

```
## `geom_smooth()` using formula = 'y ~ x'
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/02//inspect_kildisiute-1.png)<!-- -->

This tells us that a high TF-IDF metric correlates with high/specific detection rate.
However, for some clusters, there's more variability along the X-axis, so I'll select markers based on the difference in detection rates (`pct_diff`), rather than the TF-IDF metric..


```r
kildisiute_signatures <- kildisiute_markers %>%
    # filter out mitochondrial and ribosomal genes
    dplyr::filter(., !grepl("RPS|RPL|MRPS|MRPL|^MT-", gene)) %>%
    group_by(Cluster) %>%
    arrange(desc(pct_diff)) %>%
    # make sure they're unique
    distinct(gene, .keep_all = TRUE) %>%
    dplyr::slice(1:100) %>%
    ungroup()

# re-plot the same figure with th efiltered markers
kildisiute_signatures %>%
    ggplot(aes(x = pct_diff, y = tfidf)) +
    geom_point(aes(colour = Cluster), alpha = 0.3) +
    geom_smooth(aes(colour = Cluster), method = "lm", se = FALSE) +
    xlab("Diff. in detection rate inside and outside cluster") + 
    ylab("TF-IDF metric") +
    xlim(c(0, 1)) + ylim(c(0, 5)) +
    ggtitle("After filtering")
```

```
## `geom_smooth()` using formula = 'y ~ x'
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/02//kildisiute_sigs-1.png)<!-- -->

```r
table(kildisiute_signatures$Cluster)
```

```
## 
##               Human fetal adrenal Bridge 
##                                      100 
##           Human fetal adrenal Chromaffin 
##                                      100 
##               Human fetal adrenal Cortex 
##                                      100 
##        Human fetal adrenal Erythroblasts 
##                                      100 
##           Human fetal adrenal Leukocytes 
##                                      100 
##           Human fetal adrenal Mesenchyme 
##                                      100 
##            Human fetal adrenal Podocytes 
##                                      100 
##                 Human fetal adrenal SCPs 
##                                      100 
##      Human fetal adrenal Sympathoblastic 
##                                      100 
## Human fetal adrenal Vascular Endothelium 
##                                      100
```

```r
kildisiute_signatures_ens <- kildisiute_signatures %>% signatures_to_list("ens")
kildisiute_signatures_sym <- kildisiute_signatures %>% signatures_to_list("sym")

kildisiute_cell_type_anno <- kildisiute_signatures %>%
    distinct(Cluster) %>%
    mutate(Dataset = "Kildisiute et al 2021",
           Age = "Human fetal 8-21PCW",
           Cell_type = Cluster,
           Sample = "Human fetal adrenal glands")
```


## Shi et al, Science, 2021

Here we'll select markers from Shi et al, Science, 2021, which correspond to differentially expressed genes between different human fetal
ganglionic eminence progenitor clusters.


```r
shi_table_s6 <- read_xlsx(here("data/singlecell/references_normal/Shi_Science_2021/data/Shi2021_table_s6.xlsx"), skip = 1)

# how many genes per cluster?
table(shi_table_s6$cluster)
```

```
## 
## pC1 pC2 pC3 pL1 pL2 pL3 pM1 pM2 pM3 pM4 
##  49  48  48  50  49  33  48  38  48  46
```

```r
length(unique(shi_table_s6$gene))
```

```
## [1] 305
```

```r
shi_signatures_ens <- shi_table_s6 %>%
    select(cluster, gene) %>%
    mutate(cluster = paste0("Human fetal ", cluster)) %>%
    split(.$cluster) %>%
    map(~ pull(.x, gene) %>% symbols2ensembl())

shi_signatures_sym <- shi_table_s6 %>%
    select(cluster, gene) %>%
    mutate(cluster = paste0("Human fetal ", cluster)) %>%
    split(.$cluster) %>%
    map(~ pull(.x, gene))

signatures_ens2 <- c(signatures_ens, shi_signatures_ens)

shi_cell_type_anno  <- shi_table_s6 %>%
    distinct(cluster) %>%
    select(Cluster = cluster) %>%
    mutate(Cluster = paste0("Human fetal ", Cluster)) %>%
    mutate(Dataset = "Shi et al 2021", Sample = "Human fetal GE", Age = "Human fetal",
           Cell_type = case_when(
               grepl("C", Cluster) ~ "CGE progenitor",
               grepl("M", Cluster) ~ "MGE progenitor",
               grepl("L", Cluster) ~ "LGE progenitor"
           ))
```

## Van Bruggen et al, Dev Cell, 2022

This study focused on emergence of OPCs in the early human fetal brain.


```r
vanbruggen_anno <- read_tsv(here("data/singlecell/references_normal/VanBruggen_DevCell_2022/cluster_annotation.tsv"))
```

```
## Rows: 31 Columns: 3
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (1): Annotation
## dbl (1): Cluster
## lgl (1): Color
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

We noted issues with gene naming in the supplementary tables provided with the manuscript,
so we re-called markers using `Seurat::FindAllMarkers()` and load them here:


```r
load(here("data/singlecell/references_normal/VanBruggen_DevCell_2022/vanbruggen_markers2.Rda"))

scRNA_human_dev <- readRDS(here("data/singlecell/references_normal/VanBruggen_DevCell_2022/scRNA_human_dev.rds"))
vanbruggen_anno <- data.frame("Cluster" = names(sort(table(scRNA_human_dev$Clusters), decreasing = TRUE))) %>% 
    tibble::rowid_to_column(var = "Number")

vanbruggen_signatures <- vanbruggen_markers2 %>% 
    dplyr::rename(avg_logFC = avg_log2FC, Cluster = cluster) %>% 
    filter_markers(n_top = 100) %>% 
    left_join(vanbruggen_anno, by = "Cluster") %>% 
    mutate(Cluster = paste0("Human fetal ", Number, "-", Cluster))

sort(table(vanbruggen_signatures$Cluster))
```

```
## 
##                             Human fetal 23-Striatum/Cortical neurons 
##                                                                    3 
##       Human fetal 24-Radial glia/Glioblast/Forebrain progenitor EMX1 
##                                                                    3 
##              Human fetal 30-Midbrain/Hindbrain inhibitory neuroblast 
##                                                                   48 
##                                   Human fetal 14-GABAergic forebrain 
##                                                                   53 
##                                  Human fetal 6-Cortical Interneurons 
##                                                                   55 
##        Human fetal 2-Inhibitory neurons Midbrain, possibly GABAergic 
##                                                                   63 
##                              Human fetal 1-Excitatory neurons cortex 
##                                                                  100 
##                                 Human fetal 10-Radial Glia/Glioblast 
##                                                                  100 
##                                    Human fetal 11-Cortical Pyramidal 
##                                                                  100 
##                       Human fetal 12-Forebrain inhibitory neuroblast 
##                                                                  100 
##                                   Human fetal 13-Radial Glia cycling 
##                                                                  100 
##                                           Human fetal 15-Radial Glia 
##                                                                  100 
##                              Human fetal 16-Mid/Hindbrain neuroblast 
##                                                                  100 
##                                     Human fetal 17-Glioblast/Pre-OPC 
##                                                                  100 
##                     Human fetal 18-Neuroblast motorneuron/GABAergic? 
##                                                                  100 
##                                                     Human fetal 19-U 
##                                                                  100 
##                              Human fetal 20-Radial Glia VLMC primed? 
##                                                                  100 
##                                  Human fetal 21-Hindbrain neuroblast 
##                                                                  100 
## Human fetal 22-GABAergic or interneuron neuroblast probably midbrain 
##                                                                  100 
##                                                  Human fetal 25-OPCs 
##                                                                  100 
##                                                 Human fetal 26-VLMCs 
##                                                                  100 
##                        Human fetal 27-Midbrain inhibitory neuroblast 
##                                                                  100 
##                                           Human fetal 28-Endothelial 
##                                                                  100 
##                                             Human fetal 29-Microglia 
##                                                                  100 
##                        Human fetal 3-Radial Glia potential glioblast 
##                                                                  100 
##                            Human fetal 4-Excitatory neurons midbrain 
##                                                                  100 
##          Human fetal 5-Forebrain early neuroblast possibly GABAergic 
##                                                                  100 
##                   Human fetal 7-Excitatory neurons possibly midbrain 
##                                                                  100 
##                       Human fetal 8-Forebrain neural progenitor EMX1 
##                                                                  100 
##                                              Human fetal 9-Glioblast 
##                                                                  100
```

```r
vanbruggen_signatures_ens <- vanbruggen_signatures %>%
    signatures_to_list("ens")

vanbruggen_signatures_sym <- vanbruggen_signatures %>%
    signatures_to_list("sym")

vanbruggen_cell_type_anno  <- vanbruggen_signatures %>%
    distinct(Cluster) %>%
    mutate(Dataset = "VanBruggen et al 2021",
           Sample = "Human fetal forebrain",
           Age = "8-10GW",
           Cell_type = Cluster)
```


## Yu et al, Nature Neuroscience, 2021

This study also investigated human fetal ganglionic eminences and cortical interneuron development.


```r
yu_cluster_markers <- read_xlsx(here("data/singlecell/references_normal/Yu_NatNeurosci_2021/Yu_NatNeuro_2021_Sup_tables.xlsx"),
                                sheet = 3, skip = 2) %>% 
    filter(!is.na(Cluster)) %>% 
    filter(!(Cluster %in% c("MG", "EC1", "EC2"))) %>% 
    mutate(Cluster = paste0("Fetal GE ", Cluster)) %>% 
    mutate(P_val = as.numeric(P_val), P_val_adj = as.numeric(P_val_adj))

yu_MGE_markers <- read_xlsx(here("data/singlecell/references_normal/Yu_NatNeurosci_2021/Yu_NatNeuro_2021_Sup_tables.xlsx"),
                            sheet = 11, skip = 2) %>% 
    mutate(Cluster = paste0("MGE-", Cluster))

yu_LGE_markers <- read_xlsx(here("data/singlecell/references_normal/Yu_NatNeurosci_2021/Yu_NatNeuro_2021_Sup_tables.xlsx"),
                            sheet = 12, skip = 2) %>% 
    mutate(Cluster = paste0("LGE-", Cluster)) %>% 
    mutate(P_val = as.numeric(P_val), P_val_adj = as.numeric(P_val_adj))

yu_CGE_markers <- read_xlsx(here("data/singlecell/references_normal/Yu_NatNeurosci_2021/Yu_NatNeuro_2021_Sup_tables.xlsx"),
                            sheet = 15, skip = 2) %>% 
    mutate(Cluster = paste0("CGE-", Cluster))

yu_markers <- bind_rows(yu_cluster_markers,
                        yu_MGE_markers,
                        yu_LGE_markers,
                        yu_CGE_markers) %>% 
    rename(gene = Gene, avg_logFC = Avg_logFC) %>% 
    mutate(Cluster = gsub("Fetal", "Human fetal", Cluster))
```

How many clusters, and how many markers per cluster?


```r
unique(yu_markers$Cluster)
```

```
##  [1] "Human fetal GE P1"    "Human fetal GE P2"    "Human fetal GE P3"   
##  [4] "Human fetal GE P4"    "Human fetal GE P5"    "Human fetal GE ExN4" 
##  [7] "Human fetal GE ExN5"  "Human fetal GE P6"    "Human fetal GE MGE1" 
## [10] "Human fetal GE MGE2"  "Human fetal GE CGE1"  "Human fetal GE LGE1" 
## [13] "Human fetal GE LGE2"  "Human fetal GE LGE3"  "Human fetal GE OPC"  
## [16] "Human fetal GE ExN1"  "Human fetal GE ExN2"  "Human fetal GE ExN3" 
## [19] "MGE-ZEB2+,MAF+"       "MGE-POU3F2+,CNTNAP2+" "MGE-NR2F1+,MEIS2+"   
## [22] "MGE-LHX8+,NKX2-1+"    "MGE-ANGPT2+,CRABP1+"  "LGE-ISL1+,EBF1+ 4"   
## [25] "LGE-ISL1+,EBF1+ 3"    "LGE-ISL1+,EBF1+ 2"    "LGE-ISL1+,EBF1+ 1"   
## [28] "LGE-SIX3+,SOX2+ 5"    "LGE-SIX3+,SOX2+ 4"    "LGE-SIX3+,SOX2+ 3"   
## [31] "LGE-SIX3+,SOX2+ 2"    "LGE-SIX3+,SOX2+ 1"    "LGE-PAX6+,ETV1+"     
## [34] "LGE-ASPM+,TOP2A+"     "CGE-CALB2+,NPAS3+"    "CGE-NFIA+,ST18+"     
## [37] "CGE-ANK3+,ENC1+"      "CGE-GRIA1+,SST+"      "CGE-TOP2A+,HMGB2+"
```

```r
table(yu_markers$Cluster)
```

```
## 
##      CGE-ANK3+,ENC1+    CGE-CALB2+,NPAS3+      CGE-GRIA1+,SST+ 
##                  116                   50                   86 
##      CGE-NFIA+,ST18+    CGE-TOP2A+,HMGB2+  Human fetal GE CGE1 
##                  178                   59                  104 
##  Human fetal GE ExN1  Human fetal GE ExN2  Human fetal GE ExN3 
##                  287                  192                  375 
##  Human fetal GE ExN4  Human fetal GE ExN5  Human fetal GE LGE1 
##                  261                  504                  115 
##  Human fetal GE LGE2  Human fetal GE LGE3  Human fetal GE MGE1 
##                  158                  130                  130 
##  Human fetal GE MGE2   Human fetal GE OPC    Human fetal GE P1 
##                  112                  444                  504 
##    Human fetal GE P2    Human fetal GE P3    Human fetal GE P4 
##                  482                  484                  444 
##    Human fetal GE P5    Human fetal GE P6     LGE-ASPM+,TOP2A+ 
##                  361                  249                  138 
##    LGE-ISL1+,EBF1+ 1    LGE-ISL1+,EBF1+ 2    LGE-ISL1+,EBF1+ 3 
##                   67                  107                  179 
##    LGE-ISL1+,EBF1+ 4      LGE-PAX6+,ETV1+    LGE-SIX3+,SOX2+ 1 
##                  113                  194                   90 
##    LGE-SIX3+,SOX2+ 2    LGE-SIX3+,SOX2+ 3    LGE-SIX3+,SOX2+ 4 
##                  115                  246                  115 
##    LGE-SIX3+,SOX2+ 5  MGE-ANGPT2+,CRABP1+    MGE-LHX8+,NKX2-1+ 
##                  140                  104                  225 
##    MGE-NR2F1+,MEIS2+ MGE-POU3F2+,CNTNAP2+       MGE-ZEB2+,MAF+ 
##                  109                   63                   98
```

The imbalance in number of markers per cluster means that we should do some filtering to get them to roughly equal length:


```r
yu_signatures <- yu_markers %>% 
    filter_markers(n_top = 100)

sort(table(yu_signatures$Cluster))
```

```
## 
##    CGE-CALB2+,NPAS3+    CGE-TOP2A+,HMGB2+ MGE-POU3F2+,CNTNAP2+ 
##                   50                   59                   63 
##    LGE-ISL1+,EBF1+ 1      CGE-GRIA1+,SST+    LGE-SIX3+,SOX2+ 1 
##                   67                   86                   90 
##       MGE-ZEB2+,MAF+      CGE-ANK3+,ENC1+      CGE-NFIA+,ST18+ 
##                   98                  100                  100 
##  Human fetal GE CGE1  Human fetal GE ExN1  Human fetal GE ExN2 
##                  100                  100                  100 
##  Human fetal GE ExN3  Human fetal GE ExN4  Human fetal GE ExN5 
##                  100                  100                  100 
##  Human fetal GE LGE1  Human fetal GE LGE2  Human fetal GE LGE3 
##                  100                  100                  100 
##  Human fetal GE MGE1  Human fetal GE MGE2   Human fetal GE OPC 
##                  100                  100                  100 
##    Human fetal GE P1    Human fetal GE P2    Human fetal GE P3 
##                  100                  100                  100 
##    Human fetal GE P4    Human fetal GE P5    Human fetal GE P6 
##                  100                  100                  100 
##     LGE-ASPM+,TOP2A+    LGE-ISL1+,EBF1+ 2    LGE-ISL1+,EBF1+ 3 
##                  100                  100                  100 
##    LGE-ISL1+,EBF1+ 4      LGE-PAX6+,ETV1+    LGE-SIX3+,SOX2+ 2 
##                  100                  100                  100 
##    LGE-SIX3+,SOX2+ 3    LGE-SIX3+,SOX2+ 4    LGE-SIX3+,SOX2+ 5 
##                  100                  100                  100 
##  MGE-ANGPT2+,CRABP1+    MGE-LHX8+,NKX2-1+    MGE-NR2F1+,MEIS2+ 
##                  100                  100                  100
```

```r
yu_signatures_ens <- yu_signatures %>%
    signatures_to_list("ens")

yu_signatures_sym <- yu_signatures %>%
    signatures_to_list("sym")

yu_cell_type_anno  <- yu_signatures %>%
    distinct(Cluster) %>%
    mutate(Dataset = "Yu et al 2021",
           Sample = "Human fetal GE",
           Age = "Human fetal",
           Cell_type = Cluster)
```


## Compile all signatures


```r
cell_type_anno_all <- bind_rows(cell_type_anno_filt,
                                jessa_cell_type_anno,
                                shi_cell_type_anno,
                                yu_cell_type_anno,
                                vanbruggen_cell_type_anno,
                                dong_cell_type_anno,
                                kildisiute_cell_type_anno,
                                jansky_cell_type_anno,
                                aldinger_cell_type_anno) %>%
    mutate(Species = case_when(
        grepl("Velmeshev|Nowakowski|Yu|Shi|Dong|Kildisiute|Aldinger|VanBruggen|Jansky", Dataset) ~ "Human",
        TRUE ~ "Mouse"
    ))

cell_type_anno_all %>% distinct(Dataset, Species) %>% arrange(Species)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Dataset"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Species"],"name":[2],"type":["chr"],"align":["left"]}],"data":[{"1":"Nowakowski et al 2017","2":"Human"},{"1":"Velmeshev et al 2019","2":"Human"},{"1":"Shi et al 2021","2":"Human"},{"1":"Yu et al 2021","2":"Human"},{"1":"VanBruggen et al 2021","2":"Human"},{"1":"Dong et al 2020","2":"Human"},{"1":"Kildisiute et al 2021","2":"Human"},{"1":"Jansky et al 2021","2":"Human"},{"1":"Aldinger et al 2021","2":"Human"},{"1":"Mizrak et al 2019","2":"Mouse"},{"1":"Anderson et al 2020","2":"Mouse"},{"1":"Jessa et al 2022","2":"Mouse"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
signatures_ens_all <- c(signatures_ens, jessa_signatures_ens, shi_signatures_ens, yu_signatures_ens,
                        dong_signatures_ens, kildisiute_signatures_ens,
                        vanbruggen_signatures_ens, jansky_signatures_ens,
                        aldinger_signatures_ens)
signatures_sym_all <- c(signatures_sym, jessa_signatures_sym, shi_signatures_sym, yu_signatures_sym,
                        dong_signatures_sym, kildisiute_signatures_sym,
                        vanbruggen_signatures_sym, jansky_signatures_sym,
                        aldinger_signatures_sym)

# QC: check length
sig_length <- map_dbl(signatures_ens_all, length) %>%
    tibble::enframe("Signature", "Length") %>%
    arrange(Length) %>%
    mutate(Signature = factor(Signature, levels = Signature))

sig_length %>% ggplot(aes(x = Signature, y = Length)) +
    geom_col() +
    rotate_x() +
    geom_hline(yintercept = 75, colour = "red") +
    theme(axis.text.x = element_text(size = rel(0.5)))
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/02//compile_sigs-1.png)<!-- -->

```r
signatures_keep <- sig_length %>%
    filter(Length >= 75) %>%
    pull(Signature) %>%
    as.character()

# what is being dropped
sig_length %>%
    filter(Length < 75) %>%
    pull(Signature) %>% 
    as.character()
```

```
##  [1] "Human fetal 23-Striatum/Cortical neurons"                      
##  [2] "Human fetal 24-Radial glia/Glioblast/Forebrain progenitor EMX1"
##  [3] "Human fetal cerebellum 20-Choroid"                             
##  [4] "Human fetal pL3"                                               
##  [5] "Human fetal pM2"                                               
##  [6] "Human fetal 30-Midbrain/Hindbrain inhibitory neuroblast"       
##  [7] "Human fetal 14-GABAergic forebrain"                            
##  [8] "Human fetal pC1"                                               
##  [9] "Human fetal pL1"                                               
## [10] "Human fetal pM1"                                               
## [11] "CGE-CALB2+,NPAS3+"                                             
## [12] "Human fetal 2-Inhibitory neurons Midbrain, possibly GABAergic" 
## [13] "Human fetal pL2"                                               
## [14] "Human fetal 6-Cortical Interneurons"                           
## [15] "Human fetal pC2"                                               
## [16] "Human fetal pM3"                                               
## [17] "Human fetal pM4"                                               
## [18] "Human fetal pC3"                                               
## [19] "Human fetal adrenal2 connecting Progenitor cells"              
## [20] "CGE-TOP2A+,HMGB2+"                                             
## [21] "MGE-POU3F2+,CNTNAP2+"                                          
## [22] "LGE-ISL1+,EBF1+ 1"
```

```r
# total number of signatures
dim(cell_type_anno_all)
```

```
## [1] 412   6
```

Since every study uses a slightly different cell type naming scheme, we will assign a harmonized class label for each cluster, for ease of interpretation:


```r
cell_type_anno_all <- cell_type_anno_all %>% 
    mutate(Keep = ifelse(Cluster %in% signatures_keep, TRUE, FALSE)) %>% 
    mutate(Class = case_when(
        grepl("-P|prolif|TOP2A|Cycling|div", Cluster) & !grepl("Pericytes", Cluster) ~ "Proliferating progenitors",
        grepl("RG|[Rr]adial", Cluster) &
            !grepl("NRGN", Cluster) &
            !grepl("[Oo]lig", Cluster) ~ "Radial glia",
        grepl("NSC", Cluster) ~ "Neural stem cell",
        grepl("OPC|Oligodendrocyte progenitor cell|Oligodendrocyte precurso|COP", Cluster) ~ "Oligodendrocyte precursors",
        grepl("NFOL|MOL|[Oo]ligo|OL", Cluster) ~ "Oligodendrocytes",
        grepl("ASEP|EPEN|[Ee]pendymal|[Cc]horoid|CPLX", Cluster) ~ "Ependymal",
        grepl("AST|Ast|Astr|BG", Cluster) ~ "Astrocytes",
        grepl("MNG|[Mm]eninges|VLM", Cluster) ~ "Meninges",
        grepl("Mesenchym", Cluster) ~ "Mesenchyme",
        grepl("PERI|[Pp]ericyte|[Mm]ural|Vascular|ENDO|[Ee]ndo|[Ff]ibroblast|VSMC", Cluster) ~ "Mural",
        grepl("Muscle|Myo", Cluster) ~ "Muscle",
        grepl("Glia|GLIP|[Gg]lioblast", Cluster) ~ "Glial progenitors",
        grepl("ERY|Immune|MGL|MAC|Microgl|Macro|Leuko|Erythro| T ", Cluster) ~ "Immune",
        grepl("SCP|SCHW|[Cc]hromaffin|Bridge|[Ss]ympatho|Adrenal|NCC|Melano", Cluster) ~ "Neural crest lineages",
        grepl("EN|EMX|EXN|ExN|[Ee]xcit|nEN|Pyram", Cluster) & !grepl("ENC", Cluster) ~ "Excitatory neurons",
        grepl("Human fetal GE [pP]", Cluster) ~ "GE neural precursors",
        grepl("GE|GABA|nIN|MGIN|MGE|CGE|LGE|SST|PV|[Ii]nhib|CIN|PV|SST|VIP|SV2C|Somato|Parv|pC|pM|pL|Inter", Cluster) &
            !grepl("[Mm]idbrain|[Hh]indbrain", Cluster) ~ "GE/Inhibitory neurons",
        grepl("EXIP|NEURP|IP|[Ii]ntermediate|prog|Neurogenic", Cluster) &
            !grepl("VIP", Cluster) & !grepl("Human fetal progenitor:", Cluster) ~ "Neuronal progenitors",
        grepl("SPN|Str|STR|SMSN|Striatum", Cluster) ~ "Striatal neurons",
        grepl("GN|iCN|INs|MLI|Moto|MN|NRGN|[Nn]eu|MFN|CJRZN|UBC|GABAN|NEUR|THLN|SN|INH|SRN|ACHN", Cluster) & !grepl("SMSN", Cluster) ~ "Other neurons",
        grepl("Human fetal cortex:|Human fetal progenitor:|Human fetal CS12-13:", Cluster) ~ "Unlabelled",
        TRUE ~ "Other"))

save(cell_type_anno_all, signatures_keep, signatures_sym_all, signatures_ens_all, signatures_keep,
     file = glue("{out}/signatures.Rda"))
```

# Prepare supplementary table

Here, we'll prepare a supplementary table containing all the reference gene signatures
and their associated cell types/datasets from which they originate.



```r
# collapse signatures into strings
sigs_sym_df <- imap_dfr(signatures_sym_all,
                        ~ data.frame(Cluster = .y, Signature_sym = glue_collapse(.x, sep = ",")))

sigs_ens_df <- imap_dfr(signatures_ens_all,
                        ~ data.frame(Cluster = .y, Signature_ENS = glue_collapse(.x, sep = ",")))

# join with cell type annotation
TABLE_ref_sigs <- cell_type_anno_all %>% 
    left_join(sigs_sym_df) %>% 
    left_join(sigs_ens_df) %>%
    # filter to kept signatures
    filter(Keep) %>% 
    select(-Keep) %>%
    mutate(Cell_type = ifelse(is.na(Cell_type), Cluster, Cell_type))
```

```
## Joining with `by = join_by(Cluster)`
## Joining with `by = join_by(Cluster)`
```

```r
rr_write_tsv(TABLE_ref_sigs,
             glue("{out}/TABLE_reference_signatures.tsv"),
             "Cell type specific reference signatures and associated metadata and dataset info")
```

```
## ...writing description of TABLE_reference_signatures.tsv to public/output/02/TABLE_reference_signatures.desc
```


# Compare signatures

We have assembled signatures from two species, from multiple marker-identification approaches, and from multiple studies. As a sanity check, we compute pairwise Jaccard Index for overlap in genes between the signatures, in order to confirm that similar cell types do have similar gene signatures.


```r
# calculate Jaccard index as follows:
# |X and Y| / |X or Y|
jaccard <- function(x, y) length(base::intersect(x, y)) / length(base::union(x, y))

signatures_overlap <- sapply(signatures_sym_all[signatures_keep],
                             function(x) sapply(signatures_sym_all[signatures_keep],
                                                function(y) jaccard(x, y)))

anno <- cell_type_anno_all %>% select(Cluster, Dataset, Species, Class) %>% tibble::column_to_rownames(var = "Cluster")

palette_dataset <- set_qual_pal(length(unique(anno$Dataset)))
names(palette_dataset) <- unique(anno$Dataset)

hm_jaccard <- pheatmap(signatures_overlap,
                       color             = custom_magma,
                       border_color      = NA,
                       treeheight_row    = 50,
                       treeheight_col    = 0,
                       cutree_rows       = 35,
                       cutree_cols       = 35,
                       annotation_row    = anno,
                       annotation_col    = anno,
                       annotation_colors = list("Dataset" = palette_dataset,
                                                "Species" = c("Mouse" = "gray90", "Human" = "black"),
                                                "Class"   = palette_class),
                       cellwidth         = 3,
                       cellheight        = 3,
                       fontsize_row      = 4,
                       fontsize_col      = 4,
                       filename          = glue("{figout}/signature_jaccard.png"))

# save heatmap to extract clustering order
save(hm_jaccard, file = glue("{out}/hm_jaccard.Rda"))

knitr::include_graphics(glue("{figout}/signature_jaccard.png"))
```

<img src="/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/02///signature_jaccard.png" width="7333" />

- We confirm that there is coherence across datasets, with major clusters for microglia, endothelial other vascular/non-neuroetodermal cells, OPCs, astro-related cells, cycling cells, different types of neurons
- This suggests that the signatures are informative and capturing cell type-specific rather that simply dataset-specific information


# Stats


```r
cell_type_stats <- cell_type_anno_all %>%
    filter(Cluster %in% signatures_keep) %>% 
    mutate(Sample = case_when(
        Dataset == "Jessa et al 2019" & grepl("Forebrain", Sample) ~ "Mouse forebrain",
        Dataset == "Jessa et al 2019" & grepl("Pons", Sample) ~ "Mouse pons",
        TRUE ~ Sample
    )) %>% 
    mutate(Age = case_when(
        Dataset == "Jessa et al 2019" & grepl("forebrain", Sample) ~ "Mouse E12-P6",
        Dataset == "Jessa et al 2019" & grepl("pons", Sample) ~ "Mouse E12-P6",
        TRUE ~ Age
    )) %>% 
    group_by(Sample, Age, Dataset) %>% count() %>% arrange(Dataset)

sum(cell_type_stats$n)
```

```
## [1] 390
```

```r
kableExtra::kable(cell_type_stats)
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Sample </th>
   <th style="text-align:left;"> Age </th>
   <th style="text-align:left;"> Dataset </th>
   <th style="text-align:right;"> n </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Human fetal cerebellum </td>
   <td style="text-align:left;"> 9-21 PCW </td>
   <td style="text-align:left;"> Aldinger et al 2021 </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Striatum </td>
   <td style="text-align:left;"> Mouse P9 </td>
   <td style="text-align:left;"> Anderson et al 2020 </td>
   <td style="text-align:right;"> 38 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Human 4PCW embryo neural lineages </td>
   <td style="text-align:left;"> Human fetal </td>
   <td style="text-align:left;"> Dong et al 2020 </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Human 8-14PCW fetal adrenal glands </td>
   <td style="text-align:left;"> Human fetal </td>
   <td style="text-align:left;"> Dong et al 2020 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Human fetal adrenal gland </td>
   <td style="text-align:left;"> Human fetal 7-17PCW </td>
   <td style="text-align:left;"> Jansky et al 2021 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mouse forebrain </td>
   <td style="text-align:left;"> E10.5 </td>
   <td style="text-align:left;"> Jessa et al 2022 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mouse forebrain </td>
   <td style="text-align:left;"> E12.5 </td>
   <td style="text-align:left;"> Jessa et al 2022 </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mouse forebrain </td>
   <td style="text-align:left;"> E13.5 </td>
   <td style="text-align:left;"> Jessa et al 2022 </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mouse forebrain </td>
   <td style="text-align:left;"> E15.5 </td>
   <td style="text-align:left;"> Jessa et al 2022 </td>
   <td style="text-align:right;"> 17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mouse forebrain </td>
   <td style="text-align:left;"> E16.5 </td>
   <td style="text-align:left;"> Jessa et al 2022 </td>
   <td style="text-align:right;"> 19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mouse forebrain </td>
   <td style="text-align:left;"> E18.5 </td>
   <td style="text-align:left;"> Jessa et al 2022 </td>
   <td style="text-align:right;"> 23 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mouse forebrain </td>
   <td style="text-align:left;"> P0 </td>
   <td style="text-align:left;"> Jessa et al 2022 </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mouse forebrain </td>
   <td style="text-align:left;"> P3 </td>
   <td style="text-align:left;"> Jessa et al 2022 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mouse forebrain </td>
   <td style="text-align:left;"> P6 </td>
   <td style="text-align:left;"> Jessa et al 2022 </td>
   <td style="text-align:right;"> 13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Human fetal adrenal glands </td>
   <td style="text-align:left;"> Human fetal 8-21PCW </td>
   <td style="text-align:left;"> Kildisiute et al 2021 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> V-SVZ </td>
   <td style="text-align:left;"> Mouse adult </td>
   <td style="text-align:left;"> Mizrak et al 2019 </td>
   <td style="text-align:right;"> 22 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Human fetal cortex/MGE </td>
   <td style="text-align:left;"> Human fetal </td>
   <td style="text-align:left;"> Nowakowski et al 2017 </td>
   <td style="text-align:right;"> 47 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Human fetal forebrain </td>
   <td style="text-align:left;"> 8-10GW </td>
   <td style="text-align:left;"> VanBruggen et al 2021 </td>
   <td style="text-align:right;"> 24 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Human ped/adult cortex </td>
   <td style="text-align:left;"> Human ped/adult </td>
   <td style="text-align:left;"> Velmeshev et al 2019 </td>
   <td style="text-align:right;"> 13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Human fetal GE </td>
   <td style="text-align:left;"> Human fetal </td>
   <td style="text-align:left;"> Yu et al 2021 </td>
   <td style="text-align:right;"> 35 </td>
  </tr>
</tbody>
</table>


```r
cell_type_stats %>%
    tibble::rowid_to_column(var = "Order") %>% 
    ggplot(aes(x = Order, y = n)) +
    geom_col(aes(fill = Dataset), width = 0.5) +
    scale_fill_manual(values = palette_dataset) +
    coord_flip()
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/02//signature_stats_barplot-1.png)<!-- -->



<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2024-11-01 11:43:26
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
##  ! package        * version date (UTC) lib source
##  P abind            1.4-5   2016-07-21 [?] CRAN (R 4.1.2)
##  P BiocManager      1.30.15 2021-05-11 [?] CRAN (R 4.1.2)
##  P bit              4.0.4   2020-08-04 [?] CRAN (R 4.1.2)
##  P bit64            4.0.5   2020-08-30 [?] CRAN (R 4.1.2)
##  P bslib            0.3.1   2021-10-06 [?] CRAN (R 4.1.2)
##  P cachem           1.0.6   2021-08-19 [?] CRAN (R 4.1.2)
##  P callr            3.7.6   2024-03-25 [?] RSPM
##  P cellranger       1.1.0   2016-07-27 [?] CRAN (R 4.1.2)
##  P cli              3.6.1   2023-03-23 [?] RSPM (R 4.1.2)
##  P cluster          2.1.2   2021-04-17 [?] CRAN (R 4.1.2)
##  P codetools        0.2-18  2020-11-04 [?] CRAN (R 4.1.2)
##  P colorspace       2.0-2   2021-06-24 [?] CRAN (R 4.1.2)
##  P conflicted       1.2.0   2023-02-01 [?] CRAN (R 4.1.2)
##  P cowplot          1.1.1   2020-12-30 [?] CRAN (R 4.1.2)
##  P crayon           1.4.2   2021-10-29 [?] CRAN (R 4.1.2)
##  P data.table       1.14.2  2021-09-27 [?] CRAN (R 4.1.2)
##  P deldir           1.0-6   2021-10-23 [?] CRAN (R 4.1.2)
##  P devtools         2.4.5   2022-10-11 [?] CRAN (R 4.1.2)
##  P digest           0.6.35  2024-03-11 [?] CRAN (R 4.1.2)
##  P dplyr          * 1.1.1   2023-03-22 [?] CRAN (R 4.1.2)
##  P ellipsis         0.3.2   2021-04-29 [?] CRAN (R 4.1.2)
##  P evaluate         0.23    2023-11-01 [?] CRAN (R 4.1.2)
##  P fansi            1.0.2   2022-01-14 [?] CRAN (R 4.1.2)
##  P farver           2.1.0   2021-02-28 [?] CRAN (R 4.1.2)
##  P fastmap          1.1.0   2021-01-25 [?] CRAN (R 4.1.2)
##  P fitdistrplus     1.1-6   2021-09-28 [?] CRAN (R 4.1.2)
##  P fs               1.5.2   2021-12-08 [?] CRAN (R 4.1.2)
##  P future           1.25.0  2022-04-24 [?] CRAN (R 4.1.2)
##  P future.apply     1.8.1   2021-08-10 [?] CRAN (R 4.1.2)
##  P generics         0.1.3   2022-07-05 [?] CRAN (R 4.1.2)
##  P ggplot2        * 3.4.2   2023-04-03 [?] CRAN (R 4.1.2)
##  P ggrepel          0.9.1   2021-01-15 [?] CRAN (R 4.1.2)
##  P ggridges         0.5.3   2021-01-08 [?] CRAN (R 4.1.2)
##  P git2r            0.29.0  2021-11-22 [?] CRAN (R 4.1.2)
##  P globals          0.14.0  2020-11-22 [?] CRAN (R 4.1.2)
##  P glue           * 1.6.2   2022-02-24 [?] CRAN (R 4.1.2)
##  P goftest          1.2-3   2021-10-07 [?] CRAN (R 4.1.2)
##  P gridExtra        2.3     2017-09-09 [?] CRAN (R 4.1.2)
##  P gtable           0.3.0   2019-03-25 [?] CRAN (R 4.1.2)
##  P here           * 1.0.1   2020-12-13 [?] CRAN (R 4.1.2)
##  P highr            0.9     2021-04-16 [?] CRAN (R 4.1.2)
##  P hms              1.1.1   2021-09-26 [?] CRAN (R 4.1.2)
##  P htmltools        0.5.2   2021-08-25 [?] CRAN (R 4.1.2)
##  P htmlwidgets      1.5.4   2021-09-08 [?] CRAN (R 4.1.2)
##  P httpuv           1.6.5   2022-01-05 [?] CRAN (R 4.1.2)
##  P httr             1.4.2   2020-07-20 [?] CRAN (R 4.1.2)
##  P ica              1.0-2   2018-05-24 [?] CRAN (R 4.1.2)
##  P igraph           2.0.3   2024-03-13 [?] CRAN (R 4.1.2)
##  P irlba            2.3.5   2021-12-06 [?] CRAN (R 4.1.2)
##  P jquerylib        0.1.4   2021-04-26 [?] CRAN (R 4.1.2)
##  P jsonlite         1.8.8   2023-12-04 [?] CRAN (R 4.1.2)
##  P kableExtra       1.3.4   2021-02-20 [?] CRAN (R 4.1.2)
##  P KernSmooth       2.23-20 2021-05-03 [?] CRAN (R 4.1.2)
##  P knitr            1.37    2021-12-16 [?] CRAN (R 4.1.2)
##  P labeling         0.4.2   2020-10-20 [?] CRAN (R 4.1.2)
##  P later            1.3.0   2021-08-18 [?] CRAN (R 4.1.2)
##  P lattice          0.20-45 2021-09-22 [?] CRAN (R 4.1.2)
##  P lazyeval         0.2.2   2019-03-15 [?] CRAN (R 4.1.2)
##  P leiden           0.3.9   2021-07-27 [?] CRAN (R 4.1.2)
##  P lifecycle        1.0.3   2022-10-07 [?] CRAN (R 4.1.2)
##  P listenv          0.8.0   2019-12-05 [?] CRAN (R 4.1.2)
##  P lmtest           0.9-39  2021-11-07 [?] CRAN (R 4.1.2)
##  P magrittr       * 2.0.3   2022-03-30 [?] CRAN (R 4.1.2)
##  P MASS             7.3-54  2021-05-03 [?] CRAN (R 4.1.2)
##  P Matrix           1.3-4   2021-06-01 [?] CRAN (R 4.1.2)
##  P matrixStats      0.61.0  2021-09-17 [?] CRAN (R 4.1.2)
##  P memoise          2.0.1   2021-11-26 [?] CRAN (R 4.1.2)
##  P mgcv             1.8-38  2021-10-06 [?] CRAN (R 4.1.2)
##  P mime             0.12    2021-09-28 [?] CRAN (R 4.1.2)
##  P miniUI           0.1.1.1 2018-05-18 [?] CRAN (R 4.1.2)
##  P munsell          0.5.0   2018-06-12 [?] CRAN (R 4.1.2)
##  P nlme             3.1-153 2021-09-07 [?] CRAN (R 4.1.2)
##  P parallelly       1.30.0  2021-12-17 [?] CRAN (R 4.1.2)
##  P patchwork        1.1.1   2020-12-17 [?] CRAN (R 4.1.2)
##  P pbapply          1.5-0   2021-09-16 [?] CRAN (R 4.1.2)
##  P pheatmap       * 1.0.12  2019-01-04 [?] CRAN (R 4.1.2)
##  P pillar           1.9.0   2023-03-22 [?] RSPM (R 4.1.2)
##  P pkgbuild         1.4.2   2023-06-26 [?] CRAN (R 4.1.2)
##  P pkgconfig        2.0.3   2019-09-22 [?] CRAN (R 4.1.2)
##  P pkgload          1.3.3   2023-09-22 [?] CRAN (R 4.1.2)
##  P plotly           4.10.0  2021-10-09 [?] CRAN (R 4.1.2)
##  P plyr             1.8.6   2020-03-03 [?] CRAN (R 4.1.2)
##  P png              0.1-7   2013-12-03 [?] CRAN (R 4.1.2)
##  P polyclip         1.10-0  2019-03-14 [?] CRAN (R 4.1.2)
##  P prettyunits      1.1.1   2020-01-24 [?] CRAN (R 4.1.2)
##  P processx         3.8.4   2024-03-16 [?] RSPM
##  P profvis          0.3.8   2023-05-02 [?] CRAN (R 4.1.2)
##  P promises         1.2.0.1 2021-02-11 [?] CRAN (R 4.1.2)
##  P ps               1.7.6   2024-01-18 [?] RSPM
##  P purrr          * 1.0.1   2023-01-10 [?] CRAN (R 4.1.2)
##  P R6               2.5.1   2021-08-19 [?] CRAN (R 4.1.2)
##  P RANN             2.6.1   2019-01-08 [?] CRAN (R 4.1.2)
##  P RColorBrewer   * 1.1-2   2014-12-07 [?] CRAN (R 4.1.2)
##  P Rcpp             1.0.8   2022-01-13 [?] CRAN (R 4.1.2)
##  P RcppAnnoy        0.0.19  2021-07-30 [?] CRAN (R 4.1.2)
##  P readr          * 2.1.1   2021-11-30 [?] CRAN (R 4.1.2)
##  P readxl         * 1.3.1   2019-03-13 [?] CRAN (R 4.1.2)
##  P remotes          2.4.2.1 2023-07-18 [?] CRAN (R 4.1.2)
##  P renv             1.0.3   2023-09-19 [?] CRAN (R 4.1.2)
##  P reshape2         1.4.4   2020-04-09 [?] CRAN (R 4.1.2)
##  P reticulate       1.23    2022-01-14 [?] CRAN (R 4.1.2)
##  P rlang            1.1.3   2024-01-10 [?] CRAN (R 4.1.2)
##  P rmarkdown        2.11    2021-09-14 [?] CRAN (R 4.1.2)
##  P ROCR             1.0-11  2020-05-02 [?] CRAN (R 4.1.2)
##  P rpart            4.1-15  2019-04-12 [?] CRAN (R 4.1.2)
##  P rprojroot        2.0.2   2020-11-15 [?] CRAN (R 4.1.2)
##  P rstudioapi       0.13    2020-11-12 [?] CRAN (R 4.1.2)
##  P Rtsne            0.15    2018-11-10 [?] CRAN (R 4.1.2)
##  P rvest            1.0.2   2021-10-16 [?] CRAN (R 4.1.2)
##  P sass             0.4.0   2021-05-12 [?] CRAN (R 4.1.2)
##  P scales           1.2.1   2022-08-20 [?] CRAN (R 4.1.2)
##  P scattermore      0.7     2020-11-24 [?] CRAN (R 4.1.2)
##  P sctransform      0.3.3   2022-01-13 [?] CRAN (R 4.1.2)
##  P sessioninfo      1.2.2   2021-12-06 [?] CRAN (R 4.1.2)
##  P Seurat         * 4.0.0   2021-01-30 [?] CRAN (R 4.1.2)
##  P SeuratObject   * 4.0.4   2021-11-23 [?] CRAN (R 4.1.2)
##  P shiny            1.7.1   2021-10-02 [?] CRAN (R 4.1.2)
##  P spatstat         1.64-1  2020-05-12 [?] CRAN (R 4.1.2)
##  P spatstat.data    2.1-2   2021-12-17 [?] CRAN (R 4.1.2)
##  P spatstat.utils   2.3-0   2021-12-12 [?] CRAN (R 4.1.2)
##  P stringi          1.7.6   2021-11-29 [?] CRAN (R 4.1.2)
##  P stringr          1.5.0   2022-12-02 [?] CRAN (R 4.1.2)
##  P survival         3.2-13  2021-08-24 [?] CRAN (R 4.1.2)
##  P svglite          2.0.0   2021-02-20 [?] CRAN (R 4.1.2)
##  P systemfonts      1.0.3   2021-10-13 [?] CRAN (R 4.1.2)
##  P tensor           1.5     2012-05-05 [?] CRAN (R 4.1.2)
##  P tibble           3.2.1   2023-03-20 [?] RSPM (R 4.1.2)
##  P tidyr          * 1.3.0   2023-01-24 [?] CRAN (R 4.1.2)
##  P tidyselect       1.2.0   2022-10-10 [?] CRAN (R 4.1.2)
##  P tzdb             0.3.0   2022-03-28 [?] CRAN (R 4.1.2)
##  P urlchecker       1.0.1   2021-11-30 [?] CRAN (R 4.1.2)
##  P usethis          2.2.2   2023-07-06 [?] CRAN (R 4.1.2)
##  P utf8             1.2.2   2021-07-24 [?] CRAN (R 4.1.2)
##  P uwot             0.1.11  2021-12-02 [?] CRAN (R 4.1.2)
##  P vctrs            0.6.5   2023-12-01 [?] CRAN (R 4.1.2)
##  P viridis        * 0.5.1   2018-03-29 [?] RSPM (R 4.1.2)
##  P viridisLite    * 0.3.0   2018-02-01 [?] CRAN (R 4.1.2)
##  P vroom            1.5.7   2021-11-30 [?] CRAN (R 4.1.2)
##  P webshot          0.5.2   2019-11-22 [?] CRAN (R 4.1.2)
##  P withr            2.5.0   2022-03-03 [?] CRAN (R 4.1.2)
##  P xfun             0.29    2021-12-14 [?] CRAN (R 4.1.2)
##  P xml2             1.3.3   2021-11-30 [?] CRAN (R 4.1.2)
##  P xtable           1.8-4   2019-04-21 [?] CRAN (R 4.1.2)
##  P yaml             2.2.1   2020-02-01 [?] CRAN (R 4.1.2)
##  P zoo              1.8-9   2021-03-09 [?] CRAN (R 4.1.2)
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
