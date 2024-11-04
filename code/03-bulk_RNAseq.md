---
title: "03 - Prep tumor bulk RNAseq data"
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
## Document index: 03
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
## public/output/03
```

```
## public/figures/03//
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

This document performs a few preparatory steps in the analysis of bulk RNAseq data,
including checking QC, expression levels of FOXR2 and other genes, and distribution of samples in the PCA space.

For an external dataset of extra-cranial neuroblastomas, we load the raw counts, normalize the counts using DESeq2, and plot the samples in PCA space.

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
library(DESeq2)
library(scales)
library(cowplot)

source(here("include/style.R"))
source(here("code/functions/RNAseq.R"))

ggplot2::theme_set(theme_min())
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

</div>

# In-house pediatric brain tumors


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
table(info_samples$Group)
```

```
## 
##       DIPG-H3K27M DIPG-H3K27M-FOXR2            EP-PFA              ETMR 
##                19                 6                13                12 
##    HGG-H3.3G34R/V           HGG-IDH            MB-SHH            MB-WNT 
##                18                10                 8                10 
##          NB-FOXR2 
##                25
```

```r
# load counts for FOXR2
counts_foxr2 <- extract_pipeline_counts(path = here("data/RNAseq/pipeline_l3/2023-05-test_pbt/counts/Ensembl.ensGene.exon.norm.tsv.gz"),
                                        goi = "FOXR2") %>% 
    left_join(info_samples, by = c("sample" = "Nickname")) %>% 
    mutate(Group = factor(Group, levels = names(palette_groups)))
```

```
## Joining with `by = join_by(gene_symbol)`
```

## Counts


```r
counts_RNA <- read.table(here("data/RNAseq/pipeline_l3/2023-05-test_pbt/counts/Ensembl.ensGene.exon.raw.tsv.gz"),
                         header = T, sep = "\t", check.names = FALSE) %>%
    tibble::rownames_to_column(var = "ID") %>% 
    separate(ID, into = c("ENSID", "symbol"), sep = ":") %>% 
    arrange(ENSID) %>%
    .[, c("ENSID", "symbol", info_samples$Nickname)] %T>%
    rr_write_tsv(glue("{out}/TABLE_bulk_counts.tsv"),
                 "Raw counts for bulk RNAseq data")
```

```
## ...writing description of TABLE_bulk_counts.tsv to public/output/03/TABLE_bulk_counts.desc
```

## QC

Alignment stats:


```r
align_stats <- read_tsv(here("data/RNAseq/pipeline_l3/2023-05-test_pbt/alignment.statistics.tsv"))
```

```
## Rows: 121 Columns: 24
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (12): name, group, cleanReadRetentionPercentage, mappedPercentage, whole...
## dbl (12): rawReadCount, cleanReadCount, unmappedReadCount, mappedReadCount, ...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
align_stats <- align_stats %>% 
    left_join(info_samples, by = c("name" = "Nickname")) %>% 
    separate(mitochondrialPercentage, sep = "%", into = "mitochondrialPercentage") %>%
    separate(ribosomalPercentage,     sep = "%", into = "ribosomalPercentage") %>%
    mutate(mitochondrialPercentage = as.numeric(mitochondrialPercentage),
           ribosomalPercentage     = as.numeric(ribosomalPercentage))

align_stats %>% 
    select(name, group, cleanReadCount, mappedReadCount,
           mitochondrialPercentage, ribosomalPercentage) %>% 
    gather(stat, value, 3:6) %>%
    mutate(group = factor(group, levels = names(palette_groups))) %>% 
    arrange(desc(group)) %>% 
    mutate(name = factor(name, levels = unique(.$name))) %>% 
    ggplot(aes(x = name, y = value)) +
    geom_bar(aes(fill = group), stat = "identity") +
    scale_fill_manual(values = palette_groups) +
    facet_wrap(~ stat, scales = "free_y", ncol = 1) +
    rotate_x()
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/03//align_stats-1.png)<!-- -->

### TABLE: Human bulk QC

Export a supplementary table for human bulk RNA-seq QC.


```r
align_stats %>% 
    select(-ID, -group) %>% 
    dplyr::rename(ID = name) %>% 
    dplyr::relocate(Group, .after = 1) %>% 
    rr_write_tsv(glue("{out}/TABLE_bulk_stats.tsv"),
                 "QC stats for bulk RNAseq data")
```

```
## ...writing description of TABLE_bulk_stats.tsv to public/output/03/TABLE_bulk_stats.desc
```

## PCA


```r
# load the counts saved by the pipeline
load(here("data/RNAseq/pipeline_l3/2023-05-test_pbt/.counts.RData"))

# get highly-variable genes
var_idx <- apply(counts$norm$Ensembl.ensGene.exon, 1, var) %>% 
    order(decreasing = TRUE) %>% 
    .[1:1000]

pbt_counts_vst_hvg <- counts$vst$Ensembl.ensGene.exon[var_idx, ]

pca       <- prcomp(t(pbt_counts_vst_hvg))
variance  <- round(pca$sdev^2 / sum(pca$sdev^2), 3) * 100

pca_df <- data.frame(Dim1 = pca$x[, 1],
                     Dim2 = pca$x[, 2],
                     Sample = colnames(pbt_counts_vst_hvg))

p1 <- pca_df %>% 
    left_join(meta_bulk, by = c("Sample" = "ID")) %>% 
    ggplot(aes(x = Dim1, y = Dim2)) +
    geom_point(aes(colour = Group), alpha = 0.7, size = 2) +
    scale_color_manual(values = palette_groups) +
    xlab(paste0("PC1 (", variance[1], "%)")) +
    ylab(paste0("PC2 (", variance[2], "%)")) +
    color_legend_ncol(2) +
    ggtitle("PCA")

p2 <- pca_df %>% 
    left_join(meta_bulk, by = c("Sample" = "ID")) %>% 
    ggplot(aes(x = Dim1, y = Dim2)) +
    geom_point(aes(colour = Source), alpha = 0.7, size = 2) +
    xlab(paste0("PC1 (", variance[1], "%)")) +
    ylab(paste0("PC2 (", variance[2], "%)")) +
    color_legend_ncol(2) +
    ggtitle("PCA")
```

Coloured by QC stats:


```r
pca_with_align_stats <- pca_df %>% 
    left_join(info_samples, by = c("Sample" = "Nickname")) %>% 
    left_join(align_stats, by = c("Sample" = "name", "Group" = "group"))

pca_stat_plots <- list()

stats <- c("cleanReadCount",
           "mappedReadCount",
           "exonReadCount",
           "intronReadCount",
           "mitochondrialPercentage",
           "ribosomalPercentage")

for (i in seq_along(stats)) {
    
    df <- pca_with_align_stats
    df$stat <- df[, stats[i]]
    
    gg <- df %>% 
        arrange(stat) %>% 
        ggplot(aes(x = Dim1, y = Dim2)) +
        geom_point(aes(colour = stat), alpha = 0.7, size = 2) +
        scale_colour_viridis(labels = comma) +
        xlab(paste0("PC1 (", variance[1], "%)")) +
        ylab(paste0("PC2 (", variance[2], "%)")) +
        ggtitle(stats[i])
    
    pca_stat_plots[[i]] <- gg
    
}

plot_grid(plotlist = c(list(p1, p2), pca_stat_plots), align = "hv", axis = "rltb", ncol = 3)
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/03//pca_align_stats-1.png)<!-- -->

## FOXR2 expression


```r
p1 <- counts_foxr2 %>%
    ggplot(aes(x = Group, y = gene_expression)) +
    geom_boxplot(aes(fill = Group), outlier.shape = NA) +
    geom_jitter(size = 2, alpha = 0.8) +
    scale_fill_manual(values = palette_groups)  +
    ylab("Normalized expression") + xlab("Group") +
    no_legend() + rotate_x() +
    ggtitle("FOXR2 expression")

p2 <- counts_foxr2 %>% 
    filter(Group == "NB-FOXR2") %>% 
    ggplot(aes(x = sample, y = gene_expression)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    ggtitle("FOXR2 expression per sample")

plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb")
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/03//boxplot_Foxr2-1.png)<!-- -->

Check FOXR2 expression across samples:


```r
counts_foxr2 %>% 
    arrange(Group) %>% 
    mutate(sample = factor(sample, levels = unique(.$sample))) %>% 
    ggplot(aes(x = sample, y = gene_expression)) +
    geom_bar(stat = "identity", aes(fill = Group)) +
    scale_fill_manual(values = palette_groups) +
    ggtitle("FOXR2 expression per sample") +
    rotate_x()
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/03//foxr2_expression_all-1.png)<!-- -->

## MDM4 expression

The ubiquitous loss of chr1 in NB-FOXR2 may lead to overexpression of MDM2/MDM4, which phenocopies p53 loss in the tumor context ([Girish et al. *Science* 2023](https://pubmed.ncbi.nlm.nih.gov/37410869/)).

Here, we investigate this possibility by comparing expression of MDM4 across the bulk intracranial tumor cohort.


```r
# load counts for MDM genes
counts_MDM <- extract_pipeline_counts(path = here("data/RNAseq/pipeline_l3/2023-05-test_pbt/counts/Ensembl.ensGene.exon.norm.tsv.gz"),
                                        goi = c("MDM4")) %>% 
    left_join(info_samples, by = c("sample" = "Nickname")) %>% 
    mutate(Group = factor(Group, levels = names(palette_groups)))
```

```
## Joining with `by = join_by(gene_symbol)`
```

```r
palette_highlight_nbfoxr2 <- rep("grey70", length(palette_groups))
names(palette_highlight_nbfoxr2) <- names(palette_groups)
palette_highlight_nbfoxr2[which(names(palette_highlight_nbfoxr2) == "NB-FOXR2")] <- "red"

counts_MDM %>%
    filter(gene_symbol == "MDM4") %>% 
    ggplot(aes(x = Group, y = gene_expression, fill = Group)) +
    geom_violin(scale = "width") +
    scale_fill_manual(values = palette_highlight_nbfoxr2)  +
    ylab("Normalized expression") + xlab("Group") +
    no_legend() + rotate_x() +
    ggtitle("MDM4 expression") +
    stat_summary(fun.y=median, geom="crossbar", size=1, color="black", aes(width = 0.3))
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/03//mdm4_expression_all-1.png)<!-- -->

# Extra-cranial neuroblastoma

## Loading and data wrangling

Load processed bulk RNAseq counts from [Gartlgruber et al, Nature Cancer, 2020](https://www.nature.com/articles/s43018-020-00145-w#Abs1).

The data was downloaded via their shiny app (https://nbseb087.dkfz.de/project_NB_SE_viz/).


```r
gartlgruber_meta <- readRDS(here("data/RNAseq/external_data/Gartlgruber_NatCancer_2021/annotation_tumor_RNAseq.rds"))
length(unique(gartlgruber_meta$ProjectID))
```

```
## [1] 579
```

```r
# load annotation of ChIPseq data to look for overlapping samples
gartlgruber_meta_chip <- readRDS(here("data/ChIPseq/external_data/Gartlgruber_NatCancer_2021/annotation_tumor_ChIPseq.rds"))
gartlgruber_meta <- gartlgruber_meta %>% 
    left_join(gartlgruber_meta_chip, by = c("ProjectID", "MYCN", "Stage",
                                            "Age", "Risk", "Relapse"))

gartlgruber_counts <- readRDS(here("data/RNAseq/external_data/Gartlgruber_NatCancer_2021/tumor_RNAseq_Counts_Matrix.rds"))
dim(gartlgruber_counts)
```

```
## [1] 57820   579
```

```r
gartlgruber_counts[1:5, 1:5]
```

```
##                              NSP001-PT01 NSP002-PT01 NSP003-RM01 NSP005-RT01
## ENSG00000223972.4|DDX11L1              4           0           0           3
## ENSG00000227232.4|WASH7P             431        1348        1242        1117
## ENSG00000243485.2|MIR1302-11           0           0           0           0
## ENSG00000237613.2|FAM138A              0           0           0           0
## ENSG00000268020.2|OR4G4P               0           0           0           0
##                              NSP006-RM01
## ENSG00000223972.4|DDX11L1              1
## ENSG00000227232.4|WASH7P            1391
## ENSG00000243485.2|MIR1302-11           0
## ENSG00000237613.2|FAM138A              0
## ENSG00000268020.2|OR4G4P               0
```

Normalize and reformat:


```r
# normalize using DESeq2
dds <- DESeqDataSetFromMatrix(countData = gartlgruber_counts,
                              colData   = gartlgruber_meta,
                              design    = ~ 1)
dds <- DESeq(dds)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```

```
## -- replacing outliers and refitting for 6774 genes
## -- DESeq argument 'minReplicatesForReplace' = 7 
## -- original counts are preserved in counts(dds)
```

```
## estimating dispersions
```

```
## fitting model and testing
```

```r
ecnb_counts_norm <- counts(dds, normalized = TRUE)
vsd  <- vst(dds)
ecnb_counts_vst <- assay(vsd)
```

Convert to tidy format with metadata:


```r
ecnb_counts_tidy <- ecnb_counts_norm %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "id") %>% 
    separate(id, into = c("ens_id.x", "gene_symbol"), sep = "\\|") %>% 
    separate(ens_id.x, into = c("ens_id", "drop"), sep = "\\.") %>% 
    select(-drop) %>% 
    gather(sample, gene_expression, 3:ncol(.)) %>% 
    left_join(gartlgruber_meta, by = c("sample" = "ProjectID"))

save(ecnb_counts_tidy, ecnb_counts_norm, ecnb_counts_vst, file = glue("{out}/Gartlgruber_et_al_counts.Rda"))
```


<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2024-11-01 11:50:05
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
##  ! package              * version  date (UTC) lib source
##  P annotate               1.72.0   2021-10-26 [?] Bioconductor
##  P AnnotationDbi          1.56.2   2021-11-09 [?] Bioconductor
##  P Biobase              * 2.54.0   2021-10-26 [?] Bioconductor
##  P BiocGenerics         * 0.40.0   2021-10-26 [?] Bioconductor
##  P BiocManager            1.30.15  2021-05-11 [?] CRAN (R 4.1.2)
##  P BiocParallel           1.28.3   2021-12-09 [?] Bioconductor
##  P Biostrings             2.62.0   2021-10-26 [?] Bioconductor
##  P bit                    4.0.4    2020-08-04 [?] CRAN (R 4.1.2)
##  P bit64                  4.0.5    2020-08-30 [?] CRAN (R 4.1.2)
##  P bitops                 1.0-7    2021-04-24 [?] CRAN (R 4.1.2)
##  P blob                   1.2.2    2021-07-23 [?] CRAN (R 4.1.2)
##  P bslib                  0.3.1    2021-10-06 [?] CRAN (R 4.1.2)
##  P cachem                 1.0.6    2021-08-19 [?] CRAN (R 4.1.2)
##  P callr                  3.7.6    2024-03-25 [?] RSPM
##  P cellranger             1.1.0    2016-07-27 [?] CRAN (R 4.1.2)
##  P cli                    3.6.1    2023-03-23 [?] RSPM (R 4.1.2)
##  P codetools              0.2-18   2020-11-04 [?] CRAN (R 4.1.2)
##  P colorspace             2.0-2    2021-06-24 [?] CRAN (R 4.1.2)
##  P cowplot              * 1.1.1    2020-12-30 [?] CRAN (R 4.1.2)
##  P crayon                 1.4.2    2021-10-29 [?] CRAN (R 4.1.2)
##  P DBI                    1.1.2    2021-12-20 [?] CRAN (R 4.1.2)
##  P DelayedArray           0.20.0   2021-10-26 [?] Bioconductor
##  P DESeq2               * 1.34.0   2021-10-26 [?] Bioconductor
##  P devtools               2.4.5    2022-10-11 [?] CRAN (R 4.1.2)
##  P digest                 0.6.35   2024-03-11 [?] CRAN (R 4.1.2)
##  P dplyr                * 1.1.1    2023-03-22 [?] CRAN (R 4.1.2)
##  P ellipsis               0.3.2    2021-04-29 [?] CRAN (R 4.1.2)
##  P evaluate               0.23     2023-11-01 [?] CRAN (R 4.1.2)
##  P fansi                  1.0.2    2022-01-14 [?] CRAN (R 4.1.2)
##  P farver                 2.1.0    2021-02-28 [?] CRAN (R 4.1.2)
##  P fastmap                1.1.0    2021-01-25 [?] CRAN (R 4.1.2)
##  P fs                     1.5.2    2021-12-08 [?] CRAN (R 4.1.2)
##  P genefilter             1.76.0   2021-10-26 [?] Bioconductor
##  P geneplotter            1.72.0   2021-10-26 [?] Bioconductor
##  P generics               0.1.3    2022-07-05 [?] CRAN (R 4.1.2)
##  P GenomeInfoDb         * 1.30.1   2022-01-30 [?] Bioconductor
##  P GenomeInfoDbData       1.2.4    2023-11-28 [?] Bioconductor
##  P GenomicRanges        * 1.46.1   2021-11-18 [?] Bioconductor
##  P ggplot2              * 3.4.2    2023-04-03 [?] CRAN (R 4.1.2)
##  P git2r                  0.29.0   2021-11-22 [?] CRAN (R 4.1.2)
##  P glue                 * 1.6.2    2022-02-24 [?] CRAN (R 4.1.2)
##  P gridExtra              2.3      2017-09-09 [?] CRAN (R 4.1.2)
##  P gtable                 0.3.0    2019-03-25 [?] CRAN (R 4.1.2)
##  P here                 * 1.0.1    2020-12-13 [?] CRAN (R 4.1.2)
##  P highr                  0.9      2021-04-16 [?] CRAN (R 4.1.2)
##  P hms                    1.1.1    2021-09-26 [?] CRAN (R 4.1.2)
##  P htmltools              0.5.2    2021-08-25 [?] CRAN (R 4.1.2)
##  P htmlwidgets            1.5.4    2021-09-08 [?] CRAN (R 4.1.2)
##  P httpuv                 1.6.5    2022-01-05 [?] CRAN (R 4.1.2)
##  P httr                   1.4.2    2020-07-20 [?] CRAN (R 4.1.2)
##  P IRanges              * 2.28.0   2021-10-26 [?] Bioconductor
##  P jquerylib              0.1.4    2021-04-26 [?] CRAN (R 4.1.2)
##  P jsonlite               1.8.8    2023-12-04 [?] CRAN (R 4.1.2)
##  P KEGGREST               1.34.0   2021-10-26 [?] Bioconductor
##  P knitr                  1.37     2021-12-16 [?] CRAN (R 4.1.2)
##  P labeling               0.4.2    2020-10-20 [?] CRAN (R 4.1.2)
##  P later                  1.3.0    2021-08-18 [?] CRAN (R 4.1.2)
##  P lattice                0.20-45  2021-09-22 [?] CRAN (R 4.1.2)
##  P lifecycle              1.0.3    2022-10-07 [?] CRAN (R 4.1.2)
##  P locfit                 1.5-9.4  2020-03-25 [?] CRAN (R 4.1.2)
##  P magrittr             * 2.0.3    2022-03-30 [?] CRAN (R 4.1.2)
##  P Matrix                 1.3-4    2021-06-01 [?] CRAN (R 4.1.2)
##  P MatrixGenerics       * 1.6.0    2021-10-26 [?] Bioconductor
##  P matrixStats          * 0.61.0   2021-09-17 [?] CRAN (R 4.1.2)
##  P memoise                2.0.1    2021-11-26 [?] CRAN (R 4.1.2)
##  P mime                   0.12     2021-09-28 [?] CRAN (R 4.1.2)
##  P miniUI                 0.1.1.1  2018-05-18 [?] CRAN (R 4.1.2)
##  P munsell                0.5.0    2018-06-12 [?] CRAN (R 4.1.2)
##  P pillar                 1.9.0    2023-03-22 [?] RSPM (R 4.1.2)
##  P pkgbuild               1.4.2    2023-06-26 [?] CRAN (R 4.1.2)
##  P pkgconfig              2.0.3    2019-09-22 [?] CRAN (R 4.1.2)
##  P pkgload                1.3.3    2023-09-22 [?] CRAN (R 4.1.2)
##  P plyr                   1.8.6    2020-03-03 [?] CRAN (R 4.1.2)
##  P png                    0.1-7    2013-12-03 [?] CRAN (R 4.1.2)
##  P prettyunits            1.1.1    2020-01-24 [?] CRAN (R 4.1.2)
##  P processx               3.8.4    2024-03-16 [?] RSPM
##  P profvis                0.3.8    2023-05-02 [?] CRAN (R 4.1.2)
##  P promises               1.2.0.1  2021-02-11 [?] CRAN (R 4.1.2)
##  P ps                     1.7.6    2024-01-18 [?] RSPM
##  P purrr                * 1.0.1    2023-01-10 [?] CRAN (R 4.1.2)
##  P R6                     2.5.1    2021-08-19 [?] CRAN (R 4.1.2)
##  P RColorBrewer         * 1.1-2    2014-12-07 [?] CRAN (R 4.1.2)
##  P Rcpp                   1.0.8    2022-01-13 [?] CRAN (R 4.1.2)
##  P RCurl                  1.98-1.5 2021-09-17 [?] CRAN (R 4.1.2)
##  P readr                * 2.1.1    2021-11-30 [?] CRAN (R 4.1.2)
##  P readxl               * 1.3.1    2019-03-13 [?] CRAN (R 4.1.2)
##  P remotes                2.4.2.1  2023-07-18 [?] CRAN (R 4.1.2)
##  P renv                   1.0.3    2023-09-19 [?] CRAN (R 4.1.2)
##  P reshape2               1.4.4    2020-04-09 [?] CRAN (R 4.1.2)
##  P rlang                  1.1.3    2024-01-10 [?] CRAN (R 4.1.2)
##  P rmarkdown              2.11     2021-09-14 [?] CRAN (R 4.1.2)
##  P rprojroot              2.0.2    2020-11-15 [?] CRAN (R 4.1.2)
##  P RSQLite                2.2.9    2021-12-06 [?] CRAN (R 4.1.2)
##  P S4Vectors            * 0.32.4   2022-03-24 [?] Bioconductor
##  P sass                   0.4.0    2021-05-12 [?] CRAN (R 4.1.2)
##  P scales               * 1.2.1    2022-08-20 [?] CRAN (R 4.1.2)
##  P sessioninfo            1.2.2    2021-12-06 [?] CRAN (R 4.1.2)
##  P shiny                  1.7.1    2021-10-02 [?] CRAN (R 4.1.2)
##  P stringi                1.7.6    2021-11-29 [?] CRAN (R 4.1.2)
##  P stringr                1.5.0    2022-12-02 [?] CRAN (R 4.1.2)
##  P SummarizedExperiment * 1.24.0   2021-10-26 [?] Bioconductor
##  P survival               3.2-13   2021-08-24 [?] CRAN (R 4.1.2)
##  P tibble                 3.2.1    2023-03-20 [?] RSPM (R 4.1.2)
##  P tidyr                * 1.3.0    2023-01-24 [?] CRAN (R 4.1.2)
##  P tidyselect             1.2.0    2022-10-10 [?] CRAN (R 4.1.2)
##  P tzdb                   0.3.0    2022-03-28 [?] CRAN (R 4.1.2)
##  P urlchecker             1.0.1    2021-11-30 [?] CRAN (R 4.1.2)
##  P usethis                2.2.2    2023-07-06 [?] CRAN (R 4.1.2)
##  P utf8                   1.2.2    2021-07-24 [?] CRAN (R 4.1.2)
##  P vctrs                  0.6.5    2023-12-01 [?] CRAN (R 4.1.2)
##  P viridis              * 0.5.1    2018-03-29 [?] RSPM (R 4.1.2)
##  P viridisLite          * 0.3.0    2018-02-01 [?] CRAN (R 4.1.2)
##  P vroom                  1.5.7    2021-11-30 [?] CRAN (R 4.1.2)
##  P withr                  2.5.0    2022-03-03 [?] CRAN (R 4.1.2)
##  P xfun                   0.29     2021-12-14 [?] CRAN (R 4.1.2)
##  P XML                    3.99-0.8 2021-09-17 [?] CRAN (R 4.1.2)
##  P xtable                 1.8-4    2019-04-21 [?] CRAN (R 4.1.2)
##  P XVector                0.34.0   2021-10-26 [?] Bioconductor
##  P yaml                   2.2.1    2020-02-01 [?] CRAN (R 4.1.2)
##  P zlibbioc               1.40.0   2021-10-26 [?] Bioconductor
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
