---
title: "09 - FOXR2 CUT&RUN in FOXR2 p53LOF mouse models"
author: "Bhavyaa Chandarana [[bhavyaa.chandarana@mail.mcgill.ca](mailto:bhavyaa.chandarana@mail.mcgill.ca)] and Steven Hébert [[steven.hebert@ladydavis.ca](mailto:steven.hebert@ladydavis.ca)]"
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
## Document index: 09
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
## public/output/09
```

```
## public/figures/09//
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

In this document, we analyze FOXR2 CUT&RUN (similar to ChIP) performed in n=3 
cell lines derived from Foxr2 p53 LOF murine models.

Note that throughout the code we refer to the CUT&RUN data as "ChIP" due
to its similarity in the analysis stage.

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
library(ComplexHeatmap)
library(GenomeInfoDb)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(Seurat)
library(Signac)
library(gprofiler2)
library(VennDiagram)

source(here("include/style.R"))
source(here("code/functions/scRNAseq.R"))
source(here("code/functions/RNAseq.R"))

ggplot2::theme_set(theme_min())
```

# Load data

FOXR2-V5 CUT&RUN (similar to ChIP) was run with two protocols "native" and "X-link",
with n=2 replicates each. 

A FOXR2-V5 experiment was also run in a cell line lacking V5-tagged insertion (n=2 
replicates) which we use here to filter false positive peaks.


```r
chip_path <- here("data/ChIPseq/pipeline_genpipes/peaks")

# X-link FOXR2 peaks
rep1_FOXR2_peaks <- rtracklayer::import(glue("{chip_path}/AN24377_FOXR2_P53LOF_1_X-link_JAB1979A5-1.V5_peaks.narrowPeak"), format = "narrowPeak") 
rep2_FOXR2_peaks <- rtracklayer::import(glue("{chip_path}/AN24377_FOXR2_P53LOF_2_X-link_JAB1979A5-2.V5_peaks.narrowPeak"), format = "narrowPeak")

# X-link V5 peaks
rep1_V5_peaks <- rtracklayer::import(glue("{chip_path}/AH25351_G34R_1_X-link_JAB1979A5-3.V5_peaks.narrowPeak"), format = "narrowPeak")
rep2_V5_peaks <- rtracklayer::import(glue("{chip_path}/AH25351_G34R_2_X-link_JAB1979A5-4.V5_peaks.narrowPeak"), format = "narrowPeak")

# Create lists of peaks
chip_foxr2_peak_list <- list(rep1_FOXR2 = rep1_FOXR2_peaks,
               rep2_FOXR2 = rep2_FOXR2_peaks)

chip_V5_peak_list <- list(rep1_V5 = rep1_V5_peaks,
               rep2_V5 = rep2_V5_peaks)

# Remove non-standard chromosomes/contigs from each sample
standard_chrom_mm <- function(gr) { 
    gr <- gr[GenomicRanges::seqnames(gr) %in% c(1:19, 'X', 'Y')]
    GenomeInfoDb::seqlevels(gr) <- c(1:19, 'X', 'Y')
    return(gr)
}

chip_foxr2_peak_list <- map(chip_foxr2_peak_list, standard_chrom_mm)
chip_V5_peak_list <- map(chip_V5_peak_list, standard_chrom_mm)
```

# Define peaks

## Filter V5 peaks

We have control data from FOXR2-V5 CUT&RUN in a cell line lacking V5 insertion.
The peaks in this sample arise from binding regions of V5 antibody
in the genome. We use this as a background to filter out false positive peaks.


```r
# Concatenate the V5 peaks 
all_V5_peaks <- c(chip_V5_peak_list[[1]], chip_V5_peak_list[[2]]) %>%
  sortSeqlevels %>%
  sort %>%
  reduce(min.gapwidth = 0)

# Remove V5 peaks from the FOXR2 peaks
chip_foxr2_peak_list_no_V5 <- chip_foxr2_peak_list
chip_foxr2_peak_list_no_V5[[1]] <- subsetByOverlaps(chip_foxr2_peak_list_no_V5[[1]], all_V5_peaks, invert = T, minoverlap = 10)
chip_foxr2_peak_list_no_V5[[2]] <- subsetByOverlaps(chip_foxr2_peak_list_no_V5[[2]], all_V5_peaks, invert = T, minoverlap = 10)

# Check the length of the peak lists (as a QC)
map(chip_foxr2_peak_list, length)
```

```
## $rep1_FOXR2
## [1] 31861
## 
## $rep2_FOXR2
## [1] 38522
```

```r
map(chip_foxr2_peak_list_no_V5, length)
```

```
## $rep1_FOXR2
## [1] 31718
## 
## $rep2_FOXR2
## [1] 38393
```

```r
# Create a narrowpeak file for the ChIP-seq FOXR2 samples after removing the V5 peaks 
write.table(chip_foxr2_peak_list_no_V5[[1]], file = glue("{out}/rep1_FOXR2_no_V5.narrowpeak"), sep = "\t", quote = F, row.names = F, col.names = T)
write.table(chip_foxr2_peak_list_no_V5[[2]], file = glue("{out}/rep2_FOXR2_no_V5.narrowpeak"), sep = "\t", quote = F, row.names = F, col.names = T)
```

## Intersect FOXR2 CUT&RUN & ATAC-seq peaks

We expect high-quality FOXR2 CUT&RUN peaks to overlap with open chromatin
(ATAC-seq peaks) derived from scMultiome of the same mouse model.

Here, we compute the overlap between these peaks, and use ATAC peaks
to filter CUT&RUN peaks. In the end, we retain the intersect of peaks between
ATAC-seq, and the two replicates of X-link FOXR2 CUT&RUN.

Define a function to plot triple venn diagram:


```r
# Function to plot triple venn diagram from a list containing 3 GRanges objects
venn_from_gr_list <- function(gr_list, plot_title){
    
  overlap_all <- Reduce(function(x, y) subsetByOverlaps(x, y, minoverlap = 10), gr_list) %>% length()
  
  overlap_12 <- findOverlaps(gr_list[[1]], gr_list[[2]], minoverlap = 10)
  overlap_23 <- findOverlaps(gr_list[[2]], gr_list[[3]], minoverlap = 10)
  overlap_13 <- findOverlaps(gr_list[[1]], gr_list[[3]], minoverlap = 10)
  
  n_overlap_12 <- min(length(unique(queryHits(overlap_12))), length(unique(subjectHits(overlap_12))))
  n_overlap_23 <- min(length(unique(queryHits(overlap_23))), length(unique(subjectHits(overlap_23))))
  n_overlap_13 <- min(length(unique(queryHits(overlap_13))), length(unique(subjectHits(overlap_13))))
  
  venn.plot = draw.triple.venn(length(gr_list[[1]]),
                 length(gr_list[[2]]),
                 length(gr_list[[3]]),
                 n_overlap_12,
                 n_overlap_23,
                 n_overlap_13,
                 overlap_all,
                 category = c("rep1_FOXR2", "rep2_FOXR2", "ATAC"),
                 euler.d = T,
                 scaled = T,
                 fill = c("#33CC00", "#33CC00", "#009999"))

  grid.draw(venn.plot)
}
```

Load ATAC peaks:


```r
# Load ATAC peaks produced from multiome processing pipeline
Foxr2_p53_r1_ATAC_peaks <- rtracklayer::import(here("data/singlecell/pipeline_scMultiome_mm/AN24377/ATAC/macs2.peaks.standard.blacklistremoved.narrowpeak"), format = "narrowPeak")

# Note that the blacklisted peaks were removed from list of peaks
# The list of blacklisted peaks comes from this set (from Signac package)
write.table(Signac::blacklist_mm10, file = glue("{out}/blacklist_mm10.tsv"), 
            sep = "\t", quote = F, row.names = F, col.names = F)

# Modify the blacklist_mm10 file to get a sorted bed files with names
system(glue("sort -Vk1,1 -k2,2 -k3,3 {out}/blacklist_mm10.tsv | 
               cut -f1-3,6 |
               sed 's/ /_/g' > {out}/blacklist_mm10.bed"))

# Remove non-standard chromosomes
Foxr2_p53_r1_ATAC_peaks <- standard_chrom_mm(Foxr2_p53_r1_ATAC_peaks)

### Check number of peaks overlapping without any filtering
peak_list <- c(chip_foxr2_peak_list_no_V5, FOXR2_ATAC = Foxr2_p53_r1_ATAC_peaks)
venn_from_gr_list(peak_list, "venn_ChIP_ATAC_nofilter")
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/09//compare-chip-atac-1.png)<!-- -->

```r
# Add width column to the peaks
chip_foxr2_peak_list_no_V5 <- map(chip_foxr2_peak_list_no_V5, ~ { 
  .x$width <- width(.x)
  return(.x)
})

Foxr2_p53_r1_ATAC_peaks$width <- width(Foxr2_p53_r1_ATAC_peaks)
```

Filter CUT&RUN and ATAC peaks by quality metrics:


```r
# Filter peaks of ChIP-seq and ATAC-seq
# Peaks from ChIP-seq and ATAC-seq have different value distribution 
# (signalValue, qValue), so different thresholds need to be used.
filtered_2_chip_foxr2_peak_list_no_V5 <- map(chip_foxr2_peak_list_no_V5,
        ~ .x[which(.x$qValue > 6 & .x$signalValue > 6 & .x$width < 10000),])
filtered_2_Foxr2_p53_r1_ATAC_peaks <- Foxr2_p53_r1_ATAC_peaks[Foxr2_p53_r1_ATAC_peaks$qValue > 4 &
                                                            Foxr2_p53_r1_ATAC_peaks$signalValue > 2 &
                                                            Foxr2_p53_r1_ATAC_peaks$width < 10000,]

### Check number of peaks overlapping
filtered_2_peak_list <- c(filtered_2_chip_foxr2_peak_list_no_V5, FOXR2_ATAC = filtered_2_Foxr2_p53_r1_ATAC_peaks)
venn_from_gr_list(filtered_2_peak_list, "venn_ChIP_ATAC_filtered_stringent")
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/09//filter-chip-atac-1.png)<!-- -->

```r
# export the peaks into bed files
imap(filtered_2_peak_list[1:2], ~export.bed(.x, 
                                     con=glue("{out}/{.y}_qVal6_FC6_width10k.bed"), 
                                     format = "bed")) 
```

```
## $rep1_FOXR2
## BEDFile object
## resource: /project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/output/09/rep1_FOXR2_qVal6_FC6_width10k.bed 
## 
## $rep2_FOXR2
## BEDFile object
## resource: /project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/output/09/rep2_FOXR2_qVal6_FC6_width10k.bed
```

```r
export.bed(filtered_2_peak_list[[3]], 
         con=glue("{out}/FOXR2_ATAC_qVal4_FC2_width10k.bed"), 
         format = "bed")
```

Take the intersecting peaks between the n=2 FOXR2 CUT&RUN and ATAC:


```r
# Create a bed file of the overlap between the 3 sets of peaks (ChIP FOXR2 rep 1&2, ATAC FOXR2-p53 rep1) 
overlap_all <- Reduce(function(x, y) subsetByOverlaps(x, y, minoverlap = 10), filtered_2_peak_list)
export.bed(overlap_all, con = glue("{out}/intersect_ChIP_ATAC_qVal6_FC6_width10k.bed"), format = "bed")

# Create the intersect with the exact bases that intersect the ChIP FOXR2 rep 1&2 and ATAC FOXR2-p53 rep1
# Run after: module load bedtools/2.30.0
system(glue("intersectBed -a {out}/rep1_FOXR2_qVal6_FC6_width10k.bed -b {out}/rep2_FOXR2_qVal6_FC6_width10k.bed > {out}/exact_ChIP-intersect_qVal6_FC6_width10k.bed"))
system(glue("intersectBed -a {out}/exact_ChIP-intersect_qVal6_FC6_width10k.bed -b {out}/FOXR2_ATAC_qVal4_FC2_width10k.bed > {out}/exact_intersect_ChIP_ATAC_qVal6_FC6_width10k.bed"))
```

## Peak QC metric distribution

Plot histograms of metrics (signal, p-val, q-val) distributions across
CUT&RUN peaks, with QC thresholds for each metric as a vertical line.


```r
qval_threshold <- c(6, 6, 4)
fold_threshold <- c(6, 6, 2)
width_threshold1 <- 1
width_threshold2 <- 10000

# Q values
# log-log scale
map2(names(peak_list), qval_threshold,
     ~ data.frame(num = peak_list[[.x]]$name,
                  qValue = peak_list[[.x]]$qValue) %>% 
            ggplot(aes(x = log2(qValue))) +
            geom_histogram(bins = 100) +
            xlab("log2(-log10(qValue))") +
            ggtitle(glue("{.x} - q-value")) +
            geom_vline(xintercept = log2(.y), color = "red"))  %>% 
    cowplot::plot_grid(plotlist = ., ncol = 3)
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/09//histo-peaks-1.png)<!-- -->

```r
# Fold enrichment
map2(names(peak_list), fold_threshold,
     ~ data.frame(num = peak_list[[.x]]$name,
                  fold_enrichment = peak_list[[.x]]$signalValue) %>% 
            ggplot(aes(x = log2(fold_enrichment))) +
            geom_histogram(bins = 100) +
            xlab("log2(fold_enrichment)") +
            ggtitle(glue("{.x} - signal")) +
            geom_vline(xintercept = log2(.y), color = "red")) %>% 
    cowplot::plot_grid(plotlist = ., ncol = 3)
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/09//histo-peaks-2.png)<!-- -->

```r
# Peak width

width_df <- function(gr){
    
    df <- data.frame(num = gr$name)
    df$width <- width(gr)
    return(df)
    
}

imap(peak_list,
     ~ width_df(.x) %>% 
            ggplot(aes(x = width)) +
            geom_histogram(bins = 300) +
            xlab("peak width") +
            ggtitle(glue("{.y} - width")) +
            geom_vline(xintercept = width_threshold1, color = "red") +
            geom_vline(xintercept = width_threshold2, color = "red"))  %>% 
    cowplot::plot_grid(plotlist = ., ncol = 3)
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/09//histo-peaks-3.png)<!-- -->

# Motif analysis

Motif analysis with HOMER was conducted with scripts in the following directory:

```
code/scripts/ChIP
```

Originally written and run by Steven Hébert.


```r
# Load Homer motif results
motifs <- read_tsv(glue("{out}/motifs_size100_100/knownResults.txt"))
```

```
## Rows: 440 Columns: 9
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (4): Motif Name, Consensus, % of Target Sequences with Motif, % of Backg...
## dbl (5): P-value, Log P-value, q-value (Benjamini), # of Target Sequences wi...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
# Adjust p-values and count the number of significant motifs
motifs_adj <- motifs %>%
  mutate(padj_BH = p.adjust(`P-value`, method = "BH")) 
nrow(filter(motifs_adj, padj_BH < 0.05))
```

```
## [1] 193
```

Plot the top motifs and their -log10 adjusted p-value, coloring by the motif type
(ETS, SOX, other)


```r
motifs_gg <- motifs %>%
  mutate(padj_BH = p.adjust(`P-value`, method = "BH")) %>%
  slice_max(order_by = -`Log P-value`, n = 13, with_ties = T) %>%
  mutate(Motif_name = str_split(`Motif Name`, "/", simplify = T)[,1]) %>% 
  mutate(Motif_group = case_when(grepl(.$Motif_name, pattern = "ETS") ~ "ETS",
                                 grepl(.$Motif_name, pattern = "Sox") ~ "SOX",
                                 TRUE                                        ~ "Other"))

motifs_gg %>% 
  ggplot(aes(x = -log10(padj_BH), y = reorder(Motif_name, -padj_BH))) +
  geom_bar(stat = "identity", aes(fill = Motif_group)) +
  ylab("Motif") +
  xlab("-log10(padj)") +
  scale_fill_manual(values = c("ETS" = "darkgoldenrod1",
                               "SOX" = "darkolivegreen3",
                               "Other" = "grey70")) +
  theme(legend.position = "none")
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/09//plot-motifs-1.png)<!-- -->

Reorder the bars to group the motif types together, and rank by adjusted p-value
within each of the groups.


```r
df <- motifs_gg %>% 
    arrange(padj_BH) %>% 
    dplyr::select(Motif_name, Motif_group) 

order <- c(df %>% filter(Motif_group == "ETS") %>% .$Motif_name,
           df %>% filter(Motif_group == "SOX") %>% .$Motif_name,
           df %>% filter(Motif_group == "Other") %>% .$Motif_name)

motifs_gg %>% 
  mutate(Motif_name = factor(Motif_name, levels = rev(order))) %>% 
  ggplot(aes(x = -log10(padj_BH), y = Motif_name)) +
  geom_bar(stat = "identity", aes(fill = Motif_group)) +
  ylab("Motif") +
  xlab("-log10(padj)") +
  scale_fill_manual(values = c("ETS" = "darkgoldenrod1",
                               "SOX" = "darkolivegreen3",
                               "Other" = "grey70")) +
  theme(legend.position = "none")
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/09//plot-motifs-reorder-1.png)<!-- -->


<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2024-11-04 12:26:35
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
##  ! package              * version  date (UTC) lib source
##  P abind                  1.4-5    2016-07-21 [?] CRAN (R 4.1.2)
##  P AnnotationDbi        * 1.56.2   2021-11-09 [?] Bioconductor
##  P assertthat             0.2.1    2019-03-21 [?] CRAN (R 4.1.2)
##  P Biobase              * 2.54.0   2021-10-26 [?] Bioconductor
##  P BiocFileCache          2.2.1    2022-01-23 [?] Bioconductor
##  P BiocGenerics         * 0.40.0   2021-10-26 [?] Bioconductor
##  P BiocIO                 1.4.0    2021-10-26 [?] Bioconductor
##  P BiocManager            1.30.15  2021-05-11 [?] CRAN (R 4.1.2)
##  P BiocParallel           1.28.3   2021-12-09 [?] Bioconductor
##  P biomaRt                2.50.2   2022-01-13 [?] Bioconductor
##  P Biostrings             2.62.0   2021-10-26 [?] Bioconductor
##  P bit                    4.0.4    2020-08-04 [?] CRAN (R 4.1.2)
##  P bit64                  4.0.5    2020-08-30 [?] CRAN (R 4.1.2)
##  P bitops                 1.0-7    2021-04-24 [?] CRAN (R 4.1.2)
##  P blob                   1.2.2    2021-07-23 [?] CRAN (R 4.1.2)
##  P bslib                  0.3.1    2021-10-06 [?] CRAN (R 4.1.2)
##  P cachem                 1.0.6    2021-08-19 [?] CRAN (R 4.1.2)
##  P callr                  3.7.6    2024-03-25 [?] RSPM
##  P cellranger             1.1.0    2016-07-27 [?] CRAN (R 4.1.2)
##  P circlize               0.4.15   2022-05-10 [?] CRAN (R 4.1.2)
##  P cli                    3.6.1    2023-03-23 [?] RSPM (R 4.1.2)
##  P clue                   0.3-64   2023-01-31 [?] CRAN (R 4.1.2)
##  P cluster                2.1.2    2021-04-17 [?] CRAN (R 4.1.2)
##  P codetools              0.2-18   2020-11-04 [?] CRAN (R 4.1.2)
##  P colorspace             2.0-2    2021-06-24 [?] CRAN (R 4.1.2)
##  P ComplexHeatmap       * 2.10.0   2021-10-26 [?] Bioconductor
##  P cowplot              * 1.1.1    2020-12-30 [?] CRAN (R 4.1.2)
##  P crayon                 1.4.2    2021-10-29 [?] CRAN (R 4.1.2)
##  P curl                   5.2.1    2024-03-01 [?] CRAN (R 4.1.2)
##  P data.table             1.14.2   2021-09-27 [?] CRAN (R 4.1.2)
##  P DBI                    1.1.2    2021-12-20 [?] CRAN (R 4.1.2)
##  P dbplyr                 2.1.1    2021-04-06 [?] CRAN (R 4.1.2)
##  P DelayedArray           0.20.0   2021-10-26 [?] Bioconductor
##  P deldir                 1.0-6    2021-10-23 [?] CRAN (R 4.1.2)
##  P devtools               2.4.5    2022-10-11 [?] CRAN (R 4.1.2)
##  P digest                 0.6.35   2024-03-11 [?] CRAN (R 4.1.2)
##  P docopt                 0.7.1    2020-06-24 [?] CRAN (R 4.1.2)
##  P doParallel             1.0.16   2020-10-16 [?] CRAN (R 4.1.2)
##  P dplyr                * 1.1.1    2023-03-22 [?] CRAN (R 4.1.2)
##  P ellipsis               0.3.2    2021-04-29 [?] CRAN (R 4.1.2)
##  P evaluate               0.23     2023-11-01 [?] CRAN (R 4.1.2)
##  P fansi                  1.0.2    2022-01-14 [?] CRAN (R 4.1.2)
##  P farver                 2.1.0    2021-02-28 [?] CRAN (R 4.1.2)
##  P fastmap                1.1.0    2021-01-25 [?] CRAN (R 4.1.2)
##  P fastmatch              1.1-3    2021-07-23 [?] CRAN (R 4.1.2)
##  P filelock               1.0.2    2018-10-05 [?] CRAN (R 4.1.2)
##  P fitdistrplus           1.1-6    2021-09-28 [?] CRAN (R 4.1.2)
##  P foreach                1.5.1    2020-10-15 [?] CRAN (R 4.1.2)
##  P formatR                1.11     2021-06-01 [?] CRAN (R 4.1.2)
##  P fs                     1.5.2    2021-12-08 [?] CRAN (R 4.1.2)
##  P futile.logger        * 1.4.3    2016-07-10 [?] CRAN (R 4.1.2)
##  P futile.options         1.0.1    2018-04-20 [?] CRAN (R 4.1.2)
##  P future                 1.25.0   2022-04-24 [?] CRAN (R 4.1.2)
##  P future.apply           1.8.1    2021-08-10 [?] CRAN (R 4.1.2)
##  P generics               0.1.3    2022-07-05 [?] CRAN (R 4.1.2)
##  P GenomeInfoDb         * 1.30.1   2022-01-30 [?] Bioconductor
##  P GenomeInfoDbData       1.2.4    2023-11-28 [?] Bioconductor
##  P GenomicAlignments      1.30.0   2021-10-26 [?] Bioconductor
##  P GenomicFeatures      * 1.46.4   2022-01-20 [?] Bioconductor
##  P GenomicRanges        * 1.46.1   2021-11-18 [?] Bioconductor
##  P GetoptLong             1.0.5    2020-12-15 [?] CRAN (R 4.1.2)
##  P ggforce                0.3.3    2021-03-05 [?] CRAN (R 4.1.2)
##  P ggplot2              * 3.4.2    2023-04-03 [?] CRAN (R 4.1.2)
##  P ggrepel                0.9.1    2021-01-15 [?] CRAN (R 4.1.2)
##  P ggridges               0.5.3    2021-01-08 [?] CRAN (R 4.1.2)
##  P ggseqlogo              0.1      2017-07-25 [?] CRAN (R 4.1.2)
##  P git2r                  0.29.0   2021-11-22 [?] CRAN (R 4.1.2)
##  P GlobalOptions          0.1.2    2020-06-10 [?] CRAN (R 4.1.2)
##  P globals                0.14.0   2020-11-22 [?] CRAN (R 4.1.2)
##  P glue                 * 1.6.2    2022-02-24 [?] CRAN (R 4.1.2)
##  P goftest                1.2-3    2021-10-07 [?] CRAN (R 4.1.2)
##  P gprofiler2           * 0.2.3    2024-02-23 [?] CRAN (R 4.1.2)
##  P gridExtra              2.3      2017-09-09 [?] CRAN (R 4.1.2)
##  P gtable                 0.3.0    2019-03-25 [?] CRAN (R 4.1.2)
##  P here                 * 1.0.1    2020-12-13 [?] CRAN (R 4.1.2)
##  P highr                  0.9      2021-04-16 [?] CRAN (R 4.1.2)
##  P hms                    1.1.1    2021-09-26 [?] CRAN (R 4.1.2)
##  P htmltools              0.5.2    2021-08-25 [?] CRAN (R 4.1.2)
##  P htmlwidgets            1.5.4    2021-09-08 [?] CRAN (R 4.1.2)
##  P httpuv                 1.6.5    2022-01-05 [?] CRAN (R 4.1.2)
##  P httr                   1.4.2    2020-07-20 [?] CRAN (R 4.1.2)
##  P ica                    1.0-2    2018-05-24 [?] CRAN (R 4.1.2)
##  P igraph                 2.0.3    2024-03-13 [?] CRAN (R 4.1.2)
##  P IRanges              * 2.28.0   2021-10-26 [?] Bioconductor
##  P irlba                  2.3.5    2021-12-06 [?] CRAN (R 4.1.2)
##  P iterators              1.0.13   2020-10-15 [?] CRAN (R 4.1.2)
##  P jquerylib              0.1.4    2021-04-26 [?] CRAN (R 4.1.2)
##  P jsonlite               1.8.8    2023-12-04 [?] CRAN (R 4.1.2)
##  P KEGGREST               1.34.0   2021-10-26 [?] Bioconductor
##  P KernSmooth             2.23-20  2021-05-03 [?] CRAN (R 4.1.2)
##  P knitr                  1.37     2021-12-16 [?] CRAN (R 4.1.2)
##  P labeling               0.4.2    2020-10-20 [?] CRAN (R 4.1.2)
##  P lambda.r               1.2.4    2019-09-18 [?] CRAN (R 4.1.2)
##  P later                  1.3.0    2021-08-18 [?] CRAN (R 4.1.2)
##  P lattice                0.20-45  2021-09-22 [?] CRAN (R 4.1.2)
##  P lazyeval               0.2.2    2019-03-15 [?] CRAN (R 4.1.2)
##  P leiden                 0.3.9    2021-07-27 [?] CRAN (R 4.1.2)
##  P lifecycle              1.0.3    2022-10-07 [?] CRAN (R 4.1.2)
##  P listenv                0.8.0    2019-12-05 [?] CRAN (R 4.1.2)
##  P lmtest                 0.9-39   2021-11-07 [?] CRAN (R 4.1.2)
##  P lsa                    0.73.2   2020-05-04 [?] CRAN (R 4.1.2)
##  P magrittr             * 2.0.3    2022-03-30 [?] CRAN (R 4.1.2)
##  P MASS                   7.3-54   2021-05-03 [?] CRAN (R 4.1.2)
##  P Matrix                 1.3-4    2021-06-01 [?] CRAN (R 4.1.2)
##  P MatrixGenerics         1.6.0    2021-10-26 [?] Bioconductor
##  P matrixStats            0.61.0   2021-09-17 [?] CRAN (R 4.1.2)
##  P memoise                2.0.1    2021-11-26 [?] CRAN (R 4.1.2)
##  P mgcv                   1.8-38   2021-10-06 [?] CRAN (R 4.1.2)
##  P mime                   0.12     2021-09-28 [?] CRAN (R 4.1.2)
##  P miniUI                 0.1.1.1  2018-05-18 [?] CRAN (R 4.1.2)
##  P munsell                0.5.0    2018-06-12 [?] CRAN (R 4.1.2)
##  P nlme                   3.1-153  2021-09-07 [?] CRAN (R 4.1.2)
##  P parallelly             1.30.0   2021-12-17 [?] CRAN (R 4.1.2)
##  P patchwork              1.1.1    2020-12-17 [?] CRAN (R 4.1.2)
##  P pbapply                1.5-0    2021-09-16 [?] CRAN (R 4.1.2)
##  P pillar                 1.9.0    2023-03-22 [?] RSPM (R 4.1.2)
##  P pkgbuild               1.4.2    2023-06-26 [?] CRAN (R 4.1.2)
##  P pkgconfig              2.0.3    2019-09-22 [?] CRAN (R 4.1.2)
##  P pkgload                1.3.3    2023-09-22 [?] CRAN (R 4.1.2)
##  P plotly                 4.10.0   2021-10-09 [?] CRAN (R 4.1.2)
##  P plyr                   1.8.6    2020-03-03 [?] CRAN (R 4.1.2)
##  P png                    0.1-7    2013-12-03 [?] CRAN (R 4.1.2)
##  P polyclip               1.10-0   2019-03-14 [?] CRAN (R 4.1.2)
##  P prettyunits            1.1.1    2020-01-24 [?] CRAN (R 4.1.2)
##  P processx               3.8.4    2024-03-16 [?] RSPM
##  P profvis                0.3.8    2023-05-02 [?] CRAN (R 4.1.2)
##  P progress               1.2.2    2019-05-16 [?] CRAN (R 4.1.2)
##  P promises               1.2.0.1  2021-02-11 [?] CRAN (R 4.1.2)
##  P ps                     1.7.6    2024-01-18 [?] RSPM
##  P purrr                * 1.0.1    2023-01-10 [?] CRAN (R 4.1.2)
##  P qlcMatrix              0.9.7    2018-04-20 [?] CRAN (R 4.1.2)
##  P R6                     2.5.1    2021-08-19 [?] CRAN (R 4.1.2)
##  P RANN                   2.6.1    2019-01-08 [?] CRAN (R 4.1.2)
##  P rappdirs               0.3.3    2021-01-31 [?] CRAN (R 4.1.2)
##  P RColorBrewer         * 1.1-2    2014-12-07 [?] CRAN (R 4.1.2)
##  P Rcpp                   1.0.8    2022-01-13 [?] CRAN (R 4.1.2)
##  P RcppAnnoy              0.0.19   2021-07-30 [?] CRAN (R 4.1.2)
##  P RcppRoll               0.3.0    2018-06-05 [?] CRAN (R 4.1.2)
##  P RCurl                  1.98-1.5 2021-09-17 [?] CRAN (R 4.1.2)
##  P readr                * 2.1.1    2021-11-30 [?] CRAN (R 4.1.2)
##  P readxl               * 1.3.1    2019-03-13 [?] CRAN (R 4.1.2)
##  P remotes                2.4.2.1  2023-07-18 [?] CRAN (R 4.1.2)
##  P renv                   1.0.3    2023-09-19 [?] CRAN (R 4.1.2)
##  P reshape2               1.4.4    2020-04-09 [?] CRAN (R 4.1.2)
##  P restfulr               0.0.13   2017-08-06 [?] CRAN (R 4.1.2)
##  P reticulate             1.23     2022-01-14 [?] CRAN (R 4.1.2)
##  P rjson                  0.2.21   2022-01-09 [?] CRAN (R 4.1.2)
##  P rlang                  1.1.3    2024-01-10 [?] CRAN (R 4.1.2)
##  P rmarkdown              2.11     2021-09-14 [?] CRAN (R 4.1.2)
##  P ROCR                   1.0-11   2020-05-02 [?] CRAN (R 4.1.2)
##  P rpart                  4.1-15   2019-04-12 [?] CRAN (R 4.1.2)
##  P rprojroot              2.0.2    2020-11-15 [?] CRAN (R 4.1.2)
##  P Rsamtools              2.10.0   2021-10-26 [?] Bioconductor
##  P RSQLite                2.2.9    2021-12-06 [?] CRAN (R 4.1.2)
##  P rtracklayer          * 1.54.0   2021-10-26 [?] Bioconductor
##  P Rtsne                  0.15     2018-11-10 [?] CRAN (R 4.1.2)
##  P S4Vectors            * 0.32.4   2022-03-24 [?] Bioconductor
##  P sass                   0.4.0    2021-05-12 [?] CRAN (R 4.1.2)
##  P scales                 1.2.1    2022-08-20 [?] CRAN (R 4.1.2)
##  P scattermore            0.7      2020-11-24 [?] CRAN (R 4.1.2)
##  P sctransform            0.3.3    2022-01-13 [?] CRAN (R 4.1.2)
##  P sessioninfo            1.2.2    2021-12-06 [?] CRAN (R 4.1.2)
##  P Seurat               * 4.0.0    2021-01-30 [?] CRAN (R 4.1.2)
##  P SeuratObject         * 4.0.4    2021-11-23 [?] CRAN (R 4.1.2)
##  P shape                  1.4.6    2021-05-19 [?] CRAN (R 4.1.2)
##  P shiny                  1.7.1    2021-10-02 [?] CRAN (R 4.1.2)
##  P Signac               * 1.3.0    2021-07-12 [?] CRAN (R 4.1.2)
##  P slam                   0.1-50   2022-01-08 [?] CRAN (R 4.1.2)
##  P SnowballC              0.7.0    2020-04-01 [?] CRAN (R 4.1.2)
##  P sparsesvd              0.2      2019-07-15 [?] CRAN (R 4.1.2)
##  P spatstat               1.64-1   2020-05-12 [?] CRAN (R 4.1.2)
##  P spatstat.data          2.1-2    2021-12-17 [?] CRAN (R 4.1.2)
##  P spatstat.utils         2.3-0    2021-12-12 [?] CRAN (R 4.1.2)
##  P stringi                1.7.6    2021-11-29 [?] CRAN (R 4.1.2)
##  P stringr              * 1.5.0    2022-12-02 [?] CRAN (R 4.1.2)
##  P SummarizedExperiment   1.24.0   2021-10-26 [?] Bioconductor
##  P survival               3.2-13   2021-08-24 [?] CRAN (R 4.1.2)
##  P tensor                 1.5      2012-05-05 [?] CRAN (R 4.1.2)
##  P tibble                 3.2.1    2023-03-20 [?] RSPM (R 4.1.2)
##  P tidyr                * 1.3.0    2023-01-24 [?] CRAN (R 4.1.2)
##  P tidyselect             1.2.0    2022-10-10 [?] CRAN (R 4.1.2)
##  P tweenr                 1.0.2    2021-03-23 [?] CRAN (R 4.1.2)
##  P tzdb                   0.3.0    2022-03-28 [?] CRAN (R 4.1.2)
##  P urlchecker             1.0.1    2021-11-30 [?] CRAN (R 4.1.2)
##  P usethis                2.2.2    2023-07-06 [?] CRAN (R 4.1.2)
##  P utf8                   1.2.2    2021-07-24 [?] CRAN (R 4.1.2)
##  P uwot                   0.1.11   2021-12-02 [?] CRAN (R 4.1.2)
##  P vctrs                  0.6.5    2023-12-01 [?] CRAN (R 4.1.2)
##  P VennDiagram          * 1.7.3    2022-04-12 [?] CRAN (R 4.1.2)
##  P viridis              * 0.5.1    2018-03-29 [?] RSPM (R 4.1.2)
##  P viridisLite          * 0.3.0    2018-02-01 [?] CRAN (R 4.1.2)
##  P vroom                  1.5.7    2021-11-30 [?] CRAN (R 4.1.2)
##  P withr                  2.5.0    2022-03-03 [?] CRAN (R 4.1.2)
##  P xfun                   0.29     2021-12-14 [?] CRAN (R 4.1.2)
##  P XML                    3.99-0.8 2021-09-17 [?] CRAN (R 4.1.2)
##  P xml2                   1.3.3    2021-11-30 [?] CRAN (R 4.1.2)
##  P xtable                 1.8-4    2019-04-21 [?] CRAN (R 4.1.2)
##  P XVector                0.34.0   2021-10-26 [?] Bioconductor
##  P yaml                   2.2.1    2020-02-01 [?] CRAN (R 4.1.2)
##  P zlibbioc               1.40.0   2021-10-26 [?] Bioconductor
##  P zoo                    1.8-9    2021-03-09 [?] CRAN (R 4.1.2)
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
