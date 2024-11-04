---
title: "08 - Extracranial neuroblastoma"
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
## Document index: 08
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
## public/output/08
```

```
## public/figures/08//
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

This document compares NB-FOXR2 with extracranial neuroblastoma samples
stratified by stage and FOXR2 expression.

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
library(ggrepel)
library(tibble)
library(cowplot)
library(Seurat)
library(data.table)
library(MetBrewer)
library(dendextend)
library(readxl)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(ComplexHeatmap)

source(here("include/style.R"))
source(here("code/functions/RNAseq.R"))
source(here("code/functions/scRNAseq.R"))
source(here("code/functions/ssGSEA.R"))

ggplot2::theme_set(theme_min())

square_theme <- theme(aspect.ratio = 1)
large_text <- theme(text = element_text(size = 15))

conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::first)
conflicted::conflicts_prefer(dplyr::last)
conflicted::conflicts_prefer(dplyr::combine)
conflicted::conflicts_prefer(dplyr::intersect)
conflicted::conflicts_prefer(base::setdiff)
conflicted::conflicts_prefer(clusterProfiler::slice)
conflicted::conflicts_prefer(dplyr::rename)
conflicted::conflicts_prefer(dplyr::desc)
conflicted::conflicts_prefer(dplyr::filter)
```

# Load EC-NB samples

## Gartlgruber et al. 2021 (n=579)

Publication: [Gartlgruber et al. *Nature Cancer* 2021](https://pubmed.ncbi.nlm.nih.gov/35121888/)

### Load data

Load raw counts and metadata. This data was downloaded from a 
[public web app]((https://nbseb087.dkfz.de/project_NB_SE_viz/#shiny-tab-faq)) 
provided by the authors of the study.

Data was provided as un-normalized ("raw") bulk counts and per-tumor metadata.

Load DESeq normalized counts produced in in R markdown document 03, for plotting
gene expression.


```r
# Load files and remove unused ones to save space
load(here("output/03/Gartlgruber_et_al_counts.Rda"))
rm(ecnb_counts_norm)
rm(ecnb_counts_vst)
```


### Add FOXR2 positivity to metadata

Here, we label each tumor as FOXR2 positive or negative in per-tumor metadata 
based on normalized expression of FOXR2. 
We consider tumors with FOXR2 DESeq normalized expression > 2 to be FOXR2 positive. 


```r
# Set threshold of FOXR2+ expression used throughout this document
threshold <- 2

ecnb_foxr2_counts <- ecnb_counts_tidy %>%
    filter(gene_symbol == "FOXR2")

foxr2_pos_samples <- ecnb_foxr2_counts %>% 
    filter(gene_expression > threshold) %>% 
    .$sample

ecnb_foxr2_counts <- ecnb_foxr2_counts %>% mutate(FOXR2_positive = case_when(sample %in% foxr2_pos_samples ~ "Y",
                                      T ~ "N"))

# Print the number of FOXR2 positive and negative samples
ecnb_foxr2_counts %>% count(FOXR2_positive, sort = F)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["FOXR2_positive"],"name":[1],"type":["chr"],"align":["left"]},{"label":["n"],"name":[2],"type":["int"],"align":["right"]}],"data":[{"1":"N","2":"473"},{"1":"Y","2":"106"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# Print the number of FOXR2 positive and negative samples in each tumor stage
ecnb_foxr2_counts %>% count(Stage, FOXR2_positive, sort = F)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Stage"],"name":[1],"type":["fct"],"align":["left"]},{"label":["FOXR2_positive"],"name":[2],"type":["chr"],"align":["left"]},{"label":["n"],"name":[3],"type":["int"],"align":["right"]}],"data":[{"1":"1-3;4S","2":"N","3":"292"},{"1":"1-3;4S","2":"Y","3":"41"},{"1":"4","2":"N","3":"171"},{"1":"4","2":"Y","3":"58"},{"1":"NA","2":"N","3":"10"},{"1":"NA","2":"Y","3":"7"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
ecnb_counts_tidy <- ecnb_counts_tidy %>% 
    mutate(FOXR2_positive = case_when(sample %in% foxr2_pos_samples ~ "FOXR2+",
                                      T ~ "FOXR2-"))
```

## TARGET cohort (n=128)

Downloaded and processed by Steven Hébert in April 2024. 

They were downloaded as count matrices from the public Genomic
686 Data Commons (GDC) cancer portal (https://portal.gdc.cancer.gov/) 
using the Cohort Builder, with options: 

* `project = TARGET-346 NBL`
* `Experimental Strategy = RNA-seq`
* `Access = open`

### Load data

Patients were subsetted to Risk = HR (high risk) and Stage = 3 or 4 (not 4S).
This resulted in a total of 128 patient samples.

Load patient sample metadata.


```r
# Load, replace the empty value '-- with NA values and remove all-NA columns, 
# and rename id column to Sample
target_meta <- data.table::fread(here("data/RNAseq/external_data/TARGET_ECNB/clinical.tsv")) %>%
    apply(2, 
              function(col) { 
                  gsub(col, pattern = "'--", replacement = as.character(NA)) 
              }
          ) %>% 
    as.data.frame() %>% 
    select_if(~ !all(is.na(.))) %>% 
    dplyr::rename(Sample = case_submitter_id)
```

Load DESeq-normalized counts (for plotting gene expression).


```r
# load normalized counts
target_norm <- here("data/RNAseq/external_data/TARGET_ECNB/counts_ensembl_common/GENCODEv36.norm.tsv.gz") %>%
    read.table(header = T, row.names = 1, check.names = F, sep = "\t")

# Display first 5 row/col to check
target_norm[1:5, 1:5]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["TARGET-30-PAIFXV-01A"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["TARGET-30-PAIPGU-01A"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["TARGET-30-PAISNS-01A"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["TARGET-30-PAITCI-01A"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["TARGET-30-PAITEG-01A"],"name":[5],"type":["dbl"],"align":["right"]}],"data":[{"1":"1612.983582","2":"1337.289410","3":"3232.9537","4":"1649.4543","5":"1697.629906","_rn_":"ENSG00000000003:TSPAN6"},{"1":"9.170469","2":"8.985147","3":"0.0000","4":"0.0000","5":"0.607381","_rn_":"ENSG00000000005:TNMD"},{"1":"1610.945700","2":"2105.519496","3":"2272.0982","4":"4908.1899","5":"2858.335005","_rn_":"ENSG00000000419:DPM1"},{"1":"543.095546","2":"489.690523","3":"816.3660","4":"655.7644","5":"368.680269","_rn_":"ENSG00000000457:SCYL3"},{"1":"177.295732","2":"371.386085","3":"686.3254","4":"309.5681","5":"126.335249","_rn_":"ENSG00000000460:C1orf112"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

### Add FOXR2 positivity to metadata

Based on threshold = 2.


```r
# DESeq norm expression threshold used for FOXR2 +/-
threshold <- 2

target_counts_tidy <- target_norm %>% 
    as.data.frame() %>% 
    rownames_to_column("Gene") %>% 
    separate(col = "Gene", sep = ":", into = c("ENS", "gene_symbol")) %>% 
    pivot_longer(cols = -c("ENS", "gene_symbol"), names_to = "ID", values_to = "gene_expression")
    
foxr2_pos_samples_target <- target_counts_tidy %>% 
    filter(gene_symbol == "FOXR2") %>% 
    filter(gene_expression > threshold) %>% 
    .$ID 

target_meta_foxr2 <- target_counts_tidy %>% 
    mutate(FOXR2_positive = case_when(ID %in% foxr2_pos_samples_target ~ "Y",
                                      T ~ "N")) %>% 
    select(ID, FOXR2_positive) %>% distinct
```

# FOXR2 coexpression with MGE TFs in EC-NB

## Gartlgruber et al. 2021

Stage 4, high-risk tumors only
FOXR2 threshold: DESEq normalized value "2".


```r
tum_interest <- ecnb_counts_tidy %>% 
    filter(Risk == "HR") %>%
    filter(Stage == "4") %>% 
    filter(gene_symbol %in% c("FOXR2", "LHX6", "DLX5", "DLX6"))

# DLX5
foxr2_dlx5 <- tum_interest %>% 
    filter(gene_symbol %in% c("FOXR2", "DLX5")) %>% 
    pivot_wider(names_from = "gene_symbol", 
                values_from = "gene_expression", 
                id_cols = c("sample", "FOXR2_positive")) 

(cor_foxr2_dlx5 <- cor(foxr2_dlx5$FOXR2, foxr2_dlx5$DLX5, method = "pearson"))
```

```
## [1] 0.1477744
```

```r
(sig_dlx5 <- wilcox.test(foxr2_dlx5 %>% filter(FOXR2_positive == "FOXR2+") %>% .$DLX5,
                         foxr2_dlx5 %>% filter(FOXR2_positive == "FOXR2-") %>% .$DLX5))
```

```
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  foxr2_dlx5 %>% filter(FOXR2_positive == "FOXR2+") %>% .$DLX5 and foxr2_dlx5 %>% filter(FOXR2_positive == "FOXR2-") %>% .$DLX5
## W = 3549, p-value = 0.8924
## alternative hypothesis: true location shift is not equal to 0
```

```r
dlx5_scatter <- foxr2_dlx5 %>% 
    ggplot(aes(x = FOXR2, y = DLX5, alpha = 0.2)) + 
    geom_point(alpha = 0.5, size = 3) + 
    ggtitle("DLX5 - EC-NB stage 4") +
    annotate("text", label = glue("r == {round(cor_foxr2_dlx5,4)}"),
             x = 470, y = 4000, parse = T)

dlx5_vln <- foxr2_dlx5 %>% 
    ggplot(aes(x = FOXR2_positive, y = DLX5, fill = FOXR2_positive)) + 
    geom_violin(scale = "width") + 
    ggtitle("DLX5 - EC-NB stage 4") + 
    no_legend() + scale_y_log10() + 
    stat_summary(fun.y=median, geom="crossbar", size=1, 
                 color="black", aes(width = 0.3)) +
    annotate("text", label = glue("p == {round(sig_dlx5$p.value,4)}"),
             x = 2, y = 4000, parse = T)


# DLX6
foxr2_dlx6 <- tum_interest %>% 
    filter(gene_symbol %in% c("FOXR2", "DLX6")) %>% 
    pivot_wider(names_from = "gene_symbol", 
                values_from = "gene_expression", 
                id_cols = c("sample", "FOXR2_positive")) 

(cor_foxr2_dlx6 <- cor(foxr2_dlx6$FOXR2, foxr2_dlx6$DLX6, method = "pearson"))
```

```
## [1] 0.2062904
```

```r
(sig_dlx6 <- wilcox.test(foxr2_dlx6 %>% filter(FOXR2_positive == "FOXR2+") %>% .$DLX6,
                         foxr2_dlx6 %>% filter(FOXR2_positive == "FOXR2-") %>% .$DLX6))
```

```
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  foxr2_dlx6 %>% filter(FOXR2_positive == "FOXR2+") %>% .$DLX6 and foxr2_dlx6 %>% filter(FOXR2_positive == "FOXR2-") %>% .$DLX6
## W = 3461.5, p-value = 0.6946
## alternative hypothesis: true location shift is not equal to 0
```

```r
dlx6_scatter <- foxr2_dlx6 %>% 
    ggplot(aes(x = FOXR2, y = DLX6, alpha = 0.2)) + 
    geom_point(alpha = 0.5, size = 3) + ggtitle("DLX6 - EC-NB stage 4") +
    annotate("text", label = glue("r == {round(cor_foxr2_dlx6,4)}"),
             x = 470, y = 3000, parse = T)

dlx6_vln <- foxr2_dlx6 %>% 
    ggplot(aes(x = FOXR2_positive, y = DLX6, fill = FOXR2_positive)) + 
    geom_violin(scale = "width") + 
    ggtitle("DLX6 - EC-NB stage 4") + 
    no_legend() + scale_y_log10() + 
    stat_summary(fun.y=median, geom="crossbar", size=1, 
                 color="black", aes(width = 0.3))  +
    annotate("text", label = glue("p == {round(sig_dlx6$p.value,4)}"),
             x = 2, y = 4000, parse = T)



# LHX6
foxr2_lhx6 <- tum_interest %>% 
    filter(gene_symbol %in% c("FOXR2", "LHX6")) %>% 
    pivot_wider(names_from = "gene_symbol", 
                values_from = "gene_expression", 
                id_cols = c("sample", "FOXR2_positive")) 

(cor_foxr2_lhx6 <- cor(foxr2_lhx6$FOXR2, foxr2_lhx6$LHX6, method = "pearson"))
```

```
## [1] -0.04012718
```

```r
(sig_lhx6 <- wilcox.test(foxr2_lhx6 %>% filter(FOXR2_positive == "FOXR2+") %>% .$LHX6,
                         foxr2_lhx6 %>% filter(FOXR2_positive == "FOXR2-") %>% .$LHX6))
```

```
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  foxr2_lhx6 %>% filter(FOXR2_positive == "FOXR2+") %>% .$LHX6 and foxr2_lhx6 %>% filter(FOXR2_positive == "FOXR2-") %>% .$LHX6
## W = 3546, p-value = 0.8854
## alternative hypothesis: true location shift is not equal to 0
```

```r
lhx6_scatter <- foxr2_lhx6 %>% 
    ggplot(aes(x = FOXR2, y = LHX6, alpha = 0.2)) + 
    geom_point(alpha = 0.5, size = 3) + 
    ggtitle("LHX6 - EC-NB stage 4") +
    annotate("text", label = glue("r == {round(cor_foxr2_lhx6,4)}"),
             x = 470, y = 600, parse = T)

lhx6_vln <- foxr2_lhx6 %>% 
    ggplot(aes(x = FOXR2_positive, y = LHX6, fill = FOXR2_positive)) + 
    geom_violin(scale = "width") + 
    ggtitle("LHX6 - EC-NB stage 4") + 
    no_legend() + 
    stat_summary(fun.y=median, geom="crossbar", size=1, 
                 color="black", aes(width = 0.3))  +
    annotate("text", label = glue("p == {round(sig_lhx6$p.value,4)}"),
             x = 2, y = 600, parse = T)


cowplot::plot_grid(dlx5_scatter, dlx5_vln,
                                dlx6_scatter, dlx6_vln,
                                lhx6_scatter, lhx6_vln,
                   nrow = 1, align = "hv")
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/08//FOXR2-coexp-ECNB-1.png)<!-- -->

## TARGET cohort

FOXR2 + threshold: DESEq normalized value "2".


```r
tum_interest <- target_counts_tidy %>% 
    filter(gene_symbol %in% c("FOXR2", "LHX6", "DLX5", "DLX6")) %>%
    left_join(., target_meta_foxr2, by = "ID") %>% 
    dplyr::rename("sample" = "ID") %>% 
    mutate(FOXR2_positive = case_when(FOXR2_positive == "Y" ~ "FOXR2+",
                                      FOXR2_positive == "N" ~ "FOXR2-",
                                      TRUE ~ NA))

# DLX5
foxr2_dlx5 <- tum_interest %>% 
    filter(gene_symbol %in% c("FOXR2", "DLX5")) %>% 
    pivot_wider(names_from = "gene_symbol", 
                values_from = "gene_expression", 
                id_cols = c("sample", "FOXR2_positive")) 

(cor_foxr2_dlx5 <- cor(foxr2_dlx5$FOXR2, foxr2_dlx5$DLX5, method = "pearson"))
```

```
## [1] 0.07350015
```

```r
(sig_dlx5 <- wilcox.test(foxr2_dlx5 %>% filter(FOXR2_positive == "FOXR2+") %>% .$DLX5,
                         foxr2_dlx5 %>% filter(FOXR2_positive == "FOXR2-") %>% .$DLX5))
```

```
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  foxr2_dlx5 %>% filter(FOXR2_positive == "FOXR2+") %>% .$DLX5 and foxr2_dlx5 %>% filter(FOXR2_positive == "FOXR2-") %>% .$DLX5
## W = 1816, p-value = 0.2406
## alternative hypothesis: true location shift is not equal to 0
```

```r
dlx5_scatter <- foxr2_dlx5 %>% 
    ggplot(aes(x = FOXR2, y = DLX5, alpha = 0.2)) + 
    geom_point(alpha = 0.5, size = 3) + 
    ggtitle("DLX5 - EC-NB stage 4") +
    annotate("text", label = glue("r == {round(cor_foxr2_dlx5,4)}"),
            x = 900, y = 5500, parse = T)

dlx5_vln <- foxr2_dlx5 %>% 
    ggplot(aes(x = FOXR2_positive, y = DLX5, fill = FOXR2_positive)) + 
    geom_violin(scale = "width") + 
    ggtitle("DLX5 - EC-NB stage 4") + 
    no_legend() + scale_y_log10() + 
    stat_summary(fun.y=median, geom="crossbar", size=1, 
                 color="black", aes(width = 0.3)) +
    annotate("text", label = glue("p == {round(sig_dlx5$p.value,4)}"),
             x = 2, y = 5500, parse = T)


# DLX6
foxr2_dlx6 <- tum_interest %>% 
    filter(gene_symbol %in% c("FOXR2", "DLX6")) %>% 
    pivot_wider(names_from = "gene_symbol", 
                values_from = "gene_expression", 
                id_cols = c("sample", "FOXR2_positive")) 

(cor_foxr2_dlx6 <- cor(foxr2_dlx6$FOXR2, foxr2_dlx6$DLX6, method = "pearson"))
```

```
## [1] 0.1009971
```

```r
(sig_dlx6 <- wilcox.test(foxr2_dlx6 %>% filter(FOXR2_positive == "FOXR2+") %>% .$DLX6,
                         foxr2_dlx6 %>% filter(FOXR2_positive == "FOXR2-") %>% .$DLX6))
```

```
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  foxr2_dlx6 %>% filter(FOXR2_positive == "FOXR2+") %>% .$DLX6 and foxr2_dlx6 %>% filter(FOXR2_positive == "FOXR2-") %>% .$DLX6
## W = 1823.5, p-value = 0.2248
## alternative hypothesis: true location shift is not equal to 0
```

```r
dlx6_scatter <- foxr2_dlx6 %>% 
    ggplot(aes(x = FOXR2, y = DLX6, alpha = 0.2)) + 
    geom_point(alpha = 0.5, size = 3) + ggtitle("DLX6 - EC-NB stage 4") +
    annotate("text", label = glue("r == {round(cor_foxr2_dlx6,4)}"),
              x = 900, y = 3500, parse = T)

dlx6_vln <- foxr2_dlx6 %>% 
    ggplot(aes(x = FOXR2_positive, y = DLX6, fill = FOXR2_positive)) + 
    geom_violin(scale = "width") + 
    ggtitle("DLX6 - EC-NB stage 4") + 
    no_legend() + scale_y_log10() + 
    stat_summary(fun.y=median, geom="crossbar", size=1, 
                 color="black", aes(width = 0.3))  +
    annotate("text", label = glue("p == {round(sig_dlx6$p.value,4)}"),
              x = 2, y = 3500, parse = T)



# LHX6
foxr2_lhx6 <- tum_interest %>% 
    filter(gene_symbol %in% c("FOXR2", "LHX6")) %>% 
    pivot_wider(names_from = "gene_symbol", 
                values_from = "gene_expression", 
                id_cols = c("sample", "FOXR2_positive")) 

(cor_foxr2_lhx6 <- cor(foxr2_lhx6$FOXR2, foxr2_lhx6$LHX6, method = "pearson"))
```

```
## [1] -0.1011051
```

```r
(sig_lhx6 <- wilcox.test(foxr2_lhx6 %>% filter(FOXR2_positive == "FOXR2+") %>% .$LHX6,
                         foxr2_lhx6 %>% filter(FOXR2_positive == "FOXR2-") %>% .$LHX6))
```

```
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  foxr2_lhx6 %>% filter(FOXR2_positive == "FOXR2+") %>% .$LHX6 and foxr2_lhx6 %>% filter(FOXR2_positive == "FOXR2-") %>% .$LHX6
## W = 1493, p-value = 0.5729
## alternative hypothesis: true location shift is not equal to 0
```

```r
lhx6_scatter <- foxr2_lhx6 %>% 
    ggplot(aes(x = FOXR2, y = LHX6, alpha = 0.2)) + 
    geom_point(alpha = 0.5, size = 3) + 
    ggtitle("LHX6 - EC-NB stage 4") +
    annotate("text", label = glue("r == {round(cor_foxr2_lhx6,4)}"),
             x = 900, y = 900, parse = T)

lhx6_vln <- foxr2_lhx6 %>% 
    ggplot(aes(x = FOXR2_positive, y = LHX6, fill = FOXR2_positive)) + 
    geom_violin(scale = "width") + 
    ggtitle("LHX6 - EC-NB stage 4") + 
    no_legend() + 
    stat_summary(fun.y=median, geom="crossbar", size=1, 
                 color="black", aes(width = 0.3))  +
    annotate("text", label = glue("p == {round(sig_lhx6$p.value,4)}"),
             x = 2, y = 900, parse = T)


cowplot::plot_grid(dlx5_scatter, dlx5_vln,
                    dlx6_scatter, dlx6_vln,
                    lhx6_scatter, lhx6_vln,
                   nrow = 1, align = "hv")
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/08//FOXR2-coexp-TARGET-1.png)<!-- -->


<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2024-11-01 14:16:53
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
##  P abind                  1.4-5    2016-07-21 [?] CRAN (R 4.1.2)
##  P annotate               1.72.0   2021-10-26 [?] Bioconductor
##  P AnnotationDbi        * 1.56.2   2021-11-09 [?] Bioconductor
##  P ape                    5.7-1    2023-03-13 [?] RSPM (R 4.1.2)
##  P aplot                  0.2.2    2023-10-06 [?] CRAN (R 4.1.2)
##  P beachmat               2.10.0   2021-10-26 [?] Bioconductor
##  P Biobase              * 2.54.0   2021-10-26 [?] Bioconductor
##  P BiocGenerics         * 0.40.0   2021-10-26 [?] Bioconductor
##  P BiocManager            1.30.15  2021-05-11 [?] CRAN (R 4.1.2)
##  P BiocParallel           1.28.3   2021-12-09 [?] Bioconductor
##  P BiocSingular           1.10.0   2021-10-26 [?] Bioconductor
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
##  P clusterProfiler      * 4.2.2    2022-01-13 [?] Bioconductor
##  P codetools              0.2-18   2020-11-04 [?] CRAN (R 4.1.2)
##  P colorspace             2.0-2    2021-06-24 [?] CRAN (R 4.1.2)
##  P ComplexHeatmap       * 2.10.0   2021-10-26 [?] Bioconductor
##  P conflicted             1.2.0    2023-02-01 [?] CRAN (R 4.1.2)
##  P cowplot              * 1.1.1    2020-12-30 [?] CRAN (R 4.1.2)
##  P crayon                 1.4.2    2021-10-29 [?] CRAN (R 4.1.2)
##  P data.table           * 1.14.2   2021-09-27 [?] CRAN (R 4.1.2)
##  P DBI                    1.1.2    2021-12-20 [?] CRAN (R 4.1.2)
##  P DelayedArray           0.20.0   2021-10-26 [?] Bioconductor
##  P DelayedMatrixStats     1.16.0   2021-10-26 [?] Bioconductor
##  P deldir                 1.0-6    2021-10-23 [?] CRAN (R 4.1.2)
##  P dendextend           * 1.17.1   2023-03-25 [?] CRAN (R 4.1.2)
##  P devtools               2.4.5    2022-10-11 [?] CRAN (R 4.1.2)
##  P digest                 0.6.35   2024-03-11 [?] CRAN (R 4.1.2)
##  P DO.db                  2.9      2024-09-20 [?] Bioconductor
##  P doParallel             1.0.16   2020-10-16 [?] CRAN (R 4.1.2)
##  P DOSE                   3.20.1   2021-11-18 [?] Bioconductor
##  P downloader             0.4      2015-07-09 [?] CRAN (R 4.1.2)
##  P dplyr                * 1.1.1    2023-03-22 [?] CRAN (R 4.1.2)
##  P ellipsis               0.3.2    2021-04-29 [?] CRAN (R 4.1.2)
##  P enrichplot             1.14.2   2022-02-24 [?] Bioconductor
##  P evaluate               0.23     2023-11-01 [?] CRAN (R 4.1.2)
##  P fansi                  1.0.2    2022-01-14 [?] CRAN (R 4.1.2)
##  P farver                 2.1.0    2021-02-28 [?] CRAN (R 4.1.2)
##  P fastmap                1.1.0    2021-01-25 [?] CRAN (R 4.1.2)
##  P fastmatch              1.1-3    2021-07-23 [?] CRAN (R 4.1.2)
##  P fgsea                  1.20.0   2021-10-26 [?] Bioconductor
##  P fitdistrplus           1.1-6    2021-09-28 [?] CRAN (R 4.1.2)
##  P foreach                1.5.1    2020-10-15 [?] CRAN (R 4.1.2)
##  P fs                     1.5.2    2021-12-08 [?] CRAN (R 4.1.2)
##  P future                 1.25.0   2022-04-24 [?] CRAN (R 4.1.2)
##  P future.apply           1.8.1    2021-08-10 [?] CRAN (R 4.1.2)
##  P generics               0.1.3    2022-07-05 [?] CRAN (R 4.1.2)
##  P GenomeInfoDb           1.30.1   2022-01-30 [?] Bioconductor
##  P GenomeInfoDbData       1.2.4    2023-11-28 [?] Bioconductor
##  P GenomicRanges          1.46.1   2021-11-18 [?] Bioconductor
##  P GetoptLong             1.0.5    2020-12-15 [?] CRAN (R 4.1.2)
##  P ggforce                0.3.3    2021-03-05 [?] CRAN (R 4.1.2)
##  P ggfun                  0.1.4    2024-01-19 [?] CRAN (R 4.1.2)
##  P ggplot2              * 3.4.2    2023-04-03 [?] CRAN (R 4.1.2)
##  P ggplotify              0.1.2    2023-08-09 [?] CRAN (R 4.1.2)
##  P ggraph                 2.2.1    2024-03-07 [?] CRAN (R 4.1.2)
##  P ggrepel              * 0.9.1    2021-01-15 [?] CRAN (R 4.1.2)
##  P ggridges               0.5.3    2021-01-08 [?] CRAN (R 4.1.2)
##  P ggtree                 3.13.0   2024-05-19 [?] Github (YuLab-SMU/ggtree@05ef652)
##  P git2r                  0.29.0   2021-11-22 [?] CRAN (R 4.1.2)
##  P GlobalOptions          0.1.2    2020-06-10 [?] CRAN (R 4.1.2)
##  P globals                0.14.0   2020-11-22 [?] CRAN (R 4.1.2)
##  P glue                 * 1.6.2    2022-02-24 [?] CRAN (R 4.1.2)
##  P GO.db                * 3.14.0   2024-09-20 [?] Bioconductor
##  P goftest                1.2-3    2021-10-07 [?] CRAN (R 4.1.2)
##  P GOSemSim               2.20.0   2021-10-26 [?] Bioconductor
##  P graph                  1.72.0   2021-10-26 [?] Bioconductor
##  P graphlayouts           1.1.1    2024-03-09 [?] CRAN (R 4.1.2)
##  P gridExtra              2.3      2017-09-09 [?] CRAN (R 4.1.2)
##  P gridGraphics           0.5-1    2020-12-13 [?] CRAN (R 4.1.2)
##  P GSEABase               1.56.0   2021-10-26 [?] Bioconductor
##  P GSVA                 * 1.42.0   2021-10-26 [?] Bioconductor
##  P gtable                 0.3.0    2019-03-25 [?] CRAN (R 4.1.2)
##  P HDF5Array              1.22.1   2021-11-14 [?] Bioconductor
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
##  P later                  1.3.0    2021-08-18 [?] CRAN (R 4.1.2)
##  P lattice                0.20-45  2021-09-22 [?] CRAN (R 4.1.2)
##  P lazyeval               0.2.2    2019-03-15 [?] CRAN (R 4.1.2)
##  P leiden                 0.3.9    2021-07-27 [?] CRAN (R 4.1.2)
##  P lifecycle              1.0.3    2022-10-07 [?] CRAN (R 4.1.2)
##  P listenv                0.8.0    2019-12-05 [?] CRAN (R 4.1.2)
##  P lmtest                 0.9-39   2021-11-07 [?] CRAN (R 4.1.2)
##  P magrittr             * 2.0.3    2022-03-30 [?] CRAN (R 4.1.2)
##  P MASS                   7.3-54   2021-05-03 [?] CRAN (R 4.1.2)
##  P Matrix                 1.3-4    2021-06-01 [?] CRAN (R 4.1.2)
##  P MatrixGenerics         1.6.0    2021-10-26 [?] Bioconductor
##  P matrixStats            0.61.0   2021-09-17 [?] CRAN (R 4.1.2)
##  P memoise                2.0.1    2021-11-26 [?] CRAN (R 4.1.2)
##  P MetBrewer            * 0.2.0    2022-03-21 [?] CRAN (R 4.1.2)
##  P mgcv                   1.8-38   2021-10-06 [?] CRAN (R 4.1.2)
##  P mime                   0.12     2021-09-28 [?] CRAN (R 4.1.2)
##  P miniUI                 0.1.1.1  2018-05-18 [?] CRAN (R 4.1.2)
##  P munsell                0.5.0    2018-06-12 [?] CRAN (R 4.1.2)
##  P nlme                   3.1-153  2021-09-07 [?] CRAN (R 4.1.2)
##  P org.Hs.eg.db         * 3.14.0   2024-03-31 [?] Bioconductor
##  P parallelly             1.30.0   2021-12-17 [?] CRAN (R 4.1.2)
##  P patchwork              1.1.1    2020-12-17 [?] CRAN (R 4.1.2)
##  P pbapply              * 1.5-0    2021-09-16 [?] CRAN (R 4.1.2)
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
##  P promises               1.2.0.1  2021-02-11 [?] CRAN (R 4.1.2)
##  P ps                     1.7.6    2024-01-18 [?] RSPM
##  P purrr                * 1.0.1    2023-01-10 [?] CRAN (R 4.1.2)
##  P qvalue                 2.26.0   2021-10-26 [?] Bioconductor
##  P R6                     2.5.1    2021-08-19 [?] CRAN (R 4.1.2)
##  P RANN                   2.6.1    2019-01-08 [?] CRAN (R 4.1.2)
##  P RColorBrewer         * 1.1-2    2014-12-07 [?] CRAN (R 4.1.2)
##  P Rcpp                   1.0.8    2022-01-13 [?] CRAN (R 4.1.2)
##  P RcppAnnoy              0.0.19   2021-07-30 [?] CRAN (R 4.1.2)
##  P RCurl                  1.98-1.5 2021-09-17 [?] CRAN (R 4.1.2)
##  P readr                * 2.1.1    2021-11-30 [?] CRAN (R 4.1.2)
##  P readxl               * 1.3.1    2019-03-13 [?] CRAN (R 4.1.2)
##  P remotes                2.4.2.1  2023-07-18 [?] CRAN (R 4.1.2)
##  P renv                   1.0.3    2023-09-19 [?] CRAN (R 4.1.2)
##  P reshape2               1.4.4    2020-04-09 [?] CRAN (R 4.1.2)
##  P reticulate             1.23     2022-01-14 [?] CRAN (R 4.1.2)
##  P rhdf5                  2.38.1   2022-03-10 [?] Bioconductor
##  P rhdf5filters           1.6.0    2021-10-26 [?] Bioconductor
##  P Rhdf5lib               1.16.0   2021-10-26 [?] Bioconductor
##  P rjson                  0.2.21   2022-01-09 [?] CRAN (R 4.1.2)
##  P rlang                  1.1.3    2024-01-10 [?] CRAN (R 4.1.2)
##  P rmarkdown              2.11     2021-09-14 [?] CRAN (R 4.1.2)
##  P ROCR                   1.0-11   2020-05-02 [?] CRAN (R 4.1.2)
##  P rpart                  4.1-15   2019-04-12 [?] CRAN (R 4.1.2)
##  P rprojroot              2.0.2    2020-11-15 [?] CRAN (R 4.1.2)
##  P RSQLite                2.2.9    2021-12-06 [?] CRAN (R 4.1.2)
##  P rsvd                   1.0.5    2021-04-16 [?] RSPM (R 4.1.2)
##  P Rtsne                  0.15     2018-11-10 [?] CRAN (R 4.1.2)
##  P S4Vectors            * 0.32.4   2022-03-24 [?] Bioconductor
##  P sass                   0.4.0    2021-05-12 [?] CRAN (R 4.1.2)
##  P ScaledMatrix           1.2.0    2021-10-26 [?] Bioconductor
##  P scales                 1.2.1    2022-08-20 [?] CRAN (R 4.1.2)
##  P scattermore            0.7      2020-11-24 [?] CRAN (R 4.1.2)
##  P scatterpie             0.2.2    2024-04-03 [?] CRAN (R 4.1.2)
##  P sctransform            0.3.3    2022-01-13 [?] CRAN (R 4.1.2)
##  P sessioninfo            1.2.2    2021-12-06 [?] CRAN (R 4.1.2)
##  P Seurat               * 4.0.0    2021-01-30 [?] CRAN (R 4.1.2)
##  P SeuratObject         * 4.0.4    2021-11-23 [?] CRAN (R 4.1.2)
##  P shadowtext             0.1.3    2024-01-19 [?] CRAN (R 4.1.2)
##  P shape                  1.4.6    2021-05-19 [?] CRAN (R 4.1.2)
##  P shiny                  1.7.1    2021-10-02 [?] CRAN (R 4.1.2)
##  P SingleCellExperiment   1.16.0   2021-10-26 [?] Bioconductor
##  P sparseMatrixStats      1.6.0    2021-10-26 [?] Bioconductor
##  P spatstat               1.64-1   2020-05-12 [?] CRAN (R 4.1.2)
##  P spatstat.data          2.1-2    2021-12-17 [?] CRAN (R 4.1.2)
##  P spatstat.utils         2.3-0    2021-12-12 [?] CRAN (R 4.1.2)
##  P stringi                1.7.6    2021-11-29 [?] CRAN (R 4.1.2)
##  P stringr              * 1.5.0    2022-12-02 [?] CRAN (R 4.1.2)
##  P SummarizedExperiment   1.24.0   2021-10-26 [?] Bioconductor
##  P survival               3.2-13   2021-08-24 [?] CRAN (R 4.1.2)
##  P tensor                 1.5      2012-05-05 [?] CRAN (R 4.1.2)
##  P tibble               * 3.2.1    2023-03-20 [?] RSPM (R 4.1.2)
##  P tidygraph              1.3.1    2024-01-30 [?] CRAN (R 4.1.2)
##  P tidyr                * 1.3.0    2023-01-24 [?] CRAN (R 4.1.2)
##  P tidyselect             1.2.0    2022-10-10 [?] CRAN (R 4.1.2)
##  P tidytree               0.4.6    2023-12-12 [?] CRAN (R 4.1.2)
##  P treeio                 1.29.0   2024-09-20 [?] Github (GuangchuangYu/treeio@9504617)
##  P tweenr                 1.0.2    2021-03-23 [?] CRAN (R 4.1.2)
##  P tzdb                   0.3.0    2022-03-28 [?] CRAN (R 4.1.2)
##  P urlchecker             1.0.1    2021-11-30 [?] CRAN (R 4.1.2)
##  P usethis                2.2.2    2023-07-06 [?] CRAN (R 4.1.2)
##  P utf8                   1.2.2    2021-07-24 [?] CRAN (R 4.1.2)
##  P uwot                   0.1.11   2021-12-02 [?] CRAN (R 4.1.2)
##  P vctrs                  0.6.5    2023-12-01 [?] CRAN (R 4.1.2)
##  P viridis              * 0.5.1    2018-03-29 [?] RSPM (R 4.1.2)
##  P viridisLite          * 0.3.0    2018-02-01 [?] CRAN (R 4.1.2)
##  P withr                  2.5.0    2022-03-03 [?] CRAN (R 4.1.2)
##  P xfun                   0.29     2021-12-14 [?] CRAN (R 4.1.2)
##  P XML                    3.99-0.8 2021-09-17 [?] CRAN (R 4.1.2)
##  P xtable                 1.8-4    2019-04-21 [?] CRAN (R 4.1.2)
##  P XVector                0.34.0   2021-10-26 [?] Bioconductor
##  P yaml                   2.2.1    2020-02-01 [?] CRAN (R 4.1.2)
##  P yulab.utils            0.1.4    2024-01-28 [?] CRAN (R 4.1.2)
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
