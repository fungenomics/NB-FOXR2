---
title: "01 - Oncoprint"
author: "Selin Jessa [[selin.jessa@mail.mcgill.ca](mailto:selin.jessa@mail.mcgill.ca)] and Bhavyaa Chandarana [[bhavyaa.chandarana@mail.mcgill.ca](mailto:bhavyaa.chandarana@mail.mcgill.ca)]"
date: "09 October, 2024"
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
## Document index: 01
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
## public/output/01
```

```
## public/figures/01//
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

Here we produce an oncoprint as an overview of the dataset. This code
is adapted from Steven Hébert.

# Libraries


```r
knitr::opts_chunk$set(cache = FALSE)

library(here)
library(data.table)
library(magrittr)
library(tidyr)
library(dplyr)
library(readr)
library(glue)
library(ggplot2)
library(ComplexHeatmap)

source(here("include/style.R"))
source(here("code/functions/scRNAseq.R"))

ggplot2::theme_set(theme_min())
```

# Load metadata


```r
meta <- read_tsv(here("output/00/oncoprint.tsv"))
```

```
## Rows: 125 Columns: 13
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (13): Sample, ID, FOXR2_positive, Group, Source, Age, Sex, Location, RNA...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
data <- as.data.table(meta)

# Show counts per group
table(data$Group)
```

```
## 
##       DIPG-H3K27M DIPG-H3K27M-FOXR2            EP-PFA              ETMR 
##                19                 6                13                12 
##    HGG-H3.3G34R/V           HGG-IDH            MB-SHH            MB-WNT 
##                18                10                 8                10 
##          NB-FOXR2 
##                29
```

# Set heatmap annotation


```r
palette_location <- c("Cerebellum"      = "goldenrod3",
                     "Hemispheric"     = "navy",
                     "Pons"            = "gold1",
                     "Posterior Fossa" = "orange",
                     "Thalamus"        = "darkolivegreen4",
                     "Unknown"         = "gray80")

palette_bulkRNAseq   <- c("Y" = "black",
                          "N" = "gray80")

palette_scRNAseq     <- c("Y" = "black",
                          "N" = "gray80")

palette_scMultiome   <- c("Y" = "black",
                          "N" = "gray80")

palette_newSample    <- c("Y" = "black",
                          "N" = "gray80")

palette_FOXR2positive <- c("Y" = "black",
                          "N" = "gray80")

palette_sex <- c("M" = "black",
                 "F" = "gray80",
                 "Unknown" = "gray50")

# Set the order of the samples
samples_order <- data$ID

# Add additional info (row) on the bottom of the oncoprint
bot_ann <- HeatmapAnnotation(
    annotation_name_gp = gpar(fontsize = 10),
    Group = data$Group,
    Location = data$Location_broad,
    FOXR2_positive = data$FOXR2_positive,
    Sex = data$Sex,
    New_sample = data$New_sample,
    bulkRNAseq = data$RNAseq,
    scRNAseq = data$scRNAseq,
    `Age(year)` = anno_points(as.matrix(as.numeric(data$Age)), border = F, size = unit(1, "mm")),
    col = list(Group = palette_groups,
               Location = palette_location,
               FOXR2_positive = palette_FOXR2positive,
               Sex = palette_sex,
               New_sample = palette_newSample,
               bulkRNAseq = palette_bulkRNAseq,
               scRNAseq = palette_scRNAseq),
    na_col = "white",
    gp = gpar(col = "white"))
```

# Create oncoprint


```r
# plot just the heatmap annotation by including an empty 0-row heatmap
# following: https://github.com/jokergoo/ComplexHeatmap/issues/301#issuecomment-487481520
Heatmap(matrix(nc = 125, nr = 0),
        bottom_annotation = bot_ann)
```

![](/project/kleinman/bhavyaa.chandarana/from_hydra/2023-05-NB-FOXR2/public/figures/01//oncoprint-1.png)<!-- -->

<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2024-10-09 18:18:33
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
##  os       Rocky Linux 8.7 (Green Obsidian)
##  system   x86_64, linux-gnu
##  ui       X11
##  language (EN)
##  collate  en_US.UTF-8
##  ctype    en_US.UTF-8
##  tz       America/Toronto
##  date     2024-10-09
##  pandoc   1.19.2.1 @ /cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/ (via rmarkdown)
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  ! package        * version date (UTC) lib source
##  P BiocGenerics     0.40.0  2021-10-26 [?] Bioconductor
##  P BiocManager      1.30.15 2021-05-11 [?] CRAN (R 4.1.2)
##  P bit              4.0.4   2020-08-04 [?] CRAN (R 4.1.2)
##  P bit64            4.0.5   2020-08-30 [?] CRAN (R 4.1.2)
##  P bslib            0.3.1   2021-10-06 [?] CRAN (R 4.1.2)
##  P cachem           1.0.6   2021-08-19 [?] CRAN (R 4.1.2)
##  P callr            3.7.6   2024-03-25 [?] RSPM
##  P circlize         0.4.15  2022-05-10 [?] CRAN (R 4.1.2)
##  P cli              3.6.1   2023-03-23 [?] RSPM (R 4.1.2)
##  P clue             0.3-64  2023-01-31 [?] CRAN (R 4.1.2)
##  P cluster          2.1.2   2021-04-17 [?] CRAN (R 4.1.2)
##  P codetools        0.2-18  2020-11-04 [?] CRAN (R 4.1.2)
##  P colorspace       2.0-2   2021-06-24 [?] CRAN (R 4.1.2)
##  P ComplexHeatmap * 2.10.0  2021-10-26 [?] Bioconductor
##  P crayon           1.4.2   2021-10-29 [?] CRAN (R 4.1.2)
##  P data.table     * 1.14.2  2021-09-27 [?] CRAN (R 4.1.2)
##  P devtools         2.4.5   2022-10-11 [?] CRAN (R 4.1.2)
##  P digest           0.6.35  2024-03-11 [?] CRAN (R 4.1.2)
##  P doParallel       1.0.16  2020-10-16 [?] CRAN (R 4.1.2)
##  P dplyr          * 1.1.1   2023-03-22 [?] CRAN (R 4.1.2)
##  P ellipsis         0.3.2   2021-04-29 [?] CRAN (R 4.1.2)
##  P evaluate         0.23    2023-11-01 [?] CRAN (R 4.1.2)
##  P fansi            1.0.2   2022-01-14 [?] CRAN (R 4.1.2)
##  P fastmap          1.1.0   2021-01-25 [?] CRAN (R 4.1.2)
##  P foreach          1.5.1   2020-10-15 [?] CRAN (R 4.1.2)
##  P fs               1.5.2   2021-12-08 [?] CRAN (R 4.1.2)
##  P generics         0.1.3   2022-07-05 [?] CRAN (R 4.1.2)
##  P GetoptLong       1.0.5   2020-12-15 [?] CRAN (R 4.1.2)
##  P ggplot2        * 3.4.2   2023-04-03 [?] CRAN (R 4.1.2)
##  P git2r            0.29.0  2021-11-22 [?] CRAN (R 4.1.2)
##  P GlobalOptions    0.1.2   2020-06-10 [?] CRAN (R 4.1.2)
##  P glue           * 1.6.2   2022-02-24 [?] CRAN (R 4.1.2)
##  P gridExtra        2.3     2017-09-09 [?] CRAN (R 4.1.2)
##  P gtable           0.3.0   2019-03-25 [?] CRAN (R 4.1.2)
##  P here           * 1.0.1   2020-12-13 [?] CRAN (R 4.1.2)
##  P highr            0.9     2021-04-16 [?] CRAN (R 4.1.2)
##  P hms              1.1.1   2021-09-26 [?] CRAN (R 4.1.2)
##  P htmltools        0.5.2   2021-08-25 [?] CRAN (R 4.1.2)
##  P htmlwidgets      1.5.4   2021-09-08 [?] CRAN (R 4.1.2)
##  P httpuv           1.6.5   2022-01-05 [?] CRAN (R 4.1.2)
##  P IRanges          2.28.0  2021-10-26 [?] Bioconductor
##  P iterators        1.0.13  2020-10-15 [?] CRAN (R 4.1.2)
##  P jquerylib        0.1.4   2021-04-26 [?] CRAN (R 4.1.2)
##  P jsonlite         1.8.8   2023-12-04 [?] CRAN (R 4.1.2)
##  P knitr            1.37    2021-12-16 [?] CRAN (R 4.1.2)
##  P later            1.3.0   2021-08-18 [?] CRAN (R 4.1.2)
##  P lifecycle        1.0.3   2022-10-07 [?] CRAN (R 4.1.2)
##  P magrittr       * 2.0.3   2022-03-30 [?] CRAN (R 4.1.2)
##  P matrixStats      0.61.0  2021-09-17 [?] CRAN (R 4.1.2)
##  P memoise          2.0.1   2021-11-26 [?] CRAN (R 4.1.2)
##  P mime             0.12    2021-09-28 [?] CRAN (R 4.1.2)
##  P miniUI           0.1.1.1 2018-05-18 [?] CRAN (R 4.1.2)
##  P munsell          0.5.0   2018-06-12 [?] CRAN (R 4.1.2)
##  P pillar           1.9.0   2023-03-22 [?] RSPM (R 4.1.2)
##  P pkgbuild         1.4.2   2023-06-26 [?] CRAN (R 4.1.2)
##  P pkgconfig        2.0.3   2019-09-22 [?] CRAN (R 4.1.2)
##  P pkgload          1.3.3   2023-09-22 [?] CRAN (R 4.1.2)
##  P png              0.1-7   2013-12-03 [?] CRAN (R 4.1.2)
##  P prettyunits      1.1.1   2020-01-24 [?] CRAN (R 4.1.2)
##  P processx         3.8.4   2024-03-16 [?] RSPM
##  P profvis          0.3.8   2023-05-02 [?] CRAN (R 4.1.2)
##  P promises         1.2.0.1 2021-02-11 [?] CRAN (R 4.1.2)
##  P ps               1.7.6   2024-01-18 [?] RSPM
##  P purrr            1.0.1   2023-01-10 [?] CRAN (R 4.1.2)
##  P R6               2.5.1   2021-08-19 [?] CRAN (R 4.1.2)
##  P RColorBrewer   * 1.1-2   2014-12-07 [?] CRAN (R 4.1.2)
##  P Rcpp             1.0.8   2022-01-13 [?] CRAN (R 4.1.2)
##  P readr          * 2.1.1   2021-11-30 [?] CRAN (R 4.1.2)
##  P remotes          2.4.2.1 2023-07-18 [?] CRAN (R 4.1.2)
##  P renv             1.0.3   2023-09-19 [?] CRAN (R 4.1.2)
##  P rjson            0.2.21  2022-01-09 [?] CRAN (R 4.1.2)
##  P rlang            1.1.3   2024-01-10 [?] CRAN (R 4.1.2)
##  P rmarkdown        2.11    2021-09-14 [?] CRAN (R 4.1.2)
##  P rprojroot        2.0.2   2020-11-15 [?] CRAN (R 4.1.2)
##  P S4Vectors        0.32.4  2022-03-24 [?] Bioconductor
##  P sass             0.4.0   2021-05-12 [?] CRAN (R 4.1.2)
##  P scales           1.2.1   2022-08-20 [?] CRAN (R 4.1.2)
##  P sessioninfo      1.2.2   2021-12-06 [?] CRAN (R 4.1.2)
##  P shape            1.4.6   2021-05-19 [?] CRAN (R 4.1.2)
##  P shiny            1.7.1   2021-10-02 [?] CRAN (R 4.1.2)
##  P stringi          1.7.6   2021-11-29 [?] CRAN (R 4.1.2)
##  P stringr          1.5.0   2022-12-02 [?] CRAN (R 4.1.2)
##  P tibble           3.2.1   2023-03-20 [?] RSPM (R 4.1.2)
##  P tidyr          * 1.3.0   2023-01-24 [?] CRAN (R 4.1.2)
##  P tidyselect       1.2.0   2022-10-10 [?] CRAN (R 4.1.2)
##  P tzdb             0.3.0   2022-03-28 [?] CRAN (R 4.1.2)
##  P urlchecker       1.0.1   2021-11-30 [?] CRAN (R 4.1.2)
##  P usethis          2.2.2   2023-07-06 [?] CRAN (R 4.1.2)
##  P utf8             1.2.2   2021-07-24 [?] CRAN (R 4.1.2)
##  P vctrs            0.6.5   2023-12-01 [?] CRAN (R 4.1.2)
##  P viridis        * 0.5.1   2018-03-29 [?] RSPM (R 4.1.2)
##  P viridisLite    * 0.3.0   2018-02-01 [?] CRAN (R 4.1.2)
##  P vroom            1.5.7   2021-11-30 [?] CRAN (R 4.1.2)
##  P withr            2.5.0   2022-03-03 [?] CRAN (R 4.1.2)
##  P xfun             0.29    2021-12-14 [?] CRAN (R 4.1.2)
##  P xtable           1.8-4   2019-04-21 [?] CRAN (R 4.1.2)
##  P yaml             2.2.1   2020-02-01 [?] CRAN (R 4.1.2)
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
