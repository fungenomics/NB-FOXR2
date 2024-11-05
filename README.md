[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13755696.svg)](https://doi.org/10.5281/zenodo.13755696)

# NB-FOXR2 analysis code

- This is the **public** repository accompanying the study: FOXR2 targets LHX6+/DLX+ neural lineages to drive CNS neuroblastoma, Jessa*, De Cola*, Chandarana*, ..., Pathania, Jabado, Kleinman, _Cancer Research_, 2024 ([link](https://aacrjournals.org/cancerres/article/doi/10.1158/0008-5472.CAN-24-2248/749689/FOXR2-targets-LHX6-DLX-neural-lineages-to-drive))
- Link to repository: https://github.com/fungenomics/NB-FOXR2
- This repository is archived on [Zenodo](https://www.doi.org/10.5281/zenodo.13755695)
- This repository contains primarily code, see the **Data availability** section for links to data
- See the  [**Code to reproduce analyses**](https://github.com/fungenomics/NB-FOXR2#code-to-reproduce-analyses) section for links to code of each figure and table

## Contents of this repository

### Codebase overview

- This repository is meant to enhance the Materials & Methods section by providing code for the custom
analyses in the manuscript and the exact R dependencies, in order to improve reproducibility for the main results.
However, it is not a fully executable workflow.
- In general, read alignment and cell calling for single-cell sequencing data has been performed using in-house pipelines, not included here.
This repository mainly contains custom/downstream code.
- A first level of downstream analysis involves scripts applied in parallel to individual samples for specific data types.
Copies of these scripts are provided in the `code/scripts/per_sample_script_examples` (the execution of these scripts is performed in the `data/singlecell`, `data/RNAseq`, etc folders, not included here).
- A second level of downstream analysis involves custom analyses, aggregating samples and data types, use to derive the main results included in the paper. These are provided in .Rmd files in `code`, with the associated .md and rendered HTML files. 

### Codebase structure

* `renv.lock` --> lockfile containing all package versions for R analysis
* `code` --> code for R analysis, contains the .Rmd files that run the high-level analyses and produce figures included in the paper
   * `functions` --> contains .R files with custom functions used throughout the analysis
   * `scripts` --> contains .R and bash scripts for analyses that are repeated on individual samples, as well as helper scripts e.g. for formatting data
* `include` --> contains templates, palettes, etc, for this repository
* `rr_helpers.R` --> contains helper functions for working with this GitHub repository template ([`rr`](https://github.com/sjessa/rr))

### Analysis documents

Rendered HTML files for each analysis can be viewed here:

* [01: Oncoprint](https://fungenomics.github.io/NB-FOXR2/code/01-oncoprint.html)
* [02: Prepare reference gene signatures](https://fungenomics.github.io/NB-FOXR2/code/02-prep_ref_signatures.html)
* [03: Process human tumor bulk RNAseq](https://fungenomics.github.io/NB-FOXR2/code/03-bulk_RNAseq.html)
* [04: Process human tumor single-cell RNAseq](https://fungenomics.github.io/NB-FOXR2/code/04-sc_tumors.html)
* [05: Transcription factor fingerprint](https://fungenomics.github.io/NB-FOXR2/code/05-TF_fingerprint.html)
* [06.1: Project bulk human tumors to normal references](https://fungenomics.github.io/NB-FOXR2/code/06.1-bulk_tumor_projections.html)
* [06.2: Project single-cell human tumors to normal references](https://fungenomics.github.io/NB-FOXR2/code/06.2-sc_tumor_projections.html)
* [07.1: Project single-cell mouse models to normal references](https://fungenomics.github.io/NB-FOXR2/code/07.1-sc_mm_models_projections.html)
* [07.2: Validate similarity of mouse models to human disease](https://fungenomics.github.io/NB-FOXR2/code/07.2-mm_models_validation.html)
* [08: Compare to extracranial neuroblastoma](https://fungenomics.github.io/NB-FOXR2/code/08-extracranial_NB.html)
* [09: Process and analyse FOXR2 CUT&RUN chromatin binding](https://fungenomics.github.io/NB-FOXR2/code/09-ChIP.html)

## Code to reproduce analyses

Code to reproduce analyses is saved in `code/`.

### Figures

Source of each Main Figure in the R markdowns:

Figure number | Description | R markdown
-- | -- | --
1a | Oncoprint | [01](https://fungenomics.github.io/NB-FOXR2/code/01-oncoprint.html)
1b | FOXR2 expression in bulk human brain tumors | [03](https://fungenomics.github.io/NB-FOXR2/code/03-bulk_RNAseq.html)
2a | TF fingerprint expression in normal brain | NA - Custom schematic + imaging from online browsers ([link](https://db.cngb.org/stomics/mosta/spatial/), [link](https://developingmouse.brain-map.org/))
2b-c | TF fingerprint SVM and expression in normal development | [05](https://fungenomics.github.io/NB-FOXR2/code/05-TF_fingerprint.html)
3a, c-d | TF fingerprint expression in selected bulk human tumors | [05](https://fungenomics.github.io/NB-FOXR2/code/05-TF_fingerprint.html)
3b | Comparison of TF expression in bulk FOXR2+/- extracranial neuroblastoma (Gartlgruber 2021 dataset) | [08](https://fungenomics.github.io/NB-FOXR2/code/08-extracranial_NB.html)
4a-b, e-f | ssGSEA scoring in bulk human tumors | [06.1](https://fungenomics.github.io/NB-FOXR2/code/06.1-bulk_tumor_projections.html)
4c (top), 4d | Single-cell NB-FOXR2 mapping | [06.2](https://fungenomics.github.io/NB-FOXR2/code/06.2-sc_tumor_projections.html)
4c (bottom) | Single-cell NB-FOXR2 colored by sample | [04](https://fungenomics.github.io/NB-FOXR2/code/04-sc_tumors.html)
5 | Survival and imaging of murine models | NA - figures from Pathania lab
6a-b | Single-cell murine models inferCNV and UMAP colored by broad cell classes | [07.1](https://fungenomics.github.io/NB-FOXR2/code/07.1-sc_mm_models_projections.html)
6c-f | Human tumor signature generation and ssGSEA scoring in murine models (neuron-like cells) | [07.2](https://fungenomics.github.io/NB-FOXR2/code/07.2-mm_models_validation.html)
7a-c | Mouse model qPCR and survival curve | NA - figures from Pathania lab
7d-e | FOXR2 CUT&RUN peak overlap and motif enrichment | [09](https://fungenomics.github.io/NB-FOXR2/code/09-ChIP.html)
7f-i | Genomic tracks of CUT&RUN, ATAC, RNAseq | NA - figures from IGV 


Source of each Supplementary Figure in the R markdowns:

Supplementary figure number | Title | R markdown
-- | -- | --
1a-b | Single-cell tumor QC | [04](https://fungenomics.github.io/NB-FOXR2/code/04-sc_tumors.html)
1c | inferCNV in human tumors | NA - png files from `data/singlecell/pipeline_sc*/{sample}/inferCNV`
1d | MDM4 expression in human bulk tumors | [03](https://fungenomics.github.io/NB-FOXR2/code/03-bulk_RNAseq.html)
1e | MDM4 expression in human single-cell tumors | [06.2](https://fungenomics.github.io/NB-FOXR2/code/06.2-sc_tumor_projections.html)
2 | TF fingerprint expression in individual clusters of normal development | [05](https://fungenomics.github.io/NB-FOXR2/code/05-TF_fingerprint.html)
3a-b | TF fingerprint expression in all human bulk tumors | [05](https://fungenomics.github.io/NB-FOXR2/code/05-TF_fingerprint.html)
3c | Comparison of TF expression in bulk FOXR2+/- extracranial neuroblastoma (TARGET dataset) | [08](https://fungenomics.github.io/NB-FOXR2/code/08-extracranial_NB.html)
4 | Top ssGSEA scoring signatures per human bulk brain tumor, indiscriminant signature (see Methods) | [06.1](https://fungenomics.github.io/NB-FOXR2/code/06.1-bulk_tumor_projections.html)
5a | OL marker expression in human bulk brain tumors | [05](https://fungenomics.github.io/NB-FOXR2/code/05-TF_fingerprint.html)
5b | Expression of neuronal markers in normal brain | NA - from [online browser](https://celltypes.brain-map.org/rnaseq/human_ctx_smart-seq)
5c | Neuronal marker expression in human single-cell NB-FOXR2 | [06.2](https://fungenomics.github.io/NB-FOXR2/code/06.2-sc_tumor_projections.html)
6a, c-d | Single-cell murine models cell type annotation and marker expression | [07.1](https://fungenomics.github.io/NB-FOXR2/code/07.1-sc_mm_models_projections.html)
6b | inferCNV in malignant cells of murine models | NA - png files from `data/singlecell/mm_inferCNV/{sample}/malignant_only`
6e | Human tumor signature scoring in glial-like cells of murine models | [07.2](https://fungenomics.github.io/NB-FOXR2/code/07.2-mm_models_validation.html)
7 | Imaging of orthotopic engraftment | NA - figures from Pathania lab


### Tables

Source of each Supplementary Table in the R markdowns (denoted by section titled `TABLE...`):

Supplementary table number | Title | R markdown
-- | -- | --
1 | Murine models | NA – murine experimental information from Pathania lab
2 | NGS summary | NA - generated external to this repo
3 | Human tumor bulk QC | [03](https://fungenomics.github.io/NB-FOXR2/code/03-bulk_RNAseq.html)
4 | Human tumor singlecell QC | [04](https://fungenomics.github.io/NB-FOXR2/code/04-sc_tumors.html)
5 | Murine models singlecell QC | [07.1](https://fungenomics.github.io/NB-FOXR2/code/07.1-sc_mm_models_projections.html)
6 | Normal reference datasets | NA – created manually, collects references used in 02, 05, 06.2, 07.1
7 | Murine models annotation reference labels | [07.1](https://fungenomics.github.io/NB-FOXR2/code/07.1-sc_mm_models_projections.html)
8 | TF quantification | [05](https://fungenomics.github.io/NB-FOXR2/code/05-TF_fingerprint.html)
9 | Reference gene signatures | [06.1](https://fungenomics.github.io/NB-FOXR2/code/06.1-bulk_tumor_projections.html)
10 | Cell type signature enrichment | [06.1](https://fungenomics.github.io/NB-FOXR2/code/06.1-bulk_tumor_projections.html)


## Data availability

### Data from this study

- Raw and processed data for murine samples have been deposited to GEO at
[GSE270666](https://www.ncbi.xyz/geo/query/acc.cgi?acc=GSE270666).
- Raw data for human tumors have been deposited to EGA at
[EGAS00001007247](https://ega-archive.org/dacs/EGAC00001003262).
- For human samples, please see the associated [Zenodo record](https://doi.org/10.5281/zenodo.13750919)
for processed data including counts matrices, cell annotations, and bigWig files. 


### External data

- [Gartlgruber 2021](https://www.nature.com/articles/s43018-020-00145-w) dataset 
of extracranial neuroblastoma (bulk RNA-seq) analyzed in this study were obtained
as count matrices from a [public Shiny web app](https://nbseb087.dkfz.de/project_NB_SE_viz/)
provided by the authors. 

- [TARGET](https://gdc.cancer.gov/content/target-nbl-publication-summary)
dataset of extracranial neuroblastoma (bulk RNA-seq) analyzed in this 
study were obtained as count matrices from the public [Genomic Data Commons (GDC) 
cancer portal](https://portal.gdc.cancer.gov/) as described in R markdown `08`.


## GitHub / version control

The following are tracked / available on GitHub:

* `.Rmd` files, containing the code, and `.md` and rendered HTML files, containing code and outputs
* The lockfile produced by the `renv` package

The following are not tracked / available on GitHub:

* Figures in `png`/`pdf` format
* Raw data and analysis output / processed data files (see [**Data availability**](https://github.com/fungenomics/NB-FOXR2#data-availability) section above)
* The actual packages in the R library 

## Citation

If you use or modify code provided here, please cite this work as follows:

> Selin Jessa, Bhavyaa Chandarana, Steven Hébert, and Claudia L. Kleinman. (2024). NB-FOXR2 analysis code. Zenodo. https://www.doi.org/10.5281/zenodo.13755695





