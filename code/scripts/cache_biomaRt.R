# Script to download human & mouse biomaRt and cache R objects 
# Run interactively on login node, since internet is required

library(here)
require(biomaRt)
require(data.table)

# Get biomaRt
human  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse  <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

# get mm10 genes with ENS ids
mm10_gtf   <- rtracklayer::import("/lustre06/project/6004736/singlecell_pipeline/references/cellranger/refdata-cellranger-mm10-1.2.0/genes/genes.gtf")
mm10_genes <- unique(mm10_gtf$gene_name)

# Store mm10 ensembl IDs for gene symbols
genes_mm_ens_lds <- getLDS(attributes = c("mgi_symbol"),
                 filters = "mgi_symbol",
                 values = mm10_genes,
                 mart = mouse,
                 attributesL = c("ensembl_gene_id"),
                 martL = mouse,
                 uniqueRows = TRUE)
save(genes_mm_ens_lds, file = here("data/singlecell/references_genome/biomaRt_mm_symbol_to_ens_lds.Rda"))

# Store mm10 to hg19 symbol conversion
genes_lds <- getLDS(attributes = c("mgi_symbol"),
                 filters = "mgi_symbol",
                 values = mm10_genes,
                 mart = mouse,
                 attributesL = c("hgnc_symbol"),
                 martL = human,
                 uniqueRows = TRUE)
save(genes_lds, file = here("data/singlecell/references_genome/biomaRt_mm_to_hg_lds.Rda"))

# get hg19 genes with ENS ids
hg19_gtf   <- rtracklayer::import("/lustre06/project/6004736/singlecell_pipeline/references/cellranger/refdata-cellranger-mm10-1.2.0/genes/genes.gtf")
hg19_genes <- unique(hg19_gtf$gene_name)

# Store hg19 ensembl IDs for gene symbols
genes_hg_ens_lds <- getLDS(attributes = c("hgnc_symbol"),
                 filters = "hgnc_symbol",
                 values = hg19_genes,
                 mart = human,
                 attributesL = c("ensembl_gene_id"),
                 martL = human,
                 uniqueRows = TRUE)
save(genes_hg_ens_lds, file = here("data/singlecell/references_genome/biomaRt_hg_symbol_to_ens_lds.Rda"))




# Store Refseq IDs for all mm10 genes
genes_mgi_refseq <- getLDS(attributes = c("mgi_symbol"),
        filters = "mgi_symbol", 
        values = mm10_genes,
        mart = mouse,
        attributesL = c("refseq_mrna"), 
        martL = mouse,
        uniqueRows = TRUE)
save(genes_mgi_refseq, file = here("data/singlecell/references_genome/biomaRt_mm_symbol_to_refseq_lds.Rda"))

# Store hg19 to mm ensembl IDs
genes_hg_mm_ens <- getLDS(attributes = c("ensembl_gene_id"),
        filters = "hgnc_symbol", 
        values = hg19_genes,
        mart = human,
        attributesL = c("ensembl_gene_id"), 
        martL = mouse,
        uniqueRows = TRUE)
save(genes_hg_mm_ens, file = here("data/singlecell/references_genome/biomaRt_ens_mm_to_hg_lds.Rda"))