library(here)
library(gprofiler2)
require(glue)
require(data.table)

# Set output directory (for ChIP analysis)
out <- here("output/09/")

# load top genes and background genes
# "top_bound_genes" and "bg_genes"
load(glue("{out}/GO_gene_lists.Rda"))

gost <- gprofiler2::gost(query = top_bound_genes, 
                            organism = "mmusculus", ordered_query = TRUE, 
                            significant = TRUE, evcodes = TRUE, 
                            custom_bg = bg_genes)

save(gost, file = glue("{out}/gost.Rda"))

