
# bulk RNAseq pipeline helpers -------------------------------------------------

#' Based off of code from Nicolas De Jay. Extract counts for genes of interest from
#' the raptor bulk RNA-seq pipeline and wrangle them into tidy format
#'
#' @param path Character, path to counts file
#' @param goi Character, vector of one or more gene symbols for genes of interest
#' for which to retrieve counts
#'
#' @return Data frame with four columns: sample, gene_expression, gene_ensg, gene_symbol
extract_pipeline_counts <- function(path, goi, long = TRUE) {
  
  counts <- read.table(path, header = T, sep = "\t", check.names = FALSE)
  
  # Our pipeline labels genes using #{ENSG}:#{SYMBOL}
  ensg_to_symbol <- data.frame("gene_id" = rownames(counts) %>% as.character) %>%
    dplyr::mutate("gene_ensg"   = gene_id %>% as.character %>% strsplit(":") %>% sapply(`[[`, 1)) %>%
    dplyr::mutate("gene_symbol" = gene_id %>% as.character %>% strsplit(":") %>% sapply(`[[`, 2))
  
  # Convert all genes from #{SYMBOL} to #{ENSG}:#{SYMBOL}
  genes.symbol <- data.frame(genes = goi)
  colnames(genes.symbol) = "gene_symbol"
  
  genes.ensg <- genes.symbol %>% left_join(ensg_to_symbol) %>% dplyr::select(gene_id) %>% 
    filter(!is.na(gene_id))
  
  # Subset counts table
  counts.subset <- counts[genes.ensg %>% pull(gene_id) %>% as.character, , drop = F] %>%
    as.matrix()
  
  if (!long) return(counts.subset)
  
  counts.subset <- counts.subset %>%
    reshape2::melt() %>%
    setNames(c("gene_id", "sample", "gene_expression")) %>%
    left_join(ensg_to_symbol, by = "gene_id") %>%
    dplyr::select(-gene_id)
  
  return(counts.subset)
  
}


#' Convert Ensembl gene IDs to gene symbols
ensembl2symbols_safe <- function(genes, ...) {
  
  sym <- lapply(genes, ensembl2symbols, ...)
  empty <- which(unlist(map(sym, ~ identical(unlist(.x), character(0)))))
  sym[empty] <- genes[empty]
  
  return(unlist(sym))
  
}



#' Define a function to summarize the cell types into broader "classes",
#' distinguishing between different types of inhibitory neurons.
#' 
#' @param df Data frame where one column of cell type labels needs to be summarized
#' into classes
#' @param cluster_col Unquoted expression, name of the column containing
#' the cluster cell type labels
#' 
#' @example 
#' x <- data.frame(Celltype = c("Oligodendrocyte", "Astrocyte"),
#'            Score    = c(100, 150))
#' summarizeCellTypes3(x, Celltype)
summarizeCellTypes3 <- function(df, cluster_col) {
  
  cc_quo <- enquo(cluster_col)
  
  df %>%
    mutate(
      # Define some broader cell type classes
      Type = case_when(
        grepl("RG|[Rr]adial|NSC|prog", !!cc_quo) & !grepl("NRGN", !!cc_quo) & !grepl("[Oo]lig", !!cc_quo) & !grepl("GE", !!cc_quo) & !grepl("-P.{0,1}$|prolif", !!cc_quo) ~ "RGC",
        grepl("INIP|MGIN|MGE|CGE|LGE|SST|PV|[Ii]nhib|CIN", !!cc_quo) ~ "Prenatal inhib. neurons",
        grepl("PV|SST|VIP|SV2C|Somato|Parv", !!cc_quo) ~ "Pediatric/adult inhib. neurons",
        grepl("EXIP|NEURP|IP|[Ii]ntermediate|prog", !!cc_quo) & !grepl("VIP", !!cc_quo) & !grepl("[Oo]lig", !!cc_quo) ~ "Neuronal progenitors",
        grepl("CEX|PEX|[Ee]xcit|NRGN|[Nn]eu|[Cc]ortex|SPN|MFN|CJRZN|UBC|GABAN|NEUR", !!cc_quo) ~ "Other neurons",
        grepl("OPC|Oligodendrocyte progenitor cell|Oligodendrocyte precurso|COP", !!cc_quo) ~ "Oligodendrocyte precursors",
        grepl("NFOL|MOL|[Oo]ligo", !!cc_quo) ~ "Oligodendrocytes",
        grepl("ASTR|Astr", !!cc_quo) ~ "Astrocytes",
        grepl("T|Fibr|Mur|Micro|Mixed", !!cc_quo) ~ "Other",
        TRUE ~ "Other")) %>%
    mutate(Type = factor(Type, levels = c("RGC",
                                          "RGC (prolif.)",
                                          "Neuronal progenitors",
                                          "Prenatal inhib. neurons",
                                          "Pediatric/adult inhib. neurons",
                                          "Other neurons",
                                          "Glial progenitors",
                                          "Oligodendrocyte precursors",
                                          "Oligodendrocytes",
                                          "Astrocytes",
                                          "Choroid/ependymal",
                                          "Immune",
                                          "Non-neuroectoderm",
                                          "Other")))
  
}



# Helpers for GSEA analysis ----------------------------------------------------




get_gene_ranks <- function(pathway, stats, direction = "enriched") {
  
  if (direction == "enriched") rnk <- rank(-stats)
  else if (direction == "depleted") rnk <- rank(stats)
  
  all_genes_ranked <- sort(rnk)
  pathway_genes_ranked <- all_genes_ranked[pathway[pathway %in% names(stats)]]
  
  x <- data.frame(gene = names(pathway_genes_ranked),
                  rank = unname(pathway_genes_ranked)) %>%
    arrange(rank)
  
  x
  
}


prep_gsea_heatmap <- function(fgsea_df, signatures, filters, row_order) {
  
  heatmap_data_long <- fgsea_df %>%
    filter(!!! filters) %>% 
    filter(Signature %in% signatures)
  # # Subset to relevant cell types
  # summarizeCellTypes3(Cell_type) %>%
  # filter(Type %in% c("RGC",
  #                    "Prenatal inhib. neurons",
  #                    "Pediatric/adult inhib. neurons",
  #                    "Neuronal progenitors",
  #                    "Other neurons",
  #                    "Oligodendrocytes",
  #                    "Oligodendrocyte precursors",
  #                    "Astrocytes")) %>%
  # mutate(Cell_type = paste0(Age, " ", Cell_type))
  
  col_order <- heatmap_data_long %>%
    distinct(Signature, Dataset, Cell_type) %>%
    arrange(Dataset, Cell_type) %>%
    pull(Signature)
  
  col_anno <- heatmap_data_long %>%
    distinct(Signature, Cell_type, Dataset, Species) %>%
    # mutate(Species = ifelse(grepl("Mouse", Cell_type), "Mouse", "Human")) %>%
    data.frame() %>%
    select(-Cell_type) %>% 
    tibble::column_to_rownames(var = "Signature")
  
  heatmap_data_wide <- heatmap_data_long %>%
    select(Comparison, Signature, NES) %>%
    tidyr::spread(Signature, NES) %>%
    as.data.frame() %>%
    tibble::column_to_rownames(var = "Comparison") %>%
    .[row_order, col_order]
  
  signif_data_wide <- heatmap_data_long %>%
    select(Comparison, Signature, padj) %>%
    mutate(signif = case_when(
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ ""
    )) %>% 
    select(-padj) %>% 
    tidyr::spread(Signature, signif) %>%
    as.data.frame() %>%
    tibble::column_to_rownames(var = "Comparison") %>%
    .[row_order, col_order]
  
  return(list("heatmap_data_wide" = heatmap_data_wide,
              "col_anno" = col_anno,
              "signif_data_wide" = signif_data_wide))
  
}

#' Modified from the fgsea package
#'
#' Plots GSEA enrichment plot.
#' @param pathway Gene set to plot.
#' @param stats Gene-level statistics.
#' @param gseaParam GSEA parameter.
#' @param ticksSize width of vertical line corresponding to a gene (default: 0.2)
#' @return ggplot object with the enrichment plot.
#' @export
#' @examples
#' data(examplePathways)
#' data(exampleRanks)
#' \dontrun{
#' plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
#'                exampleRanks)
#' }
plotEnrichment2 <- function(pathway, stats,
                            gseaParam=1,
                            ticksSize=0.2, return_df = FALSE, colour, plot_num) {
  
  rnk <- rank(-stats)
  ord <- order(rnk)
  
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                          returnAllExtremes = TRUE)
  
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
  
  if (return_df) return(toPlot)
  
  diff <- (max(tops) - min(bottoms)) / 8
  
  # Getting rid of NOTEs
  x=y=NULL
  g <- rr_ggplot(toPlot, aes(x=x, y=y), plot_num = plot_num) +
    geom_point(color=colour, size=0.1) +
    geom_hline(yintercept=max(tops), colour="red", linetype="dashed") +
    geom_hline(yintercept=min(bottoms), colour="red", linetype="dashed") +
    geom_hline(yintercept=0, colour="black") +
    geom_line(color=colour) + theme_bw() +
    geom_segment(data=data.frame(x=pathway),
                 mapping=aes(x=x, y=-diff/2,
                             xend=x, yend=diff/2),
                 size=ticksSize,
                 colour="black") +
    
    theme(panel.border=element_blank(),
          panel.grid.minor=element_blank()) +
    
    labs(x="rank", y="enrichment score")
  g
}


#' Modified from the fgsea package
#'
#' Calculates GSEA statistics for a given query gene set
#'
#' Takes \emph{O(k log k)} time, where \emph{k} is a size of `selectedSize`.
#' @param stats Named numeric vector with gene-level statistics
#'  sorted in decreasing order (order is not checked).
#' @param selectedStats Indexes of selected genes in the `stats` array.
#' @param gseaParam GSEA weight parameter (0 is unweighted, suggested value is 1).
#' @param returnAllExtremes If TRUE return not only the most extreme point, but all of them. Can be used for enrichment plot
#' @param returnLeadingEdge If TRUE return also leading edge genes.
#' @return Value of GSEA statistic if both returnAllExtremes and returnLeadingEdge are FALSE.
#' Otherwise returns list with the folowing elements:
#' \itemize{
#' \item res -- value of GSEA statistic
#' \item tops -- vector of top peak values of cumulative enrichment statistic for each gene;
#' \item bottoms -- vector of bottom peak values of cumulative enrichment statistic for each gene;
#' \item leadingGene -- vector with indexes of leading edge genes that drive the enrichment, see \url{http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading}.
#' }
#' @export
#' @examples
#' data(exampleRanks)
#' data(examplePathways)
#' ranks <- sort(exampleRanks, decreasing=TRUE)
#' es <- calcGseaStat(ranks, na.omit(match(examplePathways[[1]], names(ranks))))
calcGseaStat <- function(stats, selectedStats, gseaParam=1,
                         returnAllExtremes=FALSE,
                         returnLeadingEdge=FALSE) {
  
  S <- selectedStats
  r <- stats
  p <- gseaParam
  
  S <- sort(S)
  
  m <- length(S)
  N <- length(r)
  if (m == N) {
    stop("GSEA statistic is not defined when all genes are selected")
  }
  NR <- (sum(abs(r[S])^p))
  rAdj <- abs(r[S])^p
  if (NR == 0) {
    # this is equivalent to rAdj being rep(eps, m)
    rCumSum <- seq_along(rAdj) / length(rAdj)
  } else {
    rCumSum <- cumsum(rAdj) / NR
  }
  
  
  tops <- rCumSum - (S - seq_along(S)) / (N - m)
  if (NR == 0) {
    # this is equivalent to rAdj being rep(eps, m)
    bottoms <- tops - 1 / m
  } else {
    bottoms <- tops - rAdj / NR
  }
  
  maxP <- max(tops)
  minP <- min(bottoms)
  
  if(maxP > -minP) {
    geneSetStatistic <- maxP
  } else if (maxP < -minP) {
    geneSetStatistic <- minP
  } else {
    geneSetStatistic <- 0
  }
  
  if (!returnAllExtremes && !returnLeadingEdge) {
    return(geneSetStatistic)
  }
  
  res <- list(res=geneSetStatistic)
  if (returnAllExtremes) {
    res <- c(res, list(tops=tops, bottoms=bottoms))
  }
  if (returnLeadingEdge) {
    leadingEdge <- if (maxP > -minP) {
      S[seq_along(S) <= which.max(bottoms)]
    } else if (maxP < -minP) {
      rev(S[seq_along(S) >= which.min(bottoms)])
    } else {
      NULL
    }
    
    res <- c(res, list(leadingEdge=leadingEdge))
  }
  res
}


