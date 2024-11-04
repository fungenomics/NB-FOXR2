
library(pbapply)
library(magrittr)


# ssGSEA -----------------------------------------------------------------------

#' integrate_cdfs
#' 
#' This code was adapted from GSVA and the GSEA code from the Broad Institute
#' 
#' @param gene_set Character vector of genes
#' @param gene_ranking Numeric vector specifying gene ranking
#' @param R Matrix of gene ranks
#' @param j Integer, specifying sample (column of expression matrix)
#' @param alpha Numeric, exponent giving the weight in the weighted ECDF $P^w_{in}$
integrate_cdfs <- function(gene_set, gene_ranking, expr_mat, R, j, alpha) {
  
  # Binary 1/0 vector indicating whether each gene (in ranked order)
  # is in the gene set or not
  indicator_in_gene_set <- match(rownames(expr_mat)[gene_ranking], gene_set)
  indicator_in_gene_set[!is.na(indicator_in_gene_set)] <- 1
  indicator_in_gene_set[ is.na(indicator_in_gene_set)] <- 0
  indicator_out <- 1 - indicator_in_gene_set
  
  # P_in (the weighted ECDF of genes in the set)
  # i.e. it will take a weighted step whenever a gene is in the set
  P_in <- cumsum( (abs(R[gene_ranking, j]) * indicator_in_gene_set)^alpha ) /
    sum( (abs(R[gene_ranking, j]) * indicator_in_gene_set)^alpha )
  
  # P_out (the ECDF of genes not in the set)
  # i.e. it will be length N - Ng, with a step whenever the gene is not in the set
  P_out <- cumsum( !indicator_in_gene_set ) /
    sum( !indicator_in_gene_set )
  
  # The difference in CDFs
  cdf_diffs <- P_in - P_out
  es <- sum(cdf_diffs)
  
  # Calculate running sum statistic, and leading edge subset
  n_g <- length(gene_set)
  n <- nrow(R)
  steps_in <- (abs(R[gene_ranking, j]) * indicator_in_gene_set)^alpha
  step_norm_in <- 1/sum(steps_in)
  step_norm_out <- 1/(n - n_g)
  
  running_sum <- cumsum(indicator_in_gene_set * steps_in * step_norm_in -  # Increment the score by a weighted step if a gene is a hit
                          indicator_out * step_norm_out)                   # Decrement by an fixed size step if a gene is a miss
  
  # The leading edge consists of all the genes in the gene set, which come up
  # in the ranked list prior to the max enrichment score (diff. between ECDFs)
  leading_edge <- names(running_sum)[1:which.max(running_sum)]
  
  # We only care about the ones in the gene set
  leading_edge_indicator <- 1 * ((rownames(R) %in% leading_edge) & (rownames(R) %in% gene_set))
  names(leading_edge_indicator) <- rownames(R)
  
  return(list(
    leading = leading_edge_indicator,
    es = sum(cdf_diffs),
    running = running_sum))
  
}

#' This code was adapted from GSVA
ssgsea_le <- function(expr_mat, gene_sets, alpha = 0.25, normalize = TRUE,
                      save_le = FALSE, n_cores = 1, verbose = FALSE) {
  
  require(pbapply)
  
  if (verbose) message("1. Converting gene expr to ranks...")
  
  # Convert gene expression data in X to ranks within each sample
  R <- apply(expr_mat, 2, function(x) as.integer(rank(x, na.last = TRUE)))
  rownames(R) <- rownames(expr_mat)
  
  # # Collect the differences between CDFs
  # diff_list <- vector("list", length = n)
  
  if (verbose) message(glue("2. Calculating enrichment, parallelizing over ", n_cores, " cores..."))
  
  # For each sample S_1 to S_n, calculate the score for each gene set
  # Parallelize over samples
  ssgsea_out <- pblapply(1:ncol(expr_mat), function(j) {
    
    # The order (by indices) by which to retrieve genes from the matrix
    # to get them from high -> low
    
    # BUG: I think that subsetting the dataframe in this way does not properly
    # deal with any NAs for that sample
    gene_ranking <- order(R[, j], decreasing = TRUE, na.last = TRUE)
    sample_out <- lapply(gene_sets, integrate_cdfs, gene_ranking, expr_mat, R, j, alpha)
    
    return(sample_out)
    
  }, cl = n_cores)
  
  # TODO: this is super unwieldy and unneeded!
  # Would be better to simply return the genes in the set which are in the leading
  # edge as a char vector, that could be saved.
  # But, that may not be the solution, since the leading edge of interest,
  # in our case, is defined per-sample, and so we'll only really care about it
  # with respect to biological groups of samples
  filter_zeros <- function(df) {
    
    df[rowSums(df) != 0, ]
    
  }
  
  if (verbose) message("3. Tidying ouput...")
  
  if (save_le) {
    
    # Tidy the binary data indicating if genes are within the leading edge subset
    leading_mat <- lapply(ssgsea_out,
                          function(sample_out) purrr::transpose(sample_out)$leading) %>% 
      # we have one element per sample, and each has one element per signature
      magrittr::set_names(colnames(expr_mat)) %>%
      # now, flip to get one element per signature
      purrr::transpose() %>%
      lapply(bind_rows) %>% 
      lapply(as.data.frame) %>% 
      lapply(magrittr::set_rownames, colnames(expr_mat)) %>% 
      lapply(filter_zeros) %>% 
      lapply(t)
    
    # Tidy the running sum data
    # Can't coerce to a dataframe because then order would be lost
    running_sum <- lapply(ssgsea_out,
                          function(sample_out) purrr::transpose(sample_out)$running) %>% 
      magrittr::set_names(colnames(expr_mat)) %>%
      purrr::transpose()
    
  } else {
    
    leading_mat <- NULL
    
  }
  
  # Tidy the enrichment scores for each set for each sample
  es <- lapply(ssgsea_out, function(sample_out) unlist(purrr::map(sample_out, "es"))) %>%
    data.frame %>%
    as.matrix()
  
  if (normalize) {
    
    # Normalize the whole dataset
    normES <- function(mat) apply(mat, 2, function(x, mat) x /
                                    (range(mat)[2] - range(mat)[1]), mat)
    
    es <- normES(es)
    
  }
  
  rownames(es) <- names(gene_sets)
  colnames(es) <- colnames(expr_mat)
  
  return(list("enrichment_scores" = es,
              "leading_edge" = leading_mat))
  
  if (verbose) message("4. Done.")
  
}


compute_scores_bin <- function(expr_mat,
                               gene_sets,
                               normalize = FALSE,
                               alpha = 0.75,
                               out_prefix = NULL,
                               bin_size = 50,
                               n_cores = 1,
                               save_le = FALSE) {
  
  # Generating the binned data
  if (ncol(expr_mat) <= bin_size) {
    binned_data <- list(expr_mat)
    message("@ skipping binning...")
  } else {
    
    binned_data <- bin_counts(counts   = expr_mat,
                              bin_size = bin_size)
    
  }
  
  
  # Computing ssGSEA for all elements
  binned_results <- lapply(1:length(binned_data), function(i) {
    
    message("@@ bin ", i, "/", length(binned_data), "...")
    
    out <- ssgsea_le(expr_mat  = binned_data[[i]],
                     gene_sets = gene_sets,
                     alpha     = alpha,
                     normalize = normalize,
                     save_le   = save_le,
                     n_cores   = n_cores)
    return(out)
  })
  
  # Putting them back together
  out_matrix <- do.call(cbind,
                        lapply(binned_results, function(res) {res$enrichment_scores}))
  
  # Tidy scores and write to file
  out_tidy <- out_matrix %>%
    as.data.frame() %>% 
    magrittr::set_rownames(NULL) %>% 
    magrittr::set_rownames(names(gene_sets)) %>% 
    tibble::rownames_to_column(var = "Signature")
  
  return(out_tidy)
  
}


#' @param genelists A named list, where each element is a character vector of
#' genes. Expects gene symbols.
#' @param summary_type Character specifying how to treat single cells for 
#' the enrichment. If "cluster_mean", computes cluster mean expression and
#' evaluates enrichment per cluster. If "metacells", assumes the seurat object
#' contains metacells (by \code{seurat2metacells}) and evaluates enrichment
#' per metacell.
#' 
#' @return A dataframe, with the first column "Cell" (where the order matches
#' the order in \code{seurat}), and following columns corresponding to signatures,
#' named as in the \code{genelists} input.
compute_scores_sc <- function(seurat,
                              genelists,
                              summary_type = "single_cells",
                              out_prefix = NULL,
                              bin_size = 100,
                              n_cores = 1,
                              save_le = FALSE) {
  
  # Generating the binned data
  binned_data <- bin_seurat(seurat = seurat,
                            bin_size = bin_size)
  
  # Computing ssGSEA for all elements
  binned_results <- lapply(1:length(binned_data), function(i) {
    
    message("@@ bin ", i, "/", length(binned_data), "...")
    
    out <- ssgsea_le(expr_mat  = binned_data[[i]],
                     gene_sets = genelists,
                     alpha     = 0.75,
                     normalize = FALSE,
                     save_le   = save_le,
                     n_cores   = n_cores)
    return(out)
  })
  
  # Putting them back together
  out_matrix <- do.call(cbind,
                        lapply(binned_results, function(my_ssGSEA_results) {my_ssGSEA_results$enrichment_scores}))
  
  # Tidy scores and write to file
  out_tidy <- out_matrix %>%
    data.frame %>% 
    magrittr::set_rownames(NULL) %>% 
    magrittr::set_rownames(names(genelists)) %>% 
    tibble::rownames_to_column(var = "Signature")
  
  if (!is.null(out_prefix)) write.table(file = glue("{out_prefix}_{summary_type}.ssgsea_scores.txt"),
                                        x = .,
                                        sep = "\t",
                                        quote = TRUE,
                                        row.names = FALSE)
  else return(out_tidy)
  
}


bin_counts <- function(counts,
                       bin_size = 50) {
  
  # Getting the amount of bins
  n_bins <- ceiling(dim(counts)[2]/bin_size)
  
  # Getting the full bins
  binned_list <- lapply(1:(n_bins-1), function(i) {
    
    lower_index <- (i-1)*bin_size + 1
    upper_index <- i*bin_size
    
    return(counts[, lower_index:upper_index])
  })
  
  # Getting the last bin
  lower_index <- (n_bins-1)*bin_size + 1
  upper_index <- dim(counts)[2]
  
  # If they are too close - could be less than 10% difference, 
  # we need to add the cells to the previous bin
  # Otherwise we can just complete the last bin
  if(upper_index - lower_index < bin_size/10) {
    binned_list[[n_bins-1]] <- cbind(binned_list[[n_bins-1]], counts[,lower_index:upper_index])
  } else {
    binned_list[[n_bins]] <- counts[,lower_index:upper_index]
  }
  
  return(binned_list)
  
}


#' For V3 Seurat objects
#' From Marie Coutelier
#' 
#' @param bin_size Numeric, number of cells for each bin. Default: 100
bin_seurat <- function(seurat,
                       bin_size = 100) {
  
  # Get the expression matrix from the Seurat object
  counts <- as.matrix(GetAssayData(object = seurat, slot = "data"))
  
  # Getting the amount of bins
  n_bins <- ceiling(dim(counts)[2]/bin_size)
  
  # Getting the full bins
  binned_list <- lapply(1:(n_bins-1), function(i) {
    
    lower_index <- (i-1)*bin_size + 1
    upper_index <- i*bin_size
    
    return(counts[, lower_index:upper_index])
  })
  
  # Getting the last bin
  lower_index <- (n_bins-1)*bin_size + 1
  upper_index <- dim(counts)[2]
  
  # If they are too close - could be less than 10% difference, 
  # we need to add the cells to the previous bin
  # Otherwise we can just complete the last bin
  if(upper_index - lower_index < bin_size/10) {
    binned_list[[n_bins-1]] <- cbind(binned_list[[n_bins-1]], counts[,lower_index:upper_index])
  } else {
    binned_list[[n_bins]] <- counts[,lower_index:upper_index]
  }
  
  return(binned_list)
  
}

#' Helper function to bin data matrix and run ssGSEA on scRNAseq seurat object
#' @param bin_by cluster or discrete column in metadata
#'
#' @export
#' @author Nisha Kabir, Selin Jessa
#' NOTE: Change by Bhavyaa read_tsv -> data.table::fread() for speed and bug fix

bin_ssgsea <- function(seurat, 
                       bin_by = "cluster",
                       summary_type = "single_cells",
                       out_prefix, 
                       genelist,
                       save_leading_edge = FALSE) { 
    
    require(data.table)
    
    if (bin_by == "cluster") {
        
        to_bin <- levels(Seurat::Idents(seurat))
        
        for (i in seq_along(to_bin)) {
            
            seurat_sub <- subset(seurat, idents = to_bin[i])
            
            if (summary_type == "cluster_mean") counts <- mean_cluster_expression(seurat_sub)
            else counts <- as.matrix(GetAssayData(seurat_sub))
            
            results <- ssgsea_le(expr_mat = counts, 
                                 gene_sets = genelist,
                                 alpha     = 0.75,
                                 normalize = FALSE,
                                 save_le   = save_leading_edge)
            results$enrichment_scores %>% 
                data.frame %>% 
                magrittr::set_rownames(NULL) %>% 
                magrittr::set_rownames(rownames(results$enrichment_scores)) %>% 
                magrittr::set_colnames(colnames(counts)) %>% 
                tibble::rownames_to_column(var = "Signature") %>%
                {write.table(file = glue("{out_prefix}.subset{i}.ssgsea_scores.txt"),
                             x = .,
                             sep = "\t",
                             quote = TRUE,
                             row.names = FALSE)}
            
            if (save_leading_edge) {
                
                leading_edge <- results$leading_edge
                save(leading_edge, file = glue("{out_prefix}.subset{i}.leading_edge.Rda"))
                
            }
            
        }
        
    } else {
        
        to_bin <- unique(seurat@meta.data[[bin_by]])
        
        for (i in seq_along(to_bin)) {
            
            cells_use <- seurat@meta.data %>% 
                rownames_to_column(var = "cell") %>% 
                filter(!!ensym(bin_by) == to_bin[i]) %>%
                pull(cell)
            seurat_sub <- seurat[ , cells_use]
            counts <- as.matrix(GetAssayData(seurat_sub))
            results <- ssgsea_le(expr_mat = counts, 
                                 gene_sets = genelist,
                                 alpha     = 0.75,
                                 normalize = FALSE,
                                 save_le   = FALSE)
            results$enrichment_scores %>% 
                data.frame %>% 
                magrittr::set_rownames(NULL) %>% 
                magrittr::set_rownames(rownames(results$enrichment_scores)) %>% 
                magrittr::set_colnames(colnames(counts)) %>% 
                tibble::rownames_to_column(var = "Signature") %>%
                {write.table(file = glue("{out_prefix}.subset{i}.ssgsea_scores.txt"),
                             x = .,
                             sep = "\t",
                             quote = TRUE,
                             row.names = FALSE)}
            
        }
        
    }
    
    # Gather output
    separate_files <- list.files(dirname(out_prefix),
                                 pattern = paste0(basename(out_prefix), ".subset.*.ssgsea_scores.txt"), full.names = TRUE)
    
    purrr::map(separate_files,
               ~ .x %>% data.table::fread() %>% as.data.frame %>% tibble::column_to_rownames(var = "Signature")) %>% 
        {do.call(cbind, .)} %>%
        as.data.frame %>%
        tibble::rownames_to_column(var = "Signature") %>%
        {write.table(file = glue("{out_prefix}.{summary_type}.ssgsea_scores.txt"),
                     x = .,
                     sep = "\t",
                     quote = TRUE,
                     row.names = FALSE)}
    
    file.remove(separate_files)
    
    if (save_leading_edge) {
        
        # Gather output
        separate_files <- list.files(dirname(out_prefix),
                                     pattern = paste0(basename(out_prefix), ".subset.*.leading_edge.Rda"), full.names = TRUE)
        
        leading_edge_list <- list()
        for (i in 1:length(genelist)) {
            
            leading_edge_list[[i]] <- map_dfc(separate_files, ~ get(load(.x)) %>% .[names(genelist)[i]] %>% unname() %>% as.data.frame())
            
        }
        
        names(leading_edge_list) <- names(genelist)
        
        save(leading_edge_list, file = glue("{out_prefix}.leading_edge.Rda"))
        
        file.remove(separate_files)
        
    }
    
}


# Leading edge 
# NOTE: I commented out most of the leading-edge related code in order to remove the cytobox dependency

#' 
#' #' Get leading edge genes
#' #'
#' #' For a given signature (cluster), return the genes which are
#' #' in the leading edge subset in some proportion of the samples. Since more, rather
#' #' than fewer, genes in the gene set tend to be in the leading edge, the default
#' #' proportion is 1, i.e. the most stringent.
#' #' 
#' #' @param leading_edge The leading edge list object produced by \code{compute_scores_le},
#' #' either a character vector with the path, or the list itself.
#' get_leading_edge <- function(leading_edge,
#'                              cluster,
#'                              samples = NULL,
#'                              proportion = 1,
#'                              return_mat = FALSE,
#'                              return_symbols = TRUE,
#'                              leading_edge_only = TRUE,
#'                              species = "hg") {
#'   
#'   if (is.character(leading_edge)) load(leading_edge)
#'   
#'   if (!(cluster %in% names(leading_edge))) stop("Cluster not found in signatures...")
#'   
#'   if (is.null(samples)) le_clust <- leading_edge[[cluster]]
#'   else if (length(samples) == 1) {
#'     
#'     # Simply return the leading edge genes, by definition, for one sample  
#'     return(rownames(leading_edge[[cluster]])[which(leading_edge[[cluster]][, samples] == 1)])
#'     
#'   } else {
#'     
#'     # Subset to the leading edge info for all relevant samples
#'     le_clust <- leading_edge[[cluster]][, samples]
#'     message("...Identifying leading edge genes for ", ncol(le_clust), " samples")
#'     
#'   }
#'   
#'   if (return_mat) return(le_clust)
#'   
#'   if (leading_edge_only) le_genes <- le_clust[rowSums(le_clust) >= proportion * ncol(le_clust), ] %>%
#'       rownames()
#'   else le_genes <- le_clust %>% rownames()
#'   
#'   if (return_symbols & grepl("^ENS", head(le_genes[[1]]))) return(ensembl2symbols(le_genes))
#'   else return(le_genes)
#'   
#'   
#' }
#' 

# 
# ensembl2symbols_safe <- function(genes, ...) {
#   
#   require(icytobox)
#   
#   sym <- lapply(genes, ensembl2symbols, ...)
#   empty <- which(unlist(map(sym, ~ identical(unlist(.x), character(0)))))
#   sym[empty] <- genes[empty]
#   
#   return(unlist(sym))
#   
# }

#' #' Rank the leading edge genes
#' #' 
#' #' For a given set of leading edge genes, rank them by their median rank expression
#' #' across a given sample group.
#' #' 
#' #' @param leading_edge The leading edge list object produced by \code{compute_scores_le},
#' #' either a character vector with the path, or the list itself.
#' #' @param leading_edge_only Logical, if TRUE, compute the ranks only for genes in 
#' #' the leading edge subset. Otherwise, compute ranks for all genes in the signature.
#' #' Default: TRUE
#' #' @param return_symbols Logical, if TRUE, return the gene symbols, otherwise return
#' #' the Ensembl IDs. Default: TRUE
#' #' @param relative_rank Logical, if TRUE, return ranks starting from 1, 2, 3... etc, 
#' #' describing how highly expressed these genes are across the samples, relative to each other.
#' #' If FALSE, return the actual median ranks of the genes across all samples. Default: FALSE
#' rank_leading_edge <- function(leading_edge = NULL, counts, cluster, gene_set = NULL,
#'                               proportion = 1,
#'                               samples = NULL,
#'                               leading_edge_only = TRUE,
#'                               return_symbols = TRUE,
#'                               relative_rank = FALSE,
#'                               return_raw_ranks = FALSE,
#'                               return_gene_set = FALSE,
#'                               debug = FALSE,
#'                               verbose = FALSE,
#'                               species = "hg") {
#'   
#'   expr_mat <- counts
#'   
#'   if (!is.null(leading_edge)) {
#'     
#'     if (verbose) message("1/3: Getting leading edge...")
#'     if (leading_edge_only) gene_set <- get_leading_edge(leading_edge,
#'                                                         samples = samples,
#'                                                         cluster = cluster,
#'                                                         proportion = proportion,
#'                                                         return_symbols = FALSE)
#'     else gene_set <- get_leading_edge(leading_edge,
#'                                       samples = samples,
#'                                       cluster = cluster,
#'                                       proportion = proportion,
#'                                       return_symbols = FALSE,
#'                                       leading_edge_only = FALSE,
#'                                       species = species)
#'     
#'     
#'     if (debug) return(gene_set)
#'     
#'   } else {
#'     
#'     if (verbose) message("1/3: Using provided gene set...")
#'     
#'   }
#'   
#'   if (verbose) message("2/3: Computing ranks...")
#'   # Convert gene expression data in expr_mat to ranks within each sample
#'   # By negating the expression values, we will rank the genes from high expression
#'   # to low, i.e. highly-expressed genes will have small ranks
#'   R <- apply(expr_mat, 2, function(x) as.integer(rank(-x)))
#'   rownames(R) <- rownames(expr_mat)
#'   R_gene_set <- R[gene_set, ]
#'   
#'   if (!is.null(samples) & length(samples) == 1) {
#'     df <- data.frame(R_gene_set[, samples]) %>%
#'       setNames("rank") %>% 
#'       tibble::rownames_to_column(var = "gene") %>% 
#'       arrange(rank) %>% 
#'       as.data.frame %>% 
#'       tibble::column_to_rownames(var = "gene")
#'     return(df)
#'   }
#'   else if (!is.null(samples) & length(samples) > 1) R_gene_set <- R_gene_set[, samples]
#'   
#'   if (return_gene_set) return(R_gene_set)
#'   
#'   if (verbose) message("3/3: Computing medians...")
#'   # Calculate the median rank of each gene across the relevant samples
#'   
#'   if (verbose) message("...Computing median ranks across ", ncol(R_gene_set), " samples")
#'   
#'   if (!return_raw_ranks) {
#'     
#'     if (return_symbols) ranks <- data.frame(gene = ensembl2symbols_safe(rownames(R_gene_set), sp = species),
#'                                             median_rank = matrixStats::rowMedians(R_gene_set))
#'     else ranks <- data.frame(gene = rownames(R_gene_set),
#'                              median_rank = matrixStats::rowMedians(R_gene_set))
#'     
#'     # Make the ranks relative, to start from 1, 2, 3, etc.
#'     if (relative_rank) ranks <- ranks %>% mutate(median_rank = min_rank(median_rank))
#'     
#'   } else {
#'     
#'     if (return_symbols) {
#'       
#'       ranks <- as.data.frame(R_gene_set) %>% 
#'         tibble::add_column(gene = ensembl2symbols_safe(rownames(R_gene_set), sp = species), .before = 1) %>% 
#'         tibble::add_column(median_rank = matrixStats::rowMedians(R_gene_set), .after = "gene")
#'       
#'     } else {
#'       
#'       ranks <- as.data.frame(R_gene_set) %>% 
#'         tibble::add_column(gene = rownames(R_gene_set), .before = 1) %>% 
#'         tibble::add_column(median_rank = matrixStats::rowMedians(R_gene_set), .after = "gene")
#'       
#'     }
#'     
#'   }
#'   
#'   return(arrange(ranks, median_rank))
#'   
#' }
#' 



#' #' TODO:
#' #' - Informative title with cluster name, sample group, enrichment score
#' #' - Make this a dashboard like representation which has a second panel with 
#' #' the gene stats within the signature
#' #' - The idea is that we would have all the information in one place, to see
#' #' how important a gene is for a signature vs. for its enrichment
#' #' 
#' #' @param le_df Data frame as returned by rank_leading_edge with return_raw_ranks = TRUE
#' plot_le <- function(le_df, title = "", colour = "gray50", n_top = NULL, gene_order = NULL,
#'                     outlier_size = -1, gene_font_size = 8) {
#'   
#'   if (ncol(le_df) == 2) stop("Did you run rank_leading_edge with 'return_raw_ranks = TRUE' ?")
#'   
#'   # Get the top N genes with smallest ranks
#'   if (!is.null(n_top)) le_df <- le_df %>% top_n(-n_top, median_rank)
#'   
#'   le_df_long <- le_df %>%
#'     gather(sample, rank, 2:ncol(.)) %>% 
#'     mutate(sample = factor(sample, levels = c("median_rank", names(le_df)[3:ncol(le_df)])))
#'   
#'   # If a gene order is not provided, then order by median rank as in input df
#'   if (is.null(gene_order)) {
#'     le_df_long <- le_df_long %>% mutate(gene = factor(gene, levels = rev(le_df$gene)))
#'   } else {
#'     le_df_long <- le_df_long %>% mutate(gene = factor(gene, levels = rev(gene_order)))
#'   }
#'   
#'   # Controlling point order with a second points layer, following the solution at
#'   # https://stackoverflow.com/questions/15706281/controlling-order-of-points-in-ggplot2-in-r
#'   p <- le_df_long %>% 
#'     ggplot(aes(x = gene, y = rank)) +
#'     geom_boxplot(width = 0.8, fill = colour, outlier.size = outlier_size) +
#'     ylab("Rank by expression") +
#'     xlab("Leading edge gene") +
#'     ggtitle(title) +
#'     coord_flip() +
#'     no_legend() +
#'     theme(axis.text.y = element_text(size = gene_font_size))
#'   
#'   return(p)
#'   
#' }




# Plotting helpers -------------------------------------------------------------

plot_dim_red <- function(dr_df, method,
                         palette = NULL,
                         title = NULL,
                         variance = NULL,
                         label = FALSE,
                         samples_to_label = NULL,
                         point_size = 1,
                         alpha = 1,
                         shape = NULL,
                         palette_shape = NULL) {
  
  p1 <- dr_df %>% 
    ggplot(aes(x = Dim1, y = Dim2))
  
  if (!is.null(shape) & !is.null(palette_shape)) {
    
    p1 <- p1 +
      geom_point(aes_string(fill = "Group", colour = "Group", shape = shape),
                 size = point_size, alpha = alpha) +
      scale_shape_manual(values = palette_shape)
    
  } else {
  
    p1 <- p1 +
      geom_point(aes(fill = Group, colour = Group),
                 size = point_size, alpha = alpha)
    
    
  }
  
    
  if (!is.null(palette)) p1 <- p1 + scale_colour_manual(values = palette)
  
  p1 <- p1 + theme_min()
  
  # Axes labels
  axes_labels <- case_when(method == "pca"  ~ c(paste0("PC1 (", variance[1], "%)"), paste0("PC2 (", variance[2], "%)")),
                           method == "tsne" ~ c("tSNE 1", "tSNE 2"))
  
  p1 <- p1 + xlab(axes_labels[1]) + ylab(axes_labels[2]) + no_ticks()
  
  # Sample labels
  if (!is.null(samples_to_label)) dr_df_lab <- dr_df %>% filter(Sample %in% samples_to_label)
  else dr_df_lab <- dr_df
  
  if (label) {
    
    require(ggrepel)
    
    p1 <- p1 + geom_text_repel(data = dr_df_lab, aes(x = Dim1, y = Dim2, label = Sample),
                               size = 3,
                               colour = "gray50",
                               segment.colour = "gray70")
    
  }
  
  if (!is.null(title)) p1 <- p1 + ggtitle(title)
  
  return(p1)
  
}


ssgsea_dim_red <- function(scores,
                           method = "pca",
                           colour_var,
                           pal = NULL,
                           return_df = FALSE,
                           return_prcomp = FALSE,
                           title = NULL,
                           samples_to_label = NULL,
                           label = ifelse(is.null(samples_to_label), FALSE, TRUE),
                           point_size = 2,
                           ...) {
  
  score_mat <- t(scores)
  
  dr <- list()
  
  if (method == "tsne") {
    
    require(Rtsne)
    
    tsne_out  <- Rtsne(score_mat, ...)
    dr$coords <- tsne_out$Y
    variance  <- NULL
    
  } else if (method == "pca") {
    
    pca       <- prcomp(score_mat)
    variance  <- round(pca$sdev^2 / sum(pca$sdev^2), 3) * 100
    dr$coords <- pca$x
    
    print(variance[1:2])
    
    if (return_prcomp) return(pca)
    
  }
  
  dr_df <- data.frame(Dim1 = dr$coords[, 1],
                      Dim2 = dr$coords[, 2],
                      Sample = rownames(score_mat))
  
  if (return_df) {
    message("@ returning data...")
    return(dr_df)
  }
  
  # dr_plot <- plotDimRed(dr_df, method, colour_var, pal, title, variance, label, samples_to_label, point_size)
  
  # return(dr_plot)
  
}

