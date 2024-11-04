
# Plotting & utility functions adapted from cytobox & from Samantha Worme

# Plotting functions -----------------------------------------------------------

#' Colour cells in reduced dimension space by gene expression
#'
#' Plot a low-dimensional embedding of the cells,
#' coloured by expression of a gene, or mean expression of a group of marker
#' genes. Defaults to UMAP space, but see the \code{reduction} argument for
#' how to plot in t-SNE or PCA space instead. This function is based on 
#' \code{Seurat::FeaturePlot} and cytobox::tsneByMeanMarkerExpression.
#'
#' @param seurat Seurat object, where dimensionality reduction has been applied,
#' i.e. (after applying Seurat::RunPCA(), Seurat::RunTSNE(), or Seurat::RunUMAP() to the object)
#' @param genes String or character vector specifying gene(s) to use
#' @param reduction String specifying the dimensionality reduction to use,
#' retrieves UMAP by default.
#' Default: "umap"
#' @param label Logical, whether to label clusters on the plot. Default: TRUE.
#' @param label_repel Logical, if \code{label} is TRUE, whether to plot cluster
#' labels repelled from the center, on a slightly transparent white background and
#' with an arrow pointing to the cluster center. If FALSE, simply plot the
#' cluster label at the cluster center. Default: TRUE.
#' @param label_size Numeric, controls the size of text labels. Default: 4.
#' @param palette String or character vector. If a string,
#' one of "viridis", "blues", or "redgrey", specifying which gradient
#' palette to use. Otherwise, a character vector of colours (from low to high)
#' to interpolate to create the scale. Default: redgrey.
#' @param point_size Numeric, size of points in scatterplot. Default: 1. (A smaller
#' value around 0.5 is better for plots which will be viewed at small scale.)
#' @param alpha Numeric, fixed alpha for points. Default: 0.6
#' @param legend Logical, whether or not to plot legend. Default: TRUE
#' @param hide_ticks Logical, whether to hide axis ticks. Default: FALSE
#' @param hide_axes Logical, whether to hide axis labels. Default: FALSE
#' @param title (Optional) String specifying the plot title
#' @param limits (Optional) A numeric vector of length two providing the limits to
#' use for the colour scale (documentation from \code{\link[ggplot2]{continous_scale}}. 
#' Default: 0 and max of the data.
#' @param label_short (Optional/Experimental!!) Logical, if TRUE, assumes clusters
#' (at levels(Ident(seurat))) consist of a prefix and a suffix separated by a non-alpha
#' numeric character (\code{"[^[:alnum:]]+"}), and tries to separate these names
#' and only plot the prefix, for shorter labels and a cleaner plot. Default: FALSE.
#' @param dim1 Numeric, index of dimension from \code{reduction} to plot on
#' the x-axis. e.g. to plot the 3rd PC on the x-axis, pass 3. Default: 1.
#' @param dim2 Numeric, like \code{dim2}, but for the y-axis. Default: 2.
#' @param return_df Logical, whether to return the mean expression dataframe.
#'
#' @export
#' @return A ggplot object
#'
#' @author Selin Jessa, Samantha Worme
plot_mean_expression <- function(seurat, genes,
                                 assay = "RNA",
                                 reduction = "umap",
                                 label = TRUE,
                                 label_repel = TRUE,
                                 label_size = 4,
                                 palette = "redgrey",
                                 point_size = 1,
                                 alpha = 0.6,
                                 legend = TRUE,
                                 hide_ticks = FALSE,
                                 hide_axes = FALSE,
                                 title = NULL,
                                 limits = NULL,
                                 label_short = FALSE,
                                 dim1 = 1,
                                 dim2 = 2,
                                 return_df = FALSE) {
    
    exp_df <- suppressWarnings(mean_gene_expression(seurat, genes, assay))
    
    # Get dimensionality reduction coordinates
    exp_df <- get_embedding(seurat, assay, exp_df, reduction, dim1, dim2) %>%
        # Order in which points will be plot, "front" points at the bottom
        dplyr::arrange(Mean_marker_expression)
    
    if (return_df) return(exp_df)
    
    # Get the variable names
    vars <- colnames(seurat[[reduction]]@cell.embeddings)[c(dim1, dim2)]
    
    # Set limits: if not provided, use default min/max
    if (is.null(limits)) limits <- c(NA, NA)
    
    # Plot
    gg <- exp_df %>%
        ggplot(aes(x = .data[[vars[1]]], y = .data[[vars[2]]])) +
        geom_point(aes(colour = Mean_marker_expression), size = point_size, alpha = alpha)
    
    if (length(palette) == 1) {
        
        if (palette == "viridis") {
            
            gg <- gg + viridis::scale_color_viridis(limits = limits)
            
        } else if (palette == "blues") {
            
            gg <- gg + scale_colour_gradientn(
                colours = RColorBrewer::brewer.pal(n = 8, name = "Blues"),
                limits = limits)
            
        } else if (palette == "redgrey") {
            
            # NOTE: palette chosen is not the default gradient from gray -> red
            # but sets a midpoint at a lighter colour
            gg <- gg + scale_color_gradientn(
                colours = grDevices::colorRampPalette(c("gray83", "#E09797", "red"))(n = 200),
                limits = limits)
            
        } else {
            
            stop("Please pass the palette as a character vector ",
                 "or specify one of: viridis, blues, redgrey")
            
        }
        
    } else if (length(palette) == 2) {
        
        gg <- gg + scale_color_gradient(low = palette[1], high = palette[2], limits = limits)
        
    } else {
        
        gg <- gg + scale_color_gradientn(colours = palette, limits = limits)
        
    }
    
    
    if (label) {
        
        if (reduction == "pca") {
            
            message("Plotting labels is currently only available for reduction = 'tsne';",
                    " returning plot without labels.")
            
        } else {
            
            centers <- suppressMessages(compute_cluster_centers(seurat, reduction = reduction))
            gg <- gg + add_labels(centers, label_repel, label_size, label_short)
            
        }
    }
    
    axes <- gsub("_", " ", vars)
    
    gg <- gg +
        xlab(axes[1]) +
        ylab(axes[2]) +
        labs(colour = case_when(assay == "RNA"                            ~ "Expression",
                                assay == "SCENIC"                         ~ "Regulon AUC",
                                assay == "chromVAR"                       ~ "Motif activity",
                                assay %in% c("peaks", "ACTIVITY", "ATAC", "promoters") ~ "Accessibility",
                                TRUE                                      ~ "Score")) +
        theme_min2()
    
    if (!is.null(title)) gg <- gg + ggtitle(title)
    if (!legend) gg <- gg + hide_legend()
    if (hide_ticks) gg <- gg + hide_ticks()
    if (hide_axes) gg <- gg + xlab(NULL) + ylab(NULL)
    
    return(gg)
    
}



#' @describeIn tsneByMeanMarkerExpression Shortcut function for plotting mean expression
#' @export
feature <- function(seurat, genes,
                    assay = "RNA",
                    per_gene = TRUE,
                    label = TRUE,
                    palette = "redgrey",
                    label_repel = FALSE,
                    label_size = 4,
                    label_short = FALSE,
                    legend = FALSE,
                    title = NULL,
                    reduction = "umap",
                    limits = c(NA, NA),
                    dim1 = 1,
                    dim2 = 2,
                    alpha = 0.6,
                    point_size = 0.5,
                    ncol = ifelse(length(genes) == 1, 1, ifelse(length(genes) %in% c(2, 4), 2, 3)),
                    hide_ticks = TRUE,
                    hide_axes = TRUE,
                    combine = TRUE) {
    
    if ((length(genes) >= 20) & per_gene) message("NOTE: you have input a lot of genes! ",
                                                  "This function by default generates ",
                                                  "one plot per gene. Set per_gene = FALSE ",
                                                  "to plot a summary statistic of all genes.")
    
    
    if (per_gene) {
        
        genes_out <- find_genes(seurat, genes, assay = assay)
        if (length(genes_out$undetected > 0)) message(paste0("NOTE: [",
                                                             paste0(genes_out$undetected, collapse = ", "),
                                                             "] undetected in the data"))
        
        if(length(genes_out$detected) == 0) stop("No genes specified were ",
                                                 "found in the data.")
        
        if ((ncol == 3) & (length(genes_out$detected) < 3)) ncol <- 2
        
        # this is a non-elegant way to ensure the supplied title is used in the one-gene case...
        if (!is.null(title) & length(genes) == 1) titles <- title
        else titles <- genes_out$detected
        
        plots <- lapply(seq_along(genes_out$detected),
                        function(i) plot_mean_expression(seurat,
                                                         genes_out$detected[i],
                                                         assay = assay,
                                                         reduction = reduction,
                                                         palette = palette,
                                                         title = titles[i],
                                                         legend = legend,
                                                         label = label,
                                                         label_short = label_short,
                                                         label_repel = label_repel,
                                                         label_size = label_size,
                                                         hide_ticks = hide_ticks,
                                                         hide_axes = hide_axes,
                                                         point_size = point_size,
                                                         limits = limits,
                                                         dim1 = dim1,
                                                         dim2 = dim2))
        
        if (combine) cowplot::plot_grid(plotlist = plots, ncol = ncol)
        else return(plots)
        
    } else {
        
        plot_mean_expression(seurat,
                             genes,
                             assay = assay,
                             reduction = reduction,
                             palette = palette,
                             title = title,
                             legend = legend,
                             label = label,
                             label_repel = label_repel,
                             label_size = label_size,
                             label_short = label_short,
                             hide_ticks = hide_ticks,
                             hide_axes = hide_axes,
                             point_size = point_size,
                             limits = limits,
                             dim1 = dim1,
                             dim2 = dim2)
    }
}


plot_dr <- function (seurat,
                     reduction         = "umap",
                     color_by          = NULL,
                     colors            = NULL,
                     color_by_type     = "discrete",
                     label             = TRUE,
                     label_repel       = TRUE,
                     label_size        = 4,
                     # Decide default point size based on density / how many points are plot
                     point_size = ifelse(length(colnames(seurat)) > 300, 0.6, 1.3),
                     alpha             = 0.8,
                     # If coloring by clusters, hide legend by default, otherwise, show it
                     legend = ifelse((is.null(color_by)) && (label), FALSE, TRUE),
                     cells             = NULL,
                     show_all_cells    = TRUE,
                     order_by          = NULL,
                     clusters_to_label = NULL,
                     hide_ticks        = TRUE,
                     title             = NULL,
                     label_prefix_only = FALSE,
                     label_suffix_only = FALSE,
                     sep               = NULL,
                     na_color          = "gray80",
                     limits            = NULL,
                     constrain_scale   = TRUE,
                     hide_axes         = FALSE,
                     dim1              = 1,
                     dim2              = 2,
                     border_color      = NULL,
                     border_size       = NULL,
                     rasterize         = FALSE) {
    
    # Error if the selected reduction has not yet been computed
    if (!(reduction %in% names(seurat@reductions))) stop(reduction, " reduction has not been computed.")
    
    # Create a dataframe holding the embedding
    embedding <- data.frame(Cell = colnames(seurat),
                            dim1 = seurat[[reduction]]@cell.embeddings[, dim1],
                            dim2 = seurat[[reduction]]@cell.embeddings[, dim2],
                            Cluster = Idents(seurat), stringsAsFactors = FALSE)
    
    # If required, order the cells in the z-axis by the specified variable
    if (!is.null(order_by)) {
        
        # We can only order by a nuemric variable, or a discrete variable with specified order
        # of levels (i.e. a factor)
        if (!is.numeric(seurat@meta.data[[order_by]]) && !is.factor(seurat@meta.data[[order_by]])) {
            
            stop("The variable specified in 'order_by' is neither numeric ",
                 "nor a factor. If the column is of type character, consider ",
                 "converting it to a factor. Otherwise, pass the name of a numeric column.")
            
        }
        
        embedding[[order_by]] <- seurat@meta.data[[order_by]]
        
        
    } else if ((!is.null(color_by)) && is.numeric(seurat@meta.data[[color_by]])) {
        # If it's not specified, by default, order by the same variable used
        # to color cells
        
        order_by <- color_by
        
    }
    
    if (!is.null(color_by)) embedding[[color_by]] <- seurat@meta.data[[color_by]]
    
    # Highlight certain cells in the plot
    if (!is.null(cells)) {
        
        # Show all cells, but only color the selected ones, by setting the color_by
        # value to NA for all other cells
        if (show_all_cells) {
            
            if (is.null(color_by)) embedding[!(embedding$Cell %in% cells), ]$Cluster <- NA
            else embedding[!(embedding$Cell %in% cells), ][[color_by]] <- NA
            
        } else { # Otherwise, only display and color the selected ones
            
            embedding <- embedding %>% filter(Cell %in% cells)
            
        }
    }
    
    # We want to sort point such that any NAs will be plot first/underneath
    color_by2 <- ifelse(is.null(color_by), "Cluster", color_by)
    
    # If the ordering variable is specified, sort NA cells at the bottom, and then order by that variable
    if (!is.null(order_by)) embedding <- embedding %>% arrange_(paste0("!is.na(", color_by2, ")"), order_by)
    else embedding <- embedding %>% arrange_(paste0("!is.na(", color_by2, ")"))
    
    # Start the plot!
    gg <- ggplot(embedding, aes(x = dim1, y = dim2))
    
    # Make sure that there are idents available to use for labeling cells / clusters
    if (label && all(is.na(Idents(seurat)))) {
        
        label <- FALSE
        message("NOTE: identity of all cells is NA, setting 'label' to FALSE.")
        
    }
    
    # Deal with the palette
    if (is.null(color_by)) {
        
        if (is.null(colors)) {
            
            # Check if there's a cluster palette stored with the object
            if (!is.null(seurat@misc$colors)) colors <- seurat@misc$colors
            else if (!is.null(seurat@misc$colours)) colors <- seurat@misc$colours
            else {
                
                # Assuming that the order of the levels is correct in the seurat object,
                # find the default colors for the clusters
                colors <- get_gg_colors(length(levels(Idents(seurat))))
                names(colors) <- levels(Idents(seurat))
                
            }
            
        }
        
        if (rasterize) {
            
            gg <- gg +
                ggrastr::rasterize(geom_point(aes(color = Cluster), size = point_size, alpha = alpha), dpi = 500) +
                scale_color_manual(values = colors, na.value = na_color)
            
        } else {
            
            gg <- gg +
                geom_point(aes(color = Cluster), size = point_size, alpha = alpha) +
                scale_color_manual(values = colors, na.value = na_color)
            
        }
        
        
    } else {
        
        # If no limits are specified for the z-axis (variable used to color cells),
        # use the ggplot2 default, which is to fit range of the data
        if (is.null(limits)) lims <- c(NA, NA)
        else lims <- limits
        
        # Add the points to the plot object
        if (rasterize) {
            
            gg <- gg +
                ggrastr::rasterize(geom_point(aes_string(color = color_by), size = point_size, alpha = alpha), dpi = 500)
            
        } else {
            
            gg <- gg +
                geom_point(aes_string(color = color_by), size = point_size, alpha = alpha)
            
        }
        
        
        # If a palette is provided, choose an appropriate scale based on whether
        # the coloring variable is discrete or continuous
        if (!is.null(colors)) {
            
            if (color_by_type == "discrete") gg <- gg + scale_color_manual(values = colors, na.value = na_color)
            else if (color_by_type == "continuous") {
                
                gg <- gg + scale_color_gradientn(colors  = colors,
                                                 na.value = na_color,
                                                 limits   = lims)
            }
            
        } else {
            
            if (color_by_type == "continuous") { # Otherwise for discrete, default ggplot2 colors are used
                
                gg <- gg + scale_color_gradientn(colors  = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "OrRd"))(n = 100),
                                                 na.value = na_color,
                                                 limits   = lims)
                
            }
        }
    }
    
    # Label clusters on the plot
    if (label) {
        
        centers <- compute_cluster_centers(seurat,
                                           reduction = reduction,
                                           dim1 = dim1, dim2 = dim2)
        
        gg <- gg + add_labels(centers           = centers,
                              label_repel       = label_repel,
                              label_size        = label_size,
                              label_prefix_only = label_prefix_only,
                              label_suffix_only = label_suffix_only,
                              sep               = sep,
                              clusters          = clusters_to_label)
    }
    
    # Add the default theme to the plot object
    gg <- gg + theme_min2(border_size = border_size)
    
    # Set the right axes titles
    if (hide_axes) gg <- gg + xlab(NULL) + ylab(NULL)
    else if (reduction == "tsne") gg <- gg + xlab(glue("tSNE {dim1}")) + ylab(glue("tSNE {dim2}"))
    else if (reduction == "umap") gg <- gg + xlab(glue("UMAP {dim1}")) + ylab(glue("UMAP {dim2}"))
    else if (reduction == "phate") gg <- gg + xlab(glue("PHATE {dim1}")) + ylab(glue("PHATE {dim2}"))
    else if (reduction == "pca") {
        
        # If the reduction is PCA, automatically compute the variance explained
        var_exp <- getVarianceExplained(seurat)
        gg <- gg +
            xlab(glue("PC{dim1} ({round(var_exp$percent.var.explained[dim1], 1)}%)")) +
            ylab(glue("PC{dim2} ({round(var_exp$percent.var.explained[dim2], 1)}%)"))
        
    } else {
        gg <- gg + xlab(glue("{reduction} {dim1}")) + ylab(glue("{reduction} {dim2}"))
    }
    
    # Remove legend
    if (!legend) gg <- gg + hide_legend()
    else if (!is.null(color_by)) {
        
        # This column typically contains the sample name in the seurat metadata;
        # use this more informative name for the legend title instead of "orig.ident"
        if (color_by == "orig.ident") gg <- gg + labs(color = "Sample")
        
    }
    
    # Plot aesthetics
    if (!is.null(title)) gg <- gg + ggtitle(title)
    
    if (hide_ticks) gg <- gg + hide_ticks()
    
    if (constrain_scale) gg <- gg + constrain_scale(seurat,
                                                    reduction = reduction,
                                                    dim1      = dim1,
                                                    dim2      = dim2)
    
    return(gg)
    
}


#' @describeIn plot_dr Plot a tSNE embedding
#' @export
tsne <- function(seurat, ...) {
    
    plot_dr(seurat, reduction = "tsne", ...)
    
}


#' @describeIn plot_dr Plot a PCA embedding
#' @export
pca <- function(seurat, ...) {
    
    plot_dr(seurat, reduction = "pca", ...)
    
}


#' @describeIn plot_dr Plot a UMAP embedding
#' @export
umap <- function(seurat, ...) {
    
    plot_dr(seurat, reduction = "umap", ...)
    
}




violin_grid <- function(seurat, genes,
                        group_col,
                        group_order = NULL,
                        order = "clusters",
                        colours = NULL,
                        scale = "width",
                        title = NULL,
                        scales = "free_x",
                        points = FALSE,
                        gene_angle = 30) {
    
    expr <- FetchData(seurat, vars = c(group_col, genes))
    
    expr <- expr %>% tidyr::gather(Marker, Expression, 2:ncol(.))
    expr$Cell <- colnames(seurat)
    expr$Cluster <- expr[[group_col]]
    expr[[group_col]] <- NULL
    
    # if cluster order is provided, subset to cells in the provided clusters
    if (!is.null(group_order)) expr <- expr %>% filter(Cluster %in% group_order)
    
    # deal with all-0 expression
    maxes <- expr %>% group_by(Marker) %>% summarise(max = max(Expression)) %>% tibble::deframe()
    
    # convert to NA b/c of the jitter
    no_expr <- names(maxes)[maxes == 0]
    if (length(no_expr) > 0) expr[expr$Marker %in% no_expr, ]$Expression <- NA
    
    # complete for non-existent genes
    expr$Marker  <- factor(expr$Marker, levels = genes)
    expr$Cluster <- factor(expr$Cluster, levels = unique(expr$Cluster))
    expr <- expr %>% 
        complete(Cluster, nesting(Marker), fill = list(Expression = NA)) %>% 
        filter(!is.na(Cluster))
    
    genes_keep <- genes[genes %in% expr$Marker]
    
    # decide on ordering of clusters
    if (!is.null(group_order) & order == "group") {
        
        # 1. use provided ordering
        
        expr$Cluster <- factor(expr$Cluster, levels = rev(group_order))
        
    } else if (length(order) == 1 & order == "sort") {
        
        # 2. if sorting by expression (to create a diagonal), sort the clusters 
        # based on expression of all the detected genes
        
        # create sorting criteria in the form of:
        # desc(gene_a), desc(gene_b), etc
        sort_criteria <- glue::glue("desc({genes_keep})")
        
        cluster_order <- expr %>%
            filter(Expression != 0) %>%
            tidyr::spread(Marker, Expression) %>%
            group_by(Cluster) %>%
            dplyr::select(-Cell) %>% 
            summarise_all(funs(median), na.rm = TRUE) %>%
            arrange_(sort_criteria) %>%
            `[[`("Cluster")
        
        expr$Cluster <- factor(expr$Cluster, levels = rev(cluster_order))
        
    } else if (is.null(group_order) & order == "group") {
        
        # 3. use factor levels of the group column
        
        seurat@meta.data[[group_col]] <- factor(seurat@meta.data[[group_col]])
        expr$Cluster <- factor(expr$Cluster, levels = rev(levels(seurat@meta.data[[group_col]])))
        
    }
    
    expr$Marker <- factor(expr$Marker, levels = genes_keep)
    
    # Set colours to the levels of the clusters, so that they are preserved
    # even if an order was specified
    
    if (!is.null(colours) && length(colours) < length(unique(expr$Cluster))) {
        
        stop("Not enough colours! Please provide as many colours ",
             "as clusters in the dataset, or one per cluster specified in the ",
             "'subset_clusters' argument.")
        
        
    } else if (is.null(colours)) {
        
        colours <- get_gg_colors(length(unique(expr$Cluster)))
        names(colours) <- unique(expr$Cluster)
        
    }
    
    # plot violin plots
    gg <- expr %>%
        ggplot(aes(x = Cluster, y = Expression)) +
        geom_violin(aes(fill = Cluster), scale = scale, size = 0.5)
    
    if (points) gg <- gg + geom_jitter(size = 0.1, alpha = 0.3, width = 0.3)
    
    gg <- gg +
        scale_fill_manual(values = colours) +
        facet_wrap(~ Marker, ncol = length(unique(expr$Marker)),
                   scales = scales) +
        theme_min2() +
        coord_flip() +
        ggplot2::theme(panel.border = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.text.x = element_blank(),
                       axis.line.x = element_blank(),
                       legend.position = "none",
                       strip.text.x = element_text(angle = gene_angle, size = 8))
    
    if (!is.null(title)) gg <- gg + ggtitle(title)
    
    return(gg)
    
}


extract_meanexp_pct <- function(seurat, genes, cluster_col, clusters = NULL, scale = TRUE) {
    
    genes <- genes[genes %in% rownames(seurat)]
    
    Idents(seurat) <- cluster_col
    
    # subset to required clusters
    if (!is.null(clusters)) seurat <- subset(seurat, idents = clusters)
    
    # get average expression
    data_meanexp <- Seurat::AverageExpression(seurat, features = genes, assays = "RNA") %>% 
        .$RNA %>% 
        t()
    
    # re-scale to [0, 1]
    # if (scale == "minmax") data_meanexp <- minmax_norm(data_meanexp)
    if (scale) data_meanexp <- apply(data_meanexp, 2, scales::rescale)
    
    data_meanexp <- data_meanexp %>%
        as.data.frame() %>% 
        tibble::rownames_to_column(var = "Cluster") %>% 
        gather(Gene, Expression, 2:ncol(.))
    
    # get detection rate
    data_pct1 <- calc_pct1(seurat, genes) %>%
        data.frame() %>%
        gather(Gene, Pct1, 2:ncol(.))
    
    # the calc_pct1 function converts "-" in gene names to ".", so here we undo that
    data_pct1$Gene <- gsub("\\.", "-", data_pct1$Gene)
    
    data_all <- left_join(data_meanexp, data_pct1, by = c("Cluster", "Gene")) %>% 
        mutate(Gene = factor(Gene, levels = rev(genes))) %>% 
        replace_na(list(Expression = 0, Pct1 = 0))
    
    return(data_all)
    
    
}

#' @param flip logical, flip x/y axes and put clusters in rows and genes in columns
bubbleplot <- function(seurat, genes, cluster_order, cluster_col,
                       clusters = NULL, scale = TRUE, max_radius = 5,
                       flip = FALSE) {
    
    data_meanexp_pct <- extract_meanexp_pct(seurat, genes,
                                            cluster_col = cluster_col,
                                            clusters = clusters,
                                            scale = scale)
    
    data_plot <- data_meanexp_pct %>%
        filter(Pct1 > 0) %>%
        mutate(Gene = factor(Gene, levels = rev(genes))) %>%
        mutate(Cluster = factor(Cluster, levels = cluster_order)) %>% 
        filter(!is.na(Cluster))
    
    p1 <- data_plot %>%
        ggplot(aes(x = Cluster, y = Gene)) +
        geom_point(aes(size = Pct1, colour = Expression), alpha = 0.8) +
        scale_radius(range = c(0, max_radius), limits = c(0, 1)) +
        scale_color_gradientn(colours = tail(rdbu, 60)) +
        theme_min2() +
        rotate_x() +
        theme(panel.grid.major.x = element_line(colour = "grey90"),
              panel.border = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_text(size = 13))
    
    if (flip) p1 + coord_flip() + theme(panel.grid.major.y = element_line(colour = "grey90"),
                                        panel.grid.major.x = element_blank())
    else p1
    
}


# Gene signatures --------------------------------------------------------------

filter_markers <- function(markers, sp = "hg", n_top = 100) {
    
    markers %>% 
        # filter out mitochondrial and ribosomal genes
        dplyr::filter(., !grepl("RPS|RPL|MRPS|MRPL|^MT-", gene)) %>% 
        group_by(Cluster) %>% 
        # get positive markers
        dplyr::filter(avg_logFC > 0) %>% 
        arrange(desc(avg_logFC)) %>% 
        # make sure they're unique
        distinct(gene, .keep_all = TRUE) %>% 
        dplyr::slice(1:n_top) %>% 
        ungroup()
    
}

# Modified by Bhavyaa 2024-03-11
signatures_to_list <- function(signatures, 
                               type = c("ens", "sym"), 
                               species = c("hg", "mm")) {
    
    sig <- signatures %>% 
        select(Cluster, gene) %>%
        split(.$Cluster) %>%
        map(~ pull(.x, gene))
    
    if (type == "ens") {
        
        if(species == "hg"){
            
            sig <- sig %>% 
                symbols2ensembl()
            
        }
        
        else if (species == "mm") {
            
            sig <- sig %>% 
                symbols2ensembl_mm()
            
        }
        
    }
        
    return(sig)
}


# Summary stats ----------------------------------------------------------------

#' mean_gene_expression
#'
#' Calculate, for each cell, the mean expression of the given set of marker genes using the
#' normalized data.
#'
#' @param seurat Seurat object
#' @param genes String or character vector specifying gene(s) to use
#'
#' @return Data frame with two columns: "Cell" which specifies the cell ID
#' and "Mean_marker_expression" which is the expression value for the marker gene, 
#' or mean if multiple genes were provided.
#'
#' @export
#' @author Selin Jessa
#' @examples
#' mean_gene_expression(pbmc, genes = c("IL32", "MS4A1"))
mean_gene_expression <- function(seurat, genes, assay = "RNA") {
    
    # Get expression data from the Seurat object
    data.frame(Cell = colnames(GetAssayData(seurat, assay = assay)),
               Mean_marker_expression = rowMeans(fetch_data(seurat, genes, assay)),
               stringsAsFactors = FALSE)
    
}



# Utility functions ------------------------------------------------------------


prop <- function(x) {
    
    sum(x)/length(x)
    
}

#' Calculate proportion of cells in each cluster where a gene is detected,
#' for each level in \code{Idents(seurat)}
#'
#' @param seurat Seurat object
#' @param genes Character, gene or genes of interest
#'
#' @return
#' @export
calc_pct1 <- function(seurat, genes) {
    
    require(data.table)
    
    # message("## Processing ", seurat@project.name)
    seurat_subset <- seurat[genes, ]
    x <- as.matrix(seurat_subset@assays$RNA@data)[genes, ]
    # Binarize
    x[x > 0] <- 1
    x <- as.data.table(t(x))
    # Add cluster info for cells
    x[, Cluster := as.character(seurat@active.ident)]
    # Get prop of cells expressing a gene, within each cluster
    x[, lapply(.SD, prop), by = Cluster]
    
}


#' Compute a min-max normalization, where for each gene, expression
#' is subtracted by the mean and divided by the range across clusters
#' 
#' Written with input from Microsoft Bing AI
#' 
#' @param mat_exp Expression matrix, cluster x gene, where the values are 
#' mean expression of each gene (or genes of interest) in each cluster
minmax_norm <- function(mat_exp) {
    
    # Compute the minimum and maximum expression values for each gene (column)
    min_values <- apply(mat_exp, 2, min)
    max_values <- apply(mat_exp, 2, max)
    
    # Compute the range of expression values for each gene
    ranges <- max_values - min_values
    
    # Normalize the expression values for each gene
    mat_norm <- sweep(mat_exp, 2, min_values, "-")
    mat_norm <- sweep(mat_norm, 2, ranges, "/")
    
    # Handle genes with zero range (same expression across all groups)
    zero_range_genes <- which(ranges == 0)
    mat_norm[, zero_range_genes] <- 0.5
    
    return(mat_norm)
}


#' find_genes
#'
#' Given a set of genes, find the ones which are detected in the sample,
#' and which are not.
#'
#' @param seurat Seurat object
#' @param genes Character vector of genes
#'
#' @return A named list with two elements: "detected" and "undetected"
#' each storing character vectors with the genes in each category
#' @export
#'
#' @author Selin Jessa
#' @examples
#' find_out <- find_genes(pbmc, c("IL32", "CD79B", "foo"))
#' find_out$detected
#' find_out$undetected
find_genes <- function(seurat, genes, assay = "RNA") {
    
    genes_detected <- genes[genes %in% rownames(GetAssayData(object = seurat, assay = assay))]
    genes_undetected <- setdiff(genes, genes_detected)
    
    list(detected = genes_detected, undetected = genes_undetected)
    
}



#' fetchData
#'
#' Subset the seurat@@data matrix by gene and cluster.
#' Similar to SeuratObject::FetchData except it doesn't thrown an error if a gene
#' is not found in the data, and accepts a specific assay argument.
#'
#' @param seurat Seurat object
#' @param genes Genes to filter
#' @param assay Assay to search for genes
#' @param clusters (Optional) Vector, include only cells with these identities
#' (e.g. cluster assignments). Searches in seurat@@active.ident.
#' @param return_cell Logical, whether or not to include a column with the cell ID.
#' Default: FALSE
#' @param return_cluster Logical, whether or not to include a column with the cluster.
#' Default: FALSE
#' @param scaled Logical, whether to fetch scaled data the assay.
#' Default: FALSE.
#'
#' @return Expression matrix for genes specified
#' @export
#'
#' @author Selin Jessa
#' @examples
#' fetchData(pbmc, c("IL32", "MS4A1"))
#' fetchData(pbmc, c("IL32"), c(1, 2))
#' fetchData(pbmc, c("IL32", "MS4A1"), c(1, 2), return_cluster = TRUE, return_cell = TRUE)
fetch_data <- function(seurat, features, assay = "RNA", clusters = NULL,
                       return_cell = FALSE, return_cluster = FALSE, scaled = FALSE) {
    
    message("@ Searching assay: ", assay)
    
    features_out <- find_genes(seurat, genes = features, assay)
    
    n_undetected <- length(features_out$undetected)
    
    if (n_undetected > 0) {
        
        if (n_undetected > 10) {
            
            message(paste0("@ NOTE: [",
                           paste0(head(features_out$undetected, 10), collapse = ", "),
                           "] and ", n_undetected - 10, " other features are undetected in ", seurat@project.name))
            
        } else {
            
            message(paste0("@ NOTE: [",
                           paste0(features_out$undetected, collapse = ", "),
                           "] undetected in ", seurat@project.name))
            
        }
    }
    
    if (length(features_out$detected) == 0) stop("No features specified were ",
                                                 "found in the data.")
    
    if (scaled) exp <- as.matrix(seurat@assays[[assay]]@scale.data)
    else exp <- as.matrix(seurat@assays[[assay]]@data)
    
    exp_filt <- as.data.frame(t(exp[which(rownames(exp) %in% features_out$detected),]))
    
    # Keep all
    if(is.null(clusters)) clusters <- levels(seurat@active.ident)
    ident_idx <- which(seurat@active.ident %in% clusters)
    
    # Handle only one gene case, and properly return a data frame
    if(nrow(exp_filt) == 1) {
        
        if(!is.null(ident)) exp_filt <- exp_filt[ident_idx]
        
        exp_filt <- as.data.frame(t(exp_filt))
        names(exp_filt) <- rownames(exp)[rownames(exp) %in% features]
        
    } else {
        if(!is.null(ident)) exp_filt <- exp_filt[ident_idx,]
    }
    
    rownames(exp_filt) <- c() # Get rid of rownames
    
    if(return_cluster) exp_filt <- tibble::add_column(exp_filt, Cluster = seurat@active.ident[ident_idx], .before = 1)
    if(return_cell) exp_filt <- tibble::add_column(exp_filt, Cell = names(seurat@active.ident)[ident_idx], .before = 1)
    
    return(exp_filt)
    
}



#' get_embedding
#' 
#' Given a Seurat object and a data frame where the rows correspond to cells and a Cell column
#' is present, the dataframe is altered to add two columns giving coordinates in a dimensionality
#' reduced space.
#' 
#' @param seurat Seurat object, where dimensionality reduction has been applied,
#' i.e. (after applying Seurat::RunPCA(), Seurat::RunTSNE(), Seurat::RunUMAP() to the object).
#' @param df Data frame with at least one column, "Cell", giving the cell ID
#' @param reduction String specifying the dimensionality reduction to use,
#' retrieves UMAP by default.
#' Default: "umap"
#' @param dim1 Numeric, dimension of embedding to use for x-axis. Default = 1.
#' @param dim2 Numeric, dimension of embedding to use for y-axis. Default = 2.
#' 
#' @export
#' @author Selin Jessa, Samantha Worme
get_embedding <- function(seurat, assay = "RNA", df, reduction = "umap", dim1 = 1, dim2 = 2) {
    
    # Get the axes for the reduced space
    # See here: http://dplyr.tidyverse.org/articles/programming.html#setting-variable-names
    vars <- colnames(seurat[[reduction]]@cell.embeddings)[c(dim1, dim2)]
    
    df$Cell <- as.character(df$Cell)
    
    embedding <- data.frame(Cell = colnames(GetAssayData(seurat, assay = assay)), stringsAsFactors = FALSE) %>%
        mutate(!!vars[1] := seurat[[reduction]]@cell.embeddings[, dim1], 
               !!vars[2] := seurat[[reduction]]@cell.embeddings[, dim2])
    
    df <- dplyr::inner_join(df, embedding, by = "Cell")
    return(df)
    
}


#' Rotate the x axis labels in a ggplot
#'
#' @param angle Integer, value in degrees to rotate labels. Default: 90.
#'
#' @return A theme element to rotate labels
#' @export
#'
#' @author Selin Jessa
#'
#' @examples
#' # gg <- mpg %>%
#' # filter(class %in% c("compact", "suv")) %>%
#' # ggplot(aes(x = displ, y = hwy)) +
#' # geom_point(aes(color = factor(year))) +
#' # facet_wrap(~ class, ncol = 2)
#'
#' # gg
#' # gg + rotateX()
rotate_x <- function(angle = 90) {
    
    theme(axis.text.x = element_text(angle = angle, hjust = 1))
    
}



#' Remove the legend in a ggplot
#'
#' @return A theme element to hide legend
#' @export
#'
#' @author Selin Jessa
#'
#' @examples
#' # gg <- mpg %>%
#' # filter(class %in% c("compact", "suv")) %>%
#' # ggplot(aes(x = displ, y = hwy)) +
#' # geom_point(aes(color = factor(year))) +
#' # facet_wrap(~ class, ncol = 2)
#'
#' # gg
#' # gg + noLegend()
hide_legend <- function() {
    
    theme(legend.position = "none")
    
}



#' Remove axis ticks and tick labels from a ggplot
#'
#' @return A theme element to remove ticks
#' @export
#'
#' @author Selin Jessa
hide_ticks <- function() {
    
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank())
    
}


#' compute_cluster_centers
#'
#' Get centers of clusters given a Seurat object, to use for labelling
#' in tSNE space. The cluster center is defined as the median X and Y coordinate
#' across cells in each cluster.
#'
#' @param seurat Seurat object, where dimensionality reduction has been applied,
#' i.e. (after applying Seurat::RunTSNE() to the object).
#'
#' @return Data frame with three columns: Cluster, mean_tSNE_1, and mean_tSNE_2
#' @export
#'
#' @author Selin Jessa
#' @examples
#'
#' compute_cluster_centers(pbmc, reduction = "pca", dim1 = 1, dim2 = 3)
compute_cluster_centers <- function(seurat, reduction = "umap", dim1 = 1, dim2 = 2) {
    
    n_clusters <- length(unique(Idents(seurat)))
    
    # Attempts at tidyeval...
    # vars <- colnames(seurat[[reduction]]@cell.embeddings)[c(1, 2)]
    # col_names <- paste0("mean_", vars)
    
    # Get the embedding
    df <- as.data.frame(seurat[[reduction]]@cell.embeddings[, c(dim1, dim2)]) %>%
        mutate(Cell = colnames(seurat),
               Cluster = Idents(seurat))
    
    # Generalize these
    colnames(df)[c(1, 2)] <- c("Dim_1", "Dim_2")
    
    # Compute cluster centers
    centers <- df %>%
        group_by(Cluster) %>%
        dplyr::summarise(mean_x = median(Dim_1),
                         mean_y = median(Dim_2))
    
    return(centers)
    
}


#' Add cluster labels to a deminsionality reduction plot
#'
#' @param centers Data frame with at least three columns: "mean_x", "mean_y",
#' and "Cluster", as returned by \code{\link{clusterCenters}}
#' @param label_repel Logical, whether to plot cluster
#' labels repelled from the center, on a slightly transparent white background and
#' with an arrow pointing to the cluster center. If FALSE, simply plot the
#' cluster label at the cluster center. Default: TRUE.
#' @param label_size Numeric, controls the size of text labels. Default: 4.
#' @param label_short (Optional/Experimental!!) Logical, if TRUE, assumes clusters
#' (at \code{seurat@@ident}) consist of a prefix and a suffix separated by a non-alpha
#' numeric character (\code{"[^[:alnum:]]+"}), and tries to separate these names
#' and only plot the prefix, for shorter labels and a cleaner plot. Default: FALSE.
#' @param clusters (Optional) Clusters for which labels should be plot (if only
#' a subset of clusters should be labelled). Default: NULL (Label all clusters).
#'
#'
#' @author Selin Jessa and Nisha Kabir
#' @export
add_labels <- function(centers, label_repel = FALSE, label_size = 4, label_prefix_only = FALSE, label_suffix_only = FALSE, sep = NULL, clusters = NULL) {
    
    if (!is.null(clusters)) centers <- filter(centers, Cluster %in% clusters)
    
    if (label_prefix_only) centers <- suppressWarnings(tidyr::separate(centers, Cluster, into = c("Cluster", "Cluster_long"), extra = "drop", sep = sep))
    else if (label_suffix_only) centers <- suppressWarnings(tidyr::separate(centers, Cluster, into = c("Cluster_long", "Cluster"), extra = "drop", sep = sep))
    
    if (label_repel) {
        
        ggrepel::geom_label_repel(data = centers,
                                  aes(x = mean_x, y = mean_y),
                                  label = centers$Cluster,
                                  size = label_size,
                                  segment.color = 'grey50',
                                  fontface = 'bold',
                                  alpha = 0.8,
                                  segment.alpha = 0.8,
                                  label.size = NA,
                                  force = 1,
                                  # Leaving these unspecified for now, since it really depends on
                                  # the dimensionality reduction
                                  # nudge_x = 5, nudge_y = 5,
                                  segment.size = 0.5,
                                  arrow = arrow(length = unit(0.01, 'npc')))
        
    } else {
        
        geom_text(data = centers,
                  aes(x = mean_x, y = mean_y, label = Cluster),
                  size = label_size)
        
    }
    
}



#' Get the limits of a the first two dimensions in a dimensionality reduction
#'
#' When plotting an embedding, we may want to plot specific cells, but
#' constrain the scale to match plots of the whole dataset. Given a dim.
#' reduction, this function extracts the x and y limits to use for plotting.
#'
#' @param seurat Seurat object for which a dimensionality reduction has been
#' computed (e.g. PCA or tSNE)
#' @param reduction String, corresponding to the dimensionality reduction to use.
#' Default: "tsne".
#' @param dim1 Numeric, dimension of embedding to use for x-axis. Default = 1.
#' @param dim2 Numeric, dimension of embedding to use for y-axis. Default = 2.
#'
#' @return A list with two elements: "xlim", which is a character vector of
#' the limits for the x-axis, and "ylim", correspondingly for the y-axis
#' @export
#' @author Selin Jessa
#'
#' @examples
#' get_dr_lims(pbmc, reduction = "tsne")
get_dr_lims <- function(seurat, reduction, dim1 = 1, dim2 = 2) {
    
    dim1 <- seurat[[reduction]]@cell.embeddings[, dim1]
    dim2 <- seurat[[reduction]]@cell.embeddings[, dim2]
    
    return(list(xlim = c(min(dim1), max(dim1)),
                ylim = c(min(dim2), max(dim2))))
    
}


#' Constrain the scale of the plot to the dimensionality reduction limits
#'
#' @inheritParams get_dr_lims
#'
#' @export
#' @author Selin Jessa
constrain_scale <- function(seurat, reduction, dim1 = 1, dim2 = 2)  {
    
    limits <- get_dr_lims(seurat = seurat, reduction = reduction, dim1 = dim1, dim2 = dim2)
    lims(x = limits$xlim, y = limits$ylim)
    
}


#' Apply a clean theme to a ggplot2 object
#'
#' @references https://github.com/sjessa/ggmin
#'
#' @importFrom ggplot2 theme_light theme
#' @author Selin Jessa
#' @export
theme_min2 <- function(base_size = 11, base_family = "",
                       border_color = "grey90",
                       border_size = 1) {
    
    theme_light(base_size = 11, base_family = "") +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, color = border_color, size = border_size),
            axis.ticks = element_line(color = border_color),
            strip.background = element_rect(fill = NA, color = NA),
            strip.text.x = element_text(color = "black", size = rel(1.2)),
            strip.text.y = element_text(color = "black", size = rel(1.2)),
            title = element_text(size = rel(0.9)),
            axis.text = element_text(color = "black", size = rel(0.8)),
            axis.title = element_text(color = "black", size = rel(0.9)),
            legend.title = element_text(color = "black", size = rel(0.9)),
            legend.key.size = unit(0.9, "lines"),
            legend.text = element_text(size = rel(0.7), color = "black"),
            legend.key = element_rect(color = NA, fill = NA),
            legend.background = element_rect(color = NA, fill = NA)
        )
}



# get top and bottom values
get_top_and_bottom <- function(df, wt, n = 10) {
    
    wt_var <- rlang::enquo(wt)
    
    df %>%
        {bind_rows(top_n(., n, !!wt_var),
                   top_n(., -n, !!wt_var))}
    
}


# a distance function for Spearman rank correlation
# from https://davetang.org/muse/2010/11/26/hierarchical-clustering-with-p-values-using-spearman/
spearman <- function(x, ...) {
    x <- as.matrix(x)
    res <- as.dist(1 - cor(x, method = "spearman", use = "everything"))
    res <- as.dist(res)
    attr(res, "method") <- "spearman"
    return(res)
}



#' Get default ggplot2/Seurat colors
#'
#' Get evenly spaced colors from around the color wheel, which are the default
#' colors assigned to clusters by Seurat. The output of this function can be
#' passed to the \code{scale_color_manual()} and \code{scale_fill_manual()} functions
#' from ggplot2, as the \code{values} argument. (\code{\link{ggColors}} points
#' to this function.)
#'
#' @param n Number of colors to return
#'
#' @return Named character vector, where names are the names of clusters, from
#' 0 to n-1, and values are the hex codes for the colors.
#' @export
#'
#' @examples
#'
#' n_clust <- 5
#' get_gg_colors(n_clust)
#'
#' @references https://stackoverflow.com/a/8197703
#' @aliases get_gg_colours
#' @importFrom grDevices hcl
get_gg_colors <- function(n) {
    
    hues <- seq(15, 375, length = n + 1)
    colors <- hcl(h = hues, l = 65, c = 100)[1:n]
    names(colors) <- seq(0, n - 1) # Since the first cluster in Seurat is 0
    
    return(colors)
    
}



symbols2ensembl <- function(genes) {
    
    load(here("data/singlecell/references_genome/biomaRt_hg_symbol_to_ens_lds.Rda"))
    genes_ens_lds %>% filter(HGNC.symbol %in% genes) %>% pull(Gene.stable.ID)
    
}

symbols2ensembl_mm <- function(genes) {
    
    load(here("data/singlecell/references_genome/biomaRt_mm_symbol_to_ens_lds.Rda"))
    genes_mm_ens_lds %>% filter(MGI.symbol %in% genes) %>% pull(Gene.stable.ID)
    
}


convert_mat_mm2hg <- function(mat, use_cache = TRUE) {
    
    require(biomaRt)
    require(data.table)
    
    if (!use_cache) {
        
        human  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
        mouse  <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
        genes2 <- getLDS(attributes = c("mgi_symbol"),
                         filters = "mgi_symbol",
                         values = rownames(mat),
                         mart = mouse,
                         attributesL = c("hgnc_symbol"),
                         martL = human,
                         uniqueRows = TRUE)
        
    } else {
        
        load(here("data/singlecell/references_genome/biomaRt_mm_to_hg_lds.Rda"))
        genes2 <- genes_lds %>% filter(MGI.symbol %in% rownames(mat))
        
    }
    
    # convert counts matrix to data.table, and put rownames in a new column
    dt <- as.data.table(mat, keep.rownames = "gene")
    
    # restrict to rows with human gene symbols
    dt <- dt[gene %in% genes2[, "MGI.symbol"], ]
    
    # get the equivalent human symbols in the same order as in the counts mat
    genes_conv_ordered <- match(dt$gene, genes2[,1])
    
    # add to the dt    
    dt$gene_hg <- genes2$HGNC.symbol[genes_conv_ordered]
    dt <- dt[,gene:=NULL]
    
    # remove duplicates by getting the index of the first instance of each gene
    idx_first_instance <- dt[, .(first_inst = min(.I)), by = gene_hg][, first_inst]
    dt <- dt[idx_first_instance]
    
    # convert back to matrix with human genes in the rownames
    mat_conv <- as.matrix(dt, rownames = "gene_hg")
    
    return(mat_conv)
    
}




# scCoAnnotate helpers ---------------------------------------------------------


#' get the mode of a vector v, from Samantha Worme
get_mode <- function(v) {
    uniqv     <- unique(v)
    matches   <- tabulate(match(v, uniqv))
    max_match <- max(matches)
    ties      <- ifelse(length(which(matches == max_match)) > 1, T, F)
    
    if (max_match == 1) {
        return("No consensus")
    } else if (ties) {
        return("No consensus")
    } else {
        return(uniqv[which.max(matches)])
    }
}


#' get consensus of labels from several automated methods, adapted from Samantha Worme
#' 
#' @param auto_labels Data frame containing predictions from Prediction_Summary.tsv file as
#' output by scCoAnnotate, or path to tsv file
#' @param tools_to_include Vector of tools to include in consensus
#' @param high_level Logical, whether to use high-level (broader labels) when calculating the consensus.
#' If TRUE, output will contain these values in the column "Consensus_broad"
#' @param ontology Data frame with the higher-level annotation, should contain two columns:
#' "Label", matching the cell types in the reference, and "Label_broad", with the high-level annotation
get_consensus_labels <- function(auto_labels,
                                 tools_to_include,
                                 high_level = FALSE,
                                 ontology = NULL,
                                 remove_timepoint = TRUE,
                                 suffix = NULL,
                                 keep_auto_labels = FALSE) {
    
    if (is.character(auto_labels)) auto_labels <- suppressMessages(read_tsv(auto_labels))
    
    predicted_labels <- auto_labels %>%
        tibble::column_to_rownames(var = "cellname") %>%
        as.matrix()
    predicted_labels <- predicted_labels[, colnames(predicted_labels) %in% tools_to_include]
    print(colnames(predicted_labels))
    
    if (remove_timepoint) {
        
        predicted_labels_no_time <- as.matrix(apply(predicted_labels, 2, FUN = function(x) gsub("^[A-Z]-[a-z][0-9]*_", "", x)))
        predicted_labels_no_time <- as.matrix(apply(predicted_labels_no_time, 2, FUN = function(x) gsub("ASTR2-2", "ASTR", x)))
        predicted_labels_no_time <- as.matrix(apply(predicted_labels_no_time, 2, FUN = function(x) gsub("[0-9]", "", x)))
        predicted_labels_no_time <- as.matrix(apply(predicted_labels_no_time, 2, FUN = function(x) gsub("^-", "", x)))
        predicted_labels <- predicted_labels_no_time
        
    }
    
    # calculate consensus
    consensus <- apply(predicted_labels, 1, get_mode) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Cell")
    colnames(consensus)[2] <- "Consensus"
    
    if (!is.null(suffix)) {
        
        colnames(consensus)[2] <- paste0(colnames(consensus)[2], "_", suffix)
        
    }
    
    # Updated by Bhavyaa to reflect new CoRAL output value
    message("@ Proportion of unresolved cells: ", round(table(consensus$Consensus)["Unresolved"]/nrow(consensus), 2))
    
    # if this option is set, then use a higher-level cell annotation provided
    # in ontology to first convert cell labels to the higher-level annotations,
    # prior to taking the consensus
    if (high_level) {
        
        tmp <- predicted_labels %>%
            as.data.frame() %>% 
            tibble::rownames_to_column(var = "Cell") %>% 
            pivot_longer(-Cell, names_to = 'Tool', values_to = 'Prediction') %>%
            left_join(ontology %>% select(Label, Label_broad), by = c('Prediction' = 'Label')) %>%
            pivot_wider(-Prediction, names_from = 'Tool', values_from = 'Label_broad')
        
        # sanity check that order was preserved
        if (!all(consensus$Cell == tmp$Cell)) stop("Cell order is not preserved")
        
        consensus_broad <- apply(select(tmp, -Cell) %>% as.matrix(), 1, get_mode)
        consensus$Consensus_broad <- consensus_broad
        
        message("@ Proportion of uncertain cells (broad): ", round(table(consensus$Consensus_broad)["Unsure"]/nrow(consensus), 2))
        
        if (!is.null(suffix)) {
            
            colnames(consensus)[3] <- paste0(colnames(consensus)[3], "_", suffix)
            
        }
        
    }
    
    if (keep_auto_labels) return(bind_cols(auto_labels, consensus))
    else return(consensus)
    
}




#' Plot a heatmap of the correlation in predictions between methods
#' adapted from Alva Annett
#' 
#' @param auto_labels Data frame as returned by \code{get_consensus_labels}
#' @param cell_col Character, name of column in auto_labels which contains cell barcodes
plot_label_cor <- function(auto_labels, tools_to_include, cell_col) {
    
    auto_labels <- auto_labels[, colnames(auto_labels) %in% c(tools_to_include, "Consensus", cell_col)]
    
    # code adapted from Alva Annett
    pred <- auto_labels %>%
        column_to_rownames(cell_col) %>%
        select(correlation, SVM, SingleCellNet, SciBet, Consensus) %>% 
        mutate(across(where(is.character), ~ifelse(. == 'unassigned', 'Unsure', .))) %>%
        mutate(across(where(is.character), ~ifelse(. == 'Unassigned', 'Unsure', .))) %>%
        mutate(across(where(is.character), ~ifelse(. == 'Unknown', 'Unsure', .))) %>%
        mutate(across(where(is.character), ~ifelse(grepl('Node', .), 'Unsure', .))) %>%
        mutate(across(where(is.character), ~ifelse(. %in% c('rand', 'Rand'), 'Unsure', .))) %>%
        mutate(across(where(is.character), ~str_replace_all(., '[.]', ''))) 
    
    cor <- pred %>%
        select(is.character) %>%
        rownames_to_column('cell') %>%
        pivot_longer(!cell) %>%
        mutate(value = factor(value)) %>%
        mutate(value = as.numeric(value)) %>%
        pivot_wider(names_from = name, values_from = value) %>%
        column_to_rownames('cell') %>%
        cor()
    
    col_fun <- circlize::colorRamp2(seq(-1, 1, length.out = 100), rdbu)
    
    count_unsure <- pred %>%
        select(is.character) %>%
        rownames_to_column('cell') %>%
        pivot_longer(!cell) %>%
        filter(value == 'Unsure') %>%
        dplyr::count(name, .drop = F) %>%
        mutate(freq = round(n/nrow(pred)*100)) %>%
        select(!n) %>%
        column_to_rownames('name')
    
    count_unsure[setdiff(names(pred %>% select(is.character)), rownames(count_unsure)),] = 0
    count_unsure = count_unsure[order(match(rownames(count_unsure), colnames(cor))), , drop = FALSE]
    
    ha <- columnAnnotation('% Unsure' = anno_barplot(count_unsure, border = F, gp = gpar(fill = '#596475', col = '#596475')))
    
    ComplexHeatmap::Heatmap(cor,
                            name = 'Correlation',
                            col = col_fun,
                            width = ncol(cor)*unit(7, "mm"),
                            height = nrow(cor)*unit(7, "mm"),
                            rect_gp = gpar(col = "white", lwd = 2), 
                            top_annotation = ha, 
                            show_column_dend = F)
    
}


