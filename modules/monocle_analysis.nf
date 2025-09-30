/*
 * Monocle3 trajectory analysis and pseudotime analysis module
 */

process MONOCLE_ANALYSIS {
    publishDir "${params.output_dir}/6.Pseudotime", mode: 'copy'
    
    input:
    path seurat_rds
    
    output:
    path "monocle3_cds_object.rds", emit: monocle_cds
    path "*.png"
    path "*.pdf"
    path "*.csv"
    path "monocle3_summary.txt"
    
    script:
    """
    #!/usr/bin/Rscript
    
    # Load required packages
    suppressMessages({
        library(Seurat)
        library(monocle3)
        library(SeuratWrappers)
        library(ggplot2)
        library(dplyr)
        library(patchwork)
        library(viridis)
    })
    
    cat("=== Monocle3 trajectory analysis and pseudotime analysis ===\\n")
    
    # Read annotated Seurat object
    cat("Reading annotated data...\\n")
    st <- readRDS("${seurat_rds}")
    
    cat("Input cell count:", ncol(st), "\\n")
    cat("Input gene count:", nrow(st), "\\n")
    
    tryCatch({
        # Check if there are enough cell types for trajectory analysis
        if("cell_type_filtered" %in% colnames(st@meta.data)) {
            celltype_counts <- table(st@meta.data\$cell_type_filtered)
            valid_celltypes <- names(celltype_counts)[celltype_counts >= ${params.monocle.min_cells_per_type}]
            valid_celltypes <- valid_celltypes[valid_celltypes != "Unknown"]
            group_column <- "cell_type_filtered"
        } else {
            # Use clustering information
            cluster_counts <- table(st@meta.data\$seurat_clusters)
            valid_celltypes <- names(cluster_counts)[cluster_counts >= ${params.monocle.min_cells_per_type}]
            group_column <- "seurat_clusters"
        }
        
        if(length(valid_celltypes) < 2) {
            cat("Insufficient cell types/clusters (<2), cannot perform trajectory analysis\\n")
            
            # Create empty output files
            empty_cds <- list(data = matrix(0, nrow = 1, ncol = 1))
            saveRDS(empty_cds, "monocle3_cds_object.rds")
            
            write.csv(data.frame(), "trajectory_dependent_genes.csv", row.names = FALSE)
            
            png("monocle_insufficient_data.png", width = 8, height = 6, units = "in", res = 300)
            plot.new()
            text(0.5, 0.5, paste("Insufficient data for trajectory analysis\\nRequires at least 2 cell types/clusters\\nEach with >=", 
                                 ${params.monocle.min_cells_per_type}, "cells"), 
                 cex = 1.2, col = "orange")
            dev.off()
            
            writeLines(c(
                "=== Monocle3 Analysis Summary ===",
                "Insufficient data for trajectory analysis",
                paste("Requires at least 2 cell types, each with >=", ${params.monocle.min_cells_per_type}, "cells")
            ), "monocle3_summary.txt")
            
            quit("no", status = 0)
        }
        
        # Filter data to keep only valid cell types
        if(group_column == "cell_type_filtered") {
            valid_cells <- st@meta.data\$cell_type_filtered %in% valid_celltypes
        } else {
            valid_cells <- st@meta.data\$seurat_clusters %in% valid_celltypes
        }
        
        st_subset <- st[, valid_cells]
        cat("Using", group_column, "for trajectory analysis\\n")
        cat("Number of valid types:", length(valid_celltypes), "\\n") 
        cat("Number of cells for analysis:", ncol(st_subset), "\\n")
        
        # Convert Seurat object to cell_data_set
        cat("Converting Seurat object to Monocle3 format...\\n")
        
        # Use SeuratWrappers for conversion
        cds <- as.cell_data_set(st_subset)
        
        # Add gene information
        fData(cds)\$gene_short_name <- rownames(cds)
        
        # Preprocess data
        cat("Preprocessing Monocle3 data...\\n")
        cds <- preprocess_cds(cds, num_dim = ${params.monocle.num_dim}, method = "PCA")
        
        # Dimensionality reduction visualization
        cat("Performing ${params.monocle.reduction_method} dimensionality reduction...\\n")
        cds <- reduce_dimension(cds, reduction_method = "${params.monocle.reduction_method}")
        
        # Cluster cells
        cat("Clustering cells...\\n")
        cds <- cluster_cells(cds)
        
        # Learn trajectory graph
        cat("Learning cell trajectory...\\n")
        cds <- learn_graph(cds)
        
        # Basic visualization
        cat("Generating basic visualization plots...\\n")
        
        # 1. UMAP partition plot
        p1 <- plot_cells(cds,
                         color_cells_by = "partition",
                         label_cell_groups = FALSE,
                         label_leaves = FALSE,
                         label_branch_points = FALSE,
                         graph_label_size = 1.5) +
            theme_bw() +
            ggtitle("Monocle3: Cell Partitions")
        
        # 2. UMAP plot colored by group
        if(group_column %in% colnames(colData(cds))) {
            p2 <- plot_cells(cds,
                           color_cells_by = group_column,
                           label_cell_groups = FALSE,
                           label_leaves = FALSE,
                           label_branch_points = FALSE,
                           graph_label_size = 1.5) +
                theme_bw() +
                ggtitle(paste("Monocle3: Cells by", tools::toTitleCase(gsub("_", " ", group_column))))
        } else {
            p2 <- plot_cells(cds,
                           color_cells_by = "cluster",
                           label_cell_groups = FALSE,
                           label_leaves = FALSE,
                           label_branch_points = FALSE,
                           graph_label_size = 1.5) +
                theme_bw() +
                ggtitle("Monocle3: Cells by Cluster")
        }
        
        # 3. Trajectory plot - colored by Cell Type
        if(group_column %in% colnames(colData(cds))) {
            p3 <- plot_cells(cds,
                             color_cells_by = group_column,
                             label_groups_by_cluster = FALSE,
                             label_leaves = TRUE,
                             label_branch_points = TRUE,
                             graph_label_size = 1.5) +
                theme_bw() +
                ggtitle("Monocle3: Trajectory Graph")
        } else {
            p3 <- plot_cells(cds,
                             color_cells_by = "cluster",
                             label_groups_by_cluster = FALSE,
                             label_leaves = TRUE,
                             label_branch_points = TRUE,
                             graph_label_size = 1.5) +
                theme_bw() +
                ggtitle("Monocle3: Trajectory Graph")
        }
        
        # Save basic trajectory plots
        p_combined <- (p1 | p2) / p3
        ggsave(plot = p_combined, 
               filename = "01_monocle3_trajectory_overview.png",
               width = 16, height = 12, dpi = 300)
        ggsave(plot = p_combined, 
               filename = "01_monocle3_trajectory_overview.pdf",
               width = 16, height = 12)
        
        # Pseudotime analysis
        cat("Calculating pseudotime...\\n")
        
        # Attempt to automatically select root node for pseudotime analysis
        pseudotime_success <- FALSE
        
        tryCatch({
            # Get leaf nodes
            leaf_nodes <- names(which(igraph::degree(principal_graph(cds)[["UMAP"]]) == 1))
            
            if(length(leaf_nodes) > 0) {
                # Select leaf node with most cells as root node
                root_node <- leaf_nodes[1]  # Simply select first leaf node
                cds <- order_cells(cds, root_cells = colnames(cds)[1:min(50, ncol(cds))])
                
                # Pseudotime visualization
                p4 <- plot_cells(cds,
                                 color_cells_by = "pseudotime",
                                 label_cell_groups = FALSE,
                                 label_leaves = TRUE,
                                 label_branch_points = TRUE,
                                 graph_label_size = 1.5) +
                    theme_bw() +
                    ggtitle("Monocle3: Pseudotime") +
                    scale_color_viridis_c()
                
                # Display pseudotime by cell type
                if(group_column %in% colnames(colData(cds))) {
                    p5 <- plot_cells(cds,
                                     color_cells_by = "pseudotime",
                                     label_cell_groups = FALSE,
                                     label_leaves = FALSE,
                                     label_branch_points = FALSE,
                                     graph_label_size = 1.5) +
                        facet_wrap(vars(!!sym(group_column))) +
                        theme_bw() +
                        ggtitle("Monocle3: Pseudotime by Cell Type") +
                        scale_color_viridis_c()
                } else {
                    p5 <- plot_cells(cds,
                                     color_cells_by = "pseudotime",
                                     label_cell_groups = FALSE,
                                     label_leaves = FALSE,
                                     label_branch_points = FALSE,
                                     graph_label_size = 1.5) +
                        facet_wrap(vars(cluster)) +
                        theme_bw() +
                        ggtitle("Monocle3: Pseudotime by Cluster") +
                        scale_color_viridis_c()
                }
                
                # Save independent pseudotime plots
                # p4: Pure Pseudotime plot
                # p5: Pseudotime by Cell Type plot
                # Save Pseudotime plot
                ggsave(plot = p4, 
                       filename = "02_monocle3_pseudotime.png",
                       width = 8, height = 6, dpi = 300)
                ggsave(plot = p4, 
                       filename = "02_monocle3_pseudotime.pdf",
                       width = 8, height = 6)
                
                # Save Pseudotime by Cell Type plot
                ggsave(plot = p5, 
                       filename = "02_monocle3_pseudotime_by_celltype.png",
                       width = 16, height = 6, dpi = 300)
                ggsave(plot = p5, 
                       filename = "02_monocle3_pseudotime_by_celltype.pdf",
                       width = 16, height = 6)
                
                pseudotime_success <- TRUE
                cat("Pseudotime analysis completed\\n")
            } else {
                cat("Cannot find leaf nodes, skipping pseudotime analysis\\n")
            }
        }, error = function(e) {
            cat("Pseudotime analysis failed:", e\$message, "\\n")
        })
        
        # Find trajectory-dependent genes
        cat("Finding trajectory-related genes...\\n")
        
        trajectory_genes_success <- FALSE
        tryCatch({
            # Select genes with highest expression variation for analysis
            cds_subset <- cds[Matrix::rowSums(exprs(cds) > 0) >= 10,]
            
            if(nrow(cds_subset) > 100) {
                # Further filter genes with highest variation
                cds_subset <- cds_subset[1:min(1000, nrow(cds_subset)),]
            }
            
            # Find genes differentially expressed along pseudotime
            if(pseudotime_success) {
                pr_test_res <- graph_test(cds_subset, neighbor_graph="knn", cores=1)
                pr_deg_ids <- row.names(subset(pr_test_res, q_value < 0.05))
                
                if(length(pr_deg_ids) > 0) {
                    cat("Found", length(pr_deg_ids), "trajectory-related genes\\n")
                    
                    # Save trajectory-related gene results
                    pr_test_res\$gene_short_name <- rownames(pr_test_res)
                    pr_test_res <- pr_test_res[order(pr_test_res\$q_value),]
                    write.csv(pr_test_res, "trajectory_dependent_genes.csv", row.names = TRUE)
                    
                    # Select top genes for visualization
                    if(length(pr_deg_ids) >= 6) {
                        top_genes <- head(pr_deg_ids, 6)
                        
                        p_genes <- plot_genes_in_pseudotime(cds_subset[top_genes,], 
                                                           color_cells_by = "pseudotime",
                                                           min_expr = 0.1) +
                            theme_bw() +
                            ggtitle("Top Trajectory-Dependent Genes")
                        
                        ggsave(plot = p_genes, 
                               filename = "03_trajectory_genes.png",
                               width = 12, height = 9, dpi = 300)
                        ggsave(plot = p_genes, 
                               filename = "03_trajectory_genes.pdf",
                               width = 12, height = 9)
                        
                        trajectory_genes_success <- TRUE
                        cat("Trajectory gene visualization completed\\n")
                    }
                } else {
                    cat("No significant trajectory-related genes found\\n")
                    write.csv(data.frame(), "trajectory_dependent_genes.csv", row.names = FALSE)
                }
            } else {
                cat("Skipping trajectory gene analysis due to pseudotime analysis failure\\n")
                write.csv(data.frame(), "trajectory_dependent_genes.csv", row.names = FALSE)
            }
        }, error = function(e) {
            cat("Trajectory gene analysis failed:", e\$message, "\\n")
            write.csv(data.frame(), "trajectory_dependent_genes.csv", row.names = FALSE)
        })
        
        # Save Monocle3 object
        saveRDS(cds, "monocle3_cds_object.rds")
        
        # Generate analysis summary
        monocle3_summary <- list(
            n_cells = ncol(cds),
            n_genes = nrow(cds),
            n_partitions = length(unique(cds@clusters[["UMAP"]]\$partitions)),
            n_clusters = length(unique(cds@clusters[["UMAP"]]\$clusters)),
            pseudotime_success = pseudotime_success,
            trajectory_genes_success = trajectory_genes_success
        )
        
        # If pseudotime information is available
        if(pseudotime_success && "pseudotime" %in% names(colData(cds))) {
            monocle3_summary\$pseudotime_range <- c(min(pseudotime(cds), na.rm = TRUE), 
                                                   max(pseudotime(cds), na.rm = TRUE))
            monocle3_summary\$cells_with_pseudotime <- sum(!is.infinite(pseudotime(cds)))
        }
        
        # Save summary
        summary_lines <- c(
            "=== Monocle3 Analysis Summary ===",
            paste("Number of analyzed cells:", monocle3_summary\$n_cells),
            paste("Number of analyzed genes:", monocle3_summary\$n_genes),
            paste("Number of cell partitions:", monocle3_summary\$n_partitions),
            paste("Number of clusters:", monocle3_summary\$n_clusters),
            paste("Grouping method:", group_column),
            paste("Pseudotime analysis:", ifelse(pseudotime_success, "Success", "Failed")),
            paste("Trajectory gene analysis:", ifelse(trajectory_genes_success, "Success", "Failed"))
        )
        
        if(pseudotime_success && "pseudotime_range" %in% names(monocle3_summary)) {
            summary_lines <- c(summary_lines,
                paste("Pseudotime range:", 
                      round(monocle3_summary\$pseudotime_range[1], 2), "-", 
                      round(monocle3_summary\$pseudotime_range[2], 2)),
                paste("Cells with pseudotime:", monocle3_summary\$cells_with_pseudotime))
        }
        
        writeLines(summary_lines, "monocle3_summary.txt")
        
        cat("\\n=== Monocle3 Analysis Summary ===\\n")
        cat("Number of analyzed cells:", monocle3_summary\$n_cells, "\\n")
        cat("Number of analyzed genes:", monocle3_summary\$n_genes, "\\n")
        cat("Number of cell partitions:", monocle3_summary\$n_partitions, "\\n")
        cat("Number of clusters:", monocle3_summary\$n_clusters, "\\n")
        
        cat("\\nMonocle3 analysis completed!\\n")
        
    }, error = function(e) {
        cat("Monocle3 analysis error:", e\$message, "\\n")
        
        # Create empty output files
        empty_cds <- list(data = matrix(0, nrow = 1, ncol = 1))
        saveRDS(empty_cds, "monocle3_cds_object.rds")
        
        write.csv(data.frame(), "trajectory_dependent_genes.csv", row.names = FALSE)
        
        png("monocle3_analysis_error.png", width = 8, height = 6, units = "in", res = 300)
        plot.new()
        text(0.5, 0.5, paste("Monocle3 analysis failed:\\n", e\$message), 
             cex = 1.2, col = "red")
        dev.off()
        
        writeLines(c(
            "=== Monocle3 Analysis Error ===",
            paste("Error message:", e\$message),
            "Possible causes:",
            "1. monocle3 or SeuratWrappers packages not properly installed",
            "2. Insufficient data dimensions for trajectory analysis",
            "3. Insufficient memory",
            "Please check error message and rerun"
        ), "monocle3_summary.txt")
    })
    """
}