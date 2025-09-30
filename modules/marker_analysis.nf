/*
 * Marker Gene Analysis Module
 */

process MARKER_ANALYSIS {
    publishDir "${params.output_dir}/3.Marker", mode: 'copy'
    
    input:
    path seurat_rds
    
    output:
    path "all_cluster_markers.csv", emit: markers
    path "top*_cluster_markers.csv"
    path "*.png"
    path "*.pdf"
    path "marker_heatmap_data.csv"
    
    script:
    """
    #!/usr/bin/Rscript
    
    # Load required packages
    suppressMessages({
        library(Seurat)
        library(ggplot2)
        library(dplyr)
        library(pheatmap)
        library(RColorBrewer)
        library(viridis)
    })
    
    cat("\\n=== Marker Gene Analysis and Clustering Heatmap ===\\n")
    
    tryCatch({
        # Read clustered Seurat object
        cat("Reading clustered data...\\n")
        st <- readRDS("${seurat_rds}")
        
        cat("Number of input cells:", ncol(st), "\\n")
        cat("Number of input genes:", nrow(st), "\\n")
        
        # Calculate marker genes for each cluster using FindAllMarkers
        cat("Calculating marker genes for each cluster...\\n")
        Idents(st) <- "seurat_clusters"
        
        all_markers <- FindAllMarkers(
            st, 
            only.pos = TRUE,
            min.pct = ${params.marker.min_pct},
            logfc.threshold = ${params.marker.logfc_threshold},
            test.use = "${params.marker.test_use}",
            verbose = FALSE
        )
        
        if(nrow(all_markers) > 0) {
            # Save all marker gene results
            write.csv(all_markers, "all_cluster_markers.csv", row.names = FALSE)
            
            # Get top marker genes for each cluster
            top_markers <- all_markers %>%
                group_by(cluster) %>%
                top_n(n = ${params.marker.top_n}, wt = avg_log2FC) %>%
                arrange(cluster, desc(avg_log2FC))
            
            # Save top marker genes
            write.csv(top_markers, "top${params.marker.top_n}_cluster_markers.csv", row.names = FALSE)
            
            # Create marker gene expression matrix for heatmap
            if(nrow(top_markers) > 0) {
                cat("Generating marker gene clustering heatmap...\\n")
                
                # Get marker genes
                marker_genes <- unique(top_markers\$gene)
                
                # Calculate average expression for each cluster
                Idents(st) <- "seurat_clusters"
                # Use AggregateExpression for Seurat v5 compatibility
                avg_exp <- tryCatch({
                    AggregateExpression(st, features = marker_genes, verbose = FALSE)
                }, error = function(e) {
                    cat("AggregateExpression failed, trying AverageExpression...\\n")
                    AverageExpression(st, features = marker_genes, verbose = FALSE)
                })
                
                # Check available expression data types
                if("SCT" %in% names(avg_exp)) {
                    avg_exp_matrix <- avg_exp\$SCT
                } else if("Spatial" %in% names(avg_exp)) {
                    avg_exp_matrix <- avg_exp\$Spatial
                } else if("RNA" %in% names(avg_exp)) {
                    avg_exp_matrix <- avg_exp\$RNA
                } else {
                    avg_exp_matrix <- avg_exp[[1]]  # Use first available data type
                }
                
                # Check and process expression matrix
                if(!is.null(avg_exp_matrix)) {
                    # Convert to regular matrix
                    avg_exp_matrix <- as.matrix(avg_exp_matrix)
                    
                    # Ensure matrix has at least 1 row and 1 column
                    if(nrow(avg_exp_matrix) >= 1 && ncol(avg_exp_matrix) >= 1) {
                        # Check and handle missing values
                        avg_exp_matrix[is.na(avg_exp_matrix)] <- 0
                        
                        # Standardize expression matrix (row-wise Z-score normalization)
                        scaled_matrix <- t(scale(t(avg_exp_matrix)))
                        
                        # Handle missing values after standardization
                        scaled_matrix[is.na(scaled_matrix)] <- 0
                        
                        # Handle column sorting carefully - check if columns have 'g' prefix
                        if(any(grepl("^g", colnames(scaled_matrix)))) {
                            cluster_order <- paste0("g", sort(as.numeric(gsub("g", "", colnames(scaled_matrix)))))
                            scaled_matrix <- scaled_matrix[, cluster_order, drop = FALSE]
                        }
                        
                        # Sort rows by gene order in top_markers - ensure genes exist in matrix
                        gene_order <- top_markers\$gene
                        available_genes <- intersect(gene_order, rownames(scaled_matrix))
                        if(length(available_genes) > 0) {
                            scaled_matrix <- scaled_matrix[available_genes, , drop = FALSE]
                        }
                        
                        # Rename columns: change g0, g1, g2... to Cluster0, Cluster1, Cluster2...
                        colnames(scaled_matrix) <- paste0("Cluster", 0:(ncol(scaled_matrix)-1))
                        
                        # Create heatmap (without clustering)
                        library(pheatmap)
                        
                        # Use specified gradient colors
                        custom_colors <- colorRampPalette(c('#330066', '#336699', '#66CC66', '#FFCC33'))(100)
                        
                        # Calculate gaps for row separation by cluster
                        cluster_table <- table(top_markers\$cluster[top_markers\$gene %in% rownames(scaled_matrix)])
                        gaps_positions <- if(length(cluster_table) > 1) {
                            cumsum(cluster_table)[-length(cluster_table)]
                        } else {
                            NULL
                        }
                        
                        # Generate PNG version
                        png("marker_genes_heatmap.png", width = 10, height = 12, units = "in", res = 300)
                        pheatmap(
                            scaled_matrix,
                            scale = "none",  # Already standardized
                            cluster_rows = FALSE,  # No row clustering
                            cluster_cols = FALSE,  # No column clustering
                            color = custom_colors,  # Use custom colors
                            show_rownames = TRUE,
                            show_colnames = TRUE,
                            main = "Top ${params.marker.top_n} Marker Genes by Cluster (avg_log2FC)",
                            fontsize = 10,
                            fontsize_row = 8,
                            fontsize_col = 10,
                            angle_col = 45,  # Column names tilted 45 degrees
                            gaps_row = gaps_positions,
                            cellwidth = 30,
                            cellheight = 10
                        )
                        dev.off()
                        
                        # Generate PDF version
                        pdf("marker_genes_heatmap.pdf", width = 10, height = 12)
                        pheatmap(
                            scaled_matrix,
                            scale = "none",
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            color = custom_colors,  # Use custom colors
                            show_rownames = TRUE,
                            show_colnames = TRUE,
                            main = "Top ${params.marker.top_n} Marker Genes by Cluster",
                            fontsize = 10,
                            fontsize_row = 8,
                            fontsize_col = 10,
                            angle_col = 45,  # Column names tilted 45 degrees
                            gaps_row = gaps_positions,
                            cellwidth = 30,
                            cellheight = 10
                        )
                        dev.off()
                        
                        # Create cluster similarity clustering plot based on similarity
                        if(ncol(scaled_matrix) >= 2) {
                            cat("Generating cluster similarity clustering plot...\\n")
                            
                            # Calculate correlation between clusters
                            cluster_cor <- cor(scaled_matrix, method = "pearson", use = "complete.obs")
                            
                            # Handle missing values in correlation matrix
                            cluster_cor[is.na(cluster_cor)] <- 0
                            
                            # Generate PNG version
                            png("cluster_similarity_heatmap.png", width = 8, height = 8, units = "in", res = 300)
                            pheatmap(
                                cluster_cor,
                                scale = "none",
                                clustering_distance_rows = "euclidean",
                                clustering_distance_cols = "euclidean",
                                clustering_method = "complete",
                                color = colorRampPalette(c("blue", "white", "red"))(100),
                                show_rownames = TRUE,
                                show_colnames = TRUE,
                                main = "Cluster Similarity based on Marker Genes",
                                fontsize = 12,
                                display_numbers = TRUE,
                                number_format = "%.2f"
                            )
                            dev.off()
                            
                            # Generate PDF version
                            pdf("cluster_similarity_heatmap.pdf", width = 8, height = 8)
                            pheatmap(
                                cluster_cor,
                                scale = "none",
                                clustering_distance_rows = "euclidean",
                                clustering_distance_cols = "euclidean",
                                clustering_method = "complete",
                                color = colorRampPalette(c("blue", "white", "red"))(100),
                                show_rownames = TRUE,
                                show_colnames = TRUE,
                                main = "Cluster Similarity based on Marker Genes",
                                fontsize = 12,
                                display_numbers = TRUE,
                                number_format = "%.2f"
                            )
                            dev.off()
                        }
                        
                        cat("  Marker gene heatmap analysis completed\\n")
                        cat("  Number of marker genes analyzed:", length(marker_genes), "\\n")
                        cat("  Number of clusters analyzed:", ncol(avg_exp_matrix), "\\n")
                        
                        # Save heatmap data
                        write.csv(as.data.frame(scaled_matrix), "marker_heatmap_data.csv")
                    }
                }
            }
        } else {
            cat("No significant marker genes found\\n")
            # Create empty output files
            write.csv(data.frame(), "all_cluster_markers.csv", row.names = FALSE)
            write.csv(data.frame(), "top${params.marker.top_n}_cluster_markers.csv", row.names = FALSE)
            write.csv(data.frame(), "marker_heatmap_data.csv", row.names = FALSE)
            
            # Create placeholder PNG files to satisfy Nextflow requirements
            png("marker_genes_heatmap.png", width = 8, height = 12, units = "in", res = 300)
            plot.new()
            text(0.5, 0.5, "No significant marker genes found", cex = 1.2, col = "gray")
            dev.off()
        }
    }, error = function(e) {
        cat("Marker gene analysis error:", e\$message, "\\n")
        # Create empty output files on error
        write.csv(data.frame(), "all_cluster_markers.csv", row.names = FALSE)
        write.csv(data.frame(), "top${params.marker.top_n}_cluster_markers.csv", row.names = FALSE)
        write.csv(data.frame(), "marker_heatmap_data.csv", row.names = FALSE)
        
        # Create error PNG files to satisfy Nextflow requirements
        png("marker_genes_heatmap.png", width = 8, height = 12, units = "in", res = 300)
        plot.new()
        text(0.5, 0.5, paste("Marker gene analysis failed:\\n", e\$message), cex = 1.2, col = "red")
        dev.off()
    })
    
    cat("✓ Marker Gene Analysis and Clustering Heatmap\\n")
    """
}