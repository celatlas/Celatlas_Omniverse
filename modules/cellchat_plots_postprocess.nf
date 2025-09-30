/*
 * CellChat Plots Post-processing Module
 * Generate heatmaps and bubble plots from saved CellChat objects at the end of the workflow
 */

process CELLCHAT_PLOTS_POSTPROCESS {
    publishDir "${params.output_dir}/5.Cellchat", mode: 'copy'
    
    input:
    path cellchat_rdata
    
    output:
    path "03_number_of_interactions_heatmap.pdf", optional: true
    path "03_number_of_interactions_heatmap.png", optional: true
    path "05_ligand_receptor_bubble.pdf", optional: true
    path "05_ligand_receptor_bubble.png", optional: true
    
    script:
    """
    #!/usr/bin/Rscript
    
    # Set R package paths
    .libPaths(c("/home/bioinfo/R/library", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library"))
    
    # Load required packages
    suppressMessages({
        library(CellChat)
        library(ggplot2)
        library(ComplexHeatmap)
        library(circlize)
        library(grDevices)
        library(grid)
        library(RColorBrewer)
    })
    
    # Configure graphics environment
    options(bitmapType = 'cairo')
    
    cat("=== CellChat Plots Post-processing ===\\n")
    cat("Starting to generate plots from saved CellChat objects...\\n")
    
    # Load CellChat results
    tryCatch({
        load("${cellchat_rdata}")
        cat("Successfully loaded CellChat object\\n")
        
        # Check if CellChat object is valid
        if(exists("cellchat") && !is.null(cellchat@net) && !is.null(cellchat@net\$count)) {
            cat("CellChat object validation passed\\n")
            
            # Check communication data
            if(sum(cellchat@net\$count) > 0) {
                cat("Found valid cell communication data\\n")
                
                # === Generate interaction heatmap ===
                cat("Generating interaction count heatmap...\\n")
                
                # Use ComplexHeatmap directly to create more stable heatmaps
                tryCatch({
                    count_mat <- cellchat@net\$count
                    weight_mat <- cellchat@net\$weight

                    plot_mat <- if(!is.null(count_mat) && sum(count_mat) > 0) count_mat else weight_mat
                    
                    if(!is.null(plot_mat) && nrow(plot_mat) > 1 && ncol(plot_mat) > 1 && sum(plot_mat) > 0) {
                        
                        cell_types <- rownames(plot_mat)
                        n_types <- length(cell_types)

                        colors <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", 
                                   "#ffff33", "#a65628", "#f781bf", "#999999", "#66c2a5",
                                   "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f",
                                   "#e5c494", "#b3b3b3", "#1b9e77", "#d95f02", "#7570b3")

                        cell_colors <- setNames(colors[1:n_types], cell_types)

                        row_sums <- rowSums(plot_mat)
                        col_sums <- colSums(plot_mat)

                        row_anno <- rowAnnotation(
                            CellType = cell_types,
                            col = list(CellType = cell_colors),
                            annotation_name_side = "bottom",
                            annotation_name_gp = gpar(fontsize = 0), 
                            simple_anno_size = unit(0.5, "cm"),
                            show_legend = FALSE
                        )

                        col_anno <- HeatmapAnnotation(
                            CellType = cell_types,
                            col = list(CellType = cell_colors),
                            annotation_name_side = "left",
                            annotation_name_gp = gpar(fontsize = 0),
                            simple_anno_size = unit(0.5, "cm"),
                            show_legend = FALSE
                        )

                        top_anno <- HeatmapAnnotation(
                            ColSums = anno_barplot(
                                col_sums,
                                bar_width = 0.6, 
                                border = FALSE,
                                gp = gpar(fill = colors[1:n_types], col = NA),
                                axis_param = list(
                                    side = "left",
                                    labels_rot = 0,
                                    gp = gpar(fontsize = 12, fontface = "bold", lwd = 2),  
                                    direction = "normal"
                                ),
                                height = unit(3, "cm")
                            ),
                            annotation_name_side = "left",
                            annotation_name_gp = gpar(fontsize = 0), 
                            gap = unit(4, "mm")
                        )

                        right_anno <- rowAnnotation(
                            RowSums = anno_barplot(
                                row_sums,
                                bar_width = 0.6, 
                                border = FALSE,
                                gp = gpar(fill = colors[1:n_types], col = NA),
                                axis_param = list(
                                    side = "bottom",
                                    labels_rot = 0,
                                    gp = gpar(fontsize = 12, fontface = "bold", lwd = 2),  
                                    direction = "normal"
                                ),
                                width = unit(3, "cm")
                            ),
                            annotation_name_side = "bottom",
                            annotation_name_gp = gpar(fontsize = 0), 
                            gap = unit(4, "mm")
                        )
                        

                        max_val <- max(plot_mat)
                        col_fun <- colorRamp2(
                            c(0, max_val * 0.3, max_val * 0.6, max_val), 
                            c("white", "#fee5d9", "#fc9272", "#de2d26")
                        )

                        ht <- Heatmap(
                            plot_mat,
                            name = "Number of interactions",
                            col = col_fun,
                            cluster_rows = FALSE,
                            cluster_columns = FALSE,
                            show_row_names = TRUE,
                            show_column_names = TRUE,
                            row_names_gp = gpar(fontsize = 14, fontface = "bold", col = "black"),
                            column_names_gp = gpar(fontsize = 14, fontface = "bold", col = "black"),
                            row_names_side = "left",
                            column_names_side = "bottom",
                            column_names_rot = 45, 
                            column_title = NULL,
                            row_title = NULL,
                            heatmap_legend_param = list(
                                title_position = "leftcenter-rot",
                                title_gp = gpar(fontsize = 14, fontface = "bold"),
                                labels_gp = gpar(fontsize = 12, fontface = "bold"),
                                legend_direction = "vertical",
                                legend_height = unit(4, "cm"),
                                legend_width = unit(1, "cm")
                            ),
                            left_annotation = row_anno,
                            bottom_annotation = col_anno,
                            top_annotation = top_anno,
                            right_annotation = right_anno
                        )

                        pdf("03_number_of_interactions_heatmap.pdf", width = 14, height = 12)
                        draw(ht, 
                             heatmap_legend_side = "right",
                             annotation_legend_side = "right",
                             gap = unit(c(6, 6), "mm")) 
                        dev.off()

                        png("03_number_of_interactions_heatmap.png", width = 14, height = 12, 
                            units = "in", res = 300, type = "cairo")
                        draw(ht, 
                             heatmap_legend_side = "right",
                             annotation_legend_side = "right", 
                             gap = unit(c(6, 6), "mm")) 
                        dev.off()
                        
                        cat("Custom ComplexHeatmap heatmap saved successfully\\n")
                        cat("Number of cell types:", n_types, "\\n")
                        
                    } else if(!is.null(plot_mat) && (nrow(plot_mat) == 1 || ncol(plot_mat) == 1)) {
                        # Matrix too small, create simple bar plot
                        library(ggplot2)
                        
                        if(nrow(plot_mat) == 1) {
                            # Single row: convert to bar plot
                            df <- data.frame(
                                Target = colnames(plot_mat),
                                Value = as.vector(plot_mat[1,]),
                                stringsAsFactors = FALSE
                            )
                            p_simple <- ggplot(df, aes(x = Target, y = Value)) +
                                geom_col(fill = "lightblue", color = "darkblue") +
                                theme_minimal() +
                                ggtitle("Cell Communication (Single Source)") +
                                theme(axis.text.x = element_text(angle = 45, hjust = 1))
                        } else {
                            # Single column: convert to bar plot
                            df <- data.frame(
                                Source = rownames(plot_mat),
                                Value = as.vector(plot_mat[,1]),
                                stringsAsFactors = FALSE
                            )
                            p_simple <- ggplot(df, aes(x = Source, y = Value)) +
                                geom_col(fill = "lightblue", color = "darkblue") +
                                theme_minimal() +
                                ggtitle("Cell Communication (Single Target)") +
                                theme(axis.text.x = element_text(angle = 45, hjust = 1))
                        }
                        
                        ggsave("03_number_of_interactions_heatmap.pdf", plot = p_simple, 
                               width = 10, height = 8, device = "pdf")
                        ggsave("03_number_of_interactions_heatmap.png", plot = p_simple, 
                               width = 10, height = 8, dpi = 300, device = "png")
                        
                        cat("Simple bar plot saved successfully\\n")
                        
                    } else {
                        # No valid data
                        library(ggplot2)
                        no_data_plot <- ggplot() + 
                            geom_text(aes(x = 0.5, y = 0.5), 
                                     label = "No interaction data available", 
                                     size = 6, color = "orange") +
                            theme_void() +
                            ggtitle("Cell Communication Heatmap")
                        
                        ggsave("03_number_of_interactions_heatmap.pdf", plot = no_data_plot, 
                               width = 10, height = 8, device = "pdf")
                        ggsave("03_number_of_interactions_heatmap.png", plot = no_data_plot, 
                               width = 10, height = 8, dpi = 300, device = "png")
                        
                        cat("No data notification plot saved successfully\\n")
                    }
                    
                }, error = function(e) {
                    cat("Heatmap generation completely failed:", e\$message, "\\n")
                    
                    # Final fallback: simple error plot
                    library(ggplot2)
                    error_plot <- ggplot() + 
                        geom_text(aes(x = 0.5, y = 0.5), 
                                 label = paste("Heatmap generation failed\\n", e\$message), 
                                 size = 5, color = "red") +
                        theme_void() +
                        ggtitle("Heatmap Generation Error")
                    
                    ggsave("03_number_of_interactions_heatmap.pdf", plot = error_plot, 
                           width = 10, height = 8, device = "pdf")
                    ggsave("03_number_of_interactions_heatmap.png", plot = error_plot, 
                           width = 10, height = 8, dpi = 300, device = "png")
                    
                    cat("Error notification plot saved successfully\\n")
                })
                
                # === Generate bubble plot ===
                cat("Generating ligand-receptor bubble plot...\\n")
                
                tryCatch({
                    # Generate bubble plot with ggplot object assignment
                    p_bubble <- netVisual_bubble(cellchat, remove.isolate = FALSE)
                    
                    # Save with ggsave
                    if(!is.null(p_bubble)) {
                        ggsave("05_ligand_receptor_bubble.pdf", plot = p_bubble, 
                               width = 12, height = 15, device = "pdf")
                        ggsave("05_ligand_receptor_bubble.png", plot = p_bubble, 
                               width = 12, height = 15, dpi = 200, device = "png")
                        
                        cat("Bubble plot saved successfully\\n")
                    }
                    
                }, error = function(e) {
                    cat("Bubble plot generation failed:", e\$message, "\\n")
                    
                    # Create fallback plot
                    tryCatch({
                        library(ggplot2)
                        
                        fallback_plot <- ggplot() + 
                            geom_text(aes(x = 0.5, y = 0.5), 
                                     label = paste("The bubble chart generation failed\\n", e\$message), 
                                     size = 5, color = "red") +
                            theme_void() +
                            ggtitle("CellChat Bubble Plot Error")
                        
                        ggsave("05_ligand_receptor_bubble.pdf", plot = fallback_plot, 
                               width = 12, height = 15, device = "pdf")
                        ggsave("05_ligand_receptor_bubble.png", plot = fallback_plot, 
                               width = 12, height = 15, dpi = 200, device = "png")
                        
                        cat("The error prompt image has been saved successfully\\n")
                        
                    }, error = function(e2) {
                        cat("Even the error prompt image failed to be saved:", e2\$message, "\\n")
                    })
                })
                
            } else {
                cat("No valid cell communication data detected\\n")
                
                # Create "no data" plots
                library(ggplot2)
                no_data_plot <- ggplot() + 
                    geom_text(aes(x = 0.5, y = 0.5), 
                             label = "There is no valid cell communication data\\nNo valid cell communication data", 
                             size = 6, color = "orange") +
                    theme_void() +
                    ggtitle("CellChat Analysis")
                
                ggsave("03_number_of_interactions_heatmap.pdf", plot = no_data_plot, 
                       width = 10, height = 8, device = "pdf")
                ggsave("03_number_of_interactions_heatmap.png", plot = no_data_plot, 
                       width = 10, height = 8, dpi = 200, device = "png")
                       
                ggsave("05_ligand_receptor_bubble.pdf", plot = no_data_plot, 
                       width = 12, height = 15, device = "pdf")
                ggsave("05_ligand_receptor_bubble.png", plot = no_data_plot, 
                       width = 12, height = 15, dpi = 200, device = "png")
            }
            
        } else {
            cat("CellChat object is invalid or empty\\n")
            
            # Create error plots
            library(ggplot2)
            error_plot <- ggplot() + 
                geom_text(aes(x = 0.5, y = 0.5), 
                         label = "The CellChat object is invalid\nInvalid CellChat object", 
                         size = 6, color = "red") +
                theme_void() +
                ggtitle("CellChat Analysis Error")
            
            ggsave("03_number_of_interactions_heatmap.pdf", plot = error_plot, 
                   width = 10, height = 8, device = "pdf")
            ggsave("03_number_of_interactions_heatmap.png", plot = error_plot, 
                   width = 10, height = 8, dpi = 200, device = "png")
                   
            ggsave("05_ligand_receptor_bubble.pdf", plot = error_plot, 
                   width = 12, height = 15, device = "pdf")
            ggsave("05_ligand_receptor_bubble.png", plot = error_plot, 
                   width = 12, height = 15, dpi = 200, device = "png")
        }
        
    }, error = function(e) {
        cat("Failed to load CellChat object:", e\$message, "\\n")
        
        # Create loading error plots
        library(ggplot2)
        load_error_plot <- ggplot() + 
            geom_text(aes(x = 0.5, y = 0.5), 
                     label = paste("加载失败:", e\$message, "\\nLoading failed"), 
                     size = 5, color = "red") +
            theme_void() +
            ggtitle("CellChat Loading Error")
        
        ggsave("03_number_of_interactions_heatmap.pdf", plot = load_error_plot, 
               width = 10, height = 8, device = "pdf")
        ggsave("03_number_of_interactions_heatmap.png", plot = load_error_plot, 
               width = 10, height = 8, dpi = 200, device = "png")
               
        ggsave("05_ligand_receptor_bubble.pdf", plot = load_error_plot, 
               width = 12, height = 15, device = "pdf")
        ggsave("05_ligand_receptor_bubble.png", plot = load_error_plot, 
               width = 12, height = 15, dpi = 200, device = "png")
    })
    
    # Verify generated files
    cat("\\n=== Verify Generated Files ===\\n")
    files_to_check <- c(
        "03_number_of_interactions_heatmap.pdf",
        "03_number_of_interactions_heatmap.png",
        "05_ligand_receptor_bubble.pdf",
        "05_ligand_receptor_bubble.png"
    )
    
    for(file in files_to_check) {
        if(file.exists(file)) {
            file_size <- file.info(file)\$size
            cat("File:", file, "Size:", file_size, "bytes\\n")
            if(file_size > 1000) {
                cat("  -> File normal\\n")
            } else {
                cat("  -> Warning: File may be abnormal (too small)\\n")
            }
        } else {
            cat("File:", file, "not found\\n")
        }
    }
    
    cat("The post-processing of the CellChat chart has been completed\\n")
    """
}