/*
 * SingleR Cell Type Annotation Module
 */

process SINGLER_ANNOTATION {
    publishDir "${params.output_dir}/4.SingleR", mode: 'copy'
    
    input:
    path seurat_rds
    
    output:
    path "seurat_annotated.rds", emit: seurat_annotated
    path "*.png"
    path "*.pdf"
    path "singler_results.csv"
    path "celltype_cluster_distribution.csv"
    path "cluster_annotation_summary.csv"
    path "complete_metadata_with_annotation.csv"
    path "02_st_final_singleR.Rdata"
    
    script:
    """
    #!/usr/bin/Rscript
    
    # Load required packages
    suppressMessages({
        library(Seurat)
        library(SingleR)
        library(celldex)
        library(BiocParallel)
        library(ggplot2)
        library(dplyr)
        library(pheatmap)
        library(RColorBrewer)
        library(patchwork)
        library(matrixStats)
    })
    
    # Use SingleR for automatic annotation - optimized version
    cat("Starting SingleR cell type annotation...\\n")
    
    tryCatch({
        # Read clustered Seurat object
        cat("Reading clustered data...\\n")
        st <- readRDS("${seurat_rds}")
        
        cat("Number of input cells:", ncol(st), "\\n")
        cat("Number of clusters:", length(unique(st\$seurat_clusters)), "\\n")
        
        # Smart reference dataset selection - Prioritize the use of specified parameters; otherwise, automatic detection will be adopted
        if("${params.species}" == "human") {
            cat("Using parameter-specified species: human\\n")
            ref <- celldex::HumanPrimaryCellAtlasData()
            species_detected <- "human"
            cat("Using HumanPrimaryCellAtlas as reference dataset\\n")
        } else if("${params.species}" == "mouse" || "${params.species}" == "mmu") {
            cat("Using parameter-specified species: mouse (mmu)\\n")
            ref <- celldex::MouseRNAseqData()
            species_detected <- "mouse"
            cat("Using MouseRNAseq as reference dataset\\n")
        } else {

            cat("Auto-detecting dataset species...\\n")
            
            # Check gene name format to determine species
            test_genes <- head(rownames(st), 100)
            human_genes <- sum(grepl("^[A-Z][A-Z0-9-]*\$", test_genes))
            mouse_genes <- sum(grepl("^[A-Z][a-z0-9-]*\$", test_genes))
            
            if(human_genes > mouse_genes) {
                cat("Detected human data, loading human reference dataset...\\n")
                ref <- celldex::HumanPrimaryCellAtlasData()
                species_detected <- "human"
                cat("Using HumanPrimaryCellAtlas as reference dataset\\n")
            } else {
                cat("Detected mouse data, loading mouse reference dataset...\\n")
                ref <- celldex::MouseRNAseqData()
                species_detected <- "mouse"
                cat("Using MouseRNAseq as reference dataset\\n")
            }
        }
        
        # Data preprocessing - get normalized data
        cat("Extracting normalized data for annotation...\\n")
        
        # Check if data slot is empty, use counts if so
        tryCatch({
            sc_data <- GetAssayData(st, assay = "Spatial", layer = "data")
            if(nrow(sc_data) == 0) {
                cat("data layer is empty, using counts data...\\n")
                sc_data <- GetAssayData(st, assay = "Spatial", layer = "counts")
            }
        }, error = function(e) {
            cat("Using fallback method for data extraction...\\n")
            sc_data <- GetAssayData(st, assay = "Spatial")
        })
        
        # Ensure data matches
        sc_data <- sc_data[, colnames(st)]
        
        # Check gene name format and find common genes
        cat("Checking gene name matching...\\n")
        test_genes <- rownames(sc_data)
        ref_genes <- rownames(ref)
        
        # Smart gene matching strategy
        common_genes <- intersect(test_genes, ref_genes)
        cat("Initial common genes count:", length(common_genes), "\\n")
        
        if(length(common_genes) < 500) {
            cat("Few common genes, trying gene name conversion...\\n")
            
            if(species_detected == "mouse") {
                # Mouse gene name conversion: Atp5a1 -> ATP5A1
                test_genes_converted <- toupper(test_genes)
            } else {
                # Human gene name conversion: ATP5A1 -> Atp5a1 or ATP5A1 -> ATP5A1
                test_genes_converted <- paste0(toupper(substr(test_genes, 1, 1)), 
                                              tolower(substr(test_genes, 2, nchar(test_genes))))
            }
            
            common_genes_converted <- intersect(test_genes_converted, ref_genes)
            cat("Common genes after conversion:", length(common_genes_converted), "\\n")
            
            if(length(common_genes_converted) > length(common_genes)) {
                rownames(sc_data) <- test_genes_converted
                common_genes <- common_genes_converted
                cat("Using converted gene names\\n")
            }
            
            # If still too few, try other conversion strategies
            if(length(common_genes) < 300) {
                # Try removing gene version numbers
                test_genes_clean <- gsub("\\\\..*\$", "", test_genes)
                common_genes_clean <- intersect(test_genes_clean, ref_genes)
                
                if(length(common_genes_clean) > length(common_genes)) {
                    rownames(sc_data) <- test_genes_clean
                    common_genes <- common_genes_clean
                    cat("Using cleaned gene names, common genes count:", length(common_genes), "\\n")
                }
            }
        }
        
        # Keep only common genes
        sc_data <- sc_data[common_genes, ]
        cat("Final genes for annotation:", nrow(sc_data), "\\n")
        
        cat("Data dimensions:", nrow(sc_data), "genes x", ncol(sc_data), "spots\\n")
        
        # Convert to regular matrix if needed
        if(inherits(sc_data, "dgCMatrix")) {
            sc_data <- as.matrix(sc_data)
            cat("Converted to regular matrix\\n")
        }
        
        # Check if data is valid
        if(is.null(sc_data) || nrow(sc_data) == 0 || ncol(sc_data) == 0) {
            stop("Cannot get valid expression data for annotation")
        }
        
        # Set random seed
        set.seed(66)
        
        # Perform SingleR annotation (based on cluster)
        cat("Starting SingleR annotation...\\n")
        clusters <- as.character(st\$seurat_clusters)
        
        # Check if clusters are valid
        if(length(clusters) == 0 || all(is.na(clusters))) {
            stop("No clustering information found, please ensure clustering analysis is completed")
        }
        
        cat("Clustering info:", length(unique(clusters)), "clusters\\n")
        
        # Ensure ref and sc_data have enough genes
        if(nrow(sc_data) < 50) {
            stop("Too few common genes for accurate annotation")
        }
        
        # Perform SingleR annotation
        cat("Starting SingleR annotation...\\n")
        
        tryCatch({
            st.hesc <- SingleR(
                test = sc_data,
                ref = ref,
                labels = ref\$label.main,
                clusters = clusters,
                assay.type.test = "logcounts",
                BPPARAM = BiocParallel::SnowParam(workers = 2)
            )
            cat("SingleR annotation completed\\n")
        }, error = function(e) {
            cat("SingleR annotation error, trying without parallel processing...\\n")
            cat("Error message:", e\$message, "\\n")
            
            # Retry without parallel processing
            st.hesc <<- SingleR(
                test = sc_data,
                ref = ref,
                labels = ref\$label.main,
                clusters = clusters,
                assay.type.test = "logcounts"
            )
            cat("SingleR retry successful\\n")
        })
        
        # Result evaluation and filtering
        # Generate SingleR heatmap with column names displayed
        p_score_heatmap <- plotScoreHeatmap(st.hesc, show_colnames = TRUE, 
                                           fontsize_col = 10, angle_col = 45)
        
        # Save PDF
        pdf("10_SingleR_score_heatmap.pdf", width = 12, height = 8)
        print(p_score_heatmap)
        dev.off()
        
        # Save PNG
        png("10_SingleR_score_heatmap.png", width = 3600, height = 2400, res = 300)
        print(p_score_heatmap)
        dev.off()
        
        # Add annotation results to Seurat object
        st\$cell_type_main <- st.hesc\$labels[match(clusters, rownames(st.hesc))]
        
        # Simplified quality control - calculate annotation scores
        annotation_scores <- st.hesc\$scores[match(clusters, rownames(st.hesc)), ]
        if(is.matrix(annotation_scores)) {
            st\$annotation_score <- rowMaxs(annotation_scores, na.rm = TRUE)
        } else {
            st\$annotation_score <- annotation_scores
        }
        
        # Use simple fixed threshold for filtering
        confidence_threshold <- 0.5
        cat("Using confidence threshold:", confidence_threshold, "\\n")
        
        st\$cell_type_filtered <- ifelse(st\$annotation_score < confidence_threshold, 
                                        "Unknown", 
                                        st\$cell_type_main)
        
        # Improved visualization of annotation results
        p1 <- DimPlot(st, group.by = "cell_type_main", reduction = "umap", 
                      label = TRUE, repel = TRUE) + 
            ggtitle("Cell Types (SingleR)") +
            theme(text = element_text(family = "sans"), legend.position = "bottom") +
            guides(color = guide_legend(ncol = 3))
        
        # Save plots
        ggsave(plot = p1, "11_UMAP_CellType.png", width = 8, height = 8)
        ggsave(plot = p1, "11_UMAP_CellType.pdf", width = 8, height = 8)
        
        # Final clustering results plot
        p_res_comparison <- DimPlot(st, group.by = "seurat_clusters", reduction = "umap", 
                                   label = TRUE) +
            ggtitle("Final Clustering (Resolution 0.2)") +
            theme(text = element_text(family = "sans"))
        
        ggsave(plot = p_res_comparison, "11_Final_Clustering.png", width = 10, height = 8)
        ggsave(plot = p_res_comparison, "11_Final_Clustering.pdf", width = 10, height = 8)
        
        # Spatial distribution visualization
        p3 <- SpatialDimPlot(st, group.by = "cell_type_filtered", 
                             pt.size.factor = 0.1, label = TRUE) +
            ggtitle("Spatial Distribution of Cell Types") +
            theme(text = element_text(family = "sans"))
        
        p4 <- SpatialFeaturePlot(st, features = "annotation_score", 
                                 pt.size.factor = 0.1) +
            ggtitle("Annotation Confidence Scores") +
            theme(text = element_text(family = "sans"))
        
        p_spatial_combined <- wrap_plots(p3, p4, ncol = 2)
        ggsave(plot = p_spatial_combined, "12_Spatial_CellType.png", width = 16, height = 8)
        ggsave(plot = p_spatial_combined, "12_Spatial_CellType.pdf", width = 16, height = 8)
        
        # Save final annotation results
        save(st, st.hesc, file = "02_st_final_singleR.Rdata")
        
        # Save annotated Seurat object
        saveRDS(st, "seurat_annotated.rds")
        
        # Output annotation statistics table
        annotation_stats <- table(st\$cell_type_filtered, st\$seurat_clusters)
        write.csv(annotation_stats, "celltype_cluster_distribution.csv")
        
        # Generate cluster-cell type correspondence table
        cluster_annotation_summary <- st@meta.data %>%
            group_by(seurat_clusters, cell_type_filtered) %>%
            summarise(
                count = n(),
                mean_score = round(mean(annotation_score, na.rm = TRUE), 3),
                .groups = 'drop'
            ) %>%
            arrange(seurat_clusters)
        
        write.csv(cluster_annotation_summary, 
                  "cluster_annotation_summary.csv", 
                  row.names = FALSE)
        
        # Save quality control information
        write.csv(st@meta.data, 
                  "complete_metadata_with_annotation.csv", 
                  row.names = TRUE)
        
        # Save SingleR results
        write.csv(as.data.frame(st.hesc), "singler_results.csv", row.names = TRUE)
        
        # Save cell type statistics
        celltype_stats <- data.frame(
            cell_type = names(table(st\$cell_type_filtered)),
            count = as.numeric(table(st\$cell_type_filtered)),
            percentage = round(as.numeric(table(st\$cell_type_filtered)) / ncol(st) * 100, 2)
        )
        write.csv(celltype_stats, "celltype_stats.csv", row.names = FALSE)
        
        # Output final annotation result summary
        cat("\\n=== Cell Type Annotation Results Summary ===\\n")
        cat("Total cells:", ncol(st), "\\n")
        cat("Number of clusters:", length(unique(st\$seurat_clusters)), "\\n")
        cat("Number of annotated cell types:", length(unique(st\$cell_type_filtered)), "\\n")
        cat("Average annotation confidence:", round(mean(st\$annotation_score, na.rm = TRUE), 3), "\\n")
        
        type_counts <- table(st\$cell_type_filtered)
        cat("\\nCell type counts:\\n")
        for(i in 1:length(type_counts)) {
            cat(names(type_counts)[i], ":", type_counts[i], "\\n")
        }
        
        # Create annotation improvement suggestions
        cat("\\n=== Annotation Suggestions ===\\n")
        unknown_ratio <- sum(st\$cell_type_filtered == "Unknown") / ncol(st)
        if(unknown_ratio > 0.2) {
            cat("- High unknown cell ratio (", round(unknown_ratio*100, 1), "%), suggestions:\\n")
            cat("  1. Try lowering confidence threshold\\n")
            cat("  2. Use more reference datasets\\n")
            cat("  3. Check data quality and preprocessing steps\\n")
        }
        
        cat("\\nAnalysis completed!\\n")
        
    }, error = function(e) {
        cat("SingleR annotation error:", e\$message, "\\n")
        
        # Create empty output files on error
        write.csv(data.frame(), "singler_results.csv", row.names = FALSE)
        write.csv(data.frame(), "celltype_cluster_distribution.csv", row.names = FALSE)
        write.csv(data.frame(), "cluster_annotation_summary.csv", row.names = FALSE)
        write.csv(data.frame(), "complete_metadata_with_annotation.csv", row.names = FALSE)
        
        # Create error plots
        png("singler_annotation_error.png", width = 8, height = 6, units = "in", res = 300)
        plot.new()
        text(0.5, 0.5, paste("SingleR annotation failed:\\n", e\$message), cex = 1.2, col = "red")
        dev.off()
        
        # Save original object as annotated (with error info)
        st <- readRDS("${seurat_rds}")
        st\$cell_type_main <- "Error"
        st\$cell_type_filtered <- "Error"
        st\$annotation_score <- 0
        saveRDS(st, "seurat_annotated.rds")
        
        # Create minimal output files
        save(st, file = "02_st_final_singleR.Rdata")
        
        celltype_stats <- data.frame(
            cell_type = "Error",
            count = ncol(st),
            percentage = 100
        )
        write.csv(celltype_stats, "celltype_stats.csv", row.names = FALSE)
    })
    
    cat("✓ SingleR Cell Type Annotation\\n")
    """
}