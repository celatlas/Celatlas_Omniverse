/*
 * QC Analysis Module
 * Quality control and data loading
 */

process QC_ANALYSIS {
    publishDir "${params.output_dir}/1.QC", mode: 'copy'
    
    input:
    tuple path(matrix_dir), path(spatial_dir)
    
    output:
    path "seurat_qc.rds", emit: seurat_object
    path "*.png"
    path "*.pdf"  
    path "qc_metrics.csv"
    
    script:
    """
    #!/usr/bin/Rscript
    
    # Set R package paths
    .libPaths(c("/home/bioinfo/R/library", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library"))
    
    # Load required packages
    suppressMessages({
        library(Seurat)
        library(ggplot2)
        library(patchwork)
        library(dplyr)
        library(png)
        cat("Successfully loaded all R packages\n")
    })
    
    # Set output directory
    output_dir <- "."
    
    # Set point size factor for spatial plots
    pt_size_factor <- if("${params.spatial_plot_pt_size}" != "null") {
        as.numeric("${params.spatial_plot_pt_size}")
    } else {
        0.1  # default value
    }
    cat("Using point size factor:", pt_size_factor, "for spatial plots\\n")
    
    cat("=== QC Analysis and Data Loading ===\\n")
    
    # Read 10X data
    cat("Reading 10X data...\\n")
    counts <- Read10X(data.dir = "${matrix_dir}")
    st <- CreateSeuratObject(counts = counts, assay = "Spatial")
    
    cat("Number of cells:", length(colnames(st)), "\\n")
    cat("Number of genes:", length(rownames(st)), "\\n")
    
    # Load spatial image information
    cat("Loading spatial information...\\n")
    image <- Read10X_Image(
        image.dir = "${spatial_dir}",
        filter.matrix = TRUE,
        slice = "slice1"
    )
    
    image <- image[colnames(st)]
    st[["slice1"]] <- image
    
    # Read tissue positions file to get coordinate information
    tissue_pos_files <- c("tissue_positions_list.csv", "tissue_positions.csv")
    tissue_pos_file <- NULL
    for(file in tissue_pos_files) {
        if(file.exists(file.path("${spatial_dir}", file))) {
            tissue_pos_file <- file.path("${spatial_dir}", file)
            break
        }
    }
    
    # Get coordinate information - compatible with different Seurat versions
    coords <- NULL
    tryCatch({
        # Try new Seurat method
        coords <- GetTissueCoordinates(st, image = "slice1")
    }, error = function(e) {
        tryCatch({
            # Try old Seurat method
            coords <- st@images\$slice1@coordinates
        }, error = function(e2) {
            cat("Unable to get coordinate information, using default scale factors\\n")
        })
    })
    
    if(!is.null(tissue_pos_file) && !is.null(coords)) {
        # Read tissue positions
        tissue_pos <- read.csv(tissue_pos_file, header = FALSE)
        
        if(ncol(tissue_pos) >= 6) {
            # Standard Visium format: V1=barcode, V2=in_tissue, V3=array_row, V4=array_col, V5=pxl_row_in_fullres, V6=pxl_col_in_fullres
            image_coords_range <- max(c(max(tissue_pos\$V5) - min(tissue_pos\$V5), 
                                        max(tissue_pos\$V6) - min(tissue_pos\$V6)))
            
            # Adjust based on coords dataframe column names
            if("imagerow" %in% colnames(coords) && "imagecol" %in% colnames(coords)) {
                spot_coords_range <- max(c(max(coords\$imagerow) - min(coords\$imagerow),
                                           max(coords\$imagecol) - min(coords\$imagecol)))
            } else if("x" %in% colnames(coords) && "y" %in% colnames(coords)) {
                spot_coords_range <- max(c(max(coords\$x) - min(coords\$x),
                                           max(coords\$y) - min(coords\$y)))
            } else {
                # Use numeric columns
                numeric_cols <- sapply(coords, is.numeric)
                if(sum(numeric_cols) >= 2) {
                    coord_matrix <- as.matrix(coords[, numeric_cols])
                    spot_coords_range <- max(apply(coord_matrix, 2, function(x) max(x) - min(x)))
                } else {
                    spot_coords_range <- 1000  # Default value
                }
            }
            
            # Calculate actual scale factors
            if(spot_coords_range > 0 && image_coords_range > 0) {
                actual_hires_factor <- image_coords_range / spot_coords_range
                actual_lowres_factor <- actual_hires_factor * 0.3
            } else {
                actual_hires_factor <- 1.0
                actual_lowres_factor <- 0.3
            }
        } else {
            # Fallback: estimate based on tissue positions
            actual_hires_factor <- 1.0
            actual_lowres_factor <- 0.3
        }
    } else if(!is.null(coords)) {
        # Only use coordinate information for dynamic adjustment
        if("imagerow" %in% colnames(coords) && "imagecol" %in% colnames(coords)) {
            coord_range <- max(c(max(coords\$imagerow) - min(coords\$imagerow),
                                 max(coords\$imagecol) - min(coords\$imagecol)))
        } else if("x" %in% colnames(coords) && "y" %in% colnames(coords)) {
            coord_range <- max(c(max(coords\$x) - min(coords\$x),
                                 max(coords\$y) - min(coords\$y)))
        } else {
            coord_range <- 1000
        }
        
        if(coord_range > 2000) {
            actual_hires_factor <- coord_range / 2000
            actual_lowres_factor <- actual_hires_factor * 0.3
        } else {
            actual_hires_factor <- 1.0
            actual_lowres_factor <- 0.3
        }
    } else {
        # Final fallback: use default values
        cat("Using default scale factors\\n")
        actual_hires_factor <- 1.0
        actual_lowres_factor <- 0.3
    }
    
    # Apply calculated scale factors
    st@images\$slice1@scale.factors\$hires <- actual_hires_factor
    st@images\$slice1@scale.factors\$lowres <- actual_lowres_factor
    
    # Load tissue slice image if available
    if(file.exists(file.path("${spatial_dir}", "tissue_lowres_image.png"))) {
        st@images\$slice1@image <- png::readPNG(file.path("${spatial_dir}", "tissue_lowres_image.png"))
    }
    
    cat("Applied scale factors - hires:", actual_hires_factor, "lowres:", actual_lowres_factor, "\\n")
    
    # Calculate QC metrics
    cat("Calculating quality control metrics...\\n")
    st[["percent.mt"]] <- PercentageFeatureSet(st, pattern = "^MT-|^mt-")
    st[["percent.rb"]] <- PercentageFeatureSet(st, pattern = "^RP[SL]|^Rp[sl]")
    
    # QC metrics statistics
    qc_stats <- data.frame(
        Cell_ID = colnames(st),
        nFeature_Spatial = st\$nFeature_Spatial,
        nCount_Spatial = st\$nCount_Spatial,
        percent.mt = st\$percent.mt,
        percent.rb = st\$percent.rb
    )
    
    write.csv(qc_stats, "qc_metrics.csv", row.names = FALSE)
    
    # Create QC plots
    cat("Generating QC plots...\\n")
    
    # 1. QC metrics violin plots with custom colors
    features <- c("nFeature_Spatial", "nCount_Spatial", "percent.mt")
    qc_colors <- c("#f38181", "#a8d8ea", "#fce38a")

    p1_list <- lapply(seq_along(features), function(i) {
    VlnPlot(st, features = features[i], pt.size = 0.1, fill = qc_colors[i]) +
        ggtitle(features[i]) +
        theme(legend.position = "none")
    })

    p1 <- p1_list[[1]] | p1_list[[2]] | p1_list[[3]]
    
    # 2. Feature correlation scatter plots
    p2 <- FeatureScatter(st, feature1 = "nCount_Spatial", feature2 = "percent.mt") + 
        geom_hline(yintercept = ${params.max_mt_percent}, color = "red", linetype = "dashed")
    p3 <- FeatureScatter(st, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") +
        geom_hline(yintercept = c(${params.min_features}, ${params.max_features}), color = "red", linetype = "dashed")
    
    # 3. Spatial QC plots
    p4 <- SpatialFeaturePlot(st, features = "nCount_Spatial", pt.size.factor = pt_size_factor) + 
        ggtitle("Spatial Distribution of Total Counts") +
        theme(legend.position = "right", text = element_text(family = "sans"))
    p5 <- SpatialFeaturePlot(st, features = "nFeature_Spatial", pt.size.factor = pt_size_factor) + 
        ggtitle("Spatial Distribution of Gene Counts") +
        theme(legend.position = "right", text = element_text(family = "sans"))
    p6 <- SpatialFeaturePlot(st, features = "percent.mt", pt.size.factor = pt_size_factor) + 
        ggtitle("Spatial Distribution of Mitochondrial %") +
        theme(legend.position = "right", text = element_text(family = "sans"))
    
    # Save plots
    png("01_QC_violin_plots.png", width = 15, height = 5, units = "in", res = 300)
    print(p1)
    dev.off()
    
    png("02_QC_scatter_plots.png", width = 12, height = 5, units = "in", res = 300)
    print(p2 | p3)
    dev.off()
    
    png("03_spatial_QC_features.png", width = 18, height = 6, units = "in", res = 300)
    print(p4 | p5 | p6)
    dev.off()
    
    # Individual PDF files
    pdf("01_QC_violin_plots.pdf", width = 15, height = 5)
    print(p1)
    dev.off()
    
    pdf("02_QC_scatter_plots.pdf", width = 12, height = 5)
    print(p2 | p3)
    dev.off()
    
    pdf("03_spatial_QC_features.pdf", width = 18, height = 6)
    print(p4 | p5 | p6)
    dev.off()
    
    # Individual spatial plots - PNG files
    png("04_spatial_total_counts.png", width = 8, height = 6, units = "in", res = 300)
    print(p4)
    dev.off()
    
    png("05_spatial_gene_counts.png", width = 8, height = 6, units = "in", res = 300)
    print(p5)
    dev.off()
    
    # Individual spatial plots - PDF files
    pdf("04_spatial_total_counts.pdf", width = 8, height = 6)
    print(p4)
    dev.off()
    
    pdf("05_spatial_gene_counts.pdf", width = 8, height = 6)
    print(p5)
    dev.off()
    
    # Data filtering
    cat("Applying quality control filtering...\\n")
    cat("Cells before filtering:", ncol(st), "\\n")
    
    # Filter low-quality cells and genes - matching original script parameters
    st <- subset(st, subset = 
        nCount_Spatial > ${params.min_counts} & 
        nCount_Spatial < ${params.max_features} & 
        nFeature_Spatial > ${params.min_features} & 
        percent.mt < ${params.max_mt_percent}
    )
    
    # Filter low-expression genes
    counts <- GetAssayData(st, assay = "Spatial", layer = "counts")
    features_keep <- rownames(counts)[Matrix::rowSums(counts > 0) >= ${params.min_cells}]
    st <- subset(st, features = features_keep)
    
    cat("Cells after filtering:", ncol(st), "\\n")
    cat("Genes after filtering:", nrow(st), "\\n")
    
    # Save filtered QC plots with custom colors
    features <- c("nFeature_Spatial", "nCount_Spatial", "percent.mt")
    qc_colors <- c("#f38181", "#a8d8ea", "#fce38a")

    p7_list <- lapply(seq_along(features), function(i) {
    VlnPlot(st, features = features[i], pt.size = 0.1, fill = qc_colors[i]) +
        ggtitle(features[i]) +
        theme(legend.position = "none")
    })

    p7 <- p7_list[[1]] | p7_list[[2]] | p7_list[[3]]
    
    png("06_QC_filtered_violin.png", width = 15, height = 5, units = "in", res = 300)
    print(p7)
    dev.off()
    pdf("06_QC_filtered_violin.pdf", width = 15, height = 5)
    print(p7)
    dev.off()

    # Save Seurat object
    cat("Saving QC-processed Seurat object...\\n")
    saveRDS(st, "seurat_qc.rds")
    
    cat("QC analysis completed!\\n")
    """
}