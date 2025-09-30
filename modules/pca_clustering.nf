/*
 * PCA and Clustering Analysis Module
 */

process PCA_CLUSTERING {
    publishDir "${params.output_dir}/2.PCA", mode: 'copy'
    
    input:
    path seurat_rds
    
    output:
    path "seurat_clustered.rds", emit: seurat_clustered
    path "*.png"
    path "*.pdf"
    path "pca_loadings.csv"
    path "cluster_stats.csv"
    
    script:
    """
    #!/usr/bin/Rscript
    
    # Load required packages
    suppressMessages({
        library(Seurat)
        library(ggplot2)
        library(patchwork)
        library(dplyr)
    })
    
    cat("=== PCA Dimensionality Reduction and Clustering Analysis ===\\n")
    
    # Read QC-processed Seurat object
    cat("Reading QC-processed data...\\n")
    st <- readRDS("${seurat_rds}")
    
    cat("Number of input cells:", ncol(st), "\\n")
    cat("Number of input genes:", nrow(st), "\\n")
    
    # Data normalization
    cat("Data normalization...\\n")
    st <- SCTransform(st, assay = "Spatial", verbose = FALSE)
    
    # PCA dimensionality reduction
    cat("PCA dimensionality reduction analysis...\\n")
    st <- RunPCA(st, assay = "SCT", verbose = FALSE, npcs = ${params.clustering.n_pcs})
    
    # Save PCA loadings
    pca_loadings <- Loadings(st, reduction = "pca")
    write.csv(as.data.frame(pca_loadings), "pca_loadings.csv")
    
    # PCA visualization
    cat("Generating PCA plots...\\n")
    
    # 1. Elbow plot - show more PCs for better evaluation
    elbow_plot_dims <- if("${params.pca_elbow_plot_dims}" != "null") {
        as.numeric("${params.pca_elbow_plot_dims}")
    } else {
        50  # default to show 50 PCs
    }
    cat("Generating elbow plot for", elbow_plot_dims, "PCs\\n")
    
    p1 <- ElbowPlot(st, ndims = elbow_plot_dims) + 
        ggtitle("PCA Elbow Plot") +
        theme(text = element_text(family = "sans"))
    
    # 2. PCA scatter plot
    p2 <- DimPlot(st, reduction = "pca", dims = c(1, 2)) + 
        theme(text = element_text(family = "sans"))
    
    # 3. PCA heatmap
    p3 <- DimHeatmap(st, dims = 1:4, cells = 500, balanced = TRUE, ncol = 2)
    
    # Save PCA plots
    png("05_PCA_elbow_plot.png", width = 8, height = 6, units = "in", res = 300)
    print(p1)
    dev.off()
    
    png("06_PCA_scatter.png", width = 8, height = 6, units = "in", res = 300)
    print(p2)
    dev.off()
    
    png("07_PCA_heatmap.png", width = 12, height = 10, units = "in", res = 300)
    print(p3)
    dev.off()
    
    # UMAP dimensionality reduction
    cat("UMAP dimensionality reduction analysis...\\n")
    st <- RunUMAP(st, reduction = "pca", dims = 1:${params.clustering.chose_pcs}, 
                  n.neighbors = ${params.clustering.n_neighbors}, 
                  min.dist = ${params.clustering.min_dist},
                  verbose = FALSE)
    
    # Find neighbors and perform clustering
    cat("Computing neighbors and clustering...\\n")
    st <- FindNeighbors(st, reduction = "pca", dims = 1:${params.clustering.chose_pcs},
                        k.param = ${params.clustering.k_param}, verbose = FALSE)
    st <- FindClusters(st, resolution = ${params.clustering.resolution}, verbose = FALSE)
    
    # Clustering statistics
    cluster_stats <- data.frame(
        Cluster = names(table(st\$seurat_clusters)),
        Cell_Count = as.numeric(table(st\$seurat_clusters))
    )
    cluster_stats\$Percentage <- round(cluster_stats\$Cell_Count / sum(cluster_stats\$Cell_Count) * 100, 2)
    write.csv(cluster_stats, "cluster_stats.csv", row.names = FALSE)
    
    cat("Clustering results:\\n")
    print(cluster_stats)
    
    # UMAP clustering visualization
    cat("Generating UMAP clustering plots...\\n")
    
    # 1. UMAP clustering plot
    p4 <- DimPlot(st, reduction = "umap", group.by = "seurat_clusters", 
                  label = TRUE, label.size = 6) + 
        ggtitle(paste("UMAP Clustering (", length(unique(st\$seurat_clusters)), " clusters)")) +
        theme(text = element_text(family = "sans"))
    
    # 2. UMAP feature plots
    p5 <- FeaturePlot(st, features = c("nCount_Spatial", "nFeature_Spatial"), 
                      reduction = "umap", ncol = 2) & 
        theme(text = element_text(family = "sans"))
    
    # Save UMAP plots
    png("08_UMAP_clusters.png", width = 10, height = 8, units = "in", res = 300)
    print(p4)
    dev.off()
    
    png("09_UMAP_features.png", width = 12, height = 6, units = "in", res = 300)
    print(p5)
    dev.off()
    
    # Spatial clustering visualization
    cat("Generating spatial clustering plots...\\n")
    
    # 1. Spatial clustering plot
    p6 <- SpatialDimPlot(st, group.by = "seurat_clusters", label = TRUE, 
                         label.size = 3, pt.size.factor = 0.3) + 
        ggtitle("Cluster Spatial Distribution") +
        theme(legend.position = "right", text = element_text(family = "sans"))
    
    # 2. Spatial feature plots
    p7 <- SpatialFeaturePlot(st, features = "nCount_Spatial", pt.size.factor = 0.2) + 
        ggtitle("Spatial Distribution of Total Counts") +
        theme(legend.position = "right", text = element_text(family = "sans"))
    
    p8 <- SpatialFeaturePlot(st, features = "nFeature_Spatial", pt.size.factor = 0.2) + 
        ggtitle("Spatial Distribution of Gene Counts") +
        theme(legend.position = "right", text = element_text(family = "sans"))
    
    # Save spatial plots
    png("10_spatial_clusters.png", width = 10, height = 8, units = "in", res = 300)
    print(p6)
    dev.off()
    
    png("11_spatial_features.png", width = 16, height = 6, units = "in", res = 300)
    print(p7 | p8)
    dev.off()
    
    # Individual PDF files
    pdf("05_PCA_elbow_plot.pdf", width = 8, height = 6)
    print(p1)
    dev.off()
    
    pdf("06_PCA_scatter.pdf", width = 8, height = 6)
    print(p2)
    dev.off()
    
    pdf("07_PCA_heatmap.pdf", width = 12, height = 10)
    print(p3)
    dev.off()
    
    pdf("08_UMAP_clusters.pdf", width = 10, height = 8)
    print(p4)
    dev.off()
    
    pdf("09_UMAP_features.pdf", width = 12, height = 6)
    print(p5)
    dev.off()
    
    pdf("10_spatial_clusters.pdf", width = 10, height = 8)
    print(p6)
    dev.off()
    
    pdf("11_spatial_features.pdf", width = 16, height = 6)
    print(p7 | p8)
    dev.off()
    
    # Save clustered Seurat object
    cat("Saving clustered Seurat object...\\n")
    saveRDS(st, "seurat_clustered.rds")
    
    cat("PCA clustering analysis completed!\\n")
    cat("- Number of clusters:", length(unique(st\$seurat_clusters)), "\\n")
    cat("- Number of cells:", ncol(st), "\\n")
    """
}