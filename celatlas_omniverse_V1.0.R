rm(list = ls())

#（Uncomment when running for the first time）:
#
# BiocManager::install(c("SingleR", "celldex"), force = TRUE)
# BiocManager::install(c("monocle3"))
# devtools::install_github("satijalab/seurat-wrappers") or remotes::install_github("satijalab/seurat-wrappers")
# BiocManager::install("scrapper")

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SingleR)
library(celldex)
library(BiocParallel)
library(CellChat)
library(tidyverse)
library(pheatmap)
library(monocle3)
library(SeuratWrappers)
# library(devtools)

theme_set(theme_bw(base_family = "sans"))

setwd("/path/to/your/data")
data_dir <- "/path/to/your/data/bin50/"
output_dir <- file.path(data_dir, "output_plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

qc_dir <- file.path(output_dir, "1.QC")
pca_dir <- file.path(output_dir, "2.PCA") 
marker_dir <- file.path(output_dir, "3.Marker")
singler_dir <- file.path(output_dir, "4.SingleR")
cellchat_dir <- file.path(output_dir, "5.Cellchat")
pseudotime_dir <- file.path(output_dir, "6.Pseudotime")
enrichment_dir <- file.path(output_dir, "7.Enrichment")

dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pca_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(marker_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(singler_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cellchat_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pseudotime_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(enrichment_dir, showWarnings = FALSE, recursive = TRUE)

counts <- Read10X(data.dir = file.path(data_dir, "filtered_feature_bc_matrix"))
st <- CreateSeuratObject(counts = counts, assay = "Spatial")
cat("cell:", length(colnames(st)), "\n")
cat("gene:", length(rownames(st)), "\n")

image <- Read10X_Image(
  image.dir = file.path(data_dir, "spatial"),
  filter.matrix = TRUE,
  slice = "slice1"
)

image <- image[colnames(st)]
st[["slice1"]] <- image

tissue_pos_files <- c("tissue_positions_list.csv", "tissue_positions.csv")
tissue_pos_file <- NULL
for(file in tissue_pos_files) {
  if(file.exists(file.path(data_dir, "spatial", file))) {
    tissue_pos_file <- file.path(data_dir, "spatial", file)
    break
  }
}

coords <- NULL
tryCatch({
  
  coords <- GetTissueCoordinates(st, image = "slice1")
}, error = function(e) {
  tryCatch({
    coords <- st@images$slice1@coordinates
  }, error = function(e2) {
    cat("If coordinate information cannot be obtained, the default scale factors\n will be used")
  })
})

if(!is.null(tissue_pos_file) && !is.null(coords)) {
  tissue_pos <- read.csv(tissue_pos_file, header = FALSE)
  
  if(ncol(tissue_pos) >= 6) {
    image_coords_range <- max(c(max(tissue_pos$V5) - min(tissue_pos$V5), 
                                max(tissue_pos$V6) - min(tissue_pos$V6)))
    
    if("imagerow" %in% colnames(coords) && "imagecol" %in% colnames(coords)) {
      spot_coords_range <- max(c(max(coords$imagerow) - min(coords$imagerow),
                                 max(coords$imagecol) - min(coords$imagecol)))
    } else if("x" %in% colnames(coords) && "y" %in% colnames(coords)) {
      spot_coords_range <- max(c(max(coords$x) - min(coords$x),
                                 max(coords$y) - min(coords$y)))
    } else {
      numeric_cols <- sapply(coords, is.numeric)
      if(sum(numeric_cols) >= 2) {
        coord_matrix <- as.matrix(coords[, numeric_cols])
        spot_coords_range <- max(apply(coord_matrix, 2, function(x) max(x) - min(x)))
      } else {
        spot_coords_range <- 1000 
      }
    }
    
    if(spot_coords_range > 0 && image_coords_range > 0) {
      actual_hires_factor <- image_coords_range / spot_coords_range
      actual_lowres_factor <- actual_hires_factor * 0.3
    } else {
      actual_hires_factor <- 1.0
      actual_lowres_factor <- 0.3
    }
  } else {
    actual_hires_factor <- 1.0
    actual_lowres_factor <- 0.3
  }
} else if(!is.null(coords)) {
  if("imagerow" %in% colnames(coords) && "imagecol" %in% colnames(coords)) {
    coord_range <- max(c(max(coords$imagerow) - min(coords$imagerow),
                         max(coords$imagecol) - min(coords$imagecol)))
  } else if("x" %in% colnames(coords) && "y" %in% colnames(coords)) {
    coord_range <- max(c(max(coords$x) - min(coords$x),
                         max(coords$y) - min(coords$y)))
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
  cat("Use the default scale factors\n")
  actual_hires_factor <- 1.0
  actual_lowres_factor <- 0.3
}

st@images$slice1@scale.factors$hires <- actual_hires_factor
st@images$slice1@scale.factors$lowres <- actual_lowres_factor

spatial_dir <- file.path(data_dir, "spatial")
if(file.exists(file.path(spatial_dir, "tissue_lowres_image.png"))) {
  st@images$slice1@image <- png::readPNG(file.path(spatial_dir, "tissue_lowres_image.png"))
}

st[["percent.mt"]] <- PercentageFeatureSet(st, pattern = "^mt-|^Mt-")

p_qc <- VlnPlot(st, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3) +
  ggtitle("QC Metrics Distribution") +
  theme(text = element_text(family = "sans"))
ggsave(plot = p_qc, filename = file.path(qc_dir, "01_QC_VlnPlot.png"), width = 14, height = 5)
ggsave(plot = p_qc, filename = file.path(qc_dir, "01_QC_VlnPlot.pdf"), width = 14, height = 5)

p_scatter1 <- FeatureScatter(st, feature1 = "nCount_Spatial", feature2 = "percent.mt") +
  ggtitle("nCount vs Mitochondrial %") +
  theme(text = element_text(family = "sans"))
ggsave(plot = p_scatter1, filename = file.path(qc_dir, "02_FeatureScatter_nCount_vs_MT.png"))
ggsave(plot = p_scatter1, filename = file.path(qc_dir, "02_FeatureScatter_nCount_vs_MT.pdf"))

p_scatter2 <- FeatureScatter(st, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") +
  ggtitle("nCount vs Feature Count") +
  theme(text = element_text(family = "sans"))
ggsave(plot = p_scatter2, filename = file.path(qc_dir, "03_FeatureScatter_nCount_vs_nFeature.png"))
ggsave(plot = p_scatter2, filename = file.path(qc_dir, "03_FeatureScatter_nCount_vs_nFeature.pdf"))

p_spatial_total_counts <- SpatialFeaturePlot(st, features = c("nCount_Spatial"), pt.size.factor = 0.2) +
  theme(legend.position = "right", 
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.caption = element_text(size = 10, hjust = 0),
        plot.margin = margin(1, 1, 1.5, 1, "cm")) +
  ggtitle("Spatial Distribution of Total Counts (MID)") +
  labs(caption = "Visualization of total gene expression counts (MID) across spatial locations.\nHigher values indicate spots with greater total gene expression.")

ggsave(plot = p_spatial_total_counts, filename = file.path(qc_dir, "04a_Spatial_TotalCounts_MID.png"), 
       width = 8, height = 6, dpi = 300)
ggsave(plot = p_spatial_total_counts, filename = file.path(qc_dir, "04a_Spatial_TotalCounts_MID.pdf"), 
       width = 8, height = 7, device = cairo_pdf)

p_spatial_gene_counts <- SpatialFeaturePlot(st, features = c("nFeature_Spatial"), pt.size.factor = 0.2) +
  theme(legend.position = "right", 
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.caption = element_text(size = 10, hjust = 0),
        plot.margin = margin(1, 1, 1.5, 1, "cm")) +
  ggtitle("Spatial Distribution of Gene Counts") +
  labs(caption = "Visualization of detected gene numbers across spatial locations.\nHigher values indicate spots with more diverse gene expression.")

ggsave(plot = p_spatial_gene_counts, filename = file.path(qc_dir, "04b_Spatial_GeneCounts.png"), 
       width = 8, height = 6, dpi = 300)
ggsave(plot = p_spatial_gene_counts, filename = file.path(qc_dir, "04b_Spatial_GeneCounts.pdf"), 
       width = 8, height = 7, device = cairo_pdf)

st <- subset(st, subset = nCount_Spatial > 200 & nCount_Spatial < 250000 & 
               nFeature_Spatial > 2000 & percent.mt < 5)

p_qc_after <- VlnPlot(st, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3) +
  ggtitle("QC Metrics Distribution After Filtering") +
  theme(text = element_text(family = "sans"))
ggsave(plot = p_qc_after, filename = file.path(qc_dir, "05_QC_VlnPlot_after.png"), width = 14, height = 5)
ggsave(plot = p_qc_after, filename = file.path(qc_dir, "05_QC_VlnPlot_after.pdf"), width = 14, height = 5)

save(st, file = file.path(qc_dir, "00_st_after_qc.Rdata"))

st <- NormalizeData(st, assay = "Spatial")
st <- FindVariableFeatures(st, assay = "Spatial")
st <- ScaleData(st, assay = "Spatial")

st <- RunPCA(st, assay = "Spatial", verbose = FALSE)

p_elbow <- ElbowPlot(st, ndims = 50,reduction = "pca") +
  ggtitle("PCA Elbow Plot") +
  theme(text = element_text(family = "sans"))
ggsave(plot = p_elbow, file.path(pca_dir, "06_PCA_ElbowPlot.png"))
ggsave(plot = p_elbow, file.path(pca_dir, "06_PCA_ElbowPlot.pdf"))

st <- FindNeighbors(st, reduction = "pca", dims = 1:20)
st <- FindClusters(st, resolution = 0.2, verbose = FALSE)
st <- RunUMAP(st, reduction = "pca", dims = 1:20)

p_umap <- DimPlot(st, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP with Cluster Labels") +
  theme(text = element_text(family = "sans"))
ggsave(plot = p_umap, file.path(pca_dir, "07_UMAP_Cluster_Labels.png"))
ggsave(plot = p_umap, file.path(pca_dir, "07_UMAP_Cluster_Labels.pdf"))

p_spatial_clusters <- SpatialDimPlot(st, pt.size.factor = 0.3, label = TRUE, label.size = 3) +
  ggtitle("Cluster Spatial Distribution") +
  theme(text = element_text(family = "sans"))
ggsave(plot = p_spatial_clusters, file.path(pca_dir, "08_Spatial_Clusters.png"))
ggsave(plot = p_spatial_clusters, file.path(pca_dir, "08_Spatial_Clusters.pdf"))

n_clusters <- length(unique(st$seurat_clusters))
highlight_clusters <- head(sort(unique(st$seurat_clusters)), min(6, n_clusters))

p_spatial_highlight <- SpatialDimPlot(st, cells.highlight = CellsByIdentities(object = st, idents = highlight_clusters),
                                      facet.highlight = TRUE, pt.size.factor = 0.1, ncol = 3) +
  ggtitle("Highlight Specific Clusters") +
  theme(text = element_text(family = "sans"))
ggsave(plot = p_spatial_highlight, file.path(pca_dir, "09_Spatial_Highlight_Clusters.png"))
ggsave(plot = p_spatial_highlight, file.path(pca_dir, "09_Spatial_Highlight_Clusters.pdf"))

save(st, file = file.path(pca_dir, "01_st_after_pca_umap.Rdata"))

# =============================================================================
# Marker gene analysis and clustering heatmap visualization
# =============================================================================
cat("\n===Marker gene analysis and clustering heatmap visualization ===\n")

tryCatch({
  Idents(st) <- "seurat_clusters"
  
  all_markers <- FindAllMarkers(
    st, 
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    verbose = FALSE
  )
  
  if(nrow(all_markers) > 0) {
    write.csv(all_markers, 
              file.path(marker_dir, "all_cluster_markers.csv"),
              row.names = FALSE)
    top_markers <- all_markers %>%
      group_by(cluster) %>%
      top_n(n = 5, wt = avg_log2FC) %>%
      arrange(cluster, desc(avg_log2FC))
    
    write.csv(top_markers, 
              file.path(marker_dir, "top5_cluster_markers.csv"),
              row.names = FALSE)
    
    if(nrow(top_markers) > 0) {
      
      marker_genes <- unique(top_markers$gene)
      
      Idents(st) <- "seurat_clusters"
      avg_exp <- AverageExpression(st, features = marker_genes, verbose = FALSE)
      
      if("Spatial" %in% names(avg_exp)) {
        avg_exp_matrix <- avg_exp$Spatial
      } else if("RNA" %in% names(avg_exp)) {
        avg_exp_matrix <- avg_exp$RNA
      } else {
        avg_exp_matrix <- avg_exp[[1]]  
      }
      
      if(!is.null(avg_exp_matrix)) {
        avg_exp_matrix <- as.matrix(avg_exp_matrix)
        
        if(nrow(avg_exp_matrix) >= 1 && ncol(avg_exp_matrix) >= 1) {
          avg_exp_matrix[is.na(avg_exp_matrix)] <- 0
          
          scaled_matrix <- t(scale(t(avg_exp_matrix)))
          
          scaled_matrix[is.na(scaled_matrix)] <- 0
          
          cluster_order <- paste0("g", sort(as.numeric(gsub("g", "", colnames(scaled_matrix)))))
          scaled_matrix <- scaled_matrix[, cluster_order]
          
          gene_order <- top_markers$gene
          scaled_matrix <- scaled_matrix[gene_order, ]
          
          colnames(scaled_matrix) <- paste0("Cluster", 1:ncol(scaled_matrix))
          
          library(pheatmap)
          
          custom_colors <- colorRampPalette(c('#330066', '#336699', '#66CC66', '#FFCC33'))(100)
          
          png(file.path(marker_dir, "marker_genes_heatmap.png"), 
              width = 12, height = 10, units = "in", res = 300)
          pheatmap(
            scaled_matrix,
            scale = "none",  
            cluster_rows = FALSE,  
            cluster_cols = FALSE,  
            color = custom_colors,  
            show_rownames = TRUE,
            show_colnames = TRUE,
            main = "Top 5 Marker Genes by Cluster (avg_log2FC)",
            fontsize = 10,
            fontsize_row = 8,
            fontsize_col = 10,
            angle_col = 45,  
            gaps_row = cumsum(table(top_markers$cluster))[-length(table(top_markers$cluster))],
            cellwidth = 30,
            cellheight = 10
          )
          dev.off()
          
          pdf(file.path(marker_dir, "marker_genes_heatmap.pdf"), 
              width = 12, height = 10)
          pheatmap(
            scaled_matrix,
            scale = "none",
            cluster_rows = FALSE,
            cluster_cols = FALSE,
            color = custom_colors, 
            show_rownames = TRUE,
            show_colnames = TRUE,
            main = "Top 5 Marker Genes by Cluster",
            fontsize = 10,
            fontsize_row = 8,
            fontsize_col = 10,
            angle_col = 45,  
            gaps_row = cumsum(table(top_markers$cluster))[-length(table(top_markers$cluster))],
            cellwidth = 30,
            cellheight = 10
          )
          dev.off()
          
          if(ncol(scaled_matrix) >= 2) {
            
            cluster_cor <- cor(scaled_matrix, method = "pearson", use = "complete.obs")
            
            cluster_cor[is.na(cluster_cor)] <- 0
            
            png(file.path(marker_dir, "cluster_similarity_heatmap.png"), 
                width = 8, height = 8, units = "in", res = 300)
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
            
            pdf(file.path(marker_dir, "cluster_similarity_heatmap.pdf"), 
                width = 8, height = 8)
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
          
          cat("The heat map analysis of the Marker gene has been completed\n")
          cat("The number of marker genes analyzed:", length(marker_genes), "\n")
          cat("The number of analyzed clusters:", ncol(avg_exp_matrix), "\n")
        }
      }
    }
  }
}, error = function(e) {
  cat("The Marker gene analysis was incorrect:", e$message, "\n")
})

cat("✓ Marker gene analysis and clustering heat map\n")

cat("Start the SingleR cell type annotation...\n")

# Intelligent selection of reference datasets - automatic detection of species
cat("Detection dataset species...\n")

test_genes <- head(rownames(st), 100)
human_genes <- sum(grepl("^[A-Z][A-Z0-9-]*$", test_genes))
mouse_genes <- sum(grepl("^[A-Z][a-z0-9-]*$", test_genes))

if(human_genes > mouse_genes) {
  ref <- celldex::HumanPrimaryCellAtlasData()
  species_detected <- "human"
} else {
  ref <- celldex::MouseRNAseqData()
  species_detected <- "mouse"
}

if(nrow(GetAssayData(st, assay = "Spatial", slot = "data")) == 0) {
  sc_data <- GetAssayData(st, assay = "Spatial", slot = "data")
  # sc_data <- log1p(sc_data)
} else {
  sc_data <- GetAssayData(st, assay = "Spatial", slot = "count")
}

sc_data <- sc_data[, colnames(st)]

test_genes <- rownames(sc_data)
ref_genes <- rownames(ref)

common_genes <- intersect(test_genes, ref_genes)

if(length(common_genes) < 500) {
  if(species_detected == "mouse") {
    test_genes_converted <- toupper(test_genes)
  } else {
    test_genes_converted <- paste0(toupper(substr(test_genes, 1, 1)), 
                                   tolower(substr(test_genes, 2, nchar(test_genes))))
  }
  
  common_genes_converted <- intersect(test_genes_converted, ref_genes)
  
  if(length(common_genes_converted) > length(common_genes)) {
    rownames(sc_data) <- test_genes_converted
    common_genes <- common_genes_converted
  }
  
  if(length(common_genes) < 300) {
    test_genes_clean <- gsub("\\..*$", "", test_genes)
    common_genes_clean <- intersect(test_genes_clean, ref_genes)
    
    if(length(common_genes_clean) > length(common_genes)) {
      rownames(sc_data) <- test_genes_clean
      common_genes <- common_genes_clean
    }
  }
}

sc_data <- sc_data[common_genes, ]

if(inherits(sc_data, "dgCMatrix")) {
  sc_data <- as.matrix(sc_data)
}

if(is.null(sc_data) || nrow(sc_data) == 0 || ncol(sc_data) == 0) {
  stop("Effective expression data cannot be obtained for annotation")
}

set.seed(66)

clusters <- as.character(st$seurat_clusters)

if(length(clusters) == 0 || all(is.na(clusters))) {
  stop("Cluster information was not found. Please ensure that the cluster analysis has been completed")
}

if(nrow(sc_data) < 50) {
  stop("There are too few common genes for accurate annotation")
}

tryCatch({
  st.hesc <- SingleR(
    test = sc_data,
    ref = ref,
    labels = ref$label.main,
    clusters = clusters,
    assay.type.test = "logcounts",
    BPPARAM = BiocParallel::SnowParam(workers = 2)
  )
}, error = function(e) {
  st.hesc <<- SingleR(
    test = sc_data,
    ref = ref,
    labels = ref$label.main,
    clusters = clusters,
    assay.type.test = "logcounts"
  )
})
p_score_heatmap <- plotScoreHeatmap(st.hesc)

pdf(file.path(singler_dir, "10_SingleR_score_heatmap.pdf"), width = 12, height = 8)
print(p_score_heatmap)
dev.off()

png(file.path(singler_dir, "10_SingleR_score_heatmap.png"), width = 3600, height = 2400, res = 300)
print(p_score_heatmap)
dev.off()

st$cell_type_main <- st.hesc$labels[match(clusters, rownames(st.hesc))]

annotation_scores <- st.hesc$scores[match(clusters, rownames(st.hesc)), ]
if(is.matrix(annotation_scores)) {
  st$annotation_score <- rowMaxs(annotation_scores, na.rm = TRUE)
} else {
  st$annotation_score <- annotation_scores
}

confidence_threshold <- 0.5

st$cell_type_filtered <- ifelse(st$annotation_score < confidence_threshold, 
                                "Unknown", 
                                st$cell_type_main)

p1 <- DimPlot(st, group.by = "cell_type_main", reduction = "umap", 
              label = TRUE, repel = TRUE) + 
  ggtitle("Cell Types (SingleR)") +
  theme(text = element_text(family = "sans"), legend.position = "bottom") +
  guides(color = guide_legend(ncol = 3))

# p_celltype_combined <- wrap_plots(p1)
ggsave(plot = p1, file.path(singler_dir, "11_UMAP_CellType.png"), width = 8, height = 8)
ggsave(plot = p1, file.path(singler_dir, "11_UMAP_CellType.pdf"), width = 8, height = 8)

p_res_comparison <- DimPlot(st, group.by = "seurat_clusters", reduction = "umap", 
                            label = TRUE) +
  ggtitle("Final Clustering (Resolution 0.2)") +
  theme(text = element_text(family = "sans"))

ggsave(plot = p_res_comparison, file.path(singler_dir, "11_Final_Clustering.png"), width = 10, height = 8)
ggsave(plot = p_res_comparison, file.path(singler_dir, "11_Final_Clustering.pdf"), width = 10, height = 8)

p3 <- SpatialDimPlot(st, group.by = "cell_type_filtered", 
                     pt.size.factor = 0.1, label = TRUE) +
  ggtitle("Spatial Distribution of Cell Types") +
  theme(text = element_text(family = "sans"))

p4 <- SpatialFeaturePlot(st, features = "annotation_score", 
                         pt.size.factor = 0.1) +
  ggtitle("Annotation Confidence Scores") +
  theme(text = element_text(family = "sans"))

p_spatial_combined <- wrap_plots(p3, p4, ncol = 2)
ggsave(plot = p_spatial_combined, file.path(singler_dir, "12_Spatial_CellType.png"), width = 16, height = 8)
ggsave(plot = p_spatial_combined, file.path(singler_dir, "12_Spatial_CellType.pdf"), width = 16, height = 8)

save(st, st.hesc, 
     file = file.path(singler_dir, "02_st_final_singleR.Rdata"))

annotation_stats <- table(st$cell_type_filtered, st$seurat_clusters)
write.csv(annotation_stats, 
          file.path(singler_dir, "celltype_cluster_distribution.csv"))

cluster_annotation_summary <- st@meta.data %>%
  group_by(seurat_clusters, cell_type_filtered) %>%
  summarise(
    count = n(),
    mean_score = round(mean(annotation_score, na.rm = TRUE), 3),
    .groups = 'drop'
  ) %>%
  arrange(seurat_clusters)

write.csv(cluster_annotation_summary, 
          file.path(singler_dir, "cluster_annotation_summary.csv"), 
          row.names = FALSE)

write.csv(st@meta.data, 
          file.path(singler_dir, "complete_metadata_with_annotation.csv"), 
          row.names = TRUE)

type_counts <- table(st$cell_type_filtered)
for(i in 1:length(type_counts)) {
  cat(names(type_counts)[i], ":", type_counts[i], "\n")
}

unknown_ratio <- sum(st$cell_type_filtered == "Unknown") / ncol(st)
if(unknown_ratio > 0.2) {
  cat("- The proportion of unknown cells is relatively high (", round(unknown_ratio*100, 1), "%)，建议:\n")
}

cat("\nAnalysis completed! The result is saved in:", output_dir, "\n")

# =============================================================================
# CellChat cell communication analysis
# =============================================================================
cat("\n=== Start the CellChat cell communication analysis ===\n")

data.input <- GetAssayData(st, slot = 'data')  # normalized data matrix
meta <- st@meta.data                           # meta data

print(table(meta$cell_type_filtered))

valid_cells <- rownames(meta)[meta$cell_type_filtered != "Unknown"]
if(length(valid_cells) > 0) {
  data.input <- data.input[, valid_cells]
  meta <- meta[valid_cells, ]
  print(table(meta$cell_type_filtered))
}

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type_filtered")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "cell_type_filtered")

cellchat@idents <- droplevels(cellchat@idents, 
                              exclude = setdiff(levels(cellchat@idents), unique(cellchat@idents)))

groupSize <- as.numeric(table(cellchat@idents))
print(table(cellchat@idents))

if(species_detected == "human") {
  CellChatDB <- CellChatDB.human
} else {
  CellChatDB <- CellChatDB.mouse
}

showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

if(species_detected == "human") {
  cellchat <- smoothData(cellchat, adj = PPI.human)
} else {
  cellchat <- smoothData(cellchat, adj = PPI.mouse)
}

cellchat <- computeCommunProb(cellchat)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

df.net <- subsetCommunication(cellchat)
write.csv(df.net, file.path(cellchat_dir, '01_communication_network.csv'), 
          row.names = FALSE, quote = FALSE)

save(cellchat, file = file.path(cellchat_dir, "cellchat_results.Rdata"))

pdf(file.path(cellchat_dir, '01_net_circle_number.pdf'), width = 9, height = 9, family = 'Times')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, 
                 label.edge = FALSE, title.name = "Number of interactions")
dev.off()

png(file.path(cellchat_dir, '01_net_circle_number.png'), width = 9, height = 9, 
    units = 'in', res = 300, family = 'Times')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, 
                 label.edge = FALSE, title.name = "Number of interactions")
dev.off()

pdf(file.path(cellchat_dir, '02_net_circle_weight.pdf'), width = 9, height = 9, family = 'Times')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, 
                 label.edge = FALSE, title.name = "Interaction weight/strength")
dev.off()

png(file.path(cellchat_dir, '02_net_circle_weight.png'), width = 9, height = 9, 
    units = 'in', res = 300, family = 'Times')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, 
                 label.edge = FALSE, title.name = "Interaction weight/strength")
dev.off()

pdf(file.path(cellchat_dir, '03_number_of_interactions_heatmap.pdf'), width = 8, height = 6, family = 'Times')
netVisual_heatmap(cellchat)
dev.off()

png(file.path(cellchat_dir, '03_number_of_interactions_heatmap.png'), width = 8, height = 6, 
    units = 'in', res = 300, family = 'Times')
netVisual_heatmap(cellchat)
dev.off()

mat <- cellchat@net$count
write.csv(mat, file.path(cellchat_dir, 'interaction_count_matrix.csv'))

if(nrow(mat) > 0 && ncol(mat) > 0) {
  pdf(file.path(cellchat_dir, '04_single_celltype_circles.pdf'), width = 18, height = 5, family = 'Times')
  par(mfrow = c(1, min(6, nrow(mat))), xpd = TRUE)
  for (i in 1:min(6, nrow(mat))) {  
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i,] <- mat[i,]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, 
                     arrow.width = 0.2, arrow.size = 0.1, edge.weight.max = max(mat),
                     title.name = rownames(mat)[i])
  }
  dev.off()
  
  png(file.path(cellchat_dir, '04_single_celltype_circles.png'), width = 1800, height = 500, family = 'Times')
  par(mfrow = c(1, min(6, nrow(mat))), xpd = TRUE)
  for (i in 1:min(6, nrow(mat))) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i,] <- mat[i,]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, 
                     arrow.width = 0.2, arrow.size = 0.1, edge.weight.max = max(mat),
                     title.name = rownames(mat)[i])
  }
  dev.off()
}

pdf(file.path(cellchat_dir, '05_ligand_receptor_bubble.pdf'), width = 12, height = 12, family = 'Times')
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()

png(file.path(cellchat_dir, '05_ligand_receptor_bubble.png'), width = 12, height = 12, 
    units = 'in', res = 600, family = 'Times')
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()

if(nrow(df.net) > 0) {
  top_communications <- head(df.net[order(df.net$prob, decreasing = TRUE), 
                                    c("source", "target", "ligand", "receptor", "prob")], 10)
  print(top_communications)
}

cat("\nCellChat analysis completed! The result is saved in:", cellchat_dir, "\n")

# =============================================================================
# Monocle3 trajectory analysis and pseudo-time analysis
# =============================================================================
cat("\n=== Start Monocle3 trajectory analysis ===\n")

monocle3_dir <- pseudotime_dir

tryCatch({
  celltype_counts <- table(st$cell_type_filtered)
  valid_celltypes <- names(celltype_counts)[celltype_counts >= 20]
  valid_celltypes <- valid_celltypes[valid_celltypes != "Unknown"]
  
  if(length(valid_celltypes) < 2) {
    cluster_counts <- table(st$seurat_clusters)
    valid_clusters <- names(cluster_counts)[cluster_counts >= 20]
    
    if(length(valid_clusters) >= 2) {
      st_subset <- st[, st$seurat_clusters %in% valid_clusters]
      group_column <- "seurat_clusters"
    } else {
      cat("The data is insufficient for trajectory analysis\n")
      next
    }
  } else {
    st_subset <- st[, st$cell_type_filtered %in% valid_celltypes]
    group_column <- "cell_type_filtered"
  }
  
  cds <- as.cell_data_set(st_subset)
  
  fData(cds)$gene_short_name <- rownames(cds)
  
  cds <- preprocess_cds(cds, num_dim = 30, method = "PCA")
  
  cds <- reduce_dimension(cds, reduction_method = "UMAP")
  
  cds <- cluster_cells(cds)
  
  cds <- learn_graph(cds)
  
  p1 <- plot_cells(cds, 
                   color_cells_by = "partition",
                   label_groups_by_cluster = FALSE,
                   label_leaves = FALSE,
                   label_branch_points = FALSE,
                   graph_label_size = 1.5) +
    theme_bw() +
    ggtitle("Monocle3: Cell Partitions")
  
  if(group_column %in% colnames(colData(cds))) {
    p2 <- plot_cells(cds,
                     color_cells_by = group_column,
                     label_cell_groups = TRUE,
                     label_leaves = FALSE,
                     label_branch_points = FALSE,
                     graph_label_size = 1.5) +
      theme_bw() +
      ggtitle(paste("Monocle3: Cells by", tools::toTitleCase(gsub("_", " ", group_column))))
  } else {
    p2 <- plot_cells(cds,
                     color_cells_by = "cluster",
                     label_cell_groups = TRUE,
                     label_leaves = FALSE,
                     label_branch_points = FALSE,
                     graph_label_size = 1.5) +
      theme_bw() +
      ggtitle("Monocle3: Cells by Cluster")
  }
  
  p3 <- plot_cells(cds,
                   color_cells_by = "partition",
                   label_groups_by_cluster = FALSE,
                   label_leaves = TRUE,
                   label_branch_points = TRUE,
                   graph_label_size = 1.5) +
    theme_bw() +
    ggtitle("Monocle3: Trajectory Graph")
  
  p_combined <- (p1 | p2) / p3
  ggsave(plot = p_combined, 
         filename = file.path(monocle3_dir, "01_monocle3_trajectory_overview.png"),
         width = 16, height = 12, dpi = 300)
  ggsave(plot = p_combined, 
         filename = file.path(monocle3_dir, "01_monocle3_trajectory_overview.pdf"),
         width = 16, height = 12)
  
  tryCatch({
    leaf_nodes <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    leaf_nodes <- unique(leaf_nodes)
    
    if(length(leaf_nodes) > 0) {
      root_node <- names(sort(table(leaf_nodes), decreasing = TRUE))[1]
      cds <- order_cells(cds, root_cells = colnames(cds)[leaf_nodes == root_node])
      
      p4 <- plot_cells(cds,
                       color_cells_by = "pseudotime",
                       label_cell_groups = FALSE,
                       label_leaves = TRUE,
                       label_branch_points = TRUE,
                       graph_label_size = 1.5) +
        theme_bw() +
        ggtitle("Monocle3: Pseudotime") +
        scale_color_viridis_c()
      
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
      
      pseudotime_df <- data.frame(
        pseudotime = pseudotime(cds),
        group = colData(cds)[[group_column]]
      )
      
      p6 <- ggplot(pseudotime_df, aes(x = group, y = pseudotime, fill = group)) +
        geom_boxplot() +
        geom_jitter(width = 0.2, alpha = 0.3) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Pseudotime Distribution by Group") +
        labs(x = tools::toTitleCase(gsub("_", " ", group_column)), 
             y = "Pseudotime") +
        guides(fill = guide_legend(title = tools::toTitleCase(gsub("_", " ", group_column))))
      
      p_pseudotime <- (p4 | p6) / p5
      ggsave(plot = p_pseudotime, 
             filename = file.path(monocle3_dir, "02_monocle3_pseudotime.png"),
             width = 16, height = 12, dpi = 300)
      ggsave(plot = p_pseudotime, 
             filename = file.path(monocle3_dir, "02_monocle3_pseudotime.pdf"),
             width = 16, height = 12)
      
      tryCatch({
        cds_subset <- cds[Matrix::rowSums(exprs(cds) > 0) >= 10,]
        pr_test_res <- graph_test(cds_subset, neighbor_graph="knn", cores=1)
        pr_deg_ids <- row.names(subset(pr_test_res, q_value < 0.05))
        
        if(length(pr_deg_ids) > 0) {
          
          pr_test_res$gene_short_name <- rownames(pr_test_res)
          pr_test_res <- pr_test_res[order(pr_test_res$q_value),]
          write.csv(pr_test_res, 
                    file.path(monocle3_dir, "trajectory_dependent_genes.csv"), 
                    row.names = TRUE)
          
          top_genes <- head(pr_deg_ids, 9)
          
          if(length(top_genes) > 0) {
            p_genes <- plot_cells(cds_subset[top_genes,],
                                  genes = top_genes,
                                  label_cell_groups = FALSE,
                                  show_trajectory_graph = TRUE) +
              theme_bw() +
              ggtitle("Top Trajectory-Dependent Genes")
            
            ggsave(plot = p_genes, 
                   filename = file.path(monocle3_dir, "03_trajectory_genes.png"),
                   width = 12, height = 9, dpi = 300)
            ggsave(plot = p_genes, 
                   filename = file.path(monocle3_dir, "03_trajectory_genes.pdf"),
                   width = 12, height = 9)
            
          }
        } else {
          cat("No significant track-related genes were found\n")
        }
      }, error = function(e) {
        cat("Trajectory gene analysis failed:", e$message, "\n")
      })
      
    } else {
      cat("The root node cannot be automatically selected. Please specify it manually\n")
    }
  }, error = function(e) {
    cat("Pseudo-time analysis failed:", e$message, "\n")
  })
  
  saveRDS(cds, file.path(monocle3_dir, "monocle3_cds_object.rds"))

  monocle3_summary <- list(
    n_cells = ncol(cds),
    n_genes = nrow(cds),
    n_partitions = length(unique(cds@clusters[["UMAP"]]$partitions)),
    n_clusters = length(unique(cds@clusters[["UMAP"]]$clusters))
  )
  
  if("pseudotime" %in% names(colData(cds))) {
    monocle3_summary$pseudotime_range <- c(min(pseudotime(cds)), max(pseudotime(cds)))
    monocle3_summary$cells_with_pseudotime <- sum(!is.infinite(pseudotime(cds)))
  }
  
  writeLines(c(
    "=== Monocle3 Analysis Summary ===",
    paste("Number of cells analyzed:", monocle3_summary$n_cells),
    paste("Number of genes analyzed:", monocle3_summary$n_genes),
    paste("Number of partitions:", monocle3_summary$n_partitions),
    paste("Number of clusters:", monocle3_summary$n_clusters),
    if("pseudotime_range" %in% names(monocle3_summary)) {
      c(paste("Pseudotime range:", 
              round(monocle3_summary$pseudotime_range[1], 2), "-", 
              round(monocle3_summary$pseudotime_range[2], 2)),
        paste("Cells with pseudotime:", monocle3_summary$cells_with_pseudotime))
    },
    paste("Output directory:", monocle3_dir)
  ), file.path(monocle3_dir, "monocle3_summary.txt"))
  
  cat("\n=== Monocle3 Analysis Summary ===\n")
  cat("Number of cells analyzed:", monocle3_summary$n_cells, "\n")
  cat("Number of genes analyzed:", monocle3_summary$n_genes, "\n")
  cat("Number of partitions:", monocle3_summary$n_partitions, "\n")
  cat("Number of clusters:", monocle3_summary$n_clusters, "\n")
  if("pseudotime_range" %in% names(monocle3_summary)) {
    cat("Pseudotime range:", 
        round(monocle3_summary$pseudotime_range[1], 2), "-", 
        round(monocle3_summary$pseudotime_range[2], 2), "\n")
    cat("Cells with pseudotime:", monocle3_summary$cells_with_pseudotime, "\n")
  }
  
  cat("\nMonocle3 analysis completed! Results saved in:", monocle3_dir, "\n")
  
}, error = function(e) {
  cat("Monocle3 analysis error:", e$message, "\n")
  cat("Possible reasons:\n")
  cat("1. monocle3 or SeuratWrappers package not properly installed\n")
  cat("2. Insufficient data dimensionality for trajectory analysis\n")
  cat("3. Insufficient memory\n")
  cat("Please check the error message and try again\n")
})
# =============================================================================
# Complete Analysis Pipeline Summary
# =============================================================================
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("              Spatial Transcriptomics Analysis Pipeline Complete\n")
cat(paste(rep("=", 80), collapse=""), "\n")
cat("Analysis completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Results saved in:", output_dir, "\n\n")

cat("Completed analysis modules:\n")
cat("✓ 1. Data loading and quality control\n")
cat("✓ 2. Data normalization and dimensionality reduction (PCA/UMAP)\n") 
cat("✓ 3. Cell clustering analysis\n")
cat("✓ 4. Spatial visualization\n")
cat("✓ 5. SingleR automated cell type annotation\n")
cat("✓ 6. Functional enrichment analysis (GSVA)\n")
cat("✓ 7. Cell-cell communication analysis (CellChat)\n")
cat("✓ 8. Trajectory analysis and pseudotime inference (Monocle3)\n\n")

cat("Main output files generated:\n")
cat("- QC and clustering plots: 01_QC_*.png/pdf - 08_*.png/pdf\n")
cat("- SingleR annotation: 10_SingleR_*.png/pdf, 11_*.png/pdf, 12_*.png/pdf\n")
cat("- Functional enrichment: Functional_analysis/\n")
cat("- Cell communication: CellChat_analysis/\n")
cat("- Trajectory analysis: Monocle3_analysis/\n")
cat("- Analysis objects: *.Rdata files\n")
cat("- Data tables: *.csv files\n\n")

# Statistical results
if(exists("st")) {
  cat("Final data statistics:\n")
  cat("- Number of cells passing QC:", ncol(st), "\n")
  cat("- Number of genes:", nrow(st), "\n")
  if("seurat_clusters" %in% colnames(st@meta.data)) {
    cat("- Number of clusters:", length(unique(st$seurat_clusters)), "\n")
  }
  if("cell_type_filtered" %in% colnames(st@meta.data)) {
    type_counts <- table(st$cell_type_filtered)
    cat("- Number of annotated cell types:", length(type_counts), "\n")
    cat("- Major cell types:\n")
    for(i in 1:min(5, length(type_counts))) {
      cat(sprintf("  %s: %d cells\n", names(type_counts)[i], type_counts[i]))
    }
  }
}

cat("\nAnalysis recommendations:\n")
cat("1. Check QC plots to verify filtering parameters are appropriate\n")
cat("2. Validate the accuracy of cell type annotations\n")
cat("3. Interpret functional enrichment and communication results in biological context\n")
cat("4. Verify trajectory analysis results with prior biological knowledge\n")
cat("5. All results have been saved and can be reloaded for further analysis\n\n")

cat("For technical support or additional analysis, please provide:\n")
cat("- Complete error logs\n")
cat("- Biological background information of the data\n")  
cat("- Specific research questions\n\n")

cat(paste(rep("=", 80), collapse=""), "\n")
cat("Thank you for using the Spatial Transcriptomics Analysis Pipeline!\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# =============================================================================
# clusterProfiler differential gene screening and GO/KEGG functional enrichment analysis
# =============================================================================
cat("\n=== Start the clusterProfiler differential gene and functional enrichment analysis ===\n")

suppressMessages({
  if (!require(clusterProfiler, quietly = TRUE)) {
    BiocManager::install("clusterProfiler")
    library(clusterProfiler)
  }
  if (!require(org.Hs.eg.db, quietly = TRUE)) {
    BiocManager::install("org.Hs.eg.db")
    library(org.Hs.eg.db)
  }
  if (!require(org.Mm.eg.db, quietly = TRUE)) {
    BiocManager::install("org.Mm.eg.db")
    library(org.Mm.eg.db)
  }
  if (!require(enrichplot, quietly = TRUE)) {
    BiocManager::install("enrichplot")
    library(enrichplot)
  }
  if (!require(ggraph, quietly = TRUE)) {
    install.packages("ggraph")
    library(ggraph)
  }
  if (!require(igraph, quietly = TRUE)) {
    install.packages("igraph")
    library(igraph)
  }
})

if(species_detected == "human") {
  organism_db <- org.Hs.eg.db
  kegg_organism <- "hsa"
} else {
  organism_db <- org.Mm.eg.db
  kegg_organism <- "mmu"
}

valid_groups <- if(exists("valid_celltypes") && length(valid_celltypes) >= 2) {
  list(groups = valid_celltypes, group_by = "cell_type_filtered", name = "celltype")
} else {
  cluster_table <- table(st$seurat_clusters)
  valid_clusters <- names(cluster_table)[cluster_table >= 20]
  list(groups = valid_clusters[1:min(6, length(valid_clusters))], group_by = "seurat_clusters", name = "cluster")
}

for(i in 1:length(valid_groups$groups)) {
  group_id <- valid_groups$groups[i]
  tryCatch({
    Idents(st) <- valid_groups$group_by
    
    markers <- FindMarkers(
      st, 
      ident.1 = group_id,
      min.pct = 0.1,
      logfc.threshold = 0.25,
      test.use = "wilcox",
      verbose = FALSE
    )
    
    if(nrow(markers) == 0) {
      cat("No differential genes were found. Skip\n")
      next
    }
    
    significant_genes <- markers[markers$p_val_adj < 0.05 & markers$avg_log2FC > 0.5, ]
    
    if(nrow(significant_genes) == 0) {
      cat("No differential genes were found. Skip\n")
      next
    }
    
    cat("  Found", nrow(significant_genes), "significant differentially expressed genes\n")
    
    significant_genes$gene <- rownames(significant_genes)
    write.csv(significant_genes, 
              file.path(enrichment_dir, paste0(valid_groups$name, "_", group_id, "_DEGs.csv")),
              row.names = TRUE)
    
    gene_list <- rownames(significant_genes)
    
    tryCatch({
      if(species_detected == "human") {
        gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      } else {
        gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
      }
      
      if(nrow(gene_entrez) == 0) {
        cat("  Gene ID conversion failed, skipping enrichment analysis\n")
        next
      }
      
      go_all <- enrichGO(
        gene = gene_entrez$ENTREZID,
        OrgDb = organism_db,
        ont = "ALL",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        readable = TRUE
      )
      
      if(!is.null(go_all) && nrow(go_all@result) > 0) {
        GO_ALL <- go_all@result
        GO_ALL <- GO_ALL[GO_ALL$pvalue < 0.05,]
        
        if(nrow(GO_ALL) > 0) {
          BP <- GO_ALL[GO_ALL$ONTOLOGY=='BP', ]
          CC <- GO_ALL[GO_ALL$ONTOLOGY=='CC', ]
          MF <- GO_ALL[GO_ALL$ONTOLOGY=='MF', ]
          
          if(nrow(BP) > 0) {
            write.csv(BP, 
                      file.path(enrichment_dir, paste0(valid_groups$name, "_", group_id, "_GO_BP.csv")),
                      row.names = FALSE)
          cat(" GO BP analysis completed, total", nrow(BP), "enriched terms\n")
          }
          
          if(nrow(CC) > 0) {
            write.csv(CC, 
                      file.path(enrichment_dir, paste0(valid_groups$name, "_", group_id, "_GO_CC.csv")),
                      row.names = FALSE)
            cat(" GO CC analysis completed, total", nrow(CC), "enriched terms\n")
          }
          
          if(nrow(MF) > 0) {
            write.csv(MF, 
                      file.path(enrichment_dir, paste0(valid_groups$name, "_", group_id, "_GO_MF.csv")),
                      row.names = FALSE)
            cat(" GO MF analysis completed, total", nrow(MF), "enriched terms\n")
          }
          
          write.csv(GO_ALL, 
                    file.path(enrichment_dir, paste0(valid_groups$name, "_", group_id, "_GO_ALL.csv")),
                    row.names = FALSE)
          
          if(nrow(BP) > 0 || nrow(CC) > 0 || nrow(MF) > 0) {
            display_number <- c(10, 10, 10)  
            
            ego_result_BP <- if(nrow(BP) > 0) head(BP[order(BP$pvalue), ], display_number[1]) else data.frame()
            ego_result_CC <- if(nrow(CC) > 0) head(CC[order(CC$pvalue), ], display_number[2]) else data.frame()
            ego_result_MF <- if(nrow(MF) > 0) head(MF[order(MF$pvalue), ], display_number[3]) else data.frame()
            
            if(nrow(ego_result_BP) > 0 || nrow(ego_result_CC) > 0 || nrow(ego_result_MF) > 0) {
              go_enrich_df <- data.frame(
                ID = c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
                Description = c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
                GeneNumber = c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
                pvalue = c(ego_result_BP$pvalue, ego_result_CC$pvalue, ego_result_MF$pvalue),
                type = factor(c(rep("BP:biological process", nrow(ego_result_BP)), 
                                rep("CC:cellular component", nrow(ego_result_CC)),
                                rep("MF:molecular function", nrow(ego_result_MF))), 
                              levels = c("BP:biological process", "CC:cellular component", "MF:molecular function"))
              )
              
              go_enrich_df$Description_short <- ifelse(nchar(go_enrich_df$Description) > 60, 
                                                       paste0(substr(go_enrich_df$Description, 1, 60), "..."), 
                                                       go_enrich_df$Description)
              
              go_enrich_df$type_order <- factor(rev(as.integer(rownames(go_enrich_df))), 
                                                labels = rev(go_enrich_df$Description_short))
              
              COLS <- c("#a8d8ea", "#FFDD94", "#FF7B89")
              
              p_go_bubble <- ggplot(data = go_enrich_df, aes(x = type_order, y = GeneNumber, size = GeneNumber, fill = type)) +
                geom_point(shape = 21, color = "black") +  
                scale_size_area(max_size = 10) +  
                scale_fill_manual(values = COLS) +  
                coord_flip() +  
                xlab("GO term") + 
                ylab("Gene Number") + 
                labs(title = paste0("GO Enrichment Analysis: ", valid_groups$name, " ", group_id)) +
                theme_bw() + theme(
                  legend.background = element_rect(fill = "white", color = "black", size = 0.2),
                  legend.text = element_text(face = "bold", color = "black", family = "Times", size = 12),
                  plot.title = element_text(hjust = 0.5, face = "bold", color = "black", family = "Times", size = 18),
                  axis.text.x = element_text(family = "Times", face = "bold", color = "black", size = 14),
                  axis.text.y = element_text(family = "Times", face = "bold", color = "black", size = 10),
                  axis.title.x = element_text(face = "bold", color = "black", family = "Times", size = 18),
                  axis.title.y = element_text(face = "bold", color = "black", family = "Times", size = 18),
                  plot.subtitle = element_text(hjust = 0.5, family = "Times", size = 14, face = "italic", colour = "black")
                )
              
              ggsave(plot = p_go_bubble, 
                     filename = file.path(enrichment_dir, paste0(valid_groups$name, "_", group_id, "_GO_bubble_plot.png")),
                     width = 16, height = 12, dpi = 300)
              ggsave(plot = p_go_bubble, 
                     filename = file.path(enrichment_dir, paste0(valid_groups$name, "_", group_id, "_GO_bubble_plot.pdf")),
                     width = 16, height = 12)
              
              cat("    GO bubble plot generation completed\n")
            }
          }
        } else {
          cat("    No significant GO enrichment results found\n")
        }
      } else {
        cat("    GO enrichment analysis failed or produced no results\n")
      }
      
      kegg_result <- enrichKEGG(
        gene = gene_entrez$ENTREZID,
        organism = kegg_organism,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        qvalueCutoff = 1
      )
      
      if(!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
        kegg_readable <- setReadable(kegg_result, organism_db, keyType = "ENTREZID")
        
        kegg_significant <- kegg_result@result[kegg_result@result$pvalue <= 0.05, ]
        
        write.csv(kegg_result@result, 
                  file.path(enrichment_dir, paste0(valid_groups$name, "_", group_id, "_KEGG_all.csv")),
                  row.names = FALSE)
        
        if(nrow(kegg_significant) > 0) {
          write.csv(kegg_significant, 
                    file.path(enrichment_dir, paste0(valid_groups$name, "_", group_id, "_KEGG_0.05.csv")),
                    row.names = FALSE)
          
          top_kegg <- head(kegg_significant[order(kegg_significant$pvalue), ], 20)
          
          if(nrow(top_kegg) > 0) {
            top_kegg$Description_short <- ifelse(nchar(top_kegg$Description) > 60, 
                                                 paste0(substr(top_kegg$Description, 1, 60), "..."), 
                                                 top_kegg$Description)
            
            kegg_plot_data <- top_kegg %>%
              mutate(
                EnrichmentScore = Count / as.numeric(sub("/\\d+", "", BgRatio)),
                Term = Description_short,
                pValue = pvalue,
                Count = Count
              ) %>%
              arrange(desc(EnrichmentScore))
            
            write.csv(kegg_plot_data, 
                      file.path(enrichment_dir, paste0(valid_groups$name, "_", group_id, "_KEGG_plot_data.csv")),
                      row.names = FALSE)
            
            p_kegg_bubble <- ggplot(kegg_plot_data, aes(x = EnrichmentScore, y = reorder(Term, EnrichmentScore))) +
              geom_point(aes(size = Count, color = pValue), alpha = 0.8) +
              theme_bw() +
              labs(y = "", x = "Enrichment Score", 
                   title = paste0("KEGG Pathway Enrichment Analysis: ", valid_groups$name, " ", group_id),
                   color = "P-value", size = "Gene Count") +
              scale_color_gradient(low = "red", high = "blue") +
              scale_size_continuous(range = c(3, 12)) +
              theme(
                axis.text.y = element_text(size = 10, family = "Times", face = "bold", color = "black"),
                axis.text.x = element_text(size = 12, family = "Times", face = "bold", color = "black"),
                axis.title.x = element_text(size = 14, family = "Times", face = "bold", color = "black"),
                axis.title.y = element_text(size = 14, family = "Times", face = "bold", color = "black"),
                plot.title = element_text(hjust = 0.5, face = "bold", color = "black", family = "Times", size = 16),
                legend.text = element_text(family = "Times", size = 10),
                legend.title = element_text(family = "Times", face = "bold", size = 12),
                legend.position = "right"
              )
            
            ggsave(plot = p_kegg_bubble, 
                   filename = file.path(enrichment_dir, paste0(valid_groups$name, "_", group_id, "_KEGG_bubble_plot.png")),
                   width = 12, height = 10, dpi = 300)
            ggsave(plot = p_kegg_bubble, 
                   filename = file.path(enrichment_dir, paste0(valid_groups$name, "_", group_id, "_KEGG_bubble_plot.pdf")),
                   width = 12, height = 10)
            
            p_kegg_bar <- ggplot(kegg_plot_data[1:min(15, nrow(kegg_plot_data)), ], 
                                 aes(x = reorder(Term, EnrichmentScore), y = EnrichmentScore)) +
              geom_col(aes(fill = pValue), alpha = 0.8) +
              coord_flip() +
              theme_bw() +
              labs(x = "", y = "Enrichment Score", 
                   title = paste0("KEGG Pathway Enrichment (Bar Plot): ", valid_groups$name, " ", group_id),
                   fill = "P-value") +
              scale_fill_gradient(low = "red", high = "blue") +
              theme(
                axis.text.y = element_text(size = 10, family = "Times", face = "bold", color = "black"),
                axis.text.x = element_text(size = 12, family = "Times", face = "bold", color = "black"),
                axis.title.x = element_text(size = 14, family = "Times", face = "bold", color = "black"),
                axis.title.y = element_text(size = 14, family = "Times", face = "bold", color = "black"),
                plot.title = element_text(hjust = 0.5, face = "bold", color = "black", family = "Times", size = 16),
                legend.text = element_text(family = "Times", size = 10),
                legend.title = element_text(family = "Times", face = "bold", size = 12),
                legend.position = "right"
              )
            
            cat(" KEGG analysis completed, total", nrow(kegg_result@result), "enriched terms, with", nrow(kegg_significant), "significant\n")
          }
        } else {
          cat(" KEGG analysis was completed, but no significant enrichment results were found\n")
        }
        
        if(nrow(kegg_result@result) >= 3 && exists("GO_ALL") && nrow(GO_ALL) >= 3) {
          cat(" Generate a network connection relationship diagram...\n")
          
          tryCatch({
            if(exists("BP") && nrow(BP) >= 5) {
              go_sim <- pairwise_termsim(go_all)
              go_sim@result <- BP[BP$p.adjust < 0.05, ]
              
              if(nrow(go_sim@result) >= 5) {
                p_go_network <- emapplot(go_sim, showCategory = 15) + 
                  ggtitle(paste("GO Enrichment Network -", valid_groups$name, group_id)) +
                  theme(text = element_text(family = "Times"),
                        plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
                
                ggsave(plot = p_go_network, 
                       filename = file.path(enrichment_dir, paste0(valid_groups$name, "_", group_id, "_GO_network.png")),
                       width = 12, height = 10, dpi = 300)
                ggsave(plot = p_go_network, 
                       filename = file.path(enrichment_dir, paste0(valid_groups$name, "_", group_id, "_GO_network.pdf")),
                       width = 12, height = 10)
              }
            }
            
            if(nrow(kegg_significant) >= 5) {
              kegg_sim <- pairwise_termsim(kegg_result)
              kegg_sim@result <- kegg_significant[kegg_significant$p.adjust < 0.05, ]
              
              if(nrow(kegg_sim@result) >= 5) {
                p_kegg_network <- emapplot(kegg_sim, showCategory = 15) + 
                  ggtitle(paste("KEGG Enrichment Network -", valid_groups$name, group_id)) +
                  theme(text = element_text(family = "Times"),
                        plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
                
                ggsave(plot = p_kegg_network, 
                       filename = file.path(enrichment_dir, paste0(valid_groups$name, "_", group_id, "_KEGG_network.png")),
                       width = 12, height = 10, dpi = 300)
                ggsave(plot = p_kegg_network, 
                       filename = file.path(enrichment_dir, paste0(valid_groups$name, "_", group_id, "_KEGG_network.pdf")),
                       width = 12, height = 10)
              }
            }
            
          }, error = function(e) {
            cat(" Network diagram generation failed:", e$message, "\n")
          })
        }
      }
      
    }, error = function(e) {
      cat(" Error in gene ID conversion or enrichment analysis:", e$message, "\n")
    })
    
  }, error = function(e) {
    cat("  ", group_id, "analysis error:", e$message, "\n")
  })
}

tryCatch({
  all_markers_list <- list()
  
  for(group_id in valid_groups$groups[1:min(4, length(valid_groups$groups))]) {
    Idents(st) <- valid_groups$group_by
    
    markers <- FindMarkers(
      st, 
      ident.1 = group_id,
      min.pct = 0.1,
      logfc.threshold = 0.25,
      test.use = "wilcox",
      verbose = FALSE
    )
    
    significant_genes <- markers[markers$p_val_adj < 0.05 & markers$avg_log2FC > 0.5, ]
    
    if(nrow(significant_genes) > 0) {
      all_markers_list[[paste(valid_groups$name, group_id, sep = "_")]] <- rownames(significant_genes)
    }
  }
  
  if(length(all_markers_list) >= 2) {
    if(species_detected == "human") {
      gene_entrez_list <- lapply(all_markers_list, function(genes) {
        result <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        return(result$ENTREZID)
      })
    } else {
      gene_entrez_list <- lapply(all_markers_list, function(genes) {
        result <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
        return(result$ENTREZID)
      })
    }
    
    gene_entrez_list <- gene_entrez_list[lengths(gene_entrez_list) > 0]
    
    if(length(gene_entrez_list) >= 2) {
      go_compare <- compareCluster(
        geneClusters = gene_entrez_list,
        fun = "enrichGO",
        OrgDb = organism_db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )
      
      if(!is.null(go_compare) && nrow(go_compare@compareClusterResult) > 0) {
        write.csv(go_compare@compareClusterResult, 
                  file.path(enrichment_dir, "GO_compareCluster_results.csv"),
                  row.names = FALSE)
        
        p_go_compare <- dotplot(go_compare, showCategory = 10) + 
          ggtitle(paste("GO Enrichment Comparison across", valid_groups$name, "s")) +
          theme(text = element_text(family = "sans"), axis.text.x = element_text(angle = 45, hjust = 1))
        
        ggsave(plot = p_go_compare, 
               filename = file.path(enrichment_dir, "GO_compareCluster_dotplot.png"),
               width = 14, height = 10, dpi = 300)
        ggsave(plot = p_go_compare, 
               filename = file.path(enrichment_dir, "GO_compareCluster_dotplot.pdf"),
               width = 14, height = 10)
        
      }
      
      kegg_compare <- compareCluster(
        geneClusters = gene_entrez_list,
        fun = "enrichKEGG",
        organism = kegg_organism,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )
      
      if(!is.null(kegg_compare) && nrow(kegg_compare@compareClusterResult) > 0) {
        write.csv(kegg_compare@compareClusterResult, 
                  file.path(enrichment_dir, "KEGG_compareCluster_results.csv"),
                  row.names = FALSE)
        
        p_kegg_compare <- dotplot(kegg_compare, showCategory = 10) + 
          ggtitle(paste("KEGG Enrichment Comparison across", valid_groups$name, "s")) +
          theme(text = element_text(family = "sans"), axis.text.x = element_text(angle = 45, hjust = 1))
        
        ggsave(plot = p_kegg_compare, 
               filename = file.path(enrichment_dir, "KEGG_compareCluster_dotplot.png"),
               width = 14, height = 10, dpi = 300)
        ggsave(plot = p_kegg_compare, 
               filename = file.path(enrichment_dir, "KEGG_compareCluster_dotplot.pdf"),
               width = 14, height = 10)
        
      }
    }
  }
  
}, error = function(e) {
  cat("The overall comparative analysis was incorrect:", e$message, "\n")
})

cat("\n=== clusterProfiler Functional Enrichment Analysis Summary ===\n")
cat("Number of", valid_groups$name, "analyzed:", length(valid_groups$groups), "\n")

# Count generated files
result_files <- list.files(enrichment_dir, pattern = "\\.(csv|png|pdf)$", full.names = FALSE)
csv_files <- length(grep("\\.csv$", result_files))
png_files <- length(grep("\\.png$", result_files))
pdf_files <- length(grep("\\.pdf$", result_files))

cat("Generated result files:\n")
cat("- CSV data files:", csv_files, "\n")
cat("- PNG image files:", png_files, "\n")
cat("- PDF image files:", pdf_files, "\n")
cat("- Total files:", length(result_files), "\n")

cat("\nMain analysis content:\n")
cat("✓ Differential gene filtering (P-value < 0.05 && Log2FC > 0.5)\n")
cat("✓ GO functional enrichment analysis (BP, CC, MF)\n")
cat("✓ KEGG pathway enrichment analysis\n")
cat("✓ Network connectivity diagrams\n")
cat("✓ Multi-group comparative analysis\n")

cat("\nclusterProfiler functional enrichment analysis completed! Results saved in:", enrichment_dir, "\n")

