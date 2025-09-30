/*
 * Functional Enrichment Analysis Module - Based on clusterProfiler
 */

process ENRICHMENT_ANALYSIS {
    publishDir "${params.output_dir}/7.Enrichment", mode: 'copy'
    
    input:
    path seurat_rds
    path markers_csv
    
    output:
    path "*.csv"
    path "*.png" 
    path "*.pdf"
    path "enrichment_summary.txt", emit: summary
    
    script:
    """
    #!/usr/bin/Rscript
    
    # Load required packages
    suppressMessages({
        library(Seurat)
        
        # Conditionally load clusterProfiler related packages
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
        if (!require(DOSE, quietly = TRUE)) {
            BiocManager::install("DOSE")
            library(DOSE)
        }
        
        library(ggplot2)
        library(dplyr)
        library(tidyr)
    })
    
    cat("=== Starting Functional Enrichment Analysis ===\\n")
    
    # Load data
    cat("Loading annotated data...\\n")
    st <- readRDS("${seurat_rds}")
    
    # Check for marker gene data
    markers_file <- "${markers_csv}"
    if(file.exists(markers_file) && file.size(markers_file) > 10) {
        all_markers <- read.csv(markers_file, stringsAsFactors = FALSE)
        cat("Loading marker gene data, total", nrow(all_markers), "genes\\n")
    } else {
        cat("No valid marker gene data found, exiting enrichment analysis\\n")
        # Create empty output files
        empty_df <- data.frame()
        write.csv(empty_df, "enrichment_results_GO.csv", row.names = FALSE)
        write.csv(empty_df, "enrichment_results_KEGG.csv", row.names = FALSE)
        
        png("enrichment_no_markers.png", width = 8, height = 6, units = "in", res = 300)
        plot.new()
        text(0.5, 0.5, "No available marker genes for enrichment analysis", cex = 1.2, col = "orange")
        dev.off()
        
        writeLines("No available marker genes for enrichment analysis", "enrichment_summary.txt")
        quit("no", status = 0)
    }
    
    tryCatch({
        # Species setting - use parameter specification only
        if("${params.species}" == "human") {
            species_detected <- "human"
            cat("Using species: human\\n")
        } else if("${params.species}" == "mouse" || "${params.species}" == "mmu") {
            species_detected <- "mouse"
            cat("Using species: mouse (mmu)\\n")
        } else {
            # Default to mouse if parameter is unclear
            species_detected <- "mouse"
            cat("Parameter unclear, defaulting to mouse\\n")
        }
        
        # Set species-related databases
        if(species_detected == "human") {
            organism_db <- org.Hs.eg.db
            kegg_organism <- "hsa"
            cat("Using human annotation database\\n")
        } else {
            organism_db <- org.Mm.eg.db
            kegg_organism <- "mmu"
            cat("Using mouse annotation database\\n")
        }
        
        # Get valid cell types or clusters (consistent with R script logic)
        if("cell_type_filtered" %in% colnames(st@meta.data)) {
            celltype_counts <- table(st@meta.data\$cell_type_filtered)
            valid_celltypes <- names(celltype_counts)[celltype_counts >= 20]
            valid_celltypes <- valid_celltypes[valid_celltypes != "Unknown"]
            
            if(length(valid_celltypes) >= 2) {
                valid_groups <- list(groups = valid_celltypes, group_by = "cell_type_filtered", name = "celltype")
            } else {
                cluster_table <- table(st@meta.data\$seurat_clusters)
                valid_clusters <- names(cluster_table)[cluster_table >= 20]
                valid_groups <- list(groups = valid_clusters[1:min(6, length(valid_clusters))], group_by = "seurat_clusters", name = "cluster")
            }
        } else {
            cluster_table <- table(st@meta.data\$seurat_clusters)
            valid_clusters <- names(cluster_table)[cluster_table >= 20]
            valid_groups <- list(groups = valid_clusters[1:min(6, length(valid_clusters))], group_by = "seurat_clusters", name = "cluster")
        }
        
        cat("Groups for differential gene analysis:", length(valid_groups\$groups), valid_groups\$name, "\\n")
        
        # Record cluster to cell type mapping (for display)
        if("cell_type_filtered" %in% colnames(st@meta.data)) {
            cluster_celltype_map <- table(st@meta.data\$seurat_clusters, st@meta.data\$cell_type_filtered)
            cat("Cluster to cell type mapping:\\\\n")
            print(cluster_celltype_map)
        }
        
        # Perform enrichment analysis for each group
        enrichment_results <- list()
        summary_info <- c("=== Functional Enrichment Analysis Summary ===", 
                         paste("Analysis species:", species_detected),
                         paste("Number of analysis groups:", length(valid_groups\$groups)),
                         paste("Grouping method:", valid_groups\$group_by), "")
        
        # Perform differential gene analysis and functional enrichment for each cell cluster
        for(i in 1:length(valid_groups\$groups)) {
            group_id <- valid_groups\$groups[i]
            cat(paste("\\nAnalyzing", valid_groups\$name, ":", group_id, "\\n"))
            
            tryCatch({
                # Set comparison groups and find marker genes
                Idents(st) <- valid_groups\$group_by
                
                # Find marker genes (compared to all other groups)
                cat("  Finding marker genes...\\n")
                markers <- FindMarkers(
                    st, 
                    ident.1 = group_id,
                    min.pct = 0.1,
                    logfc.threshold = 0.25,
                    test.use = "wilcox",
                    verbose = FALSE
                )
                
                if(nrow(markers) == 0) {
                    cat("  No differential genes found, skipping\\n")
                    next
                }
                
                # Filter significant differential genes (P-value < 0.05 && Log2FC > 0.5)
                significant_genes <- markers[markers\$p_val_adj < 0.05 & markers\$avg_log2FC > 0.5, ]
                
                if(nrow(significant_genes) == 0) {
                    cat("  No significant differential genes found, skipping\\n")
                    next
                }
                
                cat("  Found", nrow(significant_genes), "significant differential genes\\n")
                
                # Save differential gene list
                significant_genes\$gene <- rownames(significant_genes)
                write.csv(significant_genes, 
                         paste0(valid_groups\$name, "_", group_id, "_DEGs.csv"),
                         row.names = TRUE)
                
                # Gene ID conversion
                gene_list <- rownames(significant_genes)
                
                # Display first 10 gene names for debugging
                cat("  Differential gene examples:", paste(head(gene_list, 10), collapse = ", "), "\\n")
                
                # Convert gene symbols to ENTREZ ID
                gene_entrez <- tryCatch({
                    if(species_detected == "human") {
                        result <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
                        cat("  Using human annotation database to convert gene IDs\\n")
                        result
                    } else {
                        result <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
                        cat("  Using mouse annotation database to convert gene IDs\\n")
                        result
                    }
                }, error = function(e) {
                    cat("  Gene ID conversion error:", e\$message, "\\n")
                    data.frame()
                })
                
                if(is.null(gene_entrez) || nrow(gene_entrez) == 0) {
                    cat("  Gene ID conversion failed, skipping enrichment analysis\\n")
                    cat("  Conversion failure rate:", round((length(gene_list) - nrow(gene_entrez))/length(gene_list)*100, 1), "%\\n")
                    next
                }
                
                conversion_rate <- round(nrow(gene_entrez)/length(gene_list)*100, 1)
                cat("  Successfully converted", nrow(gene_entrez), "gene IDs (", conversion_rate, "%)\\n")
                
                # If conversion rate is too low, show warning but continue analysis
                if(conversion_rate < 10) {
                    cat("  Warning: Gene ID conversion rate too low (", conversion_rate, "%), enrichment results may be unreliable\\n")
                }
                
                # GO enrichment analysis (overall analysis)
                cat("  Performing GO enrichment analysis...\\n")
                
                # Analyze all GO categories together
                go_all <- enrichGO(
                    gene = gene_entrez\$ENTREZID,
                    OrgDb = organism_db,
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    readable = TRUE
                )
                
                if(!is.null(go_all) && nrow(go_all@result) > 0) {
                    # Filter significant results
                    GO_ALL <- go_all@result
                    GO_ALL <- GO_ALL[GO_ALL\$pvalue < 0.05,]
                    
                    if(nrow(GO_ALL) > 0) {
                        # 分别提取BP、CC、MF各类别结果
                        BP <- GO_ALL[GO_ALL\$ONTOLOGY=='BP', ]
                        CC <- GO_ALL[GO_ALL\$ONTOLOGY=='CC', ]
                        MF <- GO_ALL[GO_ALL\$ONTOLOGY=='MF', ]
                        
                        # 保存各类别结果
                        if(nrow(BP) > 0) {
                            write.csv(BP, 
                                     paste0(valid_groups\$name, "_", group_id, "_GO_BP.csv"),
                                     row.names = FALSE)
                            cat("    GO BP analysis completed,", nrow(BP), "enriched terms\\n")
                        }
                        
                        if(nrow(CC) > 0) {
                            write.csv(CC, 
                                     paste0(valid_groups\$name, "_", group_id, "_GO_CC.csv"),
                                     row.names = FALSE)
                            cat("    GO CC analysis completed,", nrow(CC), "enriched terms\\n")
                        }
                        
                        if(nrow(MF) > 0) {
                            write.csv(MF, 
                                     paste0(valid_groups\$name, "_", group_id, "_GO_MF.csv"),
                                     row.names = FALSE)
                            cat("    GO MF analysis completed,", nrow(MF), "enriched terms\\n")
                        }
                        
                        # 保存全部GO结果
                        write.csv(GO_ALL, 
                                 paste0(valid_groups\$name, "_", group_id, "_GO_ALL.csv"),
                                 row.names = FALSE)
                        
                        # 生成合并的GO气泡图
                        if(nrow(BP) > 0 || nrow(CC) > 0 || nrow(MF) > 0) {
                            # 设置每个类别显示的条目数
                            display_number <- c(10, 10, 10)  # BP、CC、MF各显示10条
                            
                            ego_result_BP <- if(nrow(BP) > 0) head(BP[order(BP\$pvalue), ], display_number[1]) else data.frame()
                            ego_result_CC <- if(nrow(CC) > 0) head(CC[order(CC\$pvalue), ], display_number[2]) else data.frame()
                            ego_result_MF <- if(nrow(MF) > 0) head(MF[order(MF\$pvalue), ], display_number[3]) else data.frame()
                            
                            # 构建绘图数据框
                            if(nrow(ego_result_BP) > 0 || nrow(ego_result_CC) > 0 || nrow(ego_result_MF) > 0) {
                                go_enrich_df <- data.frame(
                                    ID = c(ego_result_BP\$ID, ego_result_CC\$ID, ego_result_MF\$ID),
                                    Description = c(ego_result_BP\$Description, ego_result_CC\$Description, ego_result_MF\$Description),
                                    GeneNumber = c(ego_result_BP\$Count, ego_result_CC\$Count, ego_result_MF\$Count),
                                    pvalue = c(ego_result_BP\$pvalue, ego_result_CC\$pvalue, ego_result_MF\$pvalue),
                                    type = factor(c(rep("BP:biological process", nrow(ego_result_BP)), 
                                                   rep("CC:cellular component", nrow(ego_result_CC)),
                                                   rep("MF:molecular function", nrow(ego_result_MF))), 
                                                 levels = c("BP:biological process", "CC:cellular component", "MF:molecular function"))
                                )
                                
                                # 处理过长的条目名称，截断超过60个字符的描述
                                go_enrich_df\$Description_short <- ifelse(nchar(go_enrich_df\$Description) > 60, 
                                                                        paste0(substr(go_enrich_df\$Description, 1, 60), "..."), 
                                                                        go_enrich_df\$Description)
                                
                                # 设置绘图顺序
                                go_enrich_df\$type_order <- factor(rev(as.integer(rownames(go_enrich_df))), 
                                                                  labels = rev(go_enrich_df\$Description_short))
                                
                                # 设定颜色
                                COLS <- c("#a8d8ea", "#FFDD94", "#FF7B89")
                                
                                # 创建GO气泡图
                                p_go_bubble <- ggplot(data = go_enrich_df, aes(x = type_order, y = GeneNumber, size = GeneNumber, fill = type)) +
                                    geom_point(shape = 21, color = "black") +
                                    scale_size_area(max_size = 10) +
                                    scale_fill_manual(values = COLS) +
                                    coord_flip() +
                                    xlab("GO term") + 
                                    ylab("Gene Number") + 
                                    labs(title = paste0("GO Enrichment Analysis: ", valid_groups\$name, " ", group_id)) +
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
                                
                                # 保存GO气泡图
                                ggsave(plot = p_go_bubble, 
                                       filename = paste0(valid_groups\$name, "_", group_id, "_GO_bubble_plot.png"),
                                       width = 16, height = 12, dpi = 300)
                                ggsave(plot = p_go_bubble, 
                                       filename = paste0(valid_groups\$name, "_", group_id, "_GO_bubble_plot.pdf"),
                                       width = 16, height = 12)
                                
                                cat("    GO气泡图生成完成\\n")
                            }
                        }
                    } else {
                        cat("    No significant GO enrichment results found\\n")
                    }
                } else {
                    cat("    GO enrichment analysis failed or no results\\n")
                }
                
                # KEGG enrichment analysis
                cat("  Performing KEGG enrichment analysis...\\n")
                kegg_result <- enrichKEGG(
                    gene = gene_entrez\$ENTREZID,
                    organism = kegg_organism,
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 1
                )
                
                if(!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
                    # 添加基因符号
                    kegg_readable <- setReadable(kegg_result, organism_db, keyType = "ENTREZID")
                    
                    # Filter significant results
                    kegg_significant <- kegg_result@result[kegg_result@result\$pvalue <= 0.05, ]
                    
                    # 保存所有KEGG结果
                    write.csv(kegg_result@result, 
                             paste0(valid_groups\$name, "_", group_id, "_KEGG_all.csv"),
                             row.names = FALSE)
                    
                    # 保存显著KEGG结果
                    if(nrow(kegg_significant) > 0) {
                        write.csv(kegg_significant, 
                                 paste0(valid_groups\$name, "_", group_id, "_KEGG_0.05.csv"),
                                 row.names = FALSE)
                        
                        # 取前20个最显著的结果进行可视化
                        top_kegg <- head(kegg_significant[order(kegg_significant\$pvalue), ], 20)
                        
                        if(nrow(top_kegg) > 0) {
                            # 处理过长的条目名称
                            top_kegg\$Description_short <- ifelse(nchar(top_kegg\$Description) > 60, 
                                                                paste0(substr(top_kegg\$Description, 1, 60), "..."), 
                                                                top_kegg\$Description)
                            
                            # 创建用于气泡图的数据框
                            kegg_plot_data <- top_kegg %>%
                                mutate(
                                    EnrichmentScore = Count / as.numeric(sub("/\\\\d+", "", BgRatio)),
                                    Term = Description_short,
                                    pValue = pvalue,
                                    Count = Count
                                ) %>%
                                arrange(desc(EnrichmentScore))
                            
                            # 保存绘图数据
                            write.csv(kegg_plot_data, 
                                     paste0(valid_groups\$name, "_", group_id, "_KEGG_plot_data.csv"),
                                     row.names = FALSE)
                            
                            # 生成KEGG气泡图
                            p_kegg_bubble <- ggplot(kegg_plot_data, aes(x = EnrichmentScore, y = reorder(Term, EnrichmentScore))) +
                                geom_point(aes(size = Count, color = pValue), alpha = 0.8) +
                                theme_bw() +
                                labs(y = "", x = "Enrichment Score", 
                                     title = paste0("KEGG Pathway Enrichment Analysis: ", valid_groups\$name, " ", group_id),
                                     color = "P-value", size = "Gene Count") +
                                scale_color_gradient(low = "red", high = "blue") +
                                scale_size_continuous(range = c(3, 12)) +
                                theme(
                                    legend.background = element_rect(fill = "white", color = "black", size = 0.2),
                                    legend.text = element_text(face = "bold", color = "black", family = "Times", size = 12),
                                    plot.title = element_text(hjust = 0.5, face = "bold", color = "black", family = "Times", size = 18),
                                    axis.text.x = element_text(family = "Times", face = "bold", color = "black", size = 14),
                                    axis.text.y = element_text(family = "Times", face = "bold", color = "black", size = 10),
                                    axis.title.x = element_text(face = "bold", color = "black", family = "Times", size = 18),
                                    axis.title.y = element_text(face = "bold", color = "black", family = "Times", size = 18),
                                    plot.subtitle = element_text(hjust = 0.5, family = "Times", size = 14, face = "italic", colour = "black"),
                                    legend.position = "right"
                                )
                            
                            # 保存KEGG气泡图
                            ggsave(plot = p_kegg_bubble, 
                                   filename = paste0(valid_groups\$name, "_", group_id, "_KEGG_bubble_plot.png"),
                                   width = 12, height = 10, dpi = 300)
                            ggsave(plot = p_kegg_bubble, 
                                   filename = paste0(valid_groups\$name, "_", group_id, "_KEGG_bubble_plot.pdf"),
                                   width = 12, height = 10)
                        
                            cat("    KEGG analysis completed,", nrow(kegg_result@result), "enriched terms, of which", nrow(kegg_significant), "are significant\\n")
                        }
                    } else {
                        cat("    KEGG analysis completed, but no significant enrichment results found\\n")
                    }
                    
                    # 生成网络连接关系图
                    if(nrow(kegg_result@result) >= 3 && exists("GO_ALL") && nrow(GO_ALL) >= 3) {
                        cat("  生成网络连接关系图...\\n")
                        
                        tryCatch({
                            # GO网络图 (使用BP类别)
                            if(exists("BP") && nrow(BP) >= 5) {
                                tryCatch({
                                    # 使用pairwise_termsim计算相似性矩阵
                                    go_sim <- pairwise_termsim(go_all)
                                    # 使用所有BP条目（不再额外筛选p.adjust）
                                    go_sim@result <- BP
                                    
                                    if(nrow(go_sim@result) >= 5) {
                                        p_go_network <- emapplot(go_sim, showCategory = min(15, nrow(BP))) + 
                                            ggtitle(paste("GO Enrichment Network -", valid_groups\$name, group_id)) +
                                            theme(
                                                legend.background = element_blank(),
                                                legend.text = element_text(face = "bold", color = "black", family = "Times", size = 12),
                                                plot.title = element_text(hjust = 0.5, face = "bold", color = "black", family = "Times", size = 18),
                                                axis.text.x = element_blank(),
                                                axis.text.y = element_blank(),
                                                axis.title.x = element_blank(),
                                                axis.title.y = element_blank(),
                                                axis.ticks = element_blank(),
                                                plot.subtitle = element_text(hjust = 0.5, family = "Times", size = 14, face = "italic", colour = "black")
                                            )
                                        
                                        ggsave(plot = p_go_network, 
                                               filename = paste0(valid_groups\$name, "_", group_id, "_GO_network.png"),
                                               width = 12, height = 10, dpi = 300)
                                        ggsave(plot = p_go_network, 
                                               filename = paste0(valid_groups\$name, "_", group_id, "_GO_network.pdf"),
                                               width = 12, height = 10)
                                        
                                        cat("    GO网络图生成完成\\n")
                                    } else {
                                        cat("    GO BP条目不足，跳过网络图\\n")
                                    }
                                }, error = function(e) {
                                    cat("    GO网络图生成失败：", e\$message, "\\n")
                                })
                            } else {
                                cat("    GO BP结果不足，跳过网络图\\n")
                            }
                            
                            # KEGG网络图
                            if(nrow(kegg_significant) >= 3) {
                                tryCatch({
                                    # 使用pairwise_termsim计算相似性矩阵
                                    kegg_sim <- pairwise_termsim(kegg_result)
                                    
                                    # 使用所有显著条目（不再额外筛选p.adjust）
                                    kegg_sim@result <- kegg_significant
                                    
                                    if(nrow(kegg_sim@result) >= 3) {
                                        p_kegg_network <- emapplot(kegg_sim, showCategory = min(15, nrow(kegg_significant))) + 
                                            ggtitle(paste("KEGG Enrichment Network -", valid_groups\$name, group_id)) +
                                            theme(
                                                legend.background = element_blank(),
                                                legend.text = element_text(face = "bold", color = "black", family = "Times", size = 12),
                                                plot.title = element_text(hjust = 0.5, face = "bold", color = "black", family = "Times", size = 18),
                                                axis.text.x = element_blank(),
                                                axis.text.y = element_blank(),
                                                axis.title.x = element_blank(),
                                                axis.title.y = element_blank(),
                                                axis.ticks = element_blank(),
                                                plot.subtitle = element_text(hjust = 0.5, family = "Times", size = 14, face = "italic", colour = "black")
                                            )
                                        
                                        ggsave(plot = p_kegg_network, 
                                               filename = paste0(valid_groups\$name, "_", group_id, "_KEGG_network.png"),
                                               width = 12, height = 10, dpi = 300)
                                        ggsave(plot = p_kegg_network, 
                                               filename = paste0(valid_groups\$name, "_", group_id, "_KEGG_network.pdf"),
                                               width = 12, height = 10)
                                        
                                        cat("    KEGG网络图生成完成\\n")
                                    } else {
                                        cat("    KEGG显著条目不足，跳过网络图\\n")
                                    }
                                }, error = function(e) {
                                    cat("    KEGG网络图生成失败：", e\$message, "\\n")
                                })
                            } else {
                                cat("    KEGG富集结果不足，跳过网络图\\n")
                            }
                            
                        }, error = function(e) {
                            cat("    整体网络图生成失败：", e\$message, "\\n")
                        })
                    }
                }
                
                enrichment_results[[group_id]] <- list(go = go_all, kegg = kegg_result)
                
            }, error = function(e) {
                cat("  ", group_id, "分析出错：", e\$message, "\\n")
            })
        }
        
        # 生成整体比较分析
        cat("\\n进行整体功能富集比较分析...\\n")
        
        tryCatch({
            # 收集所有组的差异基因进行比较
            all_markers_list <- list()
            
            for(group_id in valid_groups\$groups[1:min(4, length(valid_groups\$groups))]) {
                Idents(st) <- valid_groups\$group_by
                
                markers <- FindMarkers(
                    st, 
                    ident.1 = group_id,
                    min.pct = 0.1,
                    logfc.threshold = 0.25,
                    test.use = "wilcox",
                    verbose = FALSE
                )
                
                significant_genes <- markers[markers\$p_val_adj < 0.05 & markers\$avg_log2FC > 0.5, ]
                
                if(nrow(significant_genes) > 0) {
                    all_markers_list[[paste(valid_groups\$name, group_id, sep = "_")]] <- rownames(significant_genes)
                }
            }
            
            if(length(all_markers_list) >= 2) {
                # 使用compareCluster进行比较分析
                if(species_detected == "human") {
                    gene_entrez_list <- lapply(all_markers_list, function(genes) {
                        tryCatch({
                            result <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
                            if(nrow(result) > 0) {
                                return(result\$ENTREZID)
                            } else {
                                return(character(0))
                            }
                        }, error = function(e) {
                            return(character(0))
                        })
                    })
                } else {
                    gene_entrez_list <- lapply(all_markers_list, function(genes) {
                        tryCatch({
                            result <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
                            if(nrow(result) > 0) {
                                return(result\$ENTREZID)
                            } else {
                                return(character(0))
                            }
                        }, error = function(e) {
                            return(character(0))
                        })
                    })
                }
                
                # 过滤掉空列表
                gene_entrez_list <- gene_entrez_list[lengths(gene_entrez_list) > 0]
                
                if(length(gene_entrez_list) >= 2) {
                    # GO比较分析
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
                        # 保存比较结果
                        write.csv(go_compare@compareClusterResult, 
                                 "GO_compareCluster_results.csv",
                                 row.names = FALSE)
                        
                        # 比较图
                        p_go_compare <- dotplot(go_compare, showCategory = 10) + 
                            ggtitle(paste("GO Enrichment Comparison across", valid_groups\$name)) +
                            theme(
                                legend.background = element_rect(fill = "white", color = "black", size = 0.2),
                                legend.text = element_text(face = "bold", color = "black", family = "Times", size = 12),
                                plot.title = element_text(hjust = 0.5, face = "bold", color = "black", family = "Times", size = 18),
                                axis.text.x = element_text(family = "Times", face = "bold", color = "black", size = 14, angle = 45, hjust = 1),
                                axis.text.y = element_text(family = "Times", face = "bold", color = "black", size = 10),
                                axis.title.x = element_text(face = "bold", color = "black", family = "Times", size = 18),
                                axis.title.y = element_text(face = "bold", color = "black", family = "Times", size = 18),
                                plot.subtitle = element_text(hjust = 0.5, family = "Times", size = 14, face = "italic", colour = "black")
                            )
                        
                        ggsave(plot = p_go_compare, 
                               filename = "GO_compareCluster_dotplot.png",
                               width = 14, height = 16, dpi = 300)
                        ggsave(plot = p_go_compare, 
                               filename = "GO_compareCluster_dotplot.pdf",
                               width = 14, height = 16)
                        
                        cat("GO比较分析完成\\n")
                    }
                    
                    # KEGG比较分析
                    kegg_compare <- compareCluster(
                        geneClusters = gene_entrez_list,
                        fun = "enrichKEGG",
                        organism = kegg_organism,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2
                    )
                    
                    if(!is.null(kegg_compare) && nrow(kegg_compare@compareClusterResult) > 0) {
                        # 保存比较结果
                        write.csv(kegg_compare@compareClusterResult, 
                                 "KEGG_compareCluster_results.csv",
                                 row.names = FALSE)
                        
                        # 比较图
                        p_kegg_compare <- dotplot(kegg_compare, showCategory = 10) + 
                            ggtitle(paste("KEGG Enrichment Comparison across", valid_groups\$name)) +
                            theme(
                                legend.background = element_rect(fill = "white", color = "black", size = 0.2),
                                legend.text = element_text(face = "bold", color = "black", family = "Times", size = 12),
                                plot.title = element_text(hjust = 0.5, face = "bold", color = "black", family = "Times", size = 18),
                                axis.text.x = element_text(family = "Times", face = "bold", color = "black", size = 14, angle = 45, hjust = 1),
                                axis.text.y = element_text(family = "Times", face = "bold", color = "black", size = 10),
                                axis.title.x = element_text(face = "bold", color = "black", family = "Times", size = 18),
                                axis.title.y = element_text(face = "bold", color = "black", family = "Times", size = 18),
                                plot.subtitle = element_text(hjust = 0.5, family = "Times", size = 14, face = "italic", colour = "black")
                            )
                        
                        ggsave(plot = p_kegg_compare, 
                               filename = "KEGG_compareCluster_dotplot.png",
                               width = 14, height = 16, dpi = 300)
                        ggsave(plot = p_kegg_compare, 
                               filename = "KEGG_compareCluster_dotplot.pdf",
                               width = 14, height = 16)
                        
                        cat("KEGG比较分析完成\\n")
                    }
                }
            }
            
        }, error = function(e) {
            cat("整体比较分析出错：", e\$message, "\\n")
        })
        
        # 生成分析摘要
        cat("\\n=== clusterProfiler功能富集分析摘要 ===\\n")
        cat("分析的", valid_groups\$name, "数量：", length(valid_groups\$groups), "\\n")
        
        # 统计生成的文件
        result_files <- list.files(".", pattern = "\\\\.(csv|png|pdf)\$", full.names = FALSE)
        csv_files <- length(grep("\\\\.csv\$", result_files))
        png_files <- length(grep("\\\\.png\$", result_files))
        pdf_files <- length(grep("\\\\.pdf\$", result_files))
        
        cat("生成的结果文件:\\n")
        cat("- CSV数据文件：", csv_files, "个\\n")
        cat("- PNG图片文件：", png_files, "个\\n")
        cat("- PDF图片文件：", pdf_files, "个\\n")
        cat("- 总文件数：", length(result_files), "个\\n")
        
        cat("\\n主要分析内容:\\n")
        cat("✓ 差异基因筛选 (P-value < 0.05 && Log2FC > 0.5)\\n")
        cat("✓ GO功能富集分析 (BP, CC, MF)\\n")
        cat("✓ KEGG通路富集分析\\n")
        cat("✓ 网络连接关系图\\n")
        cat("✓ 多组比较分析\\n")
        
        # 生成summary文件
        summary_info <- c(
            "=== clusterProfiler功能富集分析摘要 ===",
            paste("分析的", valid_groups\$name, "数量:", length(valid_groups\$groups)),
            "",
            "生成的结果文件:",
            paste("- CSV数据文件:", csv_files, "个"),
            paste("- PNG图片文件:", png_files, "个"),
            paste("- PDF图片文件:", pdf_files, "个"),
            paste("- 总文件数:", length(result_files), "个"),
            "",
            "主要分析内容:",
            "✓ 差异基因筛选 (P-value < 0.05 && Log2FC > 0.5)",
            "✓ GO功能富集分析 (BP, CC, MF)",
            "✓ KEGG通路富集分析",
            "✓ 网络连接关系图",
            "✓ 多组比较分析"
        )
        
        writeLines(summary_info, "enrichment_summary.txt")
        
        cat("\\nclusterProfiler functional enrichment analysis completed!\\n")
        
    }, error = function(e) {
        cat("Functional enrichment analysis error:", conditionMessage(e), "\\n")
        
        # Create empty output files
        write.csv(data.frame(), "enrichment_results_GO.csv", row.names = FALSE)
        write.csv(data.frame(), "enrichment_results_KEGG.csv", row.names = FALSE)
        
        png("enrichment_analysis_error.png", width = 8, height = 6, units = "in", res = 300)
        plot.new()
        text(0.5, 0.5, paste("Functional enrichment analysis failed:\\\\n", conditionMessage(e)), 
             cex = 1.2, col = "red")
        dev.off()
        
        writeLines(c(
            "=== Functional Enrichment Analysis Error ===",
            paste("Error message:", conditionMessage(e)),
            "Possible causes:",
            "1. clusterProfiler related packages not installed",
            "2. Annotation database issues", 
            "3. Marker gene data format issues",
            "Please check error message and rerun"
        ), "enrichment_summary.txt")
    })
    """
}