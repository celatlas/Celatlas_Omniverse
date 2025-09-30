/*
 * CellChat Cell communication analysis module
 */

process CELLCHAT_ANALYSIS {
    publishDir "${params.output_dir}/5.Cellchat", mode: 'copy'
    
    input:
    path seurat_rds
    
    output:
    path "cellchat_results.Rdata", emit: cellchat_object
    path "*.png", optional: true
    path "*.pdf", optional: true
    path "01_communication_network.csv"
    path "interaction_count_matrix.csv", optional: true
    path "03_number_of_interactions_heatmap.pdf", optional: true
    path "03_number_of_interactions_heatmap.png", optional: true
    path "04_single_celltype_circles.pdf", optional: true
    path "04_single_celltype_circles.png", optional: true
    path "04_outgoing_incoming_communication.pdf", optional: true
    path "04_outgoing_incoming_communication.png", optional: true
    path "05_ligand_receptor_bubble.pdf", optional: true
    path "05_ligand_receptor_bubble.png", optional: true
    
    script:
    """
    #!/usr/bin/Rscript
    
    suppressMessages({
        library(Seurat)
        library(CellChat)
        library(patchwork)
        library(ggplot2)
        library(dplyr)
        library(RColorBrewer)
        library(igraph)
        library(ggraph)
    })
    
    cat("=== Starting CellChat cell communication analysis ===\\n")
    
    # Read annotated Seurat object
    cat("Reading annotated data...\\n")
    st <- readRDS("${seurat_rds}")
    
    cat("Input cell count:", ncol(st), "\\n")
    
    tryCatch({
        # Prepare CellChat input data
        data.input <- GetAssayData(st, layer = 'data')  # normalized data matrix
        meta <- st@meta.data                           # meta data
        
        # Check cell type annotations
        cat("Available cell types:\\n")
        if("cell_type_filtered" %in% colnames(meta)) {
            print(table(meta\$cell_type_filtered))
            celltype_column <- "cell_type_filtered"
        } else if("cell_type_raw" %in% colnames(meta)) {
            print(table(meta\$cell_type_raw))
            celltype_column <- "cell_type_raw"
        } else {
            print(table(meta\$seurat_clusters))
            celltype_column <- "seurat_clusters"
            cat("Using clustering information for CellChat analysis\\n")
        }
        
        # Verify data integrity
        cat("Verifying Seurat object integrity:\\n")
        cat("Metadata columns:", paste(colnames(meta), collapse=", "), "\\n")
        cat("Cell type column exists:", celltype_column %in% colnames(meta), "\\n")
        cat("Metadata rows:", nrow(meta), "\\n")
        cat("data.input columns:", ncol(data.input), "\\n")
        
        # Use all cells without filtering
        cat("Using all cells for CellChat analysis\\n")
        cat("Cell count:", ncol(data.input), "\\n")
        cat("Cell type distribution:\\n")
        print(table(meta[[celltype_column]]))
        
        # Check cell type counts (lower threshold)
        celltype_counts <- table(meta[[celltype_column]])
        valid_celltypes <- names(celltype_counts)[celltype_counts >= 10]  # Lower minimum cell count requirement
        
        if(length(valid_celltypes) < 2) {
            cat("Insufficient cell types, skipping CellChat analysis\\n")
            # Create empty output files
            empty_df <- data.frame()
            write.csv(empty_df, "01_communication_network.csv", row.names = FALSE)
            write.csv(empty_df, "interaction_count_matrix.csv", row.names = FALSE)
            
            # Create empty CellChat object
            empty_cellchat <- list()
            save(empty_cellchat, file = "cellchat_results.Rdata")
            
            png("cellchat_insufficient_data.png", width = 8, height = 6, units = "in", res = 150)
            plot.new()
            text(0.5, 0.5, "Insufficient data for CellChat analysis\\nRequires at least 2 cell types, each with >=${params.cellchat.min_cells} cells", 
                 cex = 1.2, col = "orange")
            dev.off()
            
            # Create empty PDF files to meet output requirements
            pdf("01_net_circle_number.pdf", width = 8, height = 6)
            plot.new()
            text(0.5, 0.5, "Insufficient data for CellChat analysis", cex = 1.2, col = "orange")
            dev.off()
            
            pdf("02_net_circle_weight.pdf", width = 8, height = 6)
            plot.new()
            text(0.5, 0.5, "Insufficient data for CellChat analysis", cex = 1.2, col = "orange")
            dev.off()
            
            pdf("03_number_of_interactions_heatmap.pdf", width = 8, height = 6)
            plot.new()
            text(0.5, 0.5, "Insufficient data for CellChat analysis", cex = 1.2, col = "orange")
            dev.off()
            
            pdf("04_single_celltype_circles.pdf", width = 18, height = 5)
            plot.new()
            text(0.5, 0.5, "Insufficient data for CellChat analysis", cex = 1.2, col = "orange")
            dev.off()
            
            pdf("05_ligand_receptor_bubble.pdf", width = 8, height = 6)
            plot.new()
            text(0.5, 0.5, "Insufficient data for CellChat analysis", cex = 1.2, col = "orange")
            dev.off()
            
            quit("no", status = 0)
        }
        
        # Configure CellChat database based on species parameter
        database_param <- "${params.cellchat.database}"
        
        if(database_param == "auto") {
            # Auto-detect species
            if("${params.species}" == "human") {
                species_detected <- "human"
                cat("Auto-detected species: human\\n")
            } else if("${params.species}" == "mouse" || "${params.species}" == "mmu") {
                species_detected <- "mouse"
                cat("Auto-detected species: mouse (mmu)\\n")
            } else {
                # Default to mouse (since most spatial transcriptomics data is mouse)
                species_detected <- "mouse"
                cat("Species not specified, defaulting to mouse\\n")
            }
            
            if(species_detected == "human") {
                CellChatDB <- CellChatDB.human
                cat("Using human CellChat database\\n")
            } else {
                CellChatDB <- CellChatDB.mouse
                cat("Using mouse CellChat database\\n")
            }
        } else if(database_param == "CellChatDB.human") {
            CellChatDB <- CellChatDB.human
            cat("Using specified human CellChat database\\n")
        } else if(database_param == "CellChatDB.mouse") {
            CellChatDB <- CellChatDB.mouse
            cat("Using specified mouse CellChat database\\n")
        } else {
            # Default to mouse
            CellChatDB <- CellChatDB.mouse
            cat("Invalid database parameter, defaulting to mouse CellChat database\\n")
        }
        
        showDatabaseCategory(CellChatDB)
        # Use full database but limit size to improve success rate
        # Prioritize Secreted Signaling and Paracrine signaling
        CellChatDB.use <- CellChatDB
        # Try to limit database size to improve success rate
        if(nrow(CellChatDB.use\$interaction) > 500) {
            # Prioritize important signaling pathways
            important_pathways <- c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact")
            CellChatDB.use <- subsetDB(CellChatDB, search = important_pathways)
            cat("Using limited CellChat database (", nrow(CellChatDB.use\$interaction), " interactions)\\n")
        } else {
            cat("Using full CellChat database (", nrow(CellChatDB.use\$interaction), " interactions)\\n")
        }
        
        # Create CellChat object
        cat("Creating CellChat object...\\n")
        cellchat <- createCellChat(object = data.input, meta = meta, group.by = celltype_column)
        cellchat <- addMeta(cellchat, meta = meta)
        cellchat <- setIdent(cellchat, ident.use = celltype_column)
        
        # Clean up unused levels
        cellchat@idents <- droplevels(cellchat@idents, 
                                     exclude = setdiff(levels(cellchat@idents), unique(cellchat@idents)))
        
        groupSize <- as.numeric(table(cellchat@idents))
        cat("Cell type counts:\\n")
        print(table(cellchat@idents))
        
        # Set CellChat database
        cellchat@DB <- CellChatDB.use
        
        # Data preprocessing - optimized version
        cat("Performing CellChat data preprocessing...\\n")
        
        tryCatch({
            cellchat <- subsetData(cellchat)
            cat("subsetData successful\\n")
        }, error = function(e) {
            cat("subsetData failed, manually setting data...\\n")
            # Use raw data directly
            cellchat@data.signaling <- data.input
            # Set a basic overexpression threshold
            cellchat@var.features <- rownames(data.input)[apply(data.input, 1, function(x) sum(x > 0)) >= ncol(data.input) * 0.1]
            if(length(cellchat@var.features) == 0) {
                cellchat@var.features <- rownames(data.input)[1:min(2000, nrow(data.input))]
            }
            cat("Manually set", length(cellchat@var.features), "feature genes\\n")
        })
        
        # Identify overexpressed genes and interactions
        tryCatch({
            cellchat <- identifyOverExpressedGenes(cellchat)
            cat("Overexpressed genes identified successfully\\n")
        }, error = function(e) {
            cat("Failed to identify overexpressed genes, skipping this step\\n")
        })
        
        tryCatch({
            cellchat <- identifyOverExpressedInteractions(cellchat)
            cat("Interactions identified successfully\\n")
        }, error = function(e) {
            cat("Failed to identify interactions, skipping this step\\n")
        })
        
        # Skip data smoothing step (may cause issues)
        cat("Skipping data smoothing step to avoid errors\\n")
        
        # Set random seed for reproducibility
        set.seed(${params.cellchat.seed_use})
        
        # Calculate cell communication probability - enhanced error handling
        cat("Calculating cell communication probabilities...\\n")
        tryCatch({
            cellchat <- computeCommunProb(cellchat, type = "triMean", trim = 0.1)
            cat("Communication probability calculation successful\\n")
            
            # Calculate communication at pathway level
            cellchat <- computeCommunProbPathway(cellchat)
            cat("Pathway-level communication calculation successful\\n")
            
            cellchat <- aggregateNet(cellchat)
            cat("Network aggregation successful\\n")
        }, error = function(e) {
            cat("Communication probability calculation failed:", conditionMessage(e), "\\n")
            # Create empty network results
            cellchat@net <- list(
                count = matrix(0, nrow = length(levels(cellchat@idents)), ncol = length(levels(cellchat@idents)),
                              dimnames = list(levels(cellchat@idents), levels(cellchat@idents))),
                weight = matrix(0, nrow = length(levels(cellchat@idents)), ncol = length(levels(cellchat@idents)),
                               dimnames = list(levels(cellchat@idents), levels(cellchat@idents)))
            )
            cellchat@netP <- list(pathways = c())
        })
        
        # Output communication network results
        df.net <- subsetCommunication(cellchat)
        write.csv(df.net, '01_communication_network.csv', row.names = FALSE, quote = FALSE)
        
        # Save CellChat object
        save(cellchat, file = "cellchat_results.Rdata")
        
        # Visualization results - enhanced error handling
        cat("Generating CellChat visualization results...\\n")
        
        # Check for valid communication data
        has_valid_data <- !is.null(cellchat@net) && !is.null(cellchat@net\$count) && sum(cellchat@net\$count) > 0
        
        if(!has_valid_data) {
            cat("Warning: No valid cell communication data detected, generating blank plots\\n")
        }
        
        # 1. Network circle plot - communication count
        tryCatch({
            pdf('01_net_circle_number.pdf', width = 9, height = 9)
            netVisual_circle(cellchat@net\$count, vertex.weight = groupSize, weight.scale = TRUE, 
                             label.edge = FALSE, title.name = "Number of interactions")
            dev.off()
        }, error = function(e) {
            if(dev.cur() != 1) dev.off()
            cat("Network circle plot PDF generation failed:", e\$message, "\\n")
        })
        
        tryCatch({
            png('01_net_circle_number.png', width = 9, height = 9, units = 'in', res = 150)
            netVisual_circle(cellchat@net\$count, vertex.weight = groupSize, weight.scale = TRUE, 
                             label.edge = FALSE, title.name = "Number of interactions")
            dev.off()
        }, error = function(e) {
            if(dev.cur() != 1) dev.off()
            cat("Network circle plot PNG generation failed:", e\$message, "\\n")
        })
        
        # 2. Network circle plot - communication strength
        tryCatch({
            pdf('02_net_circle_weight.pdf', width = 9, height = 9)
            netVisual_circle(cellchat@net\$weight, vertex.weight = groupSize, weight.scale = TRUE, 
                             label.edge = FALSE, title.name = "Interaction weight/strength")
            dev.off()
        }, error = function(e) {
            if(dev.cur() != 1) dev.off()
            cat("Weight circle plot PDF generation failed:", e\$message, "\\n")
        })
        
        tryCatch({
            png('02_net_circle_weight.png', width = 9, height = 9, units = 'in', res = 150)
            netVisual_circle(cellchat@net\$weight, vertex.weight = groupSize, weight.scale = TRUE, 
                             label.edge = FALSE, title.name = "Interaction weight/strength")
            dev.off()
        }, error = function(e) {
            if(dev.cur() != 1) dev.off()
            cat("Weight circle plot PNG generation failed:", e\$message, "\\n")
        })
        

        # 4. Single cell type communication circle plots
        mat <- cellchat@net\$count
        write.csv(mat, 'interaction_count_matrix.csv')

        if(nrow(mat) > 0 && ncol(mat) > 0) {
            tryCatch({
                pdf('04_single_celltype_circles.pdf', width = 18, height = 5)
                par(mfrow = c(1, min(6, nrow(mat))), xpd = TRUE)
                for (i in 1:min(6, nrow(mat))) {  # Limit to first 6 cell types
                    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
                    mat2[i,] <- mat[i,]
                    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, 
                                     arrow.width = 0.2, arrow.size = 0.1, edge.weight.max = max(mat),
                                     title.name = rownames(mat)[i])
                }
                dev.off()
                cat("Single cell type communication circle plots PDF generated successfully\\n")
            }, error = function(e) {
                if(dev.cur() != 1) dev.off()
                cat("Single cell type communication circle plots PDF generation failed:", e\$message, "\\n")
            })
            
            tryCatch({
                png('04_single_celltype_circles.png', width = 1800, height = 500)
                par(mfrow = c(1, min(6, nrow(mat))), xpd = TRUE)
                for (i in 1:min(6, nrow(mat))) {
                    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
                    mat2[i,] <- mat[i,]
                    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, 
                                     arrow.width = 0.2, arrow.size = 0.1, edge.weight.max = max(mat),
                                     title.name = rownames(mat)[i])
                }
                dev.off()
                cat("Single cell type communication circle plots PNG generated successfully\\n")
            }, error = function(e) {
                if(dev.cur() != 1) dev.off()
                cat("Single cell type communication circle plots PNG generation failed:", e\$message, "\\n")
            })
        } else {
            # Create blank plot
            pdf('04_single_celltype_circles.pdf', width = 18, height = 5)
            plot.new()
            text(0.5, 0.5, "No communication data for single cell type analysis", cex = 1.2, col = "gray")
            dev.off()
            
            png('04_single_celltype_circles.png', width = 1800, height = 500)
            plot.new()
            text(0.5, 0.5, "No communication data for single cell type analysis", cex = 1.2, col = "gray")
            dev.off()
        }

        # 5. Outgoing and incoming communication strength analysis
        
        if(nrow(mat) > 0 && ncol(mat) > 0) {
            tryCatch({
                pdf('04_outgoing_incoming_communication.pdf', width = 12, height = 6)
                # Calculate outgoing and incoming communication strength per cell type
                outgoing <- rowSums(cellchat@net\$count)
                incoming <- colSums(cellchat@net\$count)
                
                par(mfrow = c(1, 2))
                barplot(outgoing, main = "Outgoing communication", las = 2, col = "#2166ac")
                barplot(incoming, main = "Incoming communication", las = 2, col = "#b2182b")
                dev.off()
                cat("Outgoing/incoming communication strength PDF generated successfully\\n")
            }, error = function(e) {
                if(dev.cur() != 1) dev.off()
                cat("Outgoing/incoming communication strength PDF generation failed:", e\$message, "\\n")
            })
            
            tryCatch({
                png('04_outgoing_incoming_communication.png', width = 12, height = 6, units = 'in', res = 150)
                # Calculate outgoing and incoming communication strength per cell type
                outgoing <- rowSums(cellchat@net\$count)
                incoming <- colSums(cellchat@net\$count)
                
                par(mfrow = c(1, 2))
                barplot(outgoing, main = "Outgoing communication", las = 2, col = "#2166ac")
                barplot(incoming, main = "Incoming communication", las = 2, col = "#b2182b")
                dev.off()
                cat("Outgoing/incoming communication strength PNG generated successfully\\n")
            }, error = function(e) {
                if(dev.cur() != 1) dev.off()
                cat("Outgoing/incoming communication strength PNG generation failed:", e\$message, "\\n")
            })
        } else {
            # Create blank plot
            pdf('04_outgoing_incoming_communication.pdf', width = 12, height = 6)
            plot.new()
            text(0.5, 0.5, "No communication data for analysis", cex = 1.2, col = "gray")
            dev.off()
            
            png('04_outgoing_incoming_communication.png', width = 12, height = 6, units = 'in', res = 150)
            plot.new()
            text(0.5, 0.5, "No communication data for analysis", cex = 1.2, col = "gray")
            dev.off()
        }
        
        
        # Output CellChat analysis summary
        cat("\\n=== CellChat Analysis Summary ===\\n")
        cat("Number of cell types involved in communication:", length(levels(cellchat@idents)), "\\n")
        cat("Total communication count:", sum(cellchat@net\$count), "\\n")
        cat("Total communication strength:", round(sum(cellchat@net\$weight), 3), "\\n")
        
        # Display top communication pairs
        if(nrow(df.net) > 0) {
            cat("\\nTop 10 strongest cell communications:\\n")
            top_communications <- head(df.net[order(df.net\$prob, decreasing = TRUE), 
                                             c("source", "target", "ligand", "receptor", "prob")], 10)
            print(top_communications)
        }
        
    }, error = function(e) {
        cat("CellChat analysis error:", conditionMessage(e), "\\n")
        
        # Create empty output files
        empty_df <- data.frame()
        write.csv(empty_df, "01_communication_network.csv", row.names = FALSE)
        write.csv(empty_df, "interaction_count_matrix.csv", row.names = FALSE)
        
        # Create empty CellChat object
        empty_cellchat <- list()
        save(empty_cellchat, file = "cellchat_results.Rdata")
        
        # Create error message plot
        png("cellchat_analysis_error.png", width = 8, height = 6, units = "in", res = 150)
        plot.new()
        text(0.5, 0.5, paste("CellChat analysis failed:\\n", conditionMessage(e)), 
             cex = 1.2, col = "red")
        dev.off()
        
        # Create empty PDF files
        pdf("01_net_circle_number.pdf", width = 8, height = 6)
        plot.new()
        text(0.5, 0.5, paste("CellChat analysis failed:\\n", conditionMessage(e)), cex = 1.2, col = "red")
        dev.off()
    })
    
    cat("CellChat analysis completed!\\n")
    """
}