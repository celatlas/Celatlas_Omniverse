#!/usr/bin/env nextflow

/*
 * Spatial Transcriptomics Analysis Pipeline - Nextflow Version
 * Author: Qingsong Li
 * Description: Modularize read_improved.R script into Nextflow workflow
 */

nextflow.enable.dsl = 2
nextflow.enable.dsl = 2

// =============================================
// Celatlas Spatial Advanced Analysis Title
// =============================================
println """
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│    ██████╗ ███████╗██╗      █████╗ ████████╗██╗      █████╗ ███████╗         │
│   ██╔════╝ ██╔════╝██║     ██╔══██╗╚══██╔══╝██║     ██╔══██╗██╔════╝         │
│   ██║      █████╗  ██║     ███████║   ██║   ██║     ███████║███████╗         │
│   ██║      ██╔══╝  ██║     ██╔══██║   ██║   ██║     ██╔══██║╚════██║         │
│   ╚██████╗ ███████╗███████╗██║  ██║   ██║   ███████╗██║  ██║███████║         │
│    ╚═════╝ ╚══════╝╚══════╝╚═╝  ╚═╝   ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝         │
│                                                                              │
│    ██████╗ ███╗   ███╗███╗   ██╗██╗██╗   ██╗███████╗██████╗ ███████╗███████╗ │
│   ██╔═══██╗████╗ ████║████╗  ██║██║██║   ██║██╔════╝██╔══██╗██╔════╝██╔════╝ │
│   ██║   ██║██╔████╔██║██╔██╗ ██║██║██║   ██║█████╗  ██████╔╝███████╗█████╗   │
│   ██║   ██║██║╚██╔╝██║██║╚██╗██║██║╚██╗ ██╔╝██╔══╝  ██╔══██╗╚════██║██╔══╝   │
│   ╚██████╔╝██║ ╚═╝ ██║██║ ╚████║██║ ╚████╔╝ ███████╗██║  ██║███████║███████╗ │
│    ╚═════╝ ╚═╝     ╚═╝╚═╝  ╚═══╝╚═╝  ╚═══╝  ╚══════╝╚═╝  ╚═╝╚══════╝╚══════╝ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
"""

// Parameter definitions
params.help = false
params.input_dir = null
params.output_dir = './results'
params.species = 'mmu'
params.modules = ['qc', 'pca', 'marker', 'singler']

// Optional module switches - consistent with configuration file
params.enable_cellchat = true
params.enable_monocle = true  
params.enable_enrichment = true

println "🚀 Starting spatial transcriptomics advanced analysis pipeline..."
println "📂 Input directory: ${params.input_dir}"
println "💾 Output directory: ${params.output_dir}"
println "📝 Analysis species: ${params.species}"

// Build complete analysis module list
def all_modules = params.modules.clone()
if (params.enable_cellchat) all_modules.add('cellchat')
if (params.enable_monocle) all_modules.add('monocle')
if (params.enable_enrichment) all_modules.add('enrichment')

println "📝 Analysis modules: ${all_modules.join(', ')}"
println "⏳ Starting data processing...\n"
params.enable_report_generation = true
params.enable_ui_optimization = true
params.report_config = null
params.qc = true
params.sample_id = "ST110250_A1"
params.tissue_type = "Brain"
params.resolution = "20 μm"
params.tissue_image_path = null  // Custom tissue slice image path

// QC analysis parameters - match original script
params.max_mt_percent = 5        // Match original script: percent.mt < 5
params.min_features = 2000       // Match original script: nFeature_Spatial > 2000  
params.max_features = 250000     // Match original script: nCount_Spatial < 250000
params.min_counts = 200          // Match original script: nCount_Spatial > 200
params.min_cells = 5

// Spatial plot display parameters
params.spatial_plot_pt_size = 0.2  // Point size factor in spatial plots

// PCA and clustering parameters - match original script
params.clustering = [
    n_pcs: 50,                   
    resolution: 0.2,            
    k_param: 20,
    n_neighbors: 30,
    min_dist: 0.3,
    chose_pcs: 20
]

// PCA visualization parameters
params.pca_elbow_plot_dims = 50    // Number of principal components shown in elbow plot

// Marker gene analysis parameters
params.marker = [
    min_pct: 0.25,
    logfc_threshold: 0.25,
    test_use: "wilcox",
    top_n: 5
]

// SingleR parameters
params.singler = [
    reference: "auto",
    fine_tune: true,
    prune_score: 0.5
]

// CellChat parameters
params.cellchat = [
    database: "auto",
    min_cells: 5,
    seed_use: 1
]

// Enrichment analysis parameters
params.enrichment = [
    organism: "auto",     // Auto-detect species
    pval_cutoff: 0.05,
    logfc_cutoff: 0.5,
    min_pct: 0.1
]

// Monocle parameters
params.monocle = [
    min_cells_per_type: 20,
    num_dim: 30,
    reduction_method: "UMAP"
]

// Import parameters
include { checkParams; helpMessage } from './lib/utils'

// Import workflow modules
include { QC_ANALYSIS } from './modules/qc_analysis'
include { PCA_CLUSTERING } from './modules/pca_clustering'
include { MARKER_ANALYSIS } from './modules/marker_analysis'
include { SINGLER_ANNOTATION } from './modules/singler_annotation'
include { CELLCHAT_ANALYSIS } from './modules/cellchat_analysis'
include { CELLCHAT_PLOTS_POSTPROCESS } from './modules/cellchat_plots_postprocess'
include { MONOCLE_ANALYSIS } from './modules/monocle_analysis'
include { ENRICHMENT_ANALYSIS } from './modules/enrichment_analysis'
include { REPORT_GENERATION; OPTIMIZE_REPORT_UI } from './modules/report_generation'

// Parameter validation
if (params.help) {
    helpMessage()
    exit 0
}

// Check required parameters
checkParams(params)


workflow {
    // Define input channels
    input_ch = Channel.fromPath("${params.input_dir}/filtered_feature_bc_matrix", type: 'dir')
    spatial_ch = Channel.fromPath("${params.input_dir}/spatial", type: 'dir')
    
    // Combine input data
    data_ch = input_ch.combine(spatial_ch)
    
    // 1. Quality control and data loading
    if ('qc' in params.modules) {
        qc_results = QC_ANALYSIS(data_ch)
        seurat_obj = qc_results.seurat_object
        
    } else {
        // If QC is skipped, create basic Seurat object
        seurat_obj = Channel.empty()
    }
    
    // 2. PCA and clustering analysis
    if ('pca' in params.modules) {
        pca_results = PCA_CLUSTERING(seurat_obj)
        clustered_obj = pca_results.seurat_clustered
        
    } else {
        clustered_obj = seurat_obj
    }
    
    // 3. Marker gene analysis
    marker_results = Channel.empty()
    if ('marker' in params.modules) {
        marker_results = MARKER_ANALYSIS(clustered_obj)
       
    }
    
    // 4. SingleR cell type annotation
    annotated_obj = Channel.empty()
    if ('singler' in params.modules) {
        singler_results = SINGLER_ANNOTATION(clustered_obj)
        annotated_obj = singler_results.seurat_annotated
       
    } else {
        annotated_obj = clustered_obj
    }
    
    // 5. CellChat cell communication analysis
    cellchat_results = Channel.empty()
    cellchat_plots_results = Channel.empty()
    if (params.enable_cellchat) {
        cellchat_results = CELLCHAT_ANALYSIS(annotated_obj)
        
        // 5.1 CellChat plot post-processing (generate plots at end of pipeline)
        cellchat_plots_results = CELLCHAT_PLOTS_POSTPROCESS(
            cellchat_results.cellchat_object
        )
    }
    
    // 6. Functional enrichment analysis
    enrichment_results = Channel.empty()
    if (params.enable_enrichment && !marker_results.isEmpty()) {
        enrichment_results = ENRICHMENT_ANALYSIS(annotated_obj, marker_results.markers)
       
    }
    
    // 7. Monocle3 trajectory analysis
    monocle_results = Channel.empty()
    if (params.enable_monocle) {
        monocle_results = MONOCLE_ANALYSIS(annotated_obj)
      
    }
    
    // 8. Generate comprehensive analysis report (if enabled)
    if (params.enable_report_generation) {
        // Collect all completed analysis modules, ensure all enabled modules are complete
        all_completed_signals = Channel.empty()
        
        // Collect completion signals from each module
        if ('qc' in params.modules) {
            all_completed_signals = all_completed_signals.mix(qc_results.seurat_object.map{ 'qc_complete' })
        }
        
        if ('pca' in params.modules) {
            all_completed_signals = all_completed_signals.mix(pca_results.seurat_clustered.map{ 'pca_complete' })
        }
        
        if ('marker' in params.modules && !marker_results.isEmpty()) {
            all_completed_signals = all_completed_signals.mix(marker_results.markers.map{ 'marker_complete' })
        }
        
        if ('singler' in params.modules) {
            all_completed_signals = all_completed_signals.mix(singler_results.seurat_annotated.map{ 'singler_complete' })
        }
        
        if (params.enable_cellchat && !cellchat_results.isEmpty()) {
            // Wait for both CellChat main analysis and plot post-processing to complete
            all_completed_signals = all_completed_signals.mix(cellchat_results.cellchat_object.map{ 'cellchat_complete' })
        }
        
        if (params.enable_enrichment && !enrichment_results.isEmpty()) {
            all_completed_signals = all_completed_signals.mix(enrichment_results.summary.map{ 'enrichment_complete' })
        }
        
        if (params.enable_monocle && !monocle_results.isEmpty()) {
            all_completed_signals = all_completed_signals.mix(monocle_results.monocle_cds.map{ 'monocle_complete' })
        }
        
        // Wait for all signals to be collected
        final_completion_signal = all_completed_signals.collect()
        
        // Generate report
        report_results = REPORT_GENERATION(
            Channel.fromPath("${params.output_dir}"),
            final_completion_signal
        )
        
        // UI optimization (if enabled)
        if (params.enable_ui_optimization) {
            optimized_report = OPTIMIZE_REPORT_UI(
                report_results.report,
                report_results.assets
            )
            
        }
        
    }
}

workflow.onComplete {
    log.info """
    ===========================================
    Pipeline execution completed!
    ===========================================
    Execution status: ${workflow.success ? 'Success' : 'Failed'}
    Work directory: ${workflow.workDir}
    Completion time: ${workflow.complete}
    Runtime duration: ${workflow.duration}
    ===========================================
    
    Result files saved in: ${params.output_dir}
    
    Module result directories:
    - QC results: ${params.output_dir}/1.QC/
    - PCA results: ${params.output_dir}/2.PCA/
    - Marker results: ${params.output_dir}/3.Marker/
    - SingleR results: ${params.output_dir}/4.SingleR/
    - CellChat results: ${params.output_dir}/5.Cellchat/
    - Monocle results: ${params.output_dir}/6.Pseudotime/
    - Enrichment results: ${params.output_dir}/7.Enrichment/
    ${params.enable_report_generation ? "    - Comprehensive report: ${params.output_dir}/reports/" : ""}
    """
    
    if (!workflow.success) {
        log.error "Pipeline execution failed, please check error messages"
    }
}

workflow.onError {
    log.error "Pipeline execution error: ${workflow.errorMessage}"
}