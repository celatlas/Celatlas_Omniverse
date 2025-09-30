/*
 * Common utility functions
 */

// Parameter checking function
def checkParams(params) {
    def errors = []
    
    // Check required input directory
    if (!params.input_dir) {
        errors.add("Please specify input directory --input_dir")
    } else if (!file(params.input_dir).exists()) {
        errors.add("Input directory does not exist: ${params.input_dir}")
    }
    
    // Check input data files
    def matrix_dir = file("${params.input_dir}/filtered_feature_bc_matrix")
    def spatial_dir = file("${params.input_dir}/spatial")
    
    if (!matrix_dir.exists()) {
        errors.add("Expression matrix directory not found: ${matrix_dir}")
    }
    
    if (!spatial_dir.exists()) {
        errors.add("Spatial information directory not found: ${spatial_dir}")
    }
    
    // Check species parameter
    def valid_species = ['auto', 'human', 'mouse', 'mmu']
    if (params.species && !(params.species in valid_species)) {
        errors.add("Invalid species parameter: ${params.species}. Supported options: ${valid_species.join(', ')}")
    }
    
    // Check core module parameters
    def valid_core_modules = ['qc', 'pca', 'marker', 'singler', 'cellchat', 'enrichment', 'monocle']
    if (params.modules) {
        def invalid_modules = params.modules.findAll { !(it in valid_core_modules) }
        if (invalid_modules) {
            errors.add("Invalid core modules: ${invalid_modules.join(', ')}. Supported core modules: ${valid_core_modules.join(', ')}")
        }
    }
    
    // Check optional module switches
    def optional_modules = ['enable_cellchat', 'enable_monocle', 'enable_enrichment']
    optional_modules.each { module ->
        if (params.containsKey(module) && !(params[module] instanceof Boolean)) {
            errors.add("${module} must be a boolean value (true/false)")
        }
    }
    
    // If errors exist, print and exit
    if (errors) {
        log.error "Parameter validation failed:"
        errors.each { log.error "  - ${it}" }
        exit 1
    }
}

// Help message function
def helpMessage() {
    log.info """
    =========================================
    Spatial Transcriptomics Analysis Pipeline v1.0
    =========================================
    
    Usage:
        nextflow run main.nf [options]
    
    Required parameters:
        --input_dir PATH        10X spatial transcriptomics data directory path
        --output_dir PATH       Results output directory (default: ./results)
    
    Optional parameters:
        --species STRING        Species type [auto|human|mouse] (default: auto)
        --modules LIST          Core modules to run (default: qc,pca,marker,singler,cellchat,enrichment,monocle)
                               Core modules: qc,pca,marker,singler,cellchat,enrichment,monocle
    
    Optional module switches:
        --enable_cellchat BOOL  Enable CellChat cell communication analysis (default: false)
        --enable_monocle BOOL   Enable Monocle3 trajectory analysis (default: false)  
        --enable_enrichment BOOL Enable functional enrichment analysis (default: false)
    
    QC parameters:
        --qc.min_features INT   Minimum genes per cell (default: 200)
        --qc.max_features INT   Maximum genes per cell (default: 5000)
        --qc.max_mt_percent INT Maximum mitochondrial gene percentage (default: 20)
        --qc.min_cells INT      Minimum cells per gene (default: 3)
    
    Clustering parameters:
        --clustering.n_pcs INT       Number of PCA components (default: 30)
        --clustering.resolution FLOAT Clustering resolution (default: 0.5)
        --clustering.k_param INT     KNN parameter (default: 20)
    
    Marker gene parameters:
        --marker.min_pct FLOAT          Minimum expression percentage (default: 0.25)
        --marker.logfc_threshold FLOAT  log2FC threshold (default: 0.25)
        --marker.test_use STRING        Statistical test method (default: wilcox)
        --marker.top_n INT              Top genes per cluster (default: 5)
    
    SingleR parameters:
        --singler.reference STRING    Reference dataset [auto|HumanPrimaryCellAtlasData|MouseRNAseqData] (default: auto)
        --singler.fine_tune BOOLEAN   Whether to perform fine tuning (default: true)
        --singler.prune_score FLOAT   Pruning score threshold (default: 0.5)
    
    CellChat parameters:
        --cellchat.database STRING    Database [auto|CellChatDB.human|CellChatDB.mouse] (default: auto)
        --cellchat.min_cells INT      Minimum cells per cell type (default: 10)
        --cellchat.thresh FLOAT       Communication significance threshold (default: 0.05)
    
    Monocle3 parameters:
        --monocle.num_dim INT             Dimensionality reduction dimensions (default: 30)
        --monocle.min_cells_per_group INT Minimum cells per group (default: 20)
        --monocle.reduction_method STRING Dimensionality reduction method (default: UMAP)
    
    Enrichment analysis parameters:
        --enrichment.padj_cutoff FLOAT    Adjusted p-value cutoff (default: 0.05)
        --enrichment.qvalue_cutoff FLOAT  Q-value cutoff (default: 0.2)
        --enrichment.min_genes INT        Minimum number of genes (default: 5)
        --enrichment.max_genes INT        Maximum number of genes (default: 500)
    
    Computational resources:
        --max_cpus INT          Maximum number of CPUs (default: 4)
        --max_memory STRING     Maximum memory (default: 8.GB)
        --max_time STRING       Maximum runtime (default: 24.h)
    
    Example usage:
        # Complete analysis pipeline
        nextflow run main.nf --input_dir /path/to/10x_data --output_dir ./results
        
        # Run only core modules
        nextflow run main.nf --input_dir /path/to/10x_data --modules qc,pca,marker,singler
        
        # Enable optional modules
        nextflow run main.nf --input_dir /path/to/10x_data --enable_cellchat true --enable_monocle true
        
        # Specify species and custom parameters
        nextflow run main.nf --input_dir /path/to/10x_data --species human \\
            --clustering.resolution 0.8 --marker.logfc_threshold 0.5
        
        # Use configuration file
        nextflow run main.nf -c custom.config
    
    For more information see: https://github.com/your-repo/spatial-transcriptomics-pipeline
    """
}

// Create output directories function
def createOutputDirs(output_dir) {
    def dirs = [
        "${output_dir}/1.QC",
        "${output_dir}/2.PCA", 
        "${output_dir}/3.Marker",
        "${output_dir}/4.SingleR",
        "${output_dir}/5.Cellchat",
        "${output_dir}/6.Pseudotime",
        "${output_dir}/7.Enrichment",
        "${output_dir}/reports"
    ]
    
    dirs.each { dir ->
        file(dir).mkdirs()
    }
}