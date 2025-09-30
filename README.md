# Celatlas Omniverse - Spatial Transcriptomics Analysis Pipeline

<div align="center">


**Advanced Spatial Transcriptomics Analysis Pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow-23.04.1-brightgreen.svg)](https://www.nextflow.io/)
[![R](https://img.shields.io/badge/R-4.5.1-blue.svg)](https://www.r-project.org/)
[![Java](https://img.shields.io/badge/Java-11-orange.svg)](https://openjdk.java.net/)

</div>

## Overview

Celatlas Omniverse is a comprehensive, modular Nextflow pipeline for spatial transcriptomics data analysis. It provides end-to-end analysis from raw 10X Visium data to publication-ready reports with interactive visualizations.

### Key Features

- **Complete Workflow**: Quality control → Clustering → Annotation → Communication analysis
- **Advanced Analytics**: Cell-cell communication, trajectory analysis, functional enrichment
- **Interactive Reports**: Professional HTML reports with embedded visualizations
- **High Performance**: Parallelized processing with resource optimization
- **Modular Design**: Enable/disable specific analysis modules as needed
- **Containerized**: Docker/Singularity support for reproducibility

## Quick Start

### Environment Installation

> **Important Note**: This pipeline is primarily R-based and requires R 4.5.1+, Nextflow 23.04.1+, and specialized bioinformatics packages for spatial transcriptomics analysis.

#### Step 1: Install R (4.5.1 recommended)

Install R 4.5.1 following the instructions for your operating system at https://cran.r-project.org/

#### Step 2: Install Java (for Nextflow)

**Using SDKMAN (Recommended - Cross-platform):**
```bash
# Install SDKMAN
curl -s https://get.sdkman.io | bash

# Open a new terminal or source the environment
source ~/.sdkman/bin/sdkman-init.sh

# Install Java 11
sdk install java 11.0.10-tem

# Confirm Java installation
java -version
```

#### Step 3: Install Nextflow

```bash
# Download and install Nextflow 23.04.1 specifically
export NXF_VER=23.04.1
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

# Verify installation
nextflow -version
```

#### Step 4: Install R Packages

```bash
# Run the R package installation script
Rscript install_no_sudo.R
```

#### What Gets Installed

**Base Environment Requirements:**
- **R 4.5.1** - Statistical analysis language
- **Java 11** - Nextflow runtime environment  
- **Nextflow 23.04.1** - Workflow engine

**R Packages:**
```
Seurat 5.2.1           # Single-cell analysis
ggplot2 3.5.2          # Data visualization
patchwork 1.3.2        # Plot composition
dplyr 1.1.4            # Data manipulation
SingleR 2.10.0         # Cell type annotation
celldex 1.18.0         # Reference databases
BiocParallel 1.42.1    # Parallel computing
CellChat 2.2.0         # Cell communication analysis
tidyverse 2.0.0        # Data science toolkit
pheatmap 1.0.13        # Heatmap plotting
monocle3 1.4.26        # Trajectory analysis
SeuratWrappers 0.4.0   # Seurat extensions
scrapper               # Additional utilities
```

#### Supported Systems
- **Ubuntu/Debian** (recommended)
- **CentOS/RHEL/Rocky Linux**
- **macOS** (requires Homebrew)

#### Installation Time
- First-time installation: 30-60 minutes
- Network speed affects download time
- Package compilation requires additional time

#### Verify Installation
After installation completes, verify with:
```bash
R --version
nextflow -version
```

Check R packages in R console:
```r
packageVersion("Seurat")
```

### Quick Run
```bash
# Run with your data
nextflow run main.nf \
    --input_dir /path/to/your/data \
    --output_dir ./results \
    --species human\
    --tissue_image_path /path/to/your/tissue_image.png
```

## Analysis Modules

| Module | Description | Output |
|--------|-------------|---------|
| **QC Analysis** | Quality control, filtering, normalization | QC metrics, filtered data |
| **PCA & Clustering** | Dimensionality reduction, cell clustering | Clusters, UMAP plots |
| **Marker Genes** | Differential expression analysis | Top marker genes per cluster |
| **Cell Annotation** | Automated cell type identification | Cell type assignments |
| **CellChat** | Cell-cell communication analysis | Communication networks |
| **Monocle3** | Pseudotime trajectory analysis | Developmental trajectories |
| **Enrichment** | GO/KEGG pathway analysis | Functional annotations |
| **Report Generation** | Interactive HTML reports | Publication-ready reports |

## Input Data Format

Your data should follow the standard 10X Visium structure:
```
input_directory/
├── filtered_feature_bc_matrix/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── spatial/
    ├── tissue_positions_list.csv
    ├── scalefactors_json.json
    └── tissue_hires_image.png
```

## Configuration

### Basic Usage
```bash
# Minimal run
nextflow run main.nf --input_dir /path/to/data

# With all modules enabled
nextflow run main.nf \
    --input_dir /path/to/data \
    --output_dir ./results \
    --species human \
    --tissue_image_path /path/to/your/tissue_image.png
    --enable_cellchat true \
    --enable_monocle true \
    --enable_enrichment true
```

### Configuration File
Create a custom config file:
```nextflow
params {
    input_dir = "/path/to/your/data"
    output_dir = "./results"
    species = "human"  // "human", "mouse", "mmu"
    tissue_image_path = "/path/to/your/tissue_image.png"

    // Sample information
    sample_id = "Sample_01"
    tissue_type = "Brain"
    resolution = "50 μm"
    
    // Module switches
    enable_cellchat = true
    enable_monocle = true
    enable_enrichment = true
    
    // Resources
    max_cpus = 16
    max_memory = '32.GB'
    max_time = '48.h'
}
```

Run with config:
```bash
nextflow run main.nf -c your_config.config
```

## Output Structure

```
results/
├── 1.QC/                    # Quality control results
│   ├── seurat_qc.rds
│   ├── qc_metrics.csv
│   └── *.png
├── 2.PCA/                   # PCA and clustering
│   ├── seurat_clustered.rds
│   ├── cluster_summary.csv
│   └── *.png
├── 3.Marker/                # Marker gene analysis
│   ├── top5_cluster_markers.csv
│   └── *.png
├── 4.SingleR/               # Cell type annotation
│   ├── cluster_annotation_summary.csv
│   └── *.png
├── 5.Cellchat/              # Cell communication (optional)
│   ├── communication_network.csv
│   └── *.png
├── 6.Pseudotime/            # Trajectory analysis (optional)
│   ├── trajectory_results.csv
│   └── *.png
├── 7.Enrichment/            # Functional enrichment (optional)
│   ├── enrichment_results.csv
│   └── *.png
└── reports/                 # HTML reports
    └── Celatlas_Omniverse_report_Sample.html
```

## Advanced Usage

### Resource Configuration
```bash
# High-memory analysis
nextflow run main.nf \
    --input_dir /path/to/data \
    --max_memory 64.GB \
    --max_cpus 32

# Time-limited run
nextflow run main.nf \
    --input_dir /path/to/data \
    --max_time 24.h
```

### Module Selection
```bash
# Run only core modules
nextflow run main.nf \
    --input_dir /path/to/data \
    --modules qc,pca,marker,singler

# Skip specific modules
nextflow run main.nf \
    --input_dir /path/to/data \
    --enable_cellchat false
```

### Resume Failed Runs
```bash
# Resume from last checkpoint
nextflow run main.nf -c config.config -resume
```

## Troubleshooting

### Common Issues

1. **Out of Memory**
   ```bash
   # Increase memory allocation
   --max_memory 64.GB
   ```

2. **Package Installation Errors**
   ```bash
   # Re-run installation script
   ./install_dependencies.sh
   ```

3. **Input Data Issues**
   - Verify data structure
   - Check file permissions
   - Ensure files are not corrupted

### Getting Help
```bash
# View help message
nextflow run main.nf --help

# Check pipeline version
nextflow run main.nf --version
```

## Documentation

- **[Installation Guide](INSTALL.md)** - Detailed installation instructions
- **[User Manual](docs/USER_MANUAL.md)** - Comprehensive usage guide
- **[API Reference](docs/API.md)** - Module documentation
- **[FAQ](docs/FAQ.md)** - Frequently asked questions

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

### Development Setup
```bash
git clone https://github.com/your-repo/celatlas_omniverse.git
cd celatlas_omniverse
git checkout develop
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

- **Issues**: [GitHub Issues](https://github.com/celatlas/Celatlas_Omniverse/issues)
- **Discussions**: [GitHub Discussions](https://github.com/celatlas/Celatlas_Omniverse/discussions)
- **Email**: <rd@celatlas.com>

---

<div align="center">
**Happy analyzing!**

[![Follow](https://img.shields.io/github/followers/celatlas?style=social)](https://github.com/celatlas)
[![Star](https://img.shields.io/github/stars/celatlas/Celatlas_Omniverse?style=social)](https://github.com/celatlas/Celatlas_Omniverse)

</div>
