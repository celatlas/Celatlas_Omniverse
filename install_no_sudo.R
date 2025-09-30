#!/usr/bin/env Rscript

# R Package Installation Script - No Sudo Required Version
# Only installs R packages, assumes R and Java are already installed

cat("=== R Package Installation (No Sudo Required) ===\n")
cat("This script only installs R packages and assumes R and Java are already installed.\n\n")

# Check R version
r_version <- R.version.string
cat("Current R version:", r_version, "\n")

if (getRversion() < "4.0.0") {
    stop("R version 4.0.0 or higher is required")
}

# Set mirror source
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("Installing BiocManager...\n")
    install.packages("BiocManager")
}

# Install remotes
if (!requireNamespace("remotes", quietly = TRUE)) {
    cat("Installing remotes...\n")
    install.packages("remotes")
}

# CRAN package list
cran_packages <- c("ggplot2", "patchwork", "dplyr", "tidyverse", "pheatmap", "Seurat")

# Install CRAN packages
cat("Installing CRAN packages...\n")
for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("Installing", pkg, "...\n")
        install.packages(pkg, dependencies = TRUE)
    } else {
        cat(pkg, "already installed\n")
    }
}

# Install Bioconductor packages
cat("Installing Bioconductor packages...\n")
BiocManager::install(c("SingleR", "celldex"), force = TRUE)
BiocManager::install("BiocParallel")
BiocManager::install("monocle3")
BiocManager::install("scrapper")

# Install GitHub packages
cat("Installing GitHub packages...\n")
if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
    cat("Installing SeuratWrappers from GitHub...\n")
    remotes::install_github("satijalab/seurat-wrappers")
}

if (!requireNamespace("CellChat", quietly = TRUE)) {
    cat("Installing CellChat from GitHub...\n")
    remotes::install_github("jinworks/CellChat")
}

# Verify installation
cat("\nVerifying installation...\n")
required_packages <- c("Seurat", "ggplot2", "patchwork", "dplyr", "SingleR",
                      "celldex", "BiocParallel", "CellChat", "tidyverse", 
                      "pheatmap", "monocle3", "SeuratWrappers", "scrapper")

success_count <- 0
for (pkg in required_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat("✓", pkg, packageVersion(pkg), "\n")
        success_count <- success_count + 1
    } else {
        cat("✗", pkg, "installation failed\n")
    }
}

cat("\nInstallation summary:", success_count, "/", length(required_packages), "packages installed\n")

if (success_count == length(required_packages)) {
    cat("🎉 All R packages installed successfully!\n")
    cat("Note: Make sure R, Java, and Nextflow are installed separately.\n")
} else {
    cat("⚠️ Some packages failed to install. Check error messages above.\n")
}