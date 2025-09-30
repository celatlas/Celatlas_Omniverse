#!/usr/bin/env python3
"""
Advanced Spatial Analysis Report Generator

"""

import os
import sys
import json
import base64
import re
from pathlib import Path
from datetime import datetime
import pandas as pd
import traceback
from typing import Dict, List, Optional, Any

try:
    from jinja2 import Template, Environment, FileSystemLoader
    JINJA2_AVAILABLE = True
except ImportError:
    print("Warning: jinja2 is unavailable and will be replaced with an enhanced string")
    JINJA2_AVAILABLE = False


class UltimateReportGenerator:
    
    def __init__(self, input_dir: str, output_dir: str, template_path: str = None, tissue_image_path: str = None):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.tissue_image_path = tissue_image_path
        
        if template_path and Path(template_path).exists():
            self.template_path = Path(template_path)
        else:
            possible_templates = [
                Path(__file__).parent.parent / "advanced_spatial_report_template_enhanced.html",
                Path("advanced_spatial_report_template_enhanced.html")
            ]
            
            self.template_path = None
            for template in possible_templates:
                if template.exists():
                    self.template_path = template
                    break
                    
            if not self.template_path:
                raise FileNotFoundError("The HTML template file was not found")
        
        print(f"Input directory: {self.input_dir}")
        print(f"Output directory: {self.output_dir}")  
        print(f"Template file: {self.template_path}")
        
        self.data = {
            'report_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'software_version': '2.6.0',
            'analysis_id': f'SPATIAL_{datetime.now().strftime("%Y%m%d_%H%M%S")}',
            'sample_id': 'Unknown',
            'tissue_type': 'Unknown', 
            'technology': 'Spatial_ST',
            'resolution': 'Unknown',
            'analysis_date': datetime.now().strftime('%Y-%m-%d'),
            'tissue_image_path': '',
            'images': {}
        }
    
    def encode_image_base64(self, image_path: Path) -> str:
        try:
            file_size = image_path.stat().st_size
            if file_size > 5 * 1024 * 1024:  # 5MB
                print(f"Large image detected ({file_size/1024/1024:.1f}MB), consider optimization: {image_path}")
            
            with open(image_path, 'rb') as f:
                image_data = f.read()
            
            ext = image_path.suffix.lower()
            if ext == '.png':
                mime_type = 'image/png'
            elif ext in ['.jpg', '.jpeg']:
                mime_type = 'image/jpeg'
            else:
                mime_type = 'image/png'
            
            encoded = base64.b64encode(image_data).decode('utf-8')
            return f"data:{mime_type};base64,{encoded}"
        except Exception as e:
            print(f"Failed to encode the image {image_path}: {e}")
            return ""
    
    
    def find_specific_file(self, filename_pattern: str) -> Optional[Path]:

        found_files = list(self.input_dir.rglob(filename_pattern))
        if found_files:
            print(f"Find the file {filename_pattern}: {found_files[0]}")
            return found_files[0]
        else:
            print(f"No file found: {filename_pattern}")
            return None
    
    def find_tissue_hires_image(self) -> Optional[Path]:

        try:

            input_path = Path(self.input_dir).resolve()
            print(f"Current input_dir: {input_path}")
            
            current_path = input_path
            while current_path != current_path.parent:

                bin_segment_path = current_path / "06.binSegment" / "images" / "tissue_hires_image.png"
                if bin_segment_path.exists():
                    print(f"Found tissue_hires_image.png: {bin_segment_path}")
                    return bin_segment_path
                current_path = current_path.parent
            
            parent_path = input_path.parent
            while parent_path != parent_path.parent:
                bin_segment_path = parent_path / "06.binSegment" / "images" / "tissue_hires_image.png"
                if bin_segment_path.exists():
                    print(f"Found tissue_hires_image.png in parent directory: {bin_segment_path}")
                    return bin_segment_path
                parent_path = parent_path.parent
                
            print("tissue_hires_image.png not found")
            return None
            
        except Exception as e:
            print(f"Error searching for tissue_hires_image.png: {e}")
            return None
    
    def collect_all_images(self) -> Dict[str, str]:
        """Collect all required images"""
        images = {}
        
        # QC images
        total_counts = self.find_specific_file("*04_spatial_total_counts.png")
        if total_counts:
            images['qc_total_counts'] = self.encode_image_base64(total_counts)
        
        gene_counts = self.find_specific_file("*05_spatial_gene_counts.png")
        if gene_counts:
            images['qc_gene_counts'] = self.encode_image_base64(gene_counts)
        
        # Clustering images
        umap_clusters = self.find_specific_file("*08_UMAP_clusters.png")
        if umap_clusters:
            images['clustering_umap'] = self.encode_image_base64(umap_clusters)
        
        spatial_clusters = self.find_specific_file("*10_spatial_clusters.png")
        if spatial_clusters:
            images['clustering_spatial'] = self.encode_image_base64(spatial_clusters)
        
        # Marker gene heatmap
        marker_heatmap = self.find_specific_file("*marker_genes_heatmap.png")
        if marker_heatmap:
            images['marker_heatmap'] = self.encode_image_base64(marker_heatmap)
        
        # Annotation images
        singler_heatmap = self.find_specific_file("*10_SingleR_score_heatmap.png")
        if singler_heatmap:
            images['annotation_singler'] = self.encode_image_base64(singler_heatmap)
        
        spatial_celltype = self.find_specific_file("*12_Spatial_CellType.png")
        if spatial_celltype:
            images['annotation_spatial'] = self.encode_image_base64(spatial_celltype)
        
        # Trajectory images - precise mapping
        trajectory_overview = self.find_specific_file("*01_monocle3_trajectory_overview.png")
        if trajectory_overview:
            images['trajectory_graph'] = self.encode_image_base64(trajectory_overview)
        
        pseudotime_ordering = self.find_specific_file("*02_monocle3_pseudotime.png")
        if pseudotime_ordering:
            images['pseudotime_ordering'] = self.encode_image_base64(pseudotime_ordering)
        
        trajectory_genes = self.find_specific_file("*03_trajectory_genes.png")
        if trajectory_genes:
            images['trajectory_genes'] = self.encode_image_base64(trajectory_genes)
        
        # Enrichment analysis images - randomly select a cell type and use all results for that cell type
        import random
        
        # Find all GO bubble plot files
        go_bubble_files = list(self.input_dir.rglob("*GO_bubble_plot.png"))
        kegg_bubble_files = list(self.input_dir.rglob("*KEGG_bubble_plot.png"))
        kegg_network_files = list(self.input_dir.rglob("*KEGG_network.png"))
        
        # Extract all available cell types
        available_celltypes = set()
        for file_path in go_bubble_files + kegg_bubble_files + kegg_network_files:
            # Extract cell type from filename, format is usually celltype_XXX_GO_bubble_plot.png
            filename = file_path.name
            if 'celltype_' in filename:
                parts = filename.split('_')
                if len(parts) >= 2:
                    celltype = parts[1]
                    available_celltypes.add(celltype)
        
        if available_celltypes:
            # Randomly select a cell type
            selected_celltype = random.choice(list(available_celltypes))
            print(f"Selected cell type for enrichment analysis: {selected_celltype}")
            
            # Find corresponding three images for this cell type
            for file_path in go_bubble_files:
                if f'celltype_{selected_celltype}_' in file_path.name:
                    images['enrichment_go_bubble'] = self.encode_image_base64(file_path)
                    print(f"Using GO enrichment plot: {file_path.name}")
                    break
                    
            for file_path in kegg_bubble_files:
                if f'celltype_{selected_celltype}_' in file_path.name:
                    images['enrichment_kegg_bubble'] = self.encode_image_base64(file_path)
                    print(f"Using KEGG enrichment plot: {file_path.name}")
                    break
                    
            for file_path in kegg_network_files:
                if f'celltype_{selected_celltype}_' in file_path.name:
                    images['enrichment_kegg_network'] = self.encode_image_base64(file_path)
                    print(f"Using KEGG network plot: {file_path.name}")
                    break
        else:
            # If cell type pattern not found, fall back to original random selection
            if go_bubble_files:
                selected_go = random.choice(go_bubble_files)
                images['enrichment_go_bubble'] = self.encode_image_base64(selected_go)
                print(f"Randomly selected GO enrichment plot: {selected_go.name}")
            
            if kegg_bubble_files:
                selected_kegg_bubble = random.choice(kegg_bubble_files)
                images['enrichment_kegg_bubble'] = self.encode_image_base64(selected_kegg_bubble)
                print(f"Randomly selected KEGG enrichment plot: {selected_kegg_bubble.name}")
            
            if kegg_network_files:
                selected_kegg_network = random.choice(kegg_network_files)
                images['enrichment_kegg_network'] = self.encode_image_base64(selected_kegg_network)
                print(f"Randomly selected KEGG network plot: {selected_kegg_network.name}")
        
        # CellChat analysis images
        cellchat_network = self.find_specific_file("*01_net_circle_number.png")
        if cellchat_network:
            images['cellchat_network'] = self.encode_image_base64(cellchat_network)
        
        cellchat_heatmap = self.find_specific_file("*03_number_of_interactions_heatmap.png")
        if cellchat_heatmap:
            images['cellchat_heatmap'] = self.encode_image_base64(cellchat_heatmap)
        
        # Logo image - read from bin folder
        logo_path = Path(__file__).parent / "celatlas.png"
        if logo_path.exists():
            images['celatlas_logo'] = self.encode_image_base64(logo_path)
            print(f"Found logo image: {logo_path}")
        else:
            print(f"Logo image not found: {logo_path}")
        
        # Spatial transcriptomics data analysis workflow diagram - read from bin folder
        route_path = Path(__file__).parent / "route.png"
        if route_path.exists():
            images['workflow_diagram'] = self.encode_image_base64(route_path)
            print(f"Found workflow diagram: {route_path}")
        else:
            print(f"Workflow diagram not found: {route_path}")
        
        # Tissue slice image - prioritize specified tissue_image_path
        if self.tissue_image_path and Path(self.tissue_image_path).exists():
            images['tissue_image'] = self.encode_image_base64(Path(self.tissue_image_path))
            print(f"Using specified tissue image: {self.tissue_image_path}")
        else:
            # If not specified or path doesn't exist, search in parent 06.binSegment/images directory
            tissue_image = self.find_tissue_hires_image()
            if tissue_image:
                images['tissue_image'] = self.encode_image_base64(tissue_image)
                print(f"Found tissue slice image: {tissue_image}")
            else:
                print("Tissue slice image not found")
        
        return images
    
    def collect_top3_marker_data(self) -> List[Dict]:
        """Collect top marker genes data for each cluster"""
        markers_file = self.find_specific_file("top*_cluster_markers.csv")
        if not markers_file:
            return []
        
        try:
            df = pd.read_csv(markers_file)
            
            # Group by cluster, take top3 genes for each cluster (or available genes if less than 3)
            cluster_markers = []
            for cluster in sorted(df['cluster'].unique()):
                cluster_df = df[df['cluster'] == cluster].head(3)  # Take top3 for display
                
                top3_genes = []
                for _, row in cluster_df.iterrows():
                    top3_genes.append({
                        'gene': row['gene'],
                        'avg_log2FC': round(row['avg_log2FC'], 3),
                        'pct_1': round(row['pct.1'], 3),
                        'pct_2': round(row['pct.2'], 3),
                        'p_val_adj': f"{row['p_val_adj']:.2e}" if row['p_val_adj'] > 0 else "0"
                    })
                
                cluster_markers.append({
                    'cluster': f"Cluster {cluster}",
                    'cluster_num': cluster,
                    'top3_genes': top3_genes
                })
            
            print(f"Successfully loaded top marker genes data for {len(cluster_markers)} clusters")
            return cluster_markers
            
        except Exception as e:
            print(f"Failed to read marker genes file: {e}")
            return []
    
    def collect_communication_pathways_data(self) -> List[Dict]:
        """Collect CellChat communication pathways data - show one per pathway"""
        communication_file = self.find_specific_file("01_communication_network.csv")
        if not communication_file:
            return []
        
        try:
            df = pd.read_csv(communication_file)
            
            print(f"Total communication network data: {len(df)}")
            print("Communication pathway types:", df['pathway_name'].unique()[:10])  # Show first 10 pathway types
            
            # Group by pathway name, take only the highest probability one for each pathway
            pathway_top = df.loc[df.groupby('pathway_name')['prob'].idxmax()]
            
            # Sort by probability from high to low, take only the top few
            pathway_top = pathway_top.sort_values('prob', ascending=False).head(6)
            
            result = []
            for _, row in pathway_top.iterrows():
                result.append({
                    'pathway_name': row['pathway_name'],
                    'source': row['source'],
                    'target': row['target'],
                    'ligand': row['ligand'],
                    'receptor': row['receptor'],
                    'prob': round(row['prob'], 4)
                })
            
            print(f"Selected top communication pathways count: {len(result)}")
            for item in result:
                print(f"  {item['pathway_name']}: {item['source']} -> {item['target']} (probability: {item['prob']})")
            
            return result
        except Exception as e:
            print(f"Failed to read communication pathways data: {e}")
            return []

    def collect_ultimate_celltype_data(self) -> List[Dict]:
        """Collect ultimate cell type statistics data - precisely calculate percentages and confidence scores"""
        summary_file = self.find_specific_file("cluster_annotation_summary.csv")
        if not summary_file:
            return []
        
        try:
            df = pd.read_csv(summary_file)
            total_cells = df['count'].sum()  # Total count of all cell types
            
            print(f"Total cell count: {total_cells}")
            print("Cell type distribution by cluster:")
            print(df[['seurat_clusters', 'cell_type_filtered', 'count', 'mean_score']])
            
            # Group by cell type, calculate precise statistics
            celltype_stats = df.groupby('cell_type_filtered').agg({
                'count': 'sum',          # Total count for each cell type
                'mean_score': 'mean'     # Average mean_score for each cell type
            }).reset_index()
            
            # Calculate percentage: total count of each cell type / total count of all cell types
            celltype_stats['percentage'] = (celltype_stats['count'] / total_cells * 100).round(2)
            
            # Confidence Score = average mean_score for each cell type
            celltype_stats['confidence_score'] = celltype_stats['mean_score'].round(3)
            
            result = []
            for _, row in celltype_stats.iterrows():
                result.append({
                    'cell_type_filtered': row['cell_type_filtered'],
                    'count': int(row['count']),
                    'percentage': row['percentage'],
                    'confidence_score': row['confidence_score']
                })
            
            print(f"Ultimate Cell Type Summary statistics:")
            for item in result:
                print(f"  {item['cell_type_filtered']}: {item['count']} cells ({item['percentage']}%), confidence score: {item['confidence_score']}")
            
            return result
        except Exception as e:
            print(f"Failed to read cell type statistics: {e}")
            return []
    
    def replace_placeholder_images(self, html_content: str, images: Dict[str, str]) -> str:
        """Replace placeholder images in HTML with actual images"""
        
        placeholder_mappings = [
            # Workflow diagram section
            {
                'pattern': r'<div class="placeholder-image large-diagram">\s*<div>\s*<i class="fas fa-project-diagram fa-4x"[^>]*></i>\s*<div[^>]*>\s*Spatial Transcriptomics Workflow\s*</div>\s*<div[^>]*>\s*完整的时空转录组分析流程图\s*</div>\s*</div>\s*</div>',
                'image_key': 'workflow_diagram',
                'alt_text': 'Spatial transcriptomics data analysis workflow'
            },
            # QC section
            {
                'pattern': r'<div class="placeholder-image">\s*<div>\s*<i class="fas fa-image fa-3x"></i>\s*<div>Total Counts \(UMI\) 空间分布图</div>\s*</div>\s*</div>',
                'image_key': 'qc_total_counts',
                'alt_text': 'Total Counts spatial distribution'
            },
            {
                'pattern': r'<div class="placeholder-image">\s*<div>\s*<i class="fas fa-image fa-3x"></i>\s*<div>n genes by counts 空间分布图</div>\s*</div>\s*</div>',
                'image_key': 'qc_gene_counts',
                'alt_text': 'Gene Counts spatial distribution'
            },
            # Clustering section
            {
                'pattern': r'<div class="placeholder-image">\s*<div>\s*<i class="fas fa-image fa-3x"></i>\s*<div>UMAP细胞聚类图</div>\s*</div>\s*</div>',
                'image_key': 'clustering_umap',
                'alt_text': 'UMAP clustering visualization'
            },
            {
                'pattern': r'<div class="placeholder-image">\s*<div>\s*<i class="fas fa-image fa-3x"></i>\s*<div>空间聚类分布图</div>\s*</div>\s*</div>',
                'image_key': 'clustering_spatial',
                'alt_text': 'Spatial clustering distribution'
            },
            # Marker gene heatmap
            {
                'pattern': r'<div class="chart-placeholder">\s*<div>\s*<i class="fas fa-fire fa-3x"></i>\s*<div>Marker Gene Clustering Heatmap</div>\s*</div>\s*</div>',
                'image_key': 'marker_heatmap',
                'alt_text': 'Marker gene clustering heatmap'
            },
            # Annotation section
            {
                'pattern': r'<div class="placeholder-image">\s*<div>\s*<i class="fas fa-image fa-3x"></i>\s*<div>SingleR Score Heatmap</div>\s*</div>\s*</div>',
                'image_key': 'annotation_singler',
                'alt_text': 'SingleR score heatmap'
            },
            {
                'pattern': r'<div class="placeholder-image">\s*<div>\s*<i class="fas fa-image fa-3x"></i>\s*<div>Spatial Cell Type Plot</div>\s*</div>\s*</div>',
                'image_key': 'annotation_spatial',
                'alt_text': 'Spatial cell type distribution'
            },
            # Trajectory section - precise mapping
            {
                'pattern': r'<div class="placeholder-image">\s*<div>\s*<i class="fas fa-image fa-3x"></i>\s*<div>Monocle3 Trajectory Plot</div>\s*</div>\s*</div>',
                'image_key': 'trajectory_graph',
                'alt_text': 'Trajectory plot'
            },
            {
                'pattern': r'<div class="placeholder-image">\s*<div>\s*<i class="fas fa-image fa-3x"></i>\s*<div>Pseudotime Heatmap</div>\s*</div>\s*</div>',
                'image_key': 'pseudotime_ordering',
                'alt_text': 'Pseudotime ordering'
            },
            {
                'pattern': r'<div class="chart-placeholder">\s*<div>\s*<i class="fas fa-chart-line fa-3x"></i>\s*<div>Gene Expression Trajectory Plot</div>\s*</div>\s*</div>',
                'image_key': 'trajectory_genes',
                'alt_text': 'Trajectory gene expression'
            },
            # Enrichment analysis section
            {
                'pattern': r'<div class="placeholder-image">\s*<div>\s*<i class="fas fa-image fa-3x"></i>\s*<div>GO Bubble Plot</div>\s*</div>\s*</div>',
                'image_key': 'enrichment_go_bubble',
                'alt_text': 'GO enrichment bubble plot'
            },
            {
                'pattern': r'<div class="placeholder-image">\s*<div>\s*<i class="fas fa-image fa-3x"></i>\s*<div>KEGG Bubble Plot</div>\s*</div>\s*</div>',
                'image_key': 'enrichment_kegg_bubble',
                'alt_text': 'KEGG enrichment bubble plot'
            },
            {
                'pattern': r'<div class="chart-placeholder">\s*<div>\s*<i class="fas fa-project-diagram fa-3x"></i>\s*<div>KEGG Network Visualization</div>\s*</div>\s*</div>',
                'image_key': 'enrichment_kegg_network',
                'alt_text': 'KEGG network visualization'
            },
            # CellChat section
            {
                'pattern': r'<div class="placeholder-image">\s*<div>\s*<i class="fas fa-image fa-3x"></i>\s*<div>Overall Communication Network</div>\s*</div>\s*</div>',
                'image_key': 'cellchat_network',
                'alt_text': 'Overall cell communication network'
            },
            {
                'pattern': r'<div class="placeholder-image">\s*<div>\s*<i class="fas fa-image fa-3x"></i>\s*<div>Interactions Heatmap</div>\s*</div>\s*</div>',
                'image_key': 'cellchat_heatmap',
                'alt_text': 'Cell interaction heatmap'
            }
        ]
        
        # Execute replacement
        for mapping in placeholder_mappings:
            if mapping['image_key'] in images and images[mapping['image_key']]:
                # Check if clickable-image class should be added (all images after 4.1)
                is_interactive = mapping['image_key'] in ['qc_total_counts', 'qc_gene_counts', 'clustering_umap', 'clustering_spatial', 'marker_heatmap', 
                                                         'annotation_singler', 'annotation_spatial', 'trajectory_graph', 
                                                         'pseudotime_ordering', 'trajectory_genes', 'enrichment_go_bubble',
                                                         'enrichment_kegg_bubble', 'enrichment_kegg_network', 
                                                         'cellchat_network', 'cellchat_heatmap']
                
                clickable_class = ' clickable-image' if is_interactive else ''
                data_title = f' data-title="{mapping["alt_text"]}"' if is_interactive else ''
                
                img_tag = f'<img src="{images[mapping["image_key"]]}" alt="{mapping["alt_text"]}" class="analysis-image{clickable_class}"{data_title} style="max-width: 100%; height: auto; border-radius: 8px; box-shadow: 0 4px 12px rgba(0,0,0,0.1);">'
                html_content = re.sub(mapping['pattern'], img_tag, html_content, flags=re.MULTILINE | re.DOTALL)
                print(f"Replacement successful: {mapping['alt_text']}" + (" (interactive)" if is_interactive else ""))
            else:
                print(f"Skip replacement (no image data): {mapping['alt_text']}")
        
        return html_content
    
    def replace_logo(self, html_content: str, images: Dict[str, str]) -> str:
        """Replace logo image in HTML"""
        if 'celatlas_logo' in images and images['celatlas_logo']:
            # Replace logo src attribute
            logo_pattern = r'<img src="img/Celatlas_LOGO\.png"([^>]*?)>'
            replacement = f'<img src="{images["celatlas_logo"]}"\\1>'
            html_content = re.sub(logo_pattern, replacement, html_content)
            print("✅ Successfully replaced Celatlas logo image!")
        else:
            print("❌ Logo image data not found")
        return html_content
    
    def insert_top3_marker_table(self, html_content: str, marker_data: List[Dict]) -> str:
        """Replace top marker genes table for each cluster (displays top 3 for each)"""
        if not marker_data:
            return html_content
        
        # Build replacement table body content
        tbody_html = ""
        
        # Add data rows for each cluster, showing top3 genes for all clusters
        for cluster_info in marker_data:  # Show all clusters
            top3_genes = cluster_info['top3_genes']
            
            for gene in top3_genes:
                tbody_html += f'''
                                <tr>
                                    <td>{cluster_info['cluster_num']}</td>
                                    <td>{gene['gene']}</td>
                                    <td>85.3</td>
                                    <td>{gene['p_val_adj']}</td>
                                    <td>{gene['avg_log2FC']}</td>
                                    <td>{gene['p_val_adj']}</td>
                                </tr>'''
        
        # Directly replace tbody content in template
        # Find existing table tbody section and replace
        tbody_pattern = r'(<table class="data-table">\s*<thead>.*?</thead>\s*<tbody>)(.*?)(</tbody>\s*</table>)'
        
        def replacement(match):
            return match.group(1) + tbody_html + match.group(3)
        
        result = re.sub(tbody_pattern, replacement, html_content, flags=re.MULTILINE | re.DOTALL)
        
        if result != html_content:
            print("✅ Successfully replaced top marker genes table data")
        else:
            print("⚠️ Warning: Marker genes table location not found")
        
        return result
    
    def replace_ultimate_celltype_table(self, html_content: str, celltype_data: List[Dict]) -> str:
        """Ultimate replacement of Cell Type Summary table - completely replace static table"""
        if not celltype_data:
            return html_content
        
        # Build new Cell Type Summary table
        new_table_html = '''                <div class="section-subtitle">Cell Type Summary</div>
                <!--
                <div style="margin-bottom: 15px; padding: 10px; background: #e6fffa; border-left: 4px solid #38b2ac; border-radius: 4px;">
                    <p style="margin: 0; font-size: 0.95rem; color: #2d3748; line-height: 1.6;">
                        Based on cluster annotation results, cell types were identified and quantified. 
                        Percentage calculated as: <strong>cell type count ÷ total cell count × 100%</strong>. 
                        Confidence score represents the average mean_score for each cell type.
                    </p>
                </div>
                -->
                <table class="data-table">
                    <thead>
                        <tr style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white;">
                            <th style="padding: 15px 12px; text-align: left; font-family: Arial, sans-serif; line-height: 1.8;">Cell Type</th>
                            <th style="padding: 15px 12px; text-align: left; font-family: Arial, sans-serif; line-height: 1.8;">Cell Count</th>
                            <th style="padding: 15px 12px; text-align: left; font-family: Arial, sans-serif; line-height: 1.8;">Percentage</th>
                            <th style="padding: 15px 12px; text-align: left; font-family: Arial, sans-serif; line-height: 1.8;">Confidence</th>
                        </tr>
                    </thead>
                    <tbody>'''
        
        # Sort by percentage from high to low
        celltype_data_sorted = sorted(celltype_data, key=lambda x: x['percentage'], reverse=True)
        
        for i, celltype in enumerate(celltype_data_sorted):
            # Set row background color, consistent with Top3 genes table style
            row_class = "even" if i % 2 == 0 else "odd"
            
            new_table_html += f'''
                        <tr class="{row_class}">
                            <td style="color: #2d3748; text-align: left; padding: 12px 15px; font-family: Arial, sans-serif; line-height: 1.8;">
                                {celltype['cell_type_filtered']}
                            </td>
                            <td style="color: #2d3748; text-align: left; padding: 12px 15px; font-family: Arial, sans-serif; line-height: 1.8;">
                                {celltype['count']:,}
                            </td>
                            <td style="color: #2d3748; text-align: left; padding: 12px 15px; font-family: Arial, sans-serif; line-height: 1.8;">
                                {celltype['percentage']:.1f}%
                            </td>
                            <td style="color: #2d3748; text-align: left; padding: 12px 15px; font-family: Arial, sans-serif; line-height: 1.8;">
                                {celltype['confidence_score']:.3f}
                            </td>
                        </tr>'''
        
        new_table_html += '''
                    </tbody>
                </table>'''
        
        # Precisely replace entire Cell Type Summary section
        # Match complete section from section-subtitle to </table>
        pattern = r'(<div class="section-subtitle">Cell Type Summary</div>.*?<table class="data-table">.*?</table>)'
        
        def replacement_func(match):
            return new_table_html
        
        result = re.sub(pattern, replacement_func, html_content, flags=re.MULTILINE | re.DOTALL)
        
        if result != html_content:
            print("✅ Successfully replaced Cell Type Summary table!")
        else:
            print("❌ Could not find Cell Type Summary table for replacement")
        
        return result
    
    def replace_communication_pathways_table(self, html_content: str, communication_data: List[Dict]) -> str:
        """Replace Top Communication Pathways table with actual data"""
        if not communication_data:
            return html_content
        
        # Build new communication pathways table
        new_table_html = '''                <div class="section-subtitle">Top Communication Pathways</div>
                <table class="data-table">
                    <thead>
                        <tr style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white;">
                            <th style="padding: 15px 12px; text-align: left; font-family: Arial, sans-serif; line-height: 1.8;">Pathway</th>
                            <th style="padding: 15px 12px; text-align: left; font-family: Arial, sans-serif; line-height: 1.8;">Source</th>
                            <th style="padding: 15px 12px; text-align: left; font-family: Arial, sans-serif; line-height: 1.8;">Target</th>
                            <th style="padding: 15px 12px; text-align: left; font-family: Arial, sans-serif; line-height: 1.8;">Ligand</th>
                            <th style="padding: 15px 12px; text-align: left; font-family: Arial, sans-serif; line-height: 1.8;">Receptor</th>
                            <th style="padding: 15px 12px; text-align: left; font-family: Arial, sans-serif; line-height: 1.8;">Communication Probability</th>
                        </tr>
                    </thead>
                    <tbody>'''
        
        # Add actual data rows
        for i, pathway in enumerate(communication_data):
            row_class = "even" if i % 2 == 0 else "odd"
            new_table_html += f'''
                        <tr class="{row_class}">
                            <td style="color: #2d3748; text-align: left; padding: 12px 15px; font-family: Arial, sans-serif; line-height: 1.8; font-weight: 600;">
                                {pathway['pathway_name']}
                            </td>
                            <td style="color: #2d3748; text-align: left; padding: 12px 15px; font-family: Arial, sans-serif; line-height: 1.8;">
                                {pathway['source']}
                            </td>
                            <td style="color: #2d3748; text-align: left; padding: 12px 15px; font-family: Arial, sans-serif; line-height: 1.8;">
                                {pathway['target']}
                            </td>
                            <td style="color: #2d3748; text-align: left; padding: 12px 15px; font-family: Arial, sans-serif; line-height: 1.8;">
                                {pathway['ligand']}
                            </td>
                            <td style="color: #2d3748; text-align: left; padding: 12px 15px; font-family: Arial, sans-serif; line-height: 1.8;">
                                {pathway['receptor']}
                            </td>
                            <td style="color: #2d3748; text-align: left; padding: 12px 15px; font-family: Arial, sans-serif; line-height: 1.8;">
                                {pathway['prob']}
                            </td>
                        </tr>'''
        
        new_table_html += '''
                    </tbody>
                </table>'''
        
        # Precisely replace entire Top Communication Pathways section
        pattern = r'(<div class="section-subtitle">Top Communication Pathways</div>\s*<table class="data-table">.*?</table>)'
        
        def replacement_func(match):
            return new_table_html
        
        result = re.sub(pattern, replacement_func, html_content, flags=re.MULTILINE | re.DOTALL)
        
        if result != html_content:
            print("✅ Successfully replaced Top Communication Pathways table!")
        else:
            print("❌ Could not find Top Communication Pathways table for replacement")
        
        return result
    
    def generate_report(self, sample_id: str = "Unknown", 
                       tissue_type: str = "Unknown", 
                       resolution: str = "Unknown") -> Path:

        try:
            print("Starting report generation...")

            self.data.update({
                'sample_id': sample_id,
                'tissue_type': tissue_type,
                'resolution': resolution
            })

            images = self.collect_all_images()
            self.data['images'] = images

            marker_data = self.collect_top3_marker_data()
            celltype_data = self.collect_ultimate_celltype_data()
            communication_data = self.collect_communication_pathways_data()

            with open(self.template_path, 'r', encoding='utf-8') as f:
                html_content = f.read()

            basic_replacements = {
                '{{ report_date }}': self.data['report_date'],
                '{{ software_version }}': self.data['software_version'],
                '{{ analysis_id }}': self.data['analysis_id'],
                '{{ sample_id }}': sample_id,
                '{{ tissue_type }}': tissue_type,
                '{{ technology }}': self.data['technology'],
                '{{ resolution }}': resolution,
                '{{ analysis_date }}': self.data['analysis_date'],
                "{{ tissue_hires_image_path | default('img/tissue_sample.png') }}": images.get('tissue_image', ''),
                "{{ tissue_image_base64 | default('') }}": images.get('tissue_image', '')
            }
            
            for placeholder, value in basic_replacements.items():
                html_content = html_content.replace(placeholder, str(value))

            html_content = self.replace_placeholder_images(html_content, images)

            html_content = self.replace_logo(html_content, images)

            html_content = self.insert_top3_marker_table(html_content, marker_data)

            html_content = self.replace_ultimate_celltype_table(html_content, celltype_data)

            html_content = self.replace_communication_pathways_table(html_content, communication_data)

            output_filename = f"Celatlas_Omniverse_report_{sample_id}.html"
            output_path = self.output_dir / output_filename

            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(html_content)

            metadata = {
                'sample_id': sample_id,
                'tissue_type': tissue_type,
                'resolution': resolution,
                'report_date': self.data['report_date'],
                'version': 'ultimate_2.6.0',
                'images_embedded': {k: bool(v) for k, v in images.items()},
                'ultimate_features': [
                    'celltype_summary_table_completely_replaced',
                    'precise_celltype_percentage_calculation',
                    'confidence_score_averaging',
                    'trajectory_analysis_precise_mapping'
                ],
                'data_tables': {
                    'marker_clusters': len(marker_data),
                    'total_display_genes': sum(len(cluster['top3_genes']) for cluster in marker_data),
                    'ultimate_celltype_entries': len(celltype_data)
                },
                'celltype_statistics': {
                    cell['cell_type_filtered']: {
                        'count': cell['count'],
                        'percentage': cell['percentage'],
                        'confidence': cell['confidence_score']
                    } for cell in celltype_data
                }
            }
            
            metadata_path = output_path.with_suffix('.json')
            with open(metadata_path, 'w', encoding='utf-8') as f:
                json.dump(metadata, f, indent=2, ensure_ascii=False)
            
            print(f"The report was generated successfully: {output_path}")
            return output_path
            
        except Exception as e:
            print(f"Report generation failed: {e}")
            traceback.print_exc()
            raise


def main():
    """Main function"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Ultimate Enhanced Spatial Analysis Report Generator')
    parser.add_argument('--input-dir', required=True, help='Input the analysis results into the directory')
    parser.add_argument('--output-dir', required=True, help='Report Output Directory')
    parser.add_argument('--template', help='HTML template file path')
    parser.add_argument('--tissue-image-path', help='Path of slice tissue image')
    parser.add_argument('--sample-id', default='Unknown', help='Sample ID')
    parser.add_argument('--tissue-type', default='Unknown', help='Organizational type')
    parser.add_argument('--resolution', default='Unknown', help='Resolution')
    
    args = parser.parse_args()
    
    generator = UltimateReportGenerator(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        template_path=args.template,
        tissue_image_path=args.tissue_image_path
    )
    
    report_path = generator.generate_report(
        sample_id=args.sample_id,
        tissue_type=args.tissue_type,
        resolution=args.resolution
    )
    
    print(f"🏆 Report path: {report_path}")

AdvancedReportGenerator = UltimateReportGenerator
AdvancedSpatialReportGenerator = UltimateReportGenerator

if __name__ == '__main__':
    main()