#!/usr/bin/env nextflow

/*
 * Report Generation Module
 * Integrate all analysis results and generate comprehensive analysis report
 */

process REPORT_GENERATION {
    tag "Generating comprehensive spatial analysis report"
    publishDir "${params.output_dir}/reports", mode: 'copy'

    input:
    path output_plots_dir
    val analysis_results

    output:
    path "Celatlas_Omniverse_report*.html", emit: report
    path "Celatlas_Omniverse_report*.json", emit: metadata, optional: true
    path "report_assets/*", emit: assets, optional: true

    script:
    """
    #!/usr/bin/env python3
    
    import os
    import sys
    from pathlib import Path
    from datetime import datetime
    
    # Add Python script path
    bin_dir = Path('${workflow.projectDir}') / 'bin'
    sys.path.insert(0, str(bin_dir))
    
    # Use new advanced report generator
    try:
        from report_generator import AdvancedReportGenerator
        print("Successfully imported enhanced report generator")
    except ImportError as e:
        print(f"Failed to import enhanced report generator, trying basic version: {e}")
        try:
            from report_generator import AdvancedSpatialReportGenerator as AdvancedReportGenerator
            print("Successfully imported basic report generator")
        except ImportError as e2:
            print(f"Failed to import report generator: {e2}")
            # Create minimal HTML report
            with open('Celatlas_Omniverse_report_error.html', 'w', encoding='utf-8') as f:
                f.write('''
                <!DOCTYPE html>
                <html>
                <head><title>Analysis Report Generation Error</title></head>
                <body>
                    <h1>Report Generation Failed</h1>
                    <p>Unable to load report generation module, please check Python environment configuration.</p>
                    <p>Error message: ''' + str(e2) + '''</p>
                </body>
                </html>
                ''')
            exit(1)
    
    # Print analysis completion information
    completion_signals = "${analysis_results}".split()
    print(f"Received completion signals: {completion_signals}")
    print(f"Number of completed modules: {len(completion_signals)}")
    print("All analysis modules completed, starting report generation...")
    
    # Wait one second to ensure all files are written
    import time
    time.sleep(1)
    
    # Create advanced report generator instance  
    print(f"Input directory: ${output_plots_dir}")
    print(f"Sample ID: ${params.sample_id ?: 'Unknown'}")
    print(f"Tissue type: ${params.tissue_type ?: 'Unknown'}")
    print(f"Resolution: ${params.resolution ?: 'Unknown'}")
    
    # Use new report generator
    tissue_image_arg = '${params.tissue_image_path}' if '${params.tissue_image_path}' != 'null' else None
    generator = AdvancedReportGenerator(
        input_dir='${output_plots_dir}',
        output_dir='.',
        template_path='${workflow.projectDir}/advanced_spatial_report_template_enhanced.html',
        tissue_image_path=tissue_image_arg
    )
    
    # Generate report
    try:
        report_path = generator.generate_report(
            sample_id='${params.sample_id ?: "Unknown"}',
            tissue_type='${params.tissue_type ?: "Unknown"}',
            resolution='${params.resolution ?: "Unknown"}'
        )
        print(f"Report generated successfully: {report_path}")
        
        # Rename report file to match Nextflow output expectations
        import shutil
        expected_filename = f"Celatlas_Omniverse_report_${params.sample_id ?: 'Unknown'}.html"
        if report_path.name != expected_filename:
            shutil.copy2(report_path, expected_filename)
            print(f"Report copied as: {expected_filename}")
        
        # Copy metadata file
        metadata_path = report_path.with_suffix('.json')
        if metadata_path.exists():
            metadata_filename = f"Celatlas_Omniverse_report_${params.sample_id ?: 'Unknown'}.json"
            shutil.copy2(metadata_path, metadata_filename)
            print(f"Metadata copied as: {metadata_filename}")
            
    except Exception as e:
        print(f"Report generation failed: {e}")
        import traceback
        traceback.print_exc()
        
        # Create error report
        with open('spatial_analysis_report_error.html', 'w', encoding='utf-8') as f:
            f.write(f'''
            <!DOCTYPE html>
            <html>
            <head>
                <title>Analysis Report Generation Error</title>
                <style>
                    body {{ font-family: Arial, sans-serif; margin: 40px; }}
                    .error {{ background: #ffe6e6; padding: 20px; border-radius: 8px; border: 1px solid #ff9999; }}
                    .sample-info {{ background: #e6f3ff; padding: 15px; border-radius: 8px; margin: 20px 0; }}
                </style>
            </head>
            <body>
                <h1>Analysis Report Generation Error</h1>
                <div class="sample-info">
                    <h3>Sample Information</h3>
                    <p><strong>Sample ID:</strong> ${params.sample_id ?: "Unknown"}</p>
                    <p><strong>Tissue Type:</strong> ${params.tissue_type ?: "Unknown"}</p>
                    <p><strong>Resolution:</strong> ${params.resolution ?: "Unknown"}</p>
                    <p><strong>Generation Time:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                </div>
                <div class="error">
                    <h3>Error Details</h3>
                    <p><strong>Error Type:</strong> {type(e).__name__}</p>
                    <p><strong>Error Message:</strong> {str(e)}</p>
                </div>
                <p>Please check input data and system configuration, then rerun the analysis.</p>
            </body>
            </html>
            ''')
    """
}

process OPTIMIZE_REPORT_UI {
    tag "Optimizing report UI and assets"
    publishDir "${params.output_dir}/reports", mode: 'copy'
    
    input:
    path report_file
    path assets_dir, stageAs: 'input_assets'
    
    output:
    path "*_optimized.html", emit: optimized_report
    path "optimized_assets/*", emit: optimized_assets, optional: true
    
    script:
    """
    #!/usr/bin/env python3
    
    import os
    import re
    import base64
    from pathlib import Path
    
    print("Starting report UI optimization...")
    
    # Create optimization output directory
    output_dir = Path('optimized_assets')
    output_dir.mkdir(exist_ok=True)
    
    # Read original report
    report_files = list(Path('.').glob('*.html'))
    if not report_files:
        print("No HTML report files found")
        exit(1)
    
    report_file = report_files[0]  # 取第一个HTML文件
    print(f"Processing report file: {report_file}")
    
    try:
        with open(report_file, 'r', encoding='utf-8') as f:
            html_content = f.read()
    except Exception as e:
        print(f"Failed to read report file: {e}")
        exit(1)
    
    # Add advanced UI enhancements
    ui_enhancements = '''
    <script>
    // Interactive feature enhancements
    document.addEventListener('DOMContentLoaded', function() {
        // Universal image click-to-zoom functionality - excluding Tissue Image
        const images = document.querySelectorAll('.result-item img, .chart-container img, .clickable-image');
        images.forEach(img => {
            // Exclude images in Tissue Image container as it has its own zoom functionality
            if (!img.closest('#tissueImageContainer')) {
                img.style.cursor = 'pointer';
                img.style.transition = 'all 0.3s ease';
                
                // Add hover effects
                img.addEventListener('mouseenter', function() {
                    this.style.transform = 'scale(1.02)';
                    this.style.boxShadow = '0 8px 25px rgba(0,0,0,0.15)';
                });
                
                img.addEventListener('mouseleave', function() {
                    this.style.transform = 'scale(1)';
                    this.style.boxShadow = '0 4px 15px rgba(0,0,0,0.1)';
                });
                
                img.addEventListener('click', function(e) {
                    e.preventDefault();
                    e.stopPropagation();
                    
                    const modal = document.createElement('div');
                    modal.className = 'image-modal';
                    modal.style.cssText = `
                        position: fixed; top: 0; left: 0; width: 100vw; height: 100vh;
                        background: rgba(0,0,0,0.95); z-index: 10000; display: flex;
                        align-items: center; justify-content: center; cursor: pointer;
                        opacity: 0; transition: opacity 0.3s ease;
                    `;
                    
                    const modalImg = document.createElement('img');
                    modalImg.src = this.src;
                    modalImg.alt = this.alt || 'Image enlarged view';
                    modalImg.style.cssText = `
                        max-width: 95vw; max-height: 95vh; border-radius: 12px;
                        box-shadow: 0 25px 80px rgba(0,0,0,0.6);
                        transition: transform 0.3s ease; cursor: pointer;
                        object-fit: contain;
                    `;
                    
                    // Add close button
                    const closeBtn = document.createElement('button');
                    closeBtn.innerHTML = '×';
                    closeBtn.style.cssText = `
                        position: absolute; top: 20px; right: 30px; background: rgba(255,255,255,0.2);
                        border: 2px solid rgba(255,255,255,0.3); color: white; font-size: 30px;
                        width: 50px; height: 50px; border-radius: 50%; cursor: pointer;
                        backdrop-filter: blur(10px); transition: all 0.3s ease;
                        display: flex; align-items: center; justify-content: center;
                    `;
                    
                    closeBtn.addEventListener('mouseenter', function() {
                        this.style.background = 'rgba(255,255,255,0.3)';
                        this.style.transform = 'scale(1.1)';
                    });
                    
                    closeBtn.addEventListener('mouseleave', function() {
                        this.style.background = 'rgba(255,255,255,0.2)';
                        this.style.transform = 'scale(1)';
                    });
                    
                    modal.appendChild(modalImg);
                    modal.appendChild(closeBtn);
                    document.body.appendChild(modal);
                    
                    // Animation entrance
                    requestAnimationFrame(() => {
                        modal.style.opacity = '1';
                    });
                    
                    // Close functionality
                    const closeModal = function() {
                        modal.style.opacity = '0';
                        setTimeout(() => {
                            if (document.body.contains(modal)) {
                                document.body.removeChild(modal);
                            }
                        }, 300);
                    };
                    
                    modal.addEventListener('click', function(e) {
                        if (e.target === modal) {
                            closeModal();
                        }
                    });
                    
                    closeBtn.addEventListener('click', closeModal);
                    
                    // Close with ESC key
                    const escHandler = function(e) {
                        if (e.key === 'Escape') {
                            closeModal();
                            document.removeEventListener('keydown', escHandler);
                        }
                    };
                    document.addEventListener('keydown', escHandler);
                });
            }
        });
        
        // Add scroll animations
        const sections = document.querySelectorAll('.section');
        const observer = new IntersectionObserver((entries) => {
            entries.forEach(entry => {
                if (entry.isIntersecting) {
                    entry.target.style.animation = 'fadeInUp 0.6s ease-out';
                }
            });
        });
        
        sections.forEach(section => observer.observe(section));
        
        // Add back-to-top button
        const backToTop = document.createElement('button');
        backToTop.innerHTML = '↑';
        backToTop.style.cssText = `
            position: fixed; bottom: 30px; right: 30px; width: 50px; height: 50px;
            background: linear-gradient(135deg, #667eea, #764ba2); color: white;
            border: none; border-radius: 50%; font-size: 20px; cursor: pointer;
            box-shadow: 0 4px 15px rgba(102, 126, 234, 0.3); z-index: 1000;
            transition: all 0.3s ease; opacity: 0; transform: translateY(20px);
        `;
        
        document.body.appendChild(backToTop);
        
        window.addEventListener('scroll', function() {
            if (window.pageYOffset > 300) {
                backToTop.style.opacity = '1';
                backToTop.style.transform = 'translateY(0)';
            } else {
                backToTop.style.opacity = '0';
                backToTop.style.transform = 'translateY(20px)';
            }
        });
        
        backToTop.addEventListener('click', function() {
            window.scrollTo({top: 0, behavior: 'smooth'});
        });
    });
    </script>
    
    <style>
    /* Animation effects */
    @keyframes fadeInUp {
        from { opacity: 0; transform: translateY(30px); }
        to { opacity: 1; transform: translateY(0); }
    }
    
    /* Enhanced hover effects */
    .result-item img {
        transition: all 0.3s ease;
    }
    
    .result-item img:hover {
        transform: scale(1.05);
        box-shadow: 0 10px 30px rgba(0,0,0,0.2);
    }
    
    /* Advanced gradient background */
    body {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        background-attachment: fixed;
    }
    
    /* Improved container styles */
    .container {
        backdrop-filter: blur(10px);
        border: 1px solid rgba(255,255,255,0.2);
    }
    
    /* Professional table styles */
    .data-table {
        border-collapse: separate;
        border-spacing: 0;
        overflow: hidden;
    }
    
    .data-table th:first-child {
        border-top-left-radius: 8px;
    }
    
    .data-table th:last-child {
        border-top-right-radius: 8px;
    }
    
    .data-table tr:last-child td:first-child {
        border-bottom-left-radius: 8px;
    }
    
    .data-table tr:last-child td:last-child {
        border-bottom-right-radius: 8px;
    }
    
    /* Statistics card 3D effects */
    .stat-card {
        position: relative;
        overflow: hidden;
    }
    
    .stat-card::before {
        content: '';
        position: absolute;
        top: -50%;
        left: -50%;
        width: 200%;
        height: 200%;
        background: linear-gradient(45deg, transparent, rgba(255,255,255,0.1), transparent);
        transform: rotate(45deg);
        transition: all 0.6s;
        opacity: 0;
    }
    
    .stat-card:hover::before {
        animation: shine 0.6s ease-in-out;
    }
    
    @keyframes shine {
        0% { opacity: 0; transform: translateX(-100%) translateY(-100%) rotate(45deg); }
        50% { opacity: 1; }
        100% { opacity: 0; transform: translateX(100%) translateY(100%) rotate(45deg); }
    }
    
    /* Loading animations */
    .section {
        opacity: 0;
        animation: fadeInUp 0.6s ease-out forwards;
    }
    
    .section:nth-child(1) { animation-delay: 0.1s; }
    .section:nth-child(2) { animation-delay: 0.2s; }
    .section:nth-child(3) { animation-delay: 0.3s; }
    .section:nth-child(4) { animation-delay: 0.4s; }
    .section:nth-child(5) { animation-delay: 0.5s; }
    .section:nth-child(6) { animation-delay: 0.6s; }
    </style>
    '''
    
    # Insert enhancements
    enhanced_content = html_content.replace('</body>', f'{ui_enhancements}</body>')
    
    # Generate optimized filename
    optimized_filename = report_file.stem + '_optimized.html'
    
    # Write optimized report
    try:
        with open(optimized_filename, 'w', encoding='utf-8') as f:
            f.write(enhanced_content)
        print(f"UI optimization completed: {optimized_filename}")
    except Exception as e:
        print(f"Failed to write optimized report: {e}")
        # Copy original file as backup
        import shutil
        shutil.copy2(report_file, optimized_filename)
        print(f"Using original report as backup: {optimized_filename}")
    """
}