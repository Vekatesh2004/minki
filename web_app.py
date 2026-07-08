#!/usr/bin/env python3
"""
Pharmacogenomics Pipeline Web Interface
A Flask web application for the pharmacogenomics analysis pipeline
"""

from flask import Flask, render_template, request, redirect, url_for, flash, jsonify, send_file
from werkzeug.utils import secure_filename
import os
import json
import asyncio
import uuid
from datetime import datetime
from pathlib import Path
import sys

# Add current directory to Python path for module imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from modules.vcf_parser import VCFParser
from modules.vep_annotator import VEPAnnotator
from modules.protein_structure import ProteinStructureAnalyzer
from modules.drug_matcher import DrugMatcher
from modules.report_generator import ReportGenerator

# Initialize Flask app
app = Flask(__name__)
app.secret_key = 'pharmacogenomics_secret_key_' + str(uuid.uuid4())

# Configuration
UPLOAD_FOLDER = 'uploads'
RESULTS_FOLDER = 'results'
ALLOWED_EXTENSIONS = {'vcf', 'vcf.gz'}
MAX_FILE_SIZE = 50 * 1024 * 1024  # 50MB

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['RESULTS_FOLDER'] = RESULTS_FOLDER
app.config['MAX_CONTENT_LENGTH'] = MAX_FILE_SIZE

# Ensure directories exist
Path(UPLOAD_FOLDER).mkdir(exist_ok=True)
Path(RESULTS_FOLDER).mkdir(exist_ok=True)

# Global pipeline instance
pipeline_components = {}

def init_pipeline():
    """Initialize pipeline components"""
    global pipeline_components
    
    try:
        # Load configuration
        with open('config.json', 'r') as f:
            config = json.load(f)
        
        # Initialize components
        pipeline_components = {
            'config': config,
            'vcf_parser': VCFParser(config),
            'vep_annotator': VEPAnnotator(config),
            'structure_analyzer': ProteinStructureAnalyzer(config),
            'drug_matcher': DrugMatcher(config),
            'report_generator': ReportGenerator(config)
        }
        
        return True
    except Exception as e:
        print(f"Error initializing pipeline: {e}")
        return False

def allowed_file(filename):
    """Check if file extension is allowed"""
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS or \
           filename.lower().endswith('.vcf.gz')

@app.route('/')
def index():
    """Main page"""
    return render_template('index.html')

@app.route('/upload', methods=['GET', 'POST'])
def upload_file():
    """File upload page"""
    if request.method == 'POST':
        # Check if file was uploaded
        if 'vcf_file' not in request.files:
            flash('No file selected', 'error')
            return redirect(request.url)
        
        file = request.files['vcf_file']
        sample_id = request.form.get('sample_id', '').strip()
        
        # Validate file
        if file.filename == '':
            flash('No file selected', 'error')
            return redirect(request.url)
        
        if not allowed_file(file.filename):
            flash('Please upload a VCF file (.vcf or .vcf.gz)', 'error')
            return redirect(request.url)
        
        # Validate sample ID
        if not sample_id:
            flash('Please provide a sample ID', 'error')
            return redirect(request.url)
        
        # Save file
        filename = secure_filename(file.filename)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        unique_filename = f"{timestamp}_{filename}"
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], unique_filename)
        
        try:
            file.save(file_path)
            flash(f'File uploaded successfully: {filename}', 'success')
            
            # Redirect to analysis page
            return redirect(url_for('analyze', filename=unique_filename, sample_id=sample_id))
            
        except Exception as e:
            flash(f'Error uploading file: {str(e)}', 'error')
            return redirect(request.url)
    
    return render_template('upload.html')

@app.route('/analyze')
def analyze():
    """Analysis page"""
    filename = request.args.get('filename')
    sample_id = request.args.get('sample_id')
    
    if not filename or not sample_id:
        flash('Missing file or sample ID', 'error')
        return redirect(url_for('upload_file'))
    
    file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    if not os.path.exists(file_path):
        flash('File not found', 'error')
        return redirect(url_for('upload_file'))
    
    return render_template('analyze.html', filename=filename, sample_id=sample_id)

@app.route('/api/analyze', methods=['POST'])
def api_analyze():
    """API endpoint for running analysis"""
    data = request.get_json()
    filename = data.get('filename')
    sample_id = data.get('sample_id')
    analysis_type = data.get('analysis_type', 'basic')  # basic, full, custom
    
    if not filename or not sample_id:
        return jsonify({'success': False, 'error': 'Missing parameters'}), 400
    
    file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    if not os.path.exists(file_path):
        return jsonify({'success': False, 'error': 'File not found'}), 404
    
    try:
        # Run analysis based on type
        if analysis_type == 'basic':
            results = asyncio.run(run_basic_analysis(file_path, sample_id))
        elif analysis_type == 'full':
            results = asyncio.run(run_full_analysis(file_path, sample_id))
        else:
            results = asyncio.run(run_custom_analysis(file_path, sample_id, data))
        
        return jsonify({'success': True, 'results': results})
        
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

async def run_basic_analysis(file_path, sample_id):
    """Run basic VCF parsing and drug matching"""
    vcf_parser = pipeline_components['vcf_parser']
    drug_matcher = pipeline_components['drug_matcher']
    
    # Parse VCF
    vcf_results = await vcf_parser.parse_vcf(file_path, sample_id)
    
    # Get unique genes for drug matching (simplified)
    genes_to_check = ['CYP2D6', 'CYP2C19', 'CYP2C9', 'TPMT', 'DPYD', 'VKORC1']
    drug_results = []
    
    for gene in genes_to_check:
        try:
            drug_result = await drug_matcher.match_drugs(gene)
            if drug_result.get('matched_drugs'):
                drug_results.append(drug_result)
        except Exception as e:
            print(f"Error matching drugs for {gene}: {e}")
    
    return {
        'analysis_type': 'basic',
        'vcf_results': vcf_results,
        'drug_results': drug_results,
        'summary': {
            'total_variants': len(vcf_results.get('variants', [])),
            'genes_checked': len(genes_to_check),
            'drug_matches': len(drug_results)
        }
    }

async def run_full_analysis(file_path, sample_id):
    """Run complete analysis pipeline (limited to avoid API rate limits)"""
    vcf_parser = pipeline_components['vcf_parser']
    drug_matcher = pipeline_components['drug_matcher']
    report_generator = pipeline_components['report_generator']
    
    # Parse VCF
    vcf_results = await vcf_parser.parse_vcf(file_path, sample_id)
    
    # Analyze pharmacogenes
    variants = vcf_results.get('variants', [])
    drug_results = []
    
    # Check known pharmacogenes
    pharmacogenes = ['CYP2D6', 'CYP2C19', 'CYP2C9', 'TPMT', 'DPYD', 'VKORC1']
    
    for gene in pharmacogenes:
        try:
            drug_result = await drug_matcher.match_drugs(gene)
            if drug_result.get('matched_drugs'):
                drug_results.append(drug_result)
        except Exception as e:
            print(f"Error matching drugs for {gene}: {e}")
    
    # Generate report
    report_result = await report_generator.generate_report(
        sample_id, vcf_results, [], [], drug_results
    )
    
    return {
        'analysis_type': 'full',
        'vcf_results': vcf_results,
        'drug_results': drug_results,
        'report_result': report_result,
        'summary': {
            'total_variants': len(variants),
            'pharmacogenes_analyzed': len(pharmacogenes),
            'drug_matches': sum(len(r.get('matched_drugs', [])) for r in drug_results),
            'report_files': report_result.get('report_files', {})
        }
    }

async def run_custom_analysis(file_path, sample_id, options):
    """Run custom analysis based on user options"""
    # For now, same as basic but could be extended
    return await run_basic_analysis(file_path, sample_id)

@app.route('/results/<sample_id>')
def view_results(sample_id):
    """View analysis results"""
    # Look for result files
    results_pattern = f"{sample_id}_pharmacogenomics_*.html"
    results_dir = Path(app.config['RESULTS_FOLDER'])
    
    html_files = list(results_dir.glob(results_pattern))
    
    if html_files:
        # Get the most recent file
        latest_file = max(html_files, key=lambda p: p.stat().st_mtime)
        return send_file(latest_file, as_attachment=False)
    else:
        flash(f'No results found for sample {sample_id}', 'error')
        return redirect(url_for('index'))

@app.route('/api/status')
def api_status():
    """API status endpoint"""
    return jsonify({
        'status': 'active',
        'pipeline_initialized': bool(pipeline_components),
        'components': list(pipeline_components.keys()) if pipeline_components else []
    })

@app.route('/demo')
def demo():
    """Demo page with example data"""
    return render_template('demo.html')

@app.route('/api/demo/analyze')
def api_demo_analyze():
    """Run demo analysis with example VCF"""
    try:
        example_vcf = 'examples/sample_pharmacogenomics.vcf'
        if os.path.exists(example_vcf):
            results = asyncio.run(run_basic_analysis(example_vcf, 'DEMO_SAMPLE'))
            return jsonify({'success': True, 'results': results})
        else:
            return jsonify({'success': False, 'error': 'Demo VCF file not found'}), 404
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

if __name__ == '__main__':
    # Initialize pipeline
    if init_pipeline():
        print("✓ Pipeline initialized successfully")
        print("🚀 Starting web interface...")
        print("📖 Open your browser to: http://localhost:5000")
        app.run(debug=True, host='0.0.0.0', port=5000)
    else:
        print("❌ Failed to initialize pipeline")
        print("Make sure config.json exists and is valid")