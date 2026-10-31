#!/usr/bin/env python3
"""
Pharmacogenomics Pipeline - Streamlit Web Interface
A Streamlit web application for pharmacogenomics analysis from VCF files.
"""

import asyncio
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import streamlit as st

# Add project directory to Python path for module imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pharmacogenomics_mcp import PharmacogenomicsMCP

st.set_page_config(
    page_title="Pharmacogenomics Pipeline",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

UPLOAD_FOLDER = "uploads"
RESULTS_FOLDER = "results"
EXAMPLE_VCF = "examples/sample_pharmacogenomics.vcf"
CONFIG_PATH = "config.json"

Path(UPLOAD_FOLDER).mkdir(exist_ok=True)
Path(RESULTS_FOLDER).mkdir(exist_ok=True)


def init_pipeline():
    """Initialize the pipeline components."""
    try:
        pipeline = PharmacogenomicsMCP(config_path=CONFIG_PATH)
        return pipeline, None
    except Exception as e:
        return None, str(e)


def allowed_file(filename):
    """Check if file extension is allowed."""
    return filename.lower().endswith((".vcf", ".vcf.gz"))


def save_uploaded_file(uploaded_file):
    """Save an uploaded file to the uploads folder."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_name = "".join(c if c.isalnum() or c in "._-" else "_" for c in uploaded_file.name)
    unique_filename = f"{timestamp}_{safe_name}"
    file_path = os.path.join(UPLOAD_FOLDER, unique_filename)
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    return file_path


def run_analysis(pipeline, vcf_file, sample_id):
    """Run the full pharmacogenomics pipeline asynchronously."""
    try:
        results = asyncio.run(
            pipeline.run_full_pipeline(vcf_file, sample_id, output_dir=RESULTS_FOLDER)
        )
        return results, None
    except Exception as e:
        return None, str(e)


def format_percentage(value, total):
    """Format a count as a percentage."""
    if total == 0:
        return "0%"
    return f"{value / total * 100:.1f}%"


def display_qc_metrics(vcf_results):
    """Display quality control metrics."""
    qc = vcf_results.get("qc_summary", {})
    total = vcf_results.get("total_variants", 0)

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Total Variants", total)
    with col2:
        st.metric("Pass Rate", format_percentage(qc.get("pass_variants", 0), total))
    with col3:
        st.metric("High Confidence", format_percentage(qc.get("high_confidence_variants", 0), total))
    with col4:
        st.metric("Ts/Tv Ratio", f"{qc.get('ts_tv_ratio', 0):.2f}")

    with st.expander("Detailed QC Statistics"):
        if qc.get("qual_stats"):
            st.write("**Quality Scores:**")
            st.json(qc["qual_stats"])
        if qc.get("depth_stats"):
            st.write("**Depth Stats:**")
            st.json(qc["depth_stats"])
        if qc.get("genotype_stats"):
            st.write("**Genotype Stats:**")
            st.json(qc["genotype_stats"])


def display_variants_table(vcf_results):
    """Display the variants table."""
    variants = vcf_results.get("variants", [])
    if not variants:
        st.info("No variants found.")
        return

    st.subheader(f"Variants ({len(variants)})")

    # Confidence filter
    confidence_filter = st.multiselect(
        "Filter by confidence",
        options=["high", "low"],
        default=["high", "low"],
    )

    filtered = [v for v in variants if v.get("confidence") in confidence_filter]

    import pandas as pd

    rows = []
    for v in filtered:
        rows.append(
            {
                "Chrom": v.get("chrom", ""),
                "Pos": v.get("pos", ""),
                "ID": v.get("id", "."),
                "Ref": v.get("ref", ""),
                "Alt": v.get("alt", ""),
                "QUAL": v.get("qual", ""),
                "FILTER": v.get("filter", ""),
                "Confidence": v.get("confidence", ""),
            }
        )

    df = pd.DataFrame(rows)
    st.dataframe(df, use_container_width=True, hide_index=True)

    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="Download Variants CSV",
        data=csv,
        file_name="variants.csv",
        mime="text/csv",
    )


def display_summary(summary):
    """Display the analysis summary."""
    st.subheader("Analysis Summary")

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Total Variants", summary.get("total_variants_analyzed", 0))
    with col2:
        st.metric("Coding Variants", summary.get("coding_variants", 0))
    with col3:
        st.metric("PGx Relevant", summary.get("pharmacogenomically_relevant", 0))
    with col4:
        st.metric("Genes Analyzed", summary.get("genes_analyzed", 0))

    st.write("**Clinical Significance Breakdown:**")
    sig = summary.get("clinical_significance_breakdown", {})
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("High", sig.get("high", 0))
    with col2:
        st.metric("Moderate", sig.get("moderate", 0))
    with col3:
        st.metric("Low", sig.get("low", 0))

    st.write("**Recommendation:**")
    st.info(summary.get("recommendation", ""))


def display_drug_results(drug_results):
    """Display drug matching results."""
    st.subheader("Drug Interactions")

    matched = [d for d in drug_results if d.get("matched_drugs")]
    if not matched:
        st.info("No drug interactions found.")
        return

    for drug_result in drug_results:
        gene = drug_result.get("gene_symbol", "Unknown")
        drugs = drug_result.get("matched_drugs", [])
        if not drugs:
            continue

        with st.expander(f"{gene} ({len(drugs)} drug matches)"):
            import pandas as pd

            rows = []
            for d in drugs:
                rows.append(
                    {
                        "Drug": d.get("drug_name", ""),
                        "Action": d.get("action", ""),
                        "Source": d.get("source", ""),
                        "Evidence Level": d.get("evidence_level", ""),
                        "Evidence Score": d.get("evidence_score", ""),
                    }
                )
            df = pd.DataFrame(rows)
            st.dataframe(df, use_container_width=True, hide_index=True)


def display_reports(report_result):
    """Display report download buttons."""
    st.subheader("Download Reports")
    report_files = report_result.get("report_files", {})
    if not report_files:
        st.warning("No report files were generated.")
        return

    for fmt, path in report_files.items():
        if path and os.path.exists(path):
            with open(path, "rb") as f:
                data = f.read()
            mime = "text/html" if fmt == "html" else "application/json"
            file_name = os.path.basename(path)
            st.download_button(
                label=f"Download {fmt.upper()} Report",
                data=data,
                file_name=file_name,
                mime=mime,
            )


def main():
    st.title("🧬 Pharmacogenomics Pipeline")
    st.markdown(
        "Analyze genetic variants from VCF files for drug interactions, "
        "protein structure impacts, and clinical significance."
    )

    st.sidebar.title("Navigation")
    page = st.sidebar.radio(
        "Go to",
        options=["Home", "Upload & Analyze", "Demo"],
        index=0,
    )

    st.sidebar.markdown("---")
    st.sidebar.markdown("### Pipeline Status")
    pipeline, init_error = init_pipeline()
    if init_error:
        st.sidebar.error(f"Pipeline init failed: {init_error}")
        st.error(f"Failed to initialize pipeline: {init_error}")
        st.stop()
    else:
        st.sidebar.success("Pipeline initialized")

    if page == "Home":
        st.header("Welcome")
        st.write(
            "Upload a VCF file to begin pharmacogenomics analysis. "
            "The pipeline performs VCF parsing, quality control, variant "
            "annotation, protein structure analysis, drug matching, and report generation."
        )
        st.markdown(
            """
            **Features:**
            - VCF file upload and parsing
            - Quality control metrics
            - Variant annotation via Ensembl VEP
            - Protein structure analysis via UniProt/AlphaFold
            - Drug interaction matching via DrugBank/PharmGKB
            - HTML and JSON report generation
            """
        )
        st.info("Go to **Upload & Analyze** to get started, or try the **Demo** with sample data.")

    elif page == "Demo":
        st.header("Demo Analysis")
        st.write("Run the analysis pipeline on the bundled example VCF file.")

        if not os.path.exists(EXAMPLE_VCF):
            st.error(f"Demo VCF file not found: {EXAMPLE_VCF}")
        else:
            st.write(f"Using example file: `{EXAMPLE_VCF}`")
            sample_id = st.text_input("Sample ID", value="DEMO_SAMPLE", key="demo_sample_id")

            if st.button("Run Demo Analysis", type="primary"):
                with st.spinner("Running analysis pipeline..."):
                    results, error = run_analysis(pipeline, EXAMPLE_VCF, sample_id)

                if error:
                    st.error(f"Analysis failed: {error}")
                else:
                    st.success("Analysis completed!")
                    st.session_state["demo_results"] = results
                    st.session_state["demo_sample_id"] = sample_id

            if "demo_results" in st.session_state:
                results = st.session_state["demo_results"]
                display_results(results)

    elif page == "Upload & Analyze":
        st.header("Upload & Analyze VCF")

        uploaded_file = st.file_uploader(
            "Choose a VCF file",
            type=["vcf", "vcf.gz"],
            help="Upload a VCF or VCF.GZ file containing genetic variants.",
        )
        sample_id = st.text_input("Sample ID", value="", key="upload_sample_id")

        if uploaded_file is not None:
            if not allowed_file(uploaded_file.name):
                st.error("Please upload a valid VCF file (.vcf or .vcf.gz).")
            elif not sample_id:
                st.warning("Please provide a sample ID.")
            else:
                st.info(f"File: {uploaded_file.name} ({uploaded_file.size / 1024:.1f} KB)")
                if st.button("Run Analysis", type="primary"):
                    file_path = save_uploaded_file(uploaded_file)
                    st.write(f"Saved to: `{file_path}`")

                    with st.spinner("Running analysis pipeline... this may take a while"):
                        results, error = run_analysis(pipeline, file_path, sample_id)

                    if error:
                        st.error(f"Analysis failed: {error}")
                    else:
                        st.success("Analysis completed!")
                        st.session_state["upload_results"] = results
                        st.session_state["upload_sample_id"] = sample_id

        if "upload_results" in st.session_state:
            results = st.session_state["upload_results"]
            display_results(results)


def display_results(results):
    """Render the full results view."""
    if results.get("pipeline_status") == "failed":
        st.error(f"Analysis failed: {results.get('error', 'Unknown error')}")
        return

    st.success(f"Analysis status: {results.get('pipeline_status', 'completed')}")

    summary = results.get("summary", {})
    vcf_results = results.get("vcf_results", {})
    annotation_results = results.get("annotation_results", [])
    structure_results = results.get("structure_results", [])
    drug_results = results.get("drug_results", [])
    report_result = results.get("report_result", {})

    tab_qc, tab_variants, tab_genes, tab_drugs, tab_reports = st.tabs(
        ["QC Metrics", "Variants", "Genes", "Drug Interactions", "Reports"]
    )

    with tab_qc:
        display_qc_metrics(vcf_results)
        st.markdown("---")
        display_summary(summary)

    with tab_variants:
        display_variants_table(vcf_results)

    with tab_genes:
        genes = summary.get("genes_with_variants", [])
        if genes:
            st.write("**Genes with variants:**")
            st.write(", ".join(genes))
        else:
            st.info("No genes with variants found.")

    with tab_drugs:
        display_drug_results(drug_results)

    with tab_reports:
        display_reports(report_result)

    st.markdown("---")
    st.caption(
        "This tool is for research purposes only and should not be used for "
        "clinical decision-making without proper validation and consultation "
        "with healthcare professionals."
    )


if __name__ == "__main__":
    main()
