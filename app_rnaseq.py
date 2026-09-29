import streamlit as st
import pandas as pd
import numpy as np

# --- PAGE CONFIGURATION ---
st.set_page_config(
    page_title="Tier 2 Live Showcase | Transcriptomics (RNA-Seq)",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# --- LIVE LINKS ---
MAIN_WEBSITE_URL = "https://bhumikasihare.github.io/genometech-studio/"
CHECKOUT_TIER2_URL = "https://genometech-checkout.streamlit.app/?tier=2"

# --- BLUISH & REDDISH THEME + SINGLE CENTER MOVING LOGO ---
st.markdown(
    '<style>'
    '.stApp {'
    '  background: radial-gradient(circle at 50% 35%, #172a54 0%, #091326 55%, #040812 100%);'
    '  color: #f8fafc;'
    '  font-family: "Inter", -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;'
    '}'
    'section[data-testid="stSidebar"] {'
    '  background-color: rgba(5, 10, 22, 0.96) !important;'
    '  border-right: 1px solid rgba(239, 68, 68, 0.35);'
    '}'
    '@keyframes floatCenterTier2 {'
    '  0% { transform: translate(-50%, -50%) translateY(0px) scale(1) rotate(0deg); opacity: 0.11; }'
    '  50% { transform: translate(-50%, -50%) translateY(-22px) scale(1.06) rotate(5deg); opacity: 0.22; }'
    '  100% { transform: translate(-50%, -50%) translateY(0px) scale(1) rotate(0deg); opacity: 0.11; }'
    '}'
    '.center-tier2-logo {'
    '  position: fixed;'
    '  top: 52%;'
    '  left: 50%;'
    '  font-size: 13rem;'
    '  pointer-events: none;'
    '  z-index: 0;'
    '  user-select: none;'
    '  filter: drop-shadow(0 0 35px rgba(239, 68, 68, 0.55)) drop-shadow(0 0 60px rgba(59, 130, 246, 0.4));'
    '  animation: floatCenterTier2 10s ease-in-out infinite;'
    '}'
    '.t2-card {'
    '  background: rgba(11, 24, 48, 0.84);'
    '  backdrop-filter: blur(14px);'
    '  border: 1px solid rgba(239, 68, 68, 0.38);'
    '  padding: 1.4rem 1.6rem;'
    '  border-radius: 16px;'
    '  box-shadow: 0 14px 34px rgba(0, 0, 0, 0.5);'
    '  margin-bottom: 1.1rem;'
    '  position: relative;'
    '  z-index: 1;'
    '}'
    '.t2-kpi {'
    '  background: linear-gradient(135deg, rgba(15, 32, 66, 0.92), rgba(127, 29, 29, 0.28));'
    '  border: 1px solid rgba(244, 63, 94, 0.45);'
    '  padding: 1.15rem;'
    '  border-radius: 14px;'
    '  text-align: center;'
    '  position: relative;'
    '  z-index: 1;'
    '}'
    '.nav-btn-back {'
    '  display: inline-flex;'
    '  align-items: center;'
    '  justify-content: center;'
    '  gap: 8px;'
    '  background: rgba(239, 68, 68, 0.15);'
    '  color: #fca5a5 !important;'
    '  border: 1px solid rgba(239, 68, 68, 0.55);'
    '  padding: 9px 18px;'
    '  border-radius: 999px;'
    '  font-size: 0.86rem;'
    '  font-weight: 800;'
    '  text-decoration: none !important;'
    '  transition: all 0.2s ease;'
    '}'
    '.nav-btn-back:hover {'
    '  background: #ef4444;'
    '  color: #ffffff !important;'
    '  box-shadow: 0 0 18px rgba(239, 68, 68, 0.55);'
    '}'
    '.nav-btn-order {'
    '  display: inline-flex;'
    '  align-items: center;'
    '  justify-content: center;'
    '  gap: 8px;'
    '  background: linear-gradient(90deg, #dc2626, #f43f5e);'
    '  color: #ffffff !important;'
    '  padding: 9px 18px;'
    '  border-radius: 999px;'
    '  font-size: 0.86rem;'
    '  font-weight: 800;'
    '  text-decoration: none !important;'
    '  box-shadow: 0 6px 20px rgba(239, 68, 68, 0.4);'
    '}'
    'div.stDownloadButton > button {'
    '  background: linear-gradient(90deg, #dc2626, #ef4444) !important;'
    '  color: #ffffff !important;'
    '  font-weight: 800 !important;'
    '  border: none !important;'
    '  border-radius: 10px !important;'
    '  padding: 0.65rem 1.2rem !important;'
    '}'
    '</style>'
    '<div class="center-tier2-logo">🧬</div>',
    unsafe_allow_html=True
)

# --- PRE-LOADED DEMO BULK RNA-SEQ DATA (NO UPLOAD NEEDED) ---
@st.cache_data
def get_tier2_demo_data():
    genes = [
        {"Gene_Symbol": "MYC", "Ensembl_ID": "ENSG00000136997", "baseMean": 4820.5, "log2FoldChange": 3.42, "pvalue": 1.2e-18, "padj": 4.5e-16, "Regulation": "Upregulated", "Pathway": "Cell Cycle & Proliferation", "CTRL_Mean_TPM": 42.1, "TUMOR_Mean_TPM": 451.8},
        {"Gene_Symbol": "VEGFA", "Ensembl_ID": "ENSG00000112715", "baseMean": 3190.2, "log2FoldChange": 2.89, "pvalue": 3.4e-15, "padj": 8.1e-13, "Regulation": "Upregulated", "Pathway": "Angiogenesis Signaling", "CTRL_Mean_TPM": 38.4, "TUMOR_Mean_TPM": 284.6},
        {"Gene_Symbol": "EGFR", "Ensembl_ID": "ENSG00000146648", "baseMean": 5640.8, "log2FoldChange": 2.64, "pvalue": 2.1e-14, "padj": 3.9e-12, "Regulation": "Upregulated", "Pathway": "RTK / MAPK Cascade", "CTRL_Mean_TPM": 64.0, "TUMOR_Mean_TPM": 398.5},
        {"Gene_Symbol": "MKI67", "Ensembl_ID": "ENSG00000148773", "baseMean": 2890.1, "log2FoldChange": 3.15, "pvalue": 8.5e-17, "padj": 2.2e-14, "Regulation": "Upregulated", "Pathway": "Mitotic Spindle Assembly", "CTRL_Mean_TPM": 29.5, "TUMOR_Mean_TPM": 261.9},
        {"Gene_Symbol": "CCND1", "Ensembl_ID": "ENSG00000110092", "baseMean": 3410.7, "log2FoldChange": 2.18, "pvalue": 1.9e-11, "padj": 1.8e-9, "Regulation": "Upregulated", "Pathway": "G1/S Transition", "CTRL_Mean_TPM": 58.2, "TUMOR_Mean_TPM": 263.4},
        {"Gene_Symbol": "MMP9", "Ensembl_ID": "ENSG00000100985", "baseMean": 1980.4, "log2FoldChange": 2.51, "pvalue": 6.2e-13, "padj": 9.4e-11, "Regulation": "Upregulated", "Pathway": "Extracellular Matrix Remodeling", "CTRL_Mean_TPM": 24.8, "TUMOR_Mean_TPM": 141.2},
        {"Gene_Symbol": "TP53", "Ensembl_ID": "ENSG00000141510", "baseMean": 4120.0, "log2FoldChange": -2.76, "pvalue": 4.0e-15, "padj": 8.9e-13, "Regulation": "Downregulated", "Pathway": "Apoptosis & DNA Repair", "CTRL_Mean_TPM": 312.4, "TUMOR_Mean_TPM": 46.1},
        {"Gene_Symbol": "PTEN", "Ensembl_ID": "ENSG00000171862", "baseMean": 3680.9, "log2FoldChange": -2.44, "pvalue": 1.5e-12, "padj": 2.1e-10, "Regulation": "Downregulated", "Pathway": "PI3K/AKT Inhibition", "CTRL_Mean_TPM": 275.0, "TUMOR_Mean_TPM": 50.8},
        {"Gene_Symbol": "CDKN1A", "Ensembl_ID": "ENSG00000124762", "baseMean": 2950.3, "log2FoldChange": -2.12, "pvalue": 3.8e-10, "padj": 3.1e-8, "Regulation": "Downregulated", "Pathway": "Cell Cycle Arrest (p21)", "CTRL_Mean_TPM": 218.6, "TUMOR_Mean_TPM": 50.2},
        {"Gene_Symbol": "BAX", "Ensembl_ID": "ENSG00000087088", "baseMean": 2140.6, "log2FoldChange": -1.89, "pvalue": 2.4e-8, "padj": 1.4e-6, "Regulation": "Downregulated", "Pathway": "Mitochondrial Apoptosis", "CTRL_Mean_TPM": 184.2, "TUMOR_Mean_TPM": 49.7},
        {"Gene_Symbol": "FOXO3", "Ensembl_ID": "ENSG00000118689", "baseMean": 2530.2, "log2FoldChange": -1.74, "pvalue": 9.1e-8, "padj": 4.8e-6, "Regulation": "Downregulated", "Pathway": "Oxidative Stress Response", "CTRL_Mean_TPM": 166.9, "TUMOR_Mean_TPM": 49.9},
        {"Gene_Symbol": "GAPDH", "Ensembl_ID": "ENSG00000111640", "baseMean": 12450.0, "log2FoldChange": 0.08, "pvalue": 0.68, "padj": 0.89, "Regulation": "Non-Significant", "Pathway": "Housekeeping Glycolysis", "CTRL_Mean_TPM": 890.2, "TUMOR_Mean_TPM": 940.5},
        {"Gene_Symbol": "ACTB", "Ensembl_ID": "ENSG00000075624", "baseMean": 14120.5, "log2FoldChange": -0.05, "pvalue": 0.79, "padj": 0.94, "Regulation": "Non-Significant", "Pathway": "Cytoskeletal Housekeeping", "CTRL_Mean_TPM": 1020.4, "TUMOR_Mean_TPM": 985.1},
        {"Gene_Symbol": "RPLP0", "Ensembl_ID": "ENSG00000089157", "baseMean": 9820.1, "log2FoldChange": 0.14, "pvalue": 0.51, "padj": 0.78, "Regulation": "Non-Significant", "Pathway": "Ribosomal Subunit", "CTRL_Mean_TPM": 640.0, "TUMOR_Mean_TPM": 705.2},
    ]
    df = pd.DataFrame(genes)
    df["neg_log10_padj"] = -np.log10(df["padj"])
    return df

df_deg = get_tier2_demo_data()

# --- TOP HEADER BAR ---
top_left, top_right = st.columns([3.2, 1.8])
with top_left:
    st.markdown(
        '<div style="display:flex; align-items:center; gap:14px;">'
        '<div style="width:56px; height:56px; border-radius:14px; background:rgba(239,68,68,0.18); border:1px solid #ef4444; display:flex; align-items:center; justify-content:center; font-size:2rem;">🧬</div>'
        '<div>'
        '<div style="font-size:0.76rem; text-transform:uppercase; letter-spacing:0.15em; color:#f87171; font-weight:800;">💎 Tier 2 Live Deliverable Showcase • Pre-Loaded Bulk RNA-Seq Cohort</div>'
        '<div style="font-size:2rem; font-weight:900; color:#ffffff; letter-spacing:-0.02em;">Transcriptomics & <span style="color:#ef4444;">Differential Expression</span></div>'
        '</div>'
        '</div>',
        unsafe_allow_html=True
    )
with top_right:
    st.markdown(
        f'<div style="display:flex; justify-content:flex-end; align-items:center; gap:10px; padding-top:10px; flex-wrap:wrap;">'
        f'<a href="{MAIN_WEBSITE_URL}" target="_blank" rel="noopener noreferrer" class="nav-btn-back">← Back to Our Website</a>'
        f'<a href="{CHECKOUT_TIER2_URL}" target="_blank" rel="noopener noreferrer" class="nav-btn-order">Select Tier 2 Modules ($410) ➔</a>'
        f'</div>',
        unsafe_allow_html=True
    )

st.markdown("<hr style='border-color: rgba(239, 68, 68, 0.3); margin: 0.85rem 0 1.2rem 0;'>", unsafe_allow_html=True)

# --- CONCISE OVERVIEW: WHAT WE DELIVER + FILES & LANGUAGES USED ---
st.markdown(
    '<div class="t2-card" style="background:linear-gradient(90deg, rgba(12,25,52,0.94), rgba(127,29,29,0.28));">'
    '<div style="font-size:1.1rem; font-weight:800; color:#f87171; margin-bottom:6px;">🔬 What You Receive in Tier 2 (Live Demo Dataset: Control vs. Tumor RNA-Seq Cohort, n=6)</div>'
    '<div style="font-size:0.9rem; color:#e2e8f0; line-height:1.55; margin-bottom:14px;">'
    'This live interactive dashboard showcases the exact deliverables provided in our <b>Tier 2 Transcriptomics Bundle</b>. '
    'Starting from raw RNA-Seq reads, we execute splice-aware genome alignment, exon-level read quantification, negative binomial differential expression modeling, and publication-grade biomarker visualizations.'
    '</div>'
    '<div style="display:grid; grid-template-columns:repeat(auto-fit, minmax(210px, 1fr)); gap:10px; font-size:0.83rem;">'
    '<div style="background:rgba(4,8,18,0.65); border:1px solid rgba(239,68,68,0.3); padding:10px 12px; border-radius:10px;">'
    '<b style="color:#fca5a5;">📂 Alignment Deliverables:</b><br><span style="color:#cbd5e1;">Raw <code>.fastq.gz</code> ➔ Splice-aligned <code>.bam</code> + splice junction logs</span>'
    '</div>'
    '<div style="background:rgba(4,8,18,0.65); border:1px solid rgba(59,130,246,0.35); padding:10px 12px; border-radius:10px;">'
    '<b style="color:#93c5fd;">🧮 Count Matrices:</b><br><span style="color:#cbd5e1;">Raw <code>featureCounts.txt</code> + normalized <code>TPM / FPKM</code> tables</span>'
    '</div>'
    '<div style="background:rgba(4,8,18,0.65); border:1px solid rgba(239,68,68,0.3); padding:10px 12px; border-radius:10px;">'
    '<b style="color:#fca5a5;">📊 Statistical & Plot Outputs:</b><br><span style="color:#cbd5e1;">Full <code>DESeq2_DEG.csv</code> table + high-res Volcano & Heatmap figures</span>'
    '</div>'
    '<div style="background:rgba(4,8,18,0.65); border:1px solid rgba(59,130,246,0.35); padding:10px 12px; border-radius:10px;">'
    '<b style="color:#93c5fd;">💻 Languages & Core Stack:</b><br><span style="color:#cbd5e1;"><code>R / Bioconductor</code> (DESeq2), <code>Linux/Bash</code> (STAR), & <code>Python</code></span>'
    '</div>'
    '</div>'
    '</div>',
    unsafe_allow_html=True
)

# --- TOP-LEVEL DEMO COHORT KPI METRICS ---
k1, k2, k3, k4 = st.columns(4)
with k1:
    st.markdown(
        '<div class="t2-kpi">'
        '<div style="font-size:0.76rem; text-transform:uppercase; color:#93c5fd; font-weight:700;">Uniquely Mapped Reads</div>'
        '<div style="font-size:1.95rem; font-weight:900; color:#ffffff; margin:4px 0;">94.6%</div>'
        '<div style="font-size:0.78rem; color:#fca5a5;">Module 1 • STAR Alignment ($150)</div>'
        '</div>',
        unsafe_allow_html=True
    )
with k2:
    st.markdown(
        '<div class="t2-kpi">'
        '<div style="font-size:0.76rem; text-transform:uppercase; color:#93c5fd; font-weight:700;">Quantified Human Genes</div>'
        '<div style="font-size:1.95rem; font-weight:900; color:#ffffff; margin:4px 0;">19,840</div>'
        '<div style="font-size:0.78rem; color:#fca5a5;">Module 2 • FeatureCounts Quant ($80)</div>'
        '</div>',
        unsafe_allow_html=True
    )
with k3:
    st.markdown(
        '<div class="t2-kpi">'
        '<div style="font-size:0.76rem; text-transform:uppercase; color:#93c5fd; font-weight:700;">Significant DEGs (padj &lt; 0.05)</div>'
        '<div style="font-size:1.95rem; font-weight:900; color:#ef4444; margin:4px 0;">1,428</div>'
        '<div style="font-size:0.78rem; color:#fca5a5;">Module 3 • DESeq2 Math ($120)</div>'
        '</div>',
        unsafe_allow_html=True
    )
with k4:
    st.markdown(
        '<div class="t2-kpi">'
        '<div style="font-size:0.76rem; text-transform:uppercase; color:#93c5fd; font-weight:700;">Top Biomarker Panel</div>'
        '<div style="font-size:1.95rem; font-weight:900; color:#ffffff; margin:4px 0;">11 Genes</div>'
        '<div style="font-size:0.78rem; color:#fca5a5;">Module 4 • Volcano & Heatmap ($60)</div>'
        '</div>',
        unsafe_allow_html=True
    )

st.write("")

# --- 4-MODULE DELIVERABLE TABS ---
tab1, tab2, tab3, tab4 = st.tabs([
    "📊 Module 3: DESeq2 Differential Math ($120)",
    "🌋 Module 4: Volcano & Expression Vis ($60)",
    "🎯 Module 1: STAR Read Alignment ($150)",
    "🧮 Module 2: FeatureCounts Quant ($80)"
])

# =====================================================================
# TAB 1: MODULE 3 — DESEQ2 DIFFERENTIAL EXPRESSION TABLE
# =====================================================================
with tab1:
    st.markdown(
        '<div class="t2-card">'
        '<div style="font-size:1.15rem; font-weight:800; color:#f87171; margin-bottom:4px;">📊 DESeq2 Differential Expression Deliverable (Control vs. Tumor)</div>'
        '<div style="font-size:0.86rem; color:#cbd5e1;">'
        'Using negative binomial generalized linear models in <code>R (DESeq2)</code>, we calculate shrinkage-corrected <code>log2FoldChange</code> and Benjamini-Hochberg adjusted p-values (<code>padj</code>). Filter the pre-loaded biomarker dataset below:'
        '</div>'
        '</div>',
        unsafe_allow_html=True
    )

    f1, f2 = st.columns(2)
    with f1:
        reg_filter = st.multiselect(
            "Filter by Regulation Direction:",
            options=["Upregulated", "Downregulated", "Non-Significant"],
            default=["Upregulated", "Downregulated", "Non-Significant"],
            key="deg_reg_filter"
        )
    with f2:
        min_lfc = st.slider("Minimum Absolute |log2FoldChange|:", min_value=0.0, max_value=3.0, value=0.0, step=0.25, key="deg_lfc_slider")

    # FIXED FILTERING LOGIC
    filtered_deg = df_deg[
        (df_deg["Regulation"].isin(reg_filter)) &
        (df_deg["log2FoldChange"].abs() >= float(min_lfc))
    ].drop(columns=["neg_log10_padj"])

    st.dataframe(filtered_deg, use_container_width=True, hide_index=True)

    d_col1, d_col2 = st.columns([1.2, 2.8])
    with d_col1:
        csv_deg = filtered_deg.to_csv(index=False).encode("utf-8")
        st.download_button(
            label="📥 Download Sample DESeq2 CSV",
            data=csv_deg,
            file_name="GenomeTech_Tier2_DESeq2_Results.csv",
            mime="text/csv",
            use_container_width=True
        )
    with d_col2:
        st.caption("💡 Generated via R (Bioconductor DESeq2) with apeglm log2FC shrinkage and False Discovery Rate (FDR) multiple-testing correction.")

# =====================================================================
# TAB 2: MODULE 4 — VOLCANO & BIOMARKER EXPRESSION VISUALIZATION
# =====================================================================
with tab2:
    v_col1, v_col2 = st.columns([1.15, 1.15], gap="large")
    with v_col1:
        st.markdown(
            '<div class="t2-card">'
            '<div style="font-size:1.1rem; font-weight:800; color:#f87171; margin-bottom:4px;">🌋 Interactive Volcano Plot (log2FC vs. -log10 padj)</div>'
            '<div style="font-size:0.84rem; color:#cbd5e1;">'
            'Visualizes statistical significance against fold-change magnitude. Oncogenes (right) are strongly upregulated in Tumor samples, while tumor suppressors (left) are downregulated.'
            '</div>'
            '</div>',
            unsafe_allow_html=True
        )
        st.scatter_chart(
            df_deg,
            x="log2FoldChange",
            y="neg_log10_padj",
            color="Regulation",
            size="baseMean",
            height=360
        )

    with v_col2:
        st.markdown(
            '<div class="t2-card">'
            '<div style="font-size:1.1rem; font-weight:800; color:#f87171; margin-bottom:4px;">🔥 Normalized Expression Profile (Mean TPM: Control vs. Tumor)</div>'
            '<div style="font-size:0.84rem; color:#cbd5e1;">'
            'Side-by-side Transcripts Per Million (TPM) expression comparison across top differentially expressed biomarkers.'
            '</div>'
            '</div>',
            unsafe_allow_html=True
        )
        sig_only = df_deg[df_deg["Regulation"] != "Non-Significant"].set_index("Gene_Symbol")[["CTRL_Mean_TPM", "TUMOR_Mean_TPM"]]
        st.bar_chart(sig_only, color=["#3b82f6", "#ef4444"], height=360)

# =====================================================================
# TAB 3: MODULE 1 — STAR SPLICE-AWARE READ ALIGNMENT
# =====================================================================
with tab3:
    s_left, s_right = st.columns([1.1, 1.4], gap="large")
    with s_left:
        st.markdown(
            '<div class="t2-card">'
            '<div style="font-size:1.15rem; font-weight:800; color:#f87171; margin-bottom:6px;">🎯 Module 1: STAR Splice-Aware Read Alignment ($150)</div>'
            '<div style="font-size:0.88rem; color:#e2e8f0; line-height:1.6;">'
            'Unlike DNA aligners, RNA-Seq requires two-pass splice-aware alignment across exon-exon junctions against the human reference transcriptome (<code>GRCh38 + GENCODE v44</code>).<br><br>'
            '• <b>Mean Uniquely Mapped Reads:</b> 94.6% across 6 samples<br>'
            '• <b>Annotated Splice Junctions:</b> 98.2% GENCODE canonical<br>'
            '• <b>Multi-Mapping Reads:</b> &lt; 3.4%<br>'
            '• <b>Deliverable Files:</b> Coordinate-sorted <code>.bam</code>, <code>SJ.out.tab</code> splice junctions, and alignment summary logs'
            '</div>'
            '</div>',
            unsafe_allow_html=True
        )
    with s_right:
        st.markdown("##### 📈 STAR Alignment Breakdown Across Cohort Replicates (%)")
        star_df = pd.DataFrame({
            "Uniquely Mapped (%)": [94.8, 94.2, 95.1, 94.5, 93.9, 95.0],
            "Multi-Mapped (%)": [3.1, 3.5, 2.9, 3.4, 3.8, 3.0],
            "Unmapped (%)": [2.1, 2.3, 2.0, 2.1, 2.3, 2.0]
        }, index=["CTRL_Rep1", "CTRL_Rep2", "CTRL_Rep3", "TUMOR_Rep1", "TUMOR_Rep2", "TUMOR_Rep3"])
        st.bar_chart(star_df, color=["#ef4444", "#3b82f6", "#334155"])

# =====================================================================
# TAB 4: MODULE 2 — FEATURECOUNTS GENE QUANTIFICATION
# =====================================================================
with tab4:
    q_left, q_right = st.columns([1.05, 1.45], gap="large")
    with q_left:
        st.markdown(
            '<div class="t2-card">'
            '<div style="font-size:1.15rem; font-weight:800; color:#f87171; margin-bottom:6px;">🧮 Module 2: FeatureCounts Gene Quantification ($80)</div>'
            '<div style="font-size:0.88rem; color:#e2e8f0; line-height:1.6;">'
            'Subread <code>featureCounts</code> assigns mapped paired-end fragments to GENCODE exon features to build your master count matrix.<br><br>'
            '• <b>Exonic Assignment Rate:</b> 88.4% of aligned fragments<br>'
            '• <b>Strand-Specificity:</b> Reverse-stranded library protocol<br>'
            '• <b>Deliverable Files:</b> Raw integer count matrix (<code>.txt</code> / <code>.csv</code>) ready for DESeq2 + normalized TPM expression matrix'
            '</div>'
            '</div>',
            unsafe_allow_html=True
        )
    with q_right:
        st.markdown("##### 📋 Sample Raw Fragment Count Matrix (6 Replicates)")
        counts_preview = pd.DataFrame({
            "Gene": ["MYC", "VEGFA", "EGFR", "MKI67", "TP53", "PTEN", "CDKN1A", "GAPDH"],
            "CTRL_1": [440, 390, 650, 305, 3210, 2810, 2240, 12300],
            "CTRL_2": [465, 410, 620, 290, 3080, 2740, 2190, 12540],
            "CTRL_3": [425, 380, 660, 315, 3190, 2790, 2150, 12410],
            "TUMOR_1": [4690, 2910, 4050, 2710, 470, 520, 510, 12610],
            "TUMOR_2": [4880, 2840, 3920, 2640, 455, 495, 490, 12390],
            "TUMOR_3": [4750, 2960, 4110, 2690, 485, 515, 525, 12580],
        })
        st.dataframe(counts_preview, use_container_width=True, hide_index=True)