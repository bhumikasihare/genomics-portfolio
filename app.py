import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

# Page Configuration
st.set_page_config(
    page_title="Genomic Variant Analysis Dashboard", page_icon="🧬", layout="wide"
)

st.title("🧬 Genomic Variant Analysis Dashboard")
st.markdown(
    "Upload your variant dataset (CSV format) to explore, filter, and"
    " visualize genomic variations, quality metrics, and chromosomal"
    " distributions."
)

# File Uploader
uploaded_file = st.file_uploader(
    "Upload Variant Dataset (CSV format)", type=["csv"]
)

if uploaded_file is not None:
  # Load Dataset
  df = pd.read_csv(uploaded_file)

  st.subheader("Raw Data Preview")
  st.dataframe(df.head())

  # Dynamic column detection for flexibility with client data
  columns = df.columns.tolist()
  chrom_col = next(
      (
          col
          for col in columns
          if "chrom" in col.lower() or "chr" in col.lower()
      ),
      columns[0] if columns else None,
  )
  qual_col = next(
      (
          col
          for col in columns
          if "qual" in col.lower() or "score" in col.lower()
      ),
      None,
  )

  # Sidebar Controls for Filtering
  st.sidebar.header("Variant Filter Parameters")

  selected_chroms = []
  if chrom_col:
    unique_chroms = df[chrom_col].unique().tolist()
    selected_chroms = st.sidebar.multiselect(
        "Select Chromosomes", unique_chroms, default=unique_chroms
    )

  min_qual = 0.0
  if qual_col:
    min_val = (
        float(df[qual_col].min()) if not df[qual_col].dropna().empty else 0.0
    )
    max_val = (
        float(df[qual_col].max()) if not df[qual_col].dropna().empty else 100.0
    )
    min_qual = st.sidebar.slider(
        "Minimum Quality Score (QUAL)",
        min_value=min_val,
        max_value=max_val,
        value=min_val,
    )

  # Apply Filters
  filtered_df = df.copy()
  if chrom_col and selected_chroms:
    filtered_df = filtered_df[filtered_df[chrom_col].isin(selected_chroms)]
  if qual_col:
    filtered_df = filtered_df[filtered_df[qual_col] >= min_qual]

  # Display Summary Metrics
  st.subheader("Dataset Summary")
  col1, col2, col3 = st.columns(3)
  col1.metric("Total Variants (Loaded)", len(df))
  col2.metric("Filtered Variants", len(filtered_df))
  if qual_col:
    avg_q = (
        round(filtered_df[qual_col].mean(), 2)
        if not filtered_df.empty
        else 0.0
    )
    col3.metric("Avg Quality Score", avg_q)
  else:
    col3.metric("Columns Analyzed", len(columns))

  # Display Filtered Data Table
  st.subheader("Filtered Variants Table")
  st.dataframe(filtered_df)

  # Visualizations
  st.subheader("Variant Visualizations")
  if chrom_col:
    fig, ax = plt.subplots(figsize=(10, 4))
    chrom_counts = filtered_df[chrom_col].value_counts().sort_index()
    chrom_counts.plot(kind="bar", ax=ax, color="skyblue", edgecolor="black")
    ax.set_title("Variant Distribution Across Chromosomes", fontsize=14)
    ax.set_xlabel("Chromosome", fontsize=12)
    ax.set_ylabel("Variant Count", fontsize=12)
    plt.xticks(rotation=45)
    st.pyplot(fig)

  if qual_col:
    fig2, ax2 = plt.subplots(figsize=(10, 4))
    ax2.hist(
        filtered_df[qual_col].dropna(),
        bins=20,
        color="salmon",
        edgecolor="black",
        alpha=0.7,
    )
    ax2.set_title("Quality Score Distribution", fontsize=14)
    ax2.set_xlabel("Quality Score", fontsize=12)
    ax2.set_ylabel("Frequency", fontsize=12)
    st.pyplot(fig2)

  # Download Filtered Data Button
  csv_data = filtered_df.to_csv(index=False).encode("utf-8")
  st.download_button(
      label="Download Filtered Variants CSV",
      data=csv_data,
      file_name="filtered_genomic_variants.csv",
      mime="text/csv",
  )
else:
  st.info("Please upload a variant dataset CSV to launch the dashboard.")