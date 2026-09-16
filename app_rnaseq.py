import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
import streamlit as st

# Page Configuration
st.set_page_config(
    page_title="RNA-Seq Differential Expression Analyzer",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 RNA-Seq Differential Expression Pipeline")
st.markdown(
    "Upload your raw count matrix CSV to perform differential expression"
    " analysis, calculate Log2 fold changes, generate publication-ready Volcano"
    " plots, and extract significant gene lists."
)

# Sidebar for User Controls
st.sidebar.header("Analysis Parameters")
p_val_threshold = st.sidebar.slider(
    "P-Value Threshold", min_value=0.01, max_value=0.10, value=0.05, step=0.01
)
log2fc_threshold = st.sidebar.slider(
    "Log2 Fold Change Threshold", min_value=0.5, max_value=3.0, value=1.0, step=0.5
)

# File Uploader
uploaded_file = st.file_uploader(
    "Upload Raw Count Matrix (CSV format)", type=["csv"]
)

if uploaded_file is not None:
  # Load Dataset
  df = pd.read_csv(uploaded_file, index_col=0)

  st.subheader("Raw Data Preview")
  st.dataframe(df.head())

  # Column Selection for Replicates
  all_columns = df.columns.tolist()
  st.sidebar.subheader("Sample Selection")
  control_cols = st.sidebar.multiselect(
      "Select Control Columns",
      all_columns,
      default=all_columns[:3] if len(all_columns) >= 6 else all_columns[:1],
  )
  treated_cols = st.sidebar.multiselect(
      "Select Treated Columns",
      all_columns,
      default=(
          all_columns[3:6]
          if len(all_columns) >= 6
          else all_columns[1:2] if len(all_columns) >= 2 else []
      ),
  )

  if control_cols and treated_cols:
    if st.button("Run Differential Expression Analysis"):
      with st.spinner(
          "Running t-tests and calculating fold changes..."
      ):
        # Calculate Means
        df["Control_Mean"] = df[control_cols].mean(axis=1)
        df["Treated_Mean"] = df[treated_cols].mean(axis=1)

        # Log2 Fold Change with pseudocount
        df["Log2_Fold_Change"] = np.log2(
            (df["Treated_Mean"] + 1) / (df["Control_Mean"] + 1)
        )

        # Statistical significance via t-test
        t_stats, p_values = ttest_ind(
            df[treated_cols], df[control_cols], axis=1
        )
        df["p_value"] = p_values
        df["p_value"] = df["p_value"].fillna(1.0)
        df["Neg_Log10_Pval"] = -np.log10(df["p_value"] + 1e-300)

        # Gene Regulation Categorization (Corrected syntax)
        df["Expression_Status"] = "Not Significant"
        df.loc[
            (df["Log2_Fold_Change"] > log2fc_threshold)
            & (df["p_value"] < p_val_threshold),
            "Expression_Status",
        ] = "Upregulated"
        df.loc[
            (df["Log2_Fold_Change"] < -log2fc_threshold)
            & (df["p_value"] < p_val_threshold),
            "Expression_Status",
        ] = "Downregulated"

        # Display Metrics Summary
        st.subheader("Analysis Summary")
        counts = df["Expression_Status"].value_counts()
        col1, col2, col3 = st.columns(3)
        col1.metric("Upregulated Genes", counts.get("Upregulated", 0))
        col2.metric("Downregulated Genes", counts.get("Downregulated", 0))
        col3.metric("Not Significant", counts.get("Not Significant", 0))

        # Generate Volcano Plot
        st.subheader("Volcano Plot Visualization")
        fig, ax = plt.subplots(figsize=(8, 6))

        colors = {
            "Not Significant": "grey",
            "Upregulated": "red",
            "Downregulated": "blue",
        }
        for status, group in df.groupby("Expression_Status"):
          ax.scatter(
              group["Log2_Fold_Change"],
              group["Neg_Log10_Pval"],
              color=colors[status],
              label=status,
              alpha=0.7,
              s=25,
          )

        ax.axhline(
            y=-np.log10(p_val_threshold),
            color="black",
            linestyle="--",
            linewidth=0.8,
            label=f"p = {p_val_threshold}",
        )
        ax.axvline(
            x=log2fc_threshold,
            color="black",
            linestyle="--",
            linewidth=0.8,
            label=f"Log2FC = ±{log2fc_threshold}",
        )
        ax.axvline(x=-log2fc_threshold, color="black", linestyle="--", linewidth=0.8)

        ax.set_title("RNA-Seq Differential Expression Volcano Plot", fontsize=14)
        ax.set_xlabel("Log2 Fold Change", fontsize=12)
        ax.set_ylabel("-Log10 (p-value)", fontsize=12)
        ax.legend()
        ax.grid(True, linestyle=":", alpha=0.5)

        st.pyplot(fig)

        # Extract and Display Significant Genes Table
        st.subheader("Flagged Significant Genes")
        sig_df = df[df["Expression_Status"] != "Not Significant"][
            [
                "Control_Mean",
                "Treated_Mean",
                "Log2_Fold_Change",
                "p_value",
                "Expression_Status",
            ]
        ].sort_values(by="p_value")

        st.dataframe(sig_df)

        # Download Button for Flagged Genes CSV
        csv_data = sig_df.to_csv().encode("utf-8")
        st.download_button(
            label="Download Flagged Genes CSV",
            data=csv_data,
            file_name="flagged_differentially_expressed_genes.csv",
            mime="text/csv",
        )
  else:
    st.warning("Please select at least one control and one treated column.")
else:
  st.info(
      "Awaiting CSV upload. You can test this using your synthetic count data"
      " or any real RNA-Seq count matrix."
  )