import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

# 1. Load your raw count matrix CSV file
file_path = "your_gene_counts.csv"

try:
  df = pd.read_csv(file_path, index_col=0)
  print(f"Successfully loaded dataset with {df.shape[0]} genes.")
except FileNotFoundError:
  print(
      f"Error: Could not find {file_path}. Please place your CSV file in the"
      " directory."
  )
  exit()

# 2. Define your control and treated sample column names based on your CSV headers
control_cols = ["Control_1", "Control_2", "Control_3"]
treated_cols = ["Treated_1", "Treated_2", "Treated_3"]

# 3. Calculate Means
df["Control_Mean"] = df[control_cols].mean(axis=1)
df["Treated_Mean"] = df[treated_cols].mean(axis=1)

# 4. Calculate Log2 Fold Change (with a pseudocount of +1 to avoid division by zero)
df["Log2_Fold_Change"] = np.log2(
    (df["Treated_Mean"] + 1) / (df["Control_Mean"] + 1)
)

# 5. Calculate statistical significance using row-wise independent t-tests
t_stats, p_values = ttest_ind(df[treated_cols], df[control_cols], axis=1)
df["p_value"] = p_values
df["p_value"] = df["p_value"].fillna(1.0)
df["Neg_Log10_Pval"] = -np.log10(df["p_value"] + 1e-300)

# 6. Categorize gene regulation significance
df["Expression_Status"] = "Not Significant"
df.loc[
    (df["Log2_Fold_Change"] > 1) & (df["p_value"] < 0.05), "Expression_Status"
] = "Upregulated"
df.loc[
    (df["Log2_Fold_Change"] < -1) & (df["p_value"] < 0.05),
    "Expression_Status",
] = "Downregulated"

print("\nDifferential Expression Summary:")
print(df["Expression_Status"].value_counts())

# 7. Generate Publication-Ready Volcano Plot
plt.figure(figsize=(8, 6))

colors = {
    "Not Significant": "grey",
    "Upregulated": "red",
    "Downregulated": "blue",
}
for status, group in df.groupby("Expression_Status"):
  plt.scatter(
      group["Log2_Fold_Change"],
      group["Neg_Log10_Pval"],
      color=colors[status],
      label=status,
      alpha=0.7,
      s=20,
  )

plt.axhline(
    y=-np.log10(0.05), color="black", linestyle="--", linewidth=0.8
)  # p = 0.05 threshold
plt.axvline(x=1, color="black", linestyle="--", linewidth=0.8)  # Up threshold
plt.axvline(x=-1, color="black", linestyle="--", linewidth=0.8)  # Down threshold

plt.title("RNA-Seq Differential Expression Volcano Plot", fontsize=14)
plt.xlabel("Log2 Fold Change", fontsize=12)
plt.ylabel("-Log10 (p-value)", fontsize=12)
plt.legend()
plt.grid(True, linestyle=":", alpha=0.5)

plt.savefig("real_volcano_plot.png", dpi=300)
print("\nVolcano plot saved as real_volcano_plot.png")

# ==========================================
# 8. EXTRACT & EXPORT FLAGGED GENE LISTS
# ==========================================

# Filter out only the significant genes
significant_genes_df = df[df["Expression_Status"] != "Not Significant"][
    [
        "Control_Mean",
        "Treated_Mean",
        "Log2_Fold_Change",
        "p_value",
        "Expression_Status",
    ]
]

# Sort by statistical significance (lowest p-value first)
significant_genes_df = significant_genes_df.sort_values(by="p_value")

# Save the flagged genes list to a CSV file for the researcher/client
output_csv_name = "flagged_differentially_expressed_genes.csv"
significant_genes_df.to_csv(output_csv_name)
print(
    f"\nExported {len(significant_genes_df)} significant genes to"
    f" '{output_csv_name}'"
)

# Print top upregulated and downregulated genes directly to the console
print("\n--- Top Upregulated Genes ---")
print(
    significant_genes_df[
        significant_genes_df["Expression_Status"] == "Upregulated"
    ]
    .head(5)
    [[
        "Log2_Fold_Change",
        "p_value",
    ]]
)

print("\n--- Top Downregulated Genes ---")
print(
    significant_genes_df[
        significant_genes_df["Expression_Status"] == "Downregulated"
    ]
    .head(5)
    [[
        "Log2_Fold_Change",
        "p_value",
    ]]
)