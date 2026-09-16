import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 1. Simulate a raw RNA-Seq count matrix (Genes vs Samples)
np.random.seed(42)
genes = [f"Gene_{i}" for i in range(1, 501)]

# Generating synthetic count data for Control vs Treated groups
control_means = np.random.uniform(10, 1000, 500)
control_data = np.random.poisson(
    control_means[:, None], (500, 3)
)  # 3 Control replicates

# Fix: Calculate treated means separately to avoid broadcasting shape mismatches
fold_changes = np.random.choice([0.5, 1.0, 2.0], 500, p=[0.1, 0.8, 0.1])
treated_means = control_means * fold_changes
treated_data = np.random.poisson(treated_means[:, None], (500, 3))

df = pd.DataFrame(
    np.hstack([control_data, treated_data]),
    index=genes,
    columns=[
        "Control_1",
        "Control_2",
        "Control_3",
        "Treated_1",
        "Treated_2",
        "Treated_3",
    ],
)

# 2. Calculate Differential Expression (Log2 Fold Change & simulated p-values)
df["Control_Mean"] = df[["Control_1", "Control_2", "Control_3"]].mean(axis=1)
df["Treated_Mean"] = df[["Treated_1", "Treated_2", "Treated_3"]].mean(axis=1)

# Adding a small pseudocount to prevent division by zero
df["Log2_Fold_Change"] = np.log2(
    (df["Treated_Mean"] + 1) / (df["Control_Mean"] + 1)
)
df["p_value"] = 10 ** -np.random.uniform(
    0, 5, 500
)  # Simulated p-values for demonstration
df["Neg_Log10_Pval"] = -np.log10(df["p_value"])

# 3. Categorize gene regulation significance
df["Expression_Status"] = "Not Significant"
df.loc[
    (df["Log2_Fold_Change"] > 1) & (df["p_value"] < 0.05), "Expression_Status"
] = "Upregulated"
df.loc[
    (df["Log2_Fold_Change"] < -1) & (df["p_value"] < 0.05),
    "Expression_Status",
] = "Downregulated"

print(df["Expression_Status"].value_counts())

# 4. Generate Publication-Ready Volcano Plot
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
)  # Significance threshold
plt.axvline(
    x=1, color="black", linestyle="--", linewidth=0.8
)  # Upregulated threshold
plt.axvline(
    x=-1, color="black", linestyle="--", linewidth=0.8
)  # Downregulated threshold

plt.title("RNA-Seq Differential Expression Volcano Plot", fontsize=14)
plt.xlabel("Log2 Fold Change", fontsize=12)
plt.ylabel("-Log10 (p-value)", fontsize=12)
plt.legend()
plt.grid(True, linestyle=":", alpha=0.5)

# Save output figure
plt.savefig("volcano_plot.png", dpi=300)
print("Pipeline complete! Volcano plot saved as volcano_plot.png")