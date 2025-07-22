import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import math

# Simulated gene expression DataFrame (TPM-like)
df = pd.DataFrame({
    'ENSG0001': [5, 10, 7, 9],
    'ENSG0002': [2, 3, 4, 6],
    'ENSG0003': [8, 6, 7, 5]
}, index=['Sample1', 'Sample2', 'Sample3', 'Sample4']).T

symbol_map = {
    'TP53': 'ENSG0001',
    'EGFR': 'ENSG0002',
    'CD274': 'ENSG0003'
}

st.title("Gene Histogram Viewer")

input_genes = st.text_input("Enter gene symbols (comma-separated):", "TP53,EGFR")

if input_genes:
    symbols = [s.strip().upper() for s in input_genes.split(',')]
    ensgs = [symbol_map.get(sym) for sym in symbols if sym in symbol_map]
    missing = [sym for sym in symbols if sym not in symbol_map]

    if missing:
        st.warning(f"Missing symbols: {', '.join(missing)}")

    if ensgs:
        n = len(ensgs)
        cols = 2
        rows = math.ceil(n / cols)
        fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 3))
        axes = axes.flatten()

        for i, ensg in enumerate(ensgs):
            vals = df.loc[ensg]
            axes[i].hist(vals, bins=10, color="skyblue", edgecolor="black")
            axes[i].set_title(f"{symbols[i]} ({ensg})")
            axes[i].set_xlabel("TPM")
            axes[i].set_ylabel("Samples")

        for j in range(i + 1, len(axes)):
            axes[j].axis("off")

        st.pyplot(fig)
