import streamlit as st
import pandas as pd
import numpy as np
import os
from glob import glob
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import MaxNLocator

# ────────────────────────────────────────────────────────────────
# Set Streamlit page config
# ────────────────────────────────────────────────────────────────
st.set_page_config(page_title="CPTAC GBM Gene Histogram Viewer", layout="wide")
st.title("CPTAC GBM Gene Expression Visualizer")

# ────────────────────────────────────────────────────────────────
# Load Data (cached)
# ────────────────────────────────────────────────────────────────
@st.cache_data
def load_cptac_data(cptac_dir):
    cptac_files = glob(os.path.join(cptac_dir, "*.tsv"))
    dfs = []
    for file in cptac_files:
        sample_name = os.path.splitext(os.path.basename(file))[0] + "_cptac"
        df = pd.read_csv(file, sep="\t", usecols=["gene_id", "tpm_unstranded"])
        df.rename(columns={"tpm_unstranded": sample_name}, inplace=True)
        df.set_index("gene_id", inplace=True)
        dfs.append(df)
    return pd.concat(dfs, axis=1)

@st.cache_data
def load_teresa_data(teresa_file):
    df = pd.read_csv(teresa_file, sep="\t", names=["gene_id", "teresa_gbm"], header=0)
    return df.set_index("gene_id")

@st.cache_data
def load_symbol_map(map_file):
    df_map = pd.read_csv(map_file, sep="\t", header=None, names=["ensg", "symbol"], dtype=str)
    sym2ensg = {s.upper(): e for e, s in df_map.values}
    ensg2sym = {e: s for e, s in df_map.values}
    return sym2ensg, ensg2sym

# ────────────────────────────────────────────────────────────────
# File Paths — Update these for your server!
# ────────────────────────────────────────────────────────────────
cptac_dir = "/Users/michael/cheng-project/cptac-gbm/rnaseq-files/clean.tpm.files"
teresa_file = "/Users/michael/cheng-project/boston-gene/rnaseq/core-rnaseq/star_rsem_gencode-v36/rsem.merged.gene_tpm.cleaned.tsv"
map_file = "/Users/michael/cheng-project/boston-gene/rnaseq/core-rnaseq/star_rsem_gencode-v36/gencode.v36.geneid.genename.map.txt"
# ────────────────────────────────────────────────────────────────
# Load and merge data
# ────────────────────────────────────────────────────────────────
cptac_df = load_cptac_data(cptac_dir)
teresa_df = load_teresa_data(teresa_file)
combined = pd.concat([cptac_df, teresa_df], axis=1)
combined_log2 = combined.map(lambda x: np.log2(x + 1))
sym2ensg, ensg2sym = load_symbol_map(map_file)

# ────────────────────────────────────────────────────────────────
# Plotting Function
# ────────────────────────────────────────────────────────────────
def plot_histograms(df, symbols, teresa_col="teresa_gbm"):
    ensg_list, title_list, missing = [], [], []
    for sym in symbols:
        key = sym.upper()
        if key in sym2ensg:
            ensg = sym2ensg[key]
            if ensg in df.index:
                ensg_list.append(ensg)
                title_list.append(sym)
            else:
                missing.append(f"{sym} (no TPM row)")
        else:
            missing.append(f"{sym} (no mapping)")

    if missing:
        st.warning("Missing or unmapped symbols: " + ", ".join(missing))
    if not ensg_list:
        return

    n = len(ensg_list)
    ncols = 4
    nrows = math.ceil(n / ncols)
    figsize = (ncols * 4, nrows * 3)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = axes.flatten()

    cohort = df.drop(columns=[teresa_col])
    teresa = df[teresa_col]

    for i, (ensg, sym) in enumerate(zip(ensg_list, title_list)):
        ax = axes[i]
        vals = cohort.loc[ensg].values
        tval = teresa.loc[ensg]
        percentile = (vals < tval).sum() / len(vals) * 100

        ax.hist(vals, bins=30, color="steelblue", edgecolor="black")
        ax.axvline(tval, color="crimson", linestyle="--", linewidth=2,
                   label=f"Teresa = {tval:.2f} ({percentile:.1f}%)")
        ax.set_title(f"{sym}")
        ax.set_xlabel("log2(TPM + 1)")
        ax.set_ylabel("CPTAC sample count")
        ax.legend()
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))  # 🔒 Force integer y-axis


    # Hide unused axes
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.tight_layout(h_pad=2.2)  # 🧱 Add space between rows
    st.pyplot(fig)  # ✅ Directly render the figure

# ────────────────────────────────────────────────────────────────
# UI Controls
# ────────────────────────────────────────────────────────────────
input_genes = st.text_input("Enter gene symbols (comma-separated):", "TP53,EGFR,CD274,IL6,IL10")
if st.button("Plot Histograms"):
    gene_list = [g.strip().upper() for g in input_genes.split(",") if g.strip()]
    plot_histograms(combined_log2, gene_list)
