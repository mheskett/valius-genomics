import streamlit as st
import pandas as pd
import yaml

with open("config.yaml") as f:
    config = yaml.safe_load(f)

somatic_table = config["somatic_mutations"]["summary_table"]

# Set Streamlit page config
st.set_page_config(page_title="Somatic Mutations Table Viewer", layout="wide")
st.title("Somatic Mutations Table Viewer")

# File path to your summary table
#file_path = "/Users/michael/cheng-project/boston-gene/exome/core-sarek/results/annotation/tumor_wes_vs_normal_wes.mutect2.filtered_snpEff.ann.table.pass.txt"

# Load and cache the data
@st.cache_data
def load_table(path):
    return pd.read_csv(path, sep="\t")

# Load data and display
df = load_table(somatic_table)
st.dataframe(df, use_container_width=True)
