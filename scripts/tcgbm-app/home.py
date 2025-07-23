import streamlit as st

st.set_page_config(page_title="BioData Explorer", layout="wide")
# Display logo
st.image("horizontal.png", width=200)  # Adjust width as needed
st.title("Welcome to the Valius BioData Explorer")
st.markdown("""
This is a multi-tool web app for exploring GBM expression data.

Use the sidebar to navigate:
- View gene expression histograms
- Explore PCA and dimensionality reduction
- Upload your own datasets
""")


st.markdown(
    """
    <hr style="margin-top: 4em;">
    <div style="text-align: center; font-size: 0.8em; color: gray;">
        © 2025 Valius Sciences — All rights reserved.<br>
        Michael Heskett, Ph. D.
    </div>
    """,
    unsafe_allow_html=True
)