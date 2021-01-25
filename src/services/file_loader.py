import pandas as pd
import scanpy as sc
import streamlit as st


@st.cache(allow_output_mutation=True)
def load_file(uploaded_file):
    return sc.read_h5ad(uploaded_file)
