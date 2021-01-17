import streamlit as st
import pandas as pd
import numpy as np
import time
from sklearn.preprocessing import StandardScaler
import plotly.graph_objects as go
from services import file_loader
import scanpy as sc


def violin_plot_cols(df):
    """ Violin plot for all columns in a data frame. """
    fig = go.Figure()

    for col in df.columns:
        fig.add_trace(go.Violin(y=df[col],
                                name=col))

    st.plotly_chart(fig)


st.title("CycIF Sample Reader")

uploaded_file = st.sidebar.file_uploader("Please select your CycIF data set")

aand_file = None
# Load the file
if uploaded_file is not None:
    aand_file = sc.read_h5ad(uploaded_file)

# Only if a file is loaded show all the interactive data
if aand_file is not None:

    df = aand_file.to_df()
    st.subheader("Raw data set")

    if st.checkbox('Show raw data set'):
        st.write(df)

    violin_plot_cols(df)
