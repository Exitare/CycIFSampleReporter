import streamlit as st
import pandas as pd
import numpy as np
import time
from sklearn.preprocessing import StandardScaler
import plotly.graph_objects as go
from services import file_loader
import scanpy as sc


def violin_plot_cols(df, selected_columns):
    """ Violin plot for all columns in a data frame. """
    fig = go.Figure()
    if len(selected_columns) == 0:
        for col in df.columns:
            fig.add_trace(go.Violin(y=df[col],
                                    name=col))
    else:
        for col in selected_columns:
            fig.add_trace(go.Violin(y=df[col],
                                    name=col))

    st.plotly_chart(fig)


st.title("CycIF Sample Reader")

# File uploader in a sidebar
uploaded_file = st.sidebar.file_uploader("Please select your CycIF data set")

aand_file = None
# Load the file or display explanation
if uploaded_file is not None:
    aand_file = sc.read_h5ad(uploaded_file)
else:
    st.subheader('FAQ')
    st.write("Supported file types:")
    st.write("- h5ad")

    st.subheader('How to')
    st.write('- Select a h5ad file using the file uploader.')
    st.write('- After you selected a file the data will be loaded automatically '
             'and you will be able to explore your data interactively.')

# Only if a file is loaded show all the interactive data
if aand_file is not None:

    df = aand_file.to_df()
    st.subheader("Raw data set")

    if st.checkbox('Show raw data set'):
        st.write(df)

    selected_columns = st.multiselect('Select columns you want to display:',
                                      df.columns)

    violin_plot_cols(df, selected_columns)
