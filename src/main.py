import streamlit as st
import pandas as pd
import phenograph
from services import plots, data_loader
import seaborn
import plotly.express as px
import numpy as np
import scanpy as sc
import math

st.title("CycIF Sample Reader")

# File uploader in a sidebar
uploaded_file = st.sidebar.file_uploader("Please select your CycIF data set")

aand_file = None
# Load the file or display explanation
if uploaded_file is None:
    st.subheader('FAQ')
    st.write("Supported file types:")
    st.write("- h5ad")

    st.subheader('How to')
    st.write('- Select a h5ad file using the file uploader.')
    st.write('- After you selected a file the data will be loaded automatically '
             'and you will be able to explore your data interactively.')


else:
    aand_file = data_loader.load_file(uploaded_file)
    df, communities, Q, k = data_loader.extract_data(aand_file)
    pg_clusters = pd.Series(communities)

    means_df = data_loader.create_means_df(communities, pg_clusters, df)

    st.subheader("Raw data set")

    if st.checkbox('Show raw data set'):
        st.write(df)

    selected_columns = st.multiselect('Select columns you want to display:',
                                      df.columns)

    st.plotly_chart(plots.create_violin_plot(df, selected_columns))
    # fig = plots.create_clustermap(means_df)
    # st.plotly_chart(fig)
    # st.plotly_chart(px.imshow(means_df, x=df.columns))

    # fig = plots.create_clustermap(means_df)
    # fig.show()

    high_low_marker_df = data_loader.create_high_low_marker_df(means_df)

    st.subheader("Mean markers per community")
    st.text("This section provides possibilites to discover the highest and \nlowest mean marker values for each phenograph community")

    threshold = st.slider('Please select the threshold:', math.floor(high_low_marker_df["Low"].min()),
                          math.ceil(high_low_marker_df["High"].max()), 0)
    community_selector = st.slider('Please select the community you want to take a closer look:',
                                   math.floor(high_low_marker_df["Community"].min()),
                                   math.ceil(high_low_marker_df["Community"].max()), -1)

    if threshold == 0 and community_selector == -1:
        st.dataframe(high_low_marker_df)
    elif threshold != 0 and community_selector == -1:
        low_df = high_low_marker_df[
            high_low_marker_df["Low"] > threshold]
        high_df = high_low_marker_df[high_low_marker_df["High"].astype('float') > threshold]
        temp_df = pd.concat([low_df, high_df]).drop_duplicates().reset_index(drop=True)
        st.dataframe(temp_df)
    elif threshold == 0 and community_selector != -1:
        temp_df = high_low_marker_df[high_low_marker_df["Community"] == community_selector]
        st.dataframe(temp_df)
    else:
        low_df = high_low_marker_df[
            (high_low_marker_df["Low"] > threshold) & (high_low_marker_df["Community"] == community_selector)]
        high_df = high_low_marker_df[(high_low_marker_df["High"].astype('float') > threshold) & (high_low_marker_df[
                                                                                                     "Community"] == community_selector)]
        temp_df = pd.concat([low_df, high_df]).drop_duplicates().reset_index(drop=True)
        st.dataframe(temp_df)
