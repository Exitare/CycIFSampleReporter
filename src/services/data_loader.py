import pandas as pd
import scanpy as sc
import streamlit as st


@st.cache
def load_file(uploaded_file):
    return sc.read_h5ad(uploaded_file)


@st.cache
def extract_data(aand_file):
    df = aand_file.to_df()
    communities = aand_file.obs['PhenoGraph_clusters']
    communities = communities.astype('object')
    Q = aand_file.uns['PhenoGraph_Q']
    k = aand_file.uns['PhenoGraph_k']

    return df, communities, Q, k


@st.cache
def create_means_df(communities, pg_clusters, df):
    means = []
    for i in range(0, pd.Series(communities).max() + 1):
        cells_index = pg_clusters[pg_clusters == i].index
        filtered_markers_df = df[df.index.isin(cells_index)]
        means.append(filtered_markers_df.mean().values)

    return pd.DataFrame(means, columns=df.columns)
