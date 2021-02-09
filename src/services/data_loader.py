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
    spatial_df = pd.DataFrame(columns=['X_centroid', 'Y_centroid'])
    spatial_df['X_centroid'] = aand_file.obs['X_centroid']
    spatial_df['Y_centroid'] = aand_file.obs['Y_centroid']
    return df, communities, Q, k, spatial_df


@st.cache
def group_biopsy_data(biopsy_data):
    """
    Returns the pre and post treatment indices for each biopsy
    """
    pre_treatment = biopsy_data[biopsy_data['biopsy'] == 1]
    post_treatment = biopsy_data[biopsy_data['biopsy'] == 2]
    return pre_treatment, post_treatment


@st.cache
def add_pre_column(df, pre_index, post_index):
    data = df.copy()
    data.at[pre_index.index, 'stage'] = 'Pre'
    data.at[post_index.index, 'stage'] = 'Post'
    return data


@st.cache
def create_means_df(communities, pg_clusters, df):
    means = []
    for i in range(0, pd.Series(communities).max() + 1):
        cells_index = pg_clusters[pg_clusters == i].index
        filtered_markers_df = df[df.index.isin(cells_index)]
        means.append(filtered_markers_df.mean().values)

    return pd.DataFrame(means, columns=df.columns)


@st.cache
def create_high_low_marker_df(means_df):
    high_low_marker_df = pd.DataFrame(columns=["Community", "Cell Name (High)", "High", "Cell Name (Low)", "Low"])
    for i, community in means_df.iterrows():
        highest = community.sort_values(ascending=False).head(3)
        lowest = community.sort_values(ascending=False).tail(3)

        for j in range(0, 3):
            high_low_marker_df = high_low_marker_df.append({
                "Community": i,
                "Cell Name (High)": highest.index[j],
                "High": highest[j],
                "Cell Name (Low)": lowest.index[j],
                "Low": lowest[j],

            }, ignore_index=True)

    return high_low_marker_df


@st.cache
def melt_data_for_pre_post_bar(data):
    means = data.groupby('stage').mean().T
    means.reset_index(level=0, inplace=True)
    means.rename(columns={"index": "marker"}, inplace=True)
    return pd.melt(means, id_vars=['marker'], value_vars=['Pre', 'Post'])
