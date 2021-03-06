import streamlit as st
import pandas as pd
import phenograph
from services import plots, data_loader
import seaborn
import plotly.express as px
import numpy as np
import scanpy as sc
import math
import matplotlib.pyplot as plt
import seaborn as sns

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
    st.write('- After selecting a file the data will be loaded automatically '
             'and one will be able to explore your data interactively.')


else:
    aand_file = data_loader.load_file(uploaded_file)
    df, communities, Q, k, spatial_df = data_loader.extract_data(aand_file)
    pg_clusters = pd.Series(communities)
    biopsy_data = pd.Series(aand_file.obs['biopsy']).to_frame()
    means_df = data_loader.create_means_df(communities, pg_clusters, df)

    st.subheader("Raw data set")

    if st.checkbox('Show raw data set'):
        st.write(df)

    st.subheader('General information')
    st.write(f'The dataset contains {communities.nunique()} communities'
             f' and {len(df.columns)} cells with {len(df.index)} rows.')
    st.write("Proportions of clusters:")

    st.subheader("Biopsy Comparison")
    # pre_index, post_index = data_loader.group_biopsy_data(biopsy_data)
    # data = data_loader.add_pre_column(df, pre_index, post_index)
    # melt = data_loader.melt_data_for_pre_post_bar(data)
    # st.pyplot(plots.create_bar_plot_pre_post(melt))

    # df_with_communities = df.copy()
    # df_with_communities['Community'] = communities

    st.subheader("Charts")
    col1, col2 = st.beta_columns(2)
    col1.header("Bar Chart")

    col1.bar_chart(communities.value_counts())
    col2.header("Area chart:")
    col2.area_chart(communities.value_counts())

    # Cell exploration
    st.subheader("Explore markers")
    st.write('Select cells to have a closer look at their distribution and to compare them with others.')

    selected_columns = st.multiselect('Select columns you want to display:',
                                      df.columns)
    st.plotly_chart(plots.create_violin_plot(df, selected_columns))

    temp_means = means_df.copy()
    temp_means['index'] = means_df.index
    test_melt = pd.melt(temp_means, id_vars=['index'], var_name='marker')
    test_melt.rename(columns={"index": "community"}, inplace=True)
    mean_community_selector = st.slider(
        'Select marker: ',
        math.floor(test_melt["community"].min()),
        math.ceil(test_melt["community"].max()), -1)

    if mean_community_selector == -1:
        g = sns.catplot(
            data=test_melt, kind="bar",
            x="marker", y="value", hue="community",
            ci="sd", palette="dark", alpha=.6, height=6
        )
        g.despine(left=True)
        g.set_axis_labels("", "Expression (Mean)")
        g.set_xticklabels(rotation=90)
        st.pyplot(g)
    else:
        temp = test_melt[test_melt['community'] == mean_community_selector]
        g = sns.catplot(
            data=temp, kind="bar",
            x="marker", y="value", hue="community",
            ci="sd", palette="dark", alpha=.6, height=6
        )
        g.despine(left=True)
        g.set_axis_labels("", "Expression (Mean)")
        g.set_xticklabels(rotation=90)
        st.pyplot(g)

    st.write(df)
    st.write(communities)
    spatial_df['communities'] = communities
    st.write(spatial_df)

    fig = px.scatter(
        x=spatial_df["X_centroid"],
        y=spatial_df["Y_centroid"],
        color=spatial_df['communities']
    )
    st.write(fig)

    high_low_marker_df = data_loader.create_high_low_marker_df(means_df)

    st.subheader("Mean markers per community")
    st.text(
        "This section provides possibilites to discover the highest and \nlowest mean marker values for each phenograph community")

    threshold = st.slider('Please select the threshold:', math.floor(high_low_marker_df["Low"].min()),
                          math.ceil(high_low_marker_df["High"].max()), 0)
    community_selector = st.slider(
        'Please select the community you want to take a closer look: (-1 selects all communities)',
        math.floor(high_low_marker_df["Community"].min()),
        math.ceil(high_low_marker_df["Community"].max()), -1)

    # TODO: Add function for that
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

    st.subheader("Clustered heatmap of phenograph communities.")

    fig = px.imshow(means_df, x=df.columns)
    st.plotly_chart(fig)
