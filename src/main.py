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
    high_low_marker_df = pd.DataFrame(columns=["Index", "Cell Name (High)", "High", "Cell Name (Low)", "Low"])
    for i, community in means_df.iterrows():
        highest = community.sort_values(ascending=False).head(3)
        lowest = community.sort_values(ascending=False).tail(3)

        for j in range(0, 3):
            high_low_marker_df = high_low_marker_df.append({
                "Index": i,
                "Cell Name (High)": highest.index[j],
                "High": highest[j],
                "Cell Name (Low)": lowest.index[j],
                "Low": lowest[j],

            }, ignore_index=True)

            # st.write(high_low_marker_df)

    threshold = st.slider('Please select your threshold:', math.floor(high_low_marker_df["Low"].min()),
                          math.ceil(high_low_marker_df["High"].max()), 0)
    st.write(threshold)
    if threshold == 0:
        st.dataframe(high_low_marker_df)
    else:
        low_df = high_low_marker_df[high_low_marker_df["Low"] > threshold]
        high_df = high_low_marker_df[high_low_marker_df["High"].astype('float') > threshold]
        st.write(high_df)
        st.write(low_df)
        temp_df = pd.concat(high_df, low_df,  join="inner")
        st.dataframe(temp_df)

    # st.write(i, "***")
    # st.write(community.sort_values(ascending=False).head(3))
    # st.write(community.sort_values(ascending=False).tail(3))
