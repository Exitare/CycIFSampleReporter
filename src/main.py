import streamlit as st
import pandas as pd
import phenograph
from services import plots, file_loader
import seaborn
import plotly.express as px


def clustering():
    communities, graph, Q = phenograph.cluster(aand_file.X, k=aand_file.uns['PhenoGraph_k'])
    return communities, graph, Q


def preprocess_data(df):
    pass


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
    aand_file = file_loader.load_file(uploaded_file)
    df = aand_file.to_df()

    communities, graph, Q = clustering()

    pg_clusters = pd.Series(communities)
    means = []
    for i in range(0, pd.Series(communities).max() + 1):
        cells_index = pg_clusters[pg_clusters == i].index
        filtered_markers_df = df[df.index.isin(cells_index)]
        means.append(filtered_markers_df.mean().values)


    st.write(means)
    means_df = pd.DataFrame(means, columns=df.columns)

    st.subheader("Raw data set")

    if st.checkbox('Show raw data set'):
        st.write(df)

    selected_columns = st.multiselect('Select columns you want to display:',
                                      df.columns)

    st.plotly_chart(plots.create_violin_plot(df, selected_columns))

    st.write(means_df)
    st.plotly_chart(plots.create_clustermap(means_df))
    st.plotly_chart(px.imshow(means_df, x=df.columns))

    # fig = plots.create_clustermap(means_df)
    # fig.show()
