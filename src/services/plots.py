import plotly.graph_objects as go
import seaborn
import streamlit as st


@st.cache
def create_violin_plot(df, selected_columns):
    """ Violin plot for all columns in a data frame. """
    fig = go.Figure()
    if len(selected_columns) == 0:
        return fig
    else:
        for col in selected_columns:
            fig.add_trace(go.Violin(y=df[col],
                                    name=col))

    return fig


@st.cache
def create_clustermap(means_df):
    return seaborn.clustermap(means_df)
