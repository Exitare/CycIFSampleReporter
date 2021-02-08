import plotly.graph_objects as go
import seaborn
import streamlit as st
import seaborn as sns

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


def create_bar_plot_pre_post(melt):
    g = sns.catplot(
        data=melt, kind="bar",
        x="marker", y="value", hue="stage",
        ci="sd", palette="dark", alpha=.6, height=6
    )
    g.despine(left=True)
    g.set_axis_labels("", "Expression (Mean)")
    g.legend.set_title("")
    g.set_xticklabels(rotation=90)
    return g
