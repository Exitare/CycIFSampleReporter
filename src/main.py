import streamlit as st
import pandas as pd
import numpy as np
import time
from sklearn.preprocessing import StandardScaler
import plotly.graph_objects as go

st.title("CycIF Sample Reader")

bx1_file = "HTAN2_HMS_0000330944_bx1_OCT_tumor_cycif.csv"

fp = st.sidebar.file_uploader("Please select your CycIF data set")


def load_data():
    return pd.read_csv(bx1_file)


def violin_plot_cols(df):
    """ Violin plot for all columns in a data frame. """
    fig = go.Figure()

    for col in df.columns:
        fig.add_trace(go.Violin(y=df[col],
                                name=col))

    st.plotly_chart(fig)


raw_df = load_data()

raw_df['bx'] = 1
raw_df['scaled'] = 'N'

if st.checkbox('Show raw data'):
    st.subheader('Raw data')
    st.dataframe(raw_df)

# Scale markers for each dataset independently. Scale approach:
# 1. log10 transformation
# 2. standard scaling
# 3. clip to maximum range
scaler = StandardScaler()
scale_func = lambda d: scaler.fit_transform(d).clip(min=-10, max=10)

marker_cols = raw_df.filter(regex="Cell Masks$").filter(regex="^(?!(Goat|DAPI))").columns
raw_df_scaled = raw_df.copy()
raw_df_scaled[marker_cols] = scale_func(raw_df[marker_cols])
raw_df_scaled['scaled'] = 'Y'

# Combine unscaled and scaled data to create complete dataset.
cycif_df = pd.concat([raw_df, raw_df_scaled])
cycif_df.reset_index(inplace=True)

if st.checkbox('Show combined dataset'):
    st.write(cycif_df)

option = st.selectbox(
    'Scaled',
    ['Y', 'N'])

st.write(option)

violin_plot_cols(cycif_df[cycif_df.scaled == option][marker_cols])
