import pandas as pd
import scanpy as sc


def load_file(uploaded_file):
    return sc.read_h5ad(uploaded_file)
