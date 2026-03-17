from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Tuple, Any, List

import pandas as pd
import numpy as np
import os


# Disable scientific notation globally
pd.options.display.float_format = '{:.10f}'.format


def add_jacobs_chesson(df):
    """
    Adds Jacobs and Chesson alpha electivity indices per lc_class
    within a single cluster (one parquet file).
    
    Input:
        df = output from compute_burn_and_available_ratios(df)
              → one row per lc_class
    """

    # ---- Jacobs Index ----
    df["jacobs"] = (df["bu"] - df["av"]) / (
        (df["bu"] + df["av"]) - (2 * df["bu"] * df["av"])
    ).replace(0, np.nan)

    # ---- Non-log Chesson α ----
    df["chesson_raw"] = (df["bu"] / df["av"]).replace([np.inf, -np.inf], np.nan)

    # Normalize within CLUSTERID (the file contains only ONE cluster)
    total_chesson = df["chesson_raw"].sum()

    df["Chesson"] = df["chesson_raw"] / total_chesson if total_chesson != 0 else np.nan

    # ---- Cleanup ----
    df.drop(columns="chesson_raw", inplace=True)

    # ---- Fixed formatting ----
    for col in ["jacobs", "Chesson"]:
        df[col] = df[col].apply(lambda x: float(f"{x:.10f}") if pd.notnull(x) else x)

    return df


def calculate_selectivity(in_dir, out_dir):
    """
    Adds 'lc_class_name' to every *.parquet file in directory
    based on the lc_class value and LC_MAP.
    """

    # Find all parquet files
    files = [f for f in os.listdir(in_dir) if f.endswith(".parquet")]

    for fname in files:
        fpath = os.path.join(out_dir, fname)

        # Load parquet
        df = pd.read_parquet(fpath)

        # compute ratios
        ratio_df = add_jacobs_chesson(df)

        # Save back (overwrite)
        ratio_df.to_parquet(fpath, index=False)

        print(f"Ratio calculated for: {fname}")

    print("All parquet files updated.")