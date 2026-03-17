import pandas as pd
import os

pd.options.display.float_format = '{:.10f}'.format

# water classes (from LC_MAP)
WATER_CLASSES = {15, 19, 20}

def compute_burn_and_available_ratios(df):
    """
    df must contain:
        - lc_class (int)
        - used (1=burned, 0=unburned)
        - area_factor_m2_per_pixel (float)
        - K_FireID, CLUSTERID, K_uniqueID, DATE (metadata columns)
    """

    # ---- 0. Extract metadata ----
    meta_cols = ["K_FireID", "CLUSTERID", "K_UniqueID", "DATE", "lc_class_name"]
    meta = df[meta_cols].iloc[0].to_dict()

    # ---- 1. Remove all water pixels BEFORE any calculations ----
    df = df[~df["lc_class"].isin(WATER_CLASSES)].copy()

    # if everything is water (rare), avoid crashes
    if df.empty:
        raise ValueError("All pixels were water — cannot compute ratios.")

    # ---- 2. Compute weighted area per lc_class ----
    grouped = (
        df.groupby("lc_class")
        .apply(lambda g: pd.Series({
            "burned_area": ((g["used"] == 1).astype(int) * g["area_factor_m2_per_pixel"]).sum(),
            "available_area": g["area_factor_m2_per_pixel"].sum()
        }))
        .reset_index()
    )

    # ---- 3. Totals across classes ----
    total_burned = grouped["burned_area"].sum()
    total_available = grouped["available_area"].sum()

    # ---- 4. Ratios ----
    grouped["bu"] = grouped["burned_area"] / total_burned if total_burned > 0 else 0
    grouped["av"] = grouped["available_area"] / total_available if total_available > 0 else 0

    # ---- 5. Add totals ----
    grouped["bu_area"] = float(f"{total_burned:.10f}")
    grouped["av_area"] = float(f"{total_available:.10f}")


    # ---- Compute burn probability ----
    grouped["bp"] = grouped.apply(
        lambda row: float(f"{(row['burned_area'] / row['available_area']) if row['available_area'] > 0 else 0:.10f}"),
        axis=1
    )


    # ---- 6. Add metadata ----
    for col, val in meta.items():
        grouped[col] = val

    # ---- 7. Format floats (no scientific notation) ----
    float_cols = ["burned_area", "available_area", "bu", "av"]
    for c in float_cols:
        grouped[c] = grouped[c].apply(lambda x: float(f"{x:.10f}"))

    # ---- 8. Move metadata to front ----
    ordered_cols = meta_cols + [c for c in grouped.columns if c not in meta_cols]
    grouped = grouped[ordered_cols]

    return grouped

 # run for all   

def calculate_ratios(in_dir, out_dir):
    """
    Adds 'lc_class_name' to every *.parquet file in directory
    based on the lc_class value and LC_MAP.
    """

    # Find all parquet files
    files = [f for f in os.listdir(in_dir) if f.endswith(".parquet")]

    for fname in files:
        fpath = os.path.join(in_dir, fname)

        # Load parquet
        df = pd.read_parquet(fpath)

        # compute ratios
        ratio_df = compute_burn_and_available_ratios(df)

        # Save back (overwrite)
        out_path = os.path.join(out_dir, fname)
        ratio_df.to_parquet(out_path, index=False)

        print(f"Ratio calculated for: {fname}")

    print("All parquet files updated.")