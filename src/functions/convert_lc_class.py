import pandas as pd
import os

# Land cover mapping dictionary
LC_MAP = {
    1: "Open Upland",
    2: "Open Poor Fen",
    3: "Open Rich Fen",
    4: "Open PPC",
    5: "Open Bog",
    6: "Open Mineral Wetland",
    7: "Forested Upland",
    8: "Forested PPC",
    9: "Forested Mineral Wetland",
    10: "Treed Mineral Wetland",
    11: "Treed Upland",
    12: "Forested Rich Fen",
    13: "Forested Bog",
    14: "Forested Poor Fen",
    15: "Water",
    16: "Treed PPC",
    17: "Treed Rich Fen",
    18: "Treed Bog",
    19: "Water",
    20: "Water",
    21: "Treed Poor Fen",
    22: "Treed Urban",
    23: "Open Urban",
    24: "Forested Urban",
    25: "Forested Agriculture",
    26: "Open Agriculture",
    27: "Treed Agriculture",
}

def add_lc_class_name_to_parquet(directory):
    """
    Adds 'lc_class_name' to every *.parquet file in directory
    based on the lc_class value and LC_MAP.
    """

    # Find all parquet files
    files = [f for f in os.listdir(directory) if f.endswith(".parquet")]

    for fname in files:
        fpath = os.path.join(directory, fname)

        # Load parquet
        df = pd.read_parquet(fpath)

        # Add name column (unknown classes become NaN)
        df["lc_class_name"] = df["lc_class"].map(LC_MAP)

        # Save back (overwrite)
        df.to_parquet(fpath, index=False)

        print(f"Updated: {fname}")

    print("All parquet files updated.")