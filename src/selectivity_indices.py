# selectivity_polars_per_fire.py

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

import yaml
import polars as pl


# ---------- Configuration helpers ----------

def load_config(filepath: Path) -> dict:
    with open(filepath, "r") as f:
        return yaml.safe_load(f)


# Load config relative to this file (adjust if needed)
SETTINGS = load_config(Path(__file__).parent.parent / "config.yml")

# Confirm keys exist
combined_path = SETTINGS["python"]["kunique_output"]  # input directory of parquet files
output_dir    = SETTINGS["python"]["kunique_selectivity_output"]       # output directory (we will force directory semantics)
print(f"combined_path={combined_path}")
print(f"output_dir={output_dir}")


# ---------- Per-fire selectivity (Polars) ----------

def compute_fire_selectivity_pl(
    parquet_path: str | Path,
    lc_col: str = "Lc_class",
    used_col: str = "Used",
) -> pl.DataFrame:
    """
    Compute burn probability, Chesson's alpha, and Jacobs' D per landcover class for one fire.

    Inputs:
      - parquet_path: path to a Parquet with columns [lc_col, used_col]
        where used_col ∈ {0,1}, and each row is a pixel within the availability domain.

    Returns (Polars DataFrame with one row per landcover class present in availability):
      [lc_col, n_available, n_burned, burn_prob, pi_avail, pi_use, chesson_alpha, jacobs_D]
    """
    parquet_path = Path(parquet_path)

    # Read only needed columns
    df = pl.read_parquet(parquet_path, columns=[lc_col, used_col, "Fire_ID", "K_UniqueID", "CLUSTERID", "AREA", "C_AREA", "P_num", "P_total", "lc_class_name"], use_pyarrow=False)

    # Keep valid rows
    df = df.filter(
        pl.col(lc_col).is_not_null() & pl.col(used_col).is_in([0, 1])
    )

    # Aggregate per landcover class
    agg = (
        df.group_by(lc_col)
          .agg(
              n_available = pl.len(), # number of pixels available in this class
              n_burned    = pl.col(used_col).sum(), # number of pixels burned in this class
              burn_prob   = pl.col(used_col).mean(), # burn probability, i.e., proportion of pixels burned in this class
              bu_area = pl.col(used_col).sum() * 30 * 30 / 10000.0, # burned area in hectares (assuming 30m pixels)
              avail_area = pl.len() * 30 * 30 / 10000.0, # available area in hectares
          )
          .sort(lc_col)
    )

    if agg.height == 0:
        # No valid rows; return empty with expected schema
        return pl.DataFrame(
            {
                lc_col: pl.Series([], dtype=df.schema.get(lc_col, pl.Int64)),
                "n_available": pl.Series([], dtype=pl.Int64),
                "n_burned": pl.Series([], dtype=pl.Int64),
                "burn_prob": pl.Series([], dtype=pl.Float64),
                "pi_avail": pl.Series([], dtype=pl.Float64),
                "pi_use": pl.Series([], dtype=pl.Float64),
                "chesson_alpha": pl.Series([], dtype=pl.Float64),
                "jacobs_D": pl.Series([], dtype=pl.Float64),
            }
        )

    # Totals across classes for this fire - including water
    A_total = int(df["P_total"][0]) # total available pixels 
    U_total = int(df["P_num"][0]) # total burned pixels 


    # get other important variables for metadata
    C_AREA = int(df["C_AREA"][0]) # total burned area in hectares
    AREA = int(df["AREA"][0]) # total area in hectares
    CLUSTERID = int(df["CLUSTERID"][0]) # cluster ID for this fire
    NFIREID = int(df["Fire_ID"][0]) # unique fire ID

    # for each lc_class, get the lc_class_name (assuming it's consistent within class)
    lc_names = (
        df.select(pl.col(lc_col), pl.col("lc_class_name"))
          .group_by(lc_col)
          .agg(pl.col("lc_class_name").first())
    )
    agg = agg.join(lc_names, on=lc_col, how="left") # add a column with the names of the Lc class


    # I need to remove the amount of pixels in the class "Water" (lc_class 15, 20, 19) from A_total, because those pixels are not available for burning.
    water_rows = agg.filter(pl.col(lc_col).is_in([15,20,19])) # get the row for Water class
    if water_rows.height > 0:
        water_available = int(water_rows.select(pl.col("n_available").sum()).item())
        A_total -= water_available
        print(f"[info] Adjusted A_total by removing Water classes: {water_available} pixels; new A_total={A_total}")
    
    # adjust burned pixels if any burned pixels were in the Water class (unlikely but for completeness)
    if water_rows.height > 0:
        water_burned = int(water_rows.select(pl.col("n_burned").sum()).item())
        U_total -= water_burned
        print(f"[info] Adjusted U_total by removing Water classes: {water_burned} pixels; new U_total={U_total}")

    # adjust A_total and U_total by removing agriculture classes (25,26,27) since those pixels are not part of the availability domain
    ag_rows = agg.filter(pl.col(lc_col).is_in([25, 26, 27]))
    if ag_rows.height > 0:
        ag_available = int(ag_rows.select(pl.col("n_available").sum()).item())
        ag_burned = int(ag_rows.select(pl.col("n_burned").sum()).item())
        A_total -= ag_available
        U_total -= ag_burned
        print(f"[info] Adjusted A_total by removing Agriculture classes: {ag_available} pixels; new A_total={A_total}")
        print(f"[info] Adjusted U_total by removing Agriculture classes: {ag_burned} pixels; new U_total={U_total}")

    # adjust A_total and U_total by removing urban (22, 23, 24)
    urban_rows = agg.filter(pl.col(lc_col).is_in([22, 23, 24]))
    if urban_rows.height > 0:
        urban_available = int(urban_rows.select(pl.col("n_available").sum()).item())
        urban_burned = int(urban_rows.select(pl.col("n_burned").sum()).item())
        A_total -= urban_available
        U_total -= urban_burned
        print(f"[info] Adjusted A_total by removing Urban classes: {urban_available} pixels; new A_total={A_total}")
        print(f"[info] Adjusted U_total by removing Urban classes: {urban_burned} pixels; new U_total={U_total}")

    # Remove agriculture, water and urban classes since those pixels are not part of the availability domain and would distort selectivity calculations.
    agg = agg.filter(~pl.col(lc_col).is_in([15, 19, 20, 22, 23, 24, 25, 26, 27]))

    # Add proportions across classes
    agg = agg.with_columns(
        pi_avail = pl.col("n_available") / pl.lit(A_total), # proportion of available pixels in this class
        pi_use   = pl.when(pl.lit(U_total) > 0) # proportion of burned pixels in this class, only if there are any burned pixels
                      .then(pl.col("n_burned") / pl.lit(U_total))
                      .otherwise(pl.lit(0.0))
    )

    # ----- Chesson's alpha -----
    # ratio = pi_use / pi_avail where pi_avail > 0
    agg = agg.with_columns(
        ratio = pl.when(pl.col("pi_avail") > 0)
                  .then(pl.col("pi_use") / pl.col("pi_avail"))
                  .otherwise(pl.lit(None))
    )
    denom = float(agg.select(pl.col("ratio").sum()).item()) if agg.height > 0 else float("nan")

    if denom and denom == denom:  # denom not 0 and not NaN
        agg = agg.with_columns(
            chesson_alpha = pl.col("ratio") / pl.lit(denom)
        )
    else:
        agg = agg.with_columns(
            chesson_alpha = pl.lit(None, dtype=pl.Float64)
        )

    # ----- Jacobs' D -----
    # D = (pi_use - pi_avail) / (pi_use + pi_avail - 2*pi_use*pi_avail)
    agg = agg.with_columns(
        _num = pl.col("pi_use") - pl.col("pi_avail"),
        _den = pl.col("pi_use") + pl.col("pi_avail") - 2.0 * pl.col("pi_use") * pl.col("pi_avail")
    ).with_columns(
        jacobs_D = pl.when(pl.col("_den") != 0.0)
                     .then(pl.col("_num") / pl.col("_den"))
                     .otherwise(pl.lit(None))
    ).drop(["ratio", "_num", "_den"])

    return agg, A_total, U_total, CLUSTERID, AREA, C_AREA, NFIREID


# ---------- Batch over a folder, save one Parquet per fire ----------

def ensure_out_dir(path_like: str | Path) -> Path:
    """
    Ensure output directory exists and return it as a Path.
    If a file path is provided, we still treat it as a directory name.
    """
    p = Path(path_like)
    # Force directory semantics
    p.mkdir(parents=True, exist_ok=True)
    return p



def detect_fire_id(
    pq_path: Path,
    prefer_column: bool = True,
    id_col: str = "K_UniqueID"
) -> str:
    """
    Determine K_UniqueID for the output filename.
    - If prefer_column=True and id_col exists in Parquet, use the unique value (if single unique).
    - Else, fallback to the file stem.
    """
    if prefer_column:
        try:
            # Read only the ID column if present
            df_id = pl.read_parquet(pq_path, columns=[id_col], use_pyarrow=False)
            if id_col in df_id.columns:
                uniques = df_id.select(pl.col(id_col).unique()).to_series().to_list()
                if len(uniques) == 1:
                    return str(int(uniques[0]))
        except Exception:
            pass
    # Fallback to filename without extension
    return pq_path.stem


def compute_and_save_all_fires(
    parquet_dir: str | Path,
    out_dir: str | Path,
    lc_col: str = "Lc_class",
    used_col: str = "Used",
    id_col: str = "K_UniqueID",
    prefer_id_column: bool = True,
    filename_pattern: str = "selectivity_{kunique_id}.parquet",
) -> None:
    """
    For each .parquet in parquet_dir:
      - compute selectivity metrics
      - save one Parquet per fire in out_dir using filename_pattern
    """
    parquet_dir = Path(parquet_dir)
    out_dir = ensure_out_dir(out_dir)

    files = sorted(parquet_dir.glob("*.parquet"))
    if not files:
        print(f"[warn] No parquet files found in: {parquet_dir}")
        return

    processed = 0
    for pq in files:
        kunique_id = detect_fire_id(pq, prefer_column=prefer_id_column, id_col=id_col)
        res, A_total, U_total, CLUSTERID, AREA, C_AREA, NFIREID = compute_fire_selectivity_pl(pq, lc_col=lc_col, used_col=used_col)

        if res.height == 0:
            print(f"[skip] No valid rows for fire_id={kunique_id} ({pq.name})")
            continue

        # Add fire_id column to saved file for traceability
        res = res.with_columns(pl.lit(kunique_id).alias("kunique_id"),
            pl.lit(A_total).alias("A_total"),
            pl.lit(U_total).alias("U_total"),
            pl.lit(CLUSTERID).alias("CLUSTERID"),
            pl.lit(AREA).alias("AREA"),
            pl.lit(C_AREA).alias("C_AREA"),
            pl.lit(NFIREID).alias("NFIREID")
        )

        res = res.select(
            ["kunique_id", lc_col, "lc_class_name", "n_available", "n_burned", "burn_prob", "bu_area", "avail_area", 
             "pi_avail", "pi_use", "chesson_alpha", "jacobs_D", "CLUSTERID", "AREA", "C_AREA", "NFIREID", "A_total", "U_total"]
        )

    

        out_file = out_dir / filename_pattern.format(kunique_id= kunique_id)
        res.write_parquet(out_file)
        processed += 1
        print(f"[ok] Wrote {out_file}")

    print(f"[summary] fires_processed={processed}, out_dir={out_dir}")


def main() -> int:
    in_dir = Path(combined_path)
    out_dir = Path(output_dir)

    if not in_dir.exists():
        print(f"[error] Input directory does not exist: {in_dir}", file=sys.stderr)
        return 2

    compute_and_save_all_fires(
        parquet_dir=in_dir,
        out_dir=out_dir,
        lc_col="Lc_class",
        used_col="Used",
        id_col="K_UniqueID",            # change if your ID column uses a different name
        prefer_id_column=True,       # set to False to always use filename stem
        filename_pattern="selectivity_{kunique_id}.parquet",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
