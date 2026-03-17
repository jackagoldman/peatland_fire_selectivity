from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Tuple, Any, List

import polars as pl


# ---------- Per-fire selectivity (Polars) ----------

def _read_unique_value(df: pl.DataFrame, col: str) -> Optional[Any]:
    if col not in df.columns:
        return None
    vals = df.select(pl.col(col).drop_nulls().unique()).to_series().to_list()
    if len(vals) == 1:
        return vals[0]
    return None


def compute_fire_selectivity_pl(
    parquet_path: str | Path,
    lc_col: str = "Lc_class",
    used_col: str = "Used",
    *,
    factor_col: str = "area_factor_m2_per_pixel",   # if present, used for areas
    pixel_size_m: float = 30.0,                      # fallback pixel size if factor_col absent
    exclude_classes: Sequence[int] = (15, 19, 20, 22, 23, 24, 25, 26, 27),  # water+urban+ag
    class_name_col: str = "lc_class_name",
    id_col_prefer: str = "K_FireID",                # for naming/metadata
) -> Tuple[pl.DataFrame, int, int, dict]:
    """
    Compute burn probability, Chesson's alpha, and Jacobs' D per landcover class for one fire.

    Inputs:
      - parquet_path: path to a Parquet with at least [lc_col, used_col], where used_col ∈ {0,1},
                      each row is a pixel within the availability domain (before exclusions).

    Returns:
      (result_df, A_total_pixels, U_total_pixels, meta_dict)

      result_df columns:
        [lc_col, class_name_col?, n_available, n_burned, burn_prob,
         avail_area_ha, bu_area_ha, pi_avail, pi_use, chesson_alpha, jacobs_D,
         plus optional passthrough metadata columns if present]
    """
    parquet_path = Path(parquet_path)

    # Build column list safely
    base_cols = [lc_col, used_col]
    optional_cols = [class_name_col, "CLUSTERID", "AREA", "C_AREA", "P_num", "P_total",
                     "Fire_ID", "K_UniqueID", "K_FireID", factor_col]
    cols_to_read = [c for c in base_cols + optional_cols if c is not None]

    # Read needed columns
    df = pl.read_parquet(parquet_path, columns=[c for c in cols_to_read if c], use_pyarrow=False)

    # sanity: used is 0/1; keep valid rows only
    df = df.filter(
        pl.col(lc_col).is_not_null() & pl.col(used_col).is_in([0, 1])
    )

    # Identify whether a per-pixel area factor is present
    has_factor = factor_col in df.columns

    # Compute per-row areas (ha) if factor available; else use constant pixel size
    if has_factor:
        # area_factor_m2_per_pixel already in m²
        df = df.with_columns(
            avail_area_ha_row = (pl.col(factor_col) / 10_000.0).cast(pl.Float64),
            bu_area_ha_row    = (pl.col(used_col) * pl.col(factor_col) / 10_000.0).cast(pl.Float64),
        )
    else:
        px_ha = (pixel_size_m * pixel_size_m) / 10_000.0
        df = df.with_columns(
            avail_area_ha_row = pl.lit(px_ha).alias("avail_area_ha_row"),
            bu_area_ha_row    = (pl.col(used_col) * px_ha).alias("bu_area_ha_row"),
        )

    # Aggregate per landcover class
    agg = (
        df.group_by(lc_col)
          .agg(
              n_available = pl.len().cast(pl.Int64),
              n_burned    = pl.col(used_col).sum().cast(pl.Int64),
              burn_prob   = pl.col(used_col).mean().cast(pl.Float64),
              avail_area_ha = pl.col("avail_area_ha_row").sum().cast(pl.Float64),
              bu_area_ha    = pl.col("bu_area_ha_row").sum().cast(pl.Float64),
          )
          .sort(lc_col)
    )

    # Early exit if no data
    if agg.height == 0:
        meta = {}
        # If available, carry some metadata (safely)
        for mcol in ("CLUSTERID", "AREA", "C_AREA", "P_num", "P_total", "Fire_ID", "K_UniqueID", "K_FireID"):
            val = _read_unique_value(df, mcol)
            if val is not None:
                meta[mcol] = val
        return agg, 0, 0, meta

    # Attach class names if present
    if class_name_col in df.columns:
        lc_names = (
            df.select(pl.col(lc_col), pl.col(class_name_col))
              .group_by(lc_col)
              .agg(pl.col(class_name_col).first())
        )
        agg = agg.join(lc_names, on=lc_col, how="left")

    # Compute totals BEFORE removing excluded classes (we will adjust below)
    # Then remove excluded classes for selectivity
    excluded = set(exclude_classes)
    excl_mask = pl.col(lc_col).is_in(list(excluded))

    # Get totals (pixels) *including* excluded to derive deltas
    totals_including = agg.select(
        pl.col("n_available").sum().alias("A_all"),
        pl.col("n_burned").sum().alias("U_all")
    )

    # Adjust totals by removing excluded classes
    totals_excluded = agg.filter(excl_mask).select(
        pl.col("n_available").sum().alias("A_exc"),
        pl.col("n_burned").sum().alias("U_exc")
    )

    A_all = int(totals_including["A_all"][0]) if totals_including.height else 0
    U_all = int(totals_including["U_all"][0]) if totals_including.height else 0
    A_exc = int(totals_excluded["A_exc"][0]) if totals_excluded.height else 0
    U_exc = int(totals_excluded["U_exc"][0]) if totals_excluded.height else 0

    # Final totals used for selectivity:
    A_total = max(A_all - A_exc, 0)
    U_total = max(U_all - U_exc, 0)

    # Remove excluded classes from the table for selectivity
    agg = agg.filter(~excl_mask)

    # Add proportions across remaining classes
    agg = agg.with_columns(
        pi_avail = (pl.col("n_available") / pl.lit(A_total)).cast(pl.Float64),
        pi_use   = pl.when(pl.lit(U_total) > 0)
                      .then((pl.col("n_burned") / pl.lit(U_total)).cast(pl.Float64))
                      .otherwise(pl.lit(0.0))
    )

    # ----- Chesson's alpha -----
    # ratio = pi_use / pi_avail where pi_avail > 0
    agg = agg.with_columns(
        ratio = pl.when(pl.col("pi_avail") > 0.0)
                  .then(pl.col("pi_use") / pl.col("pi_avail"))
                  .otherwise(pl.lit(None))
    )
    denom = float(agg.select(pl.col("ratio").sum()).item()) if agg.height > 0 else float("nan")

    if denom and denom == denom:  # not 0 and not NaN
        agg = agg.with_columns(
            chesson_alpha = (pl.col("ratio") / pl.lit(denom)).cast(pl.Float64)
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
                     .then((pl.col("_num") / pl.col("_den")).cast(pl.Float64))
                     .otherwise(pl.lit(None))
    ).drop(["ratio", "_num", "_den"])

    # Build metadata (prefer K_FireID for naming)
    meta: dict = {}
    for mcol in ("CLUSTERID", "AREA", "C_AREA", "P_num", "P_total", "Fire_ID", "K_UniqueID", "K_FireID"):
        val = _read_unique_value(df, mcol)
        if val is not None:
            meta[mcol] = val

    # Reorder/select final columns (keep class name if present)
    final_cols = [lc_col]
    if class_name_col in agg.columns:
        final_cols.append(class_name_col)
    final_cols += [
        "n_available", "n_burned", "burn_prob",
        "avail_area_ha", "bu_area_ha",
        "pi_avail", "pi_use", "chesson_alpha", "jacobs_D",
    ]
    # Keep metadata if present in agg (they usually aren't in agg; meta stays separate)
    agg = agg.select([c for c in final_cols if c in agg.columns])

    return agg, A_total, U_total, meta


# ---------- Batch over a folder, save one Parquet per fire ----------

def ensure_out_dir(path_like: str | Path) -> Path:
    """
    Ensure output directory exists and return it as a Path.
    """
    p = Path(path_like)
    p.mkdir(parents=True, exist_ok=True)
    return p


def _filename_from_k_fire_id(k_fire_id: Any) -> str:
    safe = "".join(ch if str(ch).isalnum() or ch in ("-", "_") else "_" for ch in str(k_fire_id))
    return f"selectivity_{safe}.parquet"


def compute_and_save_all_fires(
    parquet_dir: str | Path,
    out_dir: str | Path,
    lc_col: str = "Lc_class",
    used_col: str = "Used",
    *,
    factor_col: str = "area_factor_m2_per_pixel",
    pixel_size_m: float = 30.0,
    exclude_classes: Sequence[int] = (15, 19, 20, 22, 23, 24, 25, 26, 27),
    class_name_col: str = "lc_class_name",
    id_col_prefer: str = "K_FireID",               # read from data, not filename
) -> None:
    """
    For each K_FireID_*.parquet in `parquet_dir`:
      - compute selectivity metrics
      - save one Parquet per fire in `out_dir` named from the K_FireID column
    """
    parquet_dir = Path(parquet_dir)
    out_dir = ensure_out_dir(out_dir)

    files = sorted(parquet_dir.glob("K_FireID_*.parquet"))
    if not files:
        print(f"[warn] No per-fire parquet files found in: {parquet_dir}")
        return

    processed = 0
    for pq in files:
        # Compute
        res, A_total, U_total, meta = compute_fire_selectivity_pl(
            pq,
            lc_col=lc_col,
            used_col=used_col,
            factor_col=factor_col,
            pixel_size_m=pixel_size_m,
            exclude_classes=exclude_classes,
            class_name_col=class_name_col,
            id_col_prefer=id_col_prefer,
        )

        if res.height == 0:
            print(f"[skip] No valid rows for file: {pq.name}")
            continue

        # Determine K_FireID from data
        k_fire_id = meta.get("K_FireID", None)
        if k_fire_id is None:
            # Fallback: try reading the file with polars to get it
            try:
                tmp = pl.read_parquet(pq, columns=[id_col_prefer], use_pyarrow=False)
                vals = tmp.select(pl.col(id_col_prefer).drop_nulls().unique()).to_series().to_list()
                if len(vals) == 1:
                    k_fire_id = vals[0]
            except Exception:
                pass
        if k_fire_id is None:
            # Final fallback: parse from filename pattern (K_FireID_<value>.parquet)
            stem = pq.stem  # "K_FireID_<value>"
            k_fire_id = stem.replace("K_FireID_", "")

        # Attach metadata columns to result
        # If some metadata missing, compute reasonable fallbacks
        res = res.with_columns(
            pl.lit(k_fire_id).alias("K_FireID"),
            pl.lit(int(A_total)).alias("A_total_pixels"),
            pl.lit(int(U_total)).alias("U_total_pixels"),
        )
        for m in ("CLUSTERID", "AREA", "C_AREA", "P_num", "P_total", "Fire_ID", "K_UniqueID"):
            if m in meta:
                res = res.with_columns(pl.lit(meta[m]).alias(m))

        # Reorder output columns
        out_cols = ["K_FireID", lc_col]
        if class_name_col in res.columns:
            out_cols.append(class_name_col)
        out_cols += [
            "n_available", "n_burned", "burn_prob",
            "avail_area_ha", "bu_area_ha",
            "pi_avail", "pi_use", "chesson_alpha", "jacobs_D",
            "A_total_pixels", "U_total_pixels",
        ]
        # Optional metadata if present
        for m in ("CLUSTERID", "AREA", "C_AREA", "P_num", "P_total", "Fire_ID", "K_UniqueID"):
            if m in res.columns:
                out_cols.append(m)

        res = res.select([c for c in out_cols if c in res.columns])

        # Save using K_FireID from the data
        out_file = out_dir / _filename_from_k_fire_id(k_fire_id)
        res.write_parquet(out_file)
        processed += 1
        print(f"[ok] Wrote {out_file}")

    print(f"[summary] fires_processed={processed}, out_dir={out_dir}")
