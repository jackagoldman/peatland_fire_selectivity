# consolidate_parquet.py
from __future__ import annotations
import os
import glob
import math
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple, Dict

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def _list_clusters(out_dir: str | Path) -> Dict[str, Sequence[str]]:
    """
    Discover clusters under OUT_DIR structured as:
      OUT_DIR/FireID=<id>/CLUSTERID=<id>/
    Returns {fire_id: [cluster_id, ...], ...}
    """
    out_dir = Path(out_dir)
    fires: Dict[str, list] = {}
    for fire_dir in out_dir.glob("FireID=*"):
        if not fire_dir.is_dir():
            continue
        fire_id = fire_dir.name.split("=", 1)[-1]
        clusters = []
        for cdir in fire_dir.glob("CLUSTERID=*"):
            if cdir.is_dir():
                clusters.append(cdir.name.split("=", 1)[-1])
        if clusters:
            fires[fire_id] = sorted(set(clusters))
    return fires


def _get_shards_for_cluster(out_dir: str | Path, fire_id: str, cluster_id: str) -> Tuple[Path, Sequence[Path]]:
    """
    Return (cluster_dir, [part paths]) for a given FireID/CLUSTERID.
    """
    cluster_dir = Path(out_dir) / f"FireID={fire_id}" / f"CLUSTERID={cluster_id}"
    parts = sorted(cluster_dir.glob("part-*.parquet"))
    return cluster_dir, parts


def _read_area_factor_from_table(tbl: pa.Table, col: str = "area_factor_m2_per_pixel") -> Optional[float]:
    if col not in tbl.column_names:
        return None
    arr = tbl[col]
    if arr.null_count == arr.length():
        return None
    values = np.array(arr.drop_nulls().to_pylist(), dtype=float)
    if values.size == 0:
        return None
    return float(values[0])


def _validate_factor_values(factors: Sequence[float], rtol: float = 1e-9, atol: float = 1e-12) -> Tuple[bool, float, float]:
    """
    Ensure all provided factors are effectively identical within tolerance.
    Returns (ok, ref, max_abs_diff)
    """
    vals = np.array([f for f in factors if f is not None], dtype=float)
    if vals.size == 0:
        return (False, math.nan, math.nan)
    ref = float(vals[0])
    diffs = np.abs(vals - ref)
    ok = bool(np.allclose(vals, ref, rtol=rtol, atol=atol))
    return (ok, ref, float(diffs.max(initial=0.0)))


def consolidate_one_cluster(
    out_dir: str | Path,
    fire_id: str,
    cluster_id: str,
    *,
    output_name: str = "data.parquet",
    compression: str = "zstd",
    row_group_size: int = 512_000,
    factor_col: str = "area_factor_m2_per_pixel",
    validate_factor: bool = True,
    rtol: float = 1e-9,
    atol: float = 1e-12,
    delete_parts: bool = False,
    drop_duplicates_on: Optional[Sequence[str]] = None,
    dest_base_dir: Optional[str | Path] = None,   # <-- NEW: write final file to a different root
) -> Optional[Path]:
    """
    Merge all part-*.parquet shards for (FireID, CLUSTERID) into a single Parquet file.

    - Writes final to:
        dest_base_dir/FireID=<id>/CLUSTERID=<id>/output_name
      if dest_base_dir is provided; otherwise next to shards.
    - Ensures `area_factor_m2_per_pixel` is unique (within tolerance) across shards.
    - Optionally deletes part files after consolidation.
    - Optionally drops duplicates on specified subset of columns (e.g., ["x","y","lc_class","used"]).

    Returns the path to the consolidated Parquet, or None if no shards found.
    """
    cluster_dir, parts = _get_shards_for_cluster(out_dir, fire_id, cluster_id)
    if not parts:
        print(f"[consolidate] No parts found for FireID={fire_id}, CLUSTERID={cluster_id}")
        return None

    # Determine destination directory (mirror structure under dest_base_dir if provided)
    if dest_base_dir is None:
        dest_dir = cluster_dir
    else:
        dest_dir = Path(dest_base_dir) / f"FireID={fire_id}" / f"CLUSTERID={cluster_id}"
        dest_dir.mkdir(parents=True, exist_ok=True)

    out_path = dest_dir / output_name
    if out_path.exists():
        print(f"[consolidate] Already exists: {out_path}. Overwriting.")
        out_path.unlink()

    # Prepare ParquetWriter with schema from the first shard
    first_table = pq.read_table(parts[0])
    schema = first_table.schema

    # Validate factor across shards (optional)
    if validate_factor and (factor_col in schema.names):
        factors = []
        # include factor from first shard
        factors.append(_read_area_factor_from_table(first_table, factor_col))
        # check remaining shards cheaply without loading to pandas
        for p in parts[1:]:
            t = pq.read_table(p, columns=[factor_col])
            factors.append(_read_area_factor_from_table(t, factor_col))
        ok, ref, max_abs_diff = _validate_factor_values(factors, rtol=rtol, atol=atol)
        if not ok:
            print(f"[consolidate][WARN] area factor not identical across shards "
                  f"(FireID={fire_id}, CLUSTERID={cluster_id}); max abs diff={max_abs_diff:.3e}. "
                  f"Using the value present in each row as-is.")

    # Write consolidated file using ParquetWriter (streaming tables)
    with pq.ParquetWriter(out_path, schema=schema, compression=compression, use_dictionary=True) as writer:
        for p in parts:
            tbl = pq.read_table(p)  # read this shard as Arrow Table
            if drop_duplicates_on:
                pdf = tbl.to_pandas(types_mapper=pd.ArrowDtype)
                before = len(pdf)
                pdf.drop_duplicates(subset=list(drop_duplicates_on), inplace=True)
                after = len(pdf)
                if after < before:
                    print(f"[consolidate] Dropped {before - after} duplicate rows in {p.name}")
                tbl = pa.Table.from_pandas(pdf, schema=schema, preserve_index=False)
            writer.write_table(tbl, row_group_size=row_group_size)

    # Optionally delete parts in the ORIGINAL shard location
    if delete_parts:
        for p in parts:
            try:
                p.unlink()
            except Exception as e:
                print(f"[consolidate][WARN] Could not delete {p}: {e}")

    print(f"[consolidate] ✓ Wrote {out_path} (from {len(parts)} shard(s)).")
    return out_path


def consolidate_all_clusters(
    out_dir: str | Path,
    *,
    fires: Optional[Sequence[str]] = None,
    n: Optional[int] = None,                 # limit number of fires for testing
    random_sample: bool = False,
    delete_parts: bool = False,
    dest_base_dir: Optional[str | Path] = None,  # <-- NEW: new root for consolidated outputs
    **kwargs,                                # passthrough to consolidate_one_cluster
) -> None:
    """
    Consolidate shards for every (FireID, CLUSTERID) found under OUT_DIR.

    Args:
        out_dir: base directory (contains FireID=<id>/CLUSTERID=<id>/)
        fires: if provided, restrict to this set of FireIDs (strings)
        n: if provided, limit to the first n fires (or random n if random_sample=True)
        random_sample: sample fires randomly when n is specified
        delete_parts: delete part-*.parquet after consolidation (in the ORIGINAL shard directory)
        dest_base_dir: write final consolidated parquet(s) under this new root, preserving structure
        **kwargs: forwarded to consolidate_one_cluster (e.g., row_group_size, compression)
    """
    from random import sample, seed
    seed(42)

    mapping = _list_clusters(out_dir)
    if not mapping:
        print("[consolidate] No clusters found under out_dir.")
        return

    # Restrict fires if requested
    fire_ids = sorted(mapping.keys())
    if fires:
        allowed = set(map(str, fires))
        fire_ids = [fid for fid in fire_ids if fid in allowed]
    if n is not None:
        n = max(0, int(n))
        if random_sample:
            k = min(n, len(fire_ids))
            fire_ids = sample(fire_ids, k=k)
        else:
            fire_ids = fire_ids[:n]

    if not fire_ids:
        print("[consolidate] No matching FireIDs to process.")
        return

    # Consolidate
    for fid in fire_ids:
        cluster_ids = mapping[fid]
        for cid in cluster_ids:
            consolidate_one_cluster(
                out_dir=out_dir,
                fire_id=fid,
                cluster_id=cid,
                delete_parts=delete_parts,
                dest_base_dir=dest_base_dir,  # <-- pass through
                **kwargs
            )