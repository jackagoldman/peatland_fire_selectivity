from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Optional, Sequence, List, Tuple, Any

import numpy as np
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq


# ============================================================
# Filesystem discovery
# ============================================================

def list_fire_ids(out_dir: str | Path) -> List[str]:
    """
    Discover FireIDs under OUT_DIR structured as: OUT_DIR/FireID=<id>/
    Returns: ["<id1>", "<id2>", ...]
    """
    root = Path(out_dir)
    fire_ids: List[str] = []
    for fire_dir in root.glob("FireID=*"):
        if fire_dir.is_dir():
            fire_ids.append(fire_dir.name.split("=", 1)[-1])
    return sorted(set(fire_ids))


def find_cluster_parts_for_fire(out_dir: str | Path, fire_id: str) -> List[Path]:
    """
    Gather all 'part-*.parquet' under:
        OUT_DIR/FireID=<fire_id>/CLUSTERID=*/part-*.parquet
    """
    fire_root = Path(out_dir) / f"FireID={fire_id}"
    parts: List[Path] = []
    for cdir in sorted(fire_root.glob("CLUSTERID=*")):
        if not cdir.is_dir():
            continue
        parts.extend(sorted(cdir.glob("part-*.parquet")))
    return parts


# ============================================================
# Arrow helpers (robust to mixed row-group physical types)
# ============================================================

def _read_first_non_null(tbl: pa.Table, col: str) -> Optional[Any]:
    if col not in tbl.column_names:
        return None
    arr = tbl[col]
    if arr.null_count == arr.length():
        return None
    vals = arr.drop_null().to_pylist()
    return vals[0] if vals else None


def _promote_type(a: pa.DataType, b: pa.DataType) -> pa.DataType:
    """
    Pragmatic type promotion for schema union:
      * integers -> int64
      * integer + float / float + float -> float64
      * string/large_string -> large_string
      * boolean + boolean -> bool
      * timestamps (same tz) -> timestamp(us, tz)
      * anything else -> keep 'a'
    """
    if pa.types.is_integer(a) and pa.types.is_integer(b):
        return pa.int64()
    if (pa.types.is_integer(a) and pa.types.is_floating(b)) or \
       (pa.types.is_integer(b) and pa.types.is_floating(a)) or \
       (pa.types.is_floating(a) and pa.types.is_floating(b)):
        return pa.float64()
    if (pa.types.is_string(a) or pa.types.is_large_string(a)) and \
       (pa.types.is_string(b) or pa.types.is_large_string(b)):
        return pa.large_string()
    if pa.types.is_boolean(a) and pa.types.is_boolean(b):
        return pa.bool_()
    if pa.types.is_timestamp(a) and pa.types.is_timestamp(b):
        tz_a = getattr(a, "tz", None)
        tz_b = getattr(b, "tz", None)
        if tz_a == tz_b:
            return pa.timestamp("us", tz=tz_a)
    return a


def _normalize_table(tbl: pa.Table) -> pa.Table:
    """
    Decode dictionary columns to values; force CLUSTERID -> int64; combine chunks.
    """
    tbl = tbl.combine_chunks()
    cols: List[pa.Array] = []
    names: List[str] = []
    for i, field in enumerate(tbl.schema):
        col = tbl.column(i)
        if pa.types.is_dictionary(col.type):
            col = col.cast(col.type.value_type)
        if field.name == "CLUSTERID":
            col = col.cast(pa.int64())
        cols.append(col)
        names.append(field.name)
    return pa.Table.from_arrays(cols, names=names)


def _read_parquet_file_normalized(path: Path) -> pa.Table:
    """
    Read a Parquet file by row-group, normalize each RG, then concat.
    Avoids Arrow's internal schema merging errors when a single file
    has mixed row-group physical types (e.g., dict vs int).
    """
    pf = pq.ParquetFile(path)
    rg_tables: List[pa.Table] = []
    for rg_idx in range(pf.num_row_groups):
        rg_tbl = pf.read_row_group(rg_idx)
        rg_tbl = _normalize_table(rg_tbl)
        rg_tables.append(rg_tbl)
    if not rg_tables:
        return pa.Table.from_arrays([], names=[])
    return pa.concat_tables(rg_tables, promote=True)


def _union_canonical_schema(paths: Sequence[Path]) -> pa.Schema:
    """
    Scan first row-group (normalized) of each file to build a union schema.
    Ensures CLUSTERID is int64. Applies simple type promotions.
    """
    field_types: Dict[str, pa.DataType] = {}
    for p in paths:
        pf = pq.ParquetFile(p)
        if pf.num_row_groups == 0:
            continue
        sample = _normalize_table(pf.read_row_group(0))
        for field in sample.schema:
            name = field.name
            dtype = field.type
            if name == "CLUSTERID":
                dtype = pa.int64()
            if name in field_types:
                field_types[name] = _promote_type(field_types[name], dtype)
            else:
                field_types[name] = dtype

    # Order: CLUSTERID first if present; then sorted others
    names = list(field_types.keys())
    if "CLUSTERID" in names:
        names = ["CLUSTERID"] + sorted(n for n in names if n != "CLUSTERID")
    else:
        names = sorted(names)

    return pa.schema([(n, field_types[n]) for n in names])


def _align_to_schema(tbl: pa.Table, schema: pa.Schema) -> pa.Table:
    """
    Add missing columns (typed nulls) and cast/order to match target schema exactly.
    """
    aligned_cols: List[pa.Array] = []
    for field in schema:
        if field.name in tbl.column_names:
            col = tbl[field.name]
            if col.type != field.type:
                col = col.cast(field.type)
        else:
            col = pa.array([None] * tbl.num_rows, type=field.type)
        aligned_cols.append(col)
    return pa.Table.from_arrays(aligned_cols, names=[f.name for f in schema])


# ============================================================
# Output naming — from data column K_FireID (no zero padding)
# ============================================================

def _safe_id_string(x: Any) -> str:
    """
    Convert a value to a safe filename token: allow alnum, '-', '_'; replace others with '_'.
    """
    s = str(x)
    return "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in s)


def _determine_k_fire_id(parts: Sequence[Path]) -> str:
    """
    Read all parts and collect unique non-null values from K_FireID column.
    Require exactly one unique value; otherwise raise to prevent misnaming.
    """
    values: List[str] = []
    for p in parts:
        pf = pq.ParquetFile(p)
        for rg_idx in range(pf.num_row_groups):
            rg_tbl = _normalize_table(pf.read_row_group(rg_idx))
            v = _read_first_non_null(rg_tbl, "K_FireID")
            if v is not None:
                values.append(str(v))
            # Fast path: if we already saw >1 distinct, we can stop early
            if len(set(values)) > 1:
                break
        if len(set(values)) > 1:
            break

    uniq = sorted(set(values))
    if not uniq:
        raise ValueError("K_FireID column missing or all-null across shards; cannot name output.")
    if len(uniq) > 1:
        raise ValueError(f"Multiple K_FireID values found across shards: {uniq}; cannot pick a single output name.")
    return _safe_id_string(uniq[0])


def _format_output_filename_from_k_fire_id(k_fire_id: str) -> str:
    """
    Compose output filename as 'K_FireID_<value>.parquet' (no zero padding).
    """
    return f"K_FireID_{k_fire_id}.parquet"


# ============================================================
# Consolidation: one Parquet per FireID (all in ONE dest_dir)
# ============================================================

def consolidate_one_fire(
    out_dir: str | Path,
    fire_id: str,
    *,
    dest_dir: Optional[str | Path] = None,          # <-- outputs go directly here (flat)
    compression: str = "zstd",
    row_group_size: int = 512_000,
    drop_duplicates_on: Optional[Sequence[str]] = None,  # per-chunk
    factor_col: str = "area_factor_m2_per_pixel",
    validate_factor: bool = True,
    rtol: float = 1e-9,
    atol: float = 1e-12,
    delete_parts: bool = False,
) -> Optional[Path]:
    """
    Merge all cluster shards under a FireID into ONE Parquet file in `dest_dir`:
        dest_dir/K_FireID_<value>.parquet

    - Uses the *K_FireID column value* to name the file (no zero padding).
    - Robust to mixed row-group physical types (dictionary vs ints).
    - Forces CLUSTERID -> int64.
    """
    parts = find_cluster_parts_for_fire(out_dir, fire_id)
    if not parts:
        print(f"[consolidate] No parts found for FireID={fire_id}")
        return None

    # Destination: single directory for all outputs
    dest = Path(dest_dir) if dest_dir else Path(out_dir)
    dest.mkdir(parents=True, exist_ok=True)

    # Determine filename from K_FireID column
    k_fire_id = _determine_k_fire_id(parts)
    out_path = dest / _format_output_filename_from_k_fire_id(k_fire_id)

    if out_path.exists():
        print(f"[consolidate] Overwriting existing file: {out_path}")
        out_path.unlink()

    # Canonical schema across all parts for this FireID
    canonical_schema = _union_canonical_schema(parts)
    if "CLUSTERID" not in canonical_schema.names:
        raise KeyError(f"[consolidate] CLUSTERID not found in any shards for FireID={fire_id}")
    # Ensure K_FireID exists in the schema (we depend on it for naming; but keep even if absent later)
    if "K_FireID" not in canonical_schema.names:
        canonical_schema = canonical_schema.append(pa.field("K_FireID", pa.large_string()))

    # Optional factor consistency check
    if validate_factor and (factor_col in canonical_schema.names):
        vals = []
        for p in parts:
            try:
                t = _read_parquet_file_normalized(p)
                v = _read_first_non_null(t, factor_col)
                if v is not None:
                    vals.append(float(v))
            except Exception as e:
                print(f"[consolidate][WARN] Factor read failed for {p.name}: {e}")
        if len(vals) > 1:
            ref = float(vals[0])
            if not bool(np.allclose(vals, ref, rtol=rtol, atol=atol)):
                max_abs_diff = float(np.max(np.abs(np.array(vals) - ref)))
                print(f"[consolidate][WARN] area factor not identical for FireID={fire_id}; "
                      f"max abs diff={max_abs_diff:.3e}. Using per-row values as-is.")

    # Stream write per row-group
    with pq.ParquetWriter(out_path, schema=canonical_schema,
                          compression=compression, use_dictionary=True) as writer:
        for p in parts:
            pf = pq.ParquetFile(p)
            for rg_idx in range(pf.num_row_groups):
                rg_tbl = _normalize_table(pf.read_row_group(rg_idx))
                rg_tbl = _align_to_schema(rg_tbl, canonical_schema)

                # If K_FireID column is missing/null in this chunk, fill with the determined name
                if "K_FireID" in canonical_schema.names:
                    if "K_FireID" not in rg_tbl.column_names:
                        fill = pa.array([k_fire_id] * rg_tbl.num_rows, type=pa.large_string())
                        rg_tbl = rg_tbl.append_column("K_FireID", fill)
                        rg_tbl = _align_to_schema(rg_tbl, canonical_schema)
                    else:
                        # Fill nulls if any
                        col = rg_tbl["K_FireID"]
                        if col.null_count > 0:
                            as_list = [k_fire_id if x is None else str(x) for x in col.to_pylist()]
                            rg_tbl = rg_tbl.set_column(
                                rg_tbl.column_names.index("K_FireID"),
                                "K_FireID",
                                pa.array(as_list, type=pa.large_string()),
                            )

                if drop_duplicates_on:
                    df = pl.from_arrow(rg_tbl)
                    before = df.height
                    df = df.unique(subset=list(drop_duplicates_on), keep="first")
                    after = df.height
                    if after < before:
                        print(f"[consolidate] Dropped {before - after} dups in {p.name} rg#{rg_idx}")
                    rg_tbl = df.to_arrow()

                writer.write_table(rg_tbl, row_group_size=row_group_size)

    if delete_parts:
        for p in parts:
            try:
                p.unlink()
            except Exception as e:
                print(f"[consolidate][WARN] Could not delete {p}: {e}")

    print(f"[consolidate] ✓ Wrote {out_path} (from {len(parts)} file(s)).")
    return out_path


def consolidate_all_fires(
    out_dir: str | Path,
    *,
    dest_dir: Optional[str | Path] = None,          # <-- single directory for all outputs
    fires: Optional[Sequence[str]] = None,
    n: Optional[int] = None,
    random_sample: bool = False,
    **kwargs,  # passthrough to consolidate_one_fire
) -> List[Path]:
    """
    Consolidate shards into ONE file per FireID, all placed under `dest_dir`.
    File names are derived from the K_FireID column values.
    Returns a list of output Paths.
    """
    from random import sample, seed
    seed(42)

    all_ids = list_fire_ids(out_dir)
    if fires:
        allowed = set(map(str, fires))
        all_ids = [fid for fid in all_ids if fid in allowed]

    if n is not None:
        n = max(0, int(n))
        if random_sample and n < len(all_ids):
            all_ids = sample(all_ids, k=n)
        else:
            all_ids = all_ids[:n]

    if not all_ids:
        print("[consolidate] No matching FireIDs to process.")
        return []

    outputs: List[Path] = []
    for fid in all_ids:
        out = consolidate_one_fire(
            out_dir=out_dir,
            fire_id=fid,
            dest_dir=dest_dir,
            **kwargs
        )
        if out:
            outputs.append(out)

    return outputs


# ============================================================
# Concatenate all K_FireID_*.parquet (flat directory)
# ============================================================

def concatenate_fires(
    source_dir: str | Path,
    *,
    pattern: str = "K_FireID_*.parquet",     # <-- flat directory with per-fire outputs
    output_path: str | Path = "all_fires.parquet",
    compression: str = "zstd",
    row_group_size: int = 1_000_000,
) -> Optional[Path]:
    """
    Concatenate the per-Fire outputs found in a flat directory into a single Parquet file.

    - `source_dir` is globbed with `pattern` (default: K_FireID_*.parquet).
    - Assumes each file already contains a `K_FireID` column.
    """
    source_dir = Path(source_dir)
    files = sorted(source_dir.glob(pattern))
    if not files:
        print(f"[concat] No files found in {source_dir} with pattern '{pattern}'")
        return None

    # Build canonical schema across files
    canonical_schema = _union_canonical_schema(files)
    # Ensure K_FireID is in schema (should be, but guarantee it)
    if "K_FireID" not in canonical_schema.names:
        canonical_schema = canonical_schema.append(pa.field("K_FireID", pa.large_string()))

    out_path = Path(output_path)
    if out_path.exists():
        print(f"[concat] Overwriting existing file: {out_path}")
        out_path.unlink()

    with pq.ParquetWriter(out_path, schema=canonical_schema,
                          compression=compression, use_dictionary=True) as writer:
        for fp in files:
            t = _read_parquet_file_normalized(fp)
            t = _align_to_schema(t, canonical_schema)
            writer.write_table(t, row_group_size=row_group_size)

    print(f"[concat] ✓ Wrote {out_path} from {len(files)} file(s).")
    return out_path
