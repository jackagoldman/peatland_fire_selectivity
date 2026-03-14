from __future__ import annotations
import os
from pathlib import Path
from typing import Optional, Sequence, Dict, List

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def _list_fires(out_dir: str | Path) -> List[str]:
    """
    Discover FireID directories under out_dir:
      out_dir/FireID=<id>/
    Returns list of fire_id strings.
    """
    out_dir = Path(out_dir)
    fires = []
    for fire_dir in out_dir.glob("FireID=*"):
        if fire_dir.is_dir():
            fire_id = fire_dir.name.split("=", 1)[-1]
            fires.append(fire_id)
    return sorted(fires)


def _list_cluster_dirs_for_fire(out_dir: str | Path, fire_id: str) -> List[Path]:
    """
    Returns all CLUSTERID directories for a given FireID:
      out_dir/FireID=<fire_id>/CLUSTERID=<cluster_id>/
    """
    base = Path(out_dir) / f"FireID={fire_id}"
    return sorted([p for p in base.glob("CLUSTERID=*") if p.is_dir()])


def _list_shards_for_cluster(cluster_dir: Path) -> List[Path]:
    return sorted(cluster_dir.glob("part-*.parquet"))


def consolidate_one_fire(
    out_dir: str | Path,
    fire_id: str,
    *,
    dest_base_dir: Optional[str | Path] = None,
    output_name: Optional[str] = None,
    compression: str = "zstd",
    row_group_size: int = 1_000_000,
    drop_duplicates_on: Optional[Sequence[str]] = None,
    delete_parts: bool = False,
    require_nonempty: bool = True,
) -> Optional[Path]:
    """
    Merge all cluster shards for a single FireID into ONE Parquet file.

    Input layout (shards produced by your worker):
      out_dir/FireID=<id>/CLUSTERID=<cluster_id>/part-*.parquet

    Output:
      dest_base_dir/FireID=<id>.parquet      (if dest_base_dir is provided)
      or out_dir/FireID=<id>.parquet         (next to the FireID folder)

    Args:
        out_dir: root with FireID=.../CLUSTERID=.../part-*.parquet
        fire_id: the FireID to consolidate
        dest_base_dir: write final FireID file under this root if provided
        output_name: filename (default: FireID=<id>.parquet)
        compression: parquet compression (e.g., 'zstd','snappy','gzip')
        row_group_size: target rows per row group
        drop_duplicates_on: optional subset of columns to de-duplicate across shards
        delete_parts: if True, delete shards after successful write
        require_nonempty: if True, return None if no shards found

    Returns:
        Path to the consolidated Parquet, or None if nothing to consolidate.
    """
    fire_root = Path(out_dir) / f"FireID={fire_id}"
    if not fire_root.exists():
        if require_nonempty:
            print(f"[fire] No directory for FireID={fire_id}")
        return None

    cluster_dirs = _list_cluster_dirs_for_fire(out_dir, fire_id)
    parts: List[Path] = []
    for cdir in cluster_dirs:
        parts.extend(_list_shards_for_cluster(cdir))

    if not parts:
        if require_nonempty:
            print(f"[fire] No shards found for FireID={fire_id}")
        return None

    # Destination file path
    if output_name is None:
        output_name = f"FireID={fire_id}.parquet"
    if dest_base_dir:
        dest_base_dir = Path(dest_base_dir)
        dest_base_dir.mkdir(parents=True, exist_ok=True)
        out_path = dest_base_dir / output_name
    else:
        out_path = Path(out_dir) / output_name

    if out_path.exists():
        print(f"[fire] Overwriting existing file: {out_path}")
        out_path.unlink()

    # Use schema from first shard
    first_tbl = pq.read_table(parts[0])
    schema = first_tbl.schema

    # Stream-write all shards into one file
    with pq.ParquetWriter(out_path, schema=schema, compression=compression, use_dictionary=True) as writer:
        for p in parts:
            tbl = pq.read_table(p)
            if drop_duplicates_on:
                # de-dupe within this shard before writing
                pdf = tbl.to_pandas(types_mapper=pd.ArrowDtype)
                before = len(pdf)
                pdf.drop_duplicates(subset=list(drop_duplicates_on), inplace=True)
                after = len(pdf)
                if after < before:
                    print(f"[fire] Dropped {before-after} duplicate rows in {p.name}")
                tbl = pa.Table.from_pandas(pdf, schema=schema, preserve_index=False)
            writer.write_table(tbl, row_group_size=row_group_size)

    # Optionally delete shards after success
    if delete_parts:
        for p in parts:
            try:
                p.unlink()
            except Exception as e:
                print(f"[fire][WARN] Could not delete {p}: {e}")

    print(f"[fire] ✓ Wrote {out_path} from {len(parts)} shard(s).")
    return out_path


def consolidate_all_fires(
    out_dir: str | Path,
    *,
    fires: Optional[Sequence[str]] = None,
    n: Optional[int] = None,
    random_sample: bool = False,
    dest_base_dir: Optional[str | Path] = None,
    compression: str = "zstd",
    row_group_size: int = 1_000_000,
    drop_duplicates_on: Optional[Sequence[str]] = None,
    delete_parts: bool = False,
) -> None:
    """
    Consolidate all FireIDs under out_dir into one Parquet per FireID.

    Args:
        out_dir: root with FireID=<id>/CLUSTERID=<id>/part-*.parquet
        fires: explicit list of FireIDs to consolidate (strings)
        n: limit to first N fires (or random sample if random_sample=True)
        random_sample: sample N fires randomly when n is specified
        dest_base_dir: where to write FireID files (keep shards in out_dir)
        compression: parquet compression
        row_group_size: row group size for writer
        drop_duplicates_on: optional columns for de-duplication
        delete_parts: if True, delete shards after writing the FireID file
    """
    from random import seed, sample
    seed(42)

    all_fires = _list_fires(out_dir)
    if not all_fires:
        print("[fire] No FireID directories found.")
        return

    if fires is not None:
        allow = set(map(str, fires))
        fire_ids = [fid for fid in all_fires if fid in allow]
    else:
        fire_ids = all_fires

    if n is not None:
        n = max(0, int(n))
        if random_sample:
            fire_ids = sample(fire_ids, k=min(n, len(fire_ids)))
        else:
            fire_ids = fire_ids[:n]

    if not fire_ids:
        print("[fire] No FireIDs selected for consolidation.")
        return

    for fid in fire_ids:
        consolidate_one_fire(
            out_dir=out_dir,
            fire_id=fid,
            dest_base_dir=dest_base_dir,
            compression=compression,
            row_group_size=row_group_size,
            drop_duplicates_on=drop_duplicates_on,
            delete_parts=delete_parts,
        )


