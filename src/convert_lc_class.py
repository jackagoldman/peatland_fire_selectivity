#!/usr/bin/env python3
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


def resolve_paths(base_dir: Path, settings: dict) -> tuple[Path, Path]:
    in_rel = settings["python"]["kunique_input"]       # e.g., "data/ua_fire_combined/"
    out_rel = settings["python"]["kunique_output"]  # e.g., "data/full_dataset/"
    in_abs = (base_dir / in_rel).resolve()
    out_abs = (base_dir / out_rel).resolve()
    out_abs.mkdir(parents=True, exist_ok=True)
    return in_abs, out_abs


# ---------- LC mapping expression (streaming-friendly) ----------

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


def build_lc_expr(colname: str) -> pl.Expr:
    """
    Return a Polars expression that maps integer LC codes to names.
    Uses `replace` (vectorized) which is streaming-compatible.
    """
    # Ensure values are integers before replacing (non-integers become null)
    return (
        pl.col(colname)
        .cast(pl.Int64, strict=False)
        .replace(LC_MAP, default=None)
        .alias("lc_class_name")
    )


def find_lc_col(schema_names: list[str], preferred: Optional[str]) -> Optional[str]:
    """
    Choose the LC column name. Priority:
      1) explicit name from config (if present),
      2) first match among common variants (case-sensitive check).
    """
    if preferred and preferred in schema_names:
        return preferred

    candidates = ["lc_class", "Lc_class", "lc", "LC"]
    for c in candidates:
        if c in schema_names:
            return c
    return None


# ---------- Main per-file conversion ----------

def convert_one_file(in_file: Path, out_root: Path, lc_col_hint: Optional[str] = None) -> None:
    # Build output path: mirror relative directory, same filename under out_root
    # Example: in: data/ua_fire_combined/part-0001.snappy.parquet
    #          out: data/full_dataset/part-0001.snappy.parquet
    rel = in_file.relative_to(in_file.parents[0]) if in_file.parents else in_file.name
    out_file = out_root / rel  # same file name at the output root
    # Insert `_lc` before the final suffix 
    out_file = out_file.with_name(out_file.stem + "_lc" + out_file.suffix)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    # Lazy scan (no load into memory)
    lf = pl.scan_parquet(str(in_file))

    # Decide column name from the *schema* (no compute yet)
    schema_names = lf.collect_schema().names()
    lc_col = find_lc_col(schema_names, lc_col_hint)
    if lc_col is None:
        print(f"[WARN] {in_file.name}: no LC column found in {schema_names}; copying as-is")
        # Just sink the file unchanged
        lf.sink_parquet(str(out_file))
        return

    # Add mapped name column with streaming-friendly expression
    lf2 = lf.with_columns(build_lc_expr(lc_col))

    # Stream-write
    lf2.sink_parquet(
        str(out_file),
        # options you can tweak:
        compression="zstd",  # good modern default
        statistics=True,
    )

    print(f"[OK]   {in_file.name} → {out_file}")


def main() -> int:
    print("[S0] polars LC conversion (per-file, streaming)")

    # Resolve repo base and config (works in `uv run` or python)
    base_dir = Path(__file__).resolve().parent.parent if "__file__" in globals() else Path.cwd().parent
    config_path = base_dir / "config.yml"
    print(f"[S1] config: {config_path}")

    settings = load_config(config_path)
    in_dir, out_dir = resolve_paths(base_dir, settings)
    print(f"[S2] input:  {in_dir}")
    print(f"[S2] output: {out_dir}")

    # Optional: allow LC column override via config (python.lc_column)
    lc_col_hint = settings.get("python", {}).get("lc_column")

    # Enumerate parquet files
    files = sorted(in_dir.rglob("*.parquet"))
    if not files:
        print(f"[ERR] No .parquet files under {in_dir}")
        return 2

    # Process each file independently
    failures = 0
    for i, f in enumerate(files, 1):
        try:
            print(f"[S3] ({i}/{len(files)}) {f.name}")
            convert_one_file(f, out_dir, lc_col_hint=lc_col_hint)
        except Exception as e:
            failures += 1
            print(f"[FAIL] {f}: {e}", file=sys.stderr)

    print(f"[DONE] processed={len(files)} failures={failures}")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())