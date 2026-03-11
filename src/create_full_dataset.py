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

