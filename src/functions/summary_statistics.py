import polars as pl
from typing import Sequence

def summarize_lands_abundance(parquet_path: str,
                             area_col: str = "available_area",
                             burned_col: str = "burned_area",
                             class_col: str = "lc_class",
                             region_cols: Sequence[str] = ("ecozone",),
                             save_prefix: str | None = None) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Summarizes the abundance of land cover classes in terms of total available area and burned area.
    Parameters    ----------
    parquet_path : str
        Path to the input parquet file containing the data.
    area_col : str, optional
        Name of the column representing the available area, by default "available_area"
    burned_col : str, optional
        Name of the column representing the burned area, by default "burned_area"
    class_col : str, optional
        Name of the column representing the land cover class, by default "lc_class"
    region_cols : Sequence[str], optional
        List of columns representing the region (e.g., "ecozone"), by default ("ecozone",)
    save_prefix : str | None, optional
        If provided, the prefix for saving the summarized parquet file, by default None
    Returns    -------
    tuple[pl.DataFrame, pl.DataFrame]
        A tuple containing two Polars DataFrames: the first with the summarized abundance of land cover classes per class, and the second with the region-level summary.
    """
    # ---- 1. Load data ----
    df = pl.read_parquet(parquet_path)

    # ---- 2. Group by region and class, then sum areas ----
    grp = list(region_cols) + [class_col]
    per_class = (
        df.groupby(grp)
        .agg([
            pl.sum(burned_col).alias("class_burned"), 
            pl.sum(area_col).alias("class_area"),
            pl.count().alias("n_class")
        ])
    )

    # ---- 3. Compute totals per region ----
    per_class = per_class.with_columns(
        (pl.col("class_area") / pl.col("class_area").sum().over(list(region_cols))).alias("abundance")
    )

    # region-level summary of burn_frequency across classes
    region_summary = (
        per_class.groupby(list(region_cols))
                 .agg([
                     pl.count("lc_class").alias("n_classes"),
                     pl.col("bp").mean().alias("burn_freq_mean"),
                     pl.col("bp").median().alias("burn_freq_median"),
                     pl.col("bp").std().alias("burn_freq_sd"),
                     pl.col("bp").quantile(0.25).alias("burn_freq_q25"),
                     pl.col("bp").quantile(0.75).alias("burn_freq_q75"),
                     ( (pl.col("class_area") * pl.col("bp")).sum() / pl.col("class_area").sum() ).alias("burn_freq_weighted_mean"),
                     pl.col("abundance").median().alias("abundance_median")
                 ])
    )


    if save_prefix:
        per_class.write_parquet(f"{save_prefix}_per_class.parquet")
        region_summary.write_parquet(f"{save_prefix}_region_summary.parquet")

    return per_class, region_summary

    # ----- function to find fires where peatlands burned at the same rate as uplands -----
def find_peatland_fires(per_class_df: pl.DataFrame, peatland_classes: Sequence[int]) -> pl.DataFrame:
    """Identifies fires where peatland classes burned at the same rate as upland classes.
    Parameters    ----------
    per_class_df : pl.DataFrame
        DataFrame containing the summarized abundance of land cover classes per class, including burn probabilities.
    peatland_classes : Sequence[int]
        List of land cover class IDs that correspond to peatlands.
    Returns    -------
    pl.DataFrame
        A DataFrame containing the subset of fires where peatland classes burned at the same rate as upland classes.
    """
    # Filter for peatland and upland classes
    peatland = per_class_df.filter(pl.col("lc_class").is_in(peatland_classes))
    upland = per_class_df.filter(~pl.col("lc_class").is_in(peatland_classes))

    # Compute average burn probability for peatlands and uplands per fire
    peatland_burn_prob = peatland.groupby("fire_id").agg(pl.mean("bp").alias("peat_burn_prob"))
    upland_burn_prob = upland.groupby("fire_id").agg(pl.mean("bp").alias("upland_burn_prob"))

    # Join the two summaries on fire_id
    burn_comparison = peatland_burn_prob.join(upland_burn_prob, on="fire_id")

    # Identify fires where burn probabilities are approximately equal (within a certain threshold)
    threshold = 0.05  # Example threshold for "same rate"
    similar_burns = burn_comparison.filter(
        (pl.col("peat_burn_prob") - pl.col("upland_burn_prob")).abs() <= threshold
    )

    return similar_burns

    # ----- function to find fires where peatlands burned at a higher rate than uplands -----
def find_peatland_dominant_fires(per_class_df: pl.DataFrame, peatland_classes: Sequence[int]) -> pl.DataFrame:
    """Identifies fires where peatland classes burned at a higher rate than upland classes.
    Parameters    ----------
    per_class_df : pl.DataFrame
        DataFrame containing the summarized abundance of land cover classes per class, including burn probabilities.
    peatland_classes : Sequence[int]
        List of land cover class IDs that correspond to peatlands.
    Returns    -------
    pl.DataFrame
        A DataFrame containing the subset of fires where peatland classes burned at a higher rate than upland classes.
    """
    # Filter for peatland and upland classes
    peatland = per_class_df.filter(pl.col("lc_class").is_in(peatland_classes))
    upland = per_class_df.filter(~pl.col("lc_class").is_in(peatland_classes))

    # Compute average burn probability for peatlands and uplands per fire
    peatland_burn_prob = peatland.groupby("fire_id").agg(pl.mean("bp").alias("peat_burn_prob"))
    upland_burn_prob = upland.groupby("fire_id").agg(pl.mean("bp").alias("upland_burn_prob"))

    # Join the two summaries on fire_id
    burn_comparison = peatland_burn_prob.join(upland_burn_prob, on="fire_id")

    # Identify fires where peatland burn probability is significantly higher than upland burn probability
    threshold = 0.05  # Example threshold for "higher rate"
    dominant_peat_fires = burn_comparison.filter(
        (pl.col("peat_burn_prob") - pl.col("upland_burn_prob")) > threshold
    )

    return dominant_peat_fires


    # ------ function to find preference, avoidance and indifference fires for Jacobs D -----
    # Strong Preference : D > 0.5
    # Moderate Preference : 0.25 < D <= 0.5
    # Indifference : -0.25 <= D <= 0.25
    # Moderate Avoidance : -0.5 <= D < -0.25
    # Strong Avoidance : D < -0.5

def classify_fire_preference(per_class_df: pl.DataFrame, peatland_classes: Sequence[int]) -> pl.DataFrame:
    """Classifies fires based on the preference of peatland classes compared to upland classes using Jacobs D.
    Parameters    ----------
    per_class_df : pl.DataFrame
        DataFrame containing the summarized abundance of land cover classes per class, including burn probabilities.
    peatland_classes : Sequence[int]
        List of land cover class IDs that correspond to peatlands.
    Returns    -------
    pl.DataFrame
        A DataFrame containing the classification of fires based on Jacobs D.
    """

    # Sum burned area and available area for peatland and upland classes per fire
    peat_agg = (
        per_class_df
        .filter(pl.col("lc_class").is_in(peatland_classes))
        .groupby("fire_id")
        .agg([
            pl.col("class_burned").sum().alias("peat_burned"),
            pl.col("class_area").sum().alias("peat_area")
        ])
    )

    upland_agg = (
        per_class_df
        .filter(~pl.col("lc_class").is_in(peatland_classes))
        .groupby("fire_id")
        .agg([
            pl.col("class_burned").sum().alias("upland_burned"),
            pl.col("class_area").sum().alias("upland_area")
        ])
    )

    # Join peatland and upland aggregates; missing sides become 0
    burn_comparison = peat_agg.join(upland_agg, on="fire_id", how="outer")
    burn_comparison = burn_comparison.with_columns([
        pl.col("peat_burned").fill_null(0.0),
        pl.col("peat_area").fill_null(0.0),
        pl.col("upland_burned").fill_null(0.0),
        pl.col("upland_area").fill_null(0.0)
    ])

    # compute proportions: r = proportion of burned area that is peat; p = proportion of available area that is peat
    burn_comparison = burn_comparison.with_columns([
        ((pl.col("peat_burned") / (pl.col("peat_burned") + pl.col("upland_burned"))).fill_null(0.0)).alias("r"),
        ((pl.col("peat_area") / (pl.col("peat_area") + pl.col("upland_area"))).fill_null(0.0)).alias("p")
    ])

    # Jacobs D: (r - p) / (r + p - 2*r*p); handle zero denominator
    den = (pl.col("r") + pl.col("p") - 2 * pl.col("r") * pl.col("p"))
    burn_comparison = burn_comparison.with_columns(
        pl.when(den != 0).then((pl.col("r") - pl.col("p")) / den).otherwise(pl.lit(None)).alias("jacobs")
    )

    # Classify fires based on Jacobs D
    burn_comparison = burn_comparison.with_columns(
        pl.when(pl.col("jacobs") > 0.5).then("Strong Preference")
          .when((pl.col("jacobs") > 0.25) & (pl.col("jacobs") <= 0.5)).then("Moderate Preference")
          .when((pl.col("jacobs") >= -0.25) & (pl.col("jacobs") <= 0.25)).then("Indifference")
          .when((pl.col("jacobs") >= -0.5) & (pl.col("jacobs") < -0.25)).then("Moderate Avoidance")
          .when(pl.col("jacobs") < -0.5).then("Strong Avoidance")
          .otherwise("Unclassified").alias("preference_class")
    )

    return burn_comparison




def chunk_fires_three_bins(per_class_df: pl.DataFrame, peatland_classes: list[int]):
    cmp = classify_fire_preference(per_class_df, peatland_classes)
    cmp = cmp.with_columns(
        pl.when(pl.col("jacobs") > 0.25).then("preference")
          .when(pl.col("jacobs") < -0.25).then("avoidance")
          .when((pl.col("jacobs") >= -0.25) & (pl.col("jacobs") <= 0.25)).then("indifference")
          .otherwise("unclassified").alias("three_bin")
    )
    out = {}
    for label in ("preference", "indifference", "avoidance"):
        sel = cmp.filter(pl.col("three_bin") == label).select("fire_id").unique()
        out[label] = {
            "fire_ids": sel.to_series().to_list(),
            "n": sel.height
        }
    return out

# Example usage:
# per_class = pl.read_parquet("outputs/fire_selectivity_per_class.parquet")  # or result of summarize_lands_abundance
# peatland_classes = [130, 131]  # replace with your class IDs
# bins = chunk_fires_three_bins(per_class, peatland_classes)
# bins["preference"]["n"], bins["preference"]["fire_ids"][:10]