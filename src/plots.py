#####################################################################
# FIRE SELECTIVITY PLOTS
# Jack A. Goldman
# Converted from R to Python using plotnine
#####################################################################

# required libraries
import pandas as pd
import numpy as np
from plotnine import (
    ggplot, aes,
    geom_point, geom_errorbar, geom_bar, geom_text, geom_boxplot,
    geom_hline, geom_segment, geom_label, geom_violin,
    scale_x_continuous, scale_y_continuous, scale_fill_manual,
    scale_fill_brewer, scale_color_brewer, scale_alpha_continuous,
    facet_grid, facet_wrap,
    labs, xlab, ylab,
    theme, theme_classic, theme_bw,
    element_text, element_blank,
    guides, guide_legend,
    arrow, position_dodge, position_dodge2,
)
from plotnine.scales import scale_color_brewer
import warnings
warnings.filterwarnings("ignore")

# set file path for working directory
pathname = "/Users/jgoldman/Library/CloudStorage/GoogleDrive-jandrewgoldman@gmail.com/My Drive/3_POSTDOC/projects/selectivity/"

# read in selectivity data
selectivity = pd.read_parquet("data/full_dataset/fire_selectivity_full_dataset_20260318.parquet")

my_colors = ["yellow", "orange", "red", "darkred"]


####### MEAN SELECTIVITY BY CLASS WITH ARROWS ########

selectivity_mean = (
    selectivity
    .groupby("variable")["Chesson"]
    .agg(
        mean=lambda x: x.mean(),
        median=lambda x: x.median(),
        sd=lambda x: x.std(),
        p10=lambda x: x.quantile(0.10),
        p90=lambda x: x.quantile(0.90),
    )
    .reset_index()
)


##### Plot with arrows indicating increasing/decreasing selectivity from "upland" category

sel_plot = (
    selectivity_mean
    .loc[~selectivity_mean["variable"].isin(["total_peat", "mineral"])]
    .sort_values("mean")
    .reset_index(drop=True)
)
sel_plot["x"] = range(1, len(sel_plot) + 1)

# find upland index (1-based, matching R row_number)
upland_mask = sel_plot["variable"].str.lower() == "upland"
if not upland_mask.any():
    raise ValueError("No 'upland' category found in variable column.")
upland_idx = sel_plot.index[upland_mask][0] + 1  # 1-based

# compute vertical position just above the upper CI for the upland category
y_range = sel_plot["p90"].max() - sel_plot["p10"].min()
y_pos_above = sel_plot.loc[upland_mask, "p90"].values[0] + y_range * 0.05

# arrow geometry
left_start  = upland_idx - 0.25
left_end    = max(1, upland_idx - 2)
right_start = upland_idx + 0.25
right_end   = min(len(sel_plot), upland_idx + 2)

left_mid  = (left_start + left_end) / 2
right_mid = (right_start + right_end) / 2

label_offset = y_range * 0.02
label_y = y_pos_above + label_offset

# base plot using numeric x but label axis with category names
p = (
    ggplot(sel_plot, aes(x="x", y="mean"))
    + geom_point()
    + geom_errorbar(aes(ymin="p10", ymax="p90"), width=0.2)
    + scale_x_continuous(breaks=sel_plot["x"].tolist(), labels=sel_plot["variable"].tolist())
    + xlab("")
    + ylab("Chesson's index")
    + theme_classic()
    + theme(axis_text_x=element_text(angle=90, vjust=0.5))
)

# add arrows and labels
p_arrows = (
    p
    + geom_segment(
        aes(x=left_start, xend=left_end, y=y_pos_above, yend=y_pos_above),
        inherit_aes=False,
        arrow=arrow(length=0.3),
        color="black",
    )
    + geom_segment(
        aes(x=right_start, xend=right_end, y=y_pos_above, yend=y_pos_above),
        inherit_aes=False,
        arrow=arrow(length=0.3),
        color="black",
    )
    + geom_label(
        aes(x=right_mid, y=label_y, label="Increasing"),
        inherit_aes=False,
        fontweight="bold", ha="center", va="bottom",
        fill="white", size=9,
    )
    + geom_label(
        aes(x=left_mid, y=label_y, label="Decreasing"),
        inherit_aes=False,
        fontweight="bold", ha="center", va="bottom",
        fill="white", size=9,
    )
)

p_arrows.save(pathname + "figures/mean_selectivity_by_class.png", width=10, height=6, dpi=150)


## same plot but with "total_" categories removed

sel_plot2 = (
    selectivity_mean
    .loc[~selectivity_mean["variable"].str.startswith("total_")]
    .sort_values("mean")
    .reset_index(drop=True)
)
sel_plot2["x"] = range(1, len(sel_plot2) + 1)

upland_mask2 = sel_plot2["variable"].str.lower() == "upland"
if not upland_mask2.any():
    raise ValueError("No 'upland' category found in variable column.")
upland_idx2 = sel_plot2.index[upland_mask2][0] + 1

y_range2 = sel_plot2["p90"].max() - sel_plot2["p10"].min()
y_pos_above2 = sel_plot2.loc[upland_mask2, "p90"].values[0] + y_range2 * 0.05

left_start2  = upland_idx2 - 0.25
left_end2    = max(1, upland_idx2 - 2)
right_start2 = upland_idx2 + 0.25
right_end2   = min(len(sel_plot2), upland_idx2 + 2)

left_mid2  = (left_start2 + left_end2) / 2
right_mid2 = (right_start2 + right_end2) / 2

label_offset2 = y_range2 * 0.02
label_y2 = y_pos_above2 + label_offset2

p2 = (
    ggplot(sel_plot2, aes(x="x", y="mean"))
    + geom_point()
    + geom_errorbar(aes(ymin="p10", ymax="p90"), width=0.2)
    + scale_x_continuous(breaks=sel_plot2["x"].tolist(), labels=sel_plot2["variable"].tolist())
    + xlab("")
    + ylab("Chesson's index")
    + theme_classic()
    + theme(axis_text_x=element_text(angle=90, vjust=0.5))
    + geom_segment(
        aes(x=left_start2, xend=left_end2, y=y_pos_above2, yend=y_pos_above2),
        inherit_aes=False,
        arrow=arrow(length=0.3),
        color="black",
    )
    + geom_segment(
        aes(x=right_start2, xend=right_end2, y=y_pos_above2, yend=y_pos_above2),
        inherit_aes=False,
        arrow=arrow(length=0.3),
        color="black",
    )
    + geom_label(
        aes(x=right_mid2, y=label_y2, label="Increasing"),
        inherit_aes=False,
        fontweight="bold", ha="center", va="bottom",
        fill="white", size=9,
    )
    + geom_label(
        aes(x=left_mid2, y=label_y2, label="Decreasing"),
        inherit_aes=False,
        fontweight="bold", ha="center", va="bottom",
        fill="white", size=9,
    )
)

p2.save(pathname + "figures/mean_selectivity_by_class_no_total.png", width=10, height=6, dpi=150)


########### BARPLOTS OF MEAN SELECTIVITY BY FWI AND BUI QUARTILES ###########

# FWI
selectivity_mean = (
    selectivity
    .dropna(subset=["class.bin"])
    .groupby(["Quartile_Bins.FWI", "class.bin", "canopy.bin"])["Chesson"]
    .agg(
        mean=lambda x: x.mean(),
        median=lambda x: x.median(),
        sd=lambda x: x.std(),
        p05=lambda x: x.quantile(0.05),
        p95=lambda x: x.quantile(0.95),
        count="count",
    )
    .reset_index()
)

# Ensure "total" is the last canopy.bin level
canopy_order = [c for c in selectivity_mean["canopy.bin"].unique() if c != "total"] + ["total"]
selectivity_mean["canopy.bin"] = pd.Categorical(selectivity_mean["canopy.bin"], categories=canopy_order, ordered=True)

fwi_mean_barplot = (
    ggplot(
        selectivity_mean.dropna(subset=["Quartile_Bins.FWI"]),
        aes(x="Quartile_Bins.FWI", y="mean", fill="Quartile_Bins.FWI"),
    )
    + geom_bar(stat="identity", alpha=0.9)
    + geom_text(aes(label="count"), va="bottom", nudge_y=0.002, color="black", size=8)
    + labs(x="", y="Chesson's index", fill="FWI quartiles")
    + scale_fill_manual(values=my_colors)
    + theme_bw()
    + facet_grid("class.bin ~ canopy.bin", scales="fixed")
    + theme(axis_text_x=element_text(angle=90))
)

fwi_mean_barplot.save(pathname + "figures/mean_selectivity_by_fwi_.png", width=10, height=8, dpi=150)


# BUI
selectivity_mean_bui = (
    selectivity
    .dropna(subset=["class.bin"])
    .groupby(["Quartile_Bins.BUI", "class.bin", "canopy.bin"])["Chesson"]
    .agg(
        mean=lambda x: x.mean(),
        median=lambda x: x.median(),
        sd=lambda x: x.std(),
        p05=lambda x: x.quantile(0.05),
        p95=lambda x: x.quantile(0.95),
        count="count",
    )
    .reset_index()
)

canopy_order_bui = [c for c in selectivity_mean_bui["canopy.bin"].unique() if c != "total"] + ["total"]
selectivity_mean_bui["canopy.bin"] = pd.Categorical(selectivity_mean_bui["canopy.bin"], categories=canopy_order_bui, ordered=True)

bui_median_barplot = (
    ggplot(
        selectivity_mean_bui.dropna(subset=["Quartile_Bins.BUI"]),
        aes(x="Quartile_Bins.BUI", y="median", fill="Quartile_Bins.BUI"),
    )
    + geom_bar(stat="identity", alpha=0.9)
    + geom_text(aes(label="count"), va="bottom", nudge_y=0.002, color="black", size=8)
    + labs(x="", y="Chesson's index", fill="BUI quartiles")
    + scale_fill_manual(values=my_colors)
    + theme_bw()
    + facet_grid("class.bin ~ canopy.bin")
    + theme(axis_text_x=element_text(angle=90))
)

bui_median_barplot.save(pathname + "figures/median_selectivity_by_bui_.png", width=10, height=8, dpi=150)


### mean chesson for fwi quartiles — free_x facets with counts

selectivity_mean_fwi2 = (
    selectivity
    .groupby(["Quartile_Bins.FWI", "class.bin", "canopy.bin"])["Chesson"]
    .agg(
        mean=lambda x: x.mean(),
        median=lambda x: x.median(),
        sd=lambda x: x.std(),
        p05=lambda x: x.quantile(0.05),
        p95=lambda x: x.quantile(0.95),
        count="count",
    )
    .reset_index()
)

p_fwi_free = (
    ggplot(
        selectivity_mean_fwi2
        .dropna(subset=["Quartile_Bins.FWI", "class.bin"]),
        aes(x="Quartile_Bins.FWI.astype('str')", y="mean", fill="Quartile_Bins.FWI"),
    )
    + geom_bar(stat="identity", alpha=0.9)
    + geom_text(aes(label="count"), va="bottom", nudge_y=0.002, color="black", size=8)
    + labs(x="", y="Chesson's index", fill="FWI quartiles")
    + theme_bw()
    + geom_hline(yintercept=0)
    + scale_fill_brewer(palette="Set1")
    + facet_wrap("class.bin + canopy.bin", scales="free_x", ncol=4)
    + guides(fill=guide_legend(reverse=True))
    + theme(
        axis_text_x=element_blank(),
        axis_ticks_major_x=element_blank(),
    )
)

p_fwi_free.save(pathname + "figures/mean_selectivity_by_fwi_free_x_with_counts.png", width=10, height=8, dpi=150)


####### East vs West boxplots

y_min = selectivity["Chesson"].min()
y_max = selectivity["Chesson"].max()
y_pad = 0.05 * (y_max - y_min)
y_limits = (y_min - y_pad, y_max + y_pad)
y_range = y_limits[1] - y_limits[0]

count_data = (
    selectivity
    .dropna(subset=["class.bin"])
    .groupby(["class.bin", "Quartile_Bins.BUI", "region"])
    .agg(count=("Chesson", "count"))
    .reset_index()
)
count_data["y_below"] = 0 - 0.04 * y_range

chesson_east_west = (
    ggplot(
        selectivity
        .dropna(subset=["Quartile_Bins.BUI", "class.bin"]),
        aes(fill="region", x="Quartile_Bins.BUI"),
    )
    + geom_boxplot(aes(y="Chesson"), alpha=0.5, outlier_alpha=0.6)
    + geom_text(
        data=count_data.dropna(subset=["Quartile_Bins.BUI"]),
        mapping=aes(label="count", color="region", y="y_below"),
        va="top",
        position=position_dodge(width=0.75),
        angle=90,
        size=8,
    )
    + labs(x="", y="Chesson's index", fill="Region")
    + theme_bw()
    + geom_hline(yintercept=0)
    + scale_fill_brewer(palette="Set1")
    + scale_color_brewer(palette="Set1", guide=False)
    + facet_wrap("class.bin", scales="fixed", ncol=5)
    + scale_y_continuous(limits=y_limits)
    + theme(axis_text_x=element_text(angle=90))
)

chesson_east_west.save(pathname + "figures/chesson_by_class_east_west_boxplot_.png", width=12, height=6, dpi=150)


########### Quartile bin BUI count as alpha in barplot #############

selectivity_mean_bui["alpha"] = round(selectivity_mean_bui["count"] / selectivity_mean_bui["count"].max(), 3) + 0.5

y_min_bar = selectivity_mean_bui["mean"].min()
y_max_bar = 0.10
y_pad_bar = 0.05 * (y_max_bar - y_min_bar)
y_limits_bar = (y_min_bar - y_pad_bar, y_max_bar + y_pad_bar)

bui_east_west_alpha = (
    ggplot(
        selectivity_mean_bui
        .dropna(subset=["Quartile_Bins.BUI", "class.bin"]),
        aes(x="Quartile_Bins.BUI", y="mean", fill="region"),
    )
    + geom_bar(aes(alpha="count"), stat="identity", position=position_dodge())
    + labs(x="", y="Chesson's index", fill="Region")
    + scale_fill_brewer(palette="Set1")
    + scale_color_brewer(palette="Set1")
    + theme_bw()
    + facet_grid("class.bin ~ canopy.bin", scales="fixed")
    + theme(axis_text_x=element_text(angle=90))
    + scale_y_continuous(limits=y_limits_bar)
    + scale_alpha_continuous(
        range=(0.3, 1),
        limits=(0, 1250),
        breaks=list(np.linspace(0, 1250, 6).astype(int)),
        labels=[str(v) for v in np.linspace(0, 1250, 6).astype(int)],
    )
)

bui_east_west_alpha.save(pathname + "figures/bui_east_west_barplot_alpha_counts.png", width=12, height=8, dpi=150)


########## Seasonal ############

selectivity_mean_season = (
    selectivity
    .groupby(["Quartile_Bins.BUI", "class.bin", "canopy.bin", "season.bin"])["Chesson"]
    .agg(
        mean=lambda x: x.mean(),
        median=lambda x: x.median(),
        sd=lambda x: x.std(),
        p05=lambda x: x.quantile(0.05),
        p95=lambda x: x.quantile(0.95),
        count="count",
    )
    .reset_index()
)

canopy_order_s = [c for c in selectivity_mean_season["canopy.bin"].unique() if c != "total"] + ["total"]
selectivity_mean_season["canopy.bin"] = pd.Categorical(selectivity_mean_season["canopy.bin"], categories=canopy_order_s, ordered=True)

y_min_s = selectivity_mean_season["mean"].min()
y_max_s = 0.15
y_pad_s = 0.05 * (y_max_s - y_min_s)
y_limits_s = (y_min_s - y_pad_s, y_max_s + y_pad_s)

bui_spring_fall_alpha = (
    ggplot(
        selectivity_mean_season.dropna(subset=["Quartile_Bins.BUI", "class.bin"]),
        aes(x="Quartile_Bins.BUI", y="mean", fill="season.bin"),
    )
    + geom_bar(aes(alpha="count"), stat="identity", position=position_dodge())
    + labs(x="", y="Chesson's index", fill="Season")
    + scale_fill_brewer(palette="Set1")
    + theme_bw()
    + facet_grid("class.bin ~ canopy.bin", scales="fixed")
    + theme(axis_text_x=element_text(angle=90))
    + scale_y_continuous(limits=y_limits_s)
    + scale_alpha_continuous(
        range=(0.3, 1),
        limits=(0, 1250),
        breaks=list(np.linspace(0, 1250, 6).astype(int)),
        labels=[str(v) for v in np.linspace(0, 1250, 6).astype(int)],
    )
)

bui_spring_fall_alpha.save(pathname + "figures/bui_early_late_barplot_alpha_counts.png", width=12, height=8, dpi=150)


## simple barplot with counts

y_max_simple = 0.13
y_pad_simple = 0.05 * (y_max_simple - y_min_s)
y_limits_simple = (y_min_s - y_pad_simple, y_max_simple + y_pad_simple)

bui_early_late_counts = (
    ggplot(
        selectivity_mean_season.dropna(subset=["Quartile_Bins.BUI", "class.bin"]),
        aes(x="Quartile_Bins.BUI", y="mean", fill="season.bin"),
    )
    + geom_bar(stat="identity", alpha=0.9, position=position_dodge(width=0.9), width=0.8)
    + geom_text(
        aes(label="count", y="mean"),
        position=position_dodge(width=0.9),
        va="bottom",
        nudge_y=0.001,
        color="black",
        size=7,
    )
    + labs(x="", y="Chesson's index", fill="Period")
    + scale_fill_brewer(palette="Set1")
    + theme_bw()
    + facet_grid("class.bin ~ canopy.bin", scales="free")
    + theme(axis_text_x=element_text(angle=90))
    + scale_y_continuous(limits=y_limits_simple)
)

bui_early_late_counts.save(pathname + "figures/bui_early_late_barplot_simple_counts.png", width=12, height=8, dpi=150)


# filter just western ecoregion
selectivity_mean_west = (
    selectivity[selectivity["region"] == "West"]
    .groupby(["Quartile_Bins.BUI", "class.bin", "canopy.bin", "season.bin"])["Chesson"]
    .agg(
        mean=lambda x: x.mean(),
        median=lambda x: x.median(),
        sd=lambda x: x.std(),
        p05=lambda x: x.quantile(0.05),
        p95=lambda x: x.quantile(0.95),
        count="count",
    )
    .reset_index()
)

canopy_order_w = [c for c in selectivity_mean_west["canopy.bin"].unique() if c != "total"] + ["total"]
selectivity_mean_west["canopy.bin"] = pd.Categorical(selectivity_mean_west["canopy.bin"], categories=canopy_order_w, ordered=True)

bui_early_late_west = (
    ggplot(
        selectivity_mean_west.dropna(subset=["Quartile_Bins.BUI", "class.bin"]),
        aes(x="Quartile_Bins.BUI", y="mean", fill="season.bin"),
    )
    + geom_bar(stat="identity", alpha=0.9, position=position_dodge(width=0.9), width=0.8)
    + geom_text(
        aes(label="count", y="mean"),
        position=position_dodge(width=0.9),
        va="bottom",
        nudge_y=0.001,
        color="black",
        size=7,
    )
    + labs(x="", y="Chesson's index", fill="Period")
    + scale_fill_brewer(palette="Set1")
    + theme_bw()
    + facet_grid("class.bin ~ canopy.bin", scales="free")
    + theme(axis_text_x=element_text(angle=90))
    + scale_y_continuous(limits=y_limits_simple)
)

bui_early_late_west.save(pathname + "figures/bui_early_late_west.png", width=12, height=8, dpi=150)


################ mean chessons for landcover bui col isi row ##########

selectivity_mean_bui_isi = (
    selectivity
    .groupby(["class.bin", "Quartile_Bins.ISI", "Quartile_Bins.BUI"])["Chesson"]
    .agg(
        mean=lambda x: x.mean(),
        median=lambda x: x.median(),
        sd=lambda x: x.std(),
        p05=lambda x: x.quantile(0.05),
        p95=lambda x: x.quantile(0.95),
        count="count",
    )
    .reset_index()
)

# Sort class.bin by mean for the x-axis ordering
class_order = (
    selectivity_mean_bui_isi.groupby("class.bin")["mean"]
    .mean()
    .sort_values()
    .index.tolist()
)
selectivity_mean_bui_isi["class.bin"] = pd.Categorical(
    selectivity_mean_bui_isi["class.bin"], categories=class_order, ordered=True
)

mean_chesson_bui_isi = (
    ggplot(
        selectivity_mean_bui_isi
        .dropna(subset=["Quartile_Bins.ISI", "Quartile_Bins.BUI", "class.bin"])
        .query("`class.bin` != 'total'"),
        aes(x="class.bin", y="mean", fill="class.bin"),
    )
    + geom_bar(stat="identity", alpha=0.9)
    + geom_text(aes(label="count"), va="bottom", nudge_y=0.002, color="black", size=8)
    + labs(x="BUI quartile", y="Mean Chesson's Index", fill="Land cover")
    + theme_bw()
    + geom_hline(yintercept=0)
    + scale_fill_brewer(palette="Set1")
    + facet_grid("Quartile_Bins.ISI ~ Quartile_Bins.BUI")
    + theme(
        axis_text_x=element_blank(),
        axis_ticks_major_x=element_blank(),
    )
)

mean_chesson_bui_isi.save(pathname + "figures/mean_chesson_bui_isi.png", width=12, height=8, dpi=150)


### spring-fall boxplot with counts

y_min_box = selectivity["Chesson"].min()
y_max_box = 1.0
y_pad_box = 0.05 * (y_max_box - y_min_box)
y_limits_box = (y_min_box - y_pad_box, y_max_box + y_pad_box)
y_range_box = y_limits_box[1] - y_limits_box[0]

count_data_season = (
    selectivity
    .dropna(subset=["class.bin"])
    .groupby(["class.bin", "Quartile_Bins.BUI", "season.bin"])
    .agg(count=("Chesson", "count"))
    .reset_index()
)
count_data_season["y_below"] = 0 - 0.04 * y_range_box

chesson_early_late_boxplot = (
    ggplot(
        selectivity.dropna(subset=["Quartile_Bins.BUI", "class.bin"]),
        aes(fill="season.bin", x="Quartile_Bins.BUI"),
    )
    + geom_boxplot(aes(y="Chesson"), alpha=0.5, outlier_alpha=0.6)
    + geom_text(
        data=count_data_season.dropna(subset=["Quartile_Bins.BUI"]),
        mapping=aes(label="count", color="season.bin", y="y_below"),
        va="top",
        angle=90,
        position=position_dodge(width=0.75),
        size=8,
    )
    + labs(x="", y="Chesson's index", fill="Season")
    + theme_bw()
    + geom_hline(yintercept=0)
    + scale_fill_brewer(palette="Set1")
    + scale_color_brewer(palette="Set1", guide=False)
    + facet_wrap("class.bin", scales="free", ncol=5)
    + scale_y_continuous(limits=y_limits_box)
    + theme(axis_text_x=element_text(angle=90))
)

chesson_early_late_boxplot.save(pathname + "figures/chesson_early_late_boxplot_with_counts.png", width=12, height=6, dpi=150)


########## Total fires count barplot ##########

dcount = (
    selectivity
    .dropna(subset=["class.bin"])
    .groupby(["class.bin", "canopy.bin"])
    .agg(total_fires=("Chesson", "count"))
    .reset_index()
)

total_fires_plot = (
    ggplot(dcount, aes(x="class.bin", y="total_fires", fill="canopy.bin"))
    + geom_bar(stat="identity", position=position_dodge2(width=0.9, preserve="single"), width=0.8)
    + labs(x="Land cover type", y="Total Fires", fill="Forest cover type")
    + theme_bw()
    + scale_fill_brewer(palette="Set1")
    + geom_text(
        aes(label="total_fires"),
        va="bottom",
        nudge_y=1,
        color="black",
        size=8,
        position=position_dodge2(width=0.9, preserve="single"),
    )
    + theme(axis_text_x=element_text(angle=90))
)

total_fires_plot.save(pathname + "figures/total_fires_by_class_canopy.png", width=10, height=6, dpi=150)


########## Violin plot: FWI quartiles by class ##########

selectivity_mean_fwi_violin = (
    selectivity
    .groupby(["Quartile_Bins.FWI", "class.bin", "canopy.bin"])["Chesson"]
    .agg(
        mean=lambda x: x.mean(),
        median=lambda x: x.median(),
        sd=lambda x: x.std(),
        p05=lambda x: x.quantile(0.05),
        p95=lambda x: x.quantile(0.95),
        count="count",
    )
    .reset_index()
)

violin_fwi_class = (
    ggplot(
        selectivity
        .dropna(subset=["Quartile_Bins.FWI"])
        .query("variable != 'water' and variable != 'mineral'")
        .replace([np.inf, -np.inf], np.nan)
        .dropna(subset=["Chesson"]),
        aes(y="Chesson", fill="Quartile_Bins.FWI", x="Quartile_Bins.FWI"),
    )
    + geom_violin(alpha=0.5)
    + labs(x="", y="Chesson's index", fill="FWI quartiles")
    + theme_bw()
    + geom_hline(yintercept=0)
    + scale_fill_manual(values=my_colors)
    + facet_grid(". ~ class.bin", scales="free")
    + theme(
        axis_text_x=element_text(angle=90, size=20),
        axis_text_y=element_text(size=20),
        strip_text=element_text(size=20),
        legend_text=element_text(size=20),
        legend_title=element_text(size=22),
    )
)

violin_fwi_class.save(pathname + "figures/violin_fwi_class.png", width=20, height=20, dpi=150)

print("All plots saved successfully.")