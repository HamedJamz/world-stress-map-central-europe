"""
Stress maps and charts.

Maps use ETRS89-LAEA (EPSG:3035). Because grid north is not exactly geographic north
everywhere in a projected map, the two end points of each SHmax line are computed on
the WGS84 ellipsoid (pyproj.Geod.fwd) and then projected, instead of rotating the line
on the map.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from pyproj import Geod, Transformer
from shapely.geometry import LineString, box

from .data_processing import (CRS_GEOGRAPHIC, CRS_MAP, QUALITY_ORDER, REGIME_NAMES,
                              TYPE_DESCRIPTIONS, load_boundaries)
from .orientation_stats import axial_mean, azimuth_to_compass

GEOD = Geod(ellps="WGS84")
TO_MAP = Transformer.from_crs(CRS_GEOGRAPHIC, CRS_MAP, always_xy=True)

# ---- visual style ---------------------------------------------------------------------
REGIME_COLORS = {          # common WSM colour convention: NF red, SS green, TF blue
    "NF": "#d62728", "NS": "#f28e2b", "SS": "#2ca02c",
    "TS": "#5fa2dd", "TF": "#1f4e9c", "U": "#2b2b2b",
}
REGIME_ORDER = ["NF", "NS", "SS", "TS", "TF", "U"]
# symbol half-length (km, for a ~1300 km wide map) and line width per quality
QUALITY_HALF_LENGTH_KM = {"A": 21, "B": 17, "C": 13, "D": 9, "E": 6,
                          "Xmi": 6, "Xne": 6, "Xru": 6}
QUALITY_LINEWIDTH = {"A": 1.7, "B": 1.4, "C": 1.1, "D": 0.8, "E": 0.6,
                     "Xmi": 0.6, "Xne": 0.6, "Xru": 0.6}
QUALITY_COLORS = {"A": "#1b6e6a", "B": "#2a9d8f", "C": "#76c1b6",
                  "D": "#b9b9b9", "E": "#d6d6d6",
                  "Xmi": "#ece3cf", "Xne": "#ece3cf", "Xru": "#ece3cf"}
SEA, LAND, BORDER, TEXT = "#dde8f0", "#f7f5f0", "#7b7b7b", "#222222"
SOURCE_TEXT = ("Data: World Stress Map Database Release 2025 (Heidbach et al., 2025, "
               "GFZ Data Services, doi:10.5880/WSM.2025.001, CC BY 4.0). "
               "Boundaries: Natural Earth (public domain). Projection: ETRS89-LAEA (EPSG:3035).")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 11,
    "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------
def _save(fig, out_path, formats=("png",)):
    """Save a figure to one or several formats (same stem)."""
    if out_path is None:
        return
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    for ext in formats:
        fig.savefig(out_path.with_suffix("." + ext), dpi=300, facecolor="white")


def _projected_extent(extent):
    """Project a lon/lat extent (lon0, lat0, lon1, lat1) to EPSG:3035 map bounds."""
    lon0, lat0, lon1, lat1 = extent
    lons = np.r_[np.linspace(lon0, lon1, 50), np.full(50, lon1),
                 np.linspace(lon1, lon0, 50), np.full(50, lon0)]
    lats = np.r_[np.full(50, lat0), np.linspace(lat0, lat1, 50),
                 np.full(50, lat1), np.linspace(lat1, lat0, 50)]
    x, y = TO_MAP.transform(lons, lats)
    return (x.min(), y.min(), x.max(), y.max())


def shmax_segments(gdf: gpd.GeoDataFrame, scale: float = 1.0) -> np.ndarray:
    """Return projected end points of SHmax symbols, shape (n, 2, 2).

    Each symbol is centred on the observation and is 2 x half-length long; the
    half-length depends on the WSM quality (A longest) and is multiplied by ``scale``.
    End points are computed geodesically (see module docstring).
    """
    lon = gdf["LON"].to_numpy(float)
    lat = gdf["LAT"].to_numpy(float)
    az = gdf["AZI"].to_numpy(float)
    half = gdf["QUALITY"].astype(str).map(QUALITY_HALF_LENGTH_KM).fillna(6).to_numpy() \
        * 1000.0 * scale
    lon1, lat1, _ = GEOD.fwd(lon, lat, az, half)
    lon2, lat2, _ = GEOD.fwd(lon, lat, az + 180.0, half)
    x1, y1 = TO_MAP.transform(lon1, lat1)
    x2, y2 = TO_MAP.transform(lon2, lat2)
    return np.stack([np.c_[x1, y1], np.c_[x2, y2]], axis=1)


def _draw_basemap(ax, extent_xy, layers, highlight="Germany", draw_states=True):
    """Draw sea, countries, German federal states, rivers and lakes in EPSG:3035."""
    x0, y0, x1, y1 = extent_xy
    clip = gpd.GeoDataFrame(geometry=[box(x0, y0, x1, y1)], crs=CRS_MAP)
    ax.set_facecolor(SEA)
    c = layers["countries"].to_crs(CRS_MAP).clip(clip)
    c.plot(ax=ax, color=LAND, edgecolor="none", zorder=1)
    if highlight is not None:
        c[c["NAME"] == highlight].plot(ax=ax, color="#efe9dc", edgecolor="none", zorder=1.1)
    lakes = layers["lakes"].to_crs(CRS_MAP).clip(clip)
    if len(lakes):
        lakes.plot(ax=ax, color=SEA, edgecolor="#a9c3d6", linewidth=0.2, zorder=1.5)
    rivers = layers["rivers"]
    rivers = rivers[rivers["scalerank"] <= 7].to_crs(CRS_MAP).clip(clip)
    if len(rivers):
        rivers.plot(ax=ax, color="#a9c3d6", linewidth=0.45, zorder=1.6)
    if draw_states:
        st = layers["admin1_lines"]
        st = st[st["ADM0_A3"] == "DEU"].to_crs(CRS_MAP).clip(clip)
        if len(st):
            st.plot(ax=ax, color="#b3ab9a", linewidth=0.35, linestyle=(0, (3, 2)), zorder=1.8)
    c.boundary.plot(ax=ax, color=BORDER, linewidth=0.6, zorder=2)
    if highlight is not None:
        c[c["NAME"] == highlight].boundary.plot(ax=ax, color="#4a4a4a", linewidth=1.1, zorder=2.1)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(""); ax.set_ylabel("")
    for s in ax.spines.values():
        s.set_visible(True); s.set_color("#444444"); s.set_linewidth(0.8)


def _draw_graticule(ax, extent_xy, step_lon=4, step_lat=2):
    """Draw longitude/latitude lines and label them along the bottom and left frame."""
    x0, y0, x1, y1 = extent_xy
    for lon in np.arange(-12, 36, step_lon):
        lats = np.linspace(35, 65, 300)
        x, y = TO_MAP.transform(np.full_like(lats, lon), lats)
        ax.plot(x, y, color="#9aa7b0", lw=0.3, alpha=0.7, zorder=1.95)
        if y.min() < y0 < y.max():
            xb = np.interp(y0, y, x)
            if x0 < xb < x1:
                ax.text(xb, y0, f"{lon:g}°E", ha="center", va="top", fontsize=7,
                        color="#555555", transform=ax.transData, clip_on=False,
                        bbox=dict(fc="white", ec="none", pad=0.5))
    for lat in np.arange(36, 66, step_lat):
        lons = np.linspace(-12, 36, 300)
        x, y = TO_MAP.transform(lons, np.full_like(lons, lat))
        ax.plot(x, y, color="#9aa7b0", lw=0.3, alpha=0.7, zorder=1.95)
        if x.min() < x0 < x.max():
            yl = np.interp(x0, x, y)
            if y0 < yl < y1:
                ax.text(x0, yl, f"{lat:g}°N ", ha="right", va="center", fontsize=7,
                        color="#555555", clip_on=False)


def _draw_scalebar(ax, extent_xy, length_km=None):
    """Draw a simple two-segment scale bar (valid for the LAEA projection centre area)."""
    x0, y0, x1, y1 = extent_xy
    width_km = (x1 - x0) / 1000
    if length_km is None:
        length_km = 200 if width_km > 900 else 100 if width_km > 400 else 50
    L = length_km * 1000
    bx, by = x0 + 0.04 * (x1 - x0), y0 + 0.04 * (y1 - y0)
    h = 0.008 * (y1 - y0)
    ax.add_patch(Rectangle((bx - 0.02 * (x1 - x0), by - 1.6 * h), L + 0.06 * (x1 - x0),
                           6.0 * h, fc="white", ec="#bbbbbb", lw=0.4, alpha=0.9, zorder=5))
    ax.add_patch(Rectangle((bx, by), L / 2, h, fc="#222222", ec="#222222", lw=0.6, zorder=6))
    ax.add_patch(Rectangle((bx + L / 2, by), L / 2, h, fc="white", ec="#222222", lw=0.6, zorder=6))
    for frac, lab in [(0, "0"), (0.5, f"{length_km // 2:g}"), (1, f"{length_km:g} km")]:
        ax.text(bx + frac * L, by + 1.6 * h, lab, ha="center", va="bottom", fontsize=7, zorder=7)


def _draw_north_arrow(ax, extent_xy):
    """Draw a north arrow that points to geographic north at its own map position."""
    x0, y0, x1, y1 = extent_xy
    px, py = x1 - 0.06 * (x1 - x0), y1 - 0.13 * (y1 - y0)
    lon, lat = Transformer.from_crs(CRS_MAP, CRS_GEOGRAPHIC, always_xy=True).transform(px, py)
    qx, qy = TO_MAP.transform(lon, lat + 0.5)
    d = np.array([qx - px, qy - py]); d /= np.linalg.norm(d)
    L = 0.06 * (y1 - y0)
    ax.annotate("", xy=(px + d[0] * L, py + d[1] * L), xytext=(px, py),
                arrowprops=dict(arrowstyle="-|>", color="#222222", lw=1.2), zorder=7)
    ax.text(px + d[0] * L * 1.25, py + d[1] * L * 1.25, "N", ha="center", va="center",
            fontsize=9, fontweight="bold", zorder=7)


def _draw_symbols(ax, gdf, scale):
    """Draw SHmax symbols (lines) coloured by regime, length/width by quality."""
    if len(gdf) == 0:
        return
    d = gdf.copy()
    d["_q"] = pd.Categorical(d["QUALITY"].astype(str), categories=QUALITY_ORDER, ordered=True)
    # draw low quality first so that A-quality symbols end up on top
    d = d.sort_values("_q", ascending=False)
    seg = shmax_segments(d, scale)
    colors = d["REGIME"].map(REGIME_COLORS).fillna("#2b2b2b").to_numpy()
    lws = d["QUALITY"].astype(str).map(QUALITY_LINEWIDTH).fillna(0.6).to_numpy()
    ax.add_collection(LineCollection(seg, colors=colors, linewidths=lws,
                                     capstyle="round", zorder=4))
    cx, cy = TO_MAP.transform(d["LON"].to_numpy(), d["LAT"].to_numpy())
    ax.scatter(cx, cy, s=2.5, c=colors, linewidths=0, zorder=4.1)


def _legend_panel(lax, gdf, scale_ref_km, show_stats=True):
    """Side panel: regime colours, quality symbol lengths, record counts and mean SHmax."""
    lax.set_xlim(0, 1); lax.set_ylim(0, 1); lax.axis("off")
    y = 0.98
    lax.text(0, y, "Stress regime", fontweight="bold", va="top"); y -= 0.045
    present = [r for r in REGIME_ORDER if r in set(gdf["REGIME"].astype(str))] or REGIME_ORDER
    for r in REGIME_ORDER:
        alpha = 1.0 if r in present else 0.35
        lax.plot([0.02, 0.17], [y, y], color=REGIME_COLORS[r], lw=2.2, alpha=alpha,
                 solid_capstyle="round")
        n = int((gdf["REGIME"].astype(str) == r).sum())
        lax.text(0.22, y, f"{REGIME_NAMES[r]} ({r})  n={n}", va="center", fontsize=8,
                 alpha=alpha)
        y -= 0.04
    y -= 0.02
    lax.text(0, y, "WSM quality (symbol length)", fontweight="bold", va="top"); y -= 0.05
    acc = {"A": "±15°", "B": "±20°", "C": "±25°", "D": "±40°", "E": ">40°/no info"}
    qs = [q for q in QUALITY_ORDER[:5] if q in set(gdf["QUALITY"].astype(str))] or ["A", "B", "C"]
    for q in qs:
        L = QUALITY_HALF_LENGTH_KM[q] / QUALITY_HALF_LENGTH_KM["A"] * 0.16
        lax.plot([0.095 - L / 2, 0.095 + L / 2], [y, y], color="#333333",
                 lw=QUALITY_LINEWIDTH[q] * 1.2, solid_capstyle="round")
        n = int((gdf["QUALITY"].astype(str) == q).sum())
        lax.text(0.22, y, f"{q}  ({acc[q]})  n={n}", va="center", fontsize=8)
        y -= 0.04
    lax.text(0, y - 0.005, "Line = SHmax axis (orientation of max.\nhorizontal compression)",
             fontsize=7.2, va="top", color="#444444", style="italic")
    y -= 0.095
    if show_stats and len(gdf):
        st = axial_mean(gdf["AZI"])
        lax.text(0, y, "Selection summary", fontweight="bold", va="top"); y -= 0.045
        lines = [f"Records shown: {len(gdf)}",
                 f"Axial mean SHmax: {st['mean_deg']:.0f}° ({azimuth_to_compass(st['mean_deg'])})",
                 f"Mean resultant length R: {st['R']:.2f}"]
        for t in lines:
            lax.text(0, y, t, fontsize=8, va="top"); y -= 0.037
    return y


def _overview_inset(fig, rect, extent, layers):
    """Small overview map of Europe with the study-area box."""
    iax = fig.add_axes(rect)
    ext = (-11, 35.5, 32, 63)
    exy = _projected_extent(ext)
    c = layers["countries"].to_crs(CRS_MAP)
    iax.set_facecolor(SEA)
    c.plot(ax=iax, color=LAND, edgecolor="#9a9a9a", linewidth=0.25)
    c[c["NAME"] == "Germany"].plot(ax=iax, color="#d9cfb8", edgecolor="#555555", linewidth=0.4)
    lon0, lat0, lon1, lat1 = extent
    bl = LineString([(lon0, lat0), (lon1, lat0), (lon1, lat1), (lon0, lat1), (lon0, lat0)])
    bl = gpd.GeoSeries([bl.segmentize(0.25)], crs=CRS_GEOGRAPHIC).to_crs(CRS_MAP)
    bl.plot(ax=iax, color="#c0392b", linewidth=1.0)
    iax.set_xlim(exy[0] + 2.5e5, exy[2] - 4e5); iax.set_ylim(exy[1] + 2e5, exy[3] - 3e5)
    iax.set_xticks([]); iax.set_yticks([]); iax.set_aspect("equal")
    iax.set_xlabel(""); iax.set_ylabel("")
    iax.set_title("Study-area location", fontsize=8, loc="left", pad=3)
    for s in iax.spines.values():
        s.set_color("#666666"); s.set_linewidth(0.6)
    return iax


# --------------------------------------------------------------------------------------
# Main map
# --------------------------------------------------------------------------------------
def plot_stress_map(gdf: gpd.GeoDataFrame, extent=(2.5, 45.5, 19.5, 55.5), title=None,
                    subtitle=None, out_path=None, formats=("png",), highlight="Germany",
                    figsize=(11.5, 9.2), symbol_scale=None, show_inset=True):
    """Draw a stress map of SHmax orientations.

    Symbols: line through each observation, oriented along SHmax (geodesic, correct
    w.r.t. geographic north), colour = stress regime, length/width = WSM quality.
    Only records with a valid azimuth are drawn.
    """
    layers = load_boundaries()
    exy = _projected_extent(extent)
    width_km = (exy[2] - exy[0]) / 1000
    scale = symbol_scale if symbol_scale is not None else max(min(width_km / 1300, 1.3), 0.35)
    d = gdf[gdf["AZI"].notna()]

    # --- layout in inches: map keeps its true aspect ratio, legend panel on the right
    fig_w = figsize[0]
    map_w = 0.64 * fig_w
    map_h = map_w * (exy[3] - exy[1]) / (exy[2] - exy[0])
    top, bottom = 0.95, 0.75                     # inches for title / source text
    body_h = max(map_h, 6.8)                     # legend panel needs ~6.8 in of height
    fig_h = body_h + top + bottom
    fig = plt.figure(figsize=(fig_w, fig_h))
    fx = lambda inch: inch / fig_w
    fy = lambda inch: inch / fig_h
    ax = fig.add_axes([fx(0.55), fy(bottom + body_h - map_h), fx(map_w), fy(map_h)])
    panel_x = fx(0.55 + map_w + 0.45)
    lax = fig.add_axes([panel_x, fy(bottom), 1 - panel_x - fx(0.15), fy(body_h * 0.64)])
    _draw_basemap(ax, exy, layers, highlight=highlight)
    _draw_graticule(ax, exy, step_lon=4 if width_km > 600 else 2,
                    step_lat=2 if width_km > 600 else 1)
    _draw_symbols(ax, d, scale)
    _draw_scalebar(ax, exy)
    _draw_north_arrow(ax, exy)
    _legend_panel(lax, d, QUALITY_HALF_LENGTH_KM["A"] * scale)
    if show_inset:
        inset_h = body_h * 0.30
        _overview_inset(fig, [panel_x, fy(bottom + body_h - inset_h),
                              1 - panel_x - fx(0.15), fy(inset_h)], extent, layers)

    fig.text(fx(0.55), fy(fig_h - 0.18), title or "World Stress Map 2025 – SHmax orientations",
             fontsize=15, fontweight="bold", color=TEXT, va="top")
    if subtitle:
        fig.text(fx(0.55), fy(fig_h - 0.55), subtitle, fontsize=9.5, color="#444444", va="top")
    fig.text(fx(0.55), fy(0.42), SOURCE_TEXT + "\nIndependent student project – "
             "not an official WSM/GFZ product.", fontsize=6.8, color="#555555", va="top")
    _save(fig, out_path, formats)
    return fig


# --------------------------------------------------------------------------------------
# Charts
# --------------------------------------------------------------------------------------
def plot_quality_distribution(df, title=None, out_path=None, ax=None, compact=False):
    """Bar chart: number of records per WSM quality class (A-C highlighted as reliable)."""
    counts = df["QUALITY"].astype(str).value_counts().reindex(QUALITY_ORDER, fill_value=0)
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(7.2, 4.2))
    else:
        fig = ax.figure
    bars = ax.bar(counts.index, counts.values,
                  color=[QUALITY_COLORS[q] for q in counts.index],
                  edgecolor="#555555", linewidth=0.5)
    for b, v in zip(bars, counts.values):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v}", ha="center", va="bottom",
                fontsize=7.5 if compact else 8.5)
    n_ac = int(counts[["A", "B", "C"]].sum())
    ax.set_ylabel("Number of records")
    ax.set_xlabel("WSM quality class")
    ax.set_ylim(0, counts.max() * 1.15)
    ax.set_title(title or "WSM records by quality class", loc="left")
    if not compact:
        ax.text(0.99, 0.95, f"A–C (reliable): {n_ac} of {counts.sum()} "
                f"({100 * n_ac / max(counts.sum(), 1):.0f} %)",
                transform=ax.transAxes, ha="right", va="top", fontsize=9,
                bbox=dict(fc="white", ec="#cccccc"))
        ax.text(0.99, 0.80, "A ±15°  B ±20°  C ±25°  D ±40°\nE: >40° or no orientation\n"
                "X: missing info / not established /\nrarely used indicator",
                transform=ax.transAxes, ha="right", va="top", fontsize=7, color="#555555")
    ax.grid(axis="y", color="#e6e6e6", lw=0.6); ax.set_axisbelow(True)
    if own:
        fig.tight_layout(); _save(fig, out_path)
    return fig


def plot_indicator_distribution(df, title=None, out_path=None, ax=None, compact=False):
    """Horizontal stacked bars: records per stress-indicator TYPE, split A-C vs. D/E/X."""
    tab = (df.assign(rel=np.where(df["QUALITY"].astype(str).isin(["A", "B", "C"]),
                                  "A–C", "D, E, X"))
             .groupby(["TYPE", "rel"], observed=True).size().unstack(fill_value=0))
    for col in ["A–C", "D, E, X"]:
        if col not in tab:
            tab[col] = 0
    tab = tab[["A–C", "D, E, X"]]
    tab = tab.loc[tab.sum(axis=1).sort_values().index]
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(7.6, 0.33 * len(tab) + 1.6))
    else:
        fig = ax.figure
    labels = [t if compact else f"{t} – {TYPE_DESCRIPTIONS.get(t, '')}" for t in tab.index]
    ax.barh(labels, tab["A–C"], color="#2a9d8f", label="Quality A–C", edgecolor="white", lw=0.4)
    ax.barh(labels, tab["D, E, X"], left=tab["A–C"], color="#c9c9c9",
            label="Quality D, E, X", edgecolor="white", lw=0.4)
    for i, (a, b) in enumerate(zip(tab["A–C"], tab["D, E, X"])):
        ax.text(a + b, i, f" {a + b}", va="center", fontsize=7.5 if compact else 8)
    ax.set_xlabel("Number of records")
    ax.set_xlim(0, tab.sum(axis=1).max() * 1.15)
    ax.set_title(title or "WSM records by stress-indicator type", loc="left")
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    ax.grid(axis="x", color="#e6e6e6", lw=0.6); ax.set_axisbelow(True)
    ax.tick_params(axis="y", labelsize=7.5 if compact else 8)
    if own:
        fig.tight_layout(); _save(fig, out_path)
    return fig


def _rose(ax, az, bin_width=10, color="#2a9d8f", show_mean=True, title=None):
    """Axial rose diagram: each azimuth is counted at AZI and AZI+180 (symmetric petals)."""
    az = np.asarray(az, float); az = az[np.isfinite(az)]
    edges = np.arange(0, 181, bin_width)
    counts, _ = np.histogram(az % 180, bins=edges)
    theta = np.deg2rad(np.r_[edges[:-1], edges[:-1] + 180] + bin_width / 2)
    r = np.r_[counts, counts]
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.bar(theta, r, width=np.deg2rad(bin_width), color=color, edgecolor="white",
           linewidth=0.5, alpha=0.9)
    if show_mean and az.size:
        st = axial_mean(az)
        for a in (st["mean_deg"], st["mean_deg"] + 180):
            ax.plot([np.deg2rad(a)] * 2, [0, r.max() * 1.05], color="#c0392b", lw=1.8)
    ax.set_xticks(np.deg2rad(np.arange(0, 360, 45)))
    ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"], fontsize=8)
    ax.tick_params(axis="y", labelsize=6.5, colors="#666666")
    ax.set_rlabel_position(67.5)
    if title:
        ax.set_title(title, fontsize=9.5, pad=12)


def plot_orientation_distribution(df, title=None, out_path=None, bin_width=10):
    """Two panels: axial rose diagram and 0-180 deg histogram with correct axial mean."""
    az = df["AZI"].dropna().to_numpy(float)
    st = axial_mean(az)
    naive = float(np.mean(az)) if az.size else np.nan
    fig = plt.figure(figsize=(11, 4.6))
    ax1 = fig.add_subplot(1, 2, 1, projection="polar")
    _rose(ax1, az, bin_width)
    ax1.set_title(f"Rose diagram (axial, {bin_width}° bins)", fontsize=10, pad=14)
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.hist(az, bins=np.arange(0, 181, bin_width), color="#2a9d8f", edgecolor="white")
    ax2.axvline(st["mean_deg"], color="#c0392b", lw=2,
                label=f"Axial (circular) mean: {st['mean_deg']:.0f}°")
    ax2.axvline(naive, color="#555555", lw=1.4, ls="--",
                label=f"Naive arithmetic mean: {naive:.0f}° (not valid for axes)")
    ax2.set_xticks(range(0, 181, 30))
    ax2.set_xlabel("SHmax azimuth (° clockwise from north; 0° ≡ 180°)")
    ax2.set_ylabel("Number of records")
    ax2.set_title("Histogram of SHmax azimuths", loc="left", fontsize=10)
    ax2.legend(frameon=False, fontsize=8, loc="upper left")
    ax2.grid(axis="y", color="#e6e6e6", lw=0.6); ax2.set_axisbelow(True)
    fig.suptitle(title or "SHmax orientation distribution", x=0.02, ha="left",
                 fontsize=13, fontweight="bold")
    fig.text(0.02, -0.02, f"n = {st['n']}   |   axial mean = {st['mean_deg']:.1f}° "
             f"({azimuth_to_compass(st['mean_deg'])})   |   mean resultant length R = "
             f"{st['R']:.2f}   |   circular SD = {st['circ_std_deg']:.0f}°",
             fontsize=8.5, color="#444444")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    _save(fig, out_path)
    return fig


def plot_rose_by_group(df, group_col="TYPE", groups=None, out_path=None, title=None,
                       bin_width=15):
    """Small-multiple rose diagrams, one per indicator type (or other grouping)."""
    d = df[df["AZI"].notna()]
    if groups is None:
        groups = d[group_col].value_counts().index[:4].tolist()
    fig = plt.figure(figsize=(3.2 * len(groups), 3.8))
    for i, g in enumerate(groups):
        sub = d[d[group_col] == g]["AZI"]
        ax = fig.add_subplot(1, len(groups), i + 1, projection="polar")
        st = axial_mean(sub)
        _rose(ax, sub, bin_width, title=f"{g}\n n={st['n']}, mean {st['mean_deg']:.0f}°, R={st['R']:.2f}")
    fig.suptitle(title or "SHmax orientations by indicator type", x=0.02, ha="left",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    _save(fig, out_path)
    return fig


def plot_gridded_means(grid, extent=(2.5, 45.5, 19.5, 55.5), cell_deg=1.0, min_n=5,
                       out_path=None, title=None):
    """Map of axial mean SHmax per grid cell; bar colour = mean resultant length R."""
    layers = load_boundaries()
    exy = _projected_extent(extent)
    fig_w = 9.5
    map_w = 7.4
    map_h = map_w * (exy[3] - exy[1]) / (exy[2] - exy[0])
    fig = plt.figure(figsize=(fig_w, map_h + 1.6))
    ax = fig.add_axes([0.5 / fig_w, 0.7 / (map_h + 1.6), map_w / fig_w, map_h / (map_h + 1.6)])
    _draw_basemap(ax, exy, layers)
    _draw_graticule(ax, exy)
    cmap = plt.get_cmap("viridis_r")
    norm = plt.Normalize(0.2, 1.0)
    if len(grid):
        lon1, lat1, _ = GEOD.fwd(grid["cell_lon"].values, grid["cell_lat"].values,
                                 grid["mean_deg"].values, np.full(len(grid), 33000.0))
        lon2, lat2, _ = GEOD.fwd(grid["cell_lon"].values, grid["cell_lat"].values,
                                 grid["mean_deg"].values + 180, np.full(len(grid), 33000.0))
        x1, y1 = TO_MAP.transform(lon1, lat1); x2, y2 = TO_MAP.transform(lon2, lat2)
        seg = np.stack([np.c_[x1, y1], np.c_[x2, y2]], axis=1)
        lw = 1.2 + 2.2 * np.clip((grid["n"].values - min_n) / 40, 0, 1)
        ax.add_collection(LineCollection(seg, colors=cmap(norm(grid["R"].values)),
                                         linewidths=lw, capstyle="round", zorder=4))
    _draw_scalebar(ax, exy)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cax = fig.add_axes([(0.5 + map_w + 0.2) / fig_w, 0.7 / (map_h + 1.6) + 0.1,
                        0.18 / fig_w, 0.5])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("Mean resultant length R\n(1 = very consistent orientations)", fontsize=8)
    fig.text(0.5 / fig_w, 1 - 0.2 / (map_h + 1.6),
             title or f"Mean SHmax orientation per {cell_deg:g}° cell (A–C, n ≥ {min_n})",
             fontsize=13, fontweight="bold", va="top")
    fig.text(0.5 / fig_w, 1 - 0.55 / (map_h + 1.6),
             "Simple grid-cell axial means for illustration – not the official WSM smoothing method. "
             "Line width grows with the number of records.", fontsize=8, color="#444444", va="top")
    fig.text(0.5 / fig_w, 0.35 / (map_h + 1.6), SOURCE_TEXT, fontsize=6.5, color="#555555")
    _save(fig, out_path)
    return fig


def plot_quality_comparison(gdf_all, gdf_ac, extent=(2.5, 45.5, 19.5, 55.5), out_path=None):
    """Two maps side by side: all records with an azimuth vs. A-C quality only."""
    layers = load_boundaries()
    exy = _projected_extent(extent)
    scale = max(min((exy[2] - exy[0]) / 1000 / 1300, 1.3), 0.35)
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.6))
    for ax, d, lab in [(axes[0], gdf_all, "All qualities (A–E, X) with an SHmax azimuth"),
                       (axes[1], gdf_ac, "Quality A–C only (recommended)")]:
        d = d[d["AZI"].notna()]
        _draw_basemap(ax, exy, layers)
        _draw_symbols(ax, d, scale)
        st = axial_mean(d["AZI"])
        ax.set_title(f"{lab}\n n = {len(d)}, axial mean {st['mean_deg']:.0f}°, R = {st['R']:.2f}",
                     fontsize=10, loc="left")
    handles = [Line2D([0], [0], color=REGIME_COLORS[r], lw=2.2, label=f"{REGIME_NAMES[r]}")
               for r in REGIME_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Effect of quality filtering on the Central European stress pattern",
                 x=0.02, ha="left", fontsize=13, fontweight="bold")
    fig.text(0.02, -0.03, SOURCE_TEXT, fontsize=6.8, color="#555555")
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    _save(fig, out_path)
    return fig
