"""
Filters for the WSM data and a small CASMO-like workflow:

    data -> filters (region, quality, indicator, regime, depth) -> selection -> map

This only copies the basic idea of CASMO (Create A Stress Map Online, WSM/GFZ);
it is not CASMO and not official WSM software.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import geopandas as gpd
import pandas as pd
from shapely.geometry import box

from .data_processing import QUALITY_ORDER, INDICATOR_GROUPS, load_boundaries

# Named regions. Bounding boxes are (lon_min, lat_min, lon_max, lat_max) in degrees.
REGION_PRESETS: dict[str, dict] = {
    # Primary study area: Germany, selected with the Natural Earth 1:10m country polygon
    "Germany": {"kind": "polygon", "country_name": "Germany",
                "map_extent": (5.3, 46.9, 15.6, 55.3)},
    # Extended study area: Central Europe window (documented in README)
    "Central Europe": {"kind": "bbox", "bbox": (2.5, 45.5, 19.5, 55.5),
                       "map_extent": (2.5, 45.5, 19.5, 55.5)},
    "Switzerland": {"kind": "polygon", "country_name": "Switzerland",
                    "map_extent": (5.8, 45.7, 10.6, 47.9)},
    "Austria": {"kind": "polygon", "country_name": "Austria",
                "map_extent": (9.4, 46.3, 17.3, 49.1)},
    "Netherlands": {"kind": "polygon", "country_name": "Netherlands",
                    "map_extent": (3.2, 50.7, 7.3, 53.7)},
    "Upper Rhine Graben area": {"kind": "bbox", "bbox": (6.8, 47.4, 9.0, 50.2),
                                "map_extent": (6.8, 47.4, 9.0, 50.2)},
}


# --------------------------------------------------------------------------------------
# Individual filters
# --------------------------------------------------------------------------------------
def filter_by_region(gdf: gpd.GeoDataFrame, region=None) -> gpd.GeoDataFrame:
    """Select records inside a region.

    ``region`` may be
      * None                          -> no spatial filter
      * a preset name (see REGION_PRESETS), e.g. "Germany", "Central Europe"
      * a bounding box tuple (lon_min, lat_min, lon_max, lat_max)
      * a shapely Polygon / MultiPolygon in EPSG:4326
      * a list of WSM COUNTRY names, e.g. ["Germany", "Austria"] (uses the WSM field)
    """
    if region is None:
        return gdf.copy()
    if isinstance(region, str) and region in REGION_PRESETS:
        spec = REGION_PRESETS[region]
        if spec["kind"] == "bbox":
            return filter_by_region(gdf, spec["bbox"])
        countries = load_boundaries()["countries"]
        poly = countries.loc[countries["NAME"] == spec["country_name"], "geometry"].union_all()
        return gdf[gdf.within(poly)].copy()
    if isinstance(region, (tuple, list)) and len(region) == 4 and all(
            isinstance(v, (int, float)) for v in region):
        x0, y0, x1, y1 = region
        return gdf[gdf["LON"].between(x0, x1) & gdf["LAT"].between(y0, y1)].copy()
    if hasattr(region, "geom_type"):
        return gdf[gdf.within(region)].copy()
    if isinstance(region, (str, list, tuple)):
        names = [region] if isinstance(region, str) else list(region)
        return gdf[gdf["COUNTRY"].isin(names)].copy()
    raise ValueError(f"Unsupported region specification: {region!r}")


def _resolve_quality(quality) -> list[str]:
    """Translate a quality argument into an explicit list of allowed classes.

    "C" (a single A-E letter) means *minimum* quality C -> ["A", "B", "C"];
    a list such as ["A", "B"] is used as given; None means all classes.
    """
    if quality is None:
        return QUALITY_ORDER
    if isinstance(quality, str):
        if quality in "ABCDE" and len(quality) == 1:
            return QUALITY_ORDER[: "ABCDE".index(quality) + 1]
        return [quality]
    return list(quality)


def filter_by_quality(gdf: gpd.GeoDataFrame, quality=None) -> gpd.GeoDataFrame:
    """Keep only records whose WSM quality class is allowed (see ``_resolve_quality``)."""
    allowed = _resolve_quality(quality)
    return gdf[gdf["QUALITY"].astype(str).isin(allowed)].copy()


def filter_by_indicator(gdf: gpd.GeoDataFrame, indicator=None) -> gpd.GeoDataFrame:
    """Keep only records of the chosen stress-indicator type(s).

    Accepts exact TYPE codes ("BO", "FMS"), a prefix ending in "*" ("FM*", "HF*")
    or a group name from INDICATOR_GROUPS ("Focal mechanism", "Borehole breakout").
    """
    if indicator is None:
        return gdf.copy()
    items = [indicator] if isinstance(indicator, str) else list(indicator)
    types: set[str] = set()
    for it in items:
        if it.endswith("*"):
            types |= {t for t in INDICATOR_GROUPS if t.startswith(it[:-1])}
        elif it in INDICATOR_GROUPS.values():
            types |= {t for t, g in INDICATOR_GROUPS.items() if g == it}
        else:
            types.add(it)
    return gdf[gdf["TYPE"].isin(types)].copy()


def filter_by_regime(gdf: gpd.GeoDataFrame, regime=None) -> gpd.GeoDataFrame:
    """Keep only records with the given stress regime code(s): NF, NS, SS, TS, TF, U."""
    if regime is None:
        return gdf.copy()
    items = [regime] if isinstance(regime, str) else list(regime)
    return gdf[gdf["REGIME"].isin(items)].copy()


def filter_by_depth(gdf: gpd.GeoDataFrame, depth_range=None) -> gpd.GeoDataFrame:
    """Keep records with DEPTH (km) inside ``(min_km, max_km)``."""
    if depth_range is None:
        return gdf.copy()
    lo, hi = depth_range
    return gdf[gdf["DEPTH"].between(lo, hi)].copy()


# --------------------------------------------------------------------------------------
# Combined workflow
# --------------------------------------------------------------------------------------
def filter_stress_data(gdf: gpd.GeoDataFrame, region=None, quality=None, indicator=None,
                       regime=None, depth_range=None, require_azimuth: bool = True,
                       make_map: bool = False, out_path: str | Path | None = None,
                       title: str | None = None, verbose: bool = True):
    """CASMO-inspired filtering workflow: filter WSM records and (optionally) map them.

    Parameters
    ----------
    gdf            : cleaned WSM GeoDataFrame (output of clean_wsm_data + to_geodataframe)
    region         : preset name, bbox tuple, polygon or list of WSM COUNTRY names
    quality        : minimum class ("C") or list of allowed classes (["A", "B"])
    indicator      : TYPE code(s), prefix ("FM*") or group name ("Borehole breakout")
    regime         : regime code(s) such as "TF" or ["SS", "NF"]
    depth_range    : (min_km, max_km)
    require_azimuth: drop records without SHmax orientation (needed for a stress map)
    make_map       : if True, draw a stress map of the selection
    out_path       : file path for the map (PNG/PDF/SVG); shown only if None

    Returns
    -------
    (selection GeoDataFrame, step-by-step record counts DataFrame, matplotlib Figure or None)
    """
    steps = [("input", len(gdf))]
    sel = filter_by_region(gdf, region);          steps.append(("region", len(sel)))
    sel = filter_by_quality(sel, quality);        steps.append(("quality", len(sel)))
    sel = filter_by_indicator(sel, indicator);    steps.append(("indicator", len(sel)))
    sel = filter_by_regime(sel, regime);          steps.append(("regime", len(sel)))
    sel = filter_by_depth(sel, depth_range);      steps.append(("depth", len(sel)))
    if require_azimuth:
        sel = sel[sel["AZI"].notna()].copy();     steps.append(("has SHmax azimuth", len(sel)))
    counts = pd.DataFrame(steps, columns=["filter_step", "n_records"])

    if verbose:
        print(f"CASMO-inspired filter -> region={region!r}, quality={quality!r}, "
              f"indicator={indicator!r}, regime={regime!r}, depth={depth_range!r}")
        print(counts.to_string(index=False))

    fig = None
    if make_map:
        from .plotting import plot_stress_map  # local import keeps module light
        if isinstance(region, str) and region in REGION_PRESETS:
            extent = REGION_PRESETS[region]["map_extent"]
        elif isinstance(region, (tuple, list)) and len(region) == 4 and all(
                isinstance(v, (int, float)) for v in region):
            extent = tuple(region)
        else:
            b = sel.total_bounds if len(sel) else gdf.total_bounds
            pad = 0.5
            extent = (b[0] - pad, b[1] - pad, b[2] + pad, b[3] + pad)
        if title is None:
            title = f"WSM 2025 – SHmax orientations ({region if region else 'selection'})"
        subtitle = (f"Quality: {', '.join(_resolve_quality(quality))} | "
                    f"Indicator: {indicator if indicator else 'all'} | "
                    f"Regime: {regime if regime else 'all'} | n = {len(sel)}")
        highlight = "Germany"
        if isinstance(region, str) and REGION_PRESETS.get(region, {}).get("kind") == "polygon":
            highlight = REGION_PRESETS[region]["country_name"]
        fig = plot_stress_map(sel, extent=extent, title=title, subtitle=subtitle,
                              out_path=out_path, highlight=highlight)
    return sel, counts, fig


def region_polygon(name: str):
    """Return the shapely geometry of a preset region (polygon or bbox) in EPSG:4326."""
    spec = REGION_PRESETS[name]
    if spec["kind"] == "bbox":
        return box(*spec["bbox"])
    countries = load_boundaries()["countries"]
    return countries.loc[countries["NAME"] == spec["country_name"], "geometry"].union_all()
