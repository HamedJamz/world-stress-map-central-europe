"""
Load, clean and export the World Stress Map (WSM) Database Release 2025.

Data: Heidbach et al. (2025), World Stress Map Database Release 2025,
GFZ Data Services, https://doi.org/10.5880/WSM.2025.001 (CC BY 4.0).
Field definitions: WSM Technical Report 25-01, https://doi.org/10.48440/wsm.2025.001
"""

from __future__ import annotations

import hashlib
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd

# --------------------------------------------------------------------------------------
# Paths and constants
# --------------------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
OUTPUT_DIR = PROJECT_ROOT / "outputs"
NE_DIR = RAW_DIR / "natural_earth"

WSM_FILENAME = "WSM_Database_2025.csv"
# Direct file link behind the DOI landing page https://doi.org/10.5880/WSM.2025.001
# (resolved 2026-11-11; the GFZ file server path may change in future - then use the DOI).
WSM_URL = ("https://datapub.gfz.de/download/10.5880.WSM.2025.001-Scbwez/"
           "WSM_Database_2025.csv")
WSM_SHA256 = "3b24aa5b33b251b6c52ee25bcafe3d23d036deb1c38e1caae51f63b808605d1b"

# Geographic coordinate reference system of the WSM file (WGS84 lat/lon in degrees)
CRS_GEOGRAPHIC = "EPSG:4326"
# Equal-area projected CRS used for European maps (ETRS89 / LAEA Europe, metres)
CRS_MAP = "EPSG:3035"

# Quality classes as defined in the WSM 2025 quality ranking scheme (TR 25-01)
QUALITY_ORDER = ["A", "B", "C", "D", "E", "Xmi", "Xne", "Xru"]
RELIABLE_QUALITIES = ["A", "B", "C"]          # "considered reliable" (TR 25-01, Sect. 3.2)

# Stress-indicator TYPE codes (TR 25-01, Tab. 2-1) -> broader indicator group
INDICATOR_GROUPS = {
    "FMS": "Focal mechanism", "FMF": "Focal mechanism", "FMA": "Focal mechanism",
    "BO": "Borehole breakout", "BS": "Borehole slotter",
    "DIF": "Drilling-induced fracture",
    "HF": "Hydraulic fracturing", "HFH": "Hydraulic fracturing",
    "HFP": "Hydraulic fracturing", "HFS": "Hydraulic fracturing",
    "OC": "Overcoring",
    "GFI": "Geological fault-slip", "GFM": "Geological fault-slip",
    "GVA": "Volcanic vent alignment",
    "SWB": "Seismic anisotropy", "SWL": "Seismic anisotropy", "SWS": "Seismic anisotropy",
}

TYPE_DESCRIPTIONS = {
    "FMS": "Single earthquake focal mechanism",
    "FMF": "Formal inversion of several focal mechanisms",
    "FMA": "Average of focal mechanisms",
    "BO": "Borehole breakouts",
    "BS": "Borehole slotter",
    "DIF": "Drilling-induced tensile fractures",
    "HF": "Hydraulic fracturing (sub-vertical wells)",
    "HFH": "Hydraulic fracturing (sub-horizontal wells)",
    "HFP": "Hydraulic testing of pre-existing fractures",
    "HFS": "Seismicity on hydraulic fractures",
    "OC": "Overcoring / strain relief",
    "GFI": "Inversion of geological fault-slip data",
    "GFM": "Geological paleo-focal mechanism",
    "GVA": "Volcanic vent alignment",
    "SWB": "Seismic anisotropy in boreholes (sonic logs)",
    "SWL": "Seismic anisotropy in laboratory",
    "SWS": "Seismic anisotropy from seismology",
}

REGIME_NAMES = {
    "NF": "Normal faulting", "NS": "Normal / strike-slip",
    "SS": "Strike-slip", "TS": "Thrust / strike-slip",
    "TF": "Thrust faulting", "U": "Unknown",
}


# --------------------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------------------
def sha256_of(path: Path) -> str:
    """Return the SHA-256 checksum of a file (used to document the exact raw file)."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download_wsm_data(dest: Path | None = None, url: str = WSM_URL) -> Path:
    """Download the WSM 2025 CSV from GFZ Data Services if it is not already present.

    If the automatic download fails (e.g. the file-server path changed), download the
    file manually from https://doi.org/10.5880/WSM.2025.001 into ``data/raw/``.
    """
    dest = Path(dest) if dest else RAW_DIR / WSM_FILENAME
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading WSM data from {url} ...")
    urllib.request.urlretrieve(url, dest)
    return dest


def load_wsm_data(path: Path | None = None) -> pd.DataFrame:
    """Load the raw WSM CSV into a pandas DataFrame without modifying any values.

    The file is UTF-8 with a byte-order mark, so ``utf-8-sig`` is used to obtain a
    clean first column name (``ID``).
    """
    path = Path(path) if path else download_wsm_data()
    df = pd.read_csv(path, encoding="utf-8-sig", low_memory=False)
    return df


def describe_schema(df: pd.DataFrame) -> pd.DataFrame:
    """Summarise every column: dtype, number/percentage of missing values, example value."""
    rows = []
    for col in df.columns:
        s = df[col]
        example = s.dropna().iloc[0] if s.notna().any() else None
        rows.append({
            "column": col,
            "dtype": str(s.dtype),
            "n_missing": int(s.isna().sum()),
            "pct_missing": round(100 * s.isna().mean(), 1),
            "n_unique": int(s.nunique(dropna=True)),
            "example": example,
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------
# Cleaning
# --------------------------------------------------------------------------------------
def clean_wsm_data(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Clean and standardise the key WSM fields.

    Steps (all documented in the returned ``log`` dictionary):
      1. strip whitespace in text fields used for filtering;
      2. remove records with missing or out-of-range coordinates;
      3. convert AZI = 999 (WSM code for "no orientation provided") to NaN;
      4. flag any remaining azimuths outside 0-180 deg as invalid;
      5. express azimuths in the half-open interval [0, 180) - 180 deg is the same
         horizontal axis as 0 deg (axial data);
      6. add helper columns (indicator group, regime name, reliable flag, has_azimuth).

    No records are removed because of missing SHmax; they are kept and flagged so the
    quality statistics stay complete. Mapping/statistics use only records with an azimuth.
    """
    log: dict = {"n_raw": len(df)}
    out = df.copy()

    for col in ["ID", "TYPE", "QUALITY", "REGIME", "COUNTRY", "LOCALITY"]:
        if col in out.columns:
            out[col] = out[col].astype("string").str.strip()

    # --- coordinates -----------------------------------------------------------------
    missing_xy = out["LAT"].isna() | out["LON"].isna()
    bad_xy = (~out["LAT"].between(-90, 90)) | (~out["LON"].between(-180, 180))
    bad_xy &= ~missing_xy
    log["n_missing_coordinates"] = int(missing_xy.sum())
    log["n_out_of_range_coordinates"] = int(bad_xy.sum())
    log["missing_coordinates_quality"] = (
        out.loc[missing_xy, "QUALITY"].value_counts().to_dict())
    out = out.loc[~(missing_xy | bad_xy)].copy()

    # --- SHmax azimuth -------------------------------------------------------------------
    out["AZI_RAW"] = out["AZI"]
    no_azi = out["AZI"] == 999
    log["n_azi_999_no_orientation"] = int(no_azi.sum())
    log["azi_999_quality"] = out.loc[no_azi, "QUALITY"].value_counts().to_dict()
    out["AZI"] = out["AZI"].astype(float).where(~no_azi, np.nan)

    invalid_azi = out["AZI"].notna() & ~out["AZI"].between(0, 180)
    log["n_azi_out_of_range"] = int(invalid_azi.sum())
    out.loc[invalid_azi, "AZI"] = np.nan

    log["n_azi_equal_180_converted_to_0"] = int((out["AZI"] == 180).sum())
    out["AZI"] = np.mod(out["AZI"], 180.0)          # 180 -> 0, axial convention

    # --- helper columns ----------------------------------------------------------------
    out["has_azimuth"] = out["AZI"].notna()
    out["indicator_group"] = out["TYPE"].map(INDICATOR_GROUPS).fillna("Other")
    out["regime_name"] = out["REGIME"].map(REGIME_NAMES).fillna("Unknown")
    out["is_reliable_AC"] = out["QUALITY"].isin(RELIABLE_QUALITIES)
    out["QUALITY"] = pd.Categorical(out["QUALITY"], categories=QUALITY_ORDER, ordered=True)

    # --- duplicates (reported, not removed) ----------------------------------------------
    log["n_duplicate_ids"] = int(out["ID"].duplicated().sum())
    key = ["LAT", "LON", "AZI", "TYPE", "DEPTH", "QUALITY"]
    log["n_identical_location_azi_type_depth_quality"] = int(out.duplicated(key).sum())

    log["n_clean"] = len(out)
    log["n_clean_with_azimuth"] = int(out["has_azimuth"].sum())
    return out.reset_index(drop=True), log


def to_geodataframe(df: pd.DataFrame) -> gpd.GeoDataFrame:
    """Convert a cleaned DataFrame to a GeoDataFrame with point geometry in WGS84."""
    gdf = gpd.GeoDataFrame(
        df.copy(),
        geometry=gpd.points_from_xy(df["LON"], df["LAT"]),
        crs=CRS_GEOGRAPHIC,
    )
    return gdf


# --------------------------------------------------------------------------------------
# Boundaries (Natural Earth, public domain)
# --------------------------------------------------------------------------------------
def load_boundaries() -> dict[str, gpd.GeoDataFrame]:
    """Load the Natural Earth 1:10m layers (clipped to Europe) used as basemap.

    Source: Natural Earth (public domain), https://www.naturalearthdata.com,
    obtained from https://github.com/nvkelso/natural-earth-vector (GeoJSON).
    """
    layers = {
        "countries": "ne_10m_admin_0_countries_europe.gpkg",
        "admin1_lines": "ne_10m_admin_1_lines_europe.gpkg",
        "rivers": "ne_10m_rivers_europe.gpkg",
        "lakes": "ne_10m_lakes_europe.gpkg",
    }
    return {k: gpd.read_file(NE_DIR / v) for k, v in layers.items()}


# --------------------------------------------------------------------------------------
# Summaries
# --------------------------------------------------------------------------------------
def quality_counts(df: pd.DataFrame) -> pd.Series:
    """Number of records per WSM quality class, in A..E, X order."""
    return df["QUALITY"].value_counts().reindex(QUALITY_ORDER, fill_value=0)


def indicator_counts(df: pd.DataFrame, column: str = "TYPE") -> pd.Series:
    """Number of records per stress-indicator type (or group), sorted descending."""
    return df[column].value_counts()


def summary_by_indicator(df: pd.DataFrame) -> pd.DataFrame:
    """Descriptive table per indicator TYPE: counts, A-C counts and depth statistics."""
    g = df.groupby("TYPE", observed=True)
    tab = pd.DataFrame({
        "description": [TYPE_DESCRIPTIONS.get(t, "") for t in g.size().index],
        "n_records": g.size(),
        "n_A_to_C": g["is_reliable_AC"].sum().astype(int),
        "n_with_azimuth": g["has_azimuth"].sum().astype(int),
        "median_depth_km": g["DEPTH"].median().round(2),
        "max_depth_km": g["DEPTH"].max().round(2),
    })
    return tab.sort_values("n_records", ascending=False)


def export_dataset(gdf: gpd.GeoDataFrame, stem: str, out_dir: Path | None = None) -> dict:
    """Export a GeoDataFrame to CSV, GeoPackage and GeoJSON and return the file paths."""
    out_dir = Path(out_dir) if out_dir else OUTPUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    g = gdf.copy()
    g["QUALITY"] = g["QUALITY"].astype("string")
    for c in g.columns:          # nullable string dtype -> plain object for writers
        if str(g[c].dtype) == "string":
            g[c] = g[c].astype(object)
    paths = {
        "csv": out_dir / f"{stem}.csv",
        "gpkg": out_dir / f"{stem}.gpkg",
        "geojson": out_dir / f"{stem}.geojson",
    }
    g.drop(columns="geometry").to_csv(paths["csv"], index=False)
    g.to_file(paths["gpkg"], layer="wsm_records", driver="GPKG")
    g.to_file(paths["geojson"], driver="GeoJSON")
    return paths
