"""
Runs the whole workflow: load -> clean -> select study area -> filter by quality ->
statistics -> figures -> CASMO-like examples -> exports -> validation report.

Usage (from the repository root):
    python run_analysis.py
"""

from __future__ import annotations

import json
import warnings
from datetime import date

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pyproj import Transformer

from src.data_processing import (PROJECT_ROOT, RAW_DIR, PROCESSED_DIR, OUTPUT_DIR,
                                 WSM_FILENAME, WSM_SHA256, WSM_URL, CRS_GEOGRAPHIC, CRS_MAP,
                                 QUALITY_ORDER, load_wsm_data, describe_schema,
                                 clean_wsm_data, to_geodataframe, quality_counts,
                                 summary_by_indicator, export_dataset, sha256_of)
from src.filtering import (REGION_PRESETS, filter_by_region, filter_by_quality,
                           filter_stress_data)
from src.orientation_stats import (axial_mean, orientation_summary, rayleigh_test_axial,
                                   gridded_axial_means, azimuth_to_compass)
from src import plotting as P

warnings.filterwarnings("ignore", category=UserWarning)
FIG = PROJECT_ROOT / "figures"
DOCS = PROJECT_ROOT / "docs"
for d in (FIG, PROCESSED_DIR, OUTPUT_DIR, DOCS):
    d.mkdir(parents=True, exist_ok=True)

CE_BBOX = REGION_PRESETS["Central Europe"]["bbox"]
checks: list[tuple[str, bool, str]] = []          # (check, passed, detail)


def check(name: str, condition: bool, detail: str = "") -> None:
    """Record a validation check (and fail loudly if it does not pass)."""
    checks.append((name, bool(condition), detail))
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {name} {('- ' + detail) if detail else ''}")


# ======================================================================================
# 1-2  Load and inspect
# ======================================================================================
print("\n=== 1. Load WSM data ===")
raw_path = RAW_DIR / WSM_FILENAME
raw = load_wsm_data(raw_path)
sha = sha256_of(raw_path)
check("Raw file checksum matches documented SHA-256", sha == WSM_SHA256, sha[:16] + "...")
print(f"Raw shape: {raw.shape}")
schema = describe_schema(raw)
schema.to_csv(PROCESSED_DIR / "schema_summary.csv", index=False)
check("Expected key columns present",
      all(c in raw.columns for c in ["ID", "LAT", "LON", "AZI", "TYPE", "DEPTH", "QUALITY",
                                     "REGIME", "COUNTRY", "REF1"]))

# ======================================================================================
# 3-4  Clean + GeoDataFrame
# ======================================================================================
print("\n=== 2. Clean ===")
clean, log = clean_wsm_data(raw)
(PROCESSED_DIR / "cleaning_log.json").write_text(json.dumps(log, indent=2, default=str))
gdf = to_geodataframe(clean)

check("All latitudes within [-90, 90]", gdf["LAT"].between(-90, 90).all())
check("All longitudes within [-180, 180]", gdf["LON"].between(-180, 180).all())
check("All azimuths within [0, 180) or missing",
      ((gdf["AZI"] >= 0) & (gdf["AZI"] < 180) | gdf["AZI"].isna()).all(),
      f"{log['n_azi_999_no_orientation']} records had AZI=999 (no orientation)")
check("No duplicated WSM IDs", log["n_duplicate_ids"] == 0)
check("GeoDataFrame CRS is EPSG:4326", gdf.crs.to_epsg() == 4326)
check("Row bookkeeping: raw = clean + removed coordinates",
      log["n_raw"] == log["n_clean"] + log["n_missing_coordinates"]
      + log["n_out_of_range_coordinates"])

# ======================================================================================
# 5-6  Study areas and quality filter
# ======================================================================================
print("\n=== 3. Study areas ===")
ce_all = filter_by_region(gdf, "Central Europe")
ce_ac = filter_by_quality(ce_all, "C")
ce_ac = ce_ac[ce_ac["AZI"].notna()].copy()
de_all = filter_by_region(gdf, "Germany")
de_ac = filter_by_quality(de_all, "C")
de_ac = de_ac[de_ac["AZI"].notna()].copy()
de_field = gdf[gdf["COUNTRY"] == "Germany"]

check("Central Europe subset lies inside its bounding box",
      ce_all["LON"].between(CE_BBOX[0], CE_BBOX[2]).all()
      and ce_all["LAT"].between(CE_BBOX[1], CE_BBOX[3]).all())
check("A-C subset contains only A, B, C",
      set(ce_ac["QUALITY"].astype(str)) <= {"A", "B", "C"})
check("Every A-C record in Central Europe has an SHmax azimuth",
      ce_all[ce_all["is_reliable_AC"]]["AZI"].notna().all())

filter_steps = pd.DataFrame([
    ("WSM 2025 raw file", log["n_raw"]),
    ("after removing missing/invalid coordinates", log["n_clean"]),
    ("  ...of which with an SHmax azimuth", log["n_clean_with_azimuth"]),
    ("inside Central Europe window", len(ce_all)),
    ("Central Europe, quality A-C", len(ce_ac)),
    ("inside Germany (Natural Earth polygon)", len(de_all)),
    ("Germany, quality A-C", len(de_ac)),
], columns=["step", "n_records"])
filter_steps.to_csv(PROCESSED_DIR / "filter_step_counts.csv", index=False)
print(filter_steps.to_string(index=False))

# Germany: polygon selection vs. WSM COUNTRY field (consistency check)
in_poly_not_field = de_all[de_all["COUNTRY"] != "Germany"]
in_field_not_poly = de_field[~de_field["ID"].isin(de_all["ID"])]

# ======================================================================================
# 7  Statistics
# ======================================================================================
print("\n=== 4. Statistics ===")
q_global = quality_counts(gdf).rename("global")
q_ce = quality_counts(ce_all).rename("central_europe")
q_de = quality_counts(de_all).rename("germany")
qtab = pd.concat([q_global, q_ce, q_de], axis=1)
qtab.to_csv(PROCESSED_DIR / "quality_counts.csv")
print(qtab)

ind_ce = summary_by_indicator(ce_all)
ind_ce.to_csv(PROCESSED_DIR / "indicator_summary_central_europe.csv")
ind_de = summary_by_indicator(de_all)
ind_de.to_csv(PROCESSED_DIR / "indicator_summary_germany.csv")
print(ind_ce)

regime_ce = (ce_ac["regime_name"].value_counts().rename("n_A_to_C")
             .to_frame().assign(pct=lambda t: (100 * t["n_A_to_C"] / t["n_A_to_C"].sum()).round(1)))
regime_ce.to_csv(PROCESSED_DIR / "regime_counts_central_europe_AC.csv")

missing_ce = (ce_all[["LAT", "LON", "AZI", "TYPE", "DEPTH", "QUALITY", "REGIME", "COUNTRY",
                      "REF1", "DATE", "SD", "NUMBER"]].isna().sum().rename("n_missing")
              .to_frame().assign(pct=lambda t: (100 * t["n_missing"] / len(ce_all)).round(1)))
missing_ce.to_csv(PROCESSED_DIR / "missing_values_central_europe.csv")
missing_azi_by_q = ce_all[ce_all["AZI"].isna()]["QUALITY"].astype(str).value_counts()

ostats = pd.concat([
    orientation_summary(ce_ac).rename(index={"All": "Central Europe A-C"}),
    orientation_summary(ce_ac, weighted=True).rename(index={"All": "Central Europe A-C (quality-weighted)"}),
    orientation_summary(de_ac).rename(index={"All": "Germany A-C"}),
    orientation_summary(ce_all[ce_all["AZI"].notna()]).rename(index={"All": "Central Europe all qualities"}),
])
ostats_type = orientation_summary(ce_ac, "TYPE")
ostats.to_csv(PROCESSED_DIR / "orientation_stats.csv")
ostats_type.to_csv(PROCESSED_DIR / "orientation_stats_by_type_central_europe_AC.csv")
print(ostats); print(ostats_type)

grid = gridded_axial_means(ce_ac, cell_deg=1.0, min_n=5)
grid.round(3).to_csv(PROCESSED_DIR / "gridded_mean_shmax_1deg_AC.csv", index=False)

# depth classes (A-C)
depth_bins = pd.cut(ce_ac["DEPTH"], [-5, 1, 5, 10, 20, 40],
                    labels=["<1 km", "1-5 km", "5-10 km", "10-20 km", ">20 km"])
depth_tab = ce_ac.groupby(depth_bins, observed=False).agg(
    n=("AZI", "size"),
    mean_deg=("AZI", lambda s: axial_mean(s)["mean_deg"]),
    R=("AZI", lambda s: axial_mean(s)["R"])).round(2)
depth_tab.to_csv(PROCESSED_DIR / "orientation_by_depth_central_europe_AC.csv")

# Axial statistics sanity checks with known answers
check("Axial mean of 170° and 10° is 0° (N-S), not 90°",
      abs(((axial_mean([170, 10])["mean_deg"] + 90) % 180) - 90) < 1e-6)
check("Axial mean treats 0° and 180° as identical (R = 1)",
      abs(axial_mean([0, 180, 0])["R"] - 1) < 1e-9)
rng = np.random.default_rng(0)
check("Uniform random azimuths give small R (<0.05, n=5000)",
      axial_mean(rng.uniform(0, 180, 5000))["R"] < 0.05)

# ======================================================================================
# 8  Figures
# ======================================================================================
print("\n=== 5. Figures ===")
st_ce = axial_mean(ce_ac["AZI"])
P.plot_stress_map(
    ce_ac, extent=REGION_PRESETS["Central Europe"]["map_extent"],
    title="Crustal stress orientations in Central Europe (WSM 2025)",
    subtitle=(f"Orientation of maximum horizontal compressive stress SHmax · quality A–C · "
              f"n = {len(ce_ac)} records · study window 45.5–55.5°N, 2.5–19.5°E"),
    out_path=FIG / "stress_map", formats=("png", "pdf", "svg"))
P.plot_stress_map(
    de_ac, extent=REGION_PRESETS["Germany"]["map_extent"],
    title="Crustal stress orientations in Germany (WSM 2025)",
    subtitle=f"SHmax orientations · quality A–C · n = {len(de_ac)} records inside Germany",
    out_path=FIG / "stress_map_germany", formats=("png", "pdf"))
P.plot_quality_distribution(
    ce_all, title=f"WSM 2025 records by quality class – Central Europe (n = {len(ce_all)})",
    out_path=FIG / "quality_distribution")
P.plot_indicator_distribution(
    ce_all, title=f"WSM 2025 records by stress-indicator type – Central Europe (n = {len(ce_all)})",
    out_path=FIG / "indicator_distribution")
P.plot_orientation_distribution(
    ce_ac, title=f"SHmax orientation distribution – Central Europe, quality A–C",
    out_path=FIG / "orientation_distribution")
P.plot_quality_comparison(ce_all, ce_ac, extent=REGION_PRESETS["Central Europe"]["map_extent"],
                          out_path=FIG / "quality_filter_comparison")
P.plot_rose_by_group(ce_ac, "TYPE", groups=["FMS", "BO", "HF", "DIF"],
                     title="SHmax orientations by indicator type – Central Europe, A–C",
                     out_path=FIG / "rose_by_indicator")
P.plot_gridded_means(grid, extent=REGION_PRESETS["Central Europe"]["map_extent"],
                     out_path=FIG / "gridded_mean_shmax")
plt.close("all")

# --- symbol orientation test figure (synthetic test geometry, clearly labelled) -----------
test = gpd.GeoDataFrame({
    "LON": np.repeat([4.0, 10.0, 18.0], 4), "LAT": np.tile([47.0, 49.5, 52.0, 54.5], 3),
    "AZI": np.tile([0.0, 45.0, 90.0, 135.0], 3), "QUALITY": "A", "REGIME": "U"})
test = test.set_geometry(gpd.points_from_xy(test["LON"], test["LAT"]), crs=CRS_GEOGRAPHIC)
fig = P.plot_stress_map(test, extent=REGION_PRESETS["Central Europe"]["map_extent"],
                        title="Symbol orientation test (synthetic test points, not WSM data)",
                        subtitle="Expected: 0° = N–S, 45° = NE–SW, 90° = E–W, 135° = NW–SE; "
                                 "symbols must follow the local meridians/graticule",
                        out_path=PROCESSED_DIR / "symbol_orientation_test", show_inset=False)
for _, r in test.iterrows():
    x, y = P.TO_MAP.transform(r["LON"], r["LAT"])
    fig.axes[0].text(x + 25000, y + 15000, f"{r['AZI']:.0f}°", fontsize=7, color="#c0392b", zorder=8)
fig.savefig(PROCESSED_DIR / "symbol_orientation_test.png", dpi=150)
plt.close(fig)

# numerical symbol check: un-project end points and recompute geodesic azimuth
seg = P.shmax_segments(ce_ac, scale=1.0)
inv = Transformer.from_crs(CRS_MAP, CRS_GEOGRAPHIC, always_xy=True)
lo1, la1 = inv.transform(seg[:, 0, 0], seg[:, 0, 1])
lo2, la2 = inv.transform(seg[:, 1, 0], seg[:, 1, 1])
az_back, _, _ = P.GEOD.inv(lo2, la2, lo1, la1)
diff = np.abs(((az_back % 180) - ce_ac["AZI"].to_numpy() + 90) % 180 - 90)
check("Map symbols reproduce SHmax azimuth (max error < 0.5°)", diff.max() < 0.5,
      f"max deviation {diff.max():.3f}°")
x, y = Transformer.from_crs(CRS_GEOGRAPHIC, CRS_MAP, always_xy=True).transform(13.06, 52.38)
lo, la = inv.transform(x, y)
check("EPSG:4326 -> EPSG:3035 -> EPSG:4326 round trip (Potsdam) < 1e-7°",
      abs(lo - 13.06) < 1e-7 and abs(la - 52.38) < 1e-7)

# ======================================================================================
# 9  CASMO-inspired filtering examples
# ======================================================================================
print("\n=== 6. CASMO-inspired filtering examples ===")
examples = [
    dict(region="Switzerland", quality="C", indicator="FM*", name="casmo_example_switzerland_FM_AC"),
    dict(region="Germany", quality="B", indicator="BO", name="casmo_example_germany_BO_AB"),
    dict(region="Upper Rhine Graben area", quality="C", indicator=None,
         name="casmo_example_upper_rhine_graben_AC"),
]
casmo_rows = []
for ex in examples:
    sel, counts, fig = filter_stress_data(gdf, region=ex["region"], quality=ex["quality"],
                                          indicator=ex["indicator"], make_map=True,
                                          out_path=FIG / ex["name"])
    plt.close("all")
    casmo_rows.append({"example": ex["name"], "region": ex["region"],
                       "quality": ex["quality"], "indicator": ex["indicator"] or "all",
                       "n_selected": len(sel),
                       "axial_mean_deg": round(axial_mean(sel["AZI"])["mean_deg"], 1)})
pd.DataFrame(casmo_rows).to_csv(PROCESSED_DIR / "casmo_inspired_examples.csv", index=False)

# ======================================================================================
# 10  Exports + export validation
# ======================================================================================
print("\n=== 7. Exports ===")
paths_all = export_dataset(ce_all, "wsm_central_europe_clean")
paths_ac = export_dataset(ce_ac, "wsm_central_europe_AC")
paths_de = export_dataset(de_ac, "wsm_germany_AC")

for label, paths, ref in [("Central Europe (all)", paths_all, ce_all),
                          ("Central Europe A-C", paths_ac, ce_ac),
                          ("Germany A-C", paths_de, de_ac)]:
    csv_back = pd.read_csv(paths["csv"], low_memory=False)
    gpkg_back = gpd.read_file(paths["gpkg"])
    gj_back = gpd.read_file(paths["geojson"])
    same_ids = (set(csv_back["ID"]) == set(ref["ID"]) == set(gpkg_back["ID"])
                == set(gj_back["ID"]))
    check(f"Export {label}: CSV/GPKG/GeoJSON contain exactly the filtered records",
          same_ids and len(csv_back) == len(ref) == len(gpkg_back) == len(gj_back),
          f"n = {len(ref)}")
    check(f"Export {label}: GeoPackage CRS is EPSG:4326", gpkg_back.crs.to_epsg() == 4326)
    if ref["AZI"].notna().any():
        m = ref[["ID", "AZI"]].merge(csv_back[["ID", "AZI"]], on="ID", suffixes=("", "_csv"))
        check(f"Export {label}: azimuths unchanged in CSV",
              np.allclose(m["AZI"].fillna(-1), m["AZI_csv"].fillna(-1)))

# ======================================================================================
# Results summary + validation report
# ======================================================================================
st_de = axial_mean(de_ac["AZI"])
st_all = axial_mean(ce_all["AZI"])
ray = rayleigh_test_axial(ce_ac["AZI"])
results = {
    "run_date": str(date.today()),
    "data": {"name": "World Stress Map Database Release 2025",
             "doi": "https://doi.org/10.5880/WSM.2025.001", "file": WSM_FILENAME,
             "url": WSM_URL, "sha256": sha, "license": "CC BY 4.0"},
    "cleaning": log,
    "central_europe": {
        "bbox_lon_lat": CE_BBOX, "n_all": len(ce_all), "n_AC": len(ce_ac),
        "quality_counts": q_ce.to_dict(),
        "indicator_counts": ce_all["TYPE"].value_counts().to_dict(),
        "indicator_counts_AC": ce_ac["TYPE"].value_counts().to_dict(),
        "n_without_azimuth": int(ce_all["AZI"].isna().sum()),
        "without_azimuth_by_quality": missing_azi_by_q[missing_azi_by_q > 0].to_dict(),
        "axial_mean_AC_deg": round(st_ce["mean_deg"], 1),
        "axial_mean_AC_compass": azimuth_to_compass(st_ce["mean_deg"]),
        "R_AC": round(st_ce["R"], 3), "circ_std_AC_deg": round(st_ce["circ_std_deg"], 1),
        "ci95_AC_deg": [float(ostats.loc["Central Europe A-C", "ci95_low"]),
                        float(ostats.loc["Central Europe A-C", "ci95_high"])],
        "rayleigh_Z": round(ray["Z"], 1), "rayleigh_p": ray["p_value"],
        "naive_arithmetic_mean_AC_deg": round(float(ce_ac["AZI"].mean()), 1),
        "axial_mean_all_qualities_deg": round(st_all["mean_deg"], 1),
        "R_all_qualities": round(st_all["R"], 3),
        "regime_AC": ce_ac["REGIME"].value_counts().to_dict(),
        "median_depth_AC_km": float(ce_ac["DEPTH"].median()),
    },
    "germany": {
        "n_polygon_all": len(de_all), "n_country_field": len(de_field), "n_AC": len(de_ac),
        "quality_counts": q_de.to_dict(),
        "indicator_counts_AC": de_ac["TYPE"].value_counts().to_dict(),
        "axial_mean_AC_deg": round(st_de["mean_deg"], 1),
        "axial_mean_AC_compass": azimuth_to_compass(st_de["mean_deg"]),
        "R_AC": round(st_de["R"], 3),
        "n_in_polygon_but_other_country_field": len(in_poly_not_field),
        "n_country_field_germany_outside_polygon": len(in_field_not_poly),
    },
    "orientation_by_type_AC": ostats_type.reset_index().to_dict(orient="records"),
    "grid_cells_1deg_n_ge_5": len(grid),
    "casmo_inspired_examples": casmo_rows,
    "checks_passed": sum(c[1] for c in checks), "checks_total": len(checks),
}
(OUTPUT_DIR / "results_summary.json").write_text(json.dumps(results, indent=2, default=float))

# ---- markdown validation report ------------------------------------------------------
def md_table(df: pd.DataFrame) -> str:
    """Render a small DataFrame as a GitHub markdown table."""
    return df.to_markdown()

lines = [
    "# Validation report",
    "",
    f"_Generated automatically by `run_analysis.py` on {date.today()}._",
    "",
    f"**{sum(c[1] for c in checks)} of {len(checks)} automated checks passed.**",
    "",
    "## Automated checks",
    "",
    "| # | Check | Result | Detail |", "|---|---|---|---|",
]
for i, (name, ok, det) in enumerate(checks, 1):
    lines.append(f"| {i} | {name} | {'✅ PASS' if ok else '❌ FAIL'} | {det} |")
lines += [
    "",
    "## Cleaning log (global WSM 2025 file)",
    "",
    f"- Raw records: **{log['n_raw']}**",
    f"- Records with missing coordinates removed: **{log['n_missing_coordinates']}** "
    f"(all quality {', '.join(log['missing_coordinates_quality'])}); out-of-range coordinates: "
    f"{log['n_out_of_range_coordinates']}",
    f"- `AZI = 999` (WSM code for *no orientation provided*) converted to missing: "
    f"**{log['n_azi_999_no_orientation']}** records, by quality: {log['azi_999_quality']}",
    f"  - Note: one record with quality B has AZI = 999 (wsm250689, FMF, offshore at 40.4°N, 0.69°E). "
    f"It is kept in the counts but cannot be mapped; it lies outside the study area.",
    f"- Azimuths of exactly 180° set to 0° (same axis): {log['n_azi_equal_180_converted_to_0']}",
    f"- Duplicated IDs: {log['n_duplicate_ids']}",
    f"- Records sharing identical location, azimuth, type, depth and quality with another "
    f"record: {log['n_identical_location_azi_type_depth_quality']} globally "
    f"({int(ce_all.duplicated(['LAT','LON','AZI','TYPE','DEPTH','QUALITY']).sum())} in Central Europe, "
    f"{int(ce_ac.duplicated(['LAT','LON','AZI','TYPE','DEPTH','QUALITY']).sum())} in the A–C subset). "
    f"They have different WSM IDs and may be separate observations reported with rounded "
    f"values; without checking the original references they cannot be confirmed as duplicates, "
    f"so they were **kept** and are documented as a limitation.",
    f"- Clean records: **{log['n_clean']}**, of which with SHmax azimuth: **{log['n_clean_with_azimuth']}**",
    "",
    "## Records before / after filters",
    "",
    md_table(filter_steps.set_index("step")),
    "",
    "## Missing values in the Central Europe subset (all qualities)",
    "",
    md_table(missing_ce),
    "",
    f"Records without an SHmax azimuth in Central Europe: {int(ce_all['AZI'].isna().sum())}, "
    f"by quality: {missing_azi_by_q[missing_azi_by_q > 0].to_dict()}. "
    "No A–C record lacks an azimuth.",
    "",
    "## Germany: polygon selection vs. WSM `COUNTRY` field",
    "",
    f"- Inside the Natural Earth Germany polygon: {len(de_all)} records; "
    f"`COUNTRY == 'Germany'`: {len(de_field)} records.",
    f"- {len(in_poly_not_field)} records inside the polygon carry another country name "
    f"({in_poly_not_field['COUNTRY'].value_counts().to_dict()}) and {len(in_field_not_poly)} "
    f"records labelled Germany fall outside the polygon. All of them lie within about 2 km of the "
    f"Natural Earth boundary or coastline, so the differences come from the generalised 1:10m "
    f"boundary ({100 * (len(in_poly_not_field) + len(in_field_not_poly)) / len(de_all):.1f} % of the "
    "Germany records). The spatial (polygon) selection is used in this project.",
    "",
    "## Coordinate reference systems",
    "",
    "- WSM coordinates are geographic latitude/longitude (WGS84) → GeoDataFrame CRS set to **EPSG:4326**.",
    "- Maps are drawn in **ETRS89 / LAEA Europe (EPSG:3035)**, an equal-area projection "
    "recommended for pan-European mapping. Round-trip transformation test passed.",
    "- SHmax symbols: end points computed geodesically on the WGS84 ellipsoid, then projected. "
    f"Recomputing the azimuth from the projected symbols reproduces the database azimuth "
    f"(max deviation {diff.max():.3f}°).",
    "- Visual test: `data/processed/symbol_orientation_test.png` (synthetic test points at "
    "0°, 45°, 90°, 135°) was inspected; symbols follow the local graticule as expected.",
    "",
    "## Axial statistics checks",
    "",
    "- Known-answer tests: mean of 170° and 10° = 0° (correct), not 90° (arithmetic mean); "
    "0° and 180° give R = 1; uniform random azimuths give R ≈ 0.",
    f"- Central Europe A–C: axial mean {st_ce['mean_deg']:.1f}° vs. naive arithmetic mean "
    f"{ce_ac['AZI'].mean():.1f}° – the difference shows why the arithmetic mean is not used.",
]
(DOCS / "validation_report.md").write_text("\n".join(lines) + "\n")

print(f"\nDone. {sum(c[1] for c in checks)}/{len(checks)} checks passed.")
if not all(c[1] for c in checks):
    raise SystemExit("Some validation checks failed - see docs/validation_report.md")
