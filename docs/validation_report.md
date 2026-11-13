# Validation report

_Generated automatically by `run_analysis.py` on 2026-11-13._

**25 of 25 automated checks passed.**

## Automated checks

| # | Check | Result | Detail |
|---|---|---|---|
| 1 | Raw file checksum matches documented SHA-256 | ✅ PASS | 3b24aa5b33b251b6... |
| 2 | Expected key columns present | ✅ PASS |  |
| 3 | All latitudes within [-90, 90] | ✅ PASS |  |
| 4 | All longitudes within [-180, 180] | ✅ PASS |  |
| 5 | All azimuths within [0, 180) or missing | ✅ PASS | 11048 records had AZI=999 (no orientation) |
| 6 | No duplicated WSM IDs | ✅ PASS |  |
| 7 | GeoDataFrame CRS is EPSG:4326 | ✅ PASS |  |
| 8 | Row bookkeeping: raw = clean + removed coordinates | ✅ PASS |  |
| 9 | Central Europe subset lies inside its bounding box | ✅ PASS |  |
| 10 | A-C subset contains only A, B, C | ✅ PASS |  |
| 11 | Every A-C record in Central Europe has an SHmax azimuth | ✅ PASS |  |
| 12 | Axial mean of 170° and 10° is 0° (N-S), not 90° | ✅ PASS |  |
| 13 | Axial mean treats 0° and 180° as identical (R = 1) | ✅ PASS |  |
| 14 | Uniform random azimuths give small R (<0.05, n=5000) | ✅ PASS |  |
| 15 | Map symbols reproduce SHmax azimuth (max error < 0.5°) | ✅ PASS | max deviation 0.230° |
| 16 | EPSG:4326 -> EPSG:3035 -> EPSG:4326 round trip (Potsdam) < 1e-7° | ✅ PASS |  |
| 17 | Export Central Europe (all): CSV/GPKG/GeoJSON contain exactly the filtered records | ✅ PASS | n = 2216 |
| 18 | Export Central Europe (all): GeoPackage CRS is EPSG:4326 | ✅ PASS |  |
| 19 | Export Central Europe (all): azimuths unchanged in CSV | ✅ PASS |  |
| 20 | Export Central Europe A-C: CSV/GPKG/GeoJSON contain exactly the filtered records | ✅ PASS | n = 1019 |
| 21 | Export Central Europe A-C: GeoPackage CRS is EPSG:4326 | ✅ PASS |  |
| 22 | Export Central Europe A-C: azimuths unchanged in CSV | ✅ PASS |  |
| 23 | Export Germany A-C: CSV/GPKG/GeoJSON contain exactly the filtered records | ✅ PASS | n = 297 |
| 24 | Export Germany A-C: GeoPackage CRS is EPSG:4326 | ✅ PASS |  |
| 25 | Export Germany A-C: azimuths unchanged in CSV | ✅ PASS |  |

## Cleaning log (global WSM 2025 file)

- Raw records: **100842**
- Records with missing coordinates removed: **170** (all quality Xmi); out-of-range coordinates: 0
- `AZI = 999` (WSM code for *no orientation provided*) converted to missing: **11048** records, by quality: {'E': 10996, 'Xmi': 36, 'Xru': 8, 'Xne': 7, 'B': 1}
  - Note: one record with quality B has AZI = 999 (wsm250689, FMF, offshore at 40.4°N, 0.69°E). It is kept in the counts but cannot be mapped; it lies outside the study area.
- Azimuths of exactly 180° set to 0° (same axis): 243
- Duplicated IDs: 0
- Records sharing identical location, azimuth, type, depth and quality with another record: 1629 globally (32 in Central Europe, 5 in the A–C subset). They have different WSM IDs and may be separate observations reported with rounded values; without checking the original references they cannot be confirmed as duplicates, so they were **kept** and are documented as a limitation.
- Clean records: **100672**, of which with SHmax azimuth: **89624**

## Records before / after filters

| step                                       |   n_records |
|:-------------------------------------------|------------:|
| WSM 2025 raw file                          |      100842 |
| after removing missing/invalid coordinates |      100672 |
| ...of which with an SHmax azimuth          |       89624 |
| inside Central Europe window               |        2216 |
| Central Europe, quality A-C                |        1019 |
| inside Germany (Natural Earth polygon)     |         791 |
| Germany, quality A-C                       |         297 |

## Missing values in the Central Europe subset (all qualities)

|         |   n_missing |   pct |
|:--------|------------:|------:|
| LAT     |           0 |   0   |
| LON     |           0 |   0   |
| AZI     |         211 |   9.5 |
| TYPE    |           0 |   0   |
| DEPTH   |           0 |   0   |
| QUALITY |           0 |   0   |
| REGIME  |           0 |   0   |
| COUNTRY |           0 |   0   |
| REF1    |           1 |   0   |
| DATE    |         981 |  44.3 |
| SD      |        1366 |  61.6 |
| NUMBER  |        1246 |  56.2 |

Records without an SHmax azimuth in Central Europe: 211, by quality: {'E': 204, 'Xne': 7}. No A–C record lacks an azimuth.

## Germany: polygon selection vs. WSM `COUNTRY` field

- Inside the Natural Earth Germany polygon: 791 records; `COUNTRY == 'Germany'`: 786 records.
- 13 records inside the polygon carry another country name ({'Switzerland': 6, 'Netherlands': 4, 'Belgium': 1, 'Austria': 1, 'Luxembourg': 1}) and 8 records labelled Germany fall outside the polygon. All of them lie within about 2 km of the Natural Earth boundary or coastline, so the differences come from the generalised 1:10m boundary (2.7 % of the Germany records). The spatial (polygon) selection is used in this project.

## Coordinate reference systems

- WSM coordinates are geographic latitude/longitude (WGS84) → GeoDataFrame CRS set to **EPSG:4326**.
- Maps are drawn in **ETRS89 / LAEA Europe (EPSG:3035)**, an equal-area projection recommended for pan-European mapping. Round-trip transformation test passed.
- SHmax symbols: end points computed geodesically on the WGS84 ellipsoid, then projected. Recomputing the azimuth from the projected symbols reproduces the database azimuth (max deviation 0.230°).
- Visual test: `data/processed/symbol_orientation_test.png` (synthetic test points at 0°, 45°, 90°, 135°) was inspected; symbols follow the local graticule as expected.

## Axial statistics checks

- Known-answer tests: mean of 170° and 10° = 0° (correct), not 90° (arithmetic mean); 0° and 180° give R = 1; uniform random azimuths give R ≈ 0.
- Central Europe A–C: axial mean 154.9° vs. naive arithmetic mean 116.2° – the difference shows why the arithmetic mean is not used.
