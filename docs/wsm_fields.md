# Notes on the WSM 2025 fields

Definitions from WSM Technical Report 25-01 (https://doi.org/10.48440/wsm.2025.001),
plus what I found when inspecting the file.

| Field | Meaning | Notes |
|---|---|---|
| `ID` | WSM record ID | unique |
| `LAT`, `LON` | coordinates (WGS84, degrees) | 170 records without coordinates (all quality Xmi) |
| `AZI` | SHmax azimuth, clockwise from north, 0–180° | `999` = no orientation (11 048 records, mostly E) |
| `TYPE` | stress indicator | FMS, FMF, FMA, BO, DIF, HF, HFH, HFS, OC, GFI, GFM, GVA, SWB, SWL, SWS, BS |
| `QUALITY` | quality class | A, B, C, D, E, Xmi, Xne, Xru |
| `REGIME` | stress regime | NF, NS, SS, TS, TF, U |
| `DEPTH` | depth in km | |
| `COUNTRY` | country name or `Sea` | |
| `REF1`…`REF6` | references | |
| `SD`, `NUMBER`, `TOT_LEN` | standard deviation, number and total length of observations | used for the quality ranking |
| `PLATE`, `DIST` | plate (PB2002) and distance to plate boundary (km) | not used |

## Quality classes

| Class | SHmax accuracy |
|---|---|
| A | ±15° |
| B | ±20° |
| C | ±25° |
| D | ±40° |
| E | > 40° or no reliable orientation |
| Xmi / Xne / Xru | new in 2025: missing information / not-established indicator / rarely used indicator |

The technical report describes A–C as generally reliable for stress-pattern analysis, so
only A–C are used for the maps and statistics.

## Indicator types (short)

- **FMS / FMF / FMA** – earthquake focal mechanisms (single event / inversion of many events / average)
- **BO** – borehole breakouts (perpendicular to SHmax)
- **DIF** – drilling-induced tensile fractures (parallel to SHmax)
- **HF** – hydraulic fracturing
- **OC** – overcoring (shallow strain relief)
- **GFI / GFM** – geological fault-slip data, **GVA** – volcanic vent alignment
- **SW\*** – seismic anisotropy, **BS** – borehole slotter

## Stress regimes

- NF (normal faulting): S<sub>V</sub> > SHmax > Shmin
- SS (strike-slip): SHmax > S<sub>V</sub> > Shmin
- TF (thrust faulting): SHmax > Shmin > S<sub>V</sub>
- NS, TS: transitional · U: unknown
