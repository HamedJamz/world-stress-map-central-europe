"""
Statistics for SHmax orientations.

SHmax is an axis, so 0 deg and 180 deg are the same orientation (170 and 10 deg are
only 20 deg apart, but their arithmetic mean is 90 deg). To get a correct mean the
angles are doubled, averaged as unit vectors and the result is halved again.
R (0..1) is the length of the mean vector: 1 = all parallel, 0 = no preferred direction.
"""

from __future__ import annotations

import numpy as np

# Simple quality weights = 1 / angular uncertainty of the WSM quality class (deg).
QUALITY_WEIGHTS = {"A": 1 / 15, "B": 1 / 20, "C": 1 / 25}


def axial_mean(azimuths, weights=None) -> dict:
    """Return axial mean orientation, mean resultant length R and circular std.

    Parameters
    ----------
    azimuths : array-like of SHmax azimuths in degrees (0-180, NaN ignored)
    weights  : optional array-like of non-negative weights (same length)

    Returns
    -------
    dict with keys n, mean_deg (0-180), R (0-1), circ_std_deg
    """
    a = np.asarray(azimuths, dtype=float)
    w = np.ones_like(a) if weights is None else np.asarray(weights, dtype=float)
    ok = np.isfinite(a) & np.isfinite(w)
    a, w = a[ok], w[ok]
    if a.size == 0:
        return {"n": 0, "mean_deg": np.nan, "R": np.nan, "circ_std_deg": np.nan}

    doubled = np.deg2rad(2.0 * a)
    C = np.sum(w * np.cos(doubled)) / np.sum(w)
    S = np.sum(w * np.sin(doubled)) / np.sum(w)
    R = float(np.hypot(C, S))
    mean = (np.rad2deg(np.arctan2(S, C)) / 2.0) % 180.0
    # circular standard deviation of the doubled angles, halved back to axial degrees
    circ_std = np.rad2deg(np.sqrt(-2.0 * np.log(R))) / 2.0 if R > 0 else np.inf
    return {"n": int(a.size), "mean_deg": float(mean), "R": R, "circ_std_deg": float(circ_std)}


def rayleigh_test_axial(azimuths) -> dict:
    """Rayleigh test for a preferred orientation, applied to doubled angles.

    H0: orientations are uniformly distributed (no preferred SHmax orientation).
    Uses the approximation of Zar (2010) for the p-value.
    """
    a = np.asarray(azimuths, dtype=float)
    a = a[np.isfinite(a)]
    n = a.size
    if n < 2:
        return {"n": n, "Z": np.nan, "p_value": np.nan}
    R = axial_mean(a)["R"]
    Z = n * R ** 2
    p = np.exp(-Z) * (1 + (2 * Z - Z ** 2) / (4 * n)
                      - (24 * Z - 132 * Z ** 2 + 76 * Z ** 3 - 9 * Z ** 4) / (288 * n ** 2))
    return {"n": n, "Z": float(Z), "p_value": float(max(min(p, 1.0), 0.0))}


def bootstrap_axial_mean_ci(azimuths, n_boot: int = 2000, seed: int = 42) -> tuple:
    """95 % bootstrap confidence interval of the axial mean (degrees, 0-180).

    The interval is computed on deviations from the full-sample mean to avoid the
    0/180 wrap-around problem.
    """
    a = np.asarray(azimuths, dtype=float)
    a = a[np.isfinite(a)]
    if a.size < 5:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    m0 = axial_mean(a)["mean_deg"]
    devs = []
    for _ in range(n_boot):
        m = axial_mean(rng.choice(a, a.size, replace=True))["mean_deg"]
        devs.append(((m - m0 + 90) % 180) - 90)          # signed difference in (-90, 90]
    lo, hi = np.percentile(devs, [2.5, 97.5])
    return ((m0 + lo) % 180, (m0 + hi) % 180)


def orientation_summary(df, group_col: str | None = None, weighted: bool = False):
    """Axial statistics for the whole DataFrame or per group (records with azimuth only).

    Returns a pandas DataFrame with n, mean_deg, R, circ_std_deg, rayleigh p-value and
    a 95 % bootstrap CI of the mean.
    """
    import pandas as pd

    d = df[df["AZI"].notna()]
    groups = [("All", d)] if group_col is None else list(d.groupby(group_col, observed=True))
    rows = []
    for name, sub in groups:
        if len(sub) == 0:
            continue
        w = sub["QUALITY"].astype(str).map(QUALITY_WEIGHTS).values if weighted else None
        st = axial_mean(sub["AZI"].values, w)
        # bootstrap CI only for the unweighted mean (keeps the method simple)
        ci = (np.nan, np.nan) if weighted else bootstrap_axial_mean_ci(sub["AZI"].values)
        ray = rayleigh_test_axial(sub["AZI"].values)
        rows.append({"group": name, **st, "ci95_low": ci[0], "ci95_high": ci[1],
                     "rayleigh_p": ray["p_value"]})
    return pd.DataFrame(rows).set_index("group").round(3)


def gridded_axial_means(df, cell_deg: float = 1.0, min_n: int = 5):
    """Axial mean SHmax per regular lon/lat grid cell (simple binning).

    Only cells with at least ``min_n`` records are returned. This is a deliberately
    simple illustration of spatial variability; it is NOT the WSM smoothing algorithm
    (e.g. Ziegler & Heidbach) used for official mean-SHmax grids.
    """
    import pandas as pd

    d = df[df["AZI"].notna()].copy()
    d["cell_lon"] = np.floor(d["LON"] / cell_deg) * cell_deg + cell_deg / 2
    d["cell_lat"] = np.floor(d["LAT"] / cell_deg) * cell_deg + cell_deg / 2
    rows = []
    for (lon, lat), sub in d.groupby(["cell_lon", "cell_lat"]):
        if len(sub) < min_n:
            continue
        st = axial_mean(sub["AZI"].values)
        rows.append({"cell_lon": lon, "cell_lat": lat, **st})
    return pd.DataFrame(rows)


def azimuth_to_compass(az: float) -> str:
    """Describe an axial azimuth with a two-sided compass label, e.g. 147 -> 'NW-SE'."""
    if not np.isfinite(az):
        return "n/a"
    labels = ["N-S", "NNE-SSW", "NE-SW", "ENE-WSW", "E-W", "WNW-ESE", "NW-SE", "NNW-SSE"]
    return labels[int(((az % 180) + 11.25) // 22.5) % 8]
