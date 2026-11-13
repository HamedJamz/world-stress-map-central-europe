"""
Small Streamlit app to filter the WSM data and draw a stress map.

    streamlit run app/streamlit_app.py
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.data_processing import (load_wsm_data, clean_wsm_data, to_geodataframe,  # noqa: E402
                                 RAW_DIR, WSM_FILENAME, TYPE_DESCRIPTIONS)
from src.filtering import REGION_PRESETS, filter_stress_data, filter_by_region  # noqa: E402
from src.orientation_stats import axial_mean, azimuth_to_compass  # noqa: E402
from src.plotting import (plot_stress_map, plot_quality_distribution,  # noqa: E402
                          plot_indicator_distribution)

st.set_page_config(page_title="WSM stress-map explorer", layout="wide")


@st.cache_data(show_spinner="Loading and cleaning WSM 2025 data ...")
def get_data():
    """Load + clean the WSM file once and cache it for the session."""
    raw = load_wsm_data(RAW_DIR / WSM_FILENAME)
    clean, _ = clean_wsm_data(raw)
    return to_geodataframe(clean)


gdf = get_data()

st.title("CASMO-inspired WSM stress-map explorer")
st.caption("Independent student project using the public World Stress Map Database "
           "Release 2025 (Heidbach et al., 2025, doi:10.5880/WSM.2025.001, CC BY 4.0). "
           "Not official WSM/GFZ software.")

# ---------------- sidebar controls -----------------------------------------------------
with st.sidebar:
    st.header("Filters")
    region = st.selectbox("Region", list(REGION_PRESETS), index=1)
    qualities = st.multiselect("WSM quality classes", ["A", "B", "C", "D", "E"],
                               default=["A", "B", "C"],
                               help="A–C are considered reliable (±15°, ±20°, ±25°).")
    in_region = filter_by_region(gdf, region)
    types_available = in_region["TYPE"].value_counts().index.tolist()
    indicators = st.multiselect(
        "Stress-indicator types", types_available, default=types_available,
        format_func=lambda t: f"{t} – {TYPE_DESCRIPTIONS.get(t, '')}")

if not qualities or not indicators:
    st.warning("Please select at least one quality class and one indicator type.")
    st.stop()

sel, counts, _ = filter_stress_data(gdf, region=region, quality=qualities,
                                    indicator=indicators, verbose=False)

# ---------------- key numbers ----------------------------------------------------------
stats = axial_mean(sel["AZI"])
c1, c2, c3, c4 = st.columns(4)
c1.metric("Records in region", len(in_region))
c2.metric("Selected records (with SHmax)", len(sel))
c3.metric(f"Axial mean SHmax ({azimuth_to_compass(stats['mean_deg'])})" if len(sel) else "Axial mean SHmax",
          f"{stats['mean_deg']:.0f}°" if len(sel) else "–",
          help="Circular mean of doubled angles, halved (0° = 180°).")
c4.metric("Mean resultant length R", f"{stats['R']:.2f}" if len(sel) else "–")

# ---------------- map + charts ---------------------------------------------------------
if len(sel) == 0:
    st.info("No records match the current filters.")
    st.stop()

fig = plot_stress_map(sel, extent=REGION_PRESETS[region]["map_extent"],
                      title=f"SHmax orientations – {region}",
                      subtitle=f"Quality: {', '.join(qualities)} · n = {len(sel)}",
                      figsize=(11, 8))
st.pyplot(fig, width="stretch")
plt.close(fig)

left, right = st.columns(2)
with left:
    f = plot_quality_distribution(in_region, title=f"Quality classes – all records in {region}")
    st.pyplot(f); plt.close(f)
with right:
    f = plot_indicator_distribution(in_region[in_region["TYPE"].isin(indicators)],
                                    title="Indicator types (selected)", compact=True)
    st.pyplot(f); plt.close(f)

with st.expander("Filter steps and selected records"):
    st.dataframe(counts, hide_index=True)
    cols = ["ID", "LAT", "LON", "AZI", "TYPE", "QUALITY", "REGIME", "DEPTH", "COUNTRY", "REF1"]
    st.dataframe(sel[cols], hide_index=True)
    st.download_button("Download selection as CSV", sel[cols].to_csv(index=False),
                       file_name="wsm_selection.csv", mime="text/csv")
