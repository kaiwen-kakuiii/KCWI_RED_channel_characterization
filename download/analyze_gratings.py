"""
Analyze which stars have the most unique grating combinations,
and total observation counts per grating, for all IFUs combined
and for each IFU individually (Large, Medium, Small).
"""

import os

import pandas as pd

# ---- Load data ----
HERE = os.path.dirname(os.path.abspath(__file__))
INPUT_CSV = os.path.join(HERE, "outputKC", "std_red_matched.csv")

# skipinitialspace is required: `binning` is written as the quoted field "2,2", and
# without it that comma splits the row and shifts every column after it.
df = pd.read_csv(INPUT_CSV, skipinitialspace=True)
df.columns = [c.strip() for c in df.columns]

# Strip whitespace from string columns we care about
for col in ["std_name", "rgratnam", "ifunam"]:
    df[col] = df[col].astype(str).str.strip()

# std_name, not targname.  fetch_std_red.py assigns it by matching each frame's
# pointing against the starlist, so one physical star carries exactly one value and
# no name normalization is needed.  targname is what the observer typed at the
# telescope -- BD+28 4211 alone appears as seven different strings -- and is not
# used here.
STAR_COL = "std_name"
GRATING_COL = "rgratnam"
IFU_COL = "ifunam"


def analyze(subset: pd.DataFrame, label: str):
    print("=" * 70)
    print(f"IFU: {label}  (n_obs = {len(subset)})")
    print("=" * 70)

    # --- Unique grating combinations per star, and how many frames on each ---
    # gratings_used and n_frames_per_grating are parallel lists: element i of one
    # counts the frames on element i of the other.  groupby sorts the gratings, so
    # both come out in the same order.
    counts = subset.groupby([STAR_COL, GRATING_COL]).size()
    per_star = counts.groupby(level=0)

    star_summary = pd.DataFrame({
        "n_unique_gratings": per_star.size(),
        "gratings_used": per_star.agg(lambda s: list(s.index.get_level_values(1))),
        "n_frames_per_grating": per_star.agg(list),
        "n_frames": per_star.sum(),
    }).sort_values(["n_unique_gratings", "n_frames"], ascending=False)

    print("\n-- Stars ranked by number of unique gratings used --")
    print(star_summary.to_string())

    # --- Total number of observations per grating ---
    grating_counts = (
        subset[GRATING_COL]
        .value_counts()
        .rename_axis("grating")
        .reset_index(name="n_observations")
    )

    print("\n-- Total observations per grating --")
    print(grating_counts.to_string(index=False))
    print()

    return star_summary, grating_counts


results = {}

# All IFUs combined
results["ALL"] = analyze(df, "ALL (combined)")

# Each IFU individually
for ifu_value in ["Large", "Medium", "Small"]:
    subset = df[df[IFU_COL] == ifu_value]
    if subset.empty:
        print(f"(No observations found for IFU = {ifu_value})\n")
        continue
    results[ifu_value] = analyze(subset, ifu_value)

# ---- Also dump full IFU value counts, in case labels differ ----
print("=" * 70)
print("Distinct IFU values found in data (for reference):")
print(df[IFU_COL].value_counts().to_string())
