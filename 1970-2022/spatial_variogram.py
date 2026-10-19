#!/usr/bin/env python3
"""
Spatial variogram of real-estate prices.

Implements the spatial variogram analysis of Becharat, Benzaquen & Bouchaud (2025)
on the municipality-level (commune-level) price dataset, 1970-2022
(source: Cage & Piketty).

For each year, the script:
  1. samples random pairs of municipalities,
  2. computes the squared log-price difference and the distance for each pair,
  3. averages the squared differences in distance bins -> empirical variogram V(l),
  4. fits V linearly against ln(l) over a configurable distance range
     (the slope is reported as A^2 / (2*pi*D)),
  5. plots V(l) for all years on a log-x axis, colored by year.

The years and fit ranges can be selected in the code below. The current standard is to fit one year every four,
and fit ranges that span (10,200) or (10,250) km.

INPUT FILES: data.csv, data_no_inflation.csv

USAGE (CLI commands):
    python spatial_variogram.py --data data.csv
    python spatial_variogram.py --data data_no_inflation.csv --n-pairs 20_000_000


"""

import argparse
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DEFAULT_DATA_PATH = "data.csv"  # or "data_no_inflation.csv"
DEFAULT_OUTPUT_PATH = "spatial_variogram.png"

YEARS: List[int] = [
    1970, 1974, 1978, 1982, 1986, 1990, 1994,
    1998, 2002, 2006, 2010, 2014, 2018, 2022,
]


# Fitting range (l_min, l_max), in km, of the linear fit of V against ln(l)
# used by `fit_slope`. Edit it year by year; years that are not listed fall
# back to DEFAULT_FIT_RANGE.
DEFAULT_FIT_RANGE: Tuple[float, float] = (10, 200)
FIT_RANGES: Dict[int, Tuple[float, float]] = {
    1970: (10, 200),
    1974: (10, 200),
    1978: (10, 200),
    1982: (10, 200),
    1986: (10, 200),
    1990: (10, 200),
    1994: (10, 200),
    1998: (10, 200),
    2002: (10, 250),
    2006: (10, 250),
    2010: (10, 250),
    2014: (10, 250),
    2018: (10, 250),
    2022: (10, 250),
}

# Map of legacy (Italian) column names to the English names used internally.
COLUMN_ALIASES = {
    "codecommune": "commune_code",
    "anno": "year",
    "log_prezzo": "log_price",
}
REQUIRED_COLUMNS = ["year", "log_price", "x", "y"]

# Minimum number of municipalities with a valid price required to compute
# the variogram of a given year.
MIN_MUNICIPALITIES = 50


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def load_data(path: str) -> pd.DataFrame:
    """Load the municipal price dataset from a CSV file.

    Legacy Italian column names are renamed to their English equivalents.
    Raises a ValueError if a required column is missing.
    """
    data = pd.read_csv(path, dtype={"codecommune": str, "commune_code": str})
    data = data.rename(columns=COLUMN_ALIASES)

    missing = [col for col in REQUIRED_COLUMNS if col not in data.columns]
    if missing:
        raise ValueError(
            f"{path}: missing required column(s) {missing}. "
            f"Found: {list(data.columns)}"
        )
    return data


def compute_variogram(
    data: pd.DataFrame,
    year: int,
    n_pairs: int = 10_000_000,
    bin_edges: Optional[np.ndarray] = None,
    seed: int = 0,
) -> Optional[pd.DataFrame]:
    """Estimate the spatial variogram of log-prices for a given year.

    Random pairs of municipalities are drawn (with replacement), and the mean
    squared log-price difference is computed in distance bins.

    Parameters
    ----------
    data : DataFrame
        Dataset returned by `load_data`.
    year : int
        Year to analyse.
    n_pairs : int
        Number of random pairs to draw (pairs with i == j are discarded).
    bin_edges : array-like, optional
        Distance bin edges in km. Defaults to 100 equal-width (linear) bins
        between 10 and 1000 km.
    seed : int
        Seed of the random number generator, for reproducibility.

    Returns
    -------
    DataFrame with columns ``l_km`` (bin centre), ``V`` (mean squared
    difference, NaN for empty bins) and ``n_pairs`` (pairs per bin),
    or None if the year has fewer than MIN_MUNICIPALITIES valid observations.
    """
    subset = data[data["year"] == year]
    log_price = subset["log_price"].values
    x = subset["x"].values
    y = subset["y"].values
    n = len(subset)
    if n < MIN_MUNICIPALITIES:
        print(f"Year {year}: only {n} municipalities with valid data, skipping.")
        return None

    rng = np.random.default_rng(seed)
    i = rng.integers(0, n, size=n_pairs)
    j = rng.integers(0, n, size=n_pairs)
    mask = i != j
    i, j = i[mask], j[mask]

    dist_km = np.sqrt((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2) / 1000.0
    squared_diff = (log_price[i] - log_price[j]) ** 2

    if bin_edges is None:
        # Equal-width (linear) bins from 10 to 1000 km, 100 bins
        bin_edges = np.linspace(10, 1000, 101)

    bin_idx = np.digitize(dist_km, bin_edges) - 1
    n_bins = len(bin_edges) - 1
    V = np.full(n_bins, np.nan)
    counts = np.zeros(n_bins, dtype=int)
    for b in range(n_bins):
        in_bin = bin_idx == b
        counts[b] = in_bin.sum()
        if counts[b] > 0:
            V[b] = squared_diff[in_bin].mean()

    # Arithmetic centre of each bin (bins are linear)
    centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    return pd.DataFrame({"l_km": centers, "V": V, "n_pairs": counts})


def fit_slope(
    variogram: pd.DataFrame,
    l_min: float = 10,
    l_max: float = 200,
    min_pairs: int = 30,
) -> Optional[Dict[str, float]]:
    """Linear fit of V against ln(l) over the range [l_min, l_max] (km).

    Only bins with at least `min_pairs` pairs are used. Returns a dict with
    ``slope``, ``intercept`` and ``r2``, or None if fewer than 3 bins are
    available.
    """
    v = variogram.dropna(subset=["V"])
    v = v[
        (v["l_km"] >= l_min)
        & (v["l_km"] <= l_max)
        & (v["n_pairs"] >= min_pairs)
    ]
    if len(v) < 3:
        return None
    result = linregress(np.log(v["l_km"]), v["V"])
    return {
        "slope": result.slope,
        "intercept": result.intercept,
        "r2": result.rvalue ** 2,
    }


# ---------------------------------------------------------------------------
# Plotting and command-line interface
# ---------------------------------------------------------------------------

def plot_variograms(
    data: pd.DataFrame,
    years: Sequence[int],
    n_pairs: int,
    seed: int,
    show_lines: bool,
    output_path: str,
    show: bool,
) -> None:
    """Compute and plot the variogram of every year, and print the fitted slopes."""
    marker_style = "o-" if show_lines else "o"

    fig, ax = plt.subplots(figsize=(8, 5))

    # Colormap and normalization used to color each year
    cmap = plt.cm.viridis
    norm = plt.Normalize(vmin=min(years), vmax=max(years))

    for year in years:
        variogram = compute_variogram(data, year, n_pairs=n_pairs, seed=seed)
        if variogram is None:
            continue
        valid = variogram.dropna(subset=["V"])

        ax.plot(valid["l_km"], valid["V"], marker_style, ms=3, color=cmap(norm(year)))

        l_min, l_max = FIT_RANGES.get(year, DEFAULT_FIT_RANGE)
        fit = fit_slope(variogram, l_min=l_min, l_max=l_max)
        if fit:
            print(
                f"{year}: slope A^2/(2*pi*D) = {fit['slope']:.3f}  "
                f"(R^2={fit['r2']:.2f})"
            )

    ax.set_xscale("log")
    ax.set_xlabel(r"$\ell$ (km)")
    ax.set_ylabel(r"$\langle (p(r,t)-p(r',t))^2 \rangle_{|r-r'|=\ell}$")

    # Continuous colorbar for the year
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("Years")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"Figure saved to {output_path}")
    if show:
        plt.show()


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Spatial variogram of real-estate prices (1970-2022)."
    )
    parser.add_argument(
        "--data",
        default=DEFAULT_DATA_PATH,
        help="input CSV file, e.g. data.csv or data_no_inflation.csv "
        f"(default: {DEFAULT_DATA_PATH})",
    )
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT_PATH,
        help=f"output figure path (default: {DEFAULT_OUTPUT_PATH})",
    )
    parser.add_argument(
        "--n-pairs",
        type=int,
        default=10_000_000,
        help="number of random pairs sampled per year (default: 10000000)",
    )
    parser.add_argument(
        "--seed", type=int, default=0, help="random seed (default: 0)"
    )
    parser.add_argument(
        "--lines",
        action="store_true",
        help="connect the points of each year with a line (default: markers only)",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="save the figure without opening a plot window",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    data = load_data(args.data)
    plot_variograms(
        data,
        YEARS,
        n_pairs=args.n_pairs,
        seed=args.seed,
        show_lines=args.lines,
        output_path=args.output,
        show=not args.no_show,
    )


if __name__ == "__main__":
    main()