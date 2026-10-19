"""
================

Spatial and temporal analysis of French housing prices from the DVF open data
("Demandes de Valeurs Foncieres"), reproducing the analyses of

    A.-C. Becharat, M. Benzaquen and J.-P. Bouchaud,
    "The Diffusive Nature of Housing Prices",
    arXiv:2412.14624; Phys. Rev. Lett. 135, 107401 (2025).

##########################################################################################################
    
INPUT.
    
The input is the multi-year, geolocated DVF data (one file per year).
The code accepts multiple DVF files in the working directory.
In particular it accepts the files named after the pattern "DVF-XXXX.txt", where XXXX is the year of the file.
You can change the pattern accordingly with how you name your files. The standard we used in the CLI command is:
python dvf.py "/path/to/DVF-*.txt"



##########################################################################################################

INSTRUCTIONS.

Below are listed the multiple possibilites of the code. The key python functions related to every task are written
next to them (they're also commented in the code).
The last line of every section gives the CLI command to paste into the terminal to execute that specific task.

Beware: The code will take several minutes to run.


The module covers SIX analyses, all usable from the CLI command given below:

  A) Spatial variogram at ONE scale         -> run_full_pipeline(...)

    python dvf.py "/path/to/DVF-*.txt" --analysis spatial --scale france


  
  B) Price distributions (Fig. 2)           -> plot_price_distribution_by_year(...)
                                               plot_price_distribution_by_month(...)
                                               plot_monthly_volume_and_count(...)

    python dvf.py "/path/to/DVF-*.txt" --analysis distribution
                                           
    
    
  C) MULTI-scale spatial variogram          -> compute_multiscale_variogram(...) +
     (full replica of Fig. 3)                  plot_multiscale_variogram(...)

     
    python dvf.py "/path/to/DVF-*.txt" --analysis multiscale




  D) TEMPORAL variogram V1/V2 + fit         -> build_spatiotemporal_field(...) +
     (SPDE model: diffusion + noise with       compute_V1/compute_V2(...) +
      memory T, spatial cutoff S;              fit_temporal_variogram(...) +
      B=A^2/D free or fixed, nugget)           plot_temporal_variogram[_interactive](...)

    python dvf.py "/path/to/DVF-*.txt" --analysis temporal --nugget --fit-target V2


      
  E) Prices PER SQUARE METER (P^ = price/m^2) -> clean_dvf(..., single_line_only=True) +
     - price distribution (Fig. 7)              add_price_per_m2(...)
     - Paris vs Creuse (like Fig. 6, per m^2)   plot_price_m2_distribution_by_year(...)
     - multi-scale variogram (Fig. 5)           plot_paris_creuse_m2_distribution(...)
                                                compute_multiscale_variogram(..., log_price_col="log_price_m2") +
                                                plot_price_m2_multiscale_variogram(...)

                                                
    python dvf.py "/path/to/DVF-*.txt" --analysis m2
                                                
    

  F) Distribution of the log-prices PARIS vs CREUSE -> plot_paris_creuse_distribution(...)
     (replica of Fig. 6: city vs countryside);

      
    python dvf.py "/path/to/DVF-*.txt" --analysis paris_creuse

#############################################################################################################

    
METHODOLOGY (following the paper, main text and SM-2.2)
--------------------------------------------------------
1. The log-price field  p(r) := log(P)  is averaged inside hexagonal cells of fixed area:
        - 0.73 km^2 -> "city" scale        (H3 resolution 8)
        - 5    km^2 -> "department" scale  (H3 resolution 7)
        - 30   km^2 -> "region" scale      (H3 resolution 6, actual area 36.1 km^2)
        - 250  km^2 -> "France" scale      (H3 resolution 5, actual area 252.9 km^2)
   These areas coincide (almost exactly for resolutions 5, 7 and 8) with the average
   cell areas of Uber's hierarchical H3 hexagonal grid, so H3 is used here for the
   binning (this is an assumption of this implementation).

2. The field p(r) is averaged over ALL the years/files available: the transactions of
   all the DVF files are concatenated, each one is assigned to its hexagonal cell and
   the log-price is averaged per cell over the whole time window.

3. The empirical variogram

        V(l) = < (p(r) - p(r'))^2 >_{|r-r'| = l}

   is computed by averaging over pairs of cells at distance l (logarithmic distance
   binning by default, given the huge number of possible pairs for the whole of France;
   the multi-scale analysis uses linear bins, see below).

4. It is compared with the theoretical prediction of the paper (Eq. 3):

        V(l) ~ (A^2 / 2*pi*D) * log(l / l*)      for  max(a, sqrt(D*T)) << l << l*

   i.e. a LOGARITHMIC dependence on l, with slope A^2/(2*pi*D).

5. The temporal variogram uses the same fit model used in the script "temporal_variogram.py".

##########################################################################################################


DATA SELECTION: NOMINAL PRICES vs PRICES PER SQUARE METER
---------------------------------------------------------
A DVF "transaction" is the set of ALL the lines (rows) of the raw data that share the same
id_mutation. valeur_fonciere is the TOTAL price of the sale and is repeated on every line of
the transaction. The two kinds of analysis use different selections:

  * NOMINAL prices (analyses A, B, C, D, F)      ->  clean_dvf(raw)
      Keeps the transactions with exactly ONE dwelling line (Maison/Appartement) and no line
      with code_type_local = 4, and keeps only that dwelling line. Other lines of the same
      transaction (outbuildings, land) are allowed: the total price is attributed to the
      dwelling.

  * Prices per SQUARE METER (analysis E)         ->  add_price_per_m2(clean_dvf(raw, single_line_only=True))
      Keeps ONLY single-line transactions: transactions made of exactly ONE DVF line, which
      must be a Maison or an Appartement. A transaction with further lines (land,
      outbuildings, ...) has a total price that also pays for things that are not part of
      the built surface used as denominator (surface_reelle_bati), which would bias the
      price per m^2.

  The two selections are applied independently, always on the raw data. Single-line
  transactions are counted BEFORE any row-level filter (coordinates, price, ...), so a
  transaction cannot look "single-line" just because one of its other lines was dropped.

 
  
############################################################################################################    

USAGE (command line)
----------------------------
    python dvf.py "/path/to/DVF-*.txt" --analysis spatial --scale france
    python dvf.py "/path/to/DVF-*.txt" --analysis distribution
    python dvf.py "/path/to/DVF-*.txt" --analysis distribution --tail-equal-weight-bins --tail-n-bins 15
    python dvf.py "/path/to/DVF-*.txt" --analysis multiscale
    python dvf.py "/path/to/DVF-*.txt" --analysis temporal
    python dvf.py "/path/to/DVF-*.txt" --analysis temporal --a2-over-d 1.2
    python dvf.py "/path/to/DVF-*.txt" --analysis temporal --nugget --fit-target V2
    python dvf.py "/path/to/DVF-*.txt" --analysis temporal --interactive
    python dvf.py "/path/to/DVF-*.txt" --analysis distribution_m2   # Fig. 7 (price per m^2)
    python dvf.py "/path/to/DVF-*.txt" --analysis multiscale_m2     # Fig. 5 (price per m^2)
    python dvf.py "/path/to/DVF-*.txt" --analysis m2                # all the per-m^2 analyses, loading the data once
    python dvf.py "/path/to/DVF-*.txt" --analysis paris_creuse      # Fig. 6 (Paris vs Creuse)
    python dvf.py "/path/to/DVF-*.txt" --analysis paris_creuse_m2   # Paris vs Creuse, price per m^2
    python dvf.py "/path/to/DVF-*.txt" --analysis distribution --min-count-per-bin 0   # fixed bins

#################################################################################################################    






"""

from __future__ import annotations

import glob
import warnings
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
import pandas as pd

try:
    import h3
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "The 'h3' package is not installed. Run: pip install h3"
    ) from exc

try:
    from pyproj import Transformer
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "The 'pyproj' package is not installed. Run: pip install pyproj"
    ) from exc


# ---------------------------------------------------------------------------
# 1. Configuration constants
# ---------------------------------------------------------------------------

# Paper scale <-> H3 resolution. Average cell area per H3 resolution (km^2), for reference:
#   res 5: 252.90 km^2   res 6: 36.13 km^2   res 7: 5.16 km^2   res 8: 0.737 km^2
SCALE_TO_H3_RES: dict[str, int] = {
    "city": 8,          # paper: 0.73 km^2
    "departement": 7,   # paper: 5 km^2
    "region": 6,        # paper: 30 km^2 (actual H3 area: 36.1 km^2)
    "france": 5,        # paper: 250 km^2
}

# ---------------------------------------------------------------------------
# >>> MINIMUM NUMBER OF TRANSACTIONS PER HEXAGONAL CELL <<<
# Single place to change it: it is the default of every function of the module
# (build_price_field, build_spatiotemporal_field, run_full_pipeline, ...) and of the
# --min-transactions command-line option. Cells with fewer transactions are discarded (or
# set to NaN) because their mean log-price is considered too noisy. Every function also
# accepts an explicit `min_transactions` argument to override it case by case.
MIN_TRANSACTIONS_PER_CELL: int = 3

# ---------------------------------------------------------------------------
# >>> CLEANING / FIT PARAMETERS FOR THE PRICE PER SQUARE METER ANALYSES <<<
#     (P^ = valeur_fonciere / surface_reelle_bati), used as defaults by
#     add_price_per_m2 and plot_price_m2_distribution_by_year.
#   MIN_SURFACE_M2          : minimum surface (m^2) to keep a transaction (drops null or
#                             tiny surfaces, which give absurdly large prices per m^2)
#   PRICE_M2_BOUNDS         : (min, max) price per m^2 in EUR/m^2 (drops typos and
#                             extreme outliers)
#   PRICE_M2_TAIL_THRESHOLD : threshold (EUR/m^2) above which the power-law tail
#                             P^(-1-mu) is fitted (panel of the distribution plot)
MIN_SURFACE_M2: float = 9.0
PRICE_M2_BOUNDS: tuple[float, float] = (100.0, 200_000.0)
PRICE_M2_TAIL_THRESHOLD: float = 10_000.0

# ---------------------------------------------------------------------------
# >>> TAIL-FIT OPTIONS (they apply to BOTH distributions: total prices, Fig. 2, and
#     prices per m^2, Fig. 7). <<<
#   TAIL_FIT_EQUAL_WEIGHT_BINS : False -> the regression uses ALL the points of the tail:
#                                the points are much denser close to the threshold than
#                                in the extreme tail, so the last points weigh little.
#                                True (default) -> the tail is split into log-spaced bins
#                                and every bin weighs EXACTLY the same in the fit.
#   TAIL_FIT_N_BINS            : number of log-spaced bins (only if the above is True)
# Every plot function also accepts the arguments `tail_equal_weight_bins` and
# `tail_n_bins`, and the command line has --tail-equal-weight-bins /
# --no-tail-equal-weight-bins and --tail-n-bins.
TAIL_FIT_EQUAL_WEIGHT_BINS: bool = True
TAIL_FIT_N_BINS: int = 30

# ---------------------------------------------------------------------------
# >>> ADAPTIVE BINS IN THE DISTRIBUTION PLOTS (they apply to ALL of them: prices by year
#     and by month, prices per m^2, Paris vs Creuse, Paris vs Creuse per m^2). <<<
# The problem: with fixed-width bins, in the tails every bin holds 0, 1, 2, 3...
# transactions. The count is a small integer, so the points line up on HORIZONTAL STRIPES
# (y = 1, 2, 3, ... for counts; y = k/(N*width) for densities) and, in x, on a regular grid
# (the bin centres) instead of where the transactions are.
# The solution (see _adaptive_log_histogram):
#   - where the base bins hold at least DIST_MIN_COUNT_PER_BIN transactions (the "body" of
#     the distribution) nothing changes: same bins as plain fixed bins;
#   - in the tails, bins that are too poor are MERGED, from the body outwards, until every
#     point represents at least DIST_MIN_COUNT_PER_BIN transactions. The width of a merged
#     bin is its real width, so the density (count/width) is a continuous number and no
#     longer an integer;
#   - every point is drawn at the MEAN (geometric) position of the actual prices of the
#     transactions it represents, not at the bin centre.
#   DIST_MIN_COUNT_PER_BIN : integer >= 1 (higher = smoother tails but fewer points;
#                            lower = more points but noisier), or None for plain fixed
#                            bins (points at the bin centre, horizontal stripes included).
# Every plot function also accepts the `min_count_per_bin` argument, and the command line
# has --min-count-per-bin (0 = None).
DIST_MIN_COUNT_PER_BIN: Optional[int] = 1

# ---------------------------------------------------------------------------
# >>> TRANSACTION SELECTION RULES (used by clean_dvf) <<<
# A "transaction" is the set of ALL the DVF lines with the same id_mutation (a sale can
# span several lines: the dwelling, its outbuildings, land...).
#
# NOMINAL prices - select_residential_transactions (default of clean_dvf), in this order:
#   1. keep only the transactions with at least one Maison or Appartement line;
#   2. discard the WHOLE transaction if any of its lines has code_type_local = 4
#      (industrial/commercial premises);
#   3. discard the WHOLE transaction unless it contains EXACTLY ONE Maison/Appartement
#      line (so also 1 Maison + 1 Appartement, 2 Maison, ... are discarded);
#   4. of the remaining transactions keep ONLY the Maison/Appartement line: every final
#      row is thus a different transaction (the other lines of the same transaction, i.e.
#      outbuildings and land, are dropped).
#
# PRICES PER m^2 - select_single_line_transactions (clean_dvf(..., single_line_only=True)):
#   keep ONLY the transactions made of exactly ONE DVF line whose type_local is a
#   Maison or an Appartement. Land and outbuildings are not allowed even as extra lines,
#   because the total price would also pay for them while surface_reelle_bati only
#   measures the dwelling.
#
#   DWELLING_TYPES            : type_local values that count as a "dwelling"
#   EXCLUDED_CODES_TYPE_LOCAL : code_type_local values that make the transaction be
#                               discarded in the nominal-price selection
#                               (4 = "Local industriel. commercial ou assimile")
#   N_LINES_COLUMN            : column added by select_single_line_transactions with the
#                               number of DVF lines of the transaction (always 1 there);
#                               add_price_per_m2 uses it to check that it receives
#                               single-line data.
DWELLING_TYPES: tuple[str, ...] = ("Maison", "Appartement")
EXCLUDED_CODES_TYPE_LOCAL: tuple[int, ...] = (4,)
N_LINES_COLUMN: str = "n_lines_in_mutation"

# ---------------------------------------------------------------------------
# >>> PARIS vs CREUSE PLOT (replica of Fig. 6 of the paper) <<<
#   PARIS_CREUSE_DEPARTEMENTS : {legend label: code_departement}
#   PARIS_CREUSE_PRICE_RANGE  : price range P (EUR) of the x axis. It is also the price
#                               limit to clean the data with for this plot (Fig. 6 shows
#                               all the transactions up to 1e8 EUR; today this equals the
#                               default price_bounds of clean_dvf)
#   PARIS_CREUSE_N_BINS       : number of bins (log-spaced in P) inside the range
#   PARIS_CREUSE_Y_LIM        : limits of the y axis (density of the log-price)
PARIS_CREUSE_DEPARTEMENTS: dict[str, str] = {"Paris": "75", "Creuse": "23"}
PARIS_CREUSE_PRICE_RANGE: tuple[float, float] = (1_000.0, 100_000_000.0)
PARIS_CREUSE_N_BINS: int = 200
PARIS_CREUSE_Y_LIM: tuple[float, float] = (5e-4, 1.0)

# ---------------------------------------------------------------------------
# >>> PARIS vs CREUSE PLOT on the PRICES PER SQUARE METER (P^ = price/m^2) <<<
# Same plot as above, but with p^ = log(P^) instead of p = log(P) (see
# plot_paris_creuse_m2_distribution). The transactions must be prepared with
# add_price_per_m2 (on single-line data).
#   PARIS_CREUSE_M2_PRICE_RANGE : range of P^ (EUR/m^2) of the x axis and of the histogram
#                                 (default = PRICE_M2_BOUNDS, i.e. everything that
#                                 add_price_per_m2 keeps)
#   PARIS_CREUSE_M2_N_BINS      : number of base bins (log-spaced in P^); 130 on
#                                 [1e2, 2e5] gives the same bin width (~5.8% in P) as
#                                 PARIS_CREUSE_N_BINS = 200 on [1e3, 1e8]
#   PARIS_CREUSE_M2_Y_LIM       : limits of the y axis (density of p^); None = automatic.
#                                 Paris' peak exceeds 1 because the per-m^2 distribution
#                                 is much narrower than the total-price one; the lower
#                                 limit is lower than PARIS_CREUSE_Y_LIM because with
#                                 adaptive bins the tails reach lower densities
PARIS_CREUSE_M2_PRICE_RANGE: tuple[float, float] = PRICE_M2_BOUNDS
PARIS_CREUSE_M2_N_BINS: int = 130
PARIS_CREUSE_M2_Y_LIM: Optional[tuple[float, float]] = (1e-4, 3.0)

# Coordinate transformer: WGS84 (lon/lat, as in the DVF files) -> Lambert-93
# (EPSG:2154), the standard metric projected system for metropolitan France. Used to
# compute distances in km with minimal distortion.
_TO_LAMBERT93 = Transformer.from_crs("EPSG:4326", "EPSG:2154", always_xy=True)


# ---------------------------------------------------------------------------
# 2. Loading the DVF files
# ---------------------------------------------------------------------------

def load_dvf_files(
    file_pattern: str,
    sep: Optional[str] = None,
    usecols: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """
    Load and concatenate all the files matching `file_pattern`
    (e.g. "/data/DVF-*.txt"), adding a 'source_file' column and a year column 'annee'.

    The year is taken from the file name when it follows the DVF-<YYYY> convention
    (e.g. DVF-2020.txt), otherwise from date_mutation.

    sep=None -> the separator is guessed automatically for every file
    (the official DVF files use ',' or ';' depending on the source).
    """
    paths = sorted(glob.glob(file_pattern))
    if not paths:
        raise FileNotFoundError(f"No file found for the pattern: {file_pattern}")

    default_cols = usecols or [
        "id_mutation",
        "date_mutation",
        "numero_disposition",
        "nature_mutation",
        "valeur_fonciere",
        "code_commune",
        "nom_commune",
        "code_departement",
        "id_parcelle",
        "lot1_numero",
        "nombre_lots",
        "code_type_local",
        "type_local",
        "surface_reelle_bati",
        "nombre_pieces_principales",
        "surface_terrain",
        "longitude",
        "latitude",
    ]

    read_kwargs = dict(sep=sep, dtype={"code_commune": str, "code_departement": str})
    if sep is None:
        # sep=None requires the 'python' engine (separator sniffing),
        # which does not support low_memory
        read_kwargs["engine"] = "python"
    else:
        read_kwargs["low_memory"] = False

    frames = []
    for path in paths:
        try:
            df = pd.read_csv(path, usecols=lambda c: c in default_cols, **read_kwargs)
        except ValueError:
            # fallback: some optional columns may not exist in every year of the DVF
            # format -> load everything and filter afterwards
            df = pd.read_csv(path, **read_kwargs)
            df = df[[c for c in default_cols if c in df.columns]]

        df["source_file"] = path
        # year deduced from the DVF-XXXX file name or, failing that, from date_mutation
        year_from_name = "".join(ch for ch in path.split("DVF-")[-1][:4] if ch.isdigit())
        if len(year_from_name) == 4:
            df["annee"] = int(year_from_name)
        elif "date_mutation" in df.columns:
            df["annee"] = pd.to_datetime(df["date_mutation"], errors="coerce").dt.year
        frames.append(df)

    out = pd.concat(frames, ignore_index=True)
    print(f"[load_dvf_files] Loaded {len(paths)} files, {len(out):,} rows in total.")
    return out


# ---------------------------------------------------------------------------
# 3. Data cleaning
# ---------------------------------------------------------------------------

def select_residential_transactions(
    df: pd.DataFrame,
    dwelling_types: Sequence[str] = DWELLING_TYPES,
    excluded_codes_type_local: Sequence[int] = EXCLUDED_CODES_TYPE_LOCAL,
    id_col: str = "id_mutation",
) -> pd.DataFrame:
    """
    Select the TRANSACTIONS (= groups of lines with the same `id_col`) that represent the
    sale of ONE SINGLE dwelling and return, for each of them, only the dwelling line.

    This is the selection used for the NOMINAL (total) prices. For the prices per square
    meter use select_single_line_transactions instead, which is stricter.

    Rules, in order (see also the block of constants at the top of the file):
      1. keep only the transactions with at least one line with type_local in
         `dwelling_types` (Maison / Appartement);
      2. if AT LEAST ONE line of the transaction has code_type_local in
         `excluded_codes_type_local` (default 4 = industrial/commercial premises) the
         WHOLE transaction is discarded;
      3. the WHOLE transaction is discarded unless it contains EXACTLY ONE Maison/Appartement
         line: transactions with 1 Maison + 1 Appartement, 2 Maison, 2 Appartement, ...
         are therefore dropped too;
      4. of the remaining transactions keep ONLY the Maison/Appartement line (one per
         transaction): the other lines (outbuildings, land without type_local) are only
         used to decide rules 2-3 and are then discarded.

    These rules must be applied BEFORE any row-level filter (coordinates, price, ...):
    if a line were removed beforehand (e.g. premises with code_type_local 4 and no
    coordinates) the rule would not see it and the transaction would wrongly pass.
    Lines without `id_col` belong to no transaction and are discarded.
    """
    missing = [c for c in (id_col, "type_local", "code_type_local") if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing columns for the selection of the transactions: {missing}. "
            "Check that they are among the columns loaded by load_dvf_files "
            "(or use clean_dvf(..., select_transactions=False))."
        )

    # per-line flags: is it a dwelling? is it premises to exclude (code_type_local 4)?
    is_dwelling = df["type_local"].isin(list(dwelling_types)).to_numpy()
    code = pd.to_numeric(df["code_type_local"], errors="coerce")  # robust to float/str/NaN
    is_excluded = code.isin(list(excluded_codes_type_local)).to_numpy()

    # One transaction = one id_mutation. factorize gives each id an integer code
    # (-1 if the id is missing): the per-transaction counts are done with np.bincount, much
    # faster and lighter than a groupby on strings when there are millions of lines.
    codes, uniques = pd.factorize(df[id_col])
    n_tx = len(uniques)
    has_id = codes >= 0
    n_dwelling = np.bincount(codes[has_id & is_dwelling], minlength=n_tx)  # Maison/Appartement lines
    n_excluded = np.bincount(codes[has_id & is_excluded], minlength=n_tx)  # lines with code_type_local 4

    # rules 1-3, evaluated per transaction
    has_dwelling = n_dwelling > 0                                  # rule 1
    drop_code4 = has_dwelling & (n_excluded > 0)                   # rule 2
    drop_multi = has_dwelling & ~drop_code4 & (n_dwelling > 1)     # rule 3
    keep_tx = has_dwelling & ~drop_code4 & ~drop_multi             # = exactly 1 dwelling line, no code 4

    # rule 4: of every kept transaction only the dwelling line is retained
    keep_row = np.zeros(len(df), dtype=bool)
    keep_row[has_id] = keep_tx[codes[has_id]] & is_dwelling[has_id]
    out = df[keep_row]

    types_txt = "/".join(dwelling_types)
    print(
        f"[select_residential_transactions] {int(has_dwelling.sum()):,} transactions with {types_txt}: "
        f"discarded {int(drop_code4.sum()):,} (with a line having code_type_local in "
        f"{tuple(excluded_codes_type_local)}) + {int(drop_multi.sum()):,} (not exactly "
        f"one {types_txt} line) -> {int(keep_tx.sum()):,} kept "
        f"({len(out):,} rows, one per transaction)."
    )
    return out


def select_single_line_transactions(
    df: pd.DataFrame,
    dwelling_types: Sequence[str] = DWELLING_TYPES,
    id_col: str = "id_mutation",
) -> pd.DataFrame:
    """
    Select the SINGLE-LINE transactions: those made of exactly ONE DVF line (a single line
    with that `id_col` in `df`) whose type_local is in `dwelling_types` (Maison /
    Appartement). Used for the analyses of the prices per square meter.

    Why: valeur_fonciere is the total price of the whole transaction and is repeated on
    all its lines. If the sale has further lines (land, outbuildings such as a cellar or a
    garage, a second dwelling...) the price pays for all of them, while
    surface_reelle_bati only measures the dwelling, so price/surface would be biased
    (typically overestimated, e.g. for a house sold with a large plot). With a single
    line the whole price belongs to the one dwelling whose surface is the denominator.

    This is stricter than select_residential_transactions (used for the nominal prices),
    which also keeps a dwelling sold together with outbuildings or land.

    The lines of each transaction are counted on the `df` passed in, so it must be the
    complete data (after the nature_mutation filter at most), BEFORE any row-level filter
    (coordinates, price, ...): otherwise a transaction whose other lines had been dropped
    would look single-line. Lines without `id_col` belong to no transaction and are
    discarded.

    Returns a copy of the selected rows with an extra column N_LINES_COLUMN holding the
    number of DVF lines of their transaction (always 1); add_price_per_m2 checks it.
    """
    missing = [c for c in (id_col, "type_local") if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing columns for the selection of the single-line transactions: {missing}. "
            "Check that they are among the columns loaded by load_dvf_files."
        )

    is_dwelling = df["type_local"].isin(list(dwelling_types)).to_numpy()

    # One transaction = one id_mutation. factorize gives each id an integer code
    # (-1 if the id is missing); lines per transaction are counted with np.bincount.
    codes, uniques = pd.factorize(df[id_col])
    n_tx = len(uniques)
    has_id = codes >= 0
    n_lines = np.bincount(codes[has_id], minlength=n_tx)                   # DVF lines per transaction
    n_dwelling = np.bincount(codes[has_id & is_dwelling], minlength=n_tx)  # Maison/Appartement lines

    has_dwelling = n_dwelling > 0
    # a transaction with exactly one line that is a dwelling (the only line IS the dwelling)
    keep_tx = has_dwelling & (n_lines == 1)

    keep_row = np.zeros(len(df), dtype=bool)
    keep_row[has_id] = keep_tx[codes[has_id]]
    out = df[keep_row].copy()
    out[N_LINES_COLUMN] = n_lines[codes[keep_row]]

    types_txt = "/".join(dwelling_types)
    n_multi = int((has_dwelling & (n_lines > 1)).sum())
    print(
        f"[select_single_line_transactions] {int(has_dwelling.sum()):,} transactions with {types_txt}: "
        f"discarded {n_multi:,} (more than one DVF line) -> {int(keep_tx.sum()):,} kept "
        f"({len(out):,} rows, one line per transaction)."
    )
    return out


def clean_dvf(
    df: pd.DataFrame,
    type_locaux: Sequence[str] = ("Maison", "Appartement"),
    select_transactions: bool = True,
    price_bounds: tuple[float, float] = (1_000, 100_000_000),
    single_line_only: bool = False,
) -> pd.DataFrame:
    """
    Apply the standard filters for a clean price analysis:

    - nature_mutation == 'Vente' (excludes donations, exchanges, etc.)
    - (optional, default True) SELECTION OF THE TRANSACTIONS by id_mutation. In the DVF
      format valeur_fonciere is the TOTAL price of the sale, repeated on every line of the
      transaction, and a sale can include outbuildings (cellar, garage) and land besides
      the dwelling. Two rules are available (see the block of constants at the top of the
      file):

        * single_line_only=False (default) -> NOMINAL prices:
          select_residential_transactions.
            1. keep only the transactions with at least one Maison or Appartement line;
            2. if a line of the transaction has code_type_local = 4 (industrial/commercial
               premises) discard the whole transaction;
            3. discard the transaction unless it contains exactly ONE Maison/Appartement
               line (also 1 Maison + 1 Appartement);
            4. at the end keep only the Maison/Appartement lines, which all belong to
               different transactions (the outbuilding and land lines of the same
               transaction are discarded).
          The price is attributed to that single dwelling.

        * single_line_only=True -> PRICES PER SQUARE METER:
          select_single_line_transactions. Keeps ONLY the transactions made of exactly one
          DVF line, which must be a Maison or an Appartement (a dwelling sold with land or
          outbuildings is discarded, since its price would not match the surface of the
          dwelling alone). The result has an extra column N_LINES_COLUMN. Use it as
          add_price_per_m2(clean_dvf(raw, single_line_only=True)).

      With select_transactions=False this step is skipped and only the rows are filtered.
    - type_local in `type_locaux` (default Maison and Appartement; with the selection on,
      only dwellings remain anyway, so `type_locaux` serves to restrict further, e.g. to
      Maison only)
    - valid coordinates (longitude/latitude not null, inside metropolitan France)
    - valeur_fonciere within `price_bounds` (removes typos and extreme outliers)
    - exclusion of Corsica (2A, 2B)

    Adds the column 'log_price' = log(valeur_fonciere).
    """
    if single_line_only and not select_transactions:
        raise ValueError(
            "single_line_only=True needs select_transactions=True "
            "(the single-line rule is a transaction selection)."
        )

    n0 = len(df)
    out = df.copy()

    if "nature_mutation" in out.columns:
        out = out[out["nature_mutation"] == "Vente"]

    # Selection of the transactions. It must come BEFORE all the row-level filters below
    # (type_local, coordinates, price, Corsica): otherwise they would drop lines that are
    # needed to decide the rules (see the docstrings of the select_* functions).
    if select_transactions:
        if single_line_only:
            out = select_single_line_transactions(out)
        else:
            out = select_residential_transactions(out)

    if "type_local" in out.columns and type_locaux:
        out = out[out["type_local"].isin(type_locaux)]

    out = out.dropna(subset=["longitude", "latitude", "valeur_fonciere"])
    out = out[(out["longitude"].between(-6, 10)) & (out["latitude"].between(41, 52))]

    out = out[out["valeur_fonciere"].between(*price_bounds)]

    # Exclude Corse-du-Sud (2A) and Haute-Corse (2B)
    if "code_departement" in out.columns:
        out = out[~out["code_departement"].isin(["2A", "2B"])]

    out["log_price"] = np.log(out["valeur_fonciere"].astype(float))

    print(
        f"[clean_dvf] {n0:,} -> {len(out):,} rows after cleaning "
        f"({len(out) / n0:.1%} kept)."
    )
    return out


# ---------------------------------------------------------------------------
# 4. Hexagonal binning (H3) and construction of the field p(r)
# ---------------------------------------------------------------------------

def assign_hex_cells(df: pd.DataFrame, resolution: int) -> pd.DataFrame:
    """Add a 'hex_id' column with the H3 cell of every transaction."""
    out = df.copy()
    out["hex_id"] = [
        h3.latlng_to_cell(lat, lon, resolution)
        for lat, lon in zip(out["latitude"], out["longitude"])
    ]
    return out


def build_price_field(
    df: pd.DataFrame,
    resolution: int,
    min_transactions: int = MIN_TRANSACTIONS_PER_CELL,
    log_price_col: str = "log_price",
) -> pd.DataFrame:
    """
    Build the field p(r): one row per hexagonal cell, with the mean log-price computed over
    ALL the available transactions (all the years concatenated).

    min_transactions: minimum number of transactions for the mean of a cell to be
    considered reliable (cells with very few transactions are discarded so that they do
    not dominate the variogram with measurement noise, consistently with the paper's
    discussion of the "measurement noise" that shifts the curves vertically).

    log_price_col: column of `df` to average per cell (default 'log_price' = log of the
    total price; use 'log_price_m2' for the log-price per square meter, see
    add_price_per_m2).
    """
    hexed = assign_hex_cells(df, resolution)

    field = (
        hexed.groupby("hex_id")
        .agg(mean_log_price=(log_price_col, "mean"), n_transactions=(log_price_col, "size"))
        .reset_index()
    )
    field = field[field["n_transactions"] >= min_transactions].reset_index(drop=True)

    centers = [h3.cell_to_latlng(h) for h in field["hex_id"]]
    field["latitude"] = [c[0] for c in centers]
    field["longitude"] = [c[1] for c in centers]

    x, y = _TO_LAMBERT93.transform(field["longitude"].to_numpy(), field["latitude"].to_numpy())
    field["x_km"] = x / 1000.0
    field["y_km"] = y / 1000.0

    print(
        f"[build_price_field] H3 resolution={resolution} "
        f"(average area {h3.average_hexagon_area(resolution, 'km^2'):.2f} km^2): "
        f"{len(field):,} cells with >= {min_transactions} transactions."
    )
    return field


# ---------------------------------------------------------------------------
# 5. Price distributions (replica of Fig. 2 of the paper)
# ---------------------------------------------------------------------------
#
# Fig. 2 of the paper shows, for each year separately, the probability density of the
# LOG-price p := log(P), but drawn on an x axis in P (log-log plot): a "double hump" shape
# (city / countryside mixture) and a power-law right tail P^{-1-mu} with mu ~ 1.5 are
# visible. Each year is NORMALIZED separately (density, area=1) and overlaid on the others
# to show that the shape is stable in time.
#
# The monthly distributions (NOT normalized) use exactly the same binning logic but with
# density=False, so that the transaction VOLUMES of different months can be compared.


def _merge_sparse_tail(
    values: np.ndarray,
    inner: float,
    far: float,
    min_count: int,
    min_width: float,
) -> list[tuple[int, float, float, float]]:
    """
    Group in ADAPTIVE bins the log-prices `values` of ONE tail of a distribution, going from
    the body outwards (increasing values). Returns a list of (n, mean_p, low_edge,
    high_edge), one element per bin, from the innermost to the outermost.

      values    : log-prices of the tail (all >= inner)
      inner     : inner edge, i.e. where the body ends (the first bin starts here)
      far       : outer edge of the last bin (>= max(values))
      min_count : minimum number of transactions per bin
      min_width : minimum width (in p) of a bin: merged bins never become narrower than
                  the base bins, so no noise is introduced where the data are abundant

    A bin is closed as soon as it holds at least `min_count` transactions and is at least
    `min_width` wide. The inner edges between two bins are the midpoints between the last
    price of one bin and the first of the next: the widths are therefore real numbers (not
    multiples of a grid) and so are the densities (count/width), which is exactly what
    avoids the horizontal stripes. Equal prices are never separated (DVF is full of round
    prices). Any remainder at the end (fewer than `min_count` transactions) is merged into
    the last bin; if the whole tail has fewer than `min_count` transactions it becomes a
    single bin.
    """
    if values.size == 0:
        return []
    u, mult = np.unique(values, return_counts=True)  # distinct prices (increasing) and multiplicity
    cum_n = np.concatenate(([0], np.cumsum(mult)))
    q = len(u)

    bins = []  # (first distinct price, last distinct price, low edge, high edge)
    first, lo = 0, inner
    for i in range(q):
        if cum_n[i + 1] - cum_n[first] < min_count:
            continue
        hi = 0.5 * (u[i] + u[i + 1]) if i < q - 1 else far
        if hi - lo >= min_width:
            bins.append((first, i, lo, hi))
            first, lo = i + 1, hi
    if first < q:  # remainder: too few transactions (or bin too narrow) for a bin of its own
        if bins:
            f0, _, lo0, _ = bins.pop()
            bins.append((f0, q - 1, lo0, far))
        else:
            bins.append((0, q - 1, inner, far))

    out = []
    for i0, i1, lo_, hi_ in bins:
        m_ = mult[i0:i1 + 1]
        n_ = int(m_.sum())
        # mean of the prices of the bin computed relative to the first price of the group:
        # no precision is lost through cancellation and, if all the prices are equal, it
        # is exact
        mean = float(u[i0] + np.dot(u[i0:i1 + 1] - u[i0], m_) / n_)
        out.append((n_, mean, float(lo_), float(hi_)))
    return out


def _adaptive_log_histogram(
    log_prices: np.ndarray,
    bin_edges: np.ndarray,
    density: bool,
    min_count: Optional[int] = DIST_MIN_COUNT_PER_BIN,
) -> pd.DataFrame:
    """
    Histogram of p = log(P) with ADAPTIVE BINS, used by all the distribution plots instead
    of a plain fixed-bin histogram with "points at the bin centre".

    Why: with fixed-width bins, in the tails every bin holds 0, 1, 2, 3... transactions.
    The count is a small integer, so the points line up on HORIZONTAL STRIPES (y = 1, 2,
    3, ... for counts; y = k/(N*width) for densities) and on a regular grid of x (the bin
    centres), instead of where the transactions really are.

    How:
      1. start from the base bins `bin_edges` (equally spaced in p);
      2. the "body" is the contiguous series of base bins, around the most populated one,
         with at least `min_count` transactions each: here the bins remain the base ones;
      3. outside the body, to the left and to the right, the transactions are grouped in
         wider bins (see _merge_sparse_tail) so that every point represents at least
         `min_count` of them; the width is the real width of the bin, not a multiple of
         the grid;
      4. every point is placed at the geometric mean of the prices of the transactions it
         represents (their real position), not at the bin centre.

    density=True  : y = count / (N * bin width in p), with N = number of transactions
                    inside bin_edges (normalization: area 1 with respect to p).
    density=False : y = count * (width of a base bin) / bin width.
                    In the body it is the raw count; in the merged bins it is the count
                    brought back to the width of a base bin (so it stays comparable with
                    the body, but it is no longer an integer).

    min_count=None : plain fixed bins: one point for every non-empty base bin, placed at
                     the geometric centre of the bin.

    Returns a DataFrame with one row per point, sorted by increasing P:
    'P' (position, EUR), 'P_low' and 'P_high' (bin edges, EUR), 'width' (bin width in p),
    'count' (transactions represented), 'y' (plotted value).
    """
    cols = ["P", "P_low", "P_high", "width", "count", "y"]
    edges = np.asarray(bin_edges, dtype=float)
    n_base = len(edges) - 1
    w_ref = (edges[-1] - edges[0]) / n_base  # width of the base bins (equally spaced in p)

    if min_count is not None and min_count < 1:
        raise ValueError("min_count must be >= 1 (or None for plain fixed bins).")

    lp = np.asarray(log_prices, dtype=float)
    lp = lp[np.isfinite(lp) & (lp >= edges[0]) & (lp <= edges[-1])]
    if lp.size == 0:
        return pd.DataFrame({c: np.array([], dtype=float) for c in cols})
    n_in = int(lp.size)

    # same convention as np.histogram: bins closed on the left, the last one also on the right
    idx = np.minimum(np.searchsorted(edges, lp, side="right") - 1, n_base - 1)
    counts = np.bincount(idx, minlength=n_base)
    # sum of the offsets from the left edge of the own bin (they are >= 0): the mean
    # position inside the bin is computed relative to the edge, without losing precision
    # and always staying inside the bin
    offs = np.bincount(idx, weights=lp - edges[idx], minlength=n_base)

    rows = []  # (n, position in p, low edge, high edge), for increasing p
    if min_count is None:
        for j in np.flatnonzero(counts):
            rows.append((int(counts[j]), 0.5 * (edges[j] + edges[j + 1]), edges[j], edges[j + 1]))
    else:
        k = int(min_count)
        mode = int(np.argmax(counts))
        if counts[mode] >= k:
            # body: contiguous base bins, around the most populated one, with >= k transactions
            j0 = mode
            while j0 > 0 and counts[j0 - 1] >= k:
                j0 -= 1
            j1 = mode
            while j1 < n_base - 1 and counts[j1 + 1] >= k:
                j1 += 1
            core = [
                (int(counts[j]), edges[j] + offs[j] / counts[j], edges[j], edges[j + 1])
                for j in range(j0, j1 + 1)
            ]

            right_rows, left_rows = [], []
            is_right = idx > j1
            if is_right.any():
                right_rows = _merge_sparse_tail(
                    lp[is_right], inner=edges[j1 + 1], far=edges[idx[is_right].max() + 1],
                    min_count=k, min_width=w_ref,
                )
            is_left = idx < j0
            if is_left.any():
                # left tail: same algorithm on the sign-changed values (outwards = p
                # decreasing), then the result is flipped
                mirrored = _merge_sparse_tail(
                    -lp[is_left], inner=-edges[j0], far=-edges[idx[is_left].min()],
                    min_count=k, min_width=w_ref,
                )
                left_rows = [(n, -m, -hi, -lo) for n, m, lo, hi in mirrored][::-1]
            rows = left_rows + core + right_rows
        else:
            # no base bin reaches k transactions (small sample): the whole curve is "tail"
            # and is grouped from left to right
            rows = _merge_sparse_tail(
                lp, inner=edges[idx.min()], far=edges[idx.max() + 1],
                min_count=k, min_width=w_ref,
            )

    n = np.array([r[0] for r in rows], dtype=float)
    pos = np.array([r[1] for r in rows], dtype=float)
    lo = np.array([r[2] for r in rows], dtype=float)
    hi = np.array([r[3] for r in rows], dtype=float)
    width = hi - lo
    y = n / (n_in * width) if density else n * w_ref / width
    return pd.DataFrame(
        {"P": np.exp(pos), "P_low": np.exp(lo), "P_high": np.exp(hi),
         "width": width, "count": n.astype(int), "y": y}
    )


def _tail_fit_points(
    prices: np.ndarray,
    threshold: float,
    n_total: Optional[int] = None,
    equal_weight_bins: bool = False,
    n_bins: int = 15,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Prepare the tail data for the log-log regression. Returns
    (tail, survival, P_fit, S_fit):

      - tail, survival : ALL the points above the threshold (sorted) and their empirical
                         survival P(X >= P) (the black curve of the panel);
      - P_fit, S_fit   : the points on which the regression is actually done.

    equal_weight_bins=False: (P_fit, S_fit) = (tail, survival), i.e. every point of the
    tail is a point of the fit. Since the points are much denser close to the threshold
    than in the extreme tail, the fit is dominated by the region near the threshold and
    the last points (few, but the most informative about the tail) weigh very little.

    equal_weight_bins=True: the tail [first point above the threshold, P_max] is split
    into `n_bins` bins LOG-SPACED in P and every bin contributes ONE single point: the
    geometric centre of the bin, with the EXACT empirical survival evaluated there
    (fraction of transactions with price >= centre, normalized as above). All the bins,
    hence all the regions of the tail, have the same weight in the regression, regardless
    of how many points they hold (empty bins are no problem either: the survival is
    cumulative, so it is always > 0).
    """
    tail = np.sort(prices[prices > threshold])
    n = len(tail)
    if n < 20:
        raise ValueError("Too few points above the threshold for the tail fit.")
    n_ref = n if n_total is None else n_total
    # empirical survival function P(X >= P_i) ~ (n-i)/n_ref (rank-based)
    survival = (n - np.arange(n)) / n_ref

    if not equal_weight_bins:
        return tail, survival, tail, survival

    if n_bins < 3:
        raise ValueError("tail_n_bins must be >= 3 for the equal-weight-bins fit.")
    if not tail[-1] > tail[0]:
        raise ValueError("Degenerate tail (all the prices are equal): cannot create the bins.")
    edges = np.geomspace(tail[0], tail[-1], n_bins + 1)
    P_fit = np.sqrt(edges[:-1] * edges[1:])  # geometric centre of every log bin
    # number of points with price >= P_fit  (searchsorted 'left' counts those < P_fit)
    S_fit = (n - np.searchsorted(tail, P_fit, side="left")) / n_ref
    return tail, survival, P_fit, S_fit


def fit_power_law_tail(
    prices: np.ndarray,
    threshold: float = 500_000.0,
    n_total: Optional[int] = None,
    equal_weight_bins: bool = TAIL_FIT_EQUAL_WEIGHT_BINS,
    n_bins: int = TAIL_FIT_N_BINS,
) -> dict:
    """
    Fit of the power-law tail P(X>P) ~ P^{-mu} for P > threshold (in the text of the
    paper: density ~ P^{-1-mu} => the survival function goes as P^{-mu}). The fit is a
    linear regression of log(empirical survival) against log(P) in the tail.

    n_total: if None (default) the survival is normalized by the number of points of the
    tail only (it equals 1 at the threshold). If given (typically the number of ALL the
    transactions), the survival is the true P(X > P) over the whole sample (as in the
    inset of Fig. 7). The slope, hence mu, is identical in the two cases: only the
    intercept changes.

    equal_weight_bins / n_bins: see _tail_fit_points. With equal_weight_bins=False the
    regression uses all the points of the tail; with True it uses `n_bins` log-spaced
    bins of equal weight, so that the extreme tail is not overwhelmed by the many points
    close to the threshold. The default is TAIL_FIT_EQUAL_WEIGHT_BINS.
    """
    tail, _, P_fit, S_fit = _tail_fit_points(
        prices, threshold, n_total=n_total, equal_weight_bins=equal_weight_bins, n_bins=n_bins
    )
    log_p = np.log(P_fit)
    log_s = np.log(S_fit)
    mask = np.isfinite(log_s)
    slope, intercept = np.polyfit(log_p[mask], log_s[mask], 1)
    mu = -slope
    return {
        "mu": mu,
        "slope": slope,
        "intercept": intercept,
        "threshold": threshold,
        "n_tail": len(tail),
        "equal_weight_bins": bool(equal_weight_bins),
        "n_fit_points": int(mask.sum()),
    }


def plot_price_distribution_by_year(
    df: pd.DataFrame,
    year_col: str = "annee",
    n_bins: int = 150,
    tail_threshold: float = 500_000.0,
    title: str = "Distribution of log-prices by year",
    save_path: Optional[str] = None,
    log_price_col: str = "log_price",
    x_label: str = "P",
    y_label: str = r"",
    inset_y_label: str = r"$\mathbb{P}(X>P)$",
    survival_over_all: bool = False,
    tail_equal_weight_bins: bool = TAIL_FIT_EQUAL_WEIGHT_BINS,
    tail_n_bins: int = TAIL_FIT_N_BINS,
    min_count_per_bin: Optional[int] = DIST_MIN_COUNT_PER_BIN,
):
    """
    Replica of Fig. 2 of the paper: density of the log-price per year (normalized
    separately for each year), overlaid with one colour per year. The panel with the
    power-law tail fit (fit of mu) is drawn OUTSIDE the main plot, to the right of the
    colorbar, in separate axes: it can therefore never cover the data points.

    Parameters to reuse the function on other quantities (e.g. the price per square
    meter, see plot_price_m2_distribution_by_year):
      log_price_col     : column of `df` with the log-price to histogram
      x_label, y_label  : labels of the axes of the main plot
                          (x_label is also used for the x axis of the panel)
      inset_y_label     : y label of the fit panel
      survival_over_all : if True the curve of the panel is P(X>P) over the WHOLE sample
                          (not only the tail), as in Fig. 7

    Tail-fit options:
      tail_equal_weight_bins : False -> the regression uses all the points of the tail
                          (dominated by the points close to the threshold, which are the
                          most numerous). True -> the tail is split into `tail_n_bins`
                          log-spaced bins and every bin weighs the same in the fit, so
                          the last points are given an appropriate weight. The bins used
                          are shown in the panel as orange dots.
                          Default: TAIL_FIT_EQUAL_WEIGHT_BINS.
      tail_n_bins        : number of bins (only with tail_equal_weight_bins=True)

    Histogram bins (n_bins are the BASE bins):
      min_count_per_bin  : minimum number of transactions per point. In the tails the
                          base bins that are too poor are merged and every point is
                          placed at the mean position of the actual prices, so the points
                          do not line up on horizontal stripes (integer counts 1, 2,
                          3...). None = fixed bins. See DIST_MIN_COUNT_PER_BIN and
                          _adaptive_log_histogram.
    """
    import matplotlib.pyplot as plt
    from matplotlib.cm import viridis
    from matplotlib.colors import Normalize

    log_p_all = df[log_price_col].to_numpy()
    bin_edges = np.linspace(log_p_all.min(), log_p_all.max(), n_bins + 1)

    years = sorted(df[year_col].dropna().unique())
    norm = Normalize(vmin=min(years), vmax=max(years))

    # --- tail fit on the aggregated distribution (done BEFORE creating the figure, so
    # the layout depends on whether the fit succeeded) ---
    prices_all = np.exp(log_p_all)
    n_total = len(prices_all) if survival_over_all else None
    try:
        tail_fit = fit_power_law_tail(
            prices_all,
            threshold=tail_threshold,
            n_total=n_total,
            equal_weight_bins=tail_equal_weight_bins,
            n_bins=tail_n_bins,
        )
    except ValueError as exc:
        warnings.warn(f"Tail fit not performed: {exc}")
        tail_fit = None

    # --- layout: [main plot | colorbar | tail-fit panel (OUTSIDE)] ---
    if tail_fit is not None:
        fig = plt.figure(figsize=(13, 6), layout="constrained")
        gs = fig.add_gridspec(1, 3, width_ratios=[1, 0.035, 0.62])
        ax = fig.add_subplot(gs[0, 0])
        cax = fig.add_subplot(gs[0, 1])
        # the panel takes the central (vertical) part of the right column
        ax_in = fig.add_subplot(gs[0, 2].subgridspec(3, 1, height_ratios=[1, 2.4, 1])[1, 0])
    else:
        fig = plt.figure(figsize=(9, 6), layout="constrained")
        gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.035])
        ax = fig.add_subplot(gs[0, 0])
        cax = fig.add_subplot(gs[0, 1])
        ax_in = None

    for year in years:
        sub = df.loc[df[year_col] == year, log_price_col].to_numpy()
        hist = _adaptive_log_histogram(sub, bin_edges, density=True, min_count=min_count_per_bin)
        ax.scatter(hist["P"], hist["y"], s=14, color=viridis(norm(year)), label=None)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)

    sm = plt.cm.ScalarMappable(cmap=viridis, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("year")
    cbar.set_ticks([int(y) for y in years])  # integer years only (no 2021.5, ...)

    # --- tail-fit panel, OUTSIDE the main plot ---
    if ax_in is not None:
        tail, survival, P_fit, S_fit = _tail_fit_points(
            prices_all,
            tail_threshold,
            n_total=n_total,
            equal_weight_bins=tail_equal_weight_bins,
            n_bins=tail_n_bins,
        )
        ax_in.plot(tail, survival, color="black", lw=1.5, label="data (tail)")
        if tail_equal_weight_bins:
            ax_in.scatter(
                P_fit, S_fit, s=24, facecolor="tab:orange", edgecolor="black",
                linewidths=0.6, zorder=3, label=f"fitted bins (N={len(P_fit)})",
            )
        fit_line = np.exp(tail_fit["intercept"]) * tail ** tail_fit["slope"]
        ax_in.plot(
            tail, fit_line, "r--", label=rf"$\mu \approx {tail_fit['mu']:.2f}$"
        )
        ax_in.set_xscale("log")
        ax_in.set_yscale("log")
        ax_in.set_xlabel(x_label)
        ax_in.set_ylabel(inset_y_label)
        ax_in.set_title(
            "Power-law tail fit\n"
            + ("(equal-weight bins)" if tail_equal_weight_bins else "(all tail points)"),
            fontsize=10,
        )
        ax_in.grid(True, which="major", ls="-", alpha=0.25)
        ax_in.legend(fontsize=8)

    # (layout="constrained" takes care of the spacing between plot, colorbar and panel)
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"[plot_price_distribution_by_year] figure saved to: {save_path}")
    return fig, tail_fit


def plot_price_distribution_by_month(
    df: pd.DataFrame,
    date_col: str = "date_mutation",
    n_bins: int = 150,
    title: str = "Distribution of log-prices by calendar month (raw counts, NOT normalised)",
    save_path: Optional[str] = None,
    min_count_per_bin: Optional[int] = DIST_MIN_COUNT_PER_BIN,
):
    """
    Distribution of the log-prices computed per CALENDAR MONTH (January, February, ...,
    December), aggregating together the same month of ALL the available years (e.g. all
    the "Januaries" 2018-2022 end up in the same curve) - NOT normalized (density=False):
    the height of the points therefore directly reflects the VOLUME of transactions in
    that calendar month, which makes it possible to spot the seasonality of the market
    (e.g. months that are typically more/less active) that a per-curve normalization
    would hide.

    For readability every month is drawn as a cloud of POINTS (not a continuous line,
    which tends to overlap and become confusing with 12 curves), coloured with a
    continuous GRADIENT scale (cyclic colormap 'twilight', consistent with the cyclic
    nature of the months) and a colorbar with the month labels instead of a 12-entry
    legend.

    Bins (n_bins are the BASE bins): `min_count_per_bin` is the minimum number of
    transactions per point. In the body of the distribution the bins are the base ones
    and y is the raw count, as with plain fixed bins. In the tails the bins that are too
    poor are merged (otherwise the points line up on horizontal stripes, y = 1, 2, 3...)
    and y becomes the count brought back to the width of a base bin, so it stays
    comparable with the body but it is no longer an integer; every point is placed at the
    mean position of the actual prices. None = fixed bins. See DIST_MIN_COUNT_PER_BIN and
    _adaptive_log_histogram.

    Returns (fig, volumes), where `volumes` has the number of transactions per month.
    """
    import matplotlib.pyplot as plt
    from matplotlib.cm import twilight
    from matplotlib.colors import Normalize

    out = df.copy()
    out["_month"] = pd.to_datetime(out[date_col], errors="coerce").dt.month
    out = out.dropna(subset=["_month"])
    out["_month"] = out["_month"].astype(int)

    log_p_all = out["log_price"].to_numpy()
    bin_edges = np.linspace(log_p_all.min(), log_p_all.max(), n_bins + 1)

    month_names = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ]
    # CYCLIC colormap (twilight): January and December have close colours, consistent
    # with the cyclic nature of the months of the year
    norm = Normalize(vmin=1, vmax=12)

    fig, ax = plt.subplots(figsize=(10, 7))
    for month in range(1, 13):
        sub = out.loc[out["_month"] == month, "log_price"].to_numpy()
        hist = _adaptive_log_histogram(sub, bin_edges, density=False, min_count=min_count_per_bin)
        ax.scatter(
            hist["P"],
            hist["y"],
            color=twilight(norm(month)),
            s=24,
            alpha=0.85,
            edgecolors="none",
            zorder=3,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("P (EUR)")
    ax.set_ylabel("number of transactions per bin (raw count, all years summed)")
    ax.set_title(title)
    ax.grid(True, which="major", ls="-", alpha=0.25, zorder=0)
    ax.grid(True, which="minor", ls=":", alpha=0.15, zorder=0)

    # gradient scale (colorbar) instead of a 12-entry legend
    sm = plt.cm.ScalarMappable(cmap=twilight, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02, ticks=range(1, 13))
    cbar.ax.set_yticklabels(month_names)
    cbar.set_label("month")

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"[plot_price_distribution_by_month] figure saved to: {save_path}")

    volumes = out.groupby("_month").size().rename("n_transactions").reset_index()
    volumes["month"] = volumes["_month"].map(lambda m: month_names[m - 1])
    volumes = volumes[["_month", "month", "n_transactions"]].sort_values("_month")
    return fig, volumes


def compute_monthly_volume_and_count(
    df: pd.DataFrame,
    date_col: str = "date_mutation",
    value_col: str = "valeur_fonciere",
) -> pd.DataFrame:
    """
    Aggregate ALL the transactions by CALENDAR MONTH (January..December) and compute the
    AVERAGE ANNUAL values of each month by dividing by the number of unique years analysed
    in the dataset.

    Returns a DataFrame with columns month_num, month, n_transactions and total_volume
    (the last two are average annual values; total_volume is the sum of `value_col`).
    """
    out = df.copy()
    out["_month"] = pd.to_datetime(out[date_col], errors="coerce").dt.month
    out = out.dropna(subset=["_month"])
    out["_month"] = out["_month"].astype(int)

    # Number of unique years analysed. If the 'annee' column exists we use it, otherwise
    # we extract the year from the date
    if "annee" in out.columns and out["annee"].notna().sum() > 0:
        n_years = out["annee"].nunique()
    else:
        n_years = pd.to_datetime(out[date_col], errors="coerce").dt.year.nunique()

    # Safety against division by zero
    if n_years == 0:
        n_years = 1

    print(f"[compute_monthly_volume_and_count] Data divided by {n_years} analysed years.")

    # Raw totals
    agg = (
        out.groupby("_month")[value_col]
        .agg(n_transactions="size", total_volume="sum")
        .reindex(range(1, 13))
        .fillna(0.0)
    )

    # Divide by the number of years to get the annual average
    agg["n_transactions"] = agg["n_transactions"] / n_years
    agg["total_volume"] = agg["total_volume"] / n_years

    month_names = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ]
    agg = agg.reset_index().rename(columns={"_month": "month_num"})
    agg["month"] = agg["month_num"].map(lambda m: month_names[m - 1])

    return agg[["month_num", "month", "n_transactions", "total_volume"]]


def plot_monthly_volume_and_count(
    df: pd.DataFrame,
    date_col: str = "date_mutation",
    value_col: str = "valeur_fonciere",
    title: str = "Monetary volume and number of transactions - Yearly averages by calendar month",
    save_path: Optional[str] = None,
):
    """
    Histogram with every CALENDAR MONTH on the x axis and, on the y axes, the AVERAGE
    ANNUAL values of the monetary volume (bars, left axis) and of the number of
    transactions (line, right axis).

    Returns (fig, agg), where `agg` is the output of compute_monthly_volume_and_count.
    """
    import matplotlib.pyplot as plt

    agg = compute_monthly_volume_and_count(df, date_col=date_col, value_col=value_col)

    fig, ax1 = plt.subplots(figsize=(9, 6))

    color_vol = "tab:blue"
    ax1.bar(
        agg["month"], agg["total_volume"] / 1e9,
        color=color_vol, alpha=0.75, zorder=2,
        label="Average annual monetary volume",
    )
    ax1.set_xlabel("Month")
    ax1.set_ylabel("Average annual monetary volume (billion EUR)", color=color_vol)
    ax1.tick_params(axis="y", labelcolor=color_vol)
    ax1.grid(axis="y", ls="--", alpha=0.3, zorder=0)

    ax2 = ax1.twinx()
    color_cnt = "tab:red"
    line, = ax2.plot(
        agg["month"], agg["n_transactions"],
        color=color_cnt, marker="o", lw=2, markersize=7, zorder=3,
        label="Average annual number of transactions",
    )
    ax2.set_ylabel("Average annual number of transactions", color=color_cnt)
    ax2.tick_params(axis="y", labelcolor=color_cnt)
    ax2.set_ylim(bottom=0)
    ax1.set_ylim(bottom=0)

    bars_handle = plt.Rectangle((0, 0), 1, 1, color=color_vol, alpha=0.75)
    ax1.legend(
        [bars_handle, line],
        ["Average annual monetary volume", "Average annual number of transactions"],
        loc="upper left",
        fontsize=9,
    )

    ax1.set_title(title)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"[plot_monthly_volume_and_count] figure saved to: {save_path}")
    return fig, agg


# ---------------------------------------------------------------------------
# 6. Distribution of the log-prices: Paris vs Creuse (replica of Fig. 6 of the paper)
# ---------------------------------------------------------------------------
#
# Fig. 6 shows the density of the log-price p := log(P) (drawn on an x axis in P,
# log-log) of all the transactions of the 5 years of DVF data for two departments chosen
# as typical examples of a city (Paris) and of the countryside (Creuse): two very
# different shapes, a narrow hump around ~4e5 EUR versus a wide, flat distribution with a
# maximum around ~1e5 EUR. Their superposition explains the double hump of the
# distribution for the whole of France (Fig. 2).
#
# As in Fig. 2 the histogram is in p = log(P) and is normalized (area = 1 with respect to
# p), once per department. All the years are summed together: this is equivalent to
# averaging the distributions of the single years, weighted by the number of transactions
# of each year.
#
# The version on the prices per square meter (same plot with p^ = log(P^)) is
# plot_paris_creuse_m2_distribution, in section 11.


def plot_paris_creuse_distribution(
    df: pd.DataFrame,
    departements: dict[str, str] = PARIS_CREUSE_DEPARTEMENTS,
    colors: Sequence[str] = ("black", "blue"),
    dept_col: str = "code_departement",
    log_price_col: str = "log_price",
    n_bins: int = PARIS_CREUSE_N_BINS,
    price_range: tuple[float, float] = PARIS_CREUSE_PRICE_RANGE,
    y_lim: Optional[tuple[float, float]] = PARIS_CREUSE_Y_LIM,
    title: str = "",
    y_label: str = "",
    latex_style: bool = True,
    figsize: tuple[float, float] = (4.4, 3.2),
    save_path: Optional[str] = None,
    x_label: str = r"$P$",
    legend_loc: str = "upper right",
    min_count_per_bin: Optional[int] = DIST_MIN_COUNT_PER_BIN,
):
    """
    Replica of Fig. 6 of the paper: density of the log-price p = log(P) of ALL the
    transactions (all the years together) for two departments, Paris (black) and Creuse
    (blue), on log-log axes with P on the x axis.

    `df` must be the DataFrame already cleaned by clean_dvf and contain `dept_col` and
    `log_price_col`. To have the tail up to 1e8 EUR as in the figure it must be cleaned
    with price_bounds=PARIS_CREUSE_PRICE_RANGE (today this equals the default of
    clean_dvf):
        clean_pc = clean_dvf(raw, price_bounds=PARIS_CREUSE_PRICE_RANGE)

    Parameters:
      departements : {legend label: code_departement}; default Paris "75" and Creuse "23"
                     (constant PARIS_CREUSE_DEPARTEMENTS)
      colors       : one colour per department, in the same order
      n_bins       : number of base bins, equally spaced in p = log(P) (i.e. log-spaced in
                     P), between price_range[0] and price_range[1]. With few bins the
                     curve is smoother, with many it is noisier. In the tails the bins
                     that are too poor are merged (see min_count_per_bin)
      price_range  : (P_min, P_max) in EUR: range of the histogram and of the x axis
      y_lim        : limits of the y axis (default as in the figure: from 5e-4 to 1);
                     None = automatic
      title, y_label : empty by default, as in the paper's figure
      latex_style  : Computer Modern font (via mathtext, LaTeX is NOT required) to
                     resemble the paper's figure; False = matplotlib style
      figsize      : size in inches (the width/height ratio of the default axes is that
                     of the original figure)
      x_label      : label of the x axis (default P; for the prices per m^2 see
                     plot_paris_creuse_m2_distribution)
      legend_loc   : position of the legend
      min_count_per_bin : minimum number of transactions per point. In the tails the base
                     bins that are too poor are merged and every point is placed at the
                     mean position of the actual prices, so the points do not line up on
                     horizontal stripes (integer counts 1, 2, 3...). None = fixed bins.
                     See DIST_MIN_COUNT_PER_BIN and _adaptive_log_histogram.

    Returns (fig, curves). `curves` is a DataFrame with one row per plotted POINT (the x
    values differ from one curve to the other, so the format is "long"):
    'departement' (label), 'P' (position, EUR), 'P_low' and 'P_high' (bin edges, EUR),
    'width' (bin width in p), 'count' (transactions represented by the point) and
    'density' (density of p, the plotted one).
    """
    import matplotlib.pyplot as plt

    for col in (dept_col, log_price_col):
        if col not in df.columns:
            raise ValueError(
                f"Column '{col}' missing from df: it is needed for the Paris vs Creuse "
                "plot (df must be the output of clean_dvf)."
            )
    if n_bins < 2:
        raise ValueError("n_bins must be >= 2.")
    if not 0 < price_range[0] < price_range[1]:
        raise ValueError("price_range must be (P_min, P_max) with 0 < P_min < P_max.")

    # base bins equally spaced in p = log(P), as in the other distributions of the module
    bin_edges = np.linspace(np.log(price_range[0]), np.log(price_range[1]), n_bins + 1)

    dept_codes = df[dept_col].astype(str).str.strip()
    curve_parts = []  # one DataFrame of points per department (the x values change from curve to curve)
    plotted = []  # (name, colour) of the departments that have data
    for i, (name, code) in enumerate(departements.items()):
        sub = df.loc[dept_codes == str(code), log_price_col].to_numpy(dtype=float)
        print(f"[plot_paris_creuse_distribution] {name} (dept. {code}): {len(sub):,} transactions.")
        if len(sub) == 0:
            warnings.warn(
                f"No transactions for {name} (code_departement='{code}'): "
                "curve omitted. Check the DVF files loaded."
            )
            continue
        hist = _adaptive_log_histogram(sub, bin_edges, density=True, min_count=min_count_per_bin)
        print(f"[plot_paris_creuse_distribution]   -> {len(hist):,} points drawn.")
        curve_parts.append(hist.rename(columns={"y": "density"}).assign(departement=name))
        plotted.append((name, colors[i % len(colors)]))
    if not plotted:
        raise ValueError(
            "No transactions found for the requested departments: "
            f"{departements}. Check the column '{dept_col}'."
        )
    curves = pd.concat(curve_parts, ignore_index=True)[
        ["departement", "P", "P_low", "P_high", "width", "count", "density"]
    ]

    # "Paper" style: Computer Modern through mathtext/cmr10 (included in matplotlib, no
    # LaTeX needed). The tick fonts are resolved at drawing time, so the drawing
    # (draw/savefig) stays INSIDE the rc context.
    rc = {}
    if latex_style:
        rc = {
            "font.family": "serif",
            "font.serif": ["cmr10", "DejaVu Serif"],
            "mathtext.fontset": "cm",
            "axes.formatter.use_mathtext": True,
            "axes.unicode_minus": False,  # cmr10 has no glyph for the unicode minus
            "font.size": 12,
        }
    with plt.rc_context(rc):
        fig, ax = plt.subplots(figsize=figsize, layout="constrained")
        for name, color in plotted:
            pts = curves[curves["departement"] == name]
            ax.scatter(
                pts["P"], pts["density"],
                s=20, color=color, edgecolors="none", label=name,
            )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(*price_range)
        if y_lim is not None:
            ax.set_ylim(*y_lim)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        if title:
            ax.set_title(title)
        ax.legend(loc=legend_loc, scatterpoints=1)

        fig.canvas.draw()  # instantiates ticks and labels with the fonts of the rc context
        if save_path:
            fig.savefig(save_path, dpi=150)
            print(f"[plot_paris_creuse_distribution] figure saved to: {save_path}")
    return fig, curves


# ---------------------------------------------------------------------------
# 7. Spatial variogram computation
# ---------------------------------------------------------------------------

@dataclass
class VariogramResult:
    ell: np.ndarray          # centres of the distance bins (km)
    V: np.ndarray            # V(l) = <(p(r)-p(r'))^2>
    n_pairs: np.ndarray      # number of pairs used in every bin
    bin_edges: np.ndarray


def compute_variogram(
    field: pd.DataFrame,
    n_bins: int = 30,
    ell_min: Optional[float] = None,
    ell_max: Optional[float] = None,
    bin_edges: Optional[np.ndarray] = None,
    drop_empty_bins: bool = True,
    max_exact_n: int = 12_000,
    n_sampled_pairs: int = 5_000_000,
    random_state: int = 0,
    bin_center: str = "geometric",
) -> VariogramResult:
    """
    Compute the empirical spatial variogram V(l) = <(p(r)-p(r'))^2>_{|r-r'|=l} from the
    field p(r) (one row per hexagonal cell).

    For scalability:
      - if the number of cells N <= max_exact_n, ALL the C(N,2) pairs are used (exact
        computation, via scipy.spatial.distance.pdist);
      - otherwise the variogram is estimated with a Monte Carlo sampling of
        n_sampled_pairs random pairs (standard in the geostatistics literature for large
        datasets).

    Distances are computed on the Lambert-93 projected plane (x_km, y_km), which has
    negligible distortion over the whole of metropolitan France.

    bin_center controls how the centre of each distance bin is computed (used as the
    "ell" abscissa of the point): "geometric" (default, suited to LOGARITHMIC bins, such
    as those generated by default when bin_edges=None) or "arithmetic" (suited to LINEAR
    bins, such as those passed explicitly by the multi-scale variogram, see
    MULTISCALE_BIN_SPECS below).
    """
    from scipy.spatial.distance import pdist

    x = field["x_km"].to_numpy()
    y = field["y_km"].to_numpy()
    p = field["mean_log_price"].to_numpy()
    n = len(field)

    if n < 10:
        raise ValueError("Too few cells with enough data for a variogram.")

    if n <= max_exact_n:
        coords = np.column_stack([x, y])
        dist = pdist(coords, metric="euclidean")
        # squared difference: (p_i - p_j)^2 for all the pairs i<j
        diff2 = pdist(p.reshape(-1, 1), metric="sqeuclidean")
        print(f"[compute_variogram] EXACT computation on {n:,} cells ({len(dist):,} pairs).")
    else:
        rng = np.random.default_rng(random_state)
        i = rng.integers(0, n, size=n_sampled_pairs)
        j = rng.integers(0, n, size=n_sampled_pairs)
        mask = i != j
        i, j = i[mask], j[mask]
        dist = np.hypot(x[i] - x[j], y[i] - y[j])
        diff2 = (p[i] - p[j]) ** 2
        print(
            f"[compute_variogram] Monte Carlo sampling: {len(dist):,} pairs "
            f"out of {n:,} cells."
        )

    if bin_edges is None:
        ell_min = ell_min or max(dist[dist > 0].min(), 0.5)
        ell_max = ell_max or dist.max()
        bin_edges = np.logspace(np.log10(ell_min), np.log10(ell_max), n_bins + 1)
    else:
        bin_edges = np.asarray(bin_edges)
        n_bins = len(bin_edges) - 1

    bin_idx = np.digitize(dist, bin_edges) - 1
    V = np.full(n_bins, np.nan)
    npairs = np.zeros(n_bins, dtype=int)
    if bin_center == "geometric":
        ell_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])  # geometric centre (for log bins)
    elif bin_center == "arithmetic":
        ell_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])  # arithmetic centre (for linear bins)
    else:
        raise ValueError("bin_center must be 'geometric' or 'arithmetic'.")

    for b in range(n_bins):
        sel = (bin_idx == b) & (dist > 0)
        npairs[b] = sel.sum()
        if npairs[b] > 0:
            V[b] = diff2[sel].mean()

    if drop_empty_bins:
        valid = npairs > 0
        return VariogramResult(
            ell=ell_centers[valid], V=V[valid], n_pairs=npairs[valid], bin_edges=bin_edges
        )
    # the empty bins are not filtered: this is needed to align several curves (multi-scale)
    # on the same bins when they are aggregated (mean/deviation across groups)
    return VariogramResult(ell=ell_centers, V=V, n_pairs=npairs, bin_edges=bin_edges)


def fit_log_slope(
    result: VariogramResult, ell_range: tuple[float, float]
) -> dict:
    """
    Linear fit of V(l) against log(l) in the interval ell_range=(lmin,lmax), to estimate
    the slope A^2/(2*pi*D) predicted by Eq. (3) of the paper:

        V(l) ~= slope * log(l) + intercept ,   slope = A^2 / (2*pi*D)
    """
    lo, hi = ell_range
    mask = (result.ell >= lo) & (result.ell <= hi)
    if mask.sum() < 3:
        raise ValueError("Too few points in the specified fit range.")

    log_ell = np.log(result.ell[mask])
    V = result.V[mask]
    slope, intercept = np.polyfit(log_ell, V, 1)
    A2_over_D = 2 * np.pi * slope

    return {
        "slope_V_vs_logl": slope,
        "intercept": intercept,
        "A2_over_D": A2_over_D,
        "n_points_fit": int(mask.sum()),
        "ell_range_fit": ell_range,
    }


# ---------------------------------------------------------------------------
# 8. Plot (Fig. 3 style: log x axis, linear y axis)
# ---------------------------------------------------------------------------

def plot_variogram(
    result: VariogramResult,
    fit: Optional[dict] = None,
    title: str = "Spatial variogram of log-prices",
    save_path: Optional[str] = None,
):
    """Plot a spatial variogram, with the log fit (from fit_log_slope) if given."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(result.ell, result.V, s=25, color="black", label="data (DVF, all years)")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\ell$ (km)")
    ax.set_ylabel(r"$\langle (p(r,t)-p(r',t))^2\rangle_{|r-r'|=\ell}$")
    ax.set_title(title)

    if fit is not None:
        lo, hi = fit["ell_range_fit"]
        ell_line = np.logspace(np.log10(lo), np.log10(hi), 50)
        V_line = fit["slope_V_vs_logl"] * np.log(ell_line) + fit["intercept"]
        ax.plot(
            ell_line,
            V_line,
            "r--",
            label=(
                rf"log fit: $A^2/D \approx {fit['A2_over_D']:.2f}$"
            ),
        )

    ax.legend()
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"[plot_variogram] figure saved to: {save_path}")
    return fig


# ---------------------------------------------------------------------------
# 9. MULTI-scale spatial variogram (replica of Fig. 3 of the paper)
# ---------------------------------------------------------------------------
#
# The paper shows on the same plot the variogram computed at 4 scales (city / departments /
# regions / the whole of France), each with the cross-sectional variability among the
# units of the same scale drawn as a shaded band around the mean.
#
# For "France" there is a single unit (the whole country) -> no band.
# For "region"/"department"/"city" there are several units (the single regions/
# departments/cities): we compute one variogram per unit and then the mean +/- standard
# deviation across the units, on the same shared distance bins (so the curves are
# comparable point by point).

# Mapping department (2-digit INSEE code, 2A/2B for Corsica) -> French administrative
# region (2016 reform, 13 metropolitan regions + overseas departments). Used to group the
# departments into regions when computing the variogram at the "region" scale.
DEPARTEMENT_TO_REGION: dict[str, str] = {
    **{d: "Auvergne-Rhone-Alpes" for d in ["01", "03", "07", "15", "26", "38", "42", "43", "63", "69", "73", "74"]},
    **{d: "Bourgogne-Franche-Comte" for d in ["21", "25", "39", "58", "70", "71", "89", "90"]},
    **{d: "Bretagne" for d in ["22", "29", "35", "56"]},
    **{d: "Centre-Val de Loire" for d in ["18", "28", "36", "37", "41", "45"]},
    **{d: "Corse" for d in ["2A", "2B"]},
    **{d: "Grand Est" for d in ["08", "10", "51", "52", "54", "55", "57", "67", "68", "88"]},
    **{d: "Hauts-de-France" for d in ["02", "59", "60", "62", "80"]},
    **{d: "Ile-de-France" for d in ["75", "77", "78", "91", "92", "93", "94", "95"]},
    **{d: "Normandie" for d in ["14", "27", "50", "61", "76"]},
    **{d: "Nouvelle-Aquitaine" for d in ["16", "17", "19", "23", "24", "33", "40", "47", "64", "79", "86", "87"]},
    **{d: "Occitanie" for d in ["09", "11", "12", "30", "31", "32", "34", "46", "48", "65", "66", "81", "82"]},
    **{d: "Pays de la Loire" for d in ["44", "49", "53", "72", "85"]},
    **{d: "Provence-Alpes-Cote d'Azur" for d in ["04", "05", "06", "13", "83", "84"]},
    "971": "Guadeloupe", "972": "Martinique", "973": "Guyane",
    "974": "La Reunion", "976": "Mayotte",
}

# Specification of the LINEAR distance bins (constant width in km, NOT logarithmic) used to
# compute the variogram of each scale, as in Fig. 3 of the paper: every scale
# (city/department/region/France) has its own range of ell (km) and is split into
# independent linear bins, so a different resolution can be chosen for each scale (narrower
# bins where more points are needed, e.g. at the "city" scale, and wider where the range
# is large, e.g. at the "france" scale).
#
# >>> CHANGE HERE the number/width of the bins for each scale <<<
# For every scale give a dictionary with:
#   - "ell_min", "ell_max" : distance range (km) over which the bins are created
#   - "bin_width"          : width (km) of each linear bin, OR
#   - "n_bins"             : number of linear bins (alternative to bin_width)
# If both are given, "bin_width" takes precedence over "n_bins".
# The plot (plot_multiscale_variogram) keeps a LOGARITHMIC x axis as in Fig. 3: with linear
# bins the points appear denser and denser as ell grows (just like in the figure of the
# paper), because a bin of constant width in km covers a smaller and smaller fraction of a
# logarithmic decade at large ell.
MULTISCALE_BIN_SPECS: dict[str, dict] = {
    "city":        {"ell_min": 1.0,  "ell_max": 12.0,   "bin_width": 0.5},
    "departement": {"ell_min": 3.0,  "ell_max": 160.0,  "bin_width": 4.0},
    "region":      {"ell_min": 5.0,  "ell_max": 400.0,  "bin_width": 8.0},
    "france":      {"ell_min": 15.0, "ell_max": 1000.0, "bin_width": 12.0},
}

# Fallback spec used for a scale that does not appear in MULTISCALE_BIN_SPECS
# (e.g. a custom scale added by the user).
_DEFAULT_BIN_SPEC: dict = {"ell_min": 1.0, "ell_max": 1000.0, "bin_width": 10.0}


def linear_bin_edges(
    ell_min: float,
    ell_max: float,
    n_bins: Optional[int] = None,
    bin_width: Optional[float] = None,
) -> np.ndarray:
    """
    Build LINEAR bin_edges (constant width) between ell_min and ell_max.

    Specify either `n_bins` (number of bins) or `bin_width` (width in km of each bin: the
    last bin is automatically extended up to ell_max if (ell_max - ell_min) is not an
    exact multiple of bin_width). If both are given, `bin_width` takes precedence.
    """
    if ell_max <= ell_min:
        raise ValueError("ell_max must be greater than ell_min.")
    if bin_width is not None:
        if bin_width <= 0:
            raise ValueError("bin_width must be positive.")
        n = max(1, int(np.ceil((ell_max - ell_min) / bin_width)))
        edges = ell_min + np.arange(n + 1) * bin_width
        edges[-1] = max(edges[-1], ell_max)  # covers exactly up to ell_max
        return edges
    if n_bins is not None:
        if n_bins < 1:
            raise ValueError("n_bins must be >= 1.")
        return np.linspace(ell_min, ell_max, n_bins + 1)
    raise ValueError("Specify n_bins or bin_width in MULTISCALE_BIN_SPECS.")


MULTISCALE_COLORS = {
    "city": "tab:blue",
    "departement": "tab:green",
    "region": "tab:red",
    "france": "black",
}

# Names (in English) shown in the plot legends for each scale.
SCALE_DISPLAY_NAMES = {
    "city": "city",
    "departement": "department",
    "region": "region",
    "france": "France",
}

# Range of ell (km) over which to fit (log-linear fit of V against log(ell)) SEPARATELY
# for each scale, as in Fig. 3 of the paper (where the dashed reference lines cover
# different ranges for city / department / region / France). CHANGE the values
# (ell_min, ell_max in km) here to adapt them to the actual extent of your own curves: the
# useful range can be read from where the points of each scale are not NaN in the plot
# produced by plot_multiscale_variogram.
MULTISCALE_FIT_RANGES: dict[str, tuple] = {
    "city": (1.0, 8.0),
    "departement": (10.0, 75.0),
    "region": (13.0, 110.0),
    "france": (15.0, 230.0),
}


def _group_series_for_scale(df: pd.DataFrame, scale: str, top_k_cities: int = 17) -> pd.Series:
    """Return a Series aligned with df, with the group label to use for the requested
    scale (one group = one variogram that enters the average)."""
    if scale == "france":
        return pd.Series("France", index=df.index)
    if scale == "departement":
        if "code_departement" not in df.columns:
            raise ValueError("Column 'code_departement' missing from the DataFrame.")
        return df["code_departement"]
    if scale == "region":
        if "code_departement" not in df.columns:
            raise ValueError("Column 'code_departement' missing from the DataFrame.")
        return df["code_departement"].map(DEPARTEMENT_TO_REGION)
    if scale == "city":
        if "code_commune" not in df.columns:
            raise ValueError("Column 'code_commune' missing from the DataFrame.")
        counts = df["code_commune"].value_counts()
        top = counts.head(top_k_cities).index
        # NOTE: approximation. The paper uses 17 ad hoc urban areas; here the top_k_cities
        # communes (code_commune) with the most transactions are used as a proxy for the
        # "cities". For a more faithful replica, call compute_group_variograms with a
        # custom group_series that groups the communes that really make up the urban
        # area of interest (e.g. the 20 arrondissements of Paris under a single id
        # "Paris").
        return df["code_commune"].where(df["code_commune"].isin(top))
    raise ValueError(f"unknown scale: {scale}")


def compute_group_variograms(
    df: pd.DataFrame,
    resolution: int,
    group_series: pd.Series,
    bin_edges: np.ndarray,
    min_transactions_per_cell: int = MIN_TRANSACTIONS_PER_CELL,
    min_cells_per_group: int = 10,
    log_price_col: str = "log_price",
) -> dict:
    """Compute a variogram (on the same LINEAR bin_edges, shared among the groups of this
    scale) for each distinct value of group_series, returning {group: VariogramResult}."""
    if log_price_col not in df.columns:
        raise ValueError(
            f"Column '{log_price_col}' missing from df "
            "(for the prices per m^2 it is created by add_price_per_m2)."
        )
    data = df.loc[group_series.notna()].copy()
    groups = group_series.loc[data.index]

    curves = {}
    for g in sorted(groups.dropna().unique()):
        sub = data.loc[groups == g]
        try:
            field = build_price_field(
                sub,
                resolution=resolution,
                min_transactions=min_transactions_per_cell,
                log_price_col=log_price_col,
            )
        except Exception:
            continue  # a group whose field cannot be built is skipped
        if len(field) < min_cells_per_group:
            continue
        try:
            vres = compute_variogram(
                field, bin_edges=bin_edges, drop_empty_bins=False, bin_center="arithmetic"
            )
        except ValueError:
            continue
        curves[g] = vres
    print(f"[compute_group_variograms] {len(curves)}/{groups.nunique()} groups with enough data.")
    return curves


def aggregate_group_variograms(curves: dict, bin_edges: np.ndarray) -> pd.DataFrame:
    """Mean (+ standard deviation) across the curves of the various groups, bin by bin.

    bin_edges is the array of LINEAR bins used to compute all the curves in `curves` (the
    same bins for every group of the scale): the centre of each bin is therefore the
    ARITHMETIC centre (consistent with bins of constant width).
    """
    ell_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    n_bins = len(ell_centers)
    stacked = np.full((len(curves), n_bins), np.nan)
    for i, vres in enumerate(curves.values()):
        stacked[i, :] = vres.V

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        V_mean = np.nanmean(stacked, axis=0)
        V_std = np.nanstd(stacked, axis=0)
        n_groups_with_data = np.sum(~np.isnan(stacked), axis=0)

    out = pd.DataFrame(
        {"ell": ell_centers, "V_mean": V_mean, "V_std": V_std, "n_groups": n_groups_with_data}
    )
    out.attrs["bin_edges"] = bin_edges  # used by fit_multiscale_slopes
    return out


def compute_multiscale_variogram(
    df: pd.DataFrame,
    scales: Sequence[str] = ("city", "departement", "region", "france"),
    bin_specs: dict[str, dict] = MULTISCALE_BIN_SPECS,
    min_transactions_per_cell: int = MIN_TRANSACTIONS_PER_CELL,
    min_cells_per_group: int = 10,
    top_k_cities: int = 17,
    log_price_col: str = "log_price",
) -> dict[str, pd.DataFrame]:
    """Run compute_group_variograms + aggregate_group_variograms for every requested
    scale, with LINEAR distance bins that are INDEPENDENT for each scale (number/width
    configurable in `bin_specs`, by default MULTISCALE_BIN_SPECS: change that dictionary
    earlier in the code to change the resolution/range of each scale). Returns
    {scale: DataFrame with ell, V_mean, V_std, n_groups}.

    log_price_col: column with the log-price whose variogram is computed (default
    'log_price'; 'log_price_m2' for the price per square meter)."""
    out = {}
    for scale in scales:
        spec = bin_specs.get(scale, _DEFAULT_BIN_SPEC)
        if scale not in bin_specs:
            warnings.warn(
                f"No bin spec configured for the scale '{scale}' in "
                f"MULTISCALE_BIN_SPECS: using the fallback {_DEFAULT_BIN_SPEC}."
            )
        edges = linear_bin_edges(**spec)

        group_series = _group_series_for_scale(df, scale, top_k_cities=top_k_cities)
        curves = compute_group_variograms(
            df,
            resolution=SCALE_TO_H3_RES[scale],
            group_series=group_series,
            bin_edges=edges,
            min_transactions_per_cell=min_transactions_per_cell,
            min_cells_per_group=min_cells_per_group,
            log_price_col=log_price_col,
        )
        if not curves:
            warnings.warn(f"No valid group for the scale '{scale}': skipped.")
            continue
        out[scale] = aggregate_group_variograms(curves, bin_edges=edges)
    return out


def fit_multiscale_slopes(
    agg_results: dict[str, pd.DataFrame],
    fit_ranges: dict[str, tuple] = MULTISCALE_FIT_RANGES,
) -> dict[str, dict]:
    """
    Run fit_log_slope (log-linear fit of V against log(ell)) SEPARATELY for each scale
    present in agg_results, each in its own range of ell (km) defined in fit_ranges
    (default: MULTISCALE_FIT_RANGES, to be changed directly in the code to suit your
    data). Returns {scale: fit dict} (scales without a configured range, or without enough
    points in the range, are skipped with a warning).
    """
    fits = {}
    for scale, dfres in agg_results.items():
        if scale not in fit_ranges:
            warnings.warn(f"No fit range configured for the scale '{scale}': skipped.")
            continue
        valid = dfres["V_mean"].notna()
        pseudo_result = VariogramResult(
            ell=dfres.loc[valid, "ell"].to_numpy(),
            V=dfres.loc[valid, "V_mean"].to_numpy(),
            n_pairs=dfres.loc[valid, "n_groups"].to_numpy(),
            bin_edges=dfres.attrs.get("bin_edges", np.array([])),
        )
        try:
            fits[scale] = fit_log_slope(pseudo_result, ell_range=fit_ranges[scale])
        except ValueError as exc:
            warnings.warn(f"Fit failed for the scale '{scale}' in the range {fit_ranges[scale]}: {exc}")
    return fits


def plot_multiscale_variogram(
    agg_results: dict[str, pd.DataFrame],
    fit_ranges: dict[str, tuple] = MULTISCALE_FIT_RANGES,
    title: str = "Multi-scale spatial variogram",
    save_path: Optional[str] = None,
    ylabel: str = r"$\langle (p(r,t)-p(r',t))^2\rangle_{|r-r'|=\ell}$",
):
    """
    Draw, on the same plot, the mean variogram (+ cross-sectional variability band) of
    each scale, in the style of Fig. 3 of the paper: blue=city, green=departments,
    red=regions, black=France.

    The data of each scale are computed (upstream, in compute_multiscale_variogram) on
    LINEAR distance bins that are independent for each scale (see MULTISCALE_BIN_SPECS):
    the x axis of the plot stays logarithmic as in Fig. 3, so the points appear denser
    and denser at large ell, exactly as in the figure of the paper.

    For EACH scale a log-linear fit (fit_log_slope) is also run in its own range of ell,
    drawn as a dashed line of the same colour as the scale, limited to that range (like
    the reference lines of Fig. 3, which cover different ranges for each scale). The
    ranges are changed by passing a fit_ranges dictionary (default:
    MULTISCALE_FIT_RANGES, defined earlier in the code).

    Returns (fig, fits).
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 6))
    fits = fit_multiscale_slopes(agg_results, fit_ranges=fit_ranges)

    for scale, dfres in agg_results.items():
        valid = dfres["V_mean"].notna()
        color = MULTISCALE_COLORS.get(scale, None)
        ax.scatter(
            dfres.loc[valid, "ell"], dfres.loc[valid, "V_mean"], s=18, color=color,
            label=SCALE_DISPLAY_NAMES.get(scale, scale),
        )
        has_band = dfres.loc[valid, "V_std"].fillna(0).sum() > 0
        if has_band:
            ax.fill_between(
                dfres.loc[valid, "ell"],
                dfres.loc[valid, "V_mean"] - dfres.loc[valid, "V_std"],
                dfres.loc[valid, "V_mean"] + dfres.loc[valid, "V_std"],
                color=color,
                alpha=0.2,
            )

        if scale in fits:
            f = fits[scale]
            lo, hi = f["ell_range_fit"]
            ell_line = np.logspace(np.log10(lo), np.log10(hi), 50)
            V_line = f["slope_V_vs_logl"] * np.log(ell_line) + f["intercept"]
            ax.plot(
                ell_line,
                V_line,
                "--",
                color=color,
                lw=1.5,
                label=rf"fit {SCALE_DISPLAY_NAMES.get(scale, scale)}: $A^2/D \approx {f['A2_over_D']:.2f}$",
            )

    ax.set_xscale("log")
    ax.set_xlabel(r"$\ell$ (km)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"[plot_multiscale_variogram] figure saved to: {save_path}")
    return fig, fits


# ---------------------------------------------------------------------------
# 10. TEMPORAL variogram (V1, V2) + theoretical fit (Eq. 1 / SM-2.3)
# ---------------------------------------------------------------------------
#
# In the original paper the temporal variogram is NOT computed on the DVF data (which
# cover only 5 years) but on the 1970-2022 time series (reference [10] of the paper). Here
# it is computed on the DVF data anyway, choosing a "reasonable" spatial cell scale: it
# must be large enough to have a sufficient number of transactions per cell IN EVERY
# period (otherwise the monthly mean per cell is too noisy/empty), but not so large that
# very few cells remain to be averaged over in space. By default the "region" scale is
# used (H3 resolution 6, ~36 km^2), configurable.
#
# Definitions (identical to the paper):
#   V1(tau) = < (p(r,t) - p(r,t+tau))^2 >_{r,t}                    (Eq. in the text)
#   V2(tau) = < (p(r,t)-pbar(t) - (p(r,t+tau)-pbar(t+tau)))^2 >_{r,t}
# where pbar(t) = <p(r,t)>_r is the spatial mean at time t.
#
# Theoretical prediction: SPDE with diffusion + noise with memory T, spatial cutoff S,
# kappa -> 0 (see below), with T and S free in the fit and B = A^2/D free (projected
# analytically) or fixed by hand. Optional nugget c>=0.


def build_spatiotemporal_field(
    df: pd.DataFrame,
    resolution: int,
    freq: str = "M",
    min_transactions: int = MIN_TRANSACTIONS_PER_CELL,
) -> pd.DataFrame:
    """
    Build a pivot table hex_id x period, with the mean log-price per cell in each period
    (freq: 'M' monthly, 'Q' quarterly, 'Y' yearly). The cells of a (hex_id, period) pair
    with fewer than min_transactions transactions are set to NaN (too noisy).
    """
    hexed = assign_hex_cells(df, resolution)
    hexed["_period"] = pd.to_datetime(hexed["date_mutation"], errors="coerce").dt.to_period(freq)
    hexed = hexed.dropna(subset=["_period"])

    grouped = (
        hexed.groupby(["hex_id", "_period"])
        .agg(mean_log_price=("log_price", "mean"), n=("log_price", "size"))
        .reset_index()
    )
    grouped.loc[grouped["n"] < min_transactions, "mean_log_price"] = np.nan

    pivot = grouped.pivot(index="hex_id", columns="_period", values="mean_log_price")
    full_range = pd.period_range(pivot.columns.min(), pivot.columns.max(), freq=freq)
    pivot = pivot.reindex(columns=full_range)

    n_cells, n_periods = pivot.shape
    coverage = pivot.notna().mean().mean()
    print(
        f"[build_spatiotemporal_field] {n_cells:,} cells x {n_periods} periods "
        f"({freq}), average coverage {coverage:.1%}."
    )
    return pivot


@dataclass
class TemporalVariogramResult:
    tau_years: np.ndarray
    V: np.ndarray
    n_pairs: np.ndarray


def _lagged_variogram(values: np.ndarray, periods_per_year: int) -> TemporalVariogramResult:
    """Core shared by V1 and V2: given an array (n_cells, n_periods) compute
    <(x(t)-x(t+tau))^2> averaging over the available cells and times t (ignoring the NaNs
    due to cells/periods without enough data)."""
    n_periods = values.shape[1]
    max_lag = n_periods - 1
    V = np.full(max_lag, np.nan)
    npairs = np.zeros(max_lag, dtype=np.int64)

    for lag in range(1, max_lag + 1):
        diff2 = (values[:, lag:] - values[:, : n_periods - lag]) ** 2
        valid = ~np.isnan(diff2)
        n_valid = valid.sum()
        npairs[lag - 1] = n_valid
        if n_valid > 0:
            V[lag - 1] = np.nansum(diff2) / n_valid

    tau_years = np.arange(1, max_lag + 1) / periods_per_year
    ok = npairs > 0
    return TemporalVariogramResult(tau_years=tau_years[ok], V=V[ok], n_pairs=npairs[ok])


_PERIODS_PER_YEAR = {"M": 12, "Q": 4, "Y": 1}


def compute_V1(pivot: pd.DataFrame, freq: str = "M") -> TemporalVariogramResult:
    """V1(tau) = <(p(r,t)-p(r,t+tau))^2>_{r,t} computed on the whole field."""
    return _lagged_variogram(pivot.to_numpy(dtype=float), _PERIODS_PER_YEAR[freq])


def compute_V2(pivot: pd.DataFrame, freq: str = "M") -> TemporalVariogramResult:
    """
    V2(tau): like V1 but after subtracting, at every time, the spatial mean
    pbar(t) = <p(r,t)>_r (removes the common national trend, isolating the idiosyncratic
    local fluctuations).
    """
    period_mean = pivot.mean(axis=0, skipna=True)
    deviations = pivot.subtract(period_mean, axis=1)
    return _lagged_variogram(deviations.to_numpy(dtype=float), _PERIODS_PER_YEAR[freq])


# Theoretical model of the temporal variogram (Becharat, Benzaquen, Bouchaud,
# arXiv:2412.14624).
#
#     d psi(r,t)/dt = D * Laplacian psi - kappa * psi + eta(r,t) + xi(r)
#
# - xi(r) is STATIC (quenched) spatial disorder: in the difference psi(r,t) - psi(r,t+tau)
#   at a fixed point r it cancels identically, so it does NOT enter the temporal
#   variogram (it only matters for the spatial variogram).
# - eta(r,t) is a noise of intensity A^2 with temporal MEMORY T (exponential correlation
#   ~ e^{-|t-t'|/T}).
#
# For every Fourier mode k, with relaxation rate lambda = D k^2 + kappa, the stationary
# correlation function is (A^2/2) * G(lambda, tau), with
#
#     G(lam, tau) = [ T e^{-tau/T} - e^{-lam*tau}/lam ] / (T^2 lam^2 - 1)
#
# (removable singularity at T*lam = 1). Integrating over space k with an ultraviolet
# cutoff 1/a (a^2 = D*S, a = typical size of the aggregation cell) and in the limit
# kappa -> 0 (1/kappa ~ 1e3 years >> measured lags):
#
#     V(tau) = B/(4*pi) * Int_0^{1/S} dlam [ G(lam,0) - G(lam,tau) ]  (+ c)
#
# with B = A^2/D (same convention as the paper: spatial V = B/(2 pi) log l) and
# c = 2*sigma_measurement^2 an optional nugget (measurement noise, c >= 0).
#
# Parameters IDENTIFIABLE from the temporal variogram: B, T, S (and c).
# D alone is NOT identifiable.
#
# Note: kappa is sent to zero and T is the memory time of the noise (not 1/kappa), so the
# variogram does not saturate to a plateau but grows slowly (logarithmically) over the
# measured lag range.

# Composite Gauss-Legendre quadrature on u = lambda*S in [0,1], geometric panels
_GL_X, _GL_W = np.polynomial.legendre.leggauss(16)
_GL_EDGES = np.concatenate([[0.0], np.geomspace(1e-5, 1.0, 48)])
_GL_A, _GL_B = _GL_EDGES[:-1, None], _GL_EDGES[1:, None]
_GL_U = (0.5 * (_GL_B - _GL_A) * _GL_X + 0.5 * (_GL_B + _GL_A)).ravel()
_GL_WU = (0.5 * (_GL_B - _GL_A) * _GL_W).ravel()


def V_shape(tau: np.ndarray, T: float, S: float) -> np.ndarray:
    """V(tau)/B = (1/4pi) * Int_0^{1/S} [G(lam,0) - G(lam,tau)] dlam  (kappa = 0)."""
    tau = np.atleast_1d(np.asarray(tau, dtype=float))[:, None]
    lam, w = _GL_U / S, _GL_WU / S
    x = T * lam
    G0 = 1.0 / (lam * (x + 1.0))
    den = x * x - 1.0
    s = tau / T
    with np.errstate(divide="ignore", invalid="ignore"):
        Gt = (T * np.exp(-s) - np.exp(-lam * tau) / lam) / den
    near = np.abs(den) < 1e-7  # removable singularity at T*lam = 1
    if near.any():
        Gt = np.where(near[None, :], T * np.exp(-s) * (1 + s) / 2, Gt)
    return (G0 - Gt) @ w / (4 * np.pi)


def V_model(tau: np.ndarray, B: float, T: float, S: float, c: float = 0.0) -> np.ndarray:
    """Theoretical variogram: B * V_shape(tau, T, S) + c  (c = nugget = 2*sigma_measurement^2)."""
    return B * V_shape(tau, T, S) + c


def V_theory_temporal(
    tau: np.ndarray, A2_over_D: float, T: float, S: float, c: float = 0.0
) -> np.ndarray:
    """
    Alias of V_model with the legacy signature (B = A^2/D).
    Vectorized over tau (tau in years). Returns a float if tau is a scalar.
    """
    V = V_model(tau, A2_over_D, T, S, c)
    return V if V.size > 1 else float(V[0])


def _solve_linear(v: np.ndarray, V: np.ndarray, B_fixed: Optional[float], nugget: bool):
    """For fixed (T,S), B (and c>=0) are LINEAR: NNLS least squares with RELATIVE
    residuals (V_mod - V)/V. Returns (B, c, relative residuals)."""
    from scipy.optimize import nnls

    w = 1.0 / V
    if B_fixed is None:
        cols = [v * w] + ([w] if nugget else [])
        th, _ = nnls(np.column_stack(cols), np.ones_like(V))
        B, c = th[0], (th[1] if nugget else 0.0)
    else:
        B = B_fixed
        y = 1.0 - B * v * w
        c = nnls(w[:, None], y)[0][0] if nugget else 0.0
    return B, c, (B * v + c) * w - 1.0


def fit_temporal_variogram(
    tau_data: np.ndarray,
    V_data: np.ndarray,
    A2_over_D: Optional[float] = None,
    nugget: bool = False,
    T_bounds: tuple = (0.1, 100.0),
    S_bounds: tuple = (0.01, 100.0),
    grid_T: tuple = (0.2, 40.0),
    grid_S: tuple = (0.03, 40.0),
    grid: int = 80,
    n_refine: int = 6,
) -> dict:
    """
    Fit of the theoretical model (V_model, see above) to the data V(tau), tau in years.

    The optimization is over (log T, log S), with B = A^2/D (and the nugget c >= 0, if
    nugget=True) PROJECTED analytically (NNLS) at every (T,S). RELATIVE residuals
    (V_mod - V)/V, consistent with the log-log scale of the plot.

      A2_over_D = None  -> B free (estimated)
      A2_over_D = value -> B fixed (e.g. 2*pi*spatial slope, or slider)

    The cost landscape has MULTIPLE LOCAL MINIMA (even without noise): we start from all
    the local minima of a log-spaced (T,S) grid and keep the best one. The result
    includes the list of the minima found ('candidates': list of (relative rms, T, S)) to
    see the T-S degeneracy.
    """
    from scipy.optimize import least_squares
    from scipy.ndimage import minimum_filter

    tau = np.asarray(tau_data, dtype=float)
    V = np.asarray(V_data, dtype=float)
    ok = np.isfinite(tau) & np.isfinite(V) & (V > 0)
    tau, V = tau[ok], V[ok]
    if tau.size < 3:
        raise ValueError("Too few valid points (finite V and > 0) for the temporal variogram fit.")

    B_fixed = A2_over_D

    def resid(z):
        T, S = np.exp(z)
        return _solve_linear(V_shape(tau, T, S), V, B_fixed, nugget)[2]

    Tg = np.geomspace(max(grid_T[0], T_bounds[0]), min(grid_T[1], T_bounds[1]), grid)
    Sg = np.geomspace(max(grid_S[0], S_bounds[0]), min(grid_S[1], S_bounds[1]), grid)
    cmap = np.array([[np.sum(resid(np.log([T, S])) ** 2) for S in Sg] for T in Tg])
    ismin = cmap == minimum_filter(cmap, size=3, mode="nearest")
    cells = np.argwhere(ismin)[np.argsort(cmap[ismin])[:n_refine]]

    lo = np.log([T_bounds[0], S_bounds[0]])
    hi = np.log([T_bounds[1], S_bounds[1]])
    cands = []
    for i, j in cells:
        sol = least_squares(
            resid, np.log([Tg[i], Sg[j]]), bounds=(lo, hi), ftol=1e-13, xtol=1e-13, gtol=1e-13
        )
        cands.append((sol.cost, sol.x))
    cands.sort(key=lambda c_: c_[0])
    uniq = []
    for k, z in cands:  # remove duplicates (same minimum reached from several cells)
        if all(np.max(np.abs(z - z2)) > 1e-2 for _, z2 in uniq):
            uniq.append((k, z))
    cands = uniq

    T, S = np.exp(cands[0][1])
    B, c, r = _solve_linear(V_shape(tau, T, S), V, B_fixed, nugget)
    V_fit = B * V_shape(tau, T, S) + c
    return {
        "A2_over_D": float(B),  # = B (legacy key; B is fixed if given, otherwise estimated)
        "B": float(B),
        "T": float(T),
        "S": float(S),
        "c": float(c),
        "rms_rel": float(np.sqrt(np.mean(r ** 2))),
        "rmse": float(np.sqrt(np.mean((V - V_fit) ** 2))),
        "B_fixed": B_fixed is not None,
        # True if T or S ended up on a limit of the bounds: unreliable fit (model
        # unsuited to the data or bounds too narrow)
        "at_bound": bool(
            np.any(np.abs(np.log([T, S]) - lo) < 1e-3) or np.any(np.abs(np.log([T, S]) - hi) < 1e-3)
        ),
        "candidates": [
            (float(np.sqrt(2 * k / len(V))), float(np.exp(z)[0]), float(np.exp(z)[1])) for k, z in cands
        ],
    }


def plot_temporal_variogram(
    V1_result: TemporalVariogramResult,
    V2_result: Optional[TemporalVariogramResult] = None,
    fit: Optional[dict] = None,
    title: str = "Temporal variogram (V1, V2) - DVF data",
    save_path: Optional[str] = None,
):
    """STATIC version (no slider) of the plot of the two temporal variograms with the
    theoretical fit (Eq. 4) overlaid - useful for reports/non-interactive use. For manual
    adjustment of A^2/D see plot_temporal_variogram_interactive."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(V1_result.tau_years, V1_result.V, color="black", s=20, label=r"$\mathbb{V}_1(\tau)$ (data)")
    if V2_result is not None:
        ax.scatter(
            V2_result.tau_years, V2_result.V, color="tab:blue", marker="x", s=25,
            label=r"$\mathbb{V}_2(\tau)$ (data)",
        )

    if fit is not None:
        tau_line = np.geomspace(V1_result.tau_years.min(), V1_result.tau_years.max(), 200)
        c_fit = fit.get("c", 0.0)
        V_line = V_model(tau_line, fit["A2_over_D"], fit["T"], fit["S"], c_fit)
        label = (
            rf"theoretical fit: $A^2/D={fit['A2_over_D']:.2f}$, "
            rf"$T={fit['T']:.2f}$ yr, $S={fit['S']:.2f}$ yr"
        )
        if c_fit > 0:
            label += rf", $c={c_fit:.4f}$"
        ax.plot(tau_line, V_line, "r--", label=label)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\tau$ (years)")
    ax.set_ylabel(r"$\mathbb{V}(0,\tau)$")
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"[plot_temporal_variogram] figure saved to: {save_path}")
    return fig


def plot_temporal_variogram_interactive(
    V1_result: TemporalVariogramResult,
    V2_result: Optional[TemporalVariogramResult] = None,
    A2_over_D_init: float = 1.2,
    A2_over_D_range: tuple = (0.01, 5.0),
    nugget: bool = False,
    fit_on: str = "V1",
):
    """
    INTERACTIVE plot (needs an interactive matplotlib backend, e.g. running the script
    locally with `plt.show()`: it does not work in headless/non-interactive notebook
    environments) with a SLIDER for A^2/D: every time the slider changes, T and S are
    re-fitted automatically (same model and same multi-start fit as
    fit_temporal_variogram, with B = A^2/D fixed at the slider value, optional nugget)
    and the theoretical curve is updated accordingly. fit_on ('V1' or 'V2') chooses which
    variogram is fitted (V2 only if V2_result is given). For speed the interactive re-fit
    uses a coarser (T,S) grid (30x30).

    Usage:
        fig, slider = plot_temporal_variogram_interactive(v1, v2)
        plt.show()
    """
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider

    fig, ax = plt.subplots(figsize=(8, 6))
    plt.subplots_adjust(bottom=0.25)

    ax.scatter(V1_result.tau_years, V1_result.V, color="black", s=20, label=r"$\mathbb{V}_1(\tau)$ (data)")
    if V2_result is not None:
        ax.scatter(
            V2_result.tau_years, V2_result.V, color="tab:blue", marker="x", s=25,
            label=r"$\mathbb{V}_2(\tau)$ (data)",
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\tau$ (years)")
    ax.set_ylabel(r"$\mathbb{V}(0,\tau)$")
    ax.set_title("Temporal variogram - manual tuning of $A^2/D$ (exact model)")

    tau_line = np.geomspace(V1_result.tau_years.min(), V1_result.tau_years.max(), 150)

    target_result = V2_result if (fit_on == "V2" and V2_result is not None) else V1_result
    fit_name = "V2" if target_result is not V1_result else "V1"

    def _refit(a2d: float) -> dict:
        return fit_temporal_variogram(
            target_result.tau_years, target_result.V, a2d, nugget=nugget, grid=30, n_refine=3
        )

    initial_fit = _refit(A2_over_D_init)
    (line,) = ax.plot(
        tau_line,
        V_model(
            tau_line, initial_fit["A2_over_D"], initial_fit["T"], initial_fit["S"], initial_fit["c"]
        ),
        "r--",
        label=f"theoretical fit (on {fit_name})",
    )
    ax.legend(fontsize=8)

    info_text = fig.text(0.15, 0.12, "", fontsize=9)

    ax_slider = plt.axes([0.2, 0.05, 0.6, 0.03])
    slider = Slider(
        ax_slider, r"$A^2/D$", A2_over_D_range[0], A2_over_D_range[1], valinit=A2_over_D_init
    )

    def update(_val):
        a2d = slider.val
        try:
            f = _refit(a2d)
            line.set_ydata(V_model(tau_line, f["A2_over_D"], f["T"], f["S"], f["c"]))
            info_text.set_text(
                f"Current fit:  A²/D = {a2d:.3f} (fixed)   ->   T = {f['T']:.2f} yr,   "
                f"S = {f['S']:.2f} yr,   c = {f['c']:.4f}   (rel. rms = {f['rms_rel']:.3f})"
            )
        except Exception:
            info_text.set_text("Fit does not converge for this value of A²/D.")
        fig.canvas.draw_idle()

    slider.on_changed(update)
    update(A2_over_D_init)

    return fig, slider


# ---------------------------------------------------------------------------
# 11. Prices per SQUARE METER (replica of Fig. 7 and Fig. 5)
# ---------------------------------------------------------------------------
#
# Same analyses as Fig. 2 (distribution) and Fig. 3 (multi-scale variogram), but on the
# price per square meter  P^ := valeur_fonciere / surface_reelle_bati  and on its
# logarithm  p^ := log(P^).
#
#   Fig. 7: distribution of p^ per year (x axis in P^, log-log), with an inset of the
#           power-law tail P^(-1-mu) (mu ~ 1.4 in the paper).
#   Fig. 5: spatial variogram of p^ at 4 scales (city / departments / regions / France),
#           with cross-sectional variability bands.
#   Paris vs Creuse: like Fig. 6 but on p^ (plot_paris_creuse_m2_distribution).
#
# These analyses use ONLY SINGLE-LINE TRANSACTIONS (a single DVF line, which is a Maison or
# an Appartement), see select_single_line_transactions: the total price of a transaction
# with more lines (land, outbuildings) does not match the surface of the dwelling, which
# would bias the price per m^2. The nominal-price analyses (Fig. 2, 3, 6) use the broader
# selection of select_residential_transactions instead. The data of this section is
# therefore prepared with
#
#     clean_m2 = add_price_per_m2(clean_dvf(raw, single_line_only=True))
#
# and NOT from the DataFrame cleaned for the nominal prices.
#
# There is no need to duplicate the logic: the existing functions are parametric on the
# log-price column (log_price_col). Below there are only the data preparation
# (add_price_per_m2) and two thin wrappers that set the titles/labels/thresholds suited to
# the prices per m^2.


def add_price_per_m2(
    df: pd.DataFrame,
    surface_col: str = "surface_reelle_bati",
    min_surface: float = MIN_SURFACE_M2,
    price_m2_bounds: tuple[float, float] = PRICE_M2_BOUNDS,
    require_single_line: bool = True,
) -> pd.DataFrame:
    """
    Add the columns 'price_m2' = valeur_fonciere / surface_reelle_bati and 'log_price_m2'
    = log(price_m2), and discard the unusable rows:

    - missing surface or < min_surface (m^2): the price per m^2 would be undefined or
      spurious (e.g. an "apartment" of 1 m^2 is almost always a classification/recording
      error);
    - price per m^2 outside price_m2_bounds (EUR/m^2): typos and extreme outliers.

    To be called on the DataFrame cleaned with clean_dvf(raw, single_line_only=True),
    which keeps ONLY the single-line transactions: there valeur_fonciere is the whole
    price of ONE dwelling and the surface is that of the dwelling sold. The DataFrame
    cleaned for the nominal prices (clean_dvf(raw)) also contains dwellings sold together
    with land or outbuildings, whose total price would be divided by the surface of the
    dwelling alone.

    require_single_line: if True (default) a ValueError is raised unless `df` comes from
    clean_dvf(..., single_line_only=True) (i.e. it has the column N_LINES_COLUMN equal to
    1 everywhere); False skips the check (at your own risk).

    Returns a copy, `df` is not modified.
    """
    if surface_col not in df.columns:
        raise ValueError(
            f"Column '{surface_col}' missing: it is needed to compute the price per m^2. "
            "Check that it is among the columns loaded by load_dvf_files."
        )
    if require_single_line and (
        N_LINES_COLUMN not in df.columns or not (df[N_LINES_COLUMN] == 1).all()
    ):
        raise ValueError(
            "add_price_per_m2 needs single-line transactions: the total price of a "
            "transaction with more lines (land, outbuildings, ...) does not match the "
            "surface of the dwelling alone. Clean the raw data with "
            "clean_dvf(raw, single_line_only=True) first "
            "(or pass require_single_line=False to skip this check)."
        )
    n0 = len(df)
    out = df.copy()
    surface = pd.to_numeric(out[surface_col], errors="coerce")
    out["price_m2"] = out["valeur_fonciere"].astype(float) / surface

    # NaN (missing surface) and inf (surface 0) both fail the tests
    keep = (surface >= min_surface) & out["price_m2"].between(*price_m2_bounds)
    out = out[keep].copy()
    out["log_price_m2"] = np.log(out["price_m2"])

    print(
        f"[add_price_per_m2] {n0:,} -> {len(out):,} rows with a valid price per m^2 "
        f"({len(out) / max(n0, 1):.1%} kept; surface >= {min_surface:g} m^2, "
        f"price in [{price_m2_bounds[0]:g}, {price_m2_bounds[1]:g}] EUR/m^2)."
    )
    return out


def plot_price_m2_distribution_by_year(
    df: pd.DataFrame,
    year_col: str = "annee",
    n_bins: int = 150,
    tail_threshold: float = PRICE_M2_TAIL_THRESHOLD,
    title: str = "Distribution of log-prices per m$^2$ by year",
    save_path: Optional[str] = None,
    tail_equal_weight_bins: bool = TAIL_FIT_EQUAL_WEIGHT_BINS,
    tail_n_bins: int = TAIL_FIT_N_BINS,
    min_count_per_bin: Optional[int] = DIST_MIN_COUNT_PER_BIN,
):
    """
    Replica of Fig. 7: distribution of p^ = log(P^) (P^ = price per m^2) per year,
    normalized separately for each year and drawn on an x axis in P^ (log-log), with a
    panel (OUTSIDE the plot) of the power-law tail: the black curve is P(X > P^) over the
    whole sample and the dashed red line is the fit P^(-mu), with mu estimated above
    `tail_threshold` (EUR/m^2).

    tail_equal_weight_bins / tail_n_bins: if True the fit uses `tail_n_bins` log-spaced
    bins of equal weight instead of all the points of the tail (see
    plot_price_distribution_by_year / _tail_fit_points). Default:
    TAIL_FIT_EQUAL_WEIGHT_BINS.

    min_count_per_bin: minimum number of transactions per point (adaptive bins in the
    tails, points at the mean position of the actual prices; None = fixed bins): see
    plot_price_distribution_by_year and DIST_MIN_COUNT_PER_BIN.

    `df` must contain 'log_price_m2' (see add_price_per_m2). Returns (fig, tail_fit), with
    tail_fit = dict of the fit (mu, slope, ...).
    """
    return plot_price_distribution_by_year(
        df,
        year_col=year_col,
        n_bins=n_bins,
        tail_threshold=tail_threshold,
        title=title,
        save_path=save_path,
        log_price_col="log_price_m2",
        x_label=r"$\hat{P}$ (EUR/m$^2$)",
        y_label=r"density of $\hat{p} := \log \hat{P}$",
        inset_y_label=r"$\mathbb{P}(X>\hat{P})$",
        survival_over_all=True,
        tail_equal_weight_bins=tail_equal_weight_bins,
        tail_n_bins=tail_n_bins,
        min_count_per_bin=min_count_per_bin,
    )


def plot_paris_creuse_m2_distribution(
    df: pd.DataFrame,
    departements: dict[str, str] = PARIS_CREUSE_DEPARTEMENTS,
    n_bins: int = PARIS_CREUSE_M2_N_BINS,
    price_range: tuple[float, float] = PARIS_CREUSE_M2_PRICE_RANGE,
    y_lim: Optional[tuple[float, float]] = PARIS_CREUSE_M2_Y_LIM,
    title: str = "",
    y_label: str = "",
    save_path: Optional[str] = None,
    min_count_per_bin: Optional[int] = DIST_MIN_COUNT_PER_BIN,
    legend_loc: str = "best",
    **kwargs,
):
    """
    Distribution of p^ = log(P^) (P^ = price per m^2) of ALL the transactions (all the
    years together) for Paris (black) and Creuse (blue), on log-log axes with P^ on the x
    axis: the same plot as plot_paris_creuse_distribution (replica of Fig. 6) but on the
    price per square meter instead of the total price.

    `df` must contain 'log_price_m2' and 'code_departement': it is the output of
    add_price_per_m2 applied to the single-line DataFrame from clean_dvf:
        clean_m2 = add_price_per_m2(clean_dvf(raw, single_line_only=True))
        fig, curves = plot_paris_creuse_m2_distribution(clean_m2, save_path="paris_creuse_m2.png")

    Parameters (default values in the PARIS_CREUSE_M2_* constants at the top of the file):
      departements, n_bins, price_range, y_lim, title, y_label, save_path, min_count_per_bin :
          as in plot_paris_creuse_distribution, but price_range is the range of P^ in
          EUR/m^2, n_bins the number of base bins in p^ and y_lim = None gives an
          automatic y axis
      legend_loc : position of the legend; "best" (default) puts it where it covers the
          fewest points, because here the hump of Paris is on the right and could end up
          under "upper right"
      **kwargs : forwarded to plot_paris_creuse_distribution (colors, dept_col,
          latex_style, figsize)

    Returns (fig, curves) like plot_paris_creuse_distribution: `curves` has one row per
    plotted point ('departement', 'P' = P^ in EUR/m^2, 'P_low', 'P_high', 'width',
    'count', 'density').
    """
    return plot_paris_creuse_distribution(
        df,
        departements=departements,
        log_price_col="log_price_m2",
        n_bins=n_bins,
        price_range=price_range,
        y_lim=y_lim,
        title=title,
        y_label=y_label,
        x_label=r"$\hat{P}$ (EUR/m$^2$)",
        save_path=save_path,
        min_count_per_bin=min_count_per_bin,
        legend_loc=legend_loc,
        **kwargs,
    )


def plot_price_m2_multiscale_variogram(
    agg_results: dict[str, pd.DataFrame],
    fit_ranges: dict[str, tuple] = MULTISCALE_FIT_RANGES,
    title: str = "Multi-scale spatial variogram of log-prices per m$^2$",
    save_path: Optional[str] = None,
):
    """
    Replica of Fig. 5: spatial variogram of p^ = log(P^) (price per m^2) at 4 scales
    (blue=city, green=departments, red=regions, black=France). `agg_results` must be
    computed with
        compute_multiscale_variogram(clean_m2, log_price_col="log_price_m2")
    where clean_m2 = add_price_per_m2(clean_dvf(raw, single_line_only=True)). The distance
    bins, the fit ranges and the rest of the graphics are the same as in
    plot_multiscale_variogram (MULTISCALE_BIN_SPECS / MULTISCALE_FIT_RANGES).
    Returns (fig, fits).
    """
    return plot_multiscale_variogram(
        agg_results,
        fit_ranges=fit_ranges,
        title=title,
        save_path=save_path,
        ylabel=r"$\langle (\hat{p}(r,t)-\hat{p}(r',t))^2\rangle_{|r-r'|=\ell}$",
    )


# ---------------------------------------------------------------------------
# 12. Full pipeline
# ---------------------------------------------------------------------------

def run_full_pipeline(
    file_pattern: str,
    scale: str = "france",
    sep: Optional[str] = None,
    type_locaux: Sequence[str] = ("Maison", "Appartement"),
    min_transactions: int = MIN_TRANSACTIONS_PER_CELL,
    fit_ell_range: Optional[tuple[float, float]] = None,
    save_fig_path: Optional[str] = "variogram.png",
) -> dict:
    """
    Run the whole pipeline on the NOMINAL prices: load all the DVF files matching
    `file_pattern`, clean them (clean_dvf with its default selection), build the field
    p(r) averaged over ALL the years at the chosen scale, compute the spatial variogram
    and fit it with the logarithmic form of Eq. (3).

    scale: one of "city", "departement", "region", "france"
           (see SCALE_TO_H3_RES for the correspondence with the H3 resolution).
    """
    if scale not in SCALE_TO_H3_RES:
        raise ValueError(f"scale must be one of {list(SCALE_TO_H3_RES)}")

    raw = load_dvf_files(file_pattern, sep=sep)
    clean = clean_dvf(raw, type_locaux=type_locaux)
    field = build_price_field(
        clean, resolution=SCALE_TO_H3_RES[scale], min_transactions=min_transactions
    )
    variogram = compute_variogram(field)

    if fit_ell_range is None:
        # by default: exclude the first and the last ~15% of the bins (short-range
        # measurement noise, edge/plateau effects at long range)
        lo = np.percentile(variogram.ell, 15)
        hi = np.percentile(variogram.ell, 85)
        fit_ell_range = (lo, hi)

    fit = fit_log_slope(variogram, ell_range=fit_ell_range)
    fig = plot_variogram(
        variogram,
        fit=fit,
        title=f"Spatial variogram ({scale}, H3 res {SCALE_TO_H3_RES[scale]})",
        save_path=save_fig_path,
    )

    return {
        "raw_count": len(raw),
        "clean_count": len(clean),
        "field": field,
        "variogram": pd.DataFrame(
            {"ell_km": variogram.ell, "V": variogram.V, "n_pairs": variogram.n_pairs}
        ),
        "fit": fit,
        "figure": fig,
        "figure_path": save_fig_path,
    }


# ---------------------------------------------------------------------------
# 13. Command-line interface
# ---------------------------------------------------------------------------

# Analyses on the prices per m^2. They use ONLY single-line transactions (see "DATA
# SELECTION" in the module docstring); all the other analyses use the nominal prices.
M2_ANALYSES: tuple[str, ...] = ("distribution_m2", "multiscale_m2", "paris_creuse_m2", "m2")


def main(argv: Optional[Sequence[str]] = None) -> None:
    """Command-line entry point (see the module docstring for examples)."""
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Spatial and temporal variograms of French housing prices from DVF data "
            "(arXiv:2412.14624). The nominal-price analyses use the dwelling selection of "
            "clean_dvf; the analyses on the prices per m^2 use single-line transactions only."
        ),
        epilog="See the module docstring and the README for the details and more examples.",
    )
    parser.add_argument("file_pattern", help="e.g. '/data/DVF-*.txt'")
    parser.add_argument("--sep", default=None, help="CSV separator (default: auto)")
    parser.add_argument(
        "--analysis",
        default="spatial",
        choices=[
            "spatial", "distribution", "multiscale", "temporal",
            "distribution_m2", "multiscale_m2", "m2", "paris_creuse", "paris_creuse_m2",
        ],
        help=(
            "spatial: spatial variogram at one scale (replica of Fig. 3, one curve). "
            "distribution: price distribution per year (Fig. 2) + per month (volumes). "
            "multiscale: spatial variogram at several overlaid scales (full replica of Fig. 3). "
            "temporal: temporal variogram V1/V2 with theoretical fit. "
            "distribution_m2: distribution of the prices PER M^2 per year (replica of Fig. 7). "
            "multiscale_m2: multi-scale spatial variogram of the prices PER M^2 (replica of Fig. 5). "
            "m2: all the analyses on the prices per m^2 (distribution, Paris vs Creuse, "
            "multi-scale variogram), loading the data once. "
            "paris_creuse: distribution of the log-prices of Paris vs Creuse (replica of Fig. 6). "
            "paris_creuse_m2: the same Paris vs Creuse distribution on the prices per m^2. "
            "The analyses on the prices per m^2 (distribution_m2, multiscale_m2, "
            "paris_creuse_m2, m2) use ONLY single-line transactions."
        ),
    )
    parser.add_argument(
        "--scale", default="france", choices=list(SCALE_TO_H3_RES),
        help="spatial scale (only for --analysis spatial)",
    )
    parser.add_argument(
        "--min-transactions", type=int, default=MIN_TRANSACTIONS_PER_CELL,
        help=f"minimum number of transactions per cell (default: {MIN_TRANSACTIONS_PER_CELL}, "
             "can also be changed once and for all in MIN_TRANSACTIONS_PER_CELL in the code)",
    )
    parser.add_argument("--out-prefix", default="dvf_output", help="prefix of the generated image files")
    parser.add_argument(
        "--temporal-scale", default="region", choices=list(SCALE_TO_H3_RES),
        help="spatial scale of the cells used for the spatio-temporal field (only --analysis temporal)",
    )
    parser.add_argument("--temporal-freq", default="Y", choices=["M", "Q", "Y"])
    parser.add_argument(
        "--a2-over-d", type=float, default=None,
        help="A^2/D fixed for the temporal fit. If omitted, B=A^2/D is FREE "
             "(analytically projected at every (T,S))",
    )
    parser.add_argument(
        "--nugget", action="store_true",
        help="adds to the temporal model a constant c>=0 (measurement noise)",
    )
    parser.add_argument(
        "--fit-target", default="V1", choices=["V1", "V2"],
        help="temporal variogram to fit (default: V1)",
    )
    parser.add_argument(
        "--interactive", action="store_true",
        help="opens the temporal variogram plot with an interactive slider for A^2/D "
             "(needs an interactive matplotlib backend, e.g. running locally with a display)",
    )
    parser.add_argument(
        "--tail-equal-weight-bins", action=argparse.BooleanOptionalAction,
        default=TAIL_FIT_EQUAL_WEIGHT_BINS,
        help="fit of the power-law tail (distributions of the prices and of the prices per m^2) "
             "with log-spaced bins of EQUAL WEIGHT instead of all the points of the tail "
             f"(default: {TAIL_FIT_EQUAL_WEIGHT_BINS}, can be changed in TAIL_FIT_EQUAL_WEIGHT_BINS)",
    )
    parser.add_argument(
        "--tail-n-bins", type=int, default=TAIL_FIT_N_BINS,
        help=f"number of log-spaced bins for the equal-weight tail fit (default: {TAIL_FIT_N_BINS})",
    )
    parser.add_argument(
        "--min-count-per-bin", type=int, default=DIST_MIN_COUNT_PER_BIN or 0,
        help="minimum number of transactions per point in the distribution plots: in the tails "
             "the bins that are too poor are merged and the points are drawn at the mean position "
             "of the actual prices, so they do not line up on horizontal stripes "
             f"(default: {DIST_MIN_COUNT_PER_BIN}, can be changed in DIST_MIN_COUNT_PER_BIN). "
             "0 = fixed bins (points at the bin centre)",
    )
    args = parser.parse_args(argv)
    min_count_per_bin = args.min_count_per_bin if args.min_count_per_bin > 0 else None

    raw = load_dvf_files(args.file_pattern, sep=args.sep)
    if args.analysis in M2_ANALYSES:
        # Prices per m^2: ONLY single-line transactions, counted on the raw data. The
        # cleaning thresholds are in MIN_SURFACE_M2 / PRICE_M2_BOUNDS.
        clean_m2 = add_price_per_m2(clean_dvf(raw, single_line_only=True))
    elif args.analysis == "paris_creuse":
        # Fig. 6: the x axis reaches 1e8 EUR -> explicit price limit, see
        # PARIS_CREUSE_PRICE_RANGE
        clean = clean_dvf(raw, price_bounds=PARIS_CREUSE_PRICE_RANGE)
    else:
        # nominal prices: dwelling selection of clean_dvf
        clean = clean_dvf(raw)

    if args.analysis == "spatial":
        field = build_price_field(clean, resolution=SCALE_TO_H3_RES[args.scale], min_transactions=args.min_transactions)
        vres = compute_variogram(field)
        lo, hi = np.percentile(vres.ell, 15), np.percentile(vres.ell, 85)
        fit = fit_log_slope(vres, ell_range=(lo, hi))
        plot_variogram(vres, fit=fit, title=f"Spatial variogram ({args.scale})", save_path=f"{args.out_prefix}_spatial.png")
        print("\n--- SPATIAL FIT ---")
        for k, v in fit.items():
            print(f"{k}: {v}")

    elif args.analysis == "distribution":
        _, tail_fit = plot_price_distribution_by_year(
            clean,
            save_path=f"{args.out_prefix}_dist_year.png",
            tail_equal_weight_bins=args.tail_equal_weight_bins,
            tail_n_bins=args.tail_n_bins,
            min_count_per_bin=min_count_per_bin,
        )
        print("\n--- TAIL FIT (power law) ---")
        if tail_fit:
            for k, v in tail_fit.items():
                print(f"{k}: {v}")
        _, volumes = plot_price_distribution_by_month(
            clean, save_path=f"{args.out_prefix}_dist_month.png", min_count_per_bin=min_count_per_bin
        )
        print("\n--- MONTHLY VOLUMES (first rows) ---")
        print(volumes.head(12).to_string(index=False))
        _, monthly_agg = plot_monthly_volume_and_count(clean, save_path=f"{args.out_prefix}_volume_month.png")
        print("\n--- MONETARY VOLUME AND NUMBER OF TRANSACTIONS PER MONTH ---")
        print(monthly_agg.to_string(index=False))

    elif args.analysis == "multiscale":
        agg = compute_multiscale_variogram(clean, min_cells_per_group=max(10, args.min_transactions))
        # the fit ranges per scale are changed in MULTISCALE_FIT_RANGES
        # (earlier in the code), not from the command line
        fig, fits = plot_multiscale_variogram(agg, save_path=f"{args.out_prefix}_multiscale.png")
        print("\n--- FIT PER SCALE (see MULTISCALE_FIT_RANGES to change the ranges) ---")
        for scale, f in fits.items():
            print(f"{scale}: {f}")

    elif args.analysis == "temporal":
        pivot = build_spatiotemporal_field(
            clean, resolution=SCALE_TO_H3_RES[args.temporal_scale],
            freq=args.temporal_freq, min_transactions=args.min_transactions,
        )
        v1 = compute_V1(pivot, freq=args.temporal_freq)
        v2 = compute_V2(pivot, freq=args.temporal_freq)
        target = v2 if args.fit_target == "V2" else v1
        fit = fit_temporal_variogram(
            target.tau_years, target.V, A2_over_D=args.a2_over_d, nugget=args.nugget
        )
        mode = "B=A^2/D free" if args.a2_over_d is None else f"B=A^2/D fixed = {args.a2_over_d}"
        print(f"\n--- TEMPORAL FIT on {args.fit_target} ({mode}, nugget={args.nugget}) ---")
        print(f"B = A^2/D = {fit['B']:.4f}   (equivalent spatial slope B/2pi = {fit['B'] / (2 * np.pi):.4f})")
        print(f"T = {fit['T']:.3f} years   S = a^2/D = {fit['S']:.3f} years   c = {fit['c']:.5f}")
        print(f"relative RMS error = {fit['rms_rel']:.4f}   (absolute RMSE = {fit['rmse']:.5f})")
        if fit["at_bound"]:
            print("WARNING: T or S ended up on a limit of the bounds (T in [0.1, 100], S in [0.01, 100] years): "
                  "the fit is not reliable (model unsuited to these data, or lag/cell scale to be changed).")
        print("local minima found (relative rms, T, S):")
        for rms, T_, S_ in fit["candidates"]:
            print(f"    {rms:.4f}  T={T_:.3f}  S={S_:.3f}")
        if args.interactive:
            import matplotlib.pyplot as plt
            B0 = fit["B"] if args.a2_over_d is None else args.a2_over_d
            fig, slider = plot_temporal_variogram_interactive(
                v1, v2, A2_over_D_init=B0,
                A2_over_D_range=(min(0.01, 0.1 * B0), max(5.0, 3.0 * B0)),
                nugget=args.nugget, fit_on=args.fit_target,
            )
            plt.show()
        else:
            plot_temporal_variogram(v1, v2, fit=fit, save_path=f"{args.out_prefix}_temporal.png")

    elif args.analysis == "paris_creuse":
        # distribution of the log-prices of Paris vs Creuse (replica of Fig. 6)
        plot_paris_creuse_distribution(
            clean, save_path=f"{args.out_prefix}_paris_creuse.png", min_count_per_bin=min_count_per_bin
        )

    elif args.analysis in M2_ANALYSES:
        if args.analysis in ("distribution_m2", "m2"):
            _, tail_fit = plot_price_m2_distribution_by_year(
                clean_m2,
                save_path=f"{args.out_prefix}_dist_year_m2.png",
                tail_equal_weight_bins=args.tail_equal_weight_bins,
                tail_n_bins=args.tail_n_bins,
                min_count_per_bin=min_count_per_bin,
            )
            print("\n--- TAIL FIT prices per m^2 (power law) ---")
            if tail_fit:
                for k, v in tail_fit.items():
                    print(f"{k}: {v}")

        if args.analysis in ("paris_creuse_m2", "m2"):
            # distribution of the log-prices PER M^2 of Paris vs Creuse
            plot_paris_creuse_m2_distribution(
                clean_m2, save_path=f"{args.out_prefix}_paris_creuse_m2.png",
                min_count_per_bin=min_count_per_bin,
            )

        if args.analysis in ("multiscale_m2", "m2"):
            agg_m2 = compute_multiscale_variogram(
                clean_m2,
                log_price_col="log_price_m2",
                min_cells_per_group=max(10, args.min_transactions),
            )
            _, fits = plot_price_m2_multiscale_variogram(
                agg_m2, save_path=f"{args.out_prefix}_multiscale_m2.png"
            )
            print("\n--- FIT PER SCALE, prices per m^2 (see MULTISCALE_FIT_RANGES for the ranges) ---")
            for scale, f in fits.items():
                print(f"{scale}: {f}")


if __name__ == "__main__":
    main()