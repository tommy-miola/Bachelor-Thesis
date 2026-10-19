"""
Temporal variogram of real-estate prices
========================================

Fits the temporal variogram of log-prices with the diffusive model of
Becharat, Benzaquen & Bouchaud, "The Diffusive Nature of Housing Prices"
(arXiv:2412.14624).

Model (Eq. 1 of the paper):   dpsi/dt = D lap(psi) - kappa psi + eta + xi

###########################################################################################################

FIT MODEL. For each Fourier mode, with lambda = D k^2 + kappa and noise eta with memory T, the stationary
correlation function is (A^2/2) * G(lambda, tau), with

    G(lam, tau) = [ T e^{-tau/T} - e^{-lam tau}/lam ] / (T^2 lam^2 - 1)

The static field xi cancels out of the temporal variogram. Integrating over k with a UV cutoff
1/a (a^2 = D S) and kappa -> 0 (1/kappa ~ 1e3 years >> lags considered) gives

    V(tau) = B/(4 pi) * Int_0^{1/S} dlam [ G(lam,0) - G(lam,tau) ],     B = A^2 / D

(same convention as the paper, where the spatial variogram is V(l) = B/(2 pi) log l).
Identifiable parameters: B, T, S.  D alone is NOT identifiable from the variogram.

##########################################################################################################

METHOD. The temporal variogram is calculated for every commune over the whole period 1970-2022,
for the selected time lag.
For every time lag, the code calculates the variogram starting from every single year in the time period
(ex. 1970, 1971, 1972, ...) and executes the procedure for every commune. It then averages out the results to obtain
the definite value for that specific time lag.

##########################################################################################################

INPUT: The script takes "data.csv" or "data_no_inflation.csv" as input files.
Extra CLI commands can be given as shown below.

OUTPUT: The output goes to a directory called "output_temporal_variogram"

##########################################################################################################

USAGE:
    python temporal_variogram.py data.csv --target V2 --max-lag 10
    python temporal_variogram.py data_no_inflation.csv --target V1 --B-fixed 1.0 --nugget

    "nugget" adds a constant c to the fit.
    "max-lag" tells how many years the fit has to span.
    As the paper tells, and it seems to be proved by this analysis too,
    the data are reliable only for the first 10 years approximately.

"""
import argparse
import json
import os
import warnings

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from scipy.ndimage import minimum_filter
from scipy.optimize import least_squares, nnls

# ----------------------------------------------------------------------------- data


def load_panel(path, unit_col, year_col, value_col, year_min=None, year_max=None, min_coverage=1.0):
    """CSV (unit, year, log-price) -> matrix P[unit, year] of log-prices."""
    df = pd.read_csv(path, usecols=[unit_col, year_col, value_col])
    if year_min is not None:
        df = df[df[year_col] >= year_min]
    if year_max is not None:
        df = df[df[year_col] <= year_max]
    df = df[np.isfinite(df[value_col])]
    P = df.pivot_table(index=unit_col, columns=year_col, values=value_col, aggfunc="mean")
    P = P.loc[P.notna().mean(axis=1) >= min_coverage]
    if P.empty:
        raise ValueError("No unit satisfies the coverage requirement.")
    years = P.columns.to_numpy()
    if np.any(np.diff(years) != 1):
        raise ValueError("Years are not consecutive: annual lags are not defined.")
    return P.to_numpy(float), years


# ----------------------------------------------------------------------------- empirical estimators


def empirical_variograms(P, max_lag):
    """
    Four estimators of the temporal variogram for lags tau = 1..max_lag (years).
      V1     : <[p(r,t+tau) - p(r,t)]^2> over t and r                       (raw second moment)
      V1_var : mean over r of the variance over t of the increments          ('variance' reading of V1)
      V2     : like V1, after removing the spatial mean pbar(t) at each year (V2 of the paper)
      V2_var : like V2, after also removing the mean drift of each unit
    Diagnostic: common_ms = <(pbar(t+tau) - pbar(t))^2>_t.
    Exact identity (balanced panel):  V1 = V2 + common_ms.
    """
    pbar = np.nanmean(P, axis=0)
    Pt = P - pbar
    rows = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for tau in range(1, max_lag + 1):
            d = P[:, tau:] - P[:, :-tau]
            dt = Pt[:, tau:] - Pt[:, :-tau]
            db = pbar[tau:] - pbar[:-tau]
            rows.append(
                dict(
                    tau=tau,
                    V1=np.nanmean(d**2),
                    V1_var=np.nanmean(np.nanvar(d, axis=1)),
                    V2=np.nanmean(dt**2),
                    V2_var=np.nanmean(np.nanvar(dt, axis=1)),
                    common_ms=np.mean(db**2),
                    n_pairs=int(np.isfinite(d).sum()),
                )
            )
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- model

# Composite Gauss-Legendre quadrature on lambda/lambda_max in [0, 1], geometric panels
_X, _W = np.polynomial.legendre.leggauss(16)
_EDGES = np.concatenate([[0.0], np.geomspace(1e-5, 1.0, 48)])
_A, _B = _EDGES[:-1, None], _EDGES[1:, None]
_U = (0.5 * (_B - _A) * _X + 0.5 * (_B + _A)).ravel()
_WU = (0.5 * (_B - _A) * _W).ravel()


def V_shape(tau, T, S):
    """V(tau)/B = (1/4pi) Int_0^{1/S} [G(lam,0) - G(lam,tau)] dlam   (kappa = 0)."""
    tau = np.atleast_1d(np.asarray(tau, float))[:, None]
    lam, w = _U / S, _WU / S
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


def V_model(tau, B, T, S, c=0.0):
    """Exact theoretical variogram (+ optional nugget c = 2 * sigma_measurement^2)."""
    return B * V_shape(tau, T, S) + c


# ----------------------------------------------------------------------------- fit


def _solve_linear(v, V, B_fixed, nugget):
    """For fixed (T, S), B (and c >= 0) are linear: least squares on RELATIVE residuals (V_mod - V)/V."""
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


def fit_model(tau, V, B_fixed=None, nugget=False, grid=80, n_refine=6):
    """
    Fit over (log T, log S), with B (and c) projected out analytically.

    The cost landscape has MULTIPLE LOCAL MINIMA (even without noise): the optimiser is started
    from the best local minima of a (T, S) grid and the best result is kept.

    Returns (best, landscape):
      best      : dict(B, T, S, c, rms_rel)
      landscape : dict(T, S, cost, minima), with the grid axes, the cost map on that grid and the
                  distinct local minima found as (rms_rel, T, S), best first.
    """
    tau, V = np.asarray(tau, float), np.asarray(V, float)

    def resid(z):
        T, S = np.exp(z)
        return _solve_linear(V_shape(tau, T, S), V, B_fixed, nugget)[2]

    Tg, Sg = np.geomspace(0.2, 40, grid), np.geomspace(0.03, 40, grid)
    cost = np.array([[np.sum(resid(np.log([T, S])) ** 2) for S in Sg] for T in Tg])
    ismin = cost == minimum_filter(cost, size=3, mode="nearest")
    cells = np.argwhere(ismin)[np.argsort(cost[ismin])[:n_refine]]

    cands = []
    for i, j in cells:
        sol = least_squares(resid, np.log([Tg[i], Sg[j]]), bounds=(np.log([0.1, 0.01]), np.log([100, 100])),
                            ftol=1e-13, xtol=1e-13, gtol=1e-13)
        cands.append((sol.cost, sol.x))
    cands.sort(key=lambda c: c[0])

    uniq = []
    for k, z in cands:  # drop duplicates (same minimum reached from several cells)
        if all(np.max(np.abs(z - z2)) > 1e-2 for _, z2 in uniq):
            uniq.append((k, z))

    T, S = np.exp(uniq[0][1])
    B, c, r = _solve_linear(V_shape(tau, T, S), V, B_fixed, nugget)
    best = dict(B=B, T=T, S=S, c=c, rms_rel=float(np.sqrt(np.mean(r**2))))
    minima = [(float(np.sqrt(2 * k / len(V))), *np.exp(z)) for k, z in uniq]
    return best, dict(T=Tg, S=Sg, cost=cost, minima=minima)


# ----------------------------------------------------------------------------- plots


def plot_results(ev, best, landscape, target, path):
    fig = Figure(figsize=(16, 4.6))
    ax = fig.subplots(1, 3)
    t = ev.tau.values
    tt = np.linspace(t.min(), t.max(), 200)

    series = [("V1", "o"), ("V2", "^")]
    if target not in ("V1", "V2"):
        series.append((target, "s"))
    for k, m in series:
        ax[0].loglog(t, ev[k], m + "-", label=k, ms=5)
    ax[0].loglog(tt, V_model(tt, best["B"], best["T"], best["S"], best["c"]), "k--", lw=2,
                 label=f"fit on {target}")
    ax[0].set(xlabel="tau (years)", ylabel="V(tau)", title="Temporal variogram")
    ax[0].legend(fontsize=8)

    ax[1].plot(t, ev.V1 - ev.V2, "o-", label="V1 - V2")
    ax[1].plot(t, ev.common_ms, "x--", label="<(Δ pbar)^2>_t")
    ax[1].set(xlabel="tau (years)", title="Difference between V1 and V2")
    ax[1].legend(fontsize=8)

    im = ax[2].pcolormesh(landscape["S"], landscape["T"], np.log10(landscape["cost"]), shading="auto")
    ax[2].set(xscale="log", yscale="log", xlabel="S (years)", ylabel="T (years)",
              title="log10 cost (B profiled out): degeneracy valley")
    ax[2].plot(best["S"], best["T"], "r*", ms=14)
    fig.colorbar(im, ax=ax[2])

    fig.tight_layout()
    fig.savefig(path, dpi=130)


# ----------------------------------------------------------------------------- main


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", help="input CSV file")
    ap.add_argument("--unit-col", default="codecommune", help="column identifying the spatial unit")
    ap.add_argument("--year-col", default="anno", help="column with the (integer) year")
    ap.add_argument("--value-col", default="log_prezzo", help="column with the log-price")
    ap.add_argument("--target", default="V2", choices=["V1", "V1_var", "V2", "V2_var"],
                    help="empirical variogram to fit")
    ap.add_argument("--max-lag", type=int, default=10, help="largest lag tau (years)")
    ap.add_argument("--year-min", type=int)
    ap.add_argument("--year-max", type=int)
    ap.add_argument("--min-coverage", type=float, default=1.0,
                    help="minimum fraction of years present for a unit to be kept")
    ap.add_argument("--B-fixed", type=float, help="fix B = A^2/D (e.g. 2*pi*spatial slope)")
    ap.add_argument("--nugget", action="store_true", help="add a constant c >= 0 (measurement noise)")
    ap.add_argument("--D", type=float,
                    help="diffusion coefficient in km^2/year; only used to convert S and T into lengths "
                         "(D is not identifiable from the variogram)")
    ap.add_argument("--outdir", default="output_temporal_variogram")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    P, years = load_panel(a.csv, a.unit_col, a.year_col, a.value_col, a.year_min, a.year_max, a.min_coverage)
    print(f"Panel: {P.shape[0]} units x {P.shape[1]} years ({years[0]}-{years[-1]})")

    ev = empirical_variograms(P, a.max_lag)
    ev["common_share_V1ms"] = ev.common_ms / ev.V1
    ev.to_csv(os.path.join(a.outdir, "empirical_variograms.csv"), index=False)
    print("\nEmpirical variograms (log-price^2):")
    print(ev.round(4).to_string(index=False))

    best, landscape = fit_model(ev.tau.values, ev[a.target].values, B_fixed=a.B_fixed, nugget=a.nugget)
    print(f"\nFit on {a.target}  (B_fixed={a.B_fixed}, nugget={a.nugget})")
    print(f"  B = A^2/D = {best['B']:.3f}   (equivalent spatial slope B/2pi = {best['B'] / (2 * np.pi):.3f})")
    print(f"  T = {best['T']:.2f} years      S = a^2/D = {best['S']:.2f} years      c = {best['c']:.4f}")
    print(f"  RMS relative error = {best['rms_rel']:.3f}")
    if a.D is not None:
        print(f"  (for D = {a.D:g} km^2/year: a = sqrt(D S) = {np.sqrt(a.D * best['S']):.1f} km, "
              f"sqrt(D T) = {np.sqrt(a.D * best['T']):.1f} km)")

    print("  local minima found (relative RMS, T, S):")
    for rms, T_, S_ in landscape["minima"]:
        print(f"    {rms:.4f}  T={T_:.2f}  S={S_:.2f}")

    res = dict(target=a.target, **{k: float(best[k]) for k in ["B", "T", "S", "c", "rms_rel"]})
    with open(os.path.join(a.outdir, "fit.json"), "w") as f:
        json.dump(res, f, indent=2)

    plot_results(ev, best, landscape, a.target, os.path.join(a.outdir, "variogram_fit.png"))
    print(f"\nOutput in {a.outdir}/")


if __name__ == "__main__":
    main()