#!/usr/bin/env python3
"""
Diffusion of house-price increases around TGV stations.

Replication of the analysis in SM-4 and Fig. 4 of arXiv:2412.14624.

The paper declares neither t0, nor the study area, nor the fit window. Here they are explicit
parameters, definable SEPARATELY FOR EVERY EVENT (optional keys in the CITIES dictionary, or
--event on the command line):
  * t0        : reference year. Defaults in CITIES: Lyon 1981, Tours 1990, Bordeaux_1 1990, Bordeaux_2 2017,
                Rennes 2017. Set as the opening year of TGV's in every city.
  * L, shape  : area = square of half-side L ("square") or disc of radius L ("disc"), centred on the city
  * fit_years : the fit uses the years t0+1 ... t0+fit_years (truncated to the last year of the dataset)
Precedence:  --event  >  CITIES  >  global values (--half-side, --shape, --fit-years).
R(t0) = 0 identically, so sigma^2(t0) is undefined and the series starts at t0+1.


#########################################################################################################

METHOD. For every commune r, R(r, t) = p(r, t) - p(r, t0)  (p = log-price). Summing over all
communes in the study area,

    sigma^2(t) = sum_r R(r, t) |r - c|^2 / sum_r R(r, t)

where c is the FIXED city centre.
Then sigma^2(t) ~ b * x and the diffusion coefficient is D = b / 4.

########################################################################################################

INPUT: data.csv, data_no_inflation.csv. The input is given by the command line.
Example:
    python tgv_diffusion.py data.csv --outdir output_tgv

OUTPUT: the code returns a directory called "output_tgv", with the chart of sigma2(t) for every city.
It also returns the map of the studied area for every city at t0 and t0+fit_years.


########################################################################################################


FIT: straight line CONSTRAINED THROUGH THE ORIGIN,  y = b x  (least squares: b = sum(x y) / sum(x^2)).
  --origin first (default): the origin is the first point of the series (t0+1):
                            x = t - (t0+1),  y = sigma^2(t) - sigma^2(t0+1)
  --origin t0             : the origin is the physical one, sigma^2(t0) = 0 (point source):
                            x = t - t0,      y = sigma^2(t)
  R2 = 1 - SSres/SStot with SStot taken around the mean of y; with the constraint it can be negative.

--mode paper   : raw R = p(t) - p(t0) (nominal prices, as in SM-4)
--mode excess  : R = max(0, p(t) - p(t0))  (Mantains only the price increases)


#######################################################################################################

USAGE:
    python tgv_diffusion.py data.csv --outdir output_tgv --mode excess
    python tgv_diffusion.py data_no_inflation.csv --event Tours:t0=1990,L=50,fit_years=8 \\
                                        --event Bordeaux_1:L=90,shape=square,fit_years=5

    If you wanna set the studied area and fit years as the same for every city, you can go to the CLI configuration (located at
    the beginning of main) in the code and change the given default values written there.
    You can easily change the starting year t0 in the dictionary "CITIES", just below here.



"""
import argparse
import os
from itertools import cycle

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SHAPES = ("square", "disc")
MODES = ("paper", "excess")
SENSITIVITY_L = (30, 40, 50, 75, 100)   # km; area grid of the sensitivity table (event's own L is added)

# One entry per event. The part of the key before the first "_" is the city name (used to look the
# centre up by name); the suffix tells apart several events on the same city (Bordeaux_1, Bordeaux_2).
# Optional keys (otherwise the global command-line values apply):
#   L=50, shape="square", fit_years=8 (see EVENT_KEYS),  alt=[...] extra t0 values for the sensitivity table.
CITIES = {
    "Lyon":       dict(code="69123", lonlat=(4.8357, 45.7640),  t0=1981, alt=[1976]),#alt is the extra t0 value for the study of sensitivity
    "Tours":      dict(code="37261", lonlat=(0.6848, 47.3941),  t0=1990, alt=[1985]),
    "Bordeaux_1": dict(code="33063", lonlat=(-0.5792, 44.8378), t0=1990, alt=[1985]),
    "Bordeaux_2": dict(code="33063", lonlat=(-0.5792, 44.8378), t0=2017, alt=[2012]),
    "Rennes":     dict(code="35238", lonlat=(-1.6794, 48.1147), t0=2017, alt=[2012]),
}

# Per-event parameters that can be overridden from the command line (--event), with their types
EVENT_KEYS = dict(t0=int, L=float, shape=str, fit_years=int)

# Keys of analyze() results that are arrays (saved in the sigma^2 series file, not in the summary table)
SERIES_KEYS = ("tau", "sigma2", "sumR")

# Plot style (colour, marker) per event; other events get one of the extra styles in turn
EVENT_STYLES = {
    "Lyon": ("b", "^"),
    "Tours": ("g", "s"),
    "Bordeaux_1": ("r", "o"),
    "Bordeaux_2": ("orange", "p"),
    "Rennes": ("m", "d"),
}
EXTRA_STYLES = [("c", "v"), ("k", "*"), ("brown", "X"), ("olive", "h")]


def city_of(key):
    """City name of an event key: 'Bordeaux_2' -> 'Bordeaux'."""
    return key.split("_", 1)[0]


# ----------------------------------------------------------------------------- geometry
def lambert93(lon, lat):
    """lon/lat (degrees, GRS80) -> Lambert-93 (m). Only used as a fallback for the city centre."""
    a, e = 6378137.0, 0.0818191910428158
    p0, l0, p1, p2 = np.radians([46.5, 3.0, 44.0, 49.0])
    lam, phi = np.radians(lon), np.radians(lat)

    def m(p):
        return np.cos(p) / np.sqrt(1 - (e * np.sin(p)) ** 2)

    def t(p):
        return np.tan(np.pi / 4 - p / 2) / ((1 - e * np.sin(p)) / (1 + e * np.sin(p))) ** (e / 2)

    n = (np.log(m(p1)) - np.log(m(p2))) / (np.log(t(p1)) - np.log(t(p2)))
    F = m(p1) / (n * t(p1) ** n)
    rho, rho0 = a * F * t(phi) ** n, a * F * t(p0) ** n
    return 700000 + rho * np.sin(n * (lam - l0)), 6600000 + rho0 - rho * np.cos(n * (lam - l0))


def city_center(name, cfg, meta):
    """City centre in km (Lambert-93): from the dataset (INSEE code, then name), else from lon/lat."""
    if cfg["code"] in meta.index:
        r = meta.loc[cfg["code"]]
        return np.array([r.x, r.y]) / 1e3, f"dataset, code {cfg['code']}"
    m = meta[meta.nom.str.lower() == name.lower()]
    if len(m):
        return np.array([m.iloc[0].x, m.iloc[0].y]) / 1e3, f"dataset, name '{name}'"
    x, y = lambert93(*cfg["lonlat"])
    return np.array([x, y]) / 1e3, "fallback: Lambert-93 from lon/lat"


def window_mask(xy, c, L, shape):
    """Boolean mask of the points (km) inside the study area centred on c."""
    d = xy - c
    return np.abs(d).max(1) <= L if shape == "square" else np.hypot(d[:, 0], d[:, 1]) <= L


# ----------------------------------------------------------------------------- data
def load(path):
    """Read the dataset. Returns P (communes x years DataFrame of log-prices) and meta (name, x, y)."""
    df = pd.read_csv(path, dtype={"codecommune": str},
                     usecols=["codecommune", "nom", "anno", "log_prezzo", "x", "y"])
    P = df.pivot_table(index="codecommune", columns="anno", values="log_prezzo", aggfunc="mean")
    meta = df.drop_duplicates("codecommune").set_index("codecommune").loc[P.index, ["nom", "x", "y"]]
    return P, meta


# ----------------------------------------------------------------------------- per-event parameters
def parse_event_overrides(specs, cfgs):
    """Apply overrides such as 'Tours:t0=1990,L=50,fit_years=8' (one per --event)."""
    for spec in specs:
        key, _, rest = spec.partition(":")
        if key not in cfgs:
            raise ValueError(f"--event: unknown event '{key}' (available: {', '.join(cfgs)})")
        for item in filter(None, rest.split(",")):
            k, _, v = item.partition("=")
            if k not in EVENT_KEYS or not v:
                raise ValueError(f"--event {spec}: invalid parameter '{item}' (allowed: {', '.join(EVENT_KEYS)})")
            cfgs[key][k] = EVENT_KEYS[k](v)


def fit_window(years, t0, fit_years):
    """(i0, i1, tmax): column indices of t0 and of the last fit year (t0+fit_years, truncated to the dataset)."""
    if t0 not in years:
        raise ValueError(f"t0={t0} is not among the dataset years ({years.min()}-{years.max()})")
    tmax = min(t0 + fit_years, int(years.max()))
    if tmax < t0 + 2:
        raise ValueError(f"t0={t0}, fit_years={fit_years}: at least 2 years after t0 are needed "
                         f"(last dataset year: {years.max()})")
    return int(np.where(years == t0)[0][0]), int(np.where(years == tmax)[0][0]), tmax


def check_event(key, cfg, years):
    """Validate the parameters of an event and store its last fit year in cfg['t_end']."""
    if cfg["shape"] not in SHAPES:
        raise ValueError(f"{key}: invalid shape '{cfg['shape']}' (allowed: {', '.join(SHAPES)})")
    if cfg["L"] <= 0:
        raise ValueError(f"{key}: L must be > 0")
    try:
        cfg["t_end"] = fit_window(years, cfg["t0"], cfg["fit_years"])[2]
    except ValueError as e:
        raise ValueError(f"{key}: {e}") from None


# ----------------------------------------------------------------------------- sigma^2 and diffusion
def sigma2(R, xy, c):
    """sigma^2(t) about the fixed centre c, and the total weight sum_r R(r, t) for every t."""
    sw = R.sum(0)
    d2 = ((xy - c) ** 2).sum(1)
    with np.errstate(divide="ignore", invalid="ignore"):     # zero total weight -> nan/inf, handled by callers
        return (R.T @ d2) / sw, sw


def origin_xy(tau, s2, origin):
    """(x, y) coordinates of the through-the-origin fit.
    first : the first point of the series (t0+1) is the origin  ->  x = tau - tau[0],  y = s2 - s2[0]
    t0    : physical origin sigma^2(t0) = 0                     ->  x = tau,           y = s2
    """
    if origin == "first":
        return (tau - tau[0]).astype(float), s2 - s2[0]
    return tau.astype(float), np.asarray(s2, float)


def slope_fit(tau, s2, origin="first"):
    """Line y = b x through the origin: b = sum(x y) / sum(x^2) over the points with x > 0.
    Returns (b, R2), with R2 = 1 - SSres/SStot (SStot around the mean of y)."""
    x, y = origin_xy(tau, s2, origin)
    m = (x > 0) & np.isfinite(y)
    if not m.any():
        return np.nan, np.nan
    x, y = x[m], y[m]
    b = (x @ y) / (x @ x)
    ss = np.sum((y - y.mean()) ** 2)
    return b, (1 - np.sum((y - b * x) ** 2) / ss) if (ss > 0 and len(y) > 1) else np.nan


def analyze(P, years, xy_all, c, t0, L, shape, mode, fit_years, origin, n_boot, rng):
    """sigma^2(t) and its linear fit for one event.
    P: (communes x years) matrix of the log-prices of ALL communes; xy_all: their coordinates in km."""
    i0, i1, tmax = fit_window(years, t0, fit_years)
    dP_all = P[:, i0 + 1: i1 + 1] - P[:, [i0]]              # p(t) - p(t0), t = t0+1 ... tmax
    nat = np.nanmean(dP_all, axis=0)                         # national mean change
    sel = window_mask(xy_all, c, L, shape)
    dP = dP_all[sel]
    ok = np.isfinite(dP).all(1)                              # communes with no missing year in the window
    dP, xy = dP[ok], xy_all[sel][ok]
    if len(xy) < 3:
        raise ValueError(f"only {len(xy)} valid communes in the area (L={L:g}, {shape})")
    
    R = np.clip(dP, 0, None) if mode == "excess" else dP
    #R = np.clip(dP - nat, 0, None) if mode == "excess" else dP

    tau = np.arange(1, R.shape[1] + 1)
    s2, sw = sigma2(R, xy, c)
    b, r2 = slope_fit(tau, s2, origin)

    ci_lo = ci_hi = np.nan                                   # 90% bootstrap interval over communes
    if n_boot:
        bs = []
        for _ in range(n_boot):
            ix = rng.integers(0, len(xy), len(xy))
            s2b, swb = sigma2(R[ix], xy[ix], c)
            if np.all(swb > 0):
                bs.append(slope_fit(tau, s2b, origin)[0])
        if bs:
            ci_lo, ci_hi = np.nanquantile(bs, [0.05, 0.95])

    return dict(t0=t0, t_end=tmax, fit_years=tmax - t0, L=L, shape=shape, mode=mode, origin=origin,
                n_communes=len(xy), n_points=len(tau), slope=b, R2=r2, D_model=b / 4,
                neg_share=float(np.mean(dP < 0)), sumR_min=float(sw.min()), ci_lo=ci_lo, ci_hi=ci_hi,
                tau=tau, sigma2=s2, sumR=sw)


# ----------------------------------------------------------------------------- figures
def sigmoid_field(P):
    """As in Fig. 1: log-price standardized with the national mean and std of that year, then a sigmoid."""
    z = (P - np.nanmean(P, axis=0)) / np.nanstd(P, axis=0)
    return 1 / (1 + np.exp(-z))


def map_title(name, t, t0):
    """Map title: calendar year + position relative to t0 (follows t0 and the fit window automatically)."""
    return f"{name} - {t} (t0)" if t == t0 else f"{name} - {t} (t0+{t - t0})"


def draw_map(ax, xy_all, vals, c, L, shape, title):
    """Scatter map of the price field around c; the study area is outlined with a dashed line."""
    inside = window_mask(xy_all, c, L, shape)
    pad = 1.12 * L
    box = (np.abs(xy_all - c).max(1) <= pad) & np.isfinite(vals)
    out = box & ~inside
    ax.scatter(xy_all[out, 0], xy_all[out, 1], s=7, c="0.85", lw=0)
    sc = ax.scatter(xy_all[inside & box, 0], xy_all[inside & box, 1], s=13, c=vals[inside & box],
                    cmap="plasma", vmin=0.25, vmax=0.75, lw=0)
    if shape == "square":
        ax.add_patch(mpatches.Rectangle(c - L, 2 * L, 2 * L, fill=False, ls="--", lw=1.0, ec="k"))
    else:
        ax.add_patch(mpatches.Circle(c, L, fill=False, ls="--", lw=1.0, ec="k"))
    ax.plot(*c, "kx", ms=12, mew=2.5)
    bar = 50 if 2 * L >= 80 else 20
    x0, y0 = c[0] - pad + 0.04 * L, c[1] + pad - 0.10 * L
    ax.plot([x0, x0 + bar], [y0, y0], "k-", lw=2.5)
    ax.text(x0 + bar / 2, y0 + 0.02 * L, f"{bar} km", ha="center", va="bottom", fontsize=8,
            bbox=dict(fc="white", ec="none", pad=1.0))
    ax.set(xlim=(c[0] - pad, c[0] + pad), ylim=(c[1] - pad, c[1] + pad), aspect="equal", title=title)
    ax.set_xticks([])
    ax.set_yticks([])
    return sc


def map_figure(path, name, xy_all, S, years, c, L, shape, t0, t_after):
    """Before/after price-field maps of one event (years t0 and t_after)."""
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 5.2))
    for ax, t in zip(axs, (t0, t_after)):
        sc = draw_map(ax, xy_all, S[:, int(np.where(years == t)[0][0])], c, L, shape, map_title(name, t, t0))
    fig.colorbar(sc, ax=axs, fraction=0.03, pad=0.02, ticks=[0.25, 0.5, 0.75])
    fig.suptitle(rf"{name}: price field (sigmoid of standardized log-prices) - dashed line = area used for the "
                 rf"$\sigma^2$ calculation ({shape}, L = {L:g} km)", fontsize=9)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_sigma2(path, results, mode, origin):
    """Response function sigma^2(t) of all events, with the linear fits through the origin."""
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    extra = cycle(EXTRA_STYLES)
    for name, r in results.items():
        col, marker = EVENT_STYLES.get(name) or next(extra)
        x, y = origin_xy(r["tau"], r["sigma2"], r["origin"])
        ax.plot(x, y, marker, color=col,
                label=f"{name}: t0 = {r['t0']}, fit {r['t0'] + 1}-{r['t_end']}\n"
                      f"{r['shape']} L = {r['L']:g} km, slope = {r['slope']:.0f} km$^2$/yr")
        xx = np.array([0, x.max() + 0.5])
        ax.plot(xx, r["slope"] * xx, "--", color=col, lw=1.2)
    ax.axhline(0, color="0.8", lw=0.6, zorder=0)
    ax.axvline(0, color="0.8", lw=0.6, zorder=0)
    if origin == "first":
        ax.set(xlabel=r"$t-(t_0+1)$  [years]", ylabel=r"$\sigma^2(t)-\sigma^2(t_0+1)$  [km$^2$]")
    else:
        ax.set(xlabel=r"$t-t_0$  [years]", ylabel=r"$\sigma^2(t)$  [km$^2$]")
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.17), frameon=False)
    ax.set_title(f"Response function ({mode} weights), linear fit through the origin", fontsize=9)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", help="input dataset (see above for the required columns)")
    ap.add_argument("--outdir", default="output", help="output directory (default: %(default)s)")
    ap.add_argument("--half-side", type=float, default=70.0,
                    help="half-side (or radius) of the area in km: global value, can be overridden per event")
    ap.add_argument("--shape", choices=SHAPES, default="disc",
                    help="area shape: global value, can be overridden per event")
    ap.add_argument("--mode", choices=MODES, default="paper", help="weights used for the figures and maps")
    ap.add_argument("--fit-years", type=int, default=7,
                    help="years after t0 used in the fit: global value, can be overridden per event")
    ap.add_argument("--origin", choices=["first", "t0"], default="first",
                    help="origin of the fitted line: 'first' = first point of the series (t0+1); "
                         "'t0' = physical origin sigma^2(t0)=0")
    ap.add_argument("--map-after-years", type=int, default=None,
                    help="years after t0 for the 'after' map (default: last fit year of each event)")
    ap.add_argument("--event", action="append", default=[], metavar="NAME:key=val,...",
                    help="per-event override (repeatable). Keys: t0, L, shape, fit_years. "
                         "E.g.: --event Tours:t0=1990,L=50,fit_years=8 --event Bordeaux_1:shape=square")
    ap.add_argument("--boot", type=int, default=0, help="bootstrap replicas over communes (0 = no confidence interval)")
    ap.add_argument("--seed", type=int, default=0, help="random seed of the bootstrap")
    args = ap.parse_args()

    plt.switch_backend("Agg")
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # event parameters: global < CITIES < --event
    glob = dict(L=args.half_side, shape=args.shape, fit_years=args.fit_years)
    cfgs = {k: {**glob, **v} for k, v in CITIES.items()}

    Pdf, meta = load(args.csv)
    P, years = Pdf.to_numpy(float), Pdf.columns.to_numpy()
    xy_all = meta[["x", "y"]].to_numpy(float) / 1e3
    print(f"Dataset: {P.shape[0]} communes, years {years.min()}-{years.max()}")
    if not np.all(np.diff(years) == 1):
        ap.error("the dataset years are not consecutive: the analysis assumes one column per year")

    try:
        parse_event_overrides(args.event, cfgs)
        for key, cfg in cfgs.items():
            check_event(key, cfg, years)
    except ValueError as e:
        ap.error(str(e))

    rows, res_main = [], {}
    for key, cfg in cfgs.items():
        city = city_of(key)
        c, src = city_center(city, cfg, meta)
        cfg["c"] = c
        print(f"\n== {key}: centre ({c[0]:.1f}, {c[1]:.1f}) km  [{src}]")
        print(f"   t0={cfg['t0']}, fit {cfg['t0'] + 1}-{cfg['t_end']}, area {cfg['shape']} L={cfg['L']:g} km, "
              f"fit origin: {args.origin}")

        for mode in MODES:
            try:
                r = analyze(P, years, xy_all, c, cfg["t0"], cfg["L"], cfg["shape"], mode, cfg["fit_years"],
                            args.origin, args.boot, rng)
            except ValueError as e:
                ap.error(f"{key}: {e}")
            ci = f" [90% CI: {r['ci_lo']:.0f}, {r['ci_hi']:.0f}]" if args.boot else ""
            print(f"  t0={cfg['t0']} {mode:6s}: N={r['n_communes']:4d}  d(sigma^2)/dt = {r['slope']:7.1f} km^2/yr"
                  f"{ci}  D=slope/4 = {r['D_model']:.1f}  R2={r['R2']:.2f}  R<0: {100 * r['neg_share']:.0f}%")
            if r["neg_share"] > 0.05 and mode == "paper":
                print("     WARNING: many negative R weights -> sigma^2 hard to interpret (use --mode excess).")
            rows.append(dict(city=city, event=key, **{k: v for k, v in r.items() if k not in SERIES_KEYS}))
            if mode == args.mode:
                res_main[key] = r        # keyed by event: Bordeaux_1 and Bordeaux_2 share the centre, not t0

    pd.DataFrame(rows).to_csv(os.path.join(args.outdir, "tgv_diffusion.csv"), index=False)

    # ---- sensitivity: area x t0 x weights (slope of sigma^2); shape and fit years are those of the event
    print("\nSensitivity of the slope d(sigma^2)/dt [km^2/yr] to area, t0 and weights:")
    sens = []
    for key, cfg in cfgs.items():
        for t0 in dict.fromkeys([cfg["t0"], *cfg.get("alt", [])]):
            for L in sorted({*SENSITIVITY_L, cfg["L"]}):
                for mode in MODES:
                    try:
                        r = analyze(P, years, xy_all, cfg["c"], t0, L, cfg["shape"], mode, cfg["fit_years"],
                                    args.origin, 0, rng)
                    except ValueError:                       # t0 outside the dataset, empty area, ...
                        continue
                    sens.append(dict(city=city_of(key), event=key, t0=t0, L=L, mode=mode,
                                     slope=r["slope"], R2=r["R2"], N=r["n_communes"]))
    sd = pd.DataFrame(sens)
    sd.to_csv(os.path.join(args.outdir, "tgv_sensitivity.csv"), index=False)
    if len(sd):
        print(sd.set_index(["event", "t0", "mode", "L"])["slope"].unstack("L").round(0).to_string())

    # ---- sigma^2 series (with the coordinates used in the fit and the fitted line)
    def series(key, r):
        x, y = origin_xy(r["tau"], r["sigma2"], r["origin"])
        return pd.DataFrame(dict(city=city_of(key), event=key, t=r["t0"] + r["tau"], tau=r["tau"],
                                 sigma2=r["sigma2"], sumR=r["sumR"], x_fit=x, y_fit_data=y,
                                 y_fit_line=r["slope"] * x))

    pd.concat([series(k, r) for k, r in res_main.items()]).to_csv(
        os.path.join(args.outdir, "tgv_sigma2_series.csv"), index=False)

    # ---- figures
    plot_sigma2(os.path.join(args.outdir, "tgv_sigma2.png"), res_main, args.mode, args.origin)

    S = sigmoid_field(P)
    for key, cfg in cfgs.items():
        # year of the 'after' map: last fit year of the event (or t0 + --map-after-years)
        t_after = cfg["t_end"] if args.map_after_years is None \
            else min(cfg["t0"] + args.map_after_years, int(years.max()))
        map_figure(os.path.join(args.outdir, f"map_{key}_before_after.png"), key, xy_all, S, years,
                   cfg["c"], cfg["L"], cfg["shape"], cfg["t0"], t_after)
    print(f"\nOutput in {args.outdir}/")


if __name__ == "__main__":
    main()