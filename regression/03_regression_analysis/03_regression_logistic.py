#!/usr/bin/env python3
# ==============================================================================
# SCRIPT:  03_regression_logistic.py
# PURPOSE: Logistic regression comparison between Original and Synthetic Data.
#          Parallel across 18 cores via multiprocessing.Pool.
#          Each worker reads its own SD + OD pair — low per-process memory.
# ==============================================================================

import json
import os
import re
import glob
import gc
import warnings
import time
import sys
import multiprocessing as mp
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import statsmodels.api as sm

warnings.filterwarnings("ignore", category=FutureWarning)

N_WORKERS = 18

# ── Terminal Colors & Formatting ──────────────────────────────────────────────
class Colors:
    HEADER  = '\033[95m'
    BLUE    = '\033[94m'
    CYAN    = '\033[96m'
    GREEN   = '\033[92m'
    WARNING = '\033[93m'
    RED     = '\033[91m'
    ENDC    = '\033[0m'
    BOLD    = '\033[1m'
    DIM     = '\033[2m'

def _ts():
    return datetime.now().strftime('%H:%M:%S')

def log_info(msg):    print(f"{Colors.BLUE}[{_ts()}] INFO:{Colors.ENDC} {msg}", flush=True)
def log_success(msg): print(f"{Colors.GREEN}[{_ts()}] \u2713{Colors.ENDC} {msg}", flush=True)
def log_warn(msg):    print(f"{Colors.WARNING}[{_ts()}] WARN:{Colors.ENDC} {msg}", flush=True)
def log_error(msg):   print(f"{Colors.RED}[{_ts()}] ERROR:{Colors.ENDC} {msg}", flush=True)

# ── Memory Monitor (reads /proc – zero external dependencies) ─────────────────
def _read_mem():
    """Return dict with process and system memory stats in GB."""
    info = {"proc_gb": 0, "total_gb": 1, "avail_gb": 1, "used_pct": 0, "swap_gb": 0}
    try:
        with open("/proc/self/status") as f:
            for ln in f:
                if ln.startswith("VmRSS:"):
                    info["proc_gb"] = int(ln.split()[1]) / 1_048_576
    except OSError:
        pass
    try:
        with open("/proc/meminfo") as f:
            kv = {}
            for ln in f:
                parts = ln.split()
                kv[parts[0].rstrip(":")] = int(parts[1])
        info["total_gb"] = kv["MemTotal"] / 1_048_576
        info["avail_gb"] = kv["MemAvailable"] / 1_048_576
        info["used_pct"] = 100.0 * (1 - kv["MemAvailable"] / kv["MemTotal"])
        st, sf = kv.get("SwapTotal", 0), kv.get("SwapFree", 0)
        info["swap_gb"] = (st - sf) / 1_048_576
    except OSError:
        pass
    return info

def _mem_bar():
    """One-line coloured RAM + Swap status bar."""
    m = _read_mem()
    bw = 15
    filled = int(m["used_pct"] / 100 * bw)
    bar = "\u2588" * filled + "\u2591" * (bw - filled)
    clr = (Colors.GREEN if m["used_pct"] < 60
           else Colors.WARNING if m["used_pct"] < 85
           else Colors.RED)
    return (f"{clr}RAM [{bar}] {m['used_pct']:4.1f}%{Colors.ENDC}  "
            f"Proc {m['proc_gb']:.2f} GB  "
            f"Free {m['avail_gb']:.1f}/{m['total_gb']:.1f} GB  "
            f"Swap {m['swap_gb']:.1f} GB")

# ── Live Progress (overwrites 2 terminal lines in-place) ─────────────────────
_progress_started = False

def _show_progress(idx, total, filename, n_models, t0):
    global _progress_started
    elapsed = time.time() - t0
    pct  = (idx + 1) / total * 100 if total else 0
    rate = n_models / elapsed if elapsed > 1e-3 else 0.0
    if idx > 0 and elapsed > 0:
        eta = str(timedelta(seconds=int(elapsed / (idx + 1) * (total - idx - 1))))
    else:
        eta = "..."

    bw = 25
    filled = int(pct / 100 * bw)
    bar = "\u2501" * filled + ("\u2578" if filled < bw else "") + "\u2500" * max(0, bw - filled - 1)

    short = os.path.basename(filename)
    if len(short) > 52:
        short = short[:49] + "..."

    line1 = (f"  {Colors.CYAN}[{bar}]{Colors.ENDC} "
             f"{idx+1}/{total} ({pct:4.1f}%)  "
             f"{n_models:,} models  {rate:,.0f} m/s  "
             f"ETA {eta}  ({elapsed:.0f}s)")
    line2 = f"  {_mem_bar()}  {Colors.DIM}\u00bb {short}{Colors.ENDC}"

    if _progress_started:
        sys.stdout.write("\033[2A")        # move cursor up 2 lines
    _progress_started = True
    sys.stdout.write(f"\033[2K{line1}\n\033[2K{line2}\n")
    sys.stdout.flush()

# ══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*78}{Colors.ENDC}")
print(f"{Colors.HEADER}{Colors.BOLD}"
      f"   PARALLEL LOGISTIC REGRESSION PIPELINE  ({N_WORKERS} cores)                   "
      f"{Colors.ENDC}")
print(f"{Colors.HEADER}{Colors.BOLD}{'='*78}{Colors.ENDC}")
print(f"{Colors.DIM}NumPy {np.__version__}  \u2022  pandas {pd.__version__}  \u2022  "
      f"statsmodels {sm.__version__}{Colors.ENDC}")
print(f"  {_mem_bar()}\n")

# ── Configuration & Paths ─────────────────────────────────────────────────────
config_path = os.path.join("..", "config", "config.json")
with open(config_path) as f:
    config = json.load(f)

MAX_P      = 10
FEAT_NAMES = [f"X{i}" for i in range(1, MAX_P + 1)]

od_dir     = os.path.join("..", "data", "original")
sd_dir     = os.path.join("..", "data", "synthetic")
result_dir = os.path.join("..", "results")
os.makedirs(result_dir, exist_ok=True)

_OD_RE = re.compile(
    r"OD_logistic_(binary|continuous)_N(\d+)_p(\d+)_rho([\d.]+)_sig([\d.]+|NA)_p1[\w.]+")
_SD_RE = re.compile(
    r"SD_(\w+?)_logistic_(binary|continuous)_N(\d+)_p(\d+)_rho([\d.]+)_sig([\d.]+|NA)_p1[\w.]+")

# ── Helpers ───────────────────────────────────────────────────────────────────
def _pad(arr, length=MAX_P + 1):
    out = np.full(length, np.nan)
    out[:len(arr)] = arr
    return out

def _fit_logit(X, y):
    """Fit Logit and return (coefs, std_errs, pvalues, pseudo_r2, accuracy, status)."""
    Xc = sm.add_constant(X, has_constant="add")
    try:
        import warnings as _warnings
        with _warnings.catch_warnings(record=True) as _w:
            _warnings.simplefilter("always")
            r = sm.Logit(y, Xc).fit(disp=0, maxiter=200)
        warned_separation = any(
            issubclass(w.category, sm.tools.sm_exceptions.PerfectSeparationWarning)
            for w in _w
        )
        warned_hessian = any(
            issubclass(w.category, sm.tools.sm_exceptions.HessianInversionWarning)
            for w in _w
        )
        bse_nan = np.all(np.isnan(r.bse))
        if warned_separation or bse_nan:
            k = Xc.shape[1]
            nan_k = np.full(k, np.nan)
            return nan_k, nan_k.copy(), nan_k.copy(), np.nan, np.nan, "separation"
        if warned_hessian:
            k = Xc.shape[1]
            nan_k = np.full(k, np.nan)
            return nan_k, nan_k.copy(), nan_k.copy(), np.nan, np.nan, "hessian"
        if np.any(np.abs(r.params) > 20):
            k = Xc.shape[1]
            nan_k = np.full(k, np.nan)
            return nan_k, nan_k.copy(), nan_k.copy(), np.nan, np.nan, "unstable"
        accuracy = np.mean((r.predict(Xc) >= 0.5).astype(int) == y.astype(int))
        return r.params, r.bse, r.pvalues, r.prsquared, accuracy, "ok"
    except sm.tools.sm_exceptions.PerfectSeparationError:
        k = Xc.shape[1]
        nan_k = np.full(k, np.nan)
        return nan_k, nan_k.copy(), nan_k.copy(), np.nan, np.nan, "separation"
    except Exception:
        k = Xc.shape[1]
        nan_k = np.full(k, np.nan)
        return nan_k, nan_k.copy(), nan_k.copy(), np.nan, np.nan, "failed"

# ── Index OD files by (var_type, N, p, rho, sig_str) ─────────────────────────
log_info("Scanning data directories \u2026")
od_files = sorted(glob.glob(os.path.join(od_dir, "OD_logistic_*.parquet")))
sd_files = sorted(glob.glob(os.path.join(sd_dir, "SD_*_logistic_*.parquet")))

if not od_files:
    log_error(f"No OD parquet files in {od_dir}"); sys.exit(1)
if not sd_files:
    log_error(f"No SD parquet files in {sd_dir}"); sys.exit(1)

od_map = {}
for fp in od_files:
    m = _OD_RE.search(os.path.basename(fp))
    if m:
        od_map[(m.group(1), m.group(2), m.group(3),
                m.group(4), m.group(5))] = fp

log_success(f"Found {len(od_files)} OD + {len(sd_files)} SD files  "
            f"({len(od_map)} unique OD scenarios)")

# ── Parse SD files ────────────────────────────────────────────────────────────
sd_tasks = []
skipped_files = 0
for fp in sd_files:
    m = _SD_RE.search(os.path.basename(fp))
    if not m:
        log_warn(f"Regex skip: {os.path.basename(fp)}")
        skipped_files += 1
        continue
    od_key = (m.group(2), m.group(3), m.group(4), m.group(5), m.group(6))
    if od_key not in od_map:
        skipped_files += 1
        continue
    sd_tasks.append({
        "fp":       fp,
        "method":   m.group(1),
        "var_type": m.group(2),
        "N":        int(m.group(3)),
        "p":        int(m.group(4)),
        "rho":      float(m.group(5)),
        "sig_str":  m.group(6),
        "od_key":   od_key,
        "od_fp":    od_map[od_key],
    })

total_tasks = len(sd_tasks)
log_info(f"{total_tasks} SD files to process ({skipped_files} skipped)  |  "
         f"{N_WORKERS} worker processes")
if skipped_files:
    log_warn(f"Skipped {skipped_files} SD file(s) (missing OD match or bad filename)")

# ── Parallel Worker Function ──────────────────────────────────────────────────
def _process_one_sd(task):
    """Standalone worker: reads OD + SD, runs OLS per iteration, returns list[dict]."""
    import numpy as _np
    import pandas as _pd
    import statsmodels.api as _sm

    MAX_P_W = 10

    def _pad_w(arr, length=MAX_P_W + 1):
        out = _np.full(length, _np.nan)
        out[:len(arr)] = arr
        return out

    def _fit_w(X, y):
        Xc = _sm.add_constant(X, has_constant="add")
        try:
            import warnings as _warnings
            with _warnings.catch_warnings(record=True) as _w:
                _warnings.simplefilter("always")
                r = _sm.Logit(y, Xc).fit(disp=0, maxiter=200)
            # PerfectSeparation is issued as a Warning (not Error) in most statsmodels versions
            warned_separation = any(
                issubclass(w.category, _sm.tools.sm_exceptions.PerfectSeparationWarning)
                for w in _w
            )
            # HessianInversion means SEs are unreliable — coefficients may be extreme
            warned_hessian = any(
                issubclass(w.category, _sm.tools.sm_exceptions.HessianInversionWarning)
                for w in _w
            )
            # Hessian inversion failure means bse/pvalues are NaN even if fit "converged"
            bse_nan = _np.all(_np.isnan(r.bse))
            if warned_separation or bse_nan:
                k = Xc.shape[1]
                nan_k = _np.full(k, _np.nan)
                return nan_k, nan_k.copy(), nan_k.copy(), _np.nan, _np.nan, "separation"
            if warned_hessian:
                k = Xc.shape[1]
                nan_k = _np.full(k, _np.nan)
                return nan_k, nan_k.copy(), nan_k.copy(), _np.nan, _np.nan, "hessian"
            # Coefficients with extreme magnitude indicate near-separation that
            # didn't trigger an explicit warning — treat as numerically unstable.
            if _np.any(_np.abs(r.params) > 20):
                k = Xc.shape[1]
                nan_k = _np.full(k, _np.nan)
                return nan_k, nan_k.copy(), nan_k.copy(), _np.nan, _np.nan, "unstable"
            accuracy = _np.mean((r.predict(Xc) >= 0.5).astype(int) == y.astype(int))
            return r.params, r.bse, r.pvalues, r.prsquared, accuracy, "ok"
        except _sm.tools.sm_exceptions.PerfectSeparationError:
            k = Xc.shape[1]
            nan_k = _np.full(k, _np.nan)
            return nan_k, nan_k.copy(), nan_k.copy(), _np.nan, _np.nan, "separation"
        except Exception:
            k = Xc.shape[1]
            nan_k = _np.full(k, _np.nan)
            return nan_k, nan_k.copy(), nan_k.copy(), _np.nan, _np.nan, "failed"

    fp       = task["fp"]
    od_fp    = task["od_fp"]
    method   = task["method"]
    var_type = task["var_type"]
    N        = task["N"]
    p_val    = task["p"]
    rho      = task["rho"]
    sig_str  = task["sig_str"]
    sigma_out = _np.nan if sig_str == "NA" else float(sig_str)
    x_cols   = [f"X{i}" for i in range(1, p_val + 1)]

    od_df = _pd.read_parquet(od_fp)
    sd_df = _pd.read_parquet(fp)

    od_groups = {k: v for k, v in od_df.groupby("iter")}
    sd_groups = {k: v for k, v in sd_df.groupby("iter")}
    del od_df, sd_df

    rows = []
    for it, sd_grp in sd_groups.items():
        od_grp = od_groups.get(it)
        if od_grp is None:
            continue

        X_od = od_grp[x_cols].to_numpy(dtype=_np.float64)
        y_od = od_grp["y"].to_numpy(dtype=_np.float64)
        X_sd = sd_grp[x_cols].to_numpy(dtype=_np.float64)
        y_sd = sd_grp["y"].to_numpy(dtype=_np.float64)

        b_od, se_od, pv_od, r2_od, acc_od, status_od = _fit_w(X_od, y_od)
        b_sd, se_sd, pv_sd, r2_sd, acc_sd, status_sd = _fit_w(X_sd, y_sd)

        row = {
            "method": method, "var_type": var_type,
            "N": N, "p": p_val, "rho": rho,
            "sigma_2": sigma_out, "iter": int(it),
            "fit_status_od": status_od, "fit_status_sd": status_sd,
            "pseudo_r2_od": r2_od, "pseudo_r2_sd": r2_sd,
            "accuracy_od": acc_od, "accuracy_sd": acc_sd,
        }
        for i, v in enumerate(_pad_w(b_od)):  row[f"beta_od_{i}"] = v
        for i, v in enumerate(_pad_w(b_sd)):  row[f"beta_sd_{i}"] = v
        for i, v in enumerate(_pad_w(se_od)): row[f"se_od_{i}"]   = v
        for i, v in enumerate(_pad_w(se_sd)): row[f"se_sd_{i}"]   = v
        for i, v in enumerate(_pad_w(pv_od)): row[f"pval_od_{i}"] = v
        for i, v in enumerate(_pad_w(pv_sd)): row[f"pval_sd_{i}"] = v
        rows.append(row)
    return rows

# ── Dispatch with multiprocessing.Pool ────────────────────────────────────────
t_global      = time.time()
total_models  = 0
all_results   = []

log_info(f"Dispatching {total_tasks} SD files across {N_WORKERS} cores …")
print()   # blank lines for progress display
print()

with mp.Pool(processes=N_WORKERS) as pool:
    for task_idx, rows in enumerate(pool.imap_unordered(_process_one_sd, sd_tasks)):
        all_results.extend(rows)
        total_models += len(rows) * 2
        _show_progress(task_idx, total_tasks,
                       sd_tasks[task_idx]["fp"] if task_idx < len(sd_tasks) else "",
                       total_models, t_global)

print()

elapsed = time.time() - t_global
log_success(f"Completed {total_models:,} Logit fits in {elapsed:.1f}s "
            f"({total_models / max(elapsed, 1e-9):,.0f} models/sec)")

# ── Save Results ──────────────────────────────────────────────────────────────
log_info("Saving results \u2026")
df = pd.DataFrame(all_results)
del all_results
gc.collect()

nan_od = df["pseudo_r2_od"].isna().sum()
nan_sd = df["pseudo_r2_sd"].isna().sum()
if nan_od or nan_sd:
    log_warn(f"Failed fits (degenerate data): OD={nan_od}, SD={nan_sd}")

parquet_path = os.path.join(result_dir, "aggregated_logistic_metrics.parquet")
df.to_parquet(parquet_path, index=False, engine="pyarrow")

log_success(f"Saved \u2192 {parquet_path} ({os.path.getsize(parquet_path)/1e6:.1f} MB)")
print(f"\n  {_mem_bar()}")
print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*78}{Colors.ENDC}\n")