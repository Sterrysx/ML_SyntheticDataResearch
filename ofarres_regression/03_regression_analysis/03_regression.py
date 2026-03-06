#!/usr/bin/env python3
# ==============================================================================
# SCRIPT: 03_regression.py
# PURPOSE: Calculate multiple linear regressions for Original and Synthetic Data.
#          Uses SharedMemory to bypass the GIL and achieve zero-copy parallelism.
# ==============================================================================

import json
import os
import re
import glob
import gc
import warnings
import time
from datetime import datetime
from multiprocessing.shared_memory import SharedMemory

import numpy as np
import pandas as pd
import statsmodels.api as sm
from joblib import Parallel, delayed
import pyarrow as pa
import pyarrow.parquet as pq

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Terminal Colors & Formatting ──────────────────────────────────────────────
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def log_info(msg):    print(f"{Colors.BLUE}[{datetime.now().strftime('%H:%M:%S')}] INFO:{Colors.ENDC} {msg}", flush=True)
def log_success(msg): print(f"{Colors.GREEN}[{datetime.now().strftime('%H:%M:%S')}] SUCCESS:{Colors.ENDC} {msg}", flush=True)
def log_warn(msg):    print(f"{Colors.WARNING}[{datetime.now().strftime('%H:%M:%S')}] WARN:{Colors.ENDC} {msg}", flush=True)

print(f"{Colors.HEADER}{Colors.BOLD}=============================================================================={Colors.ENDC}")
print(f"{Colors.HEADER}{Colors.BOLD}                PARALLEL OLS REGRESSION EVALUATION PIPELINE                   {Colors.ENDC}")
print(f"{Colors.HEADER}{Colors.BOLD}=============================================================================={Colors.ENDC}")
print(f"Environment: NumPy {np.__version__}  •  pandas {pd.__version__}  •  statsmodels {sm.__version__}\n")

# ── Configuration & Paths ─────────────────────────────────────────────────────
config_path = os.path.join("..", "config", "config.json")
with open(config_path) as f:
    config = json.load(f)

MAX_P       = 10
FEAT_COLS   = [f"X{i}" for i in range(1, MAX_P + 1)]
SCHEMA_COLS = FEAT_COLS + ["y"]
OD_META     = ["var_type", "N", "p", "rho", "sigma_2", "iter"]
SD_META     = ["method",   "var_type", "N", "p", "rho", "sigma_2", "iter"]
ARENA_COLS  = MAX_P + 1   

od_dir         = os.path.join("..", "data", "original")
sd_dir         = os.path.join("..", "data", "synthetic")
result_dir     = os.path.join("..", "results")
od_packed_path = os.path.join("..", "data", "OD_packed.parquet")
sd_packed_path = os.path.join("..", "data", "SD_packed.parquet")
os.makedirs(result_dir, exist_ok=True)

_OD_RE = re.compile(r"OD_(binary|continuous)_N(\d+)_p(\d+)_rho([\d.]+)_sig([\d.]+|NA)_p1[\w.]+")
_SD_RE = re.compile(r"SD_(\w+?)_(binary|continuous)_N(\d+)_p(\d+)_rho([\d.]+)_sig([\d.]+|NA)_p1[\w.]+")

# ── File Parsers (Fixed NaN dictionary bug by mapping NA to -1.0) ─────────────
def _read_one_od(fp):
    m = _OD_RE.search(os.path.basename(fp))
    if not m: raise ValueError(f"Regex failed: {os.path.basename(fp)}")
    
    # Map NA to -1.0 so dictionary keys don't break on NaN comparisons
    sig = float(m.group(5)) if m.group(5) != "NA" else -1.0
    
    raw = pd.read_parquet(fp)
    iter_col = raw["iter"].astype(np.int16)
    df = raw.drop(columns=["iter"]).astype("float32")
    for c in SCHEMA_COLS:
        if c not in df.columns: df[c] = float("nan")
            
    df = df[SCHEMA_COLS].copy()
    df["var_type"], df["N"], df["p"] = m.group(1), np.int16(m.group(2)), np.int8(m.group(3))
    df["rho"], df["sigma_2"] = np.float32(m.group(4)), np.float32(sig)
    df["iter"] = iter_col
    return df

def _read_one_sd(fp):
    m = _SD_RE.search(os.path.basename(fp))
    if not m: raise ValueError(f"Regex failed: {os.path.basename(fp)}")
    
    sig = float(m.group(6)) if m.group(6) != "NA" else -1.0
    
    raw = pd.read_parquet(fp)
    iter_col = raw["iter"].astype(np.int16)
    df = raw.drop(columns=["iter"]).astype("float32")
    for c in SCHEMA_COLS:
        if c not in df.columns: df[c] = float("nan")
            
    df = df[SCHEMA_COLS].copy()
    df["method"], df["var_type"] = m.group(1), m.group(2)
    df["N"], df["p"] = np.int16(m.group(3)), np.int8(m.group(4))
    df["rho"], df["sigma_2"] = np.float32(m.group(5)), np.float32(sig)
    df["iter"] = iter_col
    return df

def _pack(file_glob, reader_fn, out_path, label, n_jobs=-1):
    files = sorted(glob.glob(file_glob))
    assert files, f"No files matched: {file_glob}"
    log_info(f"Packing {len(files):,} {label} Parquet files...")
    chunks = Parallel(n_jobs=n_jobs, backend="loky")(delayed(reader_fn)(fp) for fp in files)
    big = pd.concat(chunks, ignore_index=True)
    big.to_parquet(out_path, index=False, engine="pyarrow", compression="snappy", row_group_size=50_000)
    log_success(f"Packed {label} to {out_path} ({os.path.getsize(out_path)/1e9:.2f} GB)")

# ── Data Loading & Memory Profiling ───────────────────────────────────────────
if not os.path.exists(od_packed_path): _pack(os.path.join(od_dir, "OD_*.parquet"), _read_one_od, od_packed_path, "OD")
if not os.path.exists(sd_packed_path): _pack(os.path.join(sd_dir, "SD_*.parquet"), _read_one_sd, sd_packed_path, "SD")

log_info("Loading packed datasets into memory...")
od_df = pq.read_table(od_packed_path).to_pandas()
sd_df = pq.read_table(sd_packed_path).to_pandas()

for _df in [od_df, sd_df]: _df["var_type"] = _df["var_type"].astype("category")
sd_df["method"] = sd_df["method"].astype("category")

od_ram = od_df.memory_usage(deep=True).sum() / 1e9
sd_ram = sd_df.memory_usage(deep=True).sum() / 1e9
log_success(f"Loaded OD: {od_df.shape} ({od_ram:.2f} GB) | SD: {sd_df.shape} ({sd_ram:.2f} GB)")

# ── Shared Memory Assembly ────────────────────────────────────────────────────
# Cleanup stale memory segments from previous crashed runs
for shm_name in ["shm_od", "shm_sd"]:
    try:
        SharedMemory(name=shm_name).unlink()
    except FileNotFoundError:
        pass

def _col_to_int_key(series):
    vals = series.to_numpy()
    if vals.dtype.kind in ("U", "O"):
        return np.unique(vals, return_inverse=True)[1].astype(np.int64)
    return vals.astype(np.float64 if vals.dtype.kind == "f" else np.int64)

def _build_arena(df, meta_cols, feat_cols, max_p, arena_cols, label):
    n_rows = len(df)
    int_keys = [_col_to_int_key(df[c]) for c in meta_cols]
    sort_idx = np.lexsort(int_keys[::-1])

    arena = np.empty((n_rows, arena_cols), dtype=np.float64)
    arena[:, :max_p] = df[feat_cols].to_numpy(dtype=np.float64, na_value=np.nan)[sort_idx]
    arena[:, max_p]  = df["y"].to_numpy(dtype=np.float64)[sort_idx]

    sorted_int_keys = [k[sort_idx] for k in int_keys]
    boundary_masks = [sorted_int_keys[i][1:] != sorted_int_keys[i][:-1] for i in range(len(meta_cols))]
    boundaries = np.where(np.logical_or.reduce(boundary_masks))[0] + 1
    starts, ends = np.concatenate([[0], boundaries]), np.concatenate([boundaries, [n_rows]])

    orig_col_arrays = [df[c].to_numpy()[sort_idx] for c in meta_cols]
    idx_map = {tuple(arr[s].item() if hasattr(arr[s], "item") else arr[s] for arr in orig_col_arrays): (int(s), int(e)) for s, e in zip(starts, ends)}

    return arena, idx_map

log_info("Constructing contiguous memory arenas for zero-copy workers...")
od_arena, od_index = _build_arena(od_df, OD_META, FEAT_COLS, MAX_P, ARENA_COLS, "OD")
sd_arena, sd_index = _build_arena(sd_df, SD_META, FEAT_COLS, MAX_P, ARENA_COLS, "SD")

work_items, skipped_no_od = [], 0
for key, (sd_s, sd_e) in sd_index.items():
    method, vtype, N, p_val, rho, sigma, it = key
    od_key = (vtype, N, p_val, rho, sigma, it)
    if od_key not in od_index:
        skipped_no_od += 1
        continue
    od_s, od_e = od_index[od_key]
    work_items.append((od_s, od_e, sd_s, sd_e, int(p_val), key))

del sd_index, od_index, od_df, sd_df
gc.collect()

if skipped_no_od > 0:
    log_warn(f"Dropped {skipped_no_od} scenarios due to missing Original Data.")

shm_od_req = od_arena.nbytes / 1e9
shm_sd_req = sd_arena.nbytes / 1e9
log_info(f"Allocating System Shared Memory: OD ({shm_od_req:.2f} GB) + SD ({shm_sd_req:.2f} GB)...")

shm_od = SharedMemory(name="shm_od", create=True, size=od_arena.nbytes)
shm_sd = SharedMemory(name="shm_sd", create=True, size=sd_arena.nbytes)

np.ndarray(od_arena.shape, dtype=np.float64, buffer=shm_od.buf)[:] = od_arena
np.ndarray(sd_arena.shape, dtype=np.float64, buffer=shm_sd.buf)[:] = sd_arena

OD_SHAPE, SD_SHAPE = od_arena.shape, sd_arena.shape
del od_arena, sd_arena
gc.collect()

# ── Parallel OLS Worker ───────────────────────────────────────────────────────
def _pad(arr, length=MAX_P + 1):
    out = np.full(length, np.nan)
    out[: len(arr)] = arr
    return out

def _fit_ols_arrays(X: np.ndarray, y: np.ndarray):
    import numpy as _np
    import statsmodels.api as _sm
    X_const = _sm.add_constant(X, has_constant="add")
    try:
        res = _sm.OLS(y, X_const).fit()
        return res.params, res.bse, res.pvalues, res.rsquared_adj
    except Exception:
        return _np.full(X_const.shape[1], _np.nan), _np.full(X_const.shape[1], _np.nan), _np.full(X_const.shape[1], _np.nan), _np.nan

def process_work_item(item, od_shape, sd_shape, max_feat):
    import numpy as _np
    from multiprocessing.shared_memory import SharedMemory as _SHM

    od_s, od_e, sd_s, sd_e, p_val, meta_key = item
    method, vtype, N, p_val, rho, sigma, it = meta_key

    _shm_od = _SHM(name="shm_od", create=False, track=False)
    _shm_sd = _SHM(name="shm_sd", create=False, track=False)

    arena_od = _np.ndarray(od_shape, dtype=_np.float64, buffer=_shm_od.buf)
    arena_sd = _np.ndarray(sd_shape, dtype=_np.float64, buffer=_shm_sd.buf)

    X_od = _np.array(arena_od[od_s:od_e, :p_val], dtype=_np.float64)
    y_od = _np.array(arena_od[od_s:od_e, max_feat], dtype=_np.float64)
    X_sd = _np.array(arena_sd[sd_s:sd_e, :p_val], dtype=_np.float64)
    y_sd = _np.array(arena_sd[sd_s:sd_e, max_feat], dtype=_np.float64)

    _shm_od.close()
    _shm_sd.close()

    beta_od, se_od, pval_od, adj_r2_od = _fit_ols_arrays(X_od, y_od)
    beta_sd, se_sd, pval_sd, adj_r2_sd = _fit_ols_arrays(X_sd, y_sd)

    # Revert the -1.0 sentinel back to NaN for correct downstream logic
    sigma_out = np.nan if sigma == -1.0 else sigma

    row = {
        "method": method, "var_type": vtype, "N": N, "p": p_val,
        "rho": rho, "sigma_2": sigma_out, "iter": it,
        "adj_r2_od": adj_r2_od, "adj_r2_sd": adj_r2_sd,
    }
    for i, v in enumerate(_pad(beta_od)): row[f"beta_od_{i}"] = v
    for i, v in enumerate(_pad(beta_sd)): row[f"beta_sd_{i}"] = v
    for i, v in enumerate(_pad(se_od)):   row[f"se_od_{i}"]   = v
    for i, v in enumerate(_pad(se_sd)):   row[f"se_sd_{i}"]   = v
    for i, v in enumerate(_pad(pval_od)): row[f"pval_od_{i}"] = v
    for i, v in enumerate(_pad(pval_sd)): row[f"pval_sd_{i}"] = v
    return row

# ── Execute ───────────────────────────────────────────────────────────────────
log_info(f"Dispatching {len(work_items):,} paired OLS regressions across all cores...")
t_start = time.time()

raw_results = Parallel(n_jobs=-1, verbose=0, batch_size="auto", backend="loky")(
    delayed(process_work_item)(item, OD_SHAPE, SD_SHAPE, MAX_P) for item in work_items
)

elapsed = time.time() - t_start
log_success(f"OLS Fitting complete in {elapsed:.1f} seconds! ({(len(work_items)*2)/elapsed:.0f} models/sec)")

log_info("Releasing OS Shared Memory...")
shm_od.close(); shm_od.unlink()
shm_sd.close(); shm_sd.unlink()

# ── Save Results ──────────────────────────────────────────────────────────────
df = pd.DataFrame([r for r in raw_results if r is not None])
del raw_results
gc.collect()

nan_od = df["adj_r2_od"].isna().sum()
nan_sd = df["adj_r2_sd"].isna().sum()
if nan_od > 0 or nan_sd > 0:
    log_warn(f"Failed Matrix Inversions (Degenerate data): OD={nan_od}, SD={nan_sd}")

parquet_path = os.path.join(result_dir, "aggregated_model_metrics.parquet")
df.to_parquet(parquet_path, index=False, engine="pyarrow")

log_success(f"Final metrics saved to {parquet_path} ({os.path.getsize(parquet_path)/1e6:.1f} MB)")
print(f"{Colors.HEADER}{Colors.BOLD}=============================================================================={Colors.ENDC}\n")