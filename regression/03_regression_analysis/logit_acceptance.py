#!/usr/bin/env python3
"""Shared logistic fit acceptance logic for generation and analysis."""

from __future__ import annotations

import argparse
import json
import warnings
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.api as sm

MAXITER = 200
UNSTABLE_BETA_ABS = 20.0
UNSTABLE_BSE_ABS = 10.0


def fit_logit_with_status(X: np.ndarray, y: np.ndarray) -> dict[str, Any]:
    """Fit Logit and return coefficients, diagnostics, and acceptance status."""
    Xc = sm.add_constant(X, has_constant="add")
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = sm.Logit(y, Xc).fit(disp=0, maxiter=MAXITER)

        warned_separation = any(
            issubclass(w.category, sm.tools.sm_exceptions.PerfectSeparationWarning)
            for w in caught
        )
        warned_hessian = any(
            issubclass(w.category, sm.tools.sm_exceptions.HessianInversionWarning)
            for w in caught
        )
        bse_nan = np.all(np.isnan(result.bse))
        bse_huge = np.any(np.isfinite(result.bse) & (np.abs(result.bse) > UNSTABLE_BSE_ABS))
        if warned_separation or bse_nan:
            status = "separation"
        elif warned_hessian:
            status = "hessian"
        elif bse_huge:
            status = "unstable_se"
        elif np.any(np.abs(result.params) > UNSTABLE_BETA_ABS):
            status = "unstable"
        else:
            status = "ok"

        accepted = status == "ok"
        if accepted:
            accuracy = np.mean((result.predict(Xc) >= 0.5).astype(int) == y.astype(int))
            return {
                "status": status,
                "accepted": True,
                "params": np.asarray(result.params, dtype=float),
                "bse": np.asarray(result.bse, dtype=float),
                "pvalues": np.asarray(result.pvalues, dtype=float),
                "pseudo_r2": float(result.prsquared),
                "accuracy": float(accuracy),
            }

        k = Xc.shape[1]
        nan_k = np.full(k, np.nan)
        return {
            "status": status,
            "accepted": False,
            "params": nan_k,
            "bse": nan_k.copy(),
            "pvalues": nan_k.copy(),
            "pseudo_r2": np.nan,
            "accuracy": np.nan,
        }
    except sm.tools.sm_exceptions.PerfectSeparationError:
        status = "separation"
    except Exception:
        status = "failed"

    k = Xc.shape[1]
    nan_k = np.full(k, np.nan)
    return {
        "status": status,
        "accepted": False,
        "params": nan_k,
        "bse": nan_k.copy(),
        "pvalues": nan_k.copy(),
        "pseudo_r2": np.nan,
        "accuracy": np.nan,
    }


def evaluate_iter_groups(frame: pd.DataFrame, iter_col: str = "iter") -> list[dict[str, Any]]:
    """Run the shared fit rule for each iteration group in a frame."""
    if iter_col not in frame.columns:
        raise ValueError(f"Missing iteration column: {iter_col}")

    x_cols = sorted(
        (col for col in frame.columns if col.startswith("X")),
        key=lambda col: int(col[1:]),
    )
    if not x_cols:
        raise ValueError("No predictor columns found")
    if "y" not in frame.columns:
        raise ValueError("Missing response column 'y'")

    results: list[dict[str, Any]] = []
    for iter_value, grp in frame.groupby(iter_col, sort=True):
        fit = fit_logit_with_status(
            grp[x_cols].to_numpy(dtype=np.float64),
            grp["y"].to_numpy(dtype=np.float64),
        )
        results.append(
            {
                "iter": int(iter_value),
                "status": fit["status"],
                "accepted": bool(fit["accepted"]),
            }
        )
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to parquet file with X*, y, iter columns")
    parser.add_argument("--output", required=True, help="Path to JSON results file")
    parser.add_argument("--iter-col", default="iter", help="Iteration column name")
    args = parser.parse_args()

    frame = pd.read_parquet(args.input)
    results = evaluate_iter_groups(frame, iter_col=args.iter_col)
    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(results, handle)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
