#!/usr/bin/env python3
"""Fit runtime approximations from per-run timing history."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
HISTORY_CSV = ROOT / "results" / "pipeline_timing_history.csv"
MODEL_JSON = ROOT / "results" / "latest_runtime_model.json"
MODEL_TXT = ROOT / "results" / "latest_runtime_model.txt"

STRUCTURAL_FIELDS = [
    "num_n",
    "num_p",
    "num_rho",
    "num_sigma2",
    "num_p1",
    "od_linear_cont_scenarios",
    "od_linear_bin_scenarios",
    "od_logistic_cont_scenarios",
    "od_logistic_bin_scenarios",
    "sd_linear_expected_files",
    "sd_logistic_expected_files",
]

ANSI = {
    "reset": "\033[0m",
    "blue": "\033[38;5;39m",
    "green": "\033[38;5;40m",
    "yellow": "\033[38;5;220m",
    "orange": "\033[38;5;208m",
    "red": "\033[38;5;196m",
    "magenta": "\033[38;5;171m",
    "cyan": "\033[38;5;45m",
    "teal": "\033[38;5;44m",
    "violet": "\033[38;5;141m",
    "gray": "\033[38;5;245m",
}


def use_color() -> bool:
    return os.isatty(1)


def colorize(text: str, color: str, enabled: bool | None = None) -> str:
    if enabled is None:
        enabled = use_color()
    if not enabled:
        return text
    return f"{ANSI[color]}{text}{ANSI['reset']}"


def load_history() -> pd.DataFrame:
    if not HISTORY_CSV.exists():
        raise FileNotFoundError(
            f"Timing history not found: {HISTORY_CSV}. Run ./run_all.sh at least once first."
        )

    df = pd.read_csv(HISTORY_CSV)
    if df.empty:
        raise ValueError(f"Timing history is empty: {HISTORY_CSV}")

    numeric_cols = [c for c in df.columns if c not in {"run_timestamp_barcelona", "config_hash", "config_changed"}]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["run_timestamp_barcelona"] = pd.to_datetime(
        df["run_timestamp_barcelona"],
        format="%Y-%m-%d-%H-%M-%S",
        errors="coerce",
    )
    return df.sort_values("run_timestamp_barcelona").reset_index(drop=True)


def select_latest_family(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    latest = df.iloc[-1]
    mask = np.ones(len(df), dtype=bool)
    for field in STRUCTURAL_FIELDS:
        mask &= df[field].eq(latest[field])

    family = df.loc[mask].copy().sort_values("run_timestamp_barcelona")
    family = (
        family
        .drop_duplicates(subset=["sim_m"], keep="last")
        .sort_values("sim_m")
        .reset_index(drop=True)
    )
    signature = {field: latest[field] for field in STRUCTURAL_FIELDS}
    return family, signature


def fit_affine(x: pd.Series, y: pd.Series) -> dict:
    clean = pd.DataFrame({"x": x, "y": y}).dropna()
    n = len(clean)
    if n == 0:
        return {"status": "no_data"}

    xv = clean["x"].to_numpy(dtype=float)
    yv = clean["y"].to_numpy(dtype=float)

    if n == 1 or np.allclose(xv, xv[0]):
        intercept = float(np.mean(yv))
        slope = 0.0
        pred = np.full_like(yv, intercept)
    else:
        X = np.column_stack([np.ones(n), xv])
        coef, *_ = np.linalg.lstsq(X, yv, rcond=None)
        intercept = float(coef[0])
        slope = float(coef[1])
        pred = X @ coef

    sse = float(np.sum((yv - pred) ** 2))
    sst = float(np.sum((yv - yv.mean()) ** 2))
    r2 = None if np.isclose(sst, 0.0) else float(1.0 - sse / sst)

    return {
        "status": "ok",
        "n_rows": int(n),
        "intercept": intercept,
        "slope": slope,
        "mean_abs_error_s": float(np.mean(np.abs(yv - pred))),
        "r2": r2,
    }


def equation_str(name: str, model: dict) -> str:
    if model["status"] != "ok":
        return f"{name}: insufficient data"

    eq = f"{model['intercept']:.3f} + {model['slope']:.3f}*M"
    suffix = f"  |  n={model['n_rows']}, MAE={model['mean_abs_error_s']:.2f}s"
    if model["r2"] is not None:
        suffix += f", R²={model['r2']:.3f}"
    return f"{name}(M) ≈ {eq}{suffix}"


def predict_affine(model: dict, m_value: float) -> float | None:
    if model["status"] != "ok":
        return None
    return float(model["intercept"] + model["slope"] * m_value)


def format_hms(total_seconds: float) -> str:
    total_seconds = max(0, int(round(total_seconds)))
    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    seconds = total_seconds % 60
    return f"{hours:02d}h {minutes:02d}m {seconds:02d}s"


def build_bar(segments: list[tuple[str, float, str]], width: int = 60, colored: bool = True) -> str:
    total = sum(value for _, value, _ in segments)
    if total <= 0:
        return "-" * width

    raw = [(label, value / total * width, color) for label, value, color in segments]
    ints = [int(math.floor(amount)) for _, amount, _ in raw]
    remainder = width - sum(ints)
    order = sorted(
        range(len(raw)),
        key=lambda idx: raw[idx][1] - ints[idx],
        reverse=True,
    )
    for idx in order[:remainder]:
        ints[idx] += 1

    pieces = []
    for (label, _, color), count in zip(raw, ints):
        chunk = "█" * count if count > 0 else ""
        pieces.append(colorize(chunk, color, enabled=colored))
    return "".join(pieces)


def summarize_composition(parts: dict[str, float], colored: bool, title_prefix: str = "Latest Run") -> tuple[list[str], dict]:
    step0 = float(parts["step0_duration_s"])
    od = float(parts["od_duration_s"])
    sd = float(parts["sd_duration_s"])
    lr = float(parts["lr_duration_s"])
    pipeline = od + sd + lr
    od_linear_cont = float(parts["od_linear_cont_duration_s"])
    od_linear_bin = float(parts["od_linear_bin_duration_s"])
    od_logistic_cont = float(parts["od_logistic_cont_duration_s"])
    od_logistic_bin = float(parts["od_logistic_bin_duration_s"])
    sd_linear_cont = float(parts.get("sd_linear_cont_duration_s", 0.0))
    sd_linear_bin = float(parts.get("sd_linear_bin_duration_s", 0.0))
    sd_logistic_cont = float(parts.get("sd_logistic_cont_duration_s", 0.0))
    sd_logistic_bin = float(parts.get("sd_logistic_bin_duration_s", 0.0))
    lr_linear = float(parts["lr_linear_duration_s"])
    lr_logistic = float(parts["lr_logistic_duration_s"])

    group_bar = build_bar(
        [
            ("OD", od, "green"),
            ("SD", sd, "yellow"),
            ("LR", lr, "cyan"),
        ],
        colored=colored,
    )
    detailed_bar = build_bar(
        [
            ("OD linear cont", od_linear_cont, "green"),
            ("OD linear bin", od_linear_bin, "blue"),
            ("OD logistic cont", od_logistic_cont, "orange"),
            ("OD logistic bin", od_logistic_bin, "red"),
            ("SD linear cont", sd_linear_cont, "yellow"),
            ("SD linear bin", sd_linear_bin, "green"),
            ("SD logistic cont", sd_logistic_cont, "cyan"),
            ("SD logistic bin", sd_logistic_bin, "magenta"),
            ("LR linear", lr_linear, "teal"),
            ("LR logistic", lr_logistic, "violet"),
        ],
        colored=colored,
    )

    def pct(value: float) -> float:
        return 0.0 if pipeline <= 0 else 100.0 * value / pipeline

    lines = [
        f"{title_prefix} Composition (normalized to OD + SD + LR = 100%)",
        f"  {group_bar}",
        (
            f"  {colorize('OD', 'green', enabled=colored)}: {od:.1f}s ({pct(od):.1f}%)   "
            f"{colorize('SD', 'yellow', enabled=colored)}: {sd:.1f}s ({pct(sd):.1f}%)   "
            f"{colorize('LR', 'cyan', enabled=colored)}: {lr:.1f}s ({pct(lr):.1f}%)"
        ),
        "",
        f"{title_prefix} Detailed Composition",
        f"  {detailed_bar}",
        (
            f"  {colorize('OD linear cont', 'green', enabled=colored)}: {od_linear_cont:.1f}s ({pct(od_linear_cont):.1f}%)   "
            f"{colorize('OD linear bin', 'blue', enabled=colored)}: {od_linear_bin:.1f}s ({pct(od_linear_bin):.1f}%)"
        ),
        (
            f"  {colorize('OD logistic cont', 'orange', enabled=colored)}: {od_logistic_cont:.1f}s ({pct(od_logistic_cont):.1f}%)   "
            f"{colorize('OD logistic bin', 'red', enabled=colored)}: {od_logistic_bin:.1f}s ({pct(od_logistic_bin):.1f}%)"
        ),
        (
            f"  {colorize('SD linear cont', 'yellow', enabled=colored)}: {sd_linear_cont:.1f}s ({pct(sd_linear_cont):.1f}%)   "
            f"{colorize('SD linear bin', 'green', enabled=colored)}: {sd_linear_bin:.1f}s ({pct(sd_linear_bin):.1f}%)"
        ),
        (
            f"  {colorize('SD logistic cont', 'cyan', enabled=colored)}: {sd_logistic_cont:.1f}s ({pct(sd_logistic_cont):.1f}%)   "
            f"{colorize('SD logistic bin', 'magenta', enabled=colored)}: {sd_logistic_bin:.1f}s ({pct(sd_logistic_bin):.1f}%)"
        ),
        (
            f"  {colorize('SD total', 'yellow', enabled=colored)}: {sd:.1f}s ({pct(sd):.1f}%)   "
            f"{colorize('LR linear', 'teal', enabled=colored)}: {lr_linear:.1f}s ({pct(lr_linear):.1f}%)   "
            f"{colorize('LR logistic', 'violet', enabled=colored)}: {lr_logistic:.1f}s ({pct(lr_logistic):.1f}%)"
        ),
        f"  {colorize('Step 0', 'gray', enabled=colored)}: {step0:.1f}s (excluded from the two bars above)",
    ]

    payload = {
        "step0_duration_s": step0,
        "od_duration_s": od,
        "sd_duration_s": sd,
        "lr_duration_s": lr,
        "pipeline_without_step0_duration_s": pipeline,
    }
    return lines, payload


def summarize_latest_run(latest: pd.Series, colored: bool) -> tuple[list[str], dict]:
    parts = {
        "step0_duration_s": float(latest["step0_duration_s"]),
        "od_duration_s": float(latest["od_linear_total_duration_s"] + latest["od_logistic_total_duration_s"]),
        "sd_duration_s": float(latest["sd_duration_s"]),
        "lr_duration_s": float(latest["lr_linear_duration_s"] + latest["lr_logistic_duration_s"]),
        "od_linear_cont_duration_s": float(latest["od_linear_cont_duration_s"]),
        "od_linear_bin_duration_s": float(latest["od_linear_bin_duration_s"]),
        "od_logistic_cont_duration_s": float(latest["od_logistic_cont_duration_s"]),
        "od_logistic_bin_duration_s": float(latest["od_logistic_bin_duration_s"]),
        "sd_linear_cont_duration_s": float(latest.get("sd_linear_cont_duration_s", 0.0)),
        "sd_linear_bin_duration_s": float(latest.get("sd_linear_bin_duration_s", 0.0)),
        "sd_logistic_cont_duration_s": float(latest.get("sd_logistic_cont_duration_s", 0.0)),
        "sd_logistic_bin_duration_s": float(latest.get("sd_logistic_bin_duration_s", 0.0)),
        "lr_linear_duration_s": float(latest["lr_linear_duration_s"]),
        "lr_logistic_duration_s": float(latest["lr_logistic_duration_s"]),
    }
    return summarize_composition(parts, colored=colored, title_prefix="Latest Run")


def build_models(family: pd.DataFrame) -> dict[str, dict]:
    return {
        "pipeline_total": fit_affine(family["sim_m"], family["pipeline_total_duration_s"]),
        "od_total": fit_affine(
            family["sim_m"],
            family["od_linear_total_duration_s"] + family["od_logistic_total_duration_s"],
        ),
        "sd": fit_affine(family["sim_m"], family["sd_duration_s"]),
        "lr_total": fit_affine(
            family["sim_m"],
            family["lr_linear_duration_s"] + family["lr_logistic_duration_s"],
        ),
        "od_linear": fit_affine(family["sim_m"], family["od_linear_total_duration_s"]),
        "od_logistic": fit_affine(family["sim_m"], family["od_logistic_total_duration_s"]),
        "od_linear_cont": fit_affine(family["sim_m"], family["od_linear_cont_duration_s"]),
        "od_linear_bin": fit_affine(family["sim_m"], family["od_linear_bin_duration_s"]),
        "od_logistic_cont": fit_affine(family["sim_m"], family["od_logistic_cont_duration_s"]),
        "od_logistic_bin": fit_affine(family["sim_m"], family["od_logistic_bin_duration_s"]),
        "sd_linear_cont": fit_affine(family["sim_m"], family.get("sd_linear_cont_duration_s", pd.Series(dtype=float))),
        "sd_linear_bin": fit_affine(family["sim_m"], family.get("sd_linear_bin_duration_s", pd.Series(dtype=float))),
        "sd_logistic_cont": fit_affine(family["sim_m"], family.get("sd_logistic_cont_duration_s", pd.Series(dtype=float))),
        "sd_logistic_bin": fit_affine(family["sim_m"], family.get("sd_logistic_bin_duration_s", pd.Series(dtype=float))),
        "lr_linear": fit_affine(family["sim_m"], family["lr_linear_duration_s"]),
        "lr_logistic": fit_affine(family["sim_m"], family["lr_logistic_duration_s"]),
    }


def observed_points_lines(family: pd.DataFrame) -> list[str]:
    lines = ["Observed timing points used for the fit:"]
    for _, row in family.iterrows():
        od = float(row["od_linear_total_duration_s"] + row["od_logistic_total_duration_s"])
        sd = float(row["sd_duration_s"])
        lr = float(row["lr_linear_duration_s"] + row["lr_logistic_duration_s"])
        total = float(row["pipeline_total_duration_s"])
        lines.append(
            f"  M={int(row['sim_m'])}: total={total:.1f}s, OD={od:.1f}s, SD={sd:.1f}s, LR={lr:.1f}s"
        )
    return lines


def prompt_for_m() -> float | None:
    try:
        raw = input("Target M to predict (blank to skip): ").strip()
    except EOFError:
        return None
    if raw == "":
        return None
    return float(raw)


def prediction_lines(models: dict[str, dict], m_value: float, colored: bool) -> tuple[list[str], dict]:
    pred_parts = {
        "step0_duration_s": 0.0,
        "od_linear_cont_duration_s": max(0.0, predict_affine(models["od_linear_cont"], m_value) or 0.0),
        "od_linear_bin_duration_s": max(0.0, predict_affine(models["od_linear_bin"], m_value) or 0.0),
        "od_logistic_cont_duration_s": max(0.0, predict_affine(models["od_logistic_cont"], m_value) or 0.0),
        "od_logistic_bin_duration_s": max(0.0, predict_affine(models["od_logistic_bin"], m_value) or 0.0),
        "sd_linear_cont_duration_s": max(0.0, predict_affine(models["sd_linear_cont"], m_value) or 0.0),
        "sd_linear_bin_duration_s": max(0.0, predict_affine(models["sd_linear_bin"], m_value) or 0.0),
        "sd_logistic_cont_duration_s": max(0.0, predict_affine(models["sd_logistic_cont"], m_value) or 0.0),
        "sd_logistic_bin_duration_s": max(0.0, predict_affine(models["sd_logistic_bin"], m_value) or 0.0),
        "lr_linear_duration_s": max(0.0, predict_affine(models["lr_linear"], m_value) or 0.0),
        "lr_logistic_duration_s": max(0.0, predict_affine(models["lr_logistic"], m_value) or 0.0),
    }
    pred_parts["od_duration_s"] = (
        pred_parts["od_linear_cont_duration_s"]
        + pred_parts["od_linear_bin_duration_s"]
        + pred_parts["od_logistic_cont_duration_s"]
        + pred_parts["od_logistic_bin_duration_s"]
    )
    pred_parts["sd_duration_s"] = (
        pred_parts["sd_linear_cont_duration_s"]
        + pred_parts["sd_linear_bin_duration_s"]
        + pred_parts["sd_logistic_cont_duration_s"]
        + pred_parts["sd_logistic_bin_duration_s"]
    )
    pred_parts["lr_duration_s"] = pred_parts["lr_linear_duration_s"] + pred_parts["lr_logistic_duration_s"]
    predicted_total_direct = max(0.0, predict_affine(models["pipeline_total"], m_value) or 0.0)
    composition_lines, payload = summarize_composition(
        pred_parts,
        colored=colored,
        title_prefix=f"Predicted Runtime for M={m_value:g}",
    )
    header = [
        "",
        "=" * 80,
        f"PREDICTION FOR M={m_value:g}",
        "=" * 80,
        f"  Direct total model: {predicted_total_direct:.1f}s",
        f"  Sum of predicted components: {payload['pipeline_without_step0_duration_s']:.1f}s",
        "",
    ]
    return header + composition_lines, {
        "m": m_value,
        "direct_total_duration_s": predicted_total_direct,
        "component_sum_duration_s": payload["pipeline_without_step0_duration_s"],
        "parts": pred_parts,
    }


def compact_observed_lines(family: pd.DataFrame) -> list[str]:
    lines = ["We currently have data for:"]
    for _, row in family.iterrows():
        total = float(row["pipeline_total_duration_s"])
        lines.append(f"  M = {int(row['sim_m'])}  ->  {total:.1f}s ({format_hms(total)})")
    return lines


def compact_formula_lines(models: dict[str, dict]) -> list[str]:
    model = models["pipeline_total"]
    if model["status"] != "ok":
        return ["We do not have enough data to fit T(M) yet."]

    lines = [
        "We assume the following approximation:",
        f"  T(M) ≈ {model['intercept']:.3f} + {model['slope']:.3f}*M  seconds",
        f"  Fit quality: n={model['n_rows']}, MAE={model['mean_abs_error_s']:.2f}s"
    ]
    if model["r2"] is not None:
        lines[-1] += f", R²={model['r2']:.3f}"
    return lines


def compact_prediction_lines(models: dict[str, dict], m_value: float, colored: bool) -> list[str]:
    pred_lines, pred_payload = prediction_lines(models, m_value, colored=colored)
    total = pred_payload["direct_total_duration_s"]

    lines = [
        "",
        "=" * 80,
        f"ESTIMATE FOR M = {m_value:g}",
        "=" * 80,
        f"Estimated total runtime: {total:.1f}s ({format_hms(total)})",
        f"Estimated OD + SD + LR sum: {pred_payload['component_sum_duration_s']:.1f}s",
        "",
    ]
    lines.extend(pred_lines[7:])
    return lines


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-prompt", action="store_true", help="Do not ask for a target M")
    parser.add_argument("--predict-m", type=float, default=None, help="Predict runtime for a target M")
    parser.add_argument("--estimator", action="store_true", help="Run compact estimator mode")
    args = parser.parse_args()

    try:
        history = load_history()
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}")
        raise SystemExit(1) from exc

    family, signature = select_latest_family(history)
    latest = family.iloc[-1]

    models = build_models(family)

    interactive_mode = args.estimator or (not args.no_prompt and args.predict_m is None and os.isatty(0))

    if interactive_mode:
        compact_lines = [
            "=" * 80,
            "PIPELINE RUNTIME ESTIMATOR",
            "=" * 80,
            *compact_observed_lines(family),
            "",
            *compact_formula_lines(models),
        ]
        for line in compact_lines:
            print(line)

        predict_m = prompt_for_m()
        if predict_m is None:
            return

        for line in compact_prediction_lines(models, predict_m, colored=True):
            print(line)
        return

    head_lines = [
        "=" * 80,
        "PIPELINE RUNTIME APPROXIMATION",
        "=" * 80,
        f"History file: {HISTORY_CSV}",
        f"Stored runs: {len(history)}",
        f"Modeling family size: {len(family)} unique M values",
        f"Latest run timestamp (Barcelona): {latest['run_timestamp_barcelona']}",
        "",
        "Assumption:",
        "  This approximation treats M as the only changing variable.",
        "  The fitted family keeps the latest run for each M under the same structural grid.",
        "",
        "Structural grid used for the fit:",
        (
            f"  num_n={int(signature['num_n'])}, num_p={int(signature['num_p'])}, "
            f"num_rho={int(signature['num_rho'])}, num_sigma2={int(signature['num_sigma2'])}, "
            f"num_p1={int(signature['num_p1'])}"
        ),
        "",
        "Approximation equations:",
        f"  {equation_str('T_total', models['pipeline_total'])}",
        f"  {equation_str('T_OD', models['od_total'])}",
        f"  {equation_str('T_SD', models['sd'])}",
        f"  {equation_str('T_LR', models['lr_total'])}",
        "",
        "Detailed sub-equations:",
        f"  {equation_str('T_OD_linear', models['od_linear'])}",
        f"  {equation_str('T_OD_logistic', models['od_logistic'])}",
        f"  {equation_str('T_OD_linear_cont', models['od_linear_cont'])}",
        f"  {equation_str('T_OD_linear_bin', models['od_linear_bin'])}",
        f"  {equation_str('T_OD_logistic_cont', models['od_logistic_cont'])}",
        f"  {equation_str('T_OD_logistic_bin', models['od_logistic_bin'])}",
        f"  {equation_str('T_SD_linear_cont', models['sd_linear_cont'])}",
        f"  {equation_str('T_SD_linear_bin', models['sd_linear_bin'])}",
        f"  {equation_str('T_SD_logistic_cont', models['sd_logistic_cont'])}",
        f"  {equation_str('T_SD_logistic_bin', models['sd_logistic_bin'])}",
        f"  {equation_str('T_LR_linear', models['lr_linear'])}",
        f"  {equation_str('T_LR_logistic', models['lr_logistic'])}",
        "",
    ]

    plain_bar_lines, latest_payload = summarize_latest_run(latest, colored=False)
    console_bar_lines, _ = summarize_latest_run(latest, colored=True)
    observed_lines = observed_points_lines(family)
    tail_lines = [
        "",
        *observed_lines,
        "",
        "Observed M values used:",
        "  " + ", ".join(str(int(v)) for v in family["sim_m"].tolist()),
        "",
        f"Saved model summary to: {MODEL_JSON}",
        f"Saved text summary to: {MODEL_TXT}",
    ]
    plain_lines = head_lines + plain_bar_lines + tail_lines

    report = {
        "history_csv": str(HISTORY_CSV),
        "stored_runs": int(len(history)),
        "modeling_family_size": int(len(family)),
        "latest_run_timestamp_barcelona": latest["run_timestamp_barcelona"].strftime("%Y-%m-%d-%H-%M-%S"),
        "structural_signature": {k: int(v) for k, v in signature.items()},
        "observed_m_values": [int(v) for v in family["sim_m"].tolist()],
        "models": models,
        "latest_run": latest_payload,
    }

    predict_m = args.predict_m
    if predict_m is not None:
        plain_pred_lines, pred_payload = prediction_lines(models, predict_m, colored=False)
        console_pred_lines, _ = prediction_lines(models, predict_m, colored=True)
        plain_lines += plain_pred_lines
        report["prediction"] = pred_payload
    else:
        console_pred_lines = []

    MODEL_JSON.write_text(json.dumps(report, indent=2), encoding="utf-8")
    MODEL_TXT.write_text("\n".join(plain_lines) + "\n", encoding="utf-8")

    console_lines = head_lines + console_bar_lines + tail_lines + console_pred_lines

    for line in console_lines:
        print(line)


if __name__ == "__main__":
    main()
