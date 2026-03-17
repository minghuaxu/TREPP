import argparse
import os
import stat
import subprocess
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from tqdm import tqdm
from catboost import CatBoostClassifier
from sklearn.calibration import CalibratedClassifierCV


def build_parser() -> argparse.ArgumentParser:

    # 获取当前脚本所在目录 (trepp/src)
    src_dir = Path(__file__).resolve().parent

    p = argparse.ArgumentParser(
        prog="trepp_predict.py",
        description="TREPP.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Core
    p.add_argument("--outdir", required=True, help="Output directory (all intermediate + final outputs go here).")
    
    default_model = str(src_dir / "models/trepp_production.pkl")
    p.add_argument("-m", "--model", default=default_model, help="Path to trepp_final.pkl.")

    # Input mode A: raw BED -> feature extraction
    p.add_argument("--input-bed", required=True, help="Input BED for prediction.")
    p.add_argument("--ref-fasta", required=True, help="Reference genome FASTA.")
    p.add_argument("--db-path", required=True, help="Path to processed annotations.db.")
    default_feature_script = str(Path(__file__).resolve().parent / "extract_features.py")

    default_trf = str(src_dir / "trf-4.10.0")
    p.add_argument("--trf-path", default=default_trf, help="Path to TRF binary.")

    p.add_argument(
        "--feature-script",
        default=str(src_dir / "extract_features.py"),
        help="Feature extraction script path (used with --input-bed).",
    )

    # Output naming
    p.add_argument(
        "--pred-out",
        default=None,
        help="Prediction output TSV path. If omitted, uses <outdir>/predictions.tsv.",
    )

    p.add_argument(
        "--threshold",
        type=float,
        default=0.5,
        help="Probability threshold for pred_label. If omitted, uses meta_info.best_threshold if present, else 0.5.",
    )
    p.add_argument("--no-progress", action="store_true", help="Disable progress bar.")

    return p


def ensure_outdir(outdir: str) -> Path:
    p = Path(outdir).resolve()
    p.mkdir(parents=True, exist_ok=True)
    (p / "logs").mkdir(parents=True, exist_ok=True)
    (p / "intermediate").mkdir(parents=True, exist_ok=True)
    return p


def make_executable(path: str) -> None:
    try:
        st = os.stat(path)
        os.chmod(path, st.st_mode | stat.S_IXUSR)
    except Exception:
        # Non-fatal; extraction may still work if already executable
        pass


def run_feature_extraction(
    *,
    input_bed: str,
    ref_fasta: str,
    db_path: str,
    trf_path: str,
    feature_script: str,
    outdir: Path,
) -> Path:
    log_path = outdir / "logs" / "extract_features.log"
    features_csv = outdir / "intermediate" / "predict_features.csv"

    make_executable(trf_path)

    cmd = [
        sys.executable,
        str(Path(feature_script).resolve()),
        "--input_bed",
        str(Path(input_bed).resolve()),
        "--ref_fasta",
        str(Path(ref_fasta).resolve()),
        "--db_path",
        str(Path(db_path).resolve()),
        "--trf_path",
        str(Path(trf_path).resolve()),
        "--output_csv",
        str(features_csv),
    ]

    with open(log_path, "w", encoding="utf-8") as log_f:
        proc = subprocess.run(cmd, stdout=log_f, stderr=subprocess.STDOUT)

    if proc.returncode != 0:
        raise RuntimeError(
            f"Feature extraction failed (exit={proc.returncode}). See log: {log_path}"
        )

    if not features_csv.exists() or features_csv.stat().st_size == 0:
        raise RuntimeError(f"Feature extraction produced empty file: {features_csv}")

    return features_csv


def load_artifact(model_path: str) -> dict:
    try:
        artifact = joblib.load(model_path)
    except Exception as e:
        raise RuntimeError(f"Failed to load model artifact: {e}") from e

    required_keys = ("base_models", "meta_info", "feat_cols")
    missing = [k for k in required_keys if k not in artifact]
    if missing:
        raise ValueError(f"Model file is missing keys: {missing} (expecting V7 artifact).")

    meta_info = artifact["meta_info"]
    if "meta_model" not in meta_info:
        raise ValueError("meta_info is missing key: meta_model")

    return artifact


def run_inference(X, artifact: dict, show_progress: bool) -> np.ndarray:
    base_models = artifact["base_models"]
    meta_info = artifact["meta_info"]
    meta_model = meta_info["meta_model"]

    X_array = X.values.astype(float) if isinstance(X, pd.DataFrame) else X.astype(float)
    n_samples = X_array.shape[0]
    n_models = len(base_models)
    base_preds = np.zeros((n_samples, n_models), dtype=float)

    iterator = enumerate(base_models)
    if show_progress:
        iterator = tqdm(iterator, total=n_models, desc="Base models")

    for i, model in iterator:
        base_preds[:, i] = model.predict_proba(X_array)[:, 1]

    return meta_model.predict_proba(base_preds)[:, 1]


def ensure_id(df: pd.DataFrame) -> pd.Series:
    if "id" in df.columns:
        return df["id"].astype(str)

    required = ("chr", "start", "end", "motif")
    if all(c in df.columns for c in required):
        return (
            df["chr"].astype(str)
            + "_"
            + df["start"].astype(str)
            + "_"
            + df["end"].astype(str)
            + "_"
            + df["motif"].astype(str)
        )

    return df.index.astype(str)


def main() -> int:
    args = build_parser().parse_args()
    outdir = ensure_outdir(args.outdir)

    pred_out = Path(args.pred_out).resolve() if args.pred_out else (outdir / "trepp_predicted.tsv")


    # Need BED mode
    needed = ("input_bed", "ref_fasta", "db_path", "trf_path")
    missing = [k for k in needed if getattr(args, k) in (None, "")]
    if missing:
        print(
            "[ERROR] Either provide --features, or provide --input-bed plus "
            "--ref-fasta --db-path .",
            file=sys.stderr,
        )
        return 1

    try:
        features_path = run_feature_extraction(
            input_bed=args.input_bed,
            ref_fasta=args.ref_fasta,
            db_path=args.db_path,
            trf_path=args.trf_path,
            feature_script=args.feature_script,
            outdir=outdir,
        )
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1
    sep = ","  # extract_features.py outputs CSV by design
    log_msg = f"Extracted features to: {features_path}"

    print(log_msg)

    # Load artifact
    try:
        artifact = load_artifact(args.model)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1

    # Load features table
    try:
        df = pd.read_csv(features_path, sep=sep)
    except Exception as e:
        print(f"[ERROR] Failed to read features table: {e}", file=sys.stderr)
        return 1

    feat_cols = artifact["feat_cols"]

    # Inference
    try:
        probs = run_inference(df[feat_cols], artifact, show_progress=(not args.no_progress))
    except Exception as e:
        print(f"[ERROR] Inference failed: {e}", file=sys.stderr)
        return 1

    threshold = args.threshold

    out = pd.DataFrame(
        {
            "id": ensure_id(df),
            "prob": probs,
            "pred_label": (probs >= threshold).astype(int),
        }
    )

    pred_out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(pred_out, sep="\t", index=False)

    print(f"Saved predictions: {pred_out}")
    print(f"Samples: {len(out)}")
    print(f"Threshold: {threshold}")
    print(f"Outdir: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
