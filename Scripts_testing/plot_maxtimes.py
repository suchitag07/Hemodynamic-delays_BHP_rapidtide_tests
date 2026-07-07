import os
import json
import gzip
import argparse
import subprocess

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import kurtosis, kurtosistest
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, FuncFormatter


RAPIDTIDE_SIF = "/Software/sif/rapidtide_latest_3.1.11.sif"


def ensure_hist_exists(subj_id):
    subj = subj_id if subj_id.startswith("sub-") else f"sub-{subj_id}"

    subj_dir = os.path.join(
        "/Outputs/main/rapidtide_3.1.11/timepoint_1/spatialfilt4",
        subj,
    )

    hist_base = os.path.join(subj_dir, f"{subj}_desc-maxtime_hist")
    json_path = hist_base + ".json"
    tsv_path_gz = hist_base + ".tsv.gz"
    tsv_path = hist_base + ".tsv"

    if os.path.exists(json_path) and (os.path.exists(tsv_path) or os.path.exists(tsv_path_gz)):
        print(f"[INFO] Using existing hist files for {subj}")
        return hist_base

    print(f"[INFO] Running histnifti for {subj}")

    maxtime_nii = os.path.join(subj_dir, f"{subj}_desc-maxtime_map.nii.gz")
    mask_nii = os.path.join(subj_dir, f"{subj}_desc-corrfit_mask.nii.gz")
    out_root = hist_base
    bind_path = subj_dir

    cmd = [
        "apptainer", "exec",
        "--bind", f"{bind_path}:{bind_path}",
        RAPIDTIDE_SIF,
        "histnifti",
        maxtime_nii,
        out_root,
        "--maskfile", mask_nii,
        "--nodisplay",
    ]
    subprocess.check_call(cmd)

    return hist_base


def load_hist(hist_base):
    json_path = hist_base + ".json"
    tsv_path = hist_base + ".tsv"
    tsv_path_gz = hist_base + ".tsv.gz"

    with open(json_path, "r") as f:
        md = json.load(f)

    sfreq = md["SamplingFrequency"]
    t0 = md["StartTime"]

    if os.path.exists(tsv_path):
        y = pd.read_csv(tsv_path, sep="\t", header=None)[0].to_numpy()
    elif os.path.exists(tsv_path_gz):
        with gzip.open(tsv_path_gz, "rt") as f:
            y = pd.read_csv(f, sep="\t", header=None)[0].to_numpy()
    else:
        raise FileNotFoundError(f"Could not find {tsv_path} or {tsv_path_gz}")

    n = len(y)
    x = t0 + np.arange(n) / sfreq
    return x, y, md


def build_png(subjects, out_png, normalize=False, dpi=300):
    fig, ax = plt.subplots(figsize=(6, 4), dpi=dpi)

    for subj_id in subjects:
        hist_base = ensure_hist_exists(subj_id)
        x, y, md = load_hist(hist_base)

        y_plot = y.astype(float).copy()
        if normalize:
            total = y_plot.sum()
            if total > 0:
                y_plot = y_plot / total

        ax.plot(x, y_plot, color="steelblue", alpha=0.5, linewidth=0.8)

    ax.axvline(0.5, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Delay time (s)")
    ax.set_ylabel("Proportion of voxels" if normalize else "Voxel count")
    ax.set_title("Delay time distribution")

    if not normalize:
        ax.yaxis.set_major_locator(MultipleLocator(2000))
        ax.yaxis.set_major_formatter(
            FuncFormatter(lambda v, pos: f"{int(v):d}" if v >= 0 else "")
        )
        ax.set_ylim(bottom=0)

    legend_handles = [
        Line2D([0], [0], color="steelblue", lw=2, label=f"Subjects (n={len(subjects)})"),
    ]
    ax.legend(handles=legend_handles, loc="best", frameon=False)

    fig.tight_layout()
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(
        description="Overlay maxtime histograms for multiple subjects; run histnifti if needed."
    )
    parser.add_argument(
        "--subjects",
        nargs="+",
        required=True,
        help="Subject IDs with or without 'sub-' prefix.",
    )
    parser.add_argument(
        "--out-png",
        default="maxtime_distribution_spatialfilt4_10.png",
        help="Output PNG filename.",
    )
    parser.add_argument(
        "--normalize",
        action="store_true",
        help="Normalize each histogram to sum to 1 before plotting.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="P-value threshold for kurtosistest abnormality.",
    )
    args = parser.parse_args()

    build_png(
        subjects=args.subjects,
        out_png=args.out_png,
        normalize=args.normalize,
    )


if __name__ == "__main__":
    main()