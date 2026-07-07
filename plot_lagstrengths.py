import os
import json
import gzip
import argparse
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, FuncFormatter

RAPIDTIDE_SIF = "/Software/sif/rapidtide_latest_3.1.11.sif"


def ensure_hist_exists(subj_id):
    subj = subj_id if subj_id.startswith("sub-") else f"sub-{subj_id}"

    subj_dir = os.path.join(
        "/Outputs/main/rapidtide_3.1.11/timepoint_1/defaultfilt_5to15",
        subj,
    )

    hist_base = os.path.join(subj_dir, f"{subj}_desc-maxcorr_hist")
    json_path = hist_base + ".json"
    tsv_path_gz = hist_base + ".tsv.gz"
    tsv_path = hist_base + ".tsv"

    if os.path.exists(json_path) and (os.path.exists(tsv_path) or os.path.exists(tsv_path_gz)):
        print(f"[INFO] Using existing hist files for {subj}")
        return hist_base

    print(f"[INFO] Running histnifti for {subj}")

    maxcorr_nii = os.path.join(subj_dir, f"{subj}_desc-maxcorr_map.nii.gz")
    mask_nii = os.path.join(subj_dir, f"{subj}_desc-corrfit_mask.nii.gz")
    out_root = hist_base
    bind_path = subj_dir

    cmd = [
        "apptainer", "exec",
        "--bind", f"{bind_path}:{bind_path}",
        RAPIDTIDE_SIF,
        "histnifti",
        maxcorr_nii,
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


def weighted_skewness(x, y):
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)

    total = np.sum(y)
    if total <= 0:
        return np.nan

    mean = np.sum(y * x) / total
    var = np.sum(y * (x - mean) ** 2) / total
    if var <= 0:
        return 0.0

    std = np.sqrt(var)
    third = np.sum(y * (x - mean) ** 3) / total
    return third / (std ** 3)


def fmt_or_nan(val, ndigits=3):
    if val is None:
        return "nan"
    try:
        if np.isnan(val):
            return "nan"
    except TypeError:
        pass
    return f"{float(val):.{ndigits}f}"


def build_html(subjects, out_html, out_csv, normalize=False, skew_threshold=0.0):
    rows = []
    traces = []

    for subj_id in subjects:
        hist_base = ensure_hist_exists(subj_id)
        x, y, md = load_hist(hist_base)

        subj = subj_id if subj_id.startswith("sub-") else f"sub-{subj_id}"

        peakloc = md.get("peakloc", None)
        if peakloc is None:
            peakloc = float(x[np.argmax(y)])

        centerofmass = md.get("centerofmass", np.nan)
        pct50 = md.get("pct50", np.nan)
        skewness = weighted_skewness(x, y)

        group = "acceptable" if skewness <= skew_threshold else "review"
        color = "royalblue" if skewness <= skew_threshold else "firebrick"

        y_plot = y.astype(float).copy()
        if normalize:
            total = y_plot.sum()
            if total > 0:
                y_plot = y_plot / total

        hovertext = (
            f"<b>{subj}</b><br>"
            f"skewness={fmt_or_nan(skewness)}<br>"
            f"centerofmass={fmt_or_nan(centerofmass)}<br>"
            f"peakloc={fmt_or_nan(peakloc)}<br>"
            f"pct50={fmt_or_nan(pct50)}"
        )

        traces.append(
            {
                "subject_id": subj,
                "x": x.tolist(),
                "y": y_plot.tolist(),
                "color": color,
                "group": group,
                "hovertext": hovertext,
            }
        )

        rows.append(
            {
                "subject_id": subj,
                "peakloc": peakloc,
                "centerofmass": centerofmass,
                "pct50": pct50,
                "skewness": skewness,
                "group": group,
            }
        )

    df = pd.DataFrame(rows).sort_values("skewness")
    df.to_csv(out_csv, index=False)

    n_blue = int((df["group"] == "acceptable").sum())
    n_red = int((df["group"] == "review").sum())
    y_label = "Proportion of voxels" if normalize else "Voxel count"

    all_x = np.concatenate([np.asarray(t["x"], dtype=float) for t in traces])
    xmin = float(np.nanmin(all_x))
    xmax = float(np.nanmax(all_x))
    if xmax > xmin:
        xpad = 0.02 * (xmax - xmin)
    else:
        xpad = 0.01
    plot_xmin = max(0.0, xmin - xpad)
    plot_xmax = min(1.0, xmax + xpad)

    traces_json = json.dumps(traces)
    summary_json = json.dumps(df.to_dict(orient="records"), default=lambda x: None)

    if normalize:
        yaxis_config = f'''
    yaxis: {{
        title: "{y_label}"
    }},'''
    else:
        yaxis_config = f'''
    yaxis: {{
        title: "{y_label}",
        range: [0, 5000],
        tick0: 0,
        dtick: 1000
    }},'''

    html_template = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Maxcorr distributions</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
body {
    font-family: Arial, sans-serif;
    margin: 0;
    background: #fafafa;
    color: #222;
}
.wrap {
    max-width: 1400px;
    margin: 0 auto;
    padding: 20px;
}
h1 {
    font-size: 24px;
    margin: 0 0 8px;
}
.sub {
    color: #555;
    margin-bottom: 16px;
}
.controls {
    display: flex;
    gap: 12px;
    flex-wrap: wrap;
    align-items: center;
    margin-bottom: 12px;
}
.badge {
    display: inline-block;
    padding: 6px 10px;
    border-radius: 999px;
    background: #eee;
    font-size: 13px;
}
.panel {
    background: white;
    border: 1px solid #ddd;
    border-radius: 10px;
    padding: 12px;
    box-shadow: 0 1px 4px rgba(0,0,0,0.05);
}
#plot {
    width: 100%;
    height: 720px;
}
#info {
    margin-top: 12px;
    font-size: 14px;
}
table {
    width: 100%;
    border-collapse: collapse;
    margin-top: 14px;
    font-size: 13px;
}
th, td {
    border-bottom: 1px solid #e5e5e5;
    padding: 8px 10px;
    text-align: left;
}
th {
    background: #f5f5f5;
    position: sticky;
    top: 0;
}
.ok {
    color: royalblue;
    font-weight: 600;
}
.review {
    color: firebrick;
    font-weight: 600;
}
.hint {
    color: #666;
    font-size: 13px;
    margin-top: 8px;
}
</style>
</head>
<body>
<div class="wrap">
    <h1>Voxelwise maximum-correlation distributions</h1>
    <div class="sub">
        Blue: skewness &le; __SKEW_THRESHOLD__ (n=__N_BLUE__)
        &nbsp;&nbsp;
        Red: skewness &gt; __SKEW_THRESHOLD__ (n=__N_RED__)
    </div>

    <div class="controls">
        <span class="badge">Skewness threshold: __SKEW_THRESHOLD__</span>
        <span class="badge">Y-axis: __Y_LABEL__</span>
        <span class="badge">Click curve: highlight</span>
        <span class="badge">Double-click empty space: restore all</span>
    </div>

    <div class="panel">
        <div id="plot"></div>
        <div id="info">
            Click a curve to inspect that subject. Plotly legend click also hides/shows traces.
        </div>
    </div>

    <div class="panel" style="margin-top:14px; overflow:auto; max-height:420px;">
        <table>
            <thead>
                <tr>
                    <th>subject_id</th>
                    <th>group</th>
                    <th>skewness</th>
                    <th>centerofmass</th>
                    <th>peakloc</th>
                    <th>pct50</th>
                </tr>
            </thead>
            <tbody id="summary-body"></tbody>
        </table>
        <div class="hint">
            Negative skewness here generally means more mass at higher correlation values; positive skewness means more mass at lower values.
        </div>
    </div>
</div>

<script>
const tracesData = __TRACES_JSON__;
const summaryRows = __SUMMARY_JSON__;

const plotTraces = tracesData.map((t) => ({
    x: t.x,
    y: t.y,
    mode: "lines",
    type: "scattergl",
    name: t.subject_id,
    line: { color: t.color, width: 2 },
    opacity: 0.65,
    text: Array(t.x.length).fill(t.hovertext),
    hovertemplate: "%{text}<extra></extra>"
}));

const layout = {
    title: "Voxelwise maximum-correlation distributions",
    xaxis: {
        title: "Maximum correlation coefficient (r)",
        range: [__PLOT_XMIN__, __PLOT_XMAX__],
        tick0: 0,
        dtick: 0.1
    },
__YAXIS_CONFIG__
    shapes: [
        {
            type: "line",
            x0: 0.5,
            x1: 0.5,
            y0: 0,
            y1: 1,
            yref: "paper",
            line: {
                color: "black",
                width: 1,
                dash: "dash"
            }
        }
    ],
    hovermode: "closest",
    template: "plotly_white",
    legend: {
        orientation: "h",
        yanchor: "bottom",
        y: 1.02,
        xanchor: "left",
        x: 0
    }
};

Plotly.newPlot("plot", plotTraces, layout, { responsive: true });

function renderTable(highlightSubject = null) {
    const tbody = document.getElementById("summary-body");
    tbody.innerHTML = "";

    summaryRows.forEach(row => {
        const tr = document.createElement("tr");

        if (highlightSubject && row.subject_id === highlightSubject) {
            tr.style.background = "#fff8d6";
        }

        tr.innerHTML = `
            <td>${row.subject_id}</td>
            <td class="${row.group === "acceptable" ? "ok" : "review"}">${row.group}</td>
            <td>${row.skewness?.toFixed ? row.skewness.toFixed(4) : row.skewness}</td>
            <td>${row.centerofmass?.toFixed ? row.centerofmass.toFixed(4) : row.centerofmass}</td>
            <td>${row.peakloc?.toFixed ? row.peakloc.toFixed(4) : row.peakloc}</td>
            <td>${row.pct50?.toFixed ? row.pct50.toFixed(4) : row.pct50}</td>
        `;
        tbody.appendChild(tr);
    });
}

renderTable();

const plotDiv = document.getElementById("plot");

plotDiv.on("plotly_click", function(evt) {
    if (!evt.points || evt.points.length === 0) return;

    const clickedIndex = evt.points[0].curveNumber;
    const clickedName = plotTraces[clickedIndex].name;

    const opacities = [];
    const widths = [];

    for (let i = 0; i < plotTraces.length; i++) {
        if (i === clickedIndex) {
            opacities.push(1.0);
            widths.push(4);
        } else {
            opacities.push(0.12);
            widths.push(1.5);
        }
    }

    Plotly.restyle(plotDiv, {
        opacity: opacities,
        "line.width": widths
    });

    document.getElementById("info").textContent =
        `Selected ${clickedName}. Double-click empty plot space to restore all traces.`;

    renderTable(clickedName);
});

plotDiv.on("plotly_doubleclick", function() {
    Plotly.restyle(plotDiv, {
        opacity: plotTraces.map(() => 0.65),
        "line.width": plotTraces.map(() => 2)
    });

    document.getElementById("info").textContent =
        "Click a curve to inspect that subject. Plotly legend click also hides/shows traces.";

    renderTable();
});
</script>
</body>
</html>
"""

    html = (
        html_template
        .replace("__SKEW_THRESHOLD__", f"{skew_threshold:.2f}")
        .replace("__N_BLUE__", str(n_blue))
        .replace("__N_RED__", str(n_red))
        .replace("__Y_LABEL__", y_label)
        .replace("__TRACES_JSON__", traces_json)
        .replace("__SUMMARY_JSON__", summary_json)
        .replace("__PLOT_XMIN__", f"{plot_xmin:.6f}")
        .replace("__PLOT_XMAX__", f"{plot_xmax:.6f}")
        .replace("__YAXIS_CONFIG__", yaxis_config)
    )

    with open(out_html, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"[INFO] Saved interactive HTML to {out_html}")
    print(f"[INFO] Saved summary CSV to {out_csv}")


def build_png(subjects, out_png, normalize=False, dpi=300):
    fig, ax = plt.subplots(figsize=(6, 4), dpi=dpi)

    n_red = 0
    n_blue = 0

    for subj_id in subjects:
        hist_base = ensure_hist_exists(subj_id)
        x, y, md = load_hist(hist_base)

        skewness = weighted_skewness(x, y)

        y_plot = y.astype(float).copy()
        if normalize:
            total = y_plot.sum()
            if total > 0:
                y_plot = y_plot / total

        if skewness < 0:
            color = "darkblue"
            n_blue += 1
        else:
            color = "red"
            n_red += 1

        ax.plot(x, y_plot, color=color, alpha=0.5, linewidth=0.8)

    ax.axvline(0.5, color="black", linestyle="--", linewidth=0.8, alpha=0.8)
    ax.set_xlabel("Maximum correlation coefficient (r)")
    ax.set_ylabel("Proportion of voxels" if normalize else "Voxel count")
    ax.set_title("Voxelwise maximum-correlation distributions")

    if not normalize:
        ax.yaxis.set_major_locator(MultipleLocator(1000))
        ax.yaxis.set_major_formatter(
            FuncFormatter(lambda v, pos: f"{int(v):d}" if v >= 0 else "")
        )
        ax.set_ylim(0, 5000)

    legend_handles = [
        Line2D([0], [0], color="darkblue", lw=2, label=f"Negative-skew (n={n_blue})"),
        Line2D([0], [0], color="red", lw=2, label=f"Positive-skew (n={n_red})"),
    ]
    ax.legend(handles=legend_handles, loc="best", frameon=False, fontsize=8)

    fig.tight_layout()
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Create interactive HTML of rapidtide maxcorr distributions with skewness-based QC summary."
    )
    parser.add_argument(
        "--subjects",
        nargs="+",
        required=True,
        help="Subject IDs with or without 'sub-' prefix.",
    )
    parser.add_argument(
        "--out-html",
        default="maxcorr_distributions_skewness_n100_defaultfilt_5to15.html",
        help="Output HTML filename.",
    )
    parser.add_argument(
        "--out-csv",
        default="maxcorr_skewness_summary_n100_defaultfilt_5to15.csv",
        help="Output CSV filename.",
    )
    parser.add_argument(
        "--out-png",
        default="maxcorr_distributions_defaultfilt_5to15_100.png",
        help="Output PNG filename.",
    )
    parser.add_argument(
        "--normalize",
        action="store_true",
        help="Normalize each histogram to proportion of voxels instead of raw voxel count.",
    )
    parser.add_argument(
        "--skew-threshold",
        type=float,
        default=0.0,
        help="Skewness threshold for classifying curves. Default 0.0.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="DPI for the output PNG (default 300).",
    )

    args = parser.parse_args()

    build_html(
        args.subjects,
        args.out_html,
        args.out_csv,
        normalize=args.normalize,
        skew_threshold=args.skew_threshold,
    )

    build_png(args.subjects, args.out_png, normalize=args.normalize, dpi=args.dpi)


if __name__ == "__main__":
    main()
