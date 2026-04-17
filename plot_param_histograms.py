import numpy as np
import matplotlib.pyplot as plt
import os
import csv
from pathlib import Path

output_dir = Path('images')
os.makedirs(output_dir, exist_ok=True)

# Helper function to fix Fortran float format (missing 'E' in exponent)
def fix_fortran_float(s):
    import re
    # Insert 'E' before exponent if missing
    return float(re.sub(r'(?<![Ee])([+-]\d+)$', r'E\1', s))

# Load the binning output file
data = []
with open('output/param_binning_output.dat', 'r') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        if '********' in line:
            continue  # skip lines with overflow values
        parts = line.split()
        param = int(parts[0])
        bin_idx = int(parts[1])
        bin_left = fix_fortran_float(parts[2])
        bin_right = fix_fortran_float(parts[3])
        loglik_sum = fix_fortran_float(parts[4])
        data.append((param, bin_idx, bin_left, bin_right, loglik_sum))

data = np.array(data)
if data.size == 0:
    print('No valid data found in param_binning_output.dat. Check the file for valid numeric entries.')
    exit(1)
if data.ndim == 1:
    data = data.reshape(1, -1)

nparams = int(np.max(data[:,0]))
nparam_bins = int(np.max(data[:,1]))

# Build parameter_histograms for CSV and plotting
parameter_histograms = []
for p in range(1, nparams+1):
    param_data = data[data[:,0] == p]
    bin_edges = np.concatenate((param_data[:,2], [param_data[-1,3]]))
    bin_centers = 0.5 * (param_data[:,2] + param_data[:,3])
    accum = param_data[:,4]
    parameter_histograms.append({
        "parameter": f"param_{p}",
        "bin_edges": bin_edges,
        "bin_centers": bin_centers,
        "accumulated_loglikelihood": accum
    })


# Debug: Print bin stats for each parameter after shifting
for idx, payload in enumerate(parameter_histograms):
    edges = np.asarray(payload["bin_edges"], dtype=float)
    logliks = np.asarray(payload["accumulated_loglikelihood"], dtype=float)
    nonzero_bins = np.sum(logliks > 0)
    print(f"Parameter {idx+1}: nonzero bins = {nonzero_bins}/{len(logliks)} (after shift)")
    print(f"  loglikelihood: min={logliks.min()}, max={logliks.max()}, mean={logliks.mean()}")
    print(f"  bin edges: min={edges.min()}, max={edges.max()}")

# Find the global maximum log-likelihood across all bins for normalization
all_logliks = np.array([payload["accumulated_loglikelihood"] for payload in parameter_histograms]).flatten()
max_loglik = np.max(all_logliks)

# Write CSV (with per-parameter exp-normalized values)
histogram_csv_path = output_dir / "parameter_loglikelihood_histograms.csv"
with histogram_csv_path.open("w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["parameter", "bin_index", "bin_left", "bin_right", "bin_center", "exp_normalized_loglikelihood"])
    for payload in parameter_histograms:
        edges = np.asarray(payload["bin_edges"], dtype=float)
        centers = np.asarray(payload["bin_centers"], dtype=float)
        values = np.asarray(payload["accumulated_loglikelihood"], dtype=float)
        for idx in range(len(centers)):
            writer.writerow([
                payload["parameter"],
                idx,
                float(edges[idx]),
                float(edges[idx + 1]),
                float(centers[idx]),
                float(values[idx]),
            ])

# Plot parameter histograms in a grid (per-parameter exp-normalized)
histogram_plot_path = output_dir / "parameter_loglikelihood_histograms.png"
n_parameters = len(parameter_histograms)
if n_parameters > 0:
    import math
    ncols = min(4, n_parameters)
    nrows = int(math.ceil(n_parameters / ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4.5 * ncols, 3.5 * nrows), squeeze=False)
    for idx, payload in enumerate(parameter_histograms):
        ax = axes[idx // ncols][idx % ncols]
        edges = np.asarray(payload["bin_edges"], dtype=float)
        values = np.asarray(payload["accumulated_loglikelihood"], dtype=float)
        widths = np.diff(edges)
        ax.bar(edges[:-1], values, width=widths, align="edge", color="slateblue", edgecolor="black", linewidth=0.2)
        ax.set_title(payload["parameter"])
        ax.set_xlabel("parameter value bin")
        ax.set_ylabel("accumulated value")
    for idx in range(n_parameters, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")
    fig.tight_layout()
    fig.savefig(histogram_plot_path, dpi=150)
    plt.close(fig)

# Build 2D parameter histograms for corner plot (per-pair exp-normalized)
parameter_histograms_2d = []
for i in range(nparams):
    for j in range(nparams):
        if i <= j:
            continue
        param_i_data = data[data[:,0] == (i+1)]
        param_j_data = data[data[:,0] == (j+1)]
        if len(param_i_data) == 0 or len(param_j_data) == 0:
            continue
        xi = 0.5 * (param_i_data[:,2] + param_i_data[:,3])
        yj = 0.5 * (param_j_data[:,2] + param_j_data[:,3])
        values_i = param_i_data[:,4]
        values_j = param_j_data[:,4]
        value_grid = values_i[:,None] + values_j[None,:]
        parameter_histograms_2d.append({
            "param_i": i,
            "param_j": j,
            "accumulated_2d": value_grid
        })

# Plot corner plot (per-pair exp-normalized)
if parameter_histograms_2d and len(parameter_histograms_2d) > 0:
    max_param = max(max(p["param_i"], p["param_j"]) for p in parameter_histograms_2d)
    n_parameters = max_param + 1
    corner_plot_path = output_dir / "parameter_corner_plot.png"
    fig, axes = plt.subplots(nrows=n_parameters, ncols=n_parameters, 
                             figsize=(2.5 * n_parameters, 2.5 * n_parameters), squeeze=False)
    # Precompute bin edges and 1D histograms for all parameters
    param_bin_edges = []
    param_hist_log10 = []
    for p in range(1, n_parameters+1):
        param_data = data[data[:,0] == p]
        bin_edges = np.concatenate((param_data[:,2], [param_data[-1,3]]))
        param_bin_edges.append(bin_edges)
        logliks = param_data[:,4]
        param_hist_log10.append(logliks)
    # Plot lower triangle (i > j): 2D heatmaps (log10)
    for payload in parameter_histograms_2d:
        i = payload["param_i"]
        j = payload["param_j"]
        ax = axes[i, j]
        accumulated_2d = np.asarray(payload["accumulated_2d"], dtype=float)
        xedges = param_bin_edges[i]
        yedges = param_bin_edges[j]
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        im = ax.imshow(accumulated_2d.T, origin="lower", aspect="auto", cmap="viridis", interpolation="nearest", extent=extent)
        ax.set_xlabel(f"param_{i+1}")
        ax.set_ylabel(f"param_{j+1}")
        ax.set_title(f"param_{i+1} vs param_{j+1}")
    # Plot diagonal (i == j): 1D histograms (log10)
    for i in range(n_parameters):
        ax = axes[i, i]
        edges = param_bin_edges[i]
        values = param_hist_log10[i]
        widths = np.diff(edges)
        ax.bar(edges[:-1], values, width=widths, align="edge", color="slateblue", edgecolor="black", linewidth=0.2)
        ax.set_xlabel(f"param_{i+1}")
        ax.set_ylabel("accumulated value")
        ax.set_title(f"param_{i+1}")
    # Hide upper triangle (i < j)
    for i in range(n_parameters):
        for j in range(n_parameters):
            if i < j:
                axes[i, j].axis("off")
    fig.tight_layout()
    fig.savefig(corner_plot_path, dpi=150)
    plt.close(fig)
    print(f"Wrote {corner_plot_path}")

print(f"Wrote {histogram_csv_path}")
print(f"Wrote {histogram_plot_path}")
