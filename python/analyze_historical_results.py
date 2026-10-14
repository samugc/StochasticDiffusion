#!/usr/bin/env python3
"""Analyze only results/historical_runs.csv; never execute an optimizer.

Regenerates summaries, five parameter-specific Kruskal-Wallis tests, and figures.
Paths resolve relative to this script, so the working directory is irrelevant.
"""
import csv
import hashlib
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import scipy
from scipy.stats import kruskal

ROOT = Path(__file__).resolve().parents[1]
ALGORITHMS = ('BL', 'SA', 'ILS', 'EA', 'EA+E', 'MA+E', 'DE')
PARAMETERS = ('eta', 'beta1', 'beta2', 'beta3', 'sigma2')
REFERENCE = np.array([math.exp(-1), .1, -.009, .0002, .0001])
LABELS = (r'$\eta$', r'$\beta_1$', r'$\beta_2$', r'$\beta_3$', r'$\sigma^2$')
COLUMNS = ('algorithm', 'run', *PARAMETERS, 'fitness', 'runtime_s')
COLORS = ('#0072B2', '#E69F00', '#009E73', '#D55E00', '#CC79A7', '#56B4E9', '#6B5B45')


def load_runs(path):
    with path.open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != list(COLUMNS):
            raise ValueError('Unexpected historical_runs.csv schema')
        rows = list(reader)
    if len(rows) != 210:
        raise ValueError(f'Expected 210 runs, found {len(rows)}')
    groups, seen = {a: [] for a in ALGORITHMS}, set()
    for row in rows:
        algorithm, run = row['algorithm'], int(row['run'])
        key = algorithm, run
        if algorithm not in groups or run not in range(1, 31) or key in seen:
            raise ValueError(f'Invalid or duplicate run: {key}')
        if str(run) != row['run']:
            raise ValueError(f'Non-canonical run: {row["run"]}')
        seen.add(key)
        values = [float(row[c]) for c in (*PARAMETERS, 'fitness', 'runtime_s')]
        if not np.isfinite(values).all() or values[-1] <= 0:
            raise ValueError(f'Non-finite values or non-positive runtime: {key}')
        groups[algorithm].append((run, values))
    if seen != {(a, r) for a in ALGORITHMS for r in range(1, 31)}:
        raise ValueError('Each algorithm must have exactly runs 1..30')
    return {a: np.array([v for _, v in sorted(groups[a])]) for a in ALGORITHMS}


def holm_adjust(pvalues):
    """Holm step-down adjusted p-values, restored to original hypothesis order."""
    pvalues = np.asarray(pvalues, dtype=float)
    order = np.argsort(pvalues)
    adjusted = np.minimum(1., np.maximum.accumulate((len(order) - np.arange(len(order))) * pvalues[order]))
    result = np.empty_like(adjusted)
    result[order] = adjusted
    return result


def summarize(groups):
    output = []
    for algorithm, values in groups.items():
        absolute = np.abs(values[:, :5] - REFERENCE)
        relative = absolute / np.abs(REFERENCE)
        runtime = values[:, -1]
        row = dict(algorithm=algorithm, n=len(values), mean_runtime_s=np.mean(runtime),
                   median_runtime_s=np.median(runtime), std_runtime_s=np.std(runtime, ddof=1))
        for index, parameter in enumerate(PARAMETERS):
            row[f'mean_absolute_error_{parameter}'] = np.mean(absolute[:, index])
            row[f'median_absolute_error_{parameter}'] = np.median(absolute[:, index])
            row[f'mean_relative_absolute_error_{parameter}'] = np.mean(relative[:, index])
            row[f'median_relative_absolute_error_{parameter}'] = np.median(relative[:, index])
        output.append(row)
    return output


def statistical_report(groups, input_hash):
    results = []
    for index, parameter in enumerate(PARAMETERS):
        samples = [np.abs(v[:, index] - REFERENCE[index]) / abs(REFERENCE[index]) for v in groups.values()]
        test = kruskal(*samples)
        results.append((parameter, float(test.statistic), float(test.pvalue)))
    adjusted = holm_adjust([r[2] for r in results])
    text = [
        '# Historical statistical analysis', '',
        'Generated exclusively from `historical_runs.csv`; no optimizer was executed.', '',
        f'Input SHA-256: `{input_hash}`.', '',
        '## Outcome and method', '',
        'For each parameter separately, the outcome is relative absolute error ',
        '`abs(estimate - reference) / abs(reference)` (a fraction, not a percentage). ',
        'References: eta = exp(-1), beta1 = 0.1, beta2 = -0.009, beta3 = 0.0002, sigma2 = 0.0001.', '',
        'Kruskal–Wallis compares seven groups of 30 historical runs. SciPy applies tie correction. ',
        'Within each parameter, relative and absolute errors differ by a positive constant; ',
        'their rank tests are therefore equivalent in principle. Results below use relative error explicitly.', '',
        'Holm correction controls family-wise error at 0.05 across the five omnibus hypotheses ',
        '(one per parameter). Raw and adjusted p-values are both reported.', '',
        '| Parameter | n/group | H | Raw p | Holm p (5 tests) | Reject at 0.05 |',
        '|---|---:|---:|---:|---:|:---:|',
    ]
    for (parameter, h, p), corrected in zip(results, adjusted):
        text.append(f'| {parameter} | 30 | {h:.12g} | {p:.12g} | {corrected:.12g} | {"yes" if corrected < .05 else "no"} |')
    text += [
        '', '## Interpretation and limitations', '',
        '- These are omnibus distribution/rank tests, not a ranking or proof that any specific pair differs.',
        '- No pairwise post-hoc tests are performed. No paired Wilcoxon conclusions are published: run IDs do not establish pairing.',
        '- The tests assume independent runs. Historical documentation describes independent executions, but missing seeds and complete logs prevent auditing that assumption.',
        '- All runs use one experimental dataset. Results concern algorithm variability conditional on that dataset, not generalization across datasets.',
        '- Different distribution shapes and variances mean a rejection must not be described simply as a difference in medians.',
        '- Rounded estimates create ties, particularly for sigma2 and DE. p-values describe the recovered precision, not unobserved full-precision output.',
        '- An exact zero error at the reported precision is not evidence of zero uncertainty.',
        '- These are reconstructed historical results, not a benchmark executed in 2026. Statistical significance is not practical superiority.',
        '', '## Numerical environment', '',
        f'NumPy {np.__version__}; SciPy {scipy.__version__}; Matplotlib {matplotlib.__version__}.', '',
    ]
    return '\n'.join(text), results, adjusted


def style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.titleweight': 'bold', 'axes.labelcolor': '#273444',
                         'text.color': '#172B4D', 'axes.edgecolor': '#9BA7B4',
                         'grid.color': '#DDE3EA', 'grid.alpha': .7,
                         'savefig.facecolor': 'white'})


def boxpoints(ax, samples):
    boxes = ax.boxplot(samples, patch_artist=True, widths=.55, showfliers=False,
                       medianprops={'color': '#111111', 'linewidth': 1.7},
                       whiskerprops={'color': '#7A8793'}, capprops={'color': '#7A8793'})
    jitter = np.random.default_rng(2026).uniform(-.14, .14, size=30)
    for index, (box, sample) in enumerate(zip(boxes['boxes'], samples)):
        box.set(facecolor=COLORS[index], alpha=.28, edgecolor=COLORS[index])
        ax.scatter(index + 1 + jitter, sample, color=COLORS[index], s=14, alpha=.7, linewidths=0, zorder=3)
    ax.set_xticks(range(1, 8), ALGORITHMS)
    ax.grid(axis='y')
    ax.set_axisbelow(True)


def figures(groups, target):
    style()
    relative = {a: 100 * np.abs(v[:, :5] - REFERENCE) / np.abs(REFERENCE) for a, v in groups.items()}
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), layout='constrained')
    fig.suptitle('Historical parameter recovery · 30 runs per algorithm', fontsize=17)
    for index, ax in enumerate(axes.flat):
        if index == 5:
            ax.axis('off')
            ax.text(.03, .86, 'Each dot is one historical run.\nBoxes: median and interquartile range.\nWhiskers: 1.5 × IQR; all points shown.\n\nRelative error (%) =\n100 × |estimate − reference| / |reference|\n\nSeparate y scales for each parameter.\nσ² uses symlog: linear below 0.01%,\nlogarithmic above (zeros preserved).\n\nOne dataset · no new optimizer runs.', va='top', fontsize=10, linespacing=1.5)
            continue
        boxpoints(ax, [relative[a][:, index] for a in ALGORITHMS])
        ax.set_title(LABELS[index], fontsize=15)
        ax.set_ylabel('Relative absolute error (%)')
        if index == 4:
            ax.set_yscale('symlog', linthresh=.01)
            ax.set_ylim(-.002, max(relative[a][:, index].max() for a in ALGORITHMS) * 2)
        else:
            ax.set_ylim(bottom=-.02 * max(relative[a][:, index].max() for a in ALGORITHMS))
    fig.savefig(target / 'historical_relative_errors.png', dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5.8), layout='constrained')
    boxpoints(ax, [groups[a][:, -1] for a in ALGORITHMS])
    ax.set_title('Historical runtime · 30 runs per algorithm', fontsize=17, pad=18)
    ax.set_ylabel('Runtime (seconds)')
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Dots: individual runs · boxes: median / IQR · whiskers: 1.5 × IQR\nHistorical measurements on one dataset; no new benchmark executed')
    for index, a in enumerate(ALGORITHMS):
        ax.text(index + 1, np.max(groups[a][:, -1]) + 95,
                f'med {np.median(groups[a][:, -1]):,.1f}s', ha='center', fontsize=9)
    ax.margins(y=.12)
    fig.savefig(target / 'historical_runtime.png', dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), layout='constrained')
    fig.suptitle('Historical runtime vs parameter error · no composite accuracy metric', fontsize=16)
    for index, ax in enumerate(axes.flat):
        if index == 5:
            ax.axis('off')
            handles = [Line2D([], [], marker='o', linestyle='', color=c, label=a) for a, c in zip(ALGORITHMS, COLORS)]
            ax.legend(handles=handles, loc='upper left', ncol=2, frameon=False)
            ax.text(.03, .48, 'Each dot pairs runtime and error from\nthe same documented execution.\n30 runs per algorithm; one dataset.\n\nNo averaging across parameters.\nRuntime: logarithmic scale.\nσ² error: symlog (zeros preserved).\nOther error axes: linear.', va='top', fontsize=11, linespacing=1.4)
            continue
        for algorithm, color in zip(ALGORITHMS, COLORS):
            ax.scatter(groups[algorithm][:, -1], relative[algorithm][:, index],
                       color=color, s=20, alpha=.7, linewidths=.3, edgecolors='white')
        ax.set_title(LABELS[index], fontsize=15)
        ax.set_xlabel('Runtime (seconds; log scale)')
        ax.set_xscale('log')
        ax.set_ylabel('Relative absolute error (%)')
        if index == 4:
            ax.set_yscale('symlog', linthresh=.01)
            ax.set_ylim(-.002, max(relative[a][:, index].max() for a in ALGORITHMS) * 2)
        else:
            ax.set_ylim(bottom=-.02 * max(relative[a][:, index].max() for a in ALGORITHMS))
        ax.grid(alpha=.5)
        ax.set_axisbelow(True)
    fig.savefig(target / 'historical_tradeoff.png', dpi=180)
    plt.close(fig)


def main():
    source = ROOT / 'results/historical_runs.csv'
    groups = load_runs(source)
    rows = summarize(groups)
    report, tests, adjusted = statistical_report(groups, hashlib.sha256(source.read_bytes()).hexdigest())
    with (ROOT / 'results/historical_summary.csv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    (ROOT / 'results/statistical_tests.md').write_text(report, encoding='utf-8')
    target = ROOT / 'figures'
    target.mkdir(exist_ok=True)
    figures(groups, target)
    print('algorithm,n,mean_runtime_s,median_runtime_s,std_runtime_s')
    for row in rows:
        print(','.join(str(row[key]) for key in ('algorithm', 'n', 'mean_runtime_s', 'median_runtime_s', 'std_runtime_s')))
    for (parameter, h, p), corrected in zip(tests, adjusted):
        print(f'{parameter}: H={h:.12g}; p={p:.12g}; Holm={corrected:.12g}')


if __name__ == '__main__':
    main()
