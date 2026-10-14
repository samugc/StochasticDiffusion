# Reconstructed historical TFG results

These are **historical documentary results**, reconstructed from the individual
tables of the undergraduate thesis, not a new benchmark executed in 2026.
The dataset contains **7 algorithms × 30 runs = 210 rows**. All runs concern one
experimental dataset; seeds and complete raw logs for all executions are missing.
The row numbers identify executions within an algorithm and **do not establish
paired observations across algorithms**.

## Provenance and extraction

Primary source:

```text
/Users/samu/Documents/Uni/5º/2º Cuatri/TFG/TFG_trabajo/Plantilla-TFG-master/capitulos/aplicacion_info.tex
```

The extractor selects seven explicitly labelled individual-run `longtable`
environments. Algorithm names are normalized to BL, SA, ILS, EA, EA+E, MA+E and DE.
MA+E denotes the memetic algorithm with elitism, also called MA in some historical
figures. Original numeric text, including trailing zeros and scientific notation,
is retained in `historical_runs.csv`. No precision beyond the thesis is inferred.

`historical_provenance.json` records the source SHA-256, the seven historical error
CSV paths and SHA-256 values, the extracted CSV hash, the table label and source
line for every row, references and validation results. The hashes are computed
from the actual source files, not copied unverified from a prior audit.

Historical error CSV directory:

```text
/Users/samu/Documents/Uni/5º/2º Cuatri/TFG/TFG_trabajo/Errores/
```

The validator requires exactly 30 runs numbered 1–30 for each algorithm, no
duplicates, 210 total rows, finite numeric values and positive runtimes. It checks
all **1,050 absolute errors** against `erroresBL.csv`, `erroresSA.csv`,
`erroresILS.csv`, `erroresEA.csv`, `erroresEA+E.csv`, `erroresMA+E.csv` and
`erroresDE.csv`. A difference **greater than 1e-12 fails extraction** before any
output is written. All 1,050 errors agree exactly in the validated floating-point
calculation (maximum discrepancy 0). Agreement between transcribed sources does
not independently verify the original optimizer executions.

The aggregate tables in `conclusion.tex` are **not inputs or numerical authority**.
Summaries are regenerated exclusively from the 210 individual rows. No original
thesis PDF, third-party paper, book or old defense figure is included here.

## Schema and metric definitions

`historical_runs.csv`: `algorithm,run,eta,beta1,beta2,beta3,sigma2,fitness,runtime_s`.
Runtime is in seconds. **sigma2 is the diffusion variance σ², not σ**.
Fitness is the recorded objective/log-likelihood-related value being maximized;
it is **not parameter error or RAE**. Its precision is the precision printed in
the original tables. This analysis does not recompute fitness.

| Parameter | Reference |
|---|---:|
| eta | exp(-1) |
| beta1 | 0.1 |
| beta2 | -0.009 |
| beta3 | 0.0002 |
| sigma2 | 0.0001 |

For each parameter independently:

```text
absolute_error = abs(estimate - reference)
relative_absolute_error = absolute_error / abs(reference)
```

`historical_summary.csv` contains `n`, runtime mean/median/sample standard deviation
(`ddof=1`), and four columns per parameter: mean and median absolute error, mean
and median relative absolute error. Relative errors in the CSV are **fractions**;
multiply by 100 for percentages. No aggregate “accuracy” combines different
parameters. Mean absolute error is not the absolute error of the mean estimate.
These errors are not the historical curve-fitting quantity called RAE.

## Generated figures

- [`historical_relative_errors.png`](../figures/historical_relative_errors.png):
  separate panels for all five parameters, 30 observations per algorithm,
  median/IQR boxes and all points. Independent y scales; σ² uses a labelled
  symmetric-log scale to retain zero errors and large EA errors.
- [`historical_runtime.png`](../figures/historical_runtime.png): all runtime
  observations, median/IQR and median annotations, in seconds.
- [`historical_tradeoff.png`](../figures/historical_tradeoff.png): runtime against
  relative error of each individual parameter, using the same execution row.
  Five panels, log runtime axis, no composite metric or claimed global ranking.

Jitter only spreads dots visually and uses a fixed plotting seed; it is not an
optimizer seed or a change to the results. Boxes use 1.5 × IQR whiskers and every
observation, including outliers, is displayed. The figures describe estimation
and runtime distributions, not fitness convergence during optimization.

## Statistics and limitations

[`statistical_tests.md`](statistical_tests.md) reports Kruskal–Wallis on relative
absolute errors for each parameter. Holm correction covers the family of five
omnibus tests. No pairwise post-hoc tests or paired Wilcoxon claims are made.
Rejections establish distribution/rank differences, not which algorithm is best
or which pairs differ. Independence of historical runs cannot be fully audited
without seeds and logs. Rounded values create ties; one dataset cannot support
claims of general superiority across problems.

Known historical inconsistencies, retained as documentation rather than silently
changing the original thesis:

- Recomputed runtime means differ from old aggregate summaries. For example, ILS
  is 751.3662 s rather than 707.18 s. EA has a 4301.21 s observation, above the
  3964.60 s maximum reported in one old summary.
- SA's mean absolute eta error is 0.00818734314094, not 0.025697 from the old
  conclusion table. Some historical text mixes error of the mean and mean error.
- Constant/rounded sigma2 values produce ties and zero errors for some groups.
  They do not imply zero estimation uncertainty.
- Separate partial BL and Newton–Raphson logs found during the audit belong to
  other configurations or series and are not merged into these 210 rows.
- Historical RAE figures use a different, curve-based definition and are not
  reproduced or relabelled as parameter accuracy here.

## Reproduce the analysis

From the repository root, with Python 3.13:

```bash
python3 -m venv /tmp/stochastic-diffusion-historical-venv
/tmp/stochastic-diffusion-historical-venv/bin/python -m pip install -r python/requirements-historical.txt
MPLCONFIGDIR=/tmp/stochastic-diffusion-mpl /tmp/stochastic-diffusion-historical-venv/bin/python -B python/analyze_historical_results.py
MPLCONFIGDIR=/tmp/stochastic-diffusion-mpl /tmp/stochastic-diffusion-historical-venv/bin/python -B -m unittest discover -s python -p 'test_historical_results.py' -v
```

The analysis script reads **only `results/historical_runs.csv` as its data input**.
It does not require the original Mac directories, provenance JSON, old summaries,
historical CSVs, C++ programs or network access. It resolves paths relative to its
own location and regenerates the summary, statistical report and three figures.
Pinned top-level versions are in `python/requirements-historical.txt`; fonts use
Matplotlib's DejaVu Sans. Rendering may vary across platforms/dependency builds.

Optional re-extraction on the original Mac (Python standard library only):

```bash
python3 -B python/extract_historical_results.py \
  --source '/Users/samu/Documents/Uni/5º/2º Cuatri/TFG/TFG_trabajo/Plantilla-TFG-master/capitulos/aplicacion_info.tex' \
  --errors-dir '/Users/samu/Documents/Uni/5º/2º Cuatri/TFG/TFG_trabajo/Errores'
```

To compare the existing CSV against the original source and historical error
files **without writing anything**, append `--validate-only` to that command.
The extraction script writes only the historical CSV and its provenance JSON.
Neither script runs an optimizer or changes C++ logic.
