# Historical statistical analysis

Generated exclusively from `historical_runs.csv`; no optimizer was executed.

Input SHA-256: `4c0351472292f576509ae3d81d6664380029a810cbbc84d19b8d01a6bcdf22bb`.

## Outcome and method

For each parameter separately, the outcome is relative absolute error
`abs(estimate - reference) / abs(reference)` (a fraction, not a percentage).
References: eta = exp(-1), beta1 = 0.1, beta2 = -0.009, beta3 = 0.0002, sigma2 = 0.0001.

Kruskal–Wallis compares seven groups of 30 historical runs. SciPy applies tie correction.
Within each parameter, relative and absolute errors differ by a positive constant;
their rank tests are therefore equivalent in principle. Results below use relative error explicitly.

Holm correction controls family-wise error at 0.05 across the five omnibus hypotheses
(one per parameter). Raw and adjusted p-values are both reported.

| Parameter | n/group | H | Raw p | Holm p (5 tests) | Reject at 0.05 |
|---|---:|---:|---:|---:|:---:|
| eta | 30 | 152.655188106 | 2.12353151065e-30 | 4.2470630213e-30 | yes |
| beta1 | 30 | 162.369955076 | 1.86393782323e-32 | 7.45575129292e-32 | yes |
| beta2 | 30 | 159.228400942 | 8.62684478019e-32 | 2.58805343406e-31 | yes |
| beta3 | 30 | 169.855126385 | 4.82749414139e-34 | 2.4137470707e-33 | yes |
| sigma2 | 30 | 133.264380941 | 2.63863433119e-26 | 2.63863433119e-26 | yes |

## Interpretation and limitations

- These are omnibus distribution/rank tests, not a ranking or proof that any specific pair differs.
- No pairwise post-hoc tests are performed. No paired Wilcoxon conclusions are published: run IDs do not establish pairing.
- The tests assume independent runs. Historical documentation describes independent executions, but missing seeds and complete logs prevent auditing that assumption.
- All runs use one experimental dataset. Results concern algorithm variability conditional on that dataset, not generalization across datasets.
- Different distribution shapes and variances mean a rejection must not be described simply as a difference in medians.
- Rounded estimates create ties, particularly for sigma2 and DE. p-values describe the recovered precision, not unobserved full-precision output.
- An exact zero error at the reported precision is not evidence of zero uncertainty.
- These are reconstructed historical results, not a benchmark executed in 2026. Statistical significance is not practical superiority.

## Numerical environment

NumPy 2.2.3; SciPy 1.15.2; Matplotlib 3.10.1.
