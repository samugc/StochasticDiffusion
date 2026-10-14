"""Integrity checks and failure-path tests; no optimizer or historical logs needed."""
import csv
import hashlib
import json
import math
import tempfile
import unittest
from pathlib import Path

import analyze_historical_results as analysis
import extract_historical_results as extraction

ROOT = Path(__file__).resolve().parents[1]


class HistoricalDataTests(unittest.TestCase):
    def setUp(self):
        with (ROOT / 'results/historical_runs.csv').open(newline='') as handle:
            self.rows = list(csv.DictReader(handle))

    def write_runs(self, directory, rows):
        target = Path(directory) / 'runs.csv'
        with target.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=extraction.COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
        return target

    def test_provenance_digest_and_grid(self):
        provenance = json.loads((ROOT / 'results/historical_provenance.json').read_text())
        digest = hashlib.sha256((ROOT / 'results/historical_runs.csv').read_bytes()).hexdigest()
        self.assertEqual(digest, provenance['output']['sha256'])
        extraction.validate_rows(self.rows)
        groups = analysis.load_runs(ROOT / 'results/historical_runs.csv')
        self.assertEqual([v.shape for v in groups.values()], [(30, 7)] * 7)

    def test_duplicate_run_rejected_by_both_validators(self):
        self.rows[-1]['run'] = '29'
        with self.assertRaises(ValueError):
            extraction.validate_rows(self.rows)
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ValueError):
                analysis.load_runs(self.write_runs(directory, self.rows))

    def test_missing_run_rejected(self):
        with self.assertRaises(ValueError):
            extraction.validate_rows(self.rows[:-1])
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ValueError):
                analysis.load_runs(self.write_runs(directory, self.rows[:-1]))

    def test_nan_and_infinity_rejected(self):
        for bad in ('nan', 'inf', '-inf'):
            self.rows[0]['eta'] = bad
            with self.assertRaises(ValueError):
                extraction.validate_rows(self.rows)
            with tempfile.TemporaryDirectory() as directory:
                with self.assertRaises(ValueError):
                    analysis.load_runs(self.write_runs(directory, self.rows))

    def test_legacy_mismatch_above_tolerance_is_fatal(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            for algorithm in extraction.ALGORITHMS:
                with (path / f'errores{algorithm}.csv').open('w', newline='') as handle:
                    writer = csv.writer(handle)
                    writer.writerow(['Algoritmo', 'Ejecución', 'Parámetro', 'Error'])
                    for row in self.rows:
                        if row['algorithm'] == algorithm:
                            for parameter in extraction.PARAMETERS:
                                error = abs(float(row[parameter]) - extraction.REFERENCE[parameter])
                                writer.writerow([algorithm, row['run'], parameter, error])
            self.assertEqual(extraction.validate_legacy(self.rows, path)[0], 0)
            target = path / 'erroresBL.csv'
            with target.open(newline='') as handle:
                records = list(csv.reader(handle))
            records[1][-1] = str(float(records[1][-1]) + 2e-12)
            with target.open('w', newline='') as handle:
                csv.writer(handle).writerows(records)
            with self.assertRaisesRegex(ValueError, 'Historical error mismatch'):
                extraction.validate_legacy(self.rows, path)

    def test_documented_holm_family(self):
        actual = analysis.holm_adjust([.04, .001, .03, .2, .01])
        for result, expected in zip(actual, [.09, .005, .09, .2, .04]):
            self.assertAlmostEqual(float(result), expected)

    def test_summaries_against_independent_standard_library(self):
        import statistics
        groups = analysis.load_runs(ROOT / 'results/historical_runs.csv')
        summaries = analysis.summarize(groups)
        for summary in summaries:
            rows = [r for r in self.rows if r['algorithm'] == summary['algorithm']]
            times = [float(r['runtime_s']) for r in rows]
            for name, value in [('mean_runtime_s', statistics.mean(times)),
                                ('median_runtime_s', statistics.median(times)),
                                ('std_runtime_s', statistics.stdev(times))]:
                self.assertTrue(math.isclose(summary[name], value, rel_tol=1e-12))
            for parameter, reference in extraction.REFERENCE.items():
                absolute = [abs(float(r[parameter]) - reference) for r in rows]
                relative = [v / abs(reference) for v in absolute]
                for prefix, value in [('mean_absolute_error', statistics.mean(absolute)),
                                      ('median_absolute_error', statistics.median(absolute)),
                                      ('mean_relative_absolute_error', statistics.mean(relative)),
                                      ('median_relative_absolute_error', statistics.median(relative))]:
                    self.assertTrue(math.isclose(summary[f'{prefix}_{parameter}'], value, rel_tol=1e-12, abs_tol=1e-18))


if __name__ == '__main__':
    unittest.main()
