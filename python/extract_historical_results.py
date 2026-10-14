#!/usr/bin/env python3
"""Extract documentary runs, preserving numeric tokens; validate legacy errors.

No historical algorithm is executed. Only the seven labelled individual-run
longtables are accepted. All checks complete before any output is written.
"""
import argparse
import csv
import hashlib
import io
import json
import math
import re
from pathlib import Path

ALGORITHMS = ('BL', 'SA', 'ILS', 'EA', 'EA+E', 'MA+E', 'DE')
PARAMETERS = ('eta', 'beta1', 'beta2', 'beta3', 'sigma2')
REFERENCE = dict(zip(PARAMETERS, (math.exp(-1), .1, -.009, .0002, .0001)))
COLUMNS = ('algorithm', 'run', *PARAMETERS, 'fitness', 'runtime_s')
LABELS = dict(zip(('bl_30runs', 'sa_30runs', 'ils_30runs', 'ga_30runs',
                   'ga_elite_30runs', 'memetic_elite_30runs', 'de_30runs'), ALGORITHMS))
TOLERANCE = 1e-12


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_rows(rows):
    if len(rows) != 210:
        raise ValueError(f'Expected 210 runs, found {len(rows)}')
    seen = set()
    for row in rows:
        if set(row) != set(COLUMNS) or row['algorithm'] not in ALGORITHMS:
            raise ValueError(f'Invalid columns or algorithm: {row}')
        run = int(row['run'])
        if str(run) != row['run'] or run not in range(1, 31):
            raise ValueError(f'Invalid run: {row}')
        key = (row['algorithm'], run)
        if key in seen:
            raise ValueError(f'Duplicate run: {key}')
        seen.add(key)
        for column in (*PARAMETERS, 'fitness', 'runtime_s'):
            if not math.isfinite(float(row[column])):
                raise ValueError(f'Non-finite {column}: {key}')
        if float(row['runtime_s']) <= 0:
            raise ValueError(f'Non-positive runtime: {key}')
    expected = {(a, r) for a in ALGORITHMS for r in range(1, 31)}
    if seen != expected:
        raise ValueError(f'Missing or unexpected runs: {seen ^ expected}')


def extract(source):
    text = source.read_text(encoding='utf-8')
    rows, locations, found = [], [], set()
    blocks = list(re.finditer(r'\\begin\{longtable\}.*?\\end\{longtable\}', text, re.S))
    if len(blocks) != 7:
        raise ValueError(f'Expected seven longtables, found {len(blocks)}')
    for block in blocks:
        labels = re.findall(r'\\label\{tab:([^}]+)\}', block.group())
        matches = [label for label in labels if label in LABELS]
        if len(matches) != 1 or matches[0] in found:
            raise ValueError(f'Ambiguous or duplicate table labels: {labels}')
        label = matches[0]
        found.add(label)
        algorithm = LABELS[label]
        start = text.count('\n', 0, block.start()) + 1
        table_count = 0
        for offset, line in enumerate(block.group().splitlines()):
            if not re.match(r'\s*\d+\s*&', line):
                continue
            tokens = [s.strip() for s in line.split('\\\\', 1)[0].split('&')]
            if len(tokens) != 8:
                raise ValueError(f'Expected 8 numeric cells at line {start + offset}')
            row = dict(zip(COLUMNS, (algorithm, *tokens)))
            rows.append(row)
            locations.append({'algorithm': algorithm, 'run': int(tokens[0]),
                              'table_label': 'tab:' + label, 'source_line': start + offset})
            table_count += 1
        if table_count != 30:
            raise ValueError(f'{algorithm}: expected 30 rows, found {table_count}')
    rows.sort(key=lambda r: (ALGORITHMS.index(r['algorithm']), int(r['run'])))
    validate_rows(rows)
    return rows, locations


def validate_legacy(rows, errors_dir):
    lookup = {(r['algorithm'], int(r['run'])): r for r in rows}
    max_difference = 0.0
    files = []
    for algorithm in ALGORITHMS:
        path = errors_dir / f'errores{algorithm}.csv'
        seen = set()
        with path.open(encoding='utf-8', newline='') as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != ['Algoritmo', 'Ejecución', 'Parámetro', 'Error']:
                raise ValueError(f'Unexpected historical CSV schema: {path}')
            for record in reader:
                run, parameter = int(record['Ejecución']), record['Parámetro']
                key = (run, parameter)
                if record['Algoritmo'] != algorithm or key in seen:
                    raise ValueError(f'Wrong algorithm or duplicate historical error: {path}, {key}')
                if run not in range(1, 31) or parameter not in PARAMETERS:
                    raise ValueError(f'Unexpected historical key: {path}, {key}')
                seen.add(key)
                observed = float(record['Error'])
                reconstructed = abs(float(lookup[algorithm, run][parameter]) - REFERENCE[parameter])
                difference = abs(observed - reconstructed)
                if not math.isfinite(observed) or observed < 0 or difference > TOLERANCE:
                    raise ValueError(f'Historical error mismatch: {algorithm} run {run} {parameter}; '
                                     f'CSV={observed:.17g}, reconstructed={reconstructed:.17g}')
                max_difference = max(max_difference, difference)
        if seen != {(r, p) for r in range(1, 31) for p in PARAMETERS}:
            raise ValueError(f'Incomplete historical errors: {path}')
        files.append({'path': str(path.resolve()), 'sha256': digest(path), 'n_errors': len(seen)})
    return max_difference, files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--errors-dir', required=True, type=Path)
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'results')
    parser.add_argument('--validate-only', action='store_true', help='Check existing CSV against source and legacy errors; write nothing')
    args = parser.parse_args()
    rows, locations = extract(args.source)
    maximum, files = validate_legacy(rows, args.errors_dir)
    target = args.output_dir / 'historical_runs.csv'
    buffer = io.StringIO(newline='')
    writer = csv.DictWriter(buffer, fieldnames=COLUMNS, lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    content = buffer.getvalue().encode('utf-8')
    if args.validate_only:
        if target.read_bytes() != content:
            raise ValueError('Existing historical_runs.csv differs from exact source extraction')
    else:
        provenance = {
            'source': {'path': str(args.source.resolve()), 'sha256': digest(args.source)},
            'historical_error_sources': files,
            'output': {'path': 'results/historical_runs.csv', 'sha256': hashlib.sha256(content).hexdigest()},
            'reference_parameters': REFERENCE,
            'validation': {'runs': 210, 'runs_per_algorithm': 30, 'errors_checked': 1050,
                           'absolute_tolerance': TOLERANCE, 'max_absolute_difference': maximum},
            'numeric_policy': 'Original LaTeX numeric tokens preserved; no extra precision inferred.',
            'rows': locations,
        }
        args.output_dir.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)
        (args.output_dir / 'historical_provenance.json').write_text(
            json.dumps(provenance, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    print(f'PASS: 210 runs; 30 per algorithm; 1,050 historical errors; max difference {maximum:.17g} <= {TOLERANCE:g}')


if __name__ == '__main__':
    main()
