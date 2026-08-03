#!/usr/bin/env python3
"""Aggregate mtpole-ttj-pyconvino results.json files into one structured dataset.

Kept separate from collectCondorScans.py (scan-specific, different schema/key
structure) to avoid the same result-collision bug class that script already
hit once with its scan keys.

Primary key per fit: (exp, theory_source, order, PDF, poi_config, poly_order)
-- the shared axis vocabulary written by fit_object.py's _write_results_json
(see mtpole-ttj-pyconvino/fit_object.py). Two results.json files mapping to
the same key is a hard error by default: silently keeping one (e.g. "last
one found by os.walk") would make the aggregated dataset depend on
filesystem walk order, an easy way to quietly lose or overwrite a real fit
result.

Blinded fits (is_blinded=True) all carry blind_salt_fingerprint -- the
truncated sha256 of the salt file used, never the salt value itself (see
blinding.py). Two blinded fits with different fingerprints used different
salts, which breaks the shift-cancellation guarantee every downstream
comparison plot relies on; that is always a hard error, not something
--allow-duplicates can wave through.
"""

import argparse
import glob
import json
import os
import sys
from datetime import datetime, timezone

import pandas as pd

KEY_FIELDS = ('exp', 'theory_source', 'order', 'PDF', 'poi_config', 'poly_order')


def find_results_json(roots):
    """Yield absolute paths to every results.json under any of `roots`."""
    seen = set()
    for root in roots:
        for path in glob.glob(os.path.join(root, '**', 'results.json'), recursive=True):
            abspath = os.path.abspath(path)
            if abspath not in seen:
                seen.add(abspath)
                yield abspath


def load_result(path):
    with open(path) as fh:
        data = json.load(fh)
    missing = [f for f in KEY_FIELDS if f not in data]
    if missing:
        raise SystemExit(
            f'{path}: results.json missing required field(s) {missing} -- '
            f'produced by a pre-Phase-0 doFit.py run? Rerun the fit with the '
            f'current fit_object.py.')
    return data


def build_key(data):
    return tuple(data[f] for f in KEY_FIELDS)


def collect(roots, allow_duplicates=False):
    """Returns (rows, n_duplicates_resolved)."""
    by_key = {}
    n_dupes = 0
    for path in find_results_json(roots):
        data = load_result(path)
        key = build_key(data)
        mtime = os.path.getmtime(path)
        if key in by_key:
            other_path, other_mtime, _ = by_key[key]
            if not allow_duplicates:
                raise SystemExit(
                    f'Duplicate results for key {key}:\n'
                    f'  {other_path}\n'
                    f'  {path}\n'
                    f'Pass --allow-duplicates to keep the most recently modified '
                    f'one instead of erroring.')
            n_dupes += 1
            keep = (path, mtime, data) if mtime >= other_mtime else by_key[key]
            dropped = path if keep[0] != path else other_path
            print(f'WARNING: duplicate results for key {key}; keeping the newer '
                  f'file, dropping {dropped}', file=sys.stderr)
            by_key[key] = keep
        else:
            by_key[key] = (path, mtime, data)
    return by_key, n_dupes


def check_salt_consistency(by_key):
    """Hard error if any two blinded fits used different salts."""
    fingerprints = {}
    for key, (path, _mtime, data) in by_key.items():
        if data.get('is_blinded'):
            fp = data.get('blind_salt_fingerprint')
            if fp is None:
                raise SystemExit(f'{path}: is_blinded=True but blind_salt_fingerprint is missing.')
            fingerprints.setdefault(fp, []).append(path)
    if len(fingerprints) > 1:
        msg = ['Blinded results were produced with more than one salt -- their '
               'shifts are not comparable and must not be aggregated together:']
        for fp, paths in fingerprints.items():
            msg.append(f'  salt {fp}:')
            for p in paths:
                msg.append(f'    {p}')
        raise SystemExit('\n'.join(msg))
    return next(iter(fingerprints), None)


def to_rows(by_key):
    rows = []
    for key, (path, mtime, data) in by_key.items():
        row = dict(zip(KEY_FIELDS, key))
        for k, v in data.items():
            if k in KEY_FIELDS:
                continue
            row[k] = v
        row['source_path'] = path
        row['source_mtime'] = datetime.fromtimestamp(mtime, tz=timezone.utc).isoformat()
        rows.append(row)
    return rows


def write_outputs(df, out_prefix):
    csv_path = out_prefix + '.csv'
    df.to_csv(csv_path, index=False)
    print(f'Wrote {csv_path} ({len(df)} rows)')


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate mtpole-ttj-pyconvino results.json files into one '
                     'structured dataset (order x PDF x dataset x POI-config matrix).')
    parser.add_argument('--roots', nargs='+', required=True,
                         help='One or more directories to walk for results.json files.')
    parser.add_argument('--out-prefix', default='matrix',
                         help='Output path prefix; writes <prefix>.csv, <prefix>_meta.json.')
    parser.add_argument('--allow-duplicates', action='store_true',
                         help='On a primary-key collision, keep the most recently modified '
                              'results.json (with a warning) instead of erroring.')
    args = parser.parse_args()

    by_key, n_dupes = collect(args.roots, allow_duplicates=args.allow_duplicates)
    if not by_key:
        raise SystemExit(f'No results.json files found under {args.roots}')

    salt_fingerprint = check_salt_consistency(by_key)

    rows = to_rows(by_key)
    df = pd.DataFrame(rows)
    write_outputs(df, args.out_prefix)

    meta = {
        'generated_at': datetime.now(timezone.utc).isoformat(),
        'roots': [os.path.abspath(r) for r in args.roots],
        'fit_count': len(df),
        'duplicates_resolved': n_dupes,
        'blind_salt_fingerprint': salt_fingerprint,
    }
    meta_path = args.out_prefix + '_meta.json'
    with open(meta_path, 'w') as fh:
        json.dump(meta, fh, indent=2)
    print(f'Wrote {meta_path}')


if __name__ == '__main__':
    main()
