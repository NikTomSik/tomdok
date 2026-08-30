#!/usr/bin/env python3
"""Live monitor for the primary GNINA screen (stage 03).

Re-reads results/<target>/all_docking_results.csv while the screen runs and
renders ONE table: TOP-N by affinity among rows with CNN score > threshold,
showing drug name and clinical phase (4 = approved, 3 = Phase III).

Usage:
  python scripts/monitor_gnina.py --target 6vei
  python scripts/monitor_gnina.py --target 6vei --interval 15 --top 15 --cnn 0.7
  python scripts/monitor_gnina.py --target 6vei --once     # single render
"""
import argparse, csv, sys, time
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]


def phase_label(ph):
    if ph == 4: return 'approved'
    if ph == 3: return 'PhIII'
    if ph is None: return '?'
    return f'Ph{ph}'


def load_library():
    """chembl_id -> (pref_name, max_phase). Loaded once at start."""
    lib = {}
    p = BASE / 'data' / 'library' / 'library.tsv'
    if not p.exists():
        return lib
    with open(p, newline='', errors='replace') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            cid = (r.get('chembl_id') or '').strip()
            if not cid:
                continue
            try:
                ph = int(float(r.get('max_phase') or 0))
            except ValueError:
                ph = None
            lib[cid] = ((r.get('pref_name') or '').strip(), ph)
    return lib


def find_col(fields, keys):
    for f in fields or []:
        lf = f.lower()
        if any(k in lf for k in keys):
            return f
    return None


def read_rows(csv_path):
    """(chembl, name, affinity, cnn) rows; partial last line skipped."""
    rows = []
    try:
        with open(csv_path, newline='', errors='replace') as fh:
            rd = csv.DictReader(fh)
            c_id   = find_col(rd.fieldnames, ['chembl'])
            c_name = find_col(rd.fieldnames, ['name', 'pref'])
            c_aff  = find_col(rd.fieldnames, ['aff'])
            c_cnn  = find_col(rd.fieldnames, ['cnn'])
            if not (c_aff and c_cnn):
                return rows
            for r in rd:
                try:
                    aff = float(r[c_aff]); cnn = float(r[c_cnn])
                except (TypeError, ValueError):
                    continue
                rows.append(((r.get(c_id) or '?').strip(),
                             (r.get(c_name) or '').strip(), aff, cnn))
    except (OSError, csv.Error):
        pass
    return rows


def render(rows, lib, a, total):
    hit = [r for r in rows if r[3] > a.cnn]
    hit.sort(key=lambda r: r[2])          # most negative affinity first
    top = hit[:a.top]

    if sys.stdout.isatty() and not a.once:
        sys.stdout.write('\033[2J\033[H')  # clear screen (live mode)

    head = (f"[monitor] target={a.target}  processed={len(rows)}"
            + (f"/{total}" if total else "")
            + f"  cnn>{a.cnn}: {len(hit)}  top={len(top)}  "
            + time.strftime('%H:%M:%S'))
    print(head)
    hdr = (f"{'rank':>4}  {'CHEMBL':<14}  {'name':<18}  "
           f"{'phase':<9}  {'affinity':>8}  {'CNN':>6}")
    print(hdr)
    print('-' * len(hdr))
    if not top:
        print('   (no molecules above CNN threshold yet)')
    for i, (cid, name, aff, cnn) in enumerate(top, 1):
        lname, ph = lib.get(cid, (name, None))
        if not lname:
            lname = name or '?'
        print(f"{i:>4}  {cid:<14}  {lname[:18]:<18}  "
              f"{phase_label(ph):<9}  {aff:>8.2f}  {cnn:>6.3f}")
    print()


def main():
    ap = argparse.ArgumentParser(description='Live top-N monitor for stage 03')
    ap.add_argument('--target', default='7gqu')
    ap.add_argument('--csv', default=None, help='override CSV path')
    ap.add_argument('--interval', type=float, default=30.0, help='refresh, s')
    ap.add_argument('--cnn', type=float, default=0.7, help='CNN threshold (strictly >)')
    ap.add_argument('--top', type=int, default=15)
    ap.add_argument('--once', action='store_true', help='render once and exit')
    a = ap.parse_args()

    csv_path = Path(a.csv) if a.csv else BASE / 'results' / a.target / 'all_docking_results.csv'
    lib = load_library()
    n_pdbqt = len(list((BASE / 'pdbqt').glob('*.pdbqt')))
    total = n_pdbqt or len(lib) or None

    while True:
        render(read_rows(csv_path), lib, a, total)
        if a.once:
            break
        if (BASE / 'results' / a.target / 'summary.json').exists():
            print('[monitor] stage 03 finished - exiting.')
            break
        time.sleep(a.interval)


if __name__ == '__main__':
    main()