"""GNINA virtual screening on 2 GPUs (final version + ETA).
Stage 3 of the pipeline: reads receptor.pdbqt (with Gasteiger charges) and
docking box from receptor.json (clipped to min_box=25 A), screens the
library from data/library/library.tsv and writes:
  results/<target>/all_docking_results.csv  (affinity + cnn_score)
  results/<target>/summary.json             (for Methods)
  results/<target>/poses/                   (3D SDF poses)

Optional: --vina-top N performs Vina cross-validation of top-N hits.

Resume: on restart, skip already-successful entries; permanent-skip after 3 failures.
ETA: live estimate via average_time * remaining.

Usage:
 python scripts/03_screen_gnina.py 7gqu
 python scripts/03_screen_gnina.py 4cgz --limit 50 --gpus 0,1
 python scripts/03_screen_gnina.py 7gqu --vina-top 15 --only hit_list.txt
"""
import argparse, csv, json, multiprocessing, os, shutil, subprocess, time, sys
from pathlib import Path
from queue import Empty

MAX_FAILS = 3
G, NAMES = {}, {}


def is_float(s):
    try: float(s); return True
    except (ValueError, TypeError): return False


def fmt_eta(sec):
    sec = max(0, int(sec))
    if sec < 60: return f'{sec} s'
    if sec < 3600: return f'{sec // 60} min {sec % 60} s'
    return f'{sec // 3600} h {(sec % 3600) // 60} min'


def load_names(paths):
    """Drug names from library.tsv (single source)."""
    names = {}
    for p in paths:
        p = Path(p)
        if not p.exists(): continue
        with open(p, encoding='utf-8') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                cid = (row.get('chembl_id') or '').strip()
                nm = (row.get('pref_name') or '').strip()
                if cid and nm: names.setdefault(cid, nm)
    return names


def parse_gnina(out):
    best = None
    for line in out.splitlines():
        p = line.split()
        if len(p) >= 4 and p[0].isdigit():
            try: aff, cnn = float(p[1]), float(p[3])
            except ValueError: continue
            if not (0.0 <= cnn <= 1.0): continue
            if best is None or aff < best[0]: best = (aff, cnn)
    return best


def parse_vina(out):
    for line in out.splitlines():
        p = line.split()
        if p and p[0].isdigit():
            try: return float(p[1])
            except ValueError: continue
    return None


def worker(gpu, tq, rq):
    env = os.environ.copy(); env['CUDA_VISIBLE_DEVICES'] = str(gpu)
    while True:
        try: lig = tq.get(timeout=5)
        except Empty: break
        name = Path(lig).stem
        try:
            if G['PHASE'] == 'gnina':
                cmd = ['gnina', '-r', G['RECEPTOR'], '-l', lig, *G['BOXARGS'],
                       '--exhaustiveness', str(G['EXH']), '--num_modes', str(G['MODES']),
                       '--seed', str(G['SEED']), '--cpu', str(G['CPU']), '--device', '0',
                       '-o', str(G['POSES'] / f'{name}.sdf')]
            else:  # vina validation
                b = G['BOX']
                rec = G['RECEPTOR_VINA'] if Path(G['RECEPTOR_VINA']).exists() else G['RECEPTOR']
                cmd = ['vina', '--receptor', rec, '--ligand', lig,
                       '--center_x', str(b['center_x']), '--center_y', str(b['center_y']),
                       '--center_z', str(b['center_z']), '--size_x', str(b['size_x']),
                       '--size_y', str(b['size_y']), '--size_z', str(b['size_z']),
                       '--exhaustiveness', str(max(8, G['EXH'])), '--num_modes', str(G['MODES']),
                       '--seed', str(G['SEED']), '--cpu', str(G['CPU']),
                       '--out', str(G['VPOSES'] / f'{name}_vina.pdbqt')]
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=G['TIMEOUT'], env=env)
            if r.returncode != 0:
                raise RuntimeError((r.stderr or r.stdout).strip()[-200:] or f'exit {r.returncode}')
            if G['PHASE'] == 'gnina':
                res = parse_gnina(r.stdout)
                if not res: raise RuntimeError('no valid poses in output')
                rq.put([name, f'{res[0]:.2f}', f'{res[1]:.4f}', 'success'])
            else:
                aff = parse_vina(r.stdout)
                if aff is None: raise RuntimeError('no vina poses')
                rq.put([name, f'{aff:.2f}', '', 'success_vina'])
        except Exception as e:
            rq.put([name, '', '', 'error: ' + str(e).replace('\n', ' ')[:150]])


def run_phase(phase, tasks, gpus):
    G['PHASE'] = phase
    tq, rq = multiprocessing.Queue(), multiprocessing.Queue()
    for t in tasks: tq.put(t)
    procs = [multiprocessing.Process(target=worker, args=(g, tq, rq)) for g in gpus]
    for p in procs: p.start()
    n, t0 = 0, time.time()
    csv_path = G['CSV'] if phase == 'gnina' else G['CSV_VINA']
    sink = open(G['FAIL'], 'a') if phase == 'gnina' else None
    with open(csv_path, 'a', newline='') as f:
        w = csv.writer(f)
        while n < len(tasks):
            try: row = rq.get(timeout=2)
            except Empty:
                if all(not p.is_alive() for p in procs): break
                continue
            n += 1
            if row[3] in ('success', 'success_vina'):
                w.writerow(row)
            elif sink:
                sink.write(f'{row[0]}\t{row[3]}\n')
            f.flush()
            el = time.time() - t0
            rate = n / el * 60 if el else 0
            eta = (el / n) * (len(tasks) - n) if n else 0
            print(f'\r[{phase}] {n}/{len(tasks)} | {rate:.1f} mol/min | left ~{fmt_eta(eta)}   ',
                  end='', flush=True)
    if sink: sink.close()
    for p in procs: p.join()
    print(f'\r[{phase}] {len(tasks)}/{len(tasks)} | done in {fmt_eta(time.time() - t0)}')


def load_state():
    done = set()
    if G['CSV'].exists():
        with open(G['CSV'], newline='') as f:
            for row in csv.reader(f):
                if len(row) == 4 and row[3] == 'success' and is_float(row[1]) \
                   and is_float(row[2]) and 0.0 <= float(row[2]) <= 1.0:
                    done.add(row[0].strip())
    fails = {}
    if G['FAIL'].exists():
        for line in G['FAIL'].read_text().splitlines():
            p = line.split('\t')
            if len(p) >= 2: fails[p[0]] = fails.get(p[0], 0) + 1
    return done, fails


def load_vina_done():
    done = set()
    if G['CSV_VINA'].exists():
        with open(G['CSV_VINA'], newline='') as f:
            for row in csv.reader(f):
                if len(row) == 4 and row[3] == 'success_vina' and is_float(row[1]):
                    done.add(row[0].strip())
    return done


def finish():
    with open(G['CSV'], newline='') as f:
        rows = [r for r in csv.DictReader(f) if is_float(r['affinity'])]
    rows.sort(key=lambda r: float(r['affinity']))
    print('\nTop-15 by affinity:')
    for i, r in enumerate(rows[:15], 1):
        print(f"{i:2d}. {r['chembl_id']:<14} {NAMES.get(r['chembl_id'], '-')[:22]:<22} "
              f"{float(r['affinity']):7.2f}  CNN={float(r['cnn_score']):.4f}")
    rows.sort(key=lambda r: float(r['cnn_score']), reverse=True)
    print('\nTop-10 by CNN score:')
    for i, r in enumerate(rows[:10], 1):
        print(f"{i:2d}. {r['chembl_id']:<14} {NAMES.get(r['chembl_id'], '-')[:22]:<22} "
              f"CNN={float(r['cnn_score']):.4f}  {float(r['affinity']):7.2f}")
    summary = {'receptor': G['RECEPTOR'], 'box': G['BOX'], 'ligand_dir': G['LIGDIR'],
               'results': str(G['RESULTS']), 'exhaustiveness': G['EXH'], 'seed': G['SEED'],
               'num_modes': G['MODES'], 'n_success': len(rows), 'n_permanent_fail': G['N_PERM']}
    (G['RESULTS'] / 'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f'\nsummary.json saved to {G["RESULTS"]}')


def main():
    global NAMES
    ap = argparse.ArgumentParser(description='GNINA virtual screening (2 GPU)')
    base = Path(os.environ.get('SCREENING_ROOT',
                               str(Path(__file__).resolve().parents[1])))
    ap.add_argument('target', nargs='?', default='7gqu')
    ap.add_argument('--ligand-dir', default=str(base / 'pdbqt'))
    ap.add_argument('--results', default=None)
    ap.add_argument('--gpus', default='0,1')
    ap.add_argument('--cpu', type=int, default=4)
    ap.add_argument('--exhaustiveness', type=int, default=4)
    ap.add_argument('--num-modes', type=int, default=9)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--min-box', type=float, default=25.0,
                    help='minimum box size, A (guard against tight box)')
    ap.add_argument('--vina-top', type=int, default=0,
                    help='Vina cross-validation of top-N hits (0 = off)')
    ap.add_argument('--only', default=None,
                    help='file with CHEMBL IDs (one per line)')
    ap.add_argument('--timeout', type=int, default=600)
    ap.add_argument('--limit', type=int, default=None)
    a = ap.parse_args()
    if not shutil.which('gnina'): sys.exit('ERROR: gnina not found in PATH')

    t = a.target.lower()
    G.update(RECEPTOR=str(base / 'data' / t / 'receptor.pdbqt'),
             RECEPTOR_VINA=str(base / 'data' / t / 'receptor_vina.pdbqt'),
             AUTOBOX=str(base / 'data' / t / 'ref_ligand.pdbqt'),
             EXH=a.exhaustiveness, MODES=a.num_modes, SEED=a.seed,
             CPU=a.cpu, TIMEOUT=a.timeout, LIGDIR=a.ligand_dir)

    bj = Path(base / 'data' / t / 'receptor.json')
    if bj.exists():
        b = json.loads(bj.read_text())
        G['BOX'] = {'center_x': float(b['center_x']), 'center_y': float(b['center_y']),
                    'center_z': float(b['center_z']),
                    'size_x': max(a.min_box, float(b['size_x'])),
                    'size_y': max(a.min_box, float(b['size_y'])),
                    'size_z': max(a.min_box, float(b['size_z']))}
    else:
        G['BOX'] = None
    if G['BOX'] is not None:
        bx = G['BOX']
        G['BOXARGS'] = ['--center_x', str(bx['center_x']), '--center_y', str(bx['center_y']),
                        '--center_z', str(bx['center_z']), '--size_x', str(bx['size_x']),
                        '--size_y', str(bx['size_y']), '--size_z', str(bx['size_z'])]
    else:
        G['BOXARGS'] = ['--autobox_ligand', G['AUTOBOX']]

    res = Path(a.results) if a.results else base / 'results' / t
    G['RESULTS'] = res
    G['POSES'] = res / 'poses'
    G['VPOSES'] = res / 'vina_poses'
    G['CSV'] = res / 'all_docking_results.csv'
    G['CSV_VINA'] = res / 'vina_validation.csv'
    G['FAIL'] = res / 'failed.log'
    for d in (res, G['POSES'], G['VPOSES']): d.mkdir(parents=True, exist_ok=True)
    for p in (G['CSV'], G['CSV_VINA']):
        if not p.exists():
            with open(p, 'w', newline='') as f:
                csv.writer(f).writerow(['chembl_id', 'affinity', 'cnn_score', 'status'])

    NAMES = load_names([base / 'data' / 'library' / 'library.tsv'])

    gpus = [int(x) for x in a.gpus.split(',') if x.strip()]
    all_ligs = sorted(p for p in Path(a.ligand_dir).glob('*.pdbqt')
                      if p.stem not in ('receptor', 'ref_ligand'))
    if a.only:
        ids = {x.strip() for x in Path(a.only).read_text().splitlines() if x.strip()}
        all_ligs = [p for p in all_ligs if p.stem in ids]
    if a.limit: all_ligs = all_ligs[:a.limit]

    done, fails = load_state()
    perm = {c for c, n in fails.items() if n >= MAX_FAILS}
    todo = [p for p in all_ligs if p.stem not in done and p.stem not in perm]
    G['N_PERM'] = len(perm)

    print(f'=== Screening target={t.upper()} ===')
    print(f'Receptor : {G["RECEPTOR"]}')
    print(f'Box      : {G["BOX"] if G["BOX"] else "autobox"}')
    print(f'Ligands  : {a.ligand_dir} ({len(all_ligs)})')
    print(f'GPU {gpus} x {a.cpu} CPU | exh={a.exhaustiveness} modes={a.num_modes} seed={a.seed}')
    print(f'Done: {len(done)} | perm.fail: {len(perm)} | left: {len(todo)}')

    if todo:
        run_phase('gnina', [str(p) for p in todo], gpus)

    if a.vina_top > 0 and G['BOX'] is not None:
        with open(G['CSV'], newline='') as f:
            rows = [r for r in csv.DictReader(f) if is_float(r['affinity'])]
        rows.sort(key=lambda r: float(r['affinity']))
        by_id = {p.stem: p for p in all_ligs}
        vdone = load_vina_done()
        vtodo = [str(by_id[i]) for i in [r['chembl_id'] for r in rows[:a.vina_top]]
                 if i in by_id and i not in vdone]
        print(f'Vina cross-validation: top={a.vina_top}, left={len(vtodo)}')
        if vtodo: run_phase('vina', vtodo, gpus)

    finish()


if __name__ == '__main__':
    main()