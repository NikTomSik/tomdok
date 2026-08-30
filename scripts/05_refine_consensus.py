"""Two-stage top-hit re-evaluation + 3-engine consensus (stage 05 of 00-09).

Stage 0: primary-screen filter (aff <= aff_cut AND cnn >= cnn_cut).
         IDs without PDBQT (salts from older libraries) are silently skipped.
Stage 1: GNINA exh=32 on filtered (resume: gnina_refined_top.csv);
         re-filter on refined values.
Stage 2: Vina consensus (exh=32) + LeDock,
         resume per engine CSV; Vina cache is pulled from consensus.csv.

Engines and their receptors:
  GNINA  -> receptor.pdbqt (Gasteiger charges)
  Vina   -> receptor_vina.pdbqt (without ROOT/BRANCH/TORSDOF; auto-generated)
  LeDock -> pro.pdb (lepro output); config: Receptor / Binding pocket /
            Number of binding poses / Ligands list; scores from .dok

New (v2): if 01_redock_validation.py has written a per-target LeDock baseline
(line 'ledock_native=<x>' in validation/redock/<t>/redock_report.txt), the
consensus table gains a 'dledock' column = ledock_score - ledock_native, which
makes the systematically-soft LeDock scale interpretable per target
(negative dledock = candidate scores better than the native ligand).

Usage:
  python scripts/05_refine_consensus.py 7gqu
  python scripts/05_refine_consensus.py 7gqu --primary /path/all_docking_results.csv
"""
import argparse, csv, json, os, re, shutil, subprocess, time
import multiprocessing as mp
from pathlib import Path
from queue import Empty

BASE = Path(os.environ.get('SCREENING_ROOT', str(Path(__file__).resolve().parents[1])))
FLEX_TAGS = ('ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH', 'TORSDOF')
G, NAMES = {}, {}


def fmt_eta(sec):
    sec = max(0, int(sec))
    if sec < 60: return f'{sec} s'
    if sec < 3600: return f'{sec // 60} min {sec % 60} s'
    return f'{sec // 3600} h {(sec % 3600) // 60} min'


def box_cli():
    b = G['BOX']
    return ['--center_x', str(b['center_x']), '--center_y', str(b['center_y']),
            '--center_z', str(b['center_z']),
            '--size_x', str(max(25.0, b['size_x'])), '--size_y', str(max(25.0, b['size_y'])),
            '--size_z', str(max(25.0, b['size_z']))]


def parse_pose_table(out):
    best = None
    for line in out.splitlines():
        p = line.split()
        if len(p) >= 2 and p[0].isdigit():
            try: aff = float(p[1])
            except ValueError: continue
            cnn = None
            if len(p) >= 4:
                try: cnn = float(p[3])
                except ValueError: cnn = None
            if best is None or aff < best[0]: best = (aff, cnn)
    return best


def load_names(paths):
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


def read_ledock_native():
    """Per-target LeDock baseline written by 01_redock_validation.py
    ('ledock_native=<value>' in redock_report.txt, or ledock_baseline.json)."""
    vdir = BASE / 'validation' / 'redock' / G['TARGET']
    rpt = vdir / 'redock_report.txt'
    if rpt.exists():
        for line in rpt.read_text(errors='ignore').splitlines():
            if line.startswith('ledock_native'):
                try: return float(line.split('=')[-1].strip())
                except ValueError: return None
    js = vdir / 'ledock_baseline.json'
    if js.exists():
        try: return float(json.loads(js.read_text()).get('ledock_native'))
        except Exception: return None
    return None


# ---------------- Stage 1: GNINA exh=32 (GPU) ----------------
def gnina_worker(gpu, tq, rq):
    env = os.environ.copy(); env['CUDA_VISIBLE_DEVICES'] = str(gpu)
    gdir = G['OUT'] / 'gnina'
    while True:
        try: lig = tq.get(timeout=5)
        except Empty: break
        cid = Path(lig).stem
        out = gdir / f'{cid}_ref.sdf'
        cmd = ['gnina', '-r', str(G['RECEPTOR_PDBQT']), '-l', str(lig), *box_cli(),
               '--exhaustiveness', str(G['EXH']), '--num_modes', '9', '--seed', '42',
               '--cpu', '4', '--device', '0', '-o', str(out)]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800, env=env)
            res = parse_pose_table(r.stdout) if r.returncode == 0 else None
            rq.put((cid, res[0], res[1]) if res else (cid, None, None))
        except Exception:
            rq.put((cid, None, None))


def stage1(ids):
    (G['OUT'] / 'gnina').mkdir(parents=True, exist_ok=True)
    csv_p = G['OUT'] / 'gnina_refined_top.csv'
    done = {}
    if csv_p.exists() and csv_p.stat().st_size > 0:
        with open(csv_p, newline='') as f:
            for r in csv.DictReader(f):
                if r.get('affinity'):
                    done[r['chembl_id']] = (float(r['affinity']),
                                            float(r['cnn_score']) if r.get('cnn_score') else None)
    todo = [G['LIG_DIR'] / f'{c}.pdbqt' for c in ids
            if c not in done and (G['LIG_DIR'] / f'{c}.pdbqt').exists()]
    print(f'[Stage 1] GNINA exh={G["EXH"]}: done={len(done)} left={len(todo)}')
    hdr = (not csv_p.exists()) or csv_p.stat().st_size == 0
    if todo:
        tq, rq = mp.Queue(), mp.Queue()
        for t in todo: tq.put(t)
        procs = [mp.Process(target=gnina_worker, args=(g, tq, rq)) for g in G['GPUS']]
        for p in procs: p.start()
        n, t0 = 0, time.time()
        with open(csv_p, 'a', newline='') as f:
            w = csv.writer(f)
            if hdr: w.writerow(['chembl_id', 'affinity', 'cnn_score'])
            while n < len(todo):
                try: cid, aff, cnn = rq.get(timeout=2)
                except Empty:
                    if all(not p.is_alive() for p in procs): break
                    continue
                n += 1
                if aff is not None:
                    w.writerow([cid, f'{aff:.2f}', f'{cnn:.4f}' if cnn is not None else ''])
                    f.flush(); done[cid] = (aff, cnn)
                el = time.time() - t0
                eta = (el / n) * (len(todo) - n) if n else 0
                print(f'\r[Stage 1] {n}/{len(todo)} | left ~{fmt_eta(eta)}   ',
                      end='', flush=True)
        print(f'\r[Stage 1] {len(todo)}/{len(todo)} | done in {fmt_eta(time.time() - t0)}')
        for p in procs: p.join()
    return done


# ---------------- Stage 2: engines (CPU) ----------------
def run_vina(cid):
    rec = G['RECEPTOR_VINA'] if G['RECEPTOR_VINA'].exists() else G['RECEPTOR_PDBQT']
    out = G['OUT'] / 'vina' / f'{cid}_vina.pdbqt'
    cmd = ['vina', '--receptor', str(rec), '--ligand', str(G['LIG_DIR'] / f'{cid}.pdbqt'),
           *box_cli(), '--exhaustiveness', str(G['EXH']), '--num_modes', '9', '--seed', '42',
           '--cpu', '1', '--out', str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    res = parse_pose_table(r.stdout)
    if res is None and not G.get('VINA_DBG'):
        G['VINA_DBG'] = True
        print(f'\n[Vina DEBUG] {(r.stderr or r.stdout)[-300:]}')
    return (cid, res[0]) if res else (cid, None)


def ledock_prep():
    """Receptor for LeDock: lepro output = pro.pdb."""
    wd = G['OUT'] / 'ledock'; wd.mkdir(parents=True, exist_ok=True)
    for cand in (wd / 'receptor_lepro.pdb', wd / 'pro.pdb'):
        if cand.exists(): return cand
    if not (shutil.which('lepro') and G['RECEPTOR_PDB'].exists()): return None
    subprocess.run(['lepro', str(G['RECEPTOR_PDB'])], cwd=wd, capture_output=True, text=True)
    rec = wd / 'pro.pdb'
    return rec if rec.exists() else None


def parse_dok(path):
    """Extracts pose scores from LeDock output .dok."""
    txt = path.read_text(errors='ignore')
    for pat in (r'score\s*[:=]?\s*(-\d+\.\d+)',
                r'remark[^\n-]*(-\d+\.\d+)',
                r'energy\s*[:=]?\s*(-\d+\.\d+)',
                r'pose\s*\d+\s*(-\d+\.\d+)'):
        vals = [float(x) for x in re.findall(pat, txt, re.I)]
        if vals: return vals
    vals = []
    for line in txt.splitlines():
        if line[:1] in ('#', '/') or line.upper().startswith('REMARK'):
            vals += [float(x) for x in re.findall(r'-\d+\.\d+', line)]
    return vals


def run_ledock(cid):
    rec = G.get('LED_REC')
    if not rec: return (cid, None)
    wd = G['OUT'] / 'ledock'
    m = wd / f'{cid}.mol2'
    if not m.exists():
        subprocess.run(['obabel', str(G['LIG_DIR'] / f'{cid}.pdbqt'), '-O', str(m),
                        '--partialcharge', 'gasteiger'], capture_output=True)
    if not m.exists(): return (cid, None)
    lst = wd / f'{cid}.list'
    lst.write_text(str(m) + '\n')
    dok = wd / f'{cid}.dok'
    if dok.exists(): dok.unlink()
    b = G['BOX']
    hx, hy, hz = (max(25.0, b[f'size_{k}']) / 2 for k in 'xyz')
    cfg = wd / f'{cid}.in'
    cfg.write_text(
        f"Receptor\n{rec}\n\n"
        f"Binding pocket\n"
        f"{b['center_x'] - hx:.3f} {b['center_x'] + hx:.3f}\n"
        f"{b['center_y'] - hy:.3f} {b['center_y'] + hy:.3f}\n"
        f"{b['center_z'] - hz:.3f} {b['center_z'] + hz:.3f}\n\n"
        f"Other\nNumber of binding poses\n10\n\n"
        f"Ligands list\n{lst}\n")
    subprocess.run(['ledock', str(cfg)], capture_output=True, text=True, timeout=1800)
    if not dok.exists(): return (cid, None)
    vals = parse_dok(dok)
    return (cid, min(vals)) if vals else (cid, None)


def run_pool(label, func, ids, workers):
    res, n, t0 = {}, 0, time.time()
    with mp.Pool(workers) as pool:
        for cid, val in pool.imap_unordered(func, ids):
            res[cid] = val; n += 1
            el = time.time() - t0
            eta = (el / n) * (len(ids) - n)
            print(f'\r[{label}] {n}/{len(ids)} | left ~{fmt_eta(eta)}   ',
                  end='', flush=True)
    print(f'\r[{label}] {len(ids)}/{len(ids)} | done in {fmt_eta(time.time() - t0)}')
    return res


def load_csv(p):
    d = {}
    if p.exists():
        with open(p, newline='') as f:
            for r in csv.DictReader(f):
                if r.get('score'):
                    try: d[r['chembl_id']] = float(r['score'])
                    except ValueError: pass
    return d


def save_csv(p, d):
    with open(p, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['chembl_id', 'score'])
        for k, v in d.items(): w.writerow([k, f'{v:.2f}'])


def main():
    global NAMES
    ap = argparse.ArgumentParser()
    ap.add_argument('target', nargs='?', default='7gqu')
    ap.add_argument('--primary', default=None)
    ap.add_argument('--aff-cut', type=float, default=-9.0)
    ap.add_argument('--cnn-cut', type=float, default=0.7)
    ap.add_argument('--exh', type=int, default=32)
    ap.add_argument('--gpus', default='0,1')
    ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()
    t = a.target.lower()

    G['TARGET'] = t
    G['EXH'] = a.exh
    G['GPUS'] = [int(x) for x in a.gpus.split(',') if x.strip()]
    G['LIG_DIR'] = BASE / 'pdbqt'
    G['OUT'] = BASE / 'results' / t / 'refine'
    G['OUT'].mkdir(parents=True, exist_ok=True)
    data = BASE / 'data' / t
    G['RECEPTOR_PDBQT'] = data / 'receptor.pdbqt'
    G['RECEPTOR_VINA'] = data / 'receptor_vina.pdbqt'
    G['RECEPTOR_PDB'] = data / 'protein_clean.pdb'
    G['BOX'] = json.loads((data / 'receptor.json').read_text())
    if a.primary is None:
        a.primary = str(BASE / 'results' / t / 'all_docking_results.csv')

    # Vina receptor without ROOT/BRANCH: auto-generated if missing
    if not G['RECEPTOR_VINA'].exists() and G['RECEPTOR_PDBQT'].exists():
        lines = [l for l in G['RECEPTOR_PDBQT'].read_text().splitlines()
                 if not l.startswith(FLEX_TAGS)]
        G['RECEPTOR_VINA'].write_text('\n'.join(lines) + '\n')
        print('[prep] receptor_vina.pdbqt built from receptor.pdbqt (without ROOT/BRANCH)')

    NAMES = load_names([BASE / 'data' / 'library' / 'library.tsv',
                        BASE / 'lists' / 'approved_clean.tsv',
                        BASE / 'prep_v2' / 'data' / 'approved_phase3.tsv'])

    # ---- Stage 0: primary-screen filter ----
    with open(a.primary, newline='') as f:
        rows = [r for r in csv.DictReader(f) if r.get('affinity') and r.get('cnn_score')]
    cand = [r for r in rows if float(r['affinity']) <= a.aff_cut and float(r['cnn_score']) >= a.cnn_cut]
    cand.sort(key=lambda r: float(r['affinity']))
    top = [r['chembl_id'] for r in cand]
    missing = [c for c in top if not (G['LIG_DIR'] / f'{c}.pdbqt').exists()]
    if missing:
        print(f'[Warning] no PDBQT (salts/legacy IDs), skipped: {len(missing)}')
    top = [c for c in top if c not in set(missing)]
    print(f'[Input] passed primary-screen filter (aff<={a.aff_cut}, cnn>={a.cnn_cut}): {len(top)}')

    # ---- Stage 1: GNINA exh=32 ----
    t_start = time.time()
    refined = stage1(top)
    sel = sorted(((c, refined[c]) for c in top if c in refined), key=lambda x: x[1][0])
    with open(G['OUT'] / 'refined_filtered.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['rank', 'chembl_id', 'affinity', 'cnn_score'])
        for i, (c, v) in enumerate(sel, 1):
            w.writerow([i, c, f'{v[0]:.2f}', f'{v[1]:.4f}' if v[1] is not None else ''])
    final = [(c, v) for c, v in sel if v[0] <= a.aff_cut and v[1] is not None and v[1] >= a.cnn_cut]
    print(f'[Post-refinement filter] aff <= {a.aff_cut} and CNN >= {a.cnn_cut}: '
          f'{len(final)} of {len(sel)}')
    if not final:
        print('Nothing to validate.'); return
    ids = [c for c, _ in final]
    results = {c: {'gnina32': refined[c][0], 'cnn': refined[c][1]} for c in ids}

    # ---- Stage 2: consensus with per-engine resume ----
    vina_csv = G['OUT'] / 'vina_scores.csv'
    ledo_csv = G['OUT'] / 'ledock_scores.csv'
    # Vina cache from old consensus (avoid recomputing)
    if not vina_csv.exists() and (G['OUT'] / 'consensus.csv').exists():
        with open(G['OUT'] / 'consensus.csv', newline='') as f:
            vd = {r['chembl_id']: float(r['vina_aff']) for r in csv.DictReader(f) if r.get('vina_aff')}
        if vd:
            save_csv(vina_csv, vd); print(f'[resume] Vina cache from consensus.csv: {len(vd)}')
    vr = load_csv(vina_csv)
    if shutil.which('vina'):
        todo = [c for c in ids if c not in vr and (G['LIG_DIR'] / f'{c}.pdbqt').exists()]
        if todo:
            (G['OUT'] / 'vina').mkdir(parents=True, exist_ok=True)
            print(f'[Stage 2] Vina exh={a.exh}: {len(todo)} molecules')
            vr.update({k: v for k, v in run_pool('Vina', run_vina, todo, a.workers).items() if v is not None})
            save_csv(vina_csv, vr)
    else:
        print('SKIP Vina: conda install -c bioconda autodock-vina')
    for c in ids: results[c]['vina'] = vr.get(c)

    lr = load_csv(ledo_csv)
    if shutil.which('ledock') and shutil.which('lepro'):
        G['LED_REC'] = ledock_prep()
        todo = [c for c in ids if c not in lr and (G['LIG_DIR'] / f'{c}.pdbqt').exists()]
        if G['LED_REC'] and todo:
            print(f'[Stage 2] LeDock: {len(todo)} molecules')
            new = run_pool('LeDock', run_ledock, todo, a.workers)
            lr.update({k: v for k, v in new.items() if v is not None})
            save_csv(ledo_csv, lr)
            if not any(v is not None for v in new.values()):
                sample = G['OUT'] / 'ledock' / f'{todo[0]}.dok'
                if sample.exists():
                    print(f'[LeDock DEBUG] scores not parsed; first lines of {sample.name}:')
                    print('\n'.join(sample.read_text(errors='ignore').splitlines()[:12]))
        elif not G['LED_REC']:
            print('SKIP LeDock: lepro failed to prepare receptor')
    else:
        print('SKIP LeDock: not found (ledock/lepro)')
    for c in ids: results[c]['ledock'] = lr.get(c)

    # ---- Consensus (+ per-target LeDock baseline from stage 01) ----
    native_ld = read_ledock_native()
    if native_ld is None:
        print('[note] ledock_native not found - run 01_redock_validation.py first; '
              'dledock column left empty')
    out = G['OUT'] / 'consensus.csv'
    with open(out, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['chembl_id', 'name', 'gnina32_aff', 'cnn', 'vina_aff',
                    'ledock_score', 'n_engines_le9', 'dledock'])
        for c in ids:
            r = results[c]
            vals = [r['gnina32'], r.get('vina'), r.get('ledock')]
            n = sum(1 for v in vals if v is not None and v <= -9.0)
            dld = ''
            if r.get('ledock') is not None and native_ld is not None:
                dld = f"{r['ledock'] - native_ld:.2f}"
            w.writerow([c, NAMES.get(c, ''), f"{r['gnina32']:.2f}", f"{r['cnn']:.4f}",
                        f"{r['vina']:.2f}" if r.get('vina') is not None else '',
                        f"{r['ledock']:.2f}" if r.get('ledock') is not None else '', n, dld])
    print(f'\n[Done] Consensus: {out} | total time: {fmt_eta(time.time() - t_start)}')
    hdr = (f'{"CHEMBL":<14}{"name":<20}{"GNINA32":>8}{"CNN":>7}{"VINA":>7}'
           f'{"LeDock":>8}{"N<=-9":>6}')
    if native_ld is not None:
        hdr += f'{"dLeDock":>8}'
    print(hdr)
    for c in ids:
        r = results[c]
        vals = [r['gnina32'], r.get('vina'), r.get('ledock')]
        n = sum(1 for v in vals if v is not None and v <= -9.0)
        line = (f"{c:<14}{NAMES.get(c, '-')[:19]:<20}{r['gnina32']:8.2f}{r['cnn']:7.3f}"
                f"{r.get('vina') or 0:7.2f}{r.get('ledock') or 0:8.2f}{n:6d}")
        if native_ld is not None and r.get('ledock') is not None:
            line += f"{r['ledock'] - native_ld:8.2f}"
        print(line)


if __name__ == '__main__':
    main()