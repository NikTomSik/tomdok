"""Build a journal/Zenodo supplementary bundle (universal, target-agnostic).

Outputs (supplementary/<target>/ by default):
  S0_manifest.txt            manifest: target, funnel, parameters, validation, figures
  S1_library.tsv             full repurposing library
  S2_primary_screen_full.csv primary GNINA screen (all docked compounds)
  S3_gnina_refined.csv       refined GNINA (exh=32) top hits
  S4_consensus.csv           3-engine consensus + phase/approval + safety class
  S5_admet.csv               full ADMET profiles of consensus hits
  S6_final_candidates.csv    final CLEAN/FLAGGED prioritisation
  S7_redock_validation.txt   redocking RMSD report
  S8_parameters.json         machine-readable parameters
  S9_run_all_log.txt         run_all execution log (most complete run)
  poses/                     receptor.pdb + per-ligand pose SDFs (any viewer)
  figures/                   3D/2D figures from stage 09 (included in ZIP)

Usage:
 python scripts/08_supplementary.py --target 7gqu --zip
 python scripts/08_supplementary.py --target 4cgz --target-name "..." --zip
"""
import argparse, csv, json, re, shutil
from datetime import date
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
REPO = 'https://github.com/NikTomSik/tomdok'

def n_lines(p):
    with open(p, encoding='utf-8') as f:
        return sum(1 for _ in f) - 1

def read_rows(p):
    d = '\t' if str(p).endswith('.tsv') else ','
    with open(p, newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter=d))

def write_rows(p, rows, cols):
    with open(p, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow({c: ('' if r.get(c) is None else r.get(c)) for c in cols})

def first_existing(*paths):
    for p in paths:
        if p and Path(p).exists():
            return Path(p)
    return None

def parse_rmsd(txt):
    top = re.search(r'rmsd_top_pose=([0-9.]+)', txt) or re.search(r'Top-pose RMSD\s*:\s*([0-9.]+)', txt)
    best = re.search(r'rmsd_best_pose=([0-9.]+)', txt) or re.search(r'Best-pose RMSD\s*:\s*([0-9.]+)', txt)
    return (top.group(1) if top else None, best.group(1) if best else None)

def main():
    ap = argparse.ArgumentParser(description='Build supplementary bundle')
    ap.add_argument('--target', default='7gqu')
    ap.add_argument('--out', default=None)
    ap.add_argument('--zip', action='store_true')
    ap.add_argument('--target-name', default='',
                    help='free-text target description for the manifest, e.g. '
                         '"WRN helicase (Werner syndrome ATP-dependent helicase), Homo sapiens"')
    a = ap.parse_args()
    t = a.target.lower()
    out = Path(a.out) if a.out else BASE / 'supplementary' / t
    out.mkdir(parents=True, exist_ok=True)

    res = BASE / 'results' / t
    refine = res / 'refine'
    gdir = refine / 'gnina'
    notes = []

    def pose_of(cid):
        """Exact match first, then cid_*, then cid* (avoids prefix collisions)."""
        exact = gdir / f'{cid}.sdf'
        if exact.exists():
            return exact
        und = sorted(gdir.glob(f'{cid}_*.sdf'))
        if und:
            return und[0]
        pre = sorted(gdir.glob(f'{cid}*.sdf'))
        return pre[0] if pre else None

    # ---- S1: library ----
    lib = BASE / 'data' / 'library' / 'library.tsv'
    lib_n = None
    if lib.exists():
        shutil.copy(lib, out / 'S1_library.tsv'); lib_n = n_lines(lib)
        notes.append(('S1_library.tsv', True, f'{lib_n} compounds'))
    else:
        notes.append(('S1_library.tsv', False, 'library.tsv not found'))

    # ---- S2: primary screen ----
    s2 = res / 'all_docking_results.csv'
    dock_n = None
    if s2.exists():
        shutil.copy(s2, out / 'S2_primary_screen_full.csv'); dock_n = n_lines(s2)
        notes.append(('S2_primary_screen_full.csv', True, f'{dock_n} docked'))
    else:
        notes.append(('S2_primary_screen_full.csv', False, 'all_docking_results.csv not found'))

    # ---- S3: gnina refined ----
    s3 = refine / 'gnina_refined_top.csv'
    prim_n = None
    if s3.exists():
        shutil.copy(s3, out / 'S3_gnina_refined.csv'); prim_n = n_lines(s3)
        notes.append(('S3_gnina_refined.csv', True, f'{prim_n} primary hits'))
    else:
        notes.append(('S3_gnina_refined.csv', False, 'gnina_refined_top.csv not found'))

    # ---- S4: consensus joined with phase/approval/safety ----
    cons_p = refine / 'consensus.csv'
    fin_p = refine / 'final_candidates.csv'
    adm_p = refine / 'admet_all.csv'
    cons_n = clean_n = flag_n = None
    if cons_p.exists():
        cons = {r['chembl_id']: r for r in read_rows(cons_p)}
        fin = {r['chembl_id']: r for r in read_rows(fin_p)} if fin_p.exists() else {}
        adm = {r['chembl_id']: r for r in read_rows(adm_p)} if adm_p.exists() else {}
        joined = []
        for cid, c in cons.items():
            f, ad = fin.get(cid, {}), adm.get(cid, {})
            joined.append(dict(
                chembl_id=cid, name=c.get('name', ''),
                phase=f.get('phase', ''), approved=f.get('approved', ''),
                gnina32_aff=c.get('gnina32_aff', ''), cnn=c.get('cnn', ''),
                vina_aff=c.get('vina_aff', ''), ledock_score=c.get('ledock_score', ''),
                n_engines_le9=c.get('n_engines_le9', ''),
                qed=f.get('qed', '') or ad.get('qed', ''),
                bbb=f.get('bbb', '') or ad.get('bbb', ''),
                cls=f.get('cls', ''), flags=f.get('flags', '')))
        joined.sort(key=lambda r: float(r['gnina32_aff'] or 0))
        write_rows(out / 'S4_consensus.csv', joined,
                   ['chembl_id', 'name', 'phase', 'approved', 'gnina32_aff', 'cnn',
                    'vina_aff', 'ledock_score', 'n_engines_le9', 'qed', 'bbb', 'cls', 'flags'])
        cons_n = len(joined)
        clean_n = sum(1 for r in joined if r['cls'] == 'CLEAN')
        flag_n = cons_n - clean_n
        notes.append(('S4_consensus.csv', True,
                      f'{cons_n} hits ({clean_n} CLEAN / {flag_n} FLAGGED)'))
    else:
        notes.append(('S4_consensus.csv', False, 'consensus.csv not found'))

    # ---- S5: ADMET ----
    if adm_p.exists():
        shutil.copy(adm_p, out / 'S5_admet.csv')
        notes.append(('S5_admet.csv', True, f'{n_lines(adm_p)} profiles'))
    else:
        notes.append(('S5_admet.csv', False, 'admet_all.csv not found'))

    # ---- S6: final candidates ----
    if fin_p.exists():
        shutil.copy(fin_p, out / 'S6_final_candidates.csv')
        notes.append(('S6_final_candidates.csv', True, f'{n_lines(fin_p)} candidates'))
    else:
        notes.append(('S6_final_candidates.csv', False, 'final_candidates.csv not found'))

    # ---- S7: redock report ----
    rep = first_existing(BASE / 'validation' / 'redock' / t / 'redock_report.txt',
                         BASE / 'validation' / 'redock' / 'redock_report.txt')
    rmsd_top = rmsd_best = None
    if rep:
        shutil.copy(rep, out / 'S7_redock_validation.txt')
        rmsd_top, rmsd_best = parse_rmsd(rep.read_text(encoding='utf-8'))
        notes.append(('S7_redock_validation.txt', True,
                      f'RMSD top={rmsd_top} A, best={rmsd_best} A'))
    else:
        notes.append(('S7_redock_validation.txt', False, 'redock_report.txt not found'))

    # ---- S8: parameters (machine-readable) ----
    box_p = BASE / 'data' / t / 'receptor.json'
    box = json.loads(box_p.read_text()) if box_p.exists() else None
    params = {
        'target_pdb_id': t.upper(),
        'generated': str(date.today()),
        'pipeline': REPO,
        'box': box,
        'seed': 42,
        'primary_screen': {'engine': 'GNINA', 'exhaustiveness': 4, 'modes': 9,
                           'thresholds': {'affinity_kcal_mol_le': -9.0, 'cnn_ge': 0.7}},
        'refinement': {'gnina_exhaustiveness': 32, 'vina_exhaustiveness': 32,
                       'ledock': 'default'},
        'consensus': 'n_engines_le9 = number of engines with score <= -9.0 kcal/mol',
        'safety': {'herg': 'ChEMBL KCNH2/CHEMBL240 IC50/Ki < 10 uM',
                   'withdrawn': 'ChEMBL withdrawn_flag',
                   'pains': 'RDKit FilterCatalog',
                   'lipinski': 'RDKit, flag if violations > 1'},
    }
    (out / 'S8_parameters.json').write_text(json.dumps(params, indent=2), encoding='utf-8')
    notes.append(('S8_parameters.json', True, 'machine-readable parameters'))

    # ---- S9: run_all execution log (largest = most complete run) ----
    logs = sorted((BASE / 'logs').glob(f'run_all_{t}*.log'),
                  key=lambda q: q.stat().st_size, reverse=True)
    log_p = logs[0] if logs else None
    if log_p:
        shutil.copy(log_p, out / 'S9_run_all_log.txt')
        notes.append(('S9_run_all_log.txt', True, str(log_p.relative_to(BASE))))
    else:
        notes.append(('S9_run_all_log.txt', False, 'run run_all.py to generate the log'))

    # ---- figures (stage 09) ----
    fig_done = out / 'figures' / 'FIGURES_DONE.txt'
    if fig_done.exists():
        notes.append(('figures/ (3D+2D)', True, fig_done.read_text(encoding='utf-8').strip()))
    else:
        notes.append(('figures/ (3D+2D)', False, 'run 09_make_figures.py to generate'))

    # ---- poses: universal 3D data (open in ANY viewer) ----
    poses_out = out / 'poses'
    poses_out.mkdir(parents=True, exist_ok=True)
    rec_pdb = BASE / 'data' / t / 'protein_clean.pdb'
    if rec_pdb.exists():
        shutil.copy(rec_pdb, poses_out / 'receptor.pdb')
    n_poses = 0
    if fin_p.exists():
        for r in read_rows(fin_p):
            src = pose_of(r['chembl_id'])
            if src:
                shutil.copy(src, poses_out / f"{r['chembl_id']}.sdf"); n_poses += 1
    notes.append(('poses/ (PDB+SDF)', n_poses > 0, f'{n_poses} poses + receptor.pdb'))

    # ---- consistency warnings ----
    pdbqt_n = len(list((BASE / 'pdbqt').glob('*.pdbqt')))
    warns = []
    if lib_n and pdbqt_n and pdbqt_n < lib_n:
        warns.append(f'library ({lib_n}) -> prepared PDBQT ({pdbqt_n}): '
                     f'-{lib_n - pdbqt_n} due to stage-02 standardize/dedup and failed 3D prep '
                     f'(expected attrition, see stage-02 log)')
    if pdbqt_n and dock_n is not None and dock_n < pdbqt_n:
        warns.append(f'docked ({dock_n}) < prepared ({pdbqt_n}): '
                     f'see results/{t}/failed.log or re-run 03')

    # ---- S0: manifest ----
    verdict = ''
    if rmsd_best:
        verdict = ' -> PASS' if float(rmsd_best) < 2.0 else ' -> CHECK'
    m = ['=' * 64,
         'tomdok supplementary manifest',
         '=' * 64,
         f'Target        : {a.target_name or ("PDB " + t.upper())}',
         f'PDB ID        : {t.upper()}',
         f'Generated     : {date.today().isoformat()}',
         f'Pipeline      : tomdok ({REPO})',
         'Code license  : MIT',
         'Data license  : CC BY 4.0',
         '',
         '--- Screening funnel ---',
         f'Library (ChEMBL 37, phase >= 3, deduplicated) : {lib_n}',
         f'Prepared PDBQT (stage 02)                     : {pdbqt_n}',
         f'Docked successfully (stage 03)                : {dock_n}',
         f'Primary hits (aff <= -9.0, CNN >= 0.7)        : {prim_n}',
         f'Consensus hits (GNINA32 + Vina + LeDock)      : {cons_n}',
         f'CLEAN / FLAGGED                               : {clean_n} / {flag_n}',
         '',
         'Files:']
    for name, ok, extra in notes:
        m.append(f'  {name:<28} {"OK " if ok else "MISSING"}  {extra}')
    m += ['',
          '--- Parameters ---',
          f'Box           : {box}',
          'Primary screen: GNINA exh=4, modes=9, seed=42; thresholds aff <= -9.0 kcal/mol, CNN >= 0.7',
          'Refinement    : GNINA exh=32 + AutoDock Vina exh=32 + LeDock',
          'Consensus     : n_engines_le9 = #engines with score <= -9.0 kcal/mol',
          'Safety        : hERG/withdrawn from ChEMBL bioactivity/flags; PAINS (RDKit)',
          f'Validation    : self-docking RMSD top={rmsd_top} A, best={rmsd_best} A '
          f'(threshold 2.0 A){verdict}',
          '',
          'Column glossary: see README.md, section "Output columns".']
    if warns:
        m += ['', '--- Consistency notes ---']
        m += [f'  [!] {w}' for w in warns]
    (out / 'S0_manifest.txt').write_text('\n'.join(m) + '\n', encoding='utf-8')

    print(f'\n✅ Supplementary bundle: {out}')
    for name, ok, extra in notes:
        print(f'   {"OK " if ok else "-- "} {name:<28} {extra}')
    for w in warns:
        print(f'   [!] {w}')

    if a.zip:
        z = shutil.make_archive(str(out.parent / f'supplementary_{t}'), 'zip', out)
        print(f'   ZIP: {z} ({Path(z).stat().st_size / 1e6:.1f} MB)')

if __name__ == '__main__':
    main()