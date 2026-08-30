#!/usr/bin/env python3
"""Publication figures for ALL final candidates (universal, headless).

Stage 08 of the renumbered pipeline (00 prep, 01 redock gate, 02 library,
03 ligands, 04 screen, 05 consensus, 06 ADMET, 07 candidates, 08 figures,
09 supplementary). Runs BEFORE 09_supplementary so the bundle ZIP already
contains the figures.

For every compound that reached the final list (07_final_candidates.csv)
generates, under supplementary/<target>/figures/ (the same folder where
09_supplementary.py puts journal-ready files):

  overlay/Fig_overlay_CLEAN.png      all CLEAN candidates overlaid in the site
  overlay/Fig_overlay_CLEAN.pse      editable PyMOL session of the overlay
  overlay/Fig_overlay_FLAGGED.png    all FLAGGED candidates overlaid
  overlay/Fig_overlay_FLAGGED.pse    editable PyMOL session of the overlay
  redock/Fig_redock.png              crystal ligand vs top redock pose,
                                     captioned with RMSD from the 01 report
  per_ligand/CLEAN/<name>.png        3D: binding-site sticks + polar contacts
  per_ligand/FLAGGED/<name>.png      3D: same for flagged (liability) compounds
  interactions_2d/CLEAN/<name>.png   2D structure + PLIP interaction summary
  interactions_2d/FLAGGED/<name>.png 2D: same for flagged compounds

Writes FIGURES_DONE.txt at the end - the resume sentinel used by run_all.py
(stages 08 and 09 are considered done only when it exists).

Usage:
  python scripts/08_make_figures.py --target 7gqu
  python scripts/08_make_figures.py --target 9fza            # any target
  python scripts/08_make_figures.py --target 7gqu --skip-2d  # 3D only (no PLIP)
"""
import argparse
import csv
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]


def read_rows(p):
    d = '\t' if str(p).endswith('.tsv') else ','
    with open(p, newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter=d))


def safe(s):
    return ''.join(c if c.isalnum() else '_' for c in s)


def find_redock_poses(t):
    """New gate (01) writes redock_poses.sdf; fall back to legacy names."""
    vdir = BASE / 'validation' / 'redock' / t
    cand = vdir / 'redock_poses.sdf'
    if cand.exists():
        return cand
    for pat in ('redock_exh*.sdf', 'redock*.sdf'):
        hits = sorted(vdir.glob(pat))
        if hits:
            return hits[0]
    return None


def read_redock_report(t):
    """Parse validation/redock/<t>/redock_report.txt into a dict."""
    rpt = BASE / 'validation' / 'redock' / t / 'redock_report.txt'
    info = {}
    if rpt.exists():
        for line in rpt.read_text(errors='ignore').splitlines():
            k, _, v = line.partition('=')
            if v:
                info[k.strip()] = v.strip()
    return info


def add_caption(png, lines, font):
    """Append a white caption bar with text lines below an existing PNG."""
    from PIL import Image, ImageDraw
    img = Image.open(png).convert('RGB')
    w, h = img.size
    bar = 30 * len(lines) + 14
    canvas = Image.new('RGB', (w, h + bar), 'white')
    canvas.paste(img, (0, 0))
    d = ImageDraw.Draw(canvas)
    y = h + 6
    for ln in lines:
        d.text((12, y), ln, fill='black', font=font)
        y += 30
    canvas.save(str(png))


def main():
    ap = argparse.ArgumentParser(description='Figures for all final candidates')
    ap.add_argument('--target', default='7gqu')
    ap.add_argument('--out', default=None,
                    help='figures folder (default: supplementary/<target>/figures)')
    ap.add_argument('--skip-2d', action='store_true', help='skip 2D/PLIP diagrams')
    ap.add_argument('--dpi', type=int, default=300)
    a = ap.parse_args()
    t = a.target.lower()

    fig = Path(a.out) if a.out else BASE / 'supplementary' / t / 'figures'
    for sub in ['overlay', 'redock', 'per_ligand/CLEAN', 'per_ligand/FLAGGED',
                'interactions_2d/CLEAN', 'interactions_2d/FLAGGED']:
        (fig / sub).mkdir(parents=True, exist_ok=True)

    fin = BASE / 'results' / t / 'refine' / 'final_candidates.csv'
    if not fin.exists():
        sys.exit(f'ERROR: {fin} not found - run stages 04/06/07 first')
    rows = read_rows(fin)

    lib = {}
    lp = BASE / 'data' / 'library' / 'library.tsv'
    if lp.exists():
        for r in read_rows(lp):
            lib[r['chembl_id']] = r

    REC = BASE / 'data' / t / 'protein_clean.pdb'
    GDIR = BASE / 'results' / t / 'refine' / 'gnina'
    if not REC.exists():
        sys.exit(f'ERROR: {REC} not found - run stage 00 first')

    # ---------------- fonts / PIL ----------------
    from PIL import Image, ImageDraw, ImageFont
    try:
        import matplotlib
        ttf = str(Path(matplotlib.get_data_path()) / 'fonts' / 'ttf' / 'DejaVuSans.ttf')
        font = ImageFont.truetype(ttf, 20)
        font_small = ImageFont.truetype(ttf, 18)
    except Exception:
        font = ImageFont.load_default()
        font_small = font

    # ---------------- PyMOL (headless) ----------------
    import pymol
    from pymol import cmd, util
    pymol.finish_launching(['pymol', '-cq'])
    cmd.bg_color('white')
    cmd.set('ray_trace_mode', 1)
    cmd.set('cartoon_transparency', 0.75)
    cmd.set('dash_color', 'yellow')

    PAL = ['green', 'cyan', 'yellow', 'magenta', 'orange', 'blue',
           'salmon', 'lime', 'violet', 'red']
    CLS_COL = {'CLEAN': 'green', 'FLAGGED': 'orange'}

    def pose_of(cid):
        """Exact match first, then cid_*, then cid* (avoids prefix collisions
        like CHEMBL53 vs CHEMBL5314423 or CHEMBL1623 vs CHEMBL1623992)."""
        exact = GDIR / f'{cid}.sdf'
        if exact.exists():
            return exact
        und = sorted(GDIR.glob(f'{cid}_*.sdf'))
        if und:
            return und[0]
        pre = sorted(GDIR.glob(f'{cid}*.sdf'))
        return pre[0] if pre else None

    def load_protein():
        cmd.load(str(REC), 'rec')
        cmd.hide('everything', 'rec')
        cmd.show('cartoon', 'rec')
        cmd.color('gray70', 'rec')

    def save(png, w=1200, h=900):
        cmd.png(str(png), width=w, height=h, dpi=a.dpi, ray=1)

    # ---------------- 3D per ligand ----------------
    n3d = 0
    for r in rows:
        cid, name, cls = r['chembl_id'], r['name'], r['cls']
        pose = pose_of(cid)
        if pose is None:
            print(f'[warn] no pose for {cid}; 3D skipped')
            continue
        try:
            cmd.delete('all')
            load_protein()
            cmd.load(str(pose), 'lig')
            cmd.show('sticks', 'lig')
            cmd.color(CLS_COL.get(cls, 'cyan'), 'lig')
            util.cnc('lig')
            cmd.select('site', 'byres (rec) within 5.0 of lig')
            cmd.show('sticks', 'site')
            cmd.color('gray60', 'site')
            util.cnc('site')
            try:
                cmd.distance('polar', 'lig and not elem H', 'site and not elem H',
                             3.5, mode=2)
            except Exception:
                pass  # no polar contacts within cutoff - figure still fine
            cmd.zoom('lig', 6)
            save(fig / 'per_ligand' / cls / f'{safe(name)}.png')
            n3d += 1
            if n3d % 10 == 0:
                print(f'   3D done: {n3d}/{len(rows)}')
        except Exception as e:
            print(f'[warn] 3D failed for {cid}: {type(e).__name__}')

    # ---------------- overlays ----------------
    for tag in ['CLEAN', 'FLAGGED']:
        sub = [r for r in rows if r['cls'] == tag]
        cmd.delete('all')
        load_protein()
        loaded = []
        for i, r in enumerate(sub):
            pose = pose_of(r['chembl_id'])
            if pose is None:
                continue
            try:
                obj = f'lig_{i}'
                cmd.load(str(pose), obj)
                cmd.show('sticks', obj)
                cmd.color(PAL[i % len(PAL)], obj)
                util.cnc(obj)
                loaded.append(obj)
            except Exception as e:
                print(f'[warn] overlay load failed for {r["chembl_id"]}: '
                      f'{type(e).__name__}')
        if loaded:
            cmd.zoom(' or '.join(loaded), 6)
            save(fig / 'overlay' / f'Fig_overlay_{tag}.png', 1600, 1200)
            cmd.save(str(fig / 'overlay' / f'Fig_overlay_{tag}.pse'))
            print(f'[ok] overlay {tag}: {len(loaded)} ligands')

    # ---------------- redock overlay ----------------
    redock_ok = False
    REF = BASE / 'data' / t / 'ref_ligand.sdf'
    if not REF.exists():
        REF = BASE / 'data' / t / 'ref_ligand_prot.sdf'
    RED = find_redock_poses(t)
    rep = read_redock_report(t)
    if REF.exists() and RED is not None:
        try:
            cmd.delete('all')
            load_protein()
            cmd.load(str(REF), 'cryst')
            cmd.show('sticks', 'cryst')
            cmd.color('green', 'cryst')
            util.cnc('cryst')
            cmd.load(str(RED), 'redock')          # state 1 = top pose
            cmd.show('sticks', 'redock')
            cmd.color('magenta', 'redock')
            util.cnc('redock')
            cmd.zoom('cryst', 6)
            png = fig / 'redock' / 'Fig_redock.png'
            save(png, 1600, 1200)
            cap = [f'{t.upper()}: crystal ligand (green) vs top redock pose (magenta)']
            if rep:
                cap.append('RMSD top={} A  best={} A  rank={}  verdict={}'.format(
                    rep.get('rmsd_top_pose', 'n/a'),
                    rep.get('rmsd_best_pose', 'n/a'),
                    rep.get('native_like_rank', 'n/a'),
                    rep.get('verdict', 'n/a')))
            add_caption(png, cap, font_small)
            redock_ok = True
            print('[ok] Fig_redock.png')
        except Exception as e:
            print(f'[warn] redock figure failed: {type(e).__name__}')
    else:
        print('[warn] redock files missing; Fig_redock skipped')

    # ---------------- 2D + PLIP ----------------
    n2d = 0
    tmp = fig.parent / '_tmp_complex.pdb'
    if not a.skip_2d:
        try:
            from plip.structure.preparation import PDBComplex
            PLIP_OK = True
        except Exception as e:
            PLIP_OK = False
            print(f'[warn] PLIP unavailable ({e}); 2D without interactions')

        from rdkit import Chem
        from rdkit.Chem import Draw

        prot_lines = [l for l in REC.read_text().splitlines()
                      if l.startswith(('ATOM', 'TER'))]
        for r in rows:
            cid, name, cls = r['chembl_id'], r['name'], r['cls']
            try:
                pose = pose_of(cid)
                mol = Chem.MolFromMolFile(str(pose), removeHs=True) if pose else None
                mol_h = Chem.MolFromMolFile(str(pose), removeHs=False) if pose else None
                smi = lib.get(cid, {}).get('smiles_clean', '')
                if mol is None and smi:
                    mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    print(f'[warn] no structure for {cid}; 2D skipped')
                    continue

                inter = 'PLIP: n/a'
                if PLIP_OK and mol_h is not None:
                    try:
                        lig_block = '\n'.join(
                            l for l in Chem.MolToPDBBlock(mol_h).splitlines()
                            if not l.startswith(('END', 'TER', 'CONECT')))
                        tmp.write_text('\n'.join(prot_lines) + '\n'
                                       + lig_block + '\nEND\n')
                        p = PDBComplex()
                        p.load_pdb(str(tmp))
                        p.analyze()
                        hb = hy = pi = sa = ha = wb = 0
                        for site in p.interaction_sets.values():
                            cnt = lambda at: len(getattr(site, at, []) or [])
                            hb += cnt('hb_pdon') + cnt('hb_ldon')
                            hy += cnt('hydrophobic')
                            pi += cnt('pistacking')
                            ha += cnt('halogen')
                            wb += cnt('water_bridges')
                            sa += cnt('saltbridge_lpos') + cnt('saltbridge_pneg')
                        inter = (f'PLIP: H-bonds {hb} | hydrophobic {hy} | '
                                 f'pi-stack {pi} | salt-bridges {sa} | '
                                 f'halogen {ha} | water-bridges {wb}')
                    except Exception as e:
                        inter = f'PLIP: error ({type(e).__name__})'

                img = Draw.MolToImage(mol, size=(900, 560))
                canvas = Image.new('RGB', (900, 810), 'white')
                canvas.paste(img, (0, 0))
                d = ImageDraw.Draw(canvas)
                d.text((12, 570), f'{name}  ({cid})  [{cls}]', fill='black', font=font)
                line2 = (f"GNINA {r['gnina32']} | Vina {r['vina']} | "
                         f"LeDock {r['ledock']} kcal/mol; QED {r['qed']}")
                if r.get('dLedock'):
                    line2 += f"; dLeDock {r['dLedock']}"
                d.text((12, 600), line2, fill='black', font=font)
                d.text((12, 630), inter, fill='black', font=font)
                if r.get('flags'):
                    d.text((12, 660), f'FLAGS: {r["flags"]}', fill='red', font=font)
                canvas.save(str(fig / 'interactions_2d' / cls / f'{safe(name)}.png'))
                n2d += 1
                if n2d % 10 == 0:
                    print(f'   2D done: {n2d}/{len(rows)}')
            except Exception as e:
                print(f'[warn] 2D failed for {cid}: {type(e).__name__}')
        tmp.unlink(missing_ok=True)

    # ---------------- resume sentinel for run_all ----------------
    (fig / 'FIGURES_DONE.txt').write_text(
        f'3D={n3d}; 2D={n2d}; overlays=CLEAN,FLAGGED; '
        f'redock={"ok" if redock_ok else "skipped"}\n')
    print(f'\n✅ Figures: {fig}')
    print(f'   3D per-ligand: {n3d}; overlays: CLEAN+FLAGGED; '
          f'redock: {"ok" if redock_ok else "skipped"}; 2D: {n2d}')


if __name__ == '__main__':
    main()