#!/usr/bin/env python3
"""Publication figures for ALL final candidates (universal, headless).

For every compound that reached the final list (07_final_candidates.csv)
generates, under supplementary/<target>/figures/ (the same folder where
08_supplementary.py puts journal-ready files):

  overlay/Fig_overlay_CLEAN.png      all CLEAN candidates overlaid in the site
  overlay/Fig_overlay_CLEAN.pse      editable PyMOL session of the overlay
  overlay/Fig_overlay_FLAGGED.png    all FLAGGED candidates overlaid
  overlay/Fig_overlay_FLAGGED.pse    editable PyMOL session of the overlay
  redock/Fig_redock.png              crystal ligand vs top redock pose
  per_ligand/CLEAN/<name>.png        3D: binding-site sticks + polar contacts
  per_ligand/FLAGGED/<name>.png      3D: same for flagged (liability) compounds
  interactions_2d/CLEAN/<name>.png   2D structure + PLIP interaction summary
  interactions_2d/FLAGGED/<name>.png 2D: same for flagged compounds

Writes FIGURES_DONE.txt at the end - the resume sentinel used by run_all.py
(stages 09 and 08 are considered done only when it exists).

Usage:
  python scripts/09_make_figures.py --target 7gqu
  python scripts/09_make_figures.py --target 4cgz          # any target
  python scripts/09_make_figures.py --target 7gqu --skip-2d  # 3D only (no PLIP)
"""
import argparse, csv, sys
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]

def read_rows(p):
    d = '\t' if str(p).endswith('.tsv') else ','
    with open(p, newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter=d))

def safe(s):
    return ''.join(c if c.isalnum() else '_' for c in s)

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
            print(f'[warn] no pose for {cid}; 3D skipped'); continue
        try:
            cmd.delete('all')
            load_protein()
            cmd.load(str(pose), 'lig')
            cmd.show('sticks', 'lig')
            cmd.color(CLS_COL.get(cls, 'cyan'), 'lig'); util.cnc('lig')
            cmd.select('site', 'byres (rec) within 5.0 of lig')
            cmd.show('sticks', 'site'); cmd.color('gray60', 'site'); util.cnc('site')
            cmd.distance('polar', 'lig and not elem H', 'site and not elem H', 3.5, mode=2)
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
            if pose is None: continue
            try:
                obj = f'lig_{i}'
                cmd.load(str(pose), obj)
                cmd.show('sticks', obj)
                cmd.color(PAL[i % len(PAL)], obj); util.cnc(obj)
                loaded.append(obj)
            except Exception as e:
                print(f'[warn] overlay load failed for {r["chembl_id"]}: {type(e).__name__}')
        if loaded:
            cmd.zoom(' or '.join(loaded), 6)
            save(fig / 'overlay' / f'Fig_overlay_{tag}.png', 1600, 1200)
            cmd.save(str(fig / 'overlay' / f'Fig_overlay_{tag}.pse'))
            print(f'[ok] overlay {tag}: {len(loaded)} ligands')

    # ---------------- redock overlay ----------------
    REF = BASE / 'data' / t / 'ref_ligand.sdf'
    RED = BASE / 'validation' / 'redock' / t / 'redock_exh32.sdf'
    if not RED.exists():
        RED = BASE / 'validation' / 'redock' / 'redock_exh32.sdf'
    if REF.exists() and RED.exists():
        cmd.delete('all')
        load_protein()
        cmd.load(str(REF), 'cryst'); cmd.show('sticks', 'cryst')
        cmd.color('green', 'cryst'); util.cnc('cryst')
        cmd.load(str(RED), 'redock'); cmd.show('sticks', 'redock')  # state 1 = top pose
        cmd.color('magenta', 'redock'); util.cnc('redock')
        cmd.zoom('cryst', 6)
        save(fig / 'redock' / 'Fig_redock.png', 1600, 1200)
        print('[ok] Fig_redock.png')
    else:
        print('[warn] redock files missing; Fig_redock skipped')

    # ---------------- 2D + PLIP ----------------
    n2d = 0
    if not a.skip_2d:
        try:
            from plip.structure.preparation import PDBComplex
            PLIP_OK = True
        except Exception as e:
            PLIP_OK = False
            print(f'[warn] PLIP unavailable ({e}); 2D without interactions')
        from rdkit import Chem
        from rdkit.Chem import Draw
        from PIL import Image, ImageDraw, ImageFont
        try:
            import matplotlib
            font = ImageFont.truetype(
                str(Path(matplotlib.get_data_path()) / 'fonts' / 'ttf' / 'DejaVuSans.ttf'), 20)
        except Exception:
            font = ImageFont.load_default()

        prot_lines = [l for l in REC.read_text().splitlines()
                      if l.startswith(('ATOM', 'TER'))]
        tmp = fig.parent / '_tmp_complex.pdb'
        for r in rows:
            cid, name, cls = r['chembl_id'], r['name'], r['cls']
            try:
                pose = pose_of(cid)
                mol = Chem.MolFromMolFile(str(pose), removeHs=True) if pose else None
                smi = lib.get(cid, {}).get('smiles_clean', '')
                if mol is None and smi:
                    mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    print(f'[warn] no structure for {cid}; 2D skipped'); continue

                inter = 'PLIP: n/a'
                if PLIP_OK and pose:
                    try:
                        lig_block = Chem.MolToPDBBlock(
                            Chem.MolFromMolFile(str(pose), removeHs=False))
                        tmp.write_text('\n'.join(prot_lines) + '\n' + lig_block + '\nEND\n')
                        p = PDBComplex(); p.load_pdb(str(tmp)); p.analyze()
                        hb = hy = pi = sa = ha = wb = 0
                        for site in p.interaction_sets.values():
                            n = lambda at: len(getattr(site, at, []) or [])
                            hb += n('hb_pdon') + n('hb_ldon'); hy += n('hydrophobic')
                            pi += n('pistacking'); ha += n('halogen'); wb += n('water_bridges')
                            sa += n('saltbridge_lpos') + n('saltbridge_pneg')
                        inter = (f'PLIP: H-bonds {hb} | hydrophobic {hy} | pi-stack {pi} | '
                                 f'salt-bridges {sa} | halogen {ha} | water-bridges {wb}')
                    except Exception as e:
                        inter = f'PLIP: error ({type(e).__name__})'

                img = Draw.MolToImage(mol, size=(900, 560))
                canvas = Image.new('RGB', (900, 780), 'white')
                canvas.paste(img, (0, 0))
                d = ImageDraw.Draw(canvas)
                d.text((12, 570), f'{name}  ({cid})  [{cls}]', fill='black', font=font)
                d.text((12, 600),
                       f"GNINA {r['gnina32']} | Vina {r['vina']} | LeDock {r['ledock']} "
                       f"kcal/mol; QED {r['qed']}", fill='black', font=font)
                d.text((12, 630), inter, fill='black', font=font)
                if r['flags']:
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
        f'3D={n3d}; 2D={n2d}; overlays=CLEAN,FLAGGED; redock=ok\n')

    print(f'\n✅ Figures: {fig}')
    print(f'   3D per-ligand: {n3d}; overlays: CLEAN+FLAGGED; redock; 2D: {n2d}')

if __name__ == '__main__':
    main()