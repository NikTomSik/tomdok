#!/usr/bin/env python3
"""01_redock_validation.py — PROTOCOL GATE: self-docking RMSD + LeDock native baseline.

What it does (runs immediately after 00_prepare_target, before any GPU screen):
  1. Loads the crystallised reference ligand (ref_ligand.sdf, fallback
     ref_ligand_prot.sdf) with robust sanitization.
  2. Redocks it with GNINA (20 poses) into the prepared receptor/box.
  3. Computes top-pose and best-pose RMSD vs the crystal pose
     (symmetry-aware GetBestRMS -> AlignMol -> MCS-mapped AlignMol fallback).
  4. Runs ONE LeDock docking of the native ligand and records its best score
     as `ledock_native` — the per-target baseline that makes LeDock scores of
     candidates interpretable (dLeDock = ledock(candidate) - ledock_native).
  5. Writes validation/redock/<target>/redock_report.txt and saves
     redock_poses.sdf for the figures stage.

Gate semantics:
  exit 0  verdict PASS  (best-pose RMSD < 2.0 A by default; --basis top for top-pose)
  exit 3  verdict FAIL  (pipeline stops here unless run_all --keep-going)
  exit 1  fatal (missing receptor/box/reference, GNINA failure)

Usage:
  python scripts/01_redock_validation.py --target 9fza
  python scripts/01_redock_validation.py --target 6vei --basis top
"""
import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
import warnings
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.WARNING)   # только ERROR и FATAL
warnings.filterwarnings('ignore', message='.*molecule is tagged as 2D.*')

from rdkit import Chem
from rdkit.Chem import rdFMCS, rdMolAlign

BASE = Path(__file__).resolve().parents[1]
RMSD_PASS = 2.0          # community standard for self-docking validation
NUM_MODES = 20
SEED = 42


# --------------------------------------------------------------------------
# reference handling
# --------------------------------------------------------------------------
def load_reference(data):
    """Load the crystallographic reference; ref_ligand.sdf first, then prot."""
    for name in ('ref_ligand.sdf', 'ref_ligand_prot.sdf'):
        p = data / name
        if not p.exists():
            continue
        for m in Chem.SDMolSupplier(str(p), removeHs=True):
            if m is not None and m.GetNumConformers() > 0:
                return m, name
    return None, None


def norm(m):
    """Unified perception: sanitization + aromaticity + stereo, so that the
    GNINA output poses and the reference always have isomorphic graphs."""
    m = Chem.Mol(m)
    try:
        Chem.SanitizeMol(m)
    except Exception:
        m.UpdatePropertyCache(strict=False)
        try:
            Chem.SanitizeMol(m)
        except Exception:
            pass
    Chem.SetAromaticity(m)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    return m


def rmsd_vs_ref(pose, ref):
    """Symmetry-aware RMSD with two fallbacks (7kcc-style graph mismatches)."""
    try:
        return rdMolAlign.GetBestRMS(pose, ref)
    except Exception:
        pass
    try:
        return rdMolAlign.AlignMol(pose, ref)
    except Exception:
        pass
    mcs = rdFMCS.FindMCS((ref, pose), timeout=5,
                         ringMatchesRingOnly=True, completeRingsOnly=True)
    patt = Chem.MolFromSmarts(mcs.smartsString)
    mr, mp = ref.GetSubstructMatch(patt), pose.GetSubstructMatch(patt)
    if not (mr and mp):
        raise RuntimeError('no common substructure between pose and reference')
    return rdMolAlign.AlignMol(pose, ref, atomMap=list(zip(mp, mr)))


# --------------------------------------------------------------------------
# engines
# --------------------------------------------------------------------------
def run_gnina_redock(receptor, ref_sdf, box, out_sdf, cpu):
    cmd = ['gnina', '-r', str(receptor), '-l', str(ref_sdf),
           '--center_x', str(box['center_x']), '--center_y', str(box['center_y']),
           '--center_z', str(box['center_z']),
           '--size_x', str(max(25.0, box['size_x'])),
           '--size_y', str(max(25.0, box['size_y'])),
           '--size_z', str(max(25.0, box['size_z'])),
           '--exhaustiveness', '8', '--num_modes', str(NUM_MODES),
           '--seed', str(SEED), '--cpu', str(cpu), '-o', str(out_sdf)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if r.returncode != 0 or not out_sdf.exists():
        sys.exit('ERROR gnina redock failed: ' + (r.stderr or r.stdout)[-300:])


def ledock_native_score(data, work, ref_sdf, box):
    """One LeDock docking of the native ligand -> best score (baseline).
    Never raises: returns None if LeDock is unavailable or parsing fails."""
    if not (shutil.which('ledock') and shutil.which('lepro')):
        return None
    rec_src = data / 'protein_clean.pdb'
    if not rec_src.exists():
        return None
    try:
        work.mkdir(parents=True, exist_ok=True)
        rec = work / 'rec.pdb'
        shutil.copy(rec_src, rec)
        subprocess.run(['lepro', str(rec)], cwd=work,
                       capture_output=True, text=True, timeout=600)
        lepro_rec = work / 'pro.pdb'
        if not lepro_rec.exists() or lepro_rec.stat().st_size == 0:
            pdbs = sorted(work.glob('*.pdb'),
                          key=lambda p: p.stat().st_mtime, reverse=True)
            if not pdbs: return None
            lepro_rec = pdbs[0]

        mol2 = work / 'native.mol2'
        r = subprocess.run(['obabel', str(ref_sdf), '-O', str(mol2)],
                           capture_output=True, text=True)
        if r.returncode != 0 or not mol2.exists():
            return None

        # LeDock requires a list file for ligands
        lst = work / 'native.list'
        lst.write_text(f"{mol2}\n")

        xmin = box['center_x'] - box['size_x'] / 2
        xmax = box['center_x'] + box['size_x'] / 2
        ymin = box['center_y'] - box['size_y'] / 2
        ymax = box['center_y'] + box['size_y'] / 2
        zmin = box['center_z'] - box['size_z'] / 2
        zmax = box['center_z'] + box['size_z'] / 2
        
        # EXACT LeDock config syntax (Receptor, RMSD, Binding pocket, 
        # Number of binding poses, Ligands list, END)
        cfg = work / 'dock.in'
        cfg.write_text(
            f"Receptor\n{lepro_rec}\n\n"
            f"RMSD\n1.0\n\n"
            f"Binding pocket\n"
            f"{xmin:.3f} {xmax:.3f}\n"
            f"{ymin:.3f} {ymax:.3f}\n"
            f"{zmin:.3f} {zmax:.3f}\n\n"
            f"Number of binding poses\n10\n\n"
            f"Ligands list\n{lst}\n\n"
            f"END\n"
        )
        
        dok = work / 'native.dok'
        if dok.exists(): dok.unlink()
        
        r = subprocess.run(['ledock', str(cfg)], cwd=work,
                           capture_output=True, text=True, timeout=1800)
        
        # LeDock writes output to <ligand_stem>.dok (i.e. native.dok)
        dok_out = work / 'native.dok'
        if not dok_out.exists():
            doks = sorted(work.glob('*.dok'))
            dok_out = doks[0] if doks else None
            
        if dok_out is None or not dok_out.exists():
            return None
            
        txt = dok_out.read_text(errors='ignore')
        scores = []
        # LeDock .dok format typically has "REMARK  Score: -X.XX" or similar
        for pat in (r'score\s*[:=]?\s*(-?\d+\.\d+)',
                    r'remark[^\n]*?(-?\d+\.\d+)',
                    r'energy\s*[:=]?\s*(-?\d+\.\d+)'):
            scores.extend([float(x) for x in re.findall(pat, txt, re.I)])
            
        if not scores:
            # fallback: grab all negative floats from REMARK lines
            for line in txt.splitlines():
                if line.upper().startswith('REMARK'):
                    scores.extend([float(x) for x in re.findall(r'-?\d+\.\d+', line)])
                    
        return min(scores) if scores else None
    except Exception as e:
        print(f"[warn] LeDock native baseline failed: {e}")
        return None


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description='Redock protocol gate + LeDock baseline')
    ap.add_argument('--target', required=True)
    ap.add_argument('--basis', choices=['best', 'top'], default='best',
                    help='RMSD criterion for the PASS verdict')
    ap.add_argument('--cpu', type=int, default=4)
    a = ap.parse_args()

    t = a.target.lower()
    data = BASE / 'data' / t
    vdir = BASE / 'validation' / 'redock' / t
    vdir.mkdir(parents=True, exist_ok=True)

    receptor = data / 'receptor.pdbqt'
    box_file = data / 'receptor.json'
    if not receptor.exists():
        sys.exit(f'ERROR: missing {receptor} (run 00_prepare_target first)')
    if not box_file.exists():
        sys.exit(f'ERROR: missing {box_file} (run 00_prepare_target first)')
    box = json.loads(box_file.read_text())

    ref, ref_name = load_reference(data)
    if ref is None:
        sys.exit(f'ERROR: no readable reference ligand in {data} '
                 f'(ref_ligand.sdf / ref_ligand_prot.sdf)')
    ref = norm(ref)
    print(f'Reference loaded: {ref_name} ({ref.GetNumAtoms()} atoms, '
          f'SMILES={Chem.MolToSmiles(ref)})')

    # ---- GNINA self-dock -------------------------------------------------
    poses_sdf = vdir / 'redock_poses.sdf'
    run_gnina_redock(receptor, data / ref_name, box, poses_sdf, a.cpu)
    poses = [norm(m) for m in Chem.SDMolSupplier(str(poses_sdf), removeHs=True) if m]
    if not poses:
        sys.exit('ERROR: no poses parsed from GNINA redock output')
    print(f'Poses generated : {len(poses)}')

    rmsds = [rmsd_vs_ref(p, ref) for p in poses]
    r_top, r_best = rmsds[0], min(rmsds)
    rank = next((i + 1 for i, x in enumerate(rmsds) if x <= RMSD_PASS), 0)
    print(f'Top-pose RMSD   : {r_top:.2f} Å')
    print(f'Best-pose RMSD  : {r_best:.2f} Å')

    # ---- LeDock native baseline (optional, never fatal) ------------------
    ledock_nat = ledock_native_score(data, vdir / 'ledock', data / ref_name, box)
    if ledock_nat is not None:
        print(f'LeDock native   : {ledock_nat:.2f}')
    else:
        print('LeDock native   : n/a (LeDock unavailable or parse failed)')

    # ---- verdict + report -------------------------------------------------
    crit = r_top if a.basis == 'top' else r_best
    ok = crit <= RMSD_PASS
    print(f'Validation      : {"PASS ✅" if ok else "FAIL ❌"} '
          f'(threshold {RMSD_PASS} Å, {a.basis} pose)')

    report = vdir / 'redock_report.txt'
    report.write_text(
        f'target={t}\n'
        f'reference={ref_name}\n'
        f'ref_smiles={Chem.MolToSmiles(ref)}\n'
        f'poses={len(poses)}\n'
        f'rmsd_top_pose={r_top:.3f}\n'
        f'rmsd_best_pose={r_best:.3f}\n'
        f'native_like_rank={rank}\n'
        f'ledock_native={ledock_nat if ledock_nat is not None else "n/a"}\n'
        f'cutoff={RMSD_PASS}\n'
        f'verdict={"PASS" if ok else "FAIL"}\n')
    print(f'Report: {report}')

    if not ok:
        sys.exit(3)   # gate: stop the pipeline before expensive stages


if __name__ == '__main__':
    main()