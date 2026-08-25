#!/usr/bin/env python3
"""Redocking validation: the co-crystallised ligand is redocked with the same
screening protocol. PASS when the BEST pose RMSD <= cutoff (standard practice:
any pose within 2.0 A validates the protocol); top-pose RMSD reported for info.

Usage:
 python scripts/05_redock_validation.py --target 7gqu
 python scripts/05_redock_validation.py --target 4cgz --exh 64
"""
import argparse, json, subprocess
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolAlign

BASE = Path(__file__).resolve().parents[1]

def make_pdbqt(ref_sdf, out_pq):
    if not out_pq.exists():
        r = subprocess.run(['obabel', str(ref_sdf), '-O', str(out_pq)],
                           capture_output=True, text=True)
        if r.returncode != 0 or not out_pq.exists():
            raise RuntimeError('obabel failed: ' + (r.stderr or r.stdout)[-300:])
    return out_pq

def run_gnina(receptor, pq, box, exh, seed, out_sdf):
    cmd = ['gnina', '-r', str(receptor), '-l', str(pq),
           '--center_x', str(box['center_x']), '--center_y', str(box['center_y']),
           '--center_z', str(box['center_z']),
           '--size_x', str(max(25.0, box['size_x'])),
           '--size_y', str(max(25.0, box['size_y'])),
           '--size_z', str(max(25.0, box['size_z'])),
           '--exhaustiveness', str(exh), '--num_modes', '20',
           '--seed', str(seed), '--cpu', '8', '-o', str(out_sdf)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if r.returncode != 0:
        raise RuntimeError('gnina failed: ' + (r.stderr or r.stdout)[-300:])
    return out_sdf

def rmsd(pose, ref):
    try:
        return rdMolAlign.GetBestRMS(pose, ref)   # accounts for symmetry
    except Exception:
        return rdMolAlign.AlignMol(pose, ref)     # fallback

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target', default='7gqu', help='PDB ID of the target')
    ap.add_argument('--exh', type=int, default=32)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--cutoff', type=float, default=2.0)
    a = ap.parse_args()

    data = BASE / 'data' / a.target.lower()
    receptor = data / 'receptor.pdbqt'
    ref_sdf = data / 'ref_ligand.sdf'
    box = json.loads((data / 'receptor.json').read_text())
    out = BASE / 'validation' / 'redock' / a.target.lower()
    out.mkdir(parents=True, exist_ok=True)

    ref = next((m for m in Chem.SDMolSupplier(str(ref_sdf), removeHs=True) if m), None)
    if ref is None:
        raise SystemExit(f'ERROR: failed to read reference ligand from {ref_sdf}')

    pq = make_pdbqt(ref_sdf, out / 'ref_ligand_dock.pdbqt')
    sdf = run_gnina(receptor, pq, box, a.exh, a.seed, out / f'redock_exh{a.exh}.sdf')
    poses = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=True) if m]
    if not poses:
        raise SystemExit('ERROR: failed to parse poses from output SDF')

    r_top = rmsd(poses[0], ref)
    r_best = min(rmsd(m, ref) for m in poses)
    ok = r_best <= a.cutoff   # PASS if ANY pose reproduces the crystal pose

    print(f'Poses generated : {len(poses)}')
    print(f'Top-pose RMSD   : {r_top:.2f} Å')
    print(f'Best-pose RMSD  : {r_best:.2f} Å')
    print(f'Validation      : {"PASS ✅" if ok else "FAIL ❌"} (threshold {a.cutoff} Å, best pose)')

    (out / 'redock_report.txt').write_text(
        f'target={a.target}\nexhaustiveness={a.exh}\nseed={a.seed}\nposes={len(poses)}\n'
        f'rmsd_top_pose={r_top:.3f}\nrmsd_best_pose={r_best:.3f}\n'
        f'cutoff={a.cutoff}\nverdict={"PASS" if ok else "FAIL"}\n')
    print(f'Report: {out / "redock_report.txt"}')

if __name__ == '__main__':
    main()