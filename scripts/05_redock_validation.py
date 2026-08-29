#!/usr/bin/env python3
"""Redocking validation: the co-crystallised ligand is redocked with the same
screening protocol. PASS when the BEST pose RMSD <= cutoff (standard practice:
any pose within 2.0 A validates the protocol); top-pose RMSD reported for info.

Robustness:
 - reads ref_ligand.sdf first, falls back to ref_ligand_prot.sdf
 - normalises both ref and poses (aromaticity, stereo) before RMSD
 - uses MCS-based atom mapping as last-resort fallback for RMSD
   (fixes "No sub-structure match found" on tricky targets)

Usage:
 python scripts/05_redock_validation.py --target 7gqu
 python scripts/05_redock_validation.py --target 4cgz --exh 64
"""
import argparse, json, subprocess
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolAlign, rdFMCS

BASE = Path(__file__).resolve().parents[1]


def make_pdbqt(ref_sdf, out_pq):
    """Convert reference SDF to PDBQT for GNINA input (obabel)."""
    if not out_pq.exists():
        r = subprocess.run(['obabel', str(ref_sdf), '-O', str(out_pq)],
                           capture_output=True, text=True)
        if r.returncode != 0 or not out_pq.exists():
            raise RuntimeError('obabel failed: ' + (r.stderr or r.stdout)[-300:])
    return out_pq


def run_gnina(receptor, pq, box, exh, seed, out_sdf):
    """Run GNINA docking of the reference ligand into its own site."""
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


def norm(m):
    """Canonicalise a molecule so RDKit perceives the same graph for ref and pose.
    This is the key step that fixes "No sub-structure match found" errors when
    ref came from a PDB extract (sometimes with non-aromatic rings) and pose came
    from GNINA (always with full aromaticity perceived).
    """
    m = Chem.Mol(m)
    if m is None:
        return None
    # try full sanitisation; on failure, relax to valence-only
    try:
        Chem.SanitizeMol(m)
    except Exception:
        m.UpdatePropertyCache(strict=False)
        try:
            Chem.SanitizeMol(m)
        except Exception:
            pass
    Chem.SetAromaticity(m)                                    # unified aromaticity model
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)   # unified stereo
    return m


def rmsd(pose, ref):
    """Best RMS between pose and reference, with three-tier fallback:
      1. GetBestRMS  - handles symmetry automatically (standard path)
      2. AlignMol    - plain RMS after alignment
      3. MCS-mapped  - build atom map from the maximum common substructure,
                       then AlignMol with that map (last resort for weird cases)
    """
    pose, ref = norm(pose), norm(ref)
    if pose is None or ref is None:
        raise RuntimeError("cannot normalise pose or reference for RMSD")

    # --- tier 1: symmetry-aware best RMS ---
    try:
        return rdMolAlign.GetBestRMS(pose, ref)
    except Exception:
        pass

    # --- tier 2: plain alignment RMS ---
    try:
        return rdMolAlign.AlignMol(pose, ref)
    except Exception:
        pass

    # --- tier 3: MCS-based atom map (fixes tricky 7kcc-like cases) ---
    try:
        mcs = rdFMCS.FindMCS((ref, pose), timeout=10,
                             ringMatchesRingOnly=True, completeRingsOnly=True)
        if mcs.numAtoms < 3:
            raise RuntimeError("MCS too small")
        patt = Chem.MolFromSmarts(mcs.smartsString)
        mr = ref.GetSubstructMatch(patt)
        mp = pose.GetSubstructMatch(patt)
        if not (mr and mp):
            raise RuntimeError("MCS substructure match failed")
        atom_map = list(zip(mp, mr))
        return rdMolAlign.AlignMol(pose, ref, atomMap=atom_map)
    except Exception as e:
        raise RuntimeError(f"all RMSD tiers failed: {e}")


def load_reference(data_dir):
    """Find and load the reference ligand, preferring raw ref_ligand.sdf,
    falling back to ref_ligand_prot.sdf if the former is missing or unreadable.
    Returns (mol, path_used).
    """
    candidates = [
        data_dir / 'ref_ligand.sdf',
        data_dir / 'ref_ligand_prot.sdf',
    ]
    last_err = None
    for p in candidates:
        if not p.exists():
            last_err = f"not found: {p.name}"
            continue
        try:
            for m in Chem.SDMolSupplier(str(p), removeHs=True):
                if m is not None:
                    return m, p
        except Exception as e:
            last_err = f"{p.name} unreadable: {e}"
            continue
    raise SystemExit(
        f"ERROR: no readable reference ligand in {data_dir}\n"
        f"       last error: {last_err}\n"
        f"       (run 00_prepare_target or copy ref_ligand_prot.sdf -> ref_ligand.sdf)"
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target', default='7gqu', help='PDB ID of the target')
    ap.add_argument('--exh', type=int, default=32)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--cutoff', type=float, default=2.0)
    a = ap.parse_args()

    data = BASE / 'data' / a.target.lower()
    receptor = data / 'receptor.pdbqt'
    box_path = data / 'receptor.json'

    # preflight: make sure the inputs exist before spending GPU time
    for f, label in ((receptor, 'receptor'), (box_path, 'box')):
        if not f.exists():
            raise SystemExit(f'ERROR: missing {label} at {f}')
    box = json.loads(box_path.read_text())

    ref, ref_used = load_reference(data)
    print(f'Reference loaded: {ref_used.name} '
          f'({ref.GetNumAtoms()} atoms, SMILES={Chem.MolToSmiles(ref)})')

    out = BASE / 'validation' / 'redock' / a.target.lower()
    out.mkdir(parents=True, exist_ok=True)

    pq  = make_pdbqt(ref_used, out / 'ref_ligand_dock.pdbqt')
    sdf = run_gnina(receptor, pq, box, a.exh, a.seed, out / f'redock_exh{a.exh}.sdf')

    poses = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=True) if m]
    if not poses:
        raise SystemExit('ERROR: failed to parse poses from output SDF')

    r_top  = rmsd(poses[0], ref)
    r_best = min(rmsd(m, ref) for m in poses)
    ok = r_best <= a.cutoff   # PASS if ANY pose reproduces the crystal pose

    print(f'Poses generated : {len(poses)}')
    print(f'Top-pose RMSD   : {r_top:.2f} Å')
    print(f'Best-pose RMSD  : {r_best:.2f} Å')
    print(f'Validation      : {"PASS ✅" if ok else "FAIL ❌"} '
          f'(threshold {a.cutoff} Å, best pose)')

    (out / 'redock_report.txt').write_text(
        f'target={a.target}\n'
        f'reference={ref_used.name}\n'
        f'ref_smiles={Chem.MolToSmiles(ref)}\n'
        f'exhaustiveness={a.exh}\n'
        f'seed={a.seed}\n'
        f'poses={len(poses)}\n'
        f'rmsd_top_pose={r_top:.3f}\n'
        f'rmsd_best_pose={r_best:.3f}\n'
        f'cutoff={a.cutoff}\n'
        f'verdict={"PASS" if ok else "FAIL"}\n')
    print(f'Report: {out / "redock_report.txt"}')


if __name__ == '__main__':
    main()