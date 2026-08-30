"""Unified target preparation for virtual screening (robust version).

Fully automatic cycle for a NEW target (single command):
  python scripts/00_prepare_target.py --pdb-id 4cgz
  -> downloads 4CGZ.pdb from RCSB
  -> extracts crystallographic ligand (largest hetero-residue) to ref_ligand.sdf
  -> prepares protein, receptors and docking box under data/<pdb-id>/

Reference-ligand selection is TWO-TIER:
  * tier 1 (preferred): drug-like hetero-residues (not in SKIP_RES / COFACTOR_RES);
  * tier 2 (fallback):  cofactors / coenzymes / nucleotides (COFACTOR_RES),
    used only when no drug-like residue is present in the structure.
Bond perception: OpenBabel first (robust on PDB extracts, keeps crystal
coordinates), RDKit with permissive sanitization as fallback.

Steps:
 0. Download PDB from RCSB (if missing) + ensure ref_ligand.sdf exists in output dir
 1. Protein: PDBFixer - waters/cofactors, missing atoms, H at pH 7.4 -> protein_clean.pdb
 2. Receptor: PDBQT with Gasteiger charges (OpenBabel)               -> receptor.pdbqt
 2b. Vina receptor: without ROOT/BRANCH/TORSDOF tags                 -> receptor_vina.pdbqt
 3. Reference ligand: Dimorphite-DL pH 7.4, ETKDGv3 (seed 42), MMFF94,
    Meeko with OpenBabel fallback                                    -> ref_ligand_prot.sdf, ref_ligand.pdbqt
 4. Box: center = geometric center of crystal ligand,
    size = clip(ptp + 5 A, min_size, 27 A)                           -> receptor.json

Usage:
  python scripts/00_prepare_target.py --pdb-id 7gqu
  python scripts/00_prepare_target.py --pdb-id 6vei --reextract      # force re-extraction
  python scripts/00_prepare_target.py --pdb my.pdb --ref my_ligand.sdf --out data/my
"""
import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

PH = 7.4
RANDOM_SEED = 42
FLEX_TAGS = ('ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH', 'TORSDOF')

# waters, ions, buffers, cryo-additives, modified residues, sugars - NOT ligands
SKIP_RES = {'HOH', 'WAT', 'TIP', 'DOD', 'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE', 'MN',
            'CU', 'CO', 'NI', 'LI', 'RB', 'CS', 'SO4', 'PO4', 'EDO', 'GOL', 'ACT',
            'CIT', 'TRS', 'HEP', 'MES', 'ACE', 'NH4', 'CO3', 'SCN', 'FMT', 'IMD',
            'MPD', 'PEG', 'PG4', 'EPE', 'BME', 'DTT', 'IPA', 'BU3', 'TBU', 'IOD',
            'MSE', 'SEC', 'CSO', 'HYP', 'SEP', 'TPO', 'NAG', 'MAN', 'BMA', 'FUC'}

# cofactors / coenzymes / nucleotides - excluded from the drug-like pool;
# used ONLY as fallback when no drug-like residue is present (two-tier selection)
COFACTOR_RES = {'NDP', 'NAD', 'NAP', 'NAI', 'ATP', 'ADP', 'AMP', 'ANP', 'AGS',
                'GTP', 'GDP', 'GMP', 'GNP', 'UTP', 'UDP', 'CTP', 'CDP', 'CMP',
                'FAD', 'FMN', 'HEM', 'HEA', 'HEME', 'COA', 'ACO', 'SAM', 'SAH',
                'SIN', 'PLP', 'TPP', 'B12', 'FOL'}


def download_pdb(pdb_id, dest):
    """Downloading structure from RCSB (with mmCIF fallback for new structures)."""
    import urllib.request
    import urllib.error
    
    pdb_id_upper = pdb_id.upper()
    url = f'https://files.rcsb.org/download/{pdb_id_upper}.pdb'
    print(f'[0/4] Downloading {url} ...')
    
    try:
        urllib.request.urlretrieve(url, dest)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f'      [!] PDB format not available (404). RCSB only provides mmCIF for new entries.')
            cif_url = f'https://files.rcsb.org/download/{pdb_id_upper}.cif'
            cif_dest = dest.with_suffix('.cif')
            print(f'      Downloading {cif_url} ...')
            urllib.request.urlretrieve(cif_url, cif_dest)
            
            print(f'      Converting mmCIF to PDB using OpenMM...')
            from openmm.app import PDBxFile, PDBFile
            cif_obj = PDBxFile(str(cif_dest))
            with open(dest, 'w') as f:
                PDBFile.writeFile(cif_obj.topology, cif_obj.positions, f, keepIds=True)
            
            cif_dest.unlink(missing_ok=True) # удаляем временный .cif
            print(f'      Converted and saved: {dest.name}')
        else:
            raise


def extract_ref_ligand(pdb_path, out_sdf):
    """Extracts the largest hetero-residue from PDB -> SDF with crystallographic
    coordinates. Two-tier selection: drug-like residues preferred; cofactors
    used only as fallback. Returns atom count or None."""
    drug, cofac = {}, {}
    for line in Path(pdb_path).read_text().splitlines():
        if not line.startswith('HETATM'):
            continue
        resn = line[17:20].strip()
        if not resn or resn in SKIP_RES:
            continue
        key = (line[21:22], line[22:26].strip(), resn)
        (cofac if resn in COFACTOR_RES else drug).setdefault(key, []).append(line)

    if drug:
        pool, tier = drug, 'drug-like'
    elif cofac:
        pool, tier = cofac, 'cofactor (fallback)'
    else:
        return None

    best = max(pool.values(), key=len)          # largest residue = putative ligand
    resn = next(iter(best))[17:20].strip()
    print(f'      selected residue {resn} ({len(best)} atoms, tier: {tier})')

    tmp = out_sdf.with_suffix('.tmp.pdb')
    tmp.write_text('\n'.join(best) + '\nEND\n')

    # --- OpenBabel first: robust bond perception, keeps crystal coordinates ---
    n_atoms = None
    if shutil.which('obabel'):
        r = subprocess.run(['obabel', str(tmp), '-O', str(out_sdf)],
                           capture_output=True, text=True)
        if r.returncode == 0 and out_sdf.exists():
            m = next((x for x in Chem.SDMolSupplier(str(out_sdf), removeHs=False) if x), None)
            if m is not None and m.GetNumAtoms() >= 3:
                n_atoms = m.GetNumAtoms()

    # --- RDKit fallback with permissive sanitization ---
    if n_atoms is None:
        mol = Chem.MolFromPDBFile(str(tmp), removeHs=False, sanitize=False)
        if mol is None or mol.GetNumAtoms() < 3:
            tmp.unlink(missing_ok=True)
            return None
        try:
            mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(mol, catchErrors=True)
            Chem.SetAromaticity(mol)
        except Exception:
            pass
        w = Chem.SDWriter(str(out_sdf)); w.write(mol); w.close()
        m = next((x for x in Chem.SDMolSupplier(str(out_sdf), removeHs=False) if x), None)
        if m is None:
            tmp.unlink(missing_ok=True)
            return None
        n_atoms = m.GetNumAtoms()

    tmp.unlink(missing_ok=True)
    return n_atoms


def prep_protein(pdb_in, pdb_out):
    """Step 1: PDBFixer."""
    import pdbfixer
    from openmm.app import PDBFile
    print('[1/4] Protein preparation (PDBFixer, pH 7.4)...')
    fixer = pdbfixer.PDBFixer(filename=str(pdb_in))
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=PH)
    with open(pdb_out, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    print(f'      saved: {pdb_out.name}')


def make_receptor_pdbqt(pdb_in, pdbqt_out):
    """Step 2: receptor -> PDBQT with Gasteiger charges (GNINA/AD4/LeDock)."""
    print('[2/4] Receptor PDBQT (OpenBabel + Gasteiger)...')
    if not shutil.which('obabel'):
        sys.exit('ERROR: obabel not found (conda install -c conda-forge openbabel)')
    r = subprocess.run(['obabel', str(pdb_in), '-O', str(pdbqt_out),
                        '--partialcharge', 'gasteiger'],
                       capture_output=True, text=True)
    if r.returncode != 0 or not Path(pdbqt_out).exists():
        sys.exit('ERROR obabel: ' + (r.stderr or r.stdout)[-300:])
    print(f'      saved: {pdbqt_out.name}')


def make_vina_receptor(pdbqt_with_tags, pdbqt_clean):
    """Step 2b: strips flexible-ligand tags - Vina rejects them in rigid receptor."""
    lines = [l for l in Path(pdbqt_with_tags).read_text().splitlines()
             if not l.startswith(FLEX_TAGS)]
    Path(pdbqt_clean).write_text('\n'.join(lines) + '\n')
    print(f'      saved: {pdbqt_clean.name} (without {"/".join(FLEX_TAGS)})')


def prep_ref_ligand(sdf_in, sdf_out, pdbqt_out):
    """Step 3: protonation + 3D + Meeko (Meeko -> OpenBabel fallback)."""
    print('[3/4] Reference ligand (pH 7.4 + 3D + Meeko)...')
    mol = next((m for m in Chem.SDMolSupplier(str(sdf_in), removeHs=True) if m), None)
    if mol is None:
        sys.exit(f'ERROR: failed to read ligand from {sdf_in} '
                 f'(re-run with --reextract or provide --ref)')
    smiles = Chem.MolToSmiles(mol)
    print(f'      input SMILES: {smiles}')

    protonated = smiles
    try:
        from dimorphite_dl import protonate_smiles
        states = protonate_smiles(smiles, ph_min=PH, ph_max=PH,
                                  max_variants=1, label_states=False)
        protonated = states[0] if states else smiles
        print(f'      protonated SMILES: {protonated}')
    except Exception as e:
        print(f'      [warning] Dimorphite-DL failed: {e}, using input SMILES')

    mol = Chem.MolFromSmiles(protonated)
    if mol is None:
        sys.exit('ERROR: invalid protonated SMILES')

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = RANDOM_SEED
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=1, params=params)
    if conf_ids:
        AllChem.MMFFOptimizeMolecule(mol, confId=list(conf_ids)[0], maxIters=1000)
        w = Chem.SDWriter(str(sdf_out)); w.write(mol); w.close()
    else:
        # embedding failed (e.g. cofactor) - keep crystallographic coordinates
        print('      [warning] 3D embedding failed - keeping crystallographic coordinates')
        shutil.copy(sdf_in, sdf_out)
    print(f'      saved: {sdf_out.name}')

    # PDBQT: Meeko first, OpenBabel fallback (cofactors often break Gasteiger/Meeko)
    ok = False
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        preparator = MoleculePreparation(keep_chorded_rings=False,
                                         merge_these_atom_types=['H'],
                                         rigid_macrocycles=False)
        mol_setups = preparator.prepare(mol)
        pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(mol_setups[0])
        if is_ok:
            Path(pdbqt_out).write_text(pdbqt_string)
            ok = True
            print(f'      saved: {pdbqt_out.name} (Meeko)')
        else:
            print(f'      [warning] Meeko failed: {error_msg}')
    except Exception as e:
        print(f'      [warning] Meeko crashed: {e}')
    if not ok:
        r = subprocess.run(['obabel', str(sdf_out), '-O', str(pdbqt_out),
                            '--partialcharge', 'gasteiger'],
                           capture_output=True, text=True)
        if r.returncode != 0 or not Path(pdbqt_out).exists():
            sys.exit('ERROR: reference PDBQT failed via Meeko and obabel: '
                     + (r.stderr or r.stdout)[-300:])
        print(f'      saved: {pdbqt_out.name} (OpenBabel fallback)')


def calc_box(sdf_in, json_out, min_size=25.0):
    """Step 4: docking box from crystallographic coordinates."""
    print('[4/4] Docking box from crystallographic coordinates...')
    mol = next((m for m in Chem.SDMolSupplier(str(sdf_in), removeHs=False) if m), None)
    if mol is None:
        sys.exit(f'ERROR: failed to read ligand from {sdf_in}')
    conf = mol.GetConformer()
    pts = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    center = pts.mean(axis=0)
    size = np.clip(np.ptp(pts, axis=0) + 5.0, min_size, 27.0)
    config = {k: round(float(v), 3) for k, v in zip(
        ['center_x', 'center_y', 'center_z', 'size_x', 'size_y', 'size_z'],
        [*center, *size])}
    Path(json_out).write_text(json.dumps(config, indent=4))
    print('      ' + json.dumps(config))


def main():
    ap = argparse.ArgumentParser(description='Unified target preparation')
    ap.add_argument('--pdb-id', default='7gqu', help='PDB ID for downloading from RCSB')
    ap.add_argument('--pdb', default=None, help='input PDB (default <out>/<pdb-id>.pdb)')
    ap.add_argument('--ref', default=None,
                    help='reference ligand SDF (default <out>/ref_ligand.sdf; '
                         'if absent, extracted from PDB automatically)')
    ap.add_argument('--out', default=None, help='output folder (default: data/<pdb-id>)')
    ap.add_argument('--min-size', type=float, default=25.0,
                    help='minimum box size, A (guard against tight box)')
    ap.add_argument('--reextract', action='store_true',
                    help='force re-extraction of the reference ligand from the PDB '
                         '(ignores an existing ref_ligand.sdf)')
    a = ap.parse_args()

    out = Path(a.out) if a.out else Path(__file__).resolve().parents[1] / 'data' / a.pdb_id.lower()
    out.mkdir(parents=True, exist_ok=True)
    pdb = Path(a.pdb) if a.pdb else out / f'{a.pdb_id.lower()}.pdb'

    # raw_ref - файл, который stage 05 будет искать для redock (кристаллографические координаты)
    raw_ref = out / 'ref_ligand.sdf'
    user_ref = Path(a.ref) if a.ref else None

    if not pdb.exists():
        download_pdb(a.pdb_id, pdb)

    if user_ref:
        if not user_ref.exists():
            sys.exit(f'ERROR: provided --ref file not found: {user_ref}')
        if user_ref.resolve() != raw_ref.resolve():
            shutil.copy(user_ref, raw_ref)
            print(f'[0/4] Copied reference ligand to {raw_ref.name}')
        else:
            print(f'[0/4] Using provided reference ligand: {raw_ref.name}')
    elif a.reextract or not raw_ref.exists():
        if raw_ref.exists():
            raw_ref.unlink()
            print('[0/4] --reextract: removing old ref_ligand.sdf')
        print(f'[0/4] No reference ligand - extracting from {pdb.name}...')
        n = extract_ref_ligand(pdb, raw_ref)
        if n is None:
            sys.exit(f'ERROR: failed to extract ligand from {pdb} '
                     f'(manually: obabel {pdb} -O ref.sdf -m)')
        print(f'      saved: {raw_ref.name} ({n} atoms)')
    else:
        print(f'[0/4] Using existing reference ligand: {raw_ref.name}')

    prep_protein(pdb, out / 'protein_clean.pdb')
    make_receptor_pdbqt(out / 'protein_clean.pdb', out / 'receptor.pdbqt')
    make_vina_receptor(out / 'receptor.pdbqt', out / 'receptor_vina.pdbqt')
    prep_ref_ligand(raw_ref, out / 'ref_ligand_prot.sdf', out / 'ref_ligand.pdbqt')
    calc_box(raw_ref, out / 'receptor.json', a.min_size)

    # ---- self-check ----
    rec = (out / 'receptor.pdbqt').read_text().splitlines()
    rec_v = (out / 'receptor_vina.pdbqt').read_text().splitlines()
    n_tags = sum(1 for l in rec if l.startswith(FLEX_TAGS))
    n_tags_v = sum(1 for l in rec_v if l.startswith(FLEX_TAGS))
    charge = next((l.split()[-2] for l in rec if l.startswith('ATOM')), 'n/a')
    print(f'\n✅ Target ready: {out}')
    print(f'   receptor.pdbqt      : {len(rec)} lines, tags ROOT/BRANCH={n_tags}, '
          f'charge of first atom={charge} (GNINA/AD4/LeDock)')
    print(f'   receptor_vina.pdbqt : {len(rec_v)} lines, tags ROOT/BRANCH={n_tags_v} (Vina)')


if __name__ == '__main__':
    main()