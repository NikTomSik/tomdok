import os, csv, multiprocessing as mp
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog('rdApp.*')

BASE = Path(os.environ.get('SCREENING_ROOT', str(Path(__file__).resolve().parents[1])))
OUT = BASE / 'prep_v2'
TSVS = [BASE / 'data' / 'library' / 'library.tsv']
PDBQT_DIR = BASE / 'pdbqt'
SDF_DIR = BASE / 'sdf3d'
LOG = BASE / 'logs' / 'preparation_v3.log'
PH, SEED, NUM_CONFS = 7.4, 42, 3

NORM = rdMolStandardize.Normalizer()
TAUT = rdMolStandardize.TautomerEnumerator()
_norm = getattr(NORM, 'normalize', None) or getattr(NORM, 'Normalize')
_taut = getattr(TAUT, 'Canonicalize', None) or getattr(TAUT, 'canonicalize')

try:
    import inspect
    from dimorphite_dl import protonate_smiles
    _SIG = set(inspect.signature(protonate_smiles).parameters); HAS_DIM = True
except Exception:
    HAS_DIM, _SIG = False, set()

def protonate(smi):
    if not HAS_DIM: return smi
    kw = {}
    for k in ('min_ph', 'ph_min'):
        if k in _SIG: kw[k] = PH
    for k in ('max_ph', 'ph_max'):
        if k in _SIG: kw[k] = PH
    if 'max_variants' in _SIG: kw['max_variants'] = 1
    if 'label_states' in _SIG: kw['label_states'] = False
    try:
        st = protonate_smiles(smi, **kw)
        if st: return st[0]
    except Exception: pass
    return smi

def standardize(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None: return None, 'invalid SMILES'
    try:
        mol = _norm(mol)
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        Chem.SanitizeMol(mol)
        mol = _taut(mol)
    except Exception as e:
        return None, f'standardize: {e}'
    return mol, Chem.MolToSmiles(mol)

def make_3d(mol):
    n_ha = mol.GetNumHeavyAtoms()
    max_its = 200 if n_ha > 60 else 1000
    mol = Chem.AddHs(mol)
    cids = []
    seeds = (SEED, 7, 2026) if n_ha <= 60 else (SEED,)
    rands = (False, True) if n_ha <= 60 else (True,)
    n_confs = NUM_CONFS if n_ha <= 60 else 1
    for seed in seeds:
        for rand in rands:
            p = AllChem.ETKDGv3(); p.randomSeed = seed; p.useRandomCoords = rand
            try: ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=p)
            except Exception: ids = []
            if ids: cids = list(ids); break
        if cids: break
    if not cids: return None, 'embedding failed'
    best = None
    props = AllChem.MMFFGetMoleculeProperties(mol)
    if props is not None:
        for c in cids:
            try:
                ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=c)
                if ff is None: continue
                ff.Minimize(maxIts=max_its); e = ff.CalcEnergy()
            except Exception: continue
            if best is None or e < best[1]: best = (c, e)
    if best is None:
        for c in cids:
            try:
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=c)
                ff.Minimize(maxIts=max_its); e = ff.CalcEnergy()
            except Exception: continue
            if best is None or e < best[1]: best = (c, e)
    if best is None: return None, 'optimization failed'
    for c in list(cids):
        if c != best[0]: mol.RemoveConformer(c)
    return mol, 'ok'

def process(item):
    cid, canon = item
    out = PDBQT_DIR / f'{cid}.pdbqt'
    if out.exists(): return cid, 'SKIPPED', 'exists'
    smi_p = protonate(canon)
    mol = Chem.MolFromSmiles(smi_p) or Chem.MolFromSmiles(canon)
    if mol is None: return cid, 'ERROR', 'parse after protonation'
    if mol.GetNumHeavyAtoms() > 90: return cid, 'ERROR', 'too large (>90 HA)'
    m3, msg = make_3d(mol)
    if m3 is None: return cid, 'ERROR', msg
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        setups = MoleculePreparation().prepare(m3)
        s, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not ok: return cid, 'ERROR', f'meeko: {err}'
        out.write_text(s)
        w = Chem.SDWriter(str(SDF_DIR / f'{cid}.sdf')); w.write(m3); w.close()
        return cid, 'OK', smi_p
    except Exception as e:
        return cid, 'ERROR', f'meeko: {e}'

def main():
    for d in (PDBQT_DIR, SDF_DIR, LOG.parent): d.mkdir(parents=True, exist_ok=True)
    smi = {}
    for tsv in TSVS:
        if not tsv.exists():
            print(f'WARN: {tsv} not found'); continue
        with open(tsv, encoding='utf-8') as f:
            for r in csv.DictReader(f, delimiter='\t'):
                c = (r.get('chembl_id') or r.get('CHEMBL_ID') or '').strip()
                s = (r.get('canonical_smiles') or r.get('smiles') or '').strip()
                if c and s and c not in smi: smi[c] = s
    print(f'Merging lists: {len(smi)} unique IDs')

    seen, rows, n_filt = set(), [], 0
    for cid, s in smi.items():
        mol, canon = standardize(s)
        if mol is None or canon in seen: n_filt += 1; continue
        seen.add(canon); rows.append((cid, canon))
    print(f'After standardize+dedup: {len(rows)} (filtered out {n_filt})')

    ok = fail = skip = 0
    with open(LOG, 'w') as log, mp.Pool(min(8, os.cpu_count() or 8)) as pool:
        log.write('chembl_id\tstatus\tdetails\n')
        for i, (cid, st, det) in enumerate(pool.imap_unordered(process, rows), 1):
            log.write(f'{cid}\t{st}\t{det}\n')
            if st == 'OK': ok += 1
            elif st == 'SKIPPED': skip += 1
            else: fail += 1
            if i % 50 == 0: print(f'  {i}/{len(rows)} | ok={ok} skip={skip} fail={fail}')
    print(f'\nDone: ok={ok} skip={skip} fail={fail}')
    print(f'PDBQT total: {len(list(PDBQT_DIR.glob("*.pdbqt")))}')

if __name__ == '__main__':
    main()