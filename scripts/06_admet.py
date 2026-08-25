"""ADMET profiling of consensus hits (RDKit + rule-based, UNIVERSAL).
Reads consensus.csv, takes top N by GNINA32 (or --all), computes:
  Phys-chem    : MW, LogP, HBA, HBD, TPSA, rotatable bonds, Lipinski, QED
  Absorption   : HIA (TPSA+MW rule)
  Distribution : BBB (TPSA+LogP rule)
  Toxicity     : hERG  - ChEMBL activity vs KCNH2/CHEMBL240 (data-driven)
                 PAINS - RDKit FilterCatalog
                 withdrawn - ChEMBL withdrawn_flag (from library.tsv)
No hardcoded drug names: safety flags come from ChEMBL data, so the script
is universal for any target / library. If library.tsv already carries a
'herg' column (from 01_build_library.py), it is used preferentially;
otherwise hERG is computed here directly from chembl_37.db.
Writes admet_<all|N>.csv, prints the table + column glossary.

Usage:
 python scripts/06_admet.py --top 15
 python scripts/06_admet.py --all
 python scripts/06_admet.py --input results/7gqu/refine/consensus.csv --top 20
"""
import argparse, csv, os, sqlite3
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, QED
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

BASE = Path(os.environ.get('SCREENING_ROOT', str(Path(__file__).resolve().parents[1])))
HERG_CUTOFF_NM = 10000  # 10 uM: hERG IC50/Ki below this => RISK

LEGEND = '''
=== Column glossary (ADMET) ===
MW    - molecular weight, Da (Lipinski: <=500)
LogP  - lipophilicity (Crippen); optimum 1-5, >5 = risk
HBA   - H-bond acceptors (<=10)
HBD   - H-bond donors (<=5)
TPSA  - polar surface area, A^2 (<=140 absorption, <=90 BBB)
RB    - rotatable bonds (flexibility; >10 worse)
Lip   - Lipinski violations (0-1 ok, >1 flag)
QED   - drug-likeness 0-1 (>0.5 good)
HIA   - intestinal absorption (high/low)
BBB   - BBB penetration (yes/no)
hERG  - hERG/QT block risk from ChEMBL activity vs KCNH2 (RISK/ok)
PAINS - pan-assay interference (ALERT/ok)
WD    - withdrawn per ChEMBL withdrawn_flag (YES/no)
'''

def load_lib():
    """chembl_id -> {smi, withdrawn, herg} from library.tsv (single source)."""
    lib = {}
    p = BASE / 'data' / 'library' / 'library.tsv'
    if p.exists():
        with open(p, encoding='utf-8') as f:
            for r in csv.DictReader(f, delimiter='\t'):
                c = (r.get('chembl_id') or '').strip()
                if c:
                    lib[c] = dict(
                        smi=(r.get('canonical_smiles') or r.get('smiles_clean') or '').strip(),
                        withdrawn=(r.get('withdrawn') or '').strip(),
                        herg=(r.get('herg') or '').strip())
    return lib

def herg_from_db():
    """Universal hERG liability from ChEMBL activities (KCNH2 / CHEMBL240)."""
    db = BASE / 'data' / 'chembl_37.db'
    out = {}
    if not db.exists():
        return out
    try:
        con = sqlite3.connect(db); cur = con.cursor()
        q = """
          SELECT md.chembl_id, MIN(act.standard_value)
          FROM activities act
          JOIN assays ass  ON ass.assay_id = act.assay_id
          JOIN target_dictionary td ON td.tid = ass.tid
          JOIN molecule_dictionary md ON md.molregno = act.molregno
          WHERE (td.chembl_id = 'CHEMBL240'
                 OR lower(td.pref_name) LIKE '%herg%'
                 OR lower(td.pref_name) LIKE '%kcnh2%')
            AND act.standard_type IN ('IC50','Ki','EC50')
            AND act.standard_units = 'nM'
            AND act.standard_value IS NOT NULL
          GROUP BY md.chembl_id
        """
        out = dict(cur.execute(q))
        con.close()
    except Exception:
        pass
    return out

def build_pains():
    p = FilterCatalogParams()
    try:
        p.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)   # new RDKit
    except AttributeError:
        try:
            p.AddCatalog(FilterCatalogParams.PAINS)              # old RDKit
        except AttributeError:
            return None                                          # no catalog -> skip
    return FilterCatalog(p)

def admet_row(cid, name, mol, pains, herg_flag, wdr_flag):
    mw = Descriptors.MolWt(mol); logp = Descriptors.MolLogP(mol)
    hba = Descriptors.NumHAcceptors(mol); hbd = Descriptors.NumHDonors(mol)
    tpsa = Descriptors.TPSA(mol); rb = Descriptors.NumRotatableBonds(mol)
    lip = sum([mw > 500, logp > 5, hba > 10, hbd > 5])
    qed = QED.qed(mol)
    hia = 'high' if (tpsa <= 140 and mw <= 500) else 'low'
    bbb = 'yes' if (tpsa <= 90 and 1.0 <= logp <= 5.0) else 'no'
    herg = herg_flag or 'ok'
    wdr = 'YES' if str(wdr_flag) in ('1', 'YES', 'Yes', 'yes', 'true', 'True') else 'no'
    pa = 'ALERT' if (pains is not None and pains.HasMatch(mol)) else 'ok'
    return dict(chembl_id=cid, name=name, mw=f'{mw:.1f}', logp=f'{logp:.2f}',
                hba=hba, hbd=hbd, tpsa=f'{tpsa:.0f}', rb=rb, lipinski=lip,
                qed=f'{qed:.2f}', hia=hia, bbb=bbb, herg=herg, pains=pa, withdrawn=wdr)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default=str(BASE / 'results' / '7gqu' / 'refine' / 'consensus.csv'))
    ap.add_argument('--top', type=int, default=15)
    ap.add_argument('--all', action='store_true', help='profile entire consensus')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()

    with open(a.input, newline='') as f:
        rows = [r for r in csv.DictReader(f) if r.get('gnina32_aff')]
    rows.sort(key=lambda r: float(r['gnina32_aff']))
    if not a.all: rows = rows[:a.top]

    lib = load_lib(); pains = build_pains()
    have_lib_herg = any(e['herg'] for e in lib.values())
    herg_db = {} if have_lib_herg else herg_from_db()

    out_rows = []
    for r in rows:
        cid = r['chembl_id']; e = lib.get(cid, {})
        s = e.get('smi', '')
        if not s: continue
        mol = Chem.MolFromSmiles(s)
        if mol is None: continue
        # hERG: prefer library column, else ChEMBL activity lookup
        hf = e.get('herg', '')
        if not hf:
            v = herg_db.get(cid)
            hf = 'RISK' if (v is not None and v < HERG_CUTOFF_NM) else 'ok'
        out_rows.append(admet_row(cid, r.get('name', ''), mol, pains, hf, e.get('withdrawn', '')))

    cols = ['chembl_id', 'name', 'mw', 'logp', 'hba', 'hbd', 'tpsa', 'rb',
            'lipinski', 'qed', 'hia', 'bbb', 'herg', 'pains', 'withdrawn']
    out = Path(a.out) if a.out else Path(a.input).parent / f'admet_{"all" if a.all else a.top}.csv'
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(out_rows)

    print(f'{"CHEMBL":<13}{"name":<16}{"MW":>7}{"LogP":>6}{"HBA":>4}{"HBD":>4}{"TPSA":>6}{"RB":>3}{"Lip":>4}{"QED":>6}{"HIA":>5}{"BBB":>4}{"hERG":>6}{"PAINS":>6}{"WD":>4}')
    for r in out_rows:
        print(f"{r['chembl_id']:<13}{r['name'][:15]:<16}{r['mw']:>7}{r['logp']:>6}{r['hba']:>4}{r['hbd']:>4}"
              f"{r['tpsa']:>6}{r['rb']:>3}{r['lipinski']:>4}{r['qed']:>6}{r['hia']:>5}{r['bbb']:>4}"
              f"{r['herg']:>6}{r['pains']:>6}{r['withdrawn']:>4}")
    print(f'\n✅ Saved: {out}')
    print(LEGEND)

if __name__ == '__main__':
    main()