"""Final candidate prioritisation: 3-engine consensus + ADMET flags (no Cmax).
Adds phase (max clinical phase) and approved (YES/no).
CLEAN   = no safety flags
FLAGGED = hERG / withdrawn / PAINS / Lipinski>1
Universal: works for any target via --target.
Writes final_candidates.csv, prints both lists and glossary.

Usage:
 python scripts/07_final_candidates.py
 python scripts/07_final_candidates.py --target 4cgz
"""
import argparse, csv, os, sys
from pathlib import Path

BASE = Path(os.environ.get('SCREENING_ROOT', str(Path(__file__).resolve().parents[1])))

LEGEND = '''
=== Column glossary (final) ===
phase   - max clinical phase (4=approved, 3=PhIII, 2=PhII)
approved- regulatory approval (YES/no)
gnina32 - affinity GNINA exh=32, kcal/mol (more negative = stronger)
vina    - affinity Vina exh=32
ledock  - score LeDock (systematically softer)
n_le9   - number of engines (out of 3) with score <= -9 (consensus strength)
qed     - drug-likeness 0-1 (>0.5 good)
bbb     - BBB penetration (yes/no)
cls     - CLEAN (no safety flags) / FLAGGED (safety flags)
flags   - hERG / withdrawn / PAINS / Lip>1
'''

def read(p):
    with open(p, newline='') as f:
        return {r['chembl_id']: r for r in csv.DictReader(f)}

def load_phases():
    """Read max_phase from library.tsv (all SMILES that entered the library)."""
    ph = {}
    lib = BASE / 'data' / 'library' / 'library.tsv'
    if lib.exists():
        with open(lib, newline='', encoding='utf-8') as f:
            for r in csv.DictReader(f, delimiter='\t'):
                cid = (r.get('chembl_id') or '').strip()
                try: ph[cid] = int(r.get('max_phase') or 0)
                except (TypeError, ValueError): ph[cid] = 0
    return ph

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target', default='7gqu', help='PDB ID of the target')
    a = ap.parse_args()
    ref = BASE / 'results' / a.target.lower() / 'refine'

    for req in (ref / 'consensus.csv', ref / 'admet_all.csv'):
        if not req.exists():
            sys.exit(f'ERROR: missing {req} — run 04_refine_consensus.py and 06_admet.py first')

    cons = read(ref / 'consensus.csv')
    admet = read(ref / 'admet_all.csv')
    ph = load_phases()

    missing_admet = [cid for cid in cons if cid not in admet]
    if missing_admet:
        print(f'[Warning] {len(missing_admet)} consensus IDs absent from admet_all.csv '
              f'(kept with empty ADMET, treated as CLEAN)')

    rows = []
    for cid, c in cons.items():
        adm = admet.get(cid, {})
        flags = []
        if adm.get('withdrawn') == 'YES': flags.append('withdrawn')
        if adm.get('herg') == 'RISK': flags.append('hERG')
        if adm.get('pains') == 'ALERT': flags.append('PAINS')
        try:
            if int(adm.get('lipinski', 0)) > 1: flags.append('Lip>1')
        except (TypeError, ValueError): pass
        phase = ph.get(cid, -1)
        rows.append(dict(chembl_id=cid, name=c.get('name', ''), phase=phase,
                         approved='YES' if phase == 4 else 'no',
                         gnina32=c['gnina32_aff'], vina=c['vina_aff'], ledock=c['ledock_score'],
                         n_le9=c['n_engines_le9'], qed=adm.get('qed', ''), bbb=adm.get('bbb', ''),
                         cls='FLAGGED' if flags else 'CLEAN', flags=';'.join(flags)))
    rows.sort(key=lambda r: float(r['gnina32']))

    out = ref / 'final_candidates.csv'
    cols = ['chembl_id', 'name', 'phase', 'approved', 'gnina32', 'vina', 'ledock',
            'n_le9', 'qed', 'bbb', 'cls', 'flags']
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)

    clean = [r for r in rows if r['cls'] == 'CLEAN']
    flagged = [r for r in rows if r['cls'] == 'FLAGGED']
    print(f'=== CLEAN ({len(clean)}) — top by affinity ===')
    print(f'{"CHEMBL":<13}{"name":<16}{"ph":>3}{"app":>4}{"GNINA32":>8}{"VINA":>7}{"LeDock":>8}{"N<=-9":>6}{"QED":>6}{"BBB":>5}')
    for r in clean[:15]:
        print(f"{r['chembl_id']:<13}{r['name'][:15]:<16}{int(r['phase']):>3d}{r['approved']:>4}"
              f"{float(r['gnina32']):8.2f}{float(r['vina']):7.2f}"
              f"{float(r['ledock']):8.2f}{int(r['n_le9']):6d}{r['qed']:>6}{r['bbb']:>5}")
    print(f'\n=== FLAGGED ({len(flagged)}) — strong affinity but flagged ===')
    for r in flagged[:10]:
        print(f"{r['chembl_id']:<13}{r['name'][:15]:<16}{int(r['phase']):>3d}{r['approved']:>4}"
              f"{float(r['gnina32']):8.2f}  [{r['flags']}]")
    print(f'\n✅ Saved: {out}')
    print(LEGEND)

if __name__ == '__main__':
    main()