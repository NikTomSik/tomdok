"""Unified library from ChEMBL SQLite (chembl_37.db).

SQL filters: max_phase >= --min-phase (parameterised), Small molecule, has SMILES,
not inorganic and not polymer (*_flag columns). Optionally (via flags):
--exclude-withdrawn, --exclude-veterinary.
Deduplication by parent (largest fragment, canonical SMILES).
Adds a data-driven 'herg' column: RISK if the compound shows hERG (KCNH2 /
CHEMBL240) IC50/Ki/EC50 < 10 uM in ChEMBL activities, else 'ok' ('' if no data).

Depth: --min-phase 3 = approved+PhIII (default, repurposing),
2 = +PhII, 0 = all. Output: data/library/library.tsv (columns stable
for 02/03/04/06/07).

Usage:
 python scripts/01_build_library.py --db data/chembl_37.db
 python scripts/01_build_library.py --min-phase 2 --out /tmp/lib_ph2.tsv  # preview
"""
import argparse, csv, sqlite3, sys
from pathlib import Path
from collections import Counter

BASE = Path(__file__).resolve().parents[1]
HERG_CUTOFF_NM = 10000  # 10 uM: hERG IC50/Ki below this => RISK

OUT_COLS = ['chembl_id', 'pref_name', 'max_phase', 'first_approval', 'molecule_type',
            'oral', 'parenteral', 'topical', 'natural_product', 'prodrug', 'inorganic',
            'polymer', 'veterinary', 'withdrawn', 'canonical_smiles', 'inchi_key',
            'mw', 'alogp', 'hba', 'hbd', 'psa', 'rtb', 'ro5_viol', 'qed',
            'heavy_atoms', 'aromatic_rings', 'smiles_clean', 'herg']


def table_cols(cur, t):
    return {r[1].lower() for r in cur.execute(f'PRAGMA table_info({t})')}


def clean_smiles(smi):
    """Parent: largest fragment, canonical SMILES (desalt)."""
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(smi)
        if m is None: return ''
        frags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
        m = max(frags, key=lambda x: x.GetNumHeavyAtoms())
        return Chem.MolToSmiles(m)
    except Exception:
        return ''


def herg_map(cur):
    """molregno -> min hERG standard_value (nM); {} if tables missing."""
    try:
        qh = """
          SELECT act.molregno, MIN(act.standard_value)
          FROM activities act
          JOIN assays ass  ON ass.assay_id = act.assay_id
          JOIN target_dictionary td ON td.tid = ass.tid
          WHERE (td.chembl_id = 'CHEMBL240'
                 OR lower(td.pref_name) LIKE '%herg%'
                 OR lower(td.pref_name) LIKE '%kcnh2%')
            AND act.standard_type IN ('IC50','Ki','EC50')
            AND act.standard_units = 'nM'
            AND act.standard_value IS NOT NULL
          GROUP BY act.molregno
        """
        return dict(cur.execute(qh)), True
    except Exception:
        return {}, False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--db', default=str(BASE / 'data' / 'chembl_37.db'))
    ap.add_argument('--out', default=None)
    ap.add_argument('--min-phase', type=int, default=3, choices=[0, 1, 2, 3, 4])
    ap.add_argument('--exclude-withdrawn', action='store_true',
                    help='exclude withdrawn (default: KEEP — they are valid hits)')
    ap.add_argument('--exclude-veterinary', action='store_true')
    ap.add_argument('--limit', type=int, default=0)
    a = ap.parse_args()

    db = Path(a.db)
    if not db.exists():
        sys.exit(f'Missing {db}.\nDownload: https://ftp.ebi.ac.uk/pub/databases/chembl/'
                 f'ChEMBLdb/releases/chembl_37/chembl_37_sqlite.tar.gz')
    out = Path(a.out) if a.out else db.parent / 'library' / 'library.tsv'
    out.parent.mkdir(parents=True, exist_ok=True)

    con = sqlite3.connect(db); cur = con.cursor()
    md = table_cols(cur, 'molecule_dictionary')
    cp = table_cols(cur, 'compound_properties')

    def mcol(n, *alts):
        for x in (n, *alts):
            if x in md: return f'md.{x}'
        return 'NULL'
    def ccol(n, *alts):
        for x in (n, *alts):
            if x in cp: return f'cp.{x}'
        return 'NULL'

    inorg, poly = mcol('inorganic_flag', 'inorganic'), mcol('polymer_flag', 'polymer')
    withdr, vet = mcol('withdrawn_flag', 'withdrawn'), mcol('veterinary')

    where = ['md.max_phase >= ?',
             'cs.canonical_smiles IS NOT NULL',
             "md.molecule_type = 'Small molecule'",
             f'COALESCE({inorg},0) = 0', f'COALESCE({poly},0) = 0']
    if a.exclude_withdrawn: where.append(f'COALESCE({withdr},0) = 0')
    if a.exclude_veterinary: where.append(f'COALESCE({vet},0) = 0')

    q = f"""
    SELECT md.molregno, md.chembl_id, md.pref_name, md.max_phase,
           {mcol('first_approval')} AS first_approval, md.molecule_type,
           {mcol('oral')},{mcol('parenteral')},{mcol('topical')},{mcol('natural_product')},
           {mcol('prodrug')},{inorg},{poly},{vet},{withdr},
           cs.canonical_smiles, cs.standard_inchi_key,
           {ccol('mw_freebase','full_mwt')},{ccol('alogp')},{ccol('hba')},{ccol('hbd')},
           {ccol('psa')},{ccol('rtb')},{ccol('ro5_violations')},{ccol('qed_weighted')},
           {ccol('heavy_atoms')},{ccol('aromatic_rings')}
    FROM molecule_dictionary md
    JOIN compound_structures cs ON cs.molregno=md.molregno
    LEFT JOIN compound_properties cp ON cp.molregno=md.molregno
    WHERE {' AND '.join(where)}
    """
    cur.execute(q, (a.min_phase,))
    names = [d[0] for d in cur.description]
    rows = [dict(zip(names, r)) for r in cur.fetchall()]
    print(f'From DB after SQL filters: {len(rows)}')
    if not rows:
        with open(out, 'w', newline='', encoding='utf-8') as f:
            csv.writer(f, delimiter='\t').writerow(OUT_COLS)
        print(f'No molecules matched; wrote empty {out}'); return

    # ---- data-driven hERG liability ----
    herg, herg_ok = herg_map(cur)
    for r in rows:
        if herg_ok:
            v = herg.get(r.get('molregno'))
            r['herg'] = 'RISK' if (v is not None and v < HERG_CUTOFF_NM) else 'ok'
        else:
            r['herg'] = ''

    # ---- parent + dedup ----
    for r in rows: r['smiles_clean'] = clean_smiles(r['canonical_smiles'])
    rows = [r for r in rows if r['smiles_clean']]
    rows.sort(key=lambda r: (r['smiles_clean'] != r['canonical_smiles'],
                             -int(r['max_phase'] or 0), r['chembl_id']))
    seen, ded = set(), []
    for r in rows:
        if r['smiles_clean'] in seen: continue
        seen.add(r['smiles_clean']); ded.append(r)
        if a.limit and len(ded) >= a.limit: break
    print(f'After parent-based deduplication: {len(ded)}')

    ded.sort(key=lambda r: r['chembl_id'])
    with open(out, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(OUT_COLS)
        for r in ded:
            w.writerow([r.get(c) if r.get(c) is not None else '' for c in OUT_COLS])

    c = Counter(int(r['max_phase'] or 0) for r in ded)
    n_herg = sum(1 for r in ded if r['herg'] == 'RISK')
    print(f'✅ Saved: {out} (total {len(ded)}, hERG RISK={n_herg})')
    for ph in sorted(c, reverse=True):
        print(f'   phase {ph}: {c[ph]}')

if __name__ == '__main__':
    main()