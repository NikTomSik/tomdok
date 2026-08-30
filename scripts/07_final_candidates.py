"""Final candidate prioritisation: 3-engine consensus + ADMET flags (no Cmax).
Stage 07 of the renumbered pipeline (00 prep, 01 redock gate, 02 library,
03 ligands, 04 screen, 05 consensus, 06 ADMET, 07 candidates, 08 figures,
09 supplementary).

Adds phase (max clinical phase), approved (YES/no) and dledock:
  dledock = ledock_score - ledock_native,
where ledock_native is the per-target LeDock baseline written by
01_redock_validation.py ('ledock_native=<x>' in
validation/redock/<t>/redock_report.txt, or ledock_baseline.json).
The value is propagated from consensus.csv when 05_refine_consensus.py
already computed it, otherwise recomputed here. Negative dledock = the
candidate scores better than the native co-crystallised ligand in LeDock.

CLEAN   = no safety flags
FLAGGED = hERG / withdrawn / PAINS / Lipinski>1

Universal: works for any target via --target.
Writes final_candidates.csv, prints both lists and glossary.

Usage:
  python scripts/07_final_candidates.py --target 7gqu
  python scripts/07_final_candidates.py --target 9fza --top 20
"""
import argparse, csv, json, os, sys
from pathlib import Path

BASE = Path(os.environ.get('SCREENING_ROOT', str(Path(__file__).resolve().parents[1])))

LEGEND = '''
=== Column glossary (final) ===
phase   - max clinical phase (4=approved, 3=PhIII, 2=PhII)
approved- regulatory approval (YES/no)
gnina32 - affinity GNINA exh=32, kcal/mol (more negative = stronger)
vina    - affinity Vina exh=32
ledock  - score LeDock (systematically softer scale than GNINA/Vina)
dledock - ledock - ledock_native (per-target baseline from 01 redock);
          negative = better than the native co-crystallised ligand
n_le9   - number of engines (out of 3) with score <= -9 (consensus strength)
qed     - drug-likeness 0-1 (>0.5 good)
bbb     - BBB penetration (yes/no)
cls     - CLEAN (no safety flags) / FLAGGED (safety flags)
flags   - hERG / withdrawn / PAINS / Lip>1
'''


def read(p):
    with open(p, newline='') as f:
        return {r['chembl_id']: r for r in csv.DictReader(f)}


def fs(x, w=8, nd=2):
    """Safe float formatting: '     n/a' when missing/empty."""
    try:
        return f"{float(x):{w}.{nd}f}"
    except (TypeError, ValueError):
        return ' ' * (w - 3) + 'n/a'


def load_phases():
    """Read max_phase from library.tsv (all SMILES that entered the library)."""
    ph = {}
    lib = BASE / 'data' / 'library' / 'library.tsv'
    if lib.exists():
        with open(lib, newline='', encoding='utf-8') as f:
            for r in csv.DictReader(f, delimiter='\t'):
                cid = (r.get('chembl_id') or '').strip()
                try:
                    ph[cid] = int(r.get('max_phase') or 0)
                except (TypeError, ValueError):
                    ph[cid] = 0
    return ph


def load_ledock_native(t):
    """Per-target LeDock baseline written by 01_redock_validation.py.

    Looks for 'ledock_native=<value>' in
    validation/redock/<t>/redock_report.txt, then in ledock_baseline.json.
    Returns float or None.
    """
    vdir = BASE / 'validation' / 'redock' / t
    rpt = vdir / 'redock_report.txt'
    if rpt.exists():
        try:
            for line in rpt.read_text(errors='ignore').splitlines():
                if line.strip().startswith('ledock_native'):
                    return float(line.split('=')[-1].strip())
        except Exception:
            pass
    js = vdir / 'ledock_baseline.json'
    if js.exists():
        try:
            return float(json.loads(js.read_text()).get('ledock_native'))
        except Exception:
            pass
    return None


def main():
    ap = argparse.ArgumentParser(description='Final CLEAN/FLAGGED prioritisation')
    ap.add_argument('--target', default='7gqu', help='PDB ID of the target')
    ap.add_argument('--top', type=int, default=15, help='rows to print per list')
    a = ap.parse_args()
    t = a.target.lower()

    ref = BASE / 'results' / t / 'refine'
    for req in (ref / 'consensus.csv', ref / 'admet_all.csv'):
        if not req.exists():
            sys.exit(f'ERROR: missing {req} — run 05_refine_consensus.py and '
                     f'06_admet.py first')

    cons = read(ref / 'consensus.csv')
    admet = read(ref / 'admet_all.csv')
    ph = load_phases()
    native_ld = load_ledock_native(t)

    missing_admet = [cid for cid in cons if cid not in admet]
    if missing_admet:
        print(f'[Warning] {len(missing_admet)} consensus IDs absent from '
              f'admet_all.csv (kept with empty ADMET, treated as CLEAN)')
    if native_ld is None:
        print('[Warning] ledock_native baseline not found '
              '(run 01_redock_validation.py); dledock left empty')

    rows = []
    for cid, c in cons.items():
        adm = admet.get(cid, {})
        flags = []
        if adm.get('withdrawn') == 'YES': flags.append('withdrawn')
        if adm.get('herg') == 'RISK':     flags.append('hERG')
        if adm.get('pains') == 'ALERT':   flags.append('PAINS')
        try:
            if int(adm.get('lipinski', 0)) > 1: flags.append('Lip>1')
        except (TypeError, ValueError):
            pass

        # dledock: propagate from consensus when present, else recompute here
        dld = c.get('dledock', '')
        if dld in (None, '') and native_ld is not None:
            try:
                dld = f"{float(c['ledock_score']) - native_ld:.2f}"
            except (TypeError, ValueError):
                dld = ''

        phase = ph.get(cid, -1)
        rows.append(dict(chembl_id=cid, name=c.get('name', ''), phase=phase,
                         approved='YES' if phase == 4 else 'no',
                         gnina32=c['gnina32_aff'], vina=c['vina_aff'],
                         ledock=c['ledock_score'], dledock=dld,
                         n_le9=c['n_engines_le9'], qed=adm.get('qed', ''),
                         bbb=adm.get('bbb', ''),
                         cls='FLAGGED' if flags else 'CLEAN',
                         flags=';'.join(flags)))

    rows.sort(key=lambda r: float(r['gnina32']))

    out = ref / 'final_candidates.csv'
    cols = ['chembl_id', 'name', 'phase', 'approved', 'gnina32', 'vina', 'ledock',
            'dledock', 'n_le9', 'qed', 'bbb', 'cls', 'flags']
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    clean = [r for r in rows if r['cls'] == 'CLEAN']
    flagged = [r for r in rows if r['cls'] == 'FLAGGED']

    print(f'=== CLEAN ({len(clean)}) — top by affinity ===')
    print(f'{"CHEMBL":<13}{"name":<16}{"ph":>3}{"app":>4}{"GNINA32":>8}{"VINA":>7}'
          f'{"LeDock":>8}{"dLeDock":>8}{"N<=-9":>6}{"QED":>6}{"BBB":>5}')
    for r in clean[:a.top]:
        print(f"{r['chembl_id']:<13}{r['name'][:15]:<16}{int(r['phase']):>3d}"
              f"{r['approved']:>4}{fs(r['gnina32'])}{fs(r['vina'], 7)}"
              f"{fs(r['ledock'])}{fs(r['dledock'])}{int(r['n_le9']):6d}"
              f"{r['qed']:>6}{r['bbb']:>5}")

    print(f'\n=== FLAGGED ({len(flagged)}) — strong affinity but flagged ===')
    for r in flagged[:a.top]:
        print(f"{r['chembl_id']:<13}{r['name'][:15]:<16}{int(r['phase']):>3d}"
              f"{r['approved']:>4}{fs(r['gnina32'])}  [{r['flags']}]"
              f"{fs(r['dledock'])}")

    print(f'\n✅ Saved: {out}')
    print(LEGEND)


if __name__ == '__main__':
    main()