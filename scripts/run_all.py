"""Full 00->09 pipeline in one unattended run (universal, resume-aware).

The target PDB ID is propagated to every stage (00/03/04/05/06/07), so the
same orchestrator works for any target. Heavy stages are auto-skipped when
their outputs look complete (unless --force). Stops on first error;
--keep-going continues onward.

Robustness features (added):
- PREFLIGHT check: before stage 03, verifies receptor + readable reference
  ligand exist. Saves ~4 GPU-hours if target preparation is broken.
- redock_done(): safe check of redock_report.txt (tolerates missing/corrupt files).

Stage order: 00..07, then 09 (figures), then 08 (supplementary bundle + ZIP);
08 runs LAST so the ZIP it produces already contains the figures.

Every run is mirrored to logs/run_all_<target>.log (the logs/ folder is
created automatically); this file later becomes S9_run_all_log.txt in the
supplementary bundle (08_supplementary.py). Before each run the previous log
is auto-archived as logs/run_all_<target>_<YYYYmmdd_HHMMSS>.log, so a short
resume run can never destroy a full --force log; 08 picks the largest (most
complete) log for S9.

Usage:
 python scripts/run_all.py                      # auto-skip ready (7gqu)
 python scripts/run_all.py --target 4cgz        # any target
 python scripts/run_all.py --force              # recompute EVERYTHING
 python scripts/run_all.py --keep-going         # do not stop on errors
 python scripts/run_all.py --start 4            # start from stage N
 python scripts/run_all.py --min-phase 2        # library depth
 python scripts/run_all.py --gpus 0,1,2,3 --workers 24   # hardware scaling
"""
import argparse, os, subprocess, sys, time
from datetime import datetime
from pathlib import Path
from rdkit import Chem

BASE = Path(__file__).resolve().parents[1]
COMPLETION = 0.97   # resume markers tolerate up to ~3% permanent failures


class Logger:
    """Mirrors every byte to the terminal AND to the log file (live).

    On start, the previous non-empty log is archived with a timestamp suffix
    (rotation), so history is never lost."""
    def __init__(self, path):
        self.term = sys.stdout
        if path.exists() and path.stat().st_size > 0:
            ts = datetime.fromtimestamp(path.stat().st_mtime).strftime('%Y%m%d_%H%M%S')
            arch = path.with_name(f'{path.stem}_{ts}{path.suffix}')
            n = 1
            while arch.exists():   # never overwrite an existing archive
                arch = path.with_name(f'{path.stem}_{ts}_{n}{path.suffix}')
                n += 1
            path.rename(arch)
        self.f = open(path, 'wb')
    def say(self, msg=''):
        b = (msg + '\n').encode('utf-8')
        self.term.buffer.write(b); self.term.buffer.flush()
        self.f.write(b); self.f.flush()
    def pipe(self, cmd):
        """Run a stage, streaming stdout+stderr live to terminal and log."""
        p = subprocess.Popen(cmd, cwd=BASE, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        fd = p.stdout.fileno()
        while True:
            chunk = os.read(fd, 4096)
            if not chunk:
                break
            self.term.buffer.write(chunk); self.term.buffer.flush()
            self.f.write(chunk); self.f.flush()
        return p.wait()
    def close(self):
        self.f.close()


def n_rows(p):
    try:
        with open(p) as f: return sum(1 for _ in f) - 1
    except Exception:
        return 0


def n_pdbqt(): return len(list((BASE / 'pdbqt').glob('*.pdbqt')))
def n_lib(): return n_rows(BASE / 'data' / 'library' / 'library.tsv')


def preflight(target, log):
    """Check critical inputs BEFORE expensive GPU stages (03 screen_gnina).
    
    This is the key robustness fix: if 00_prepare_target produced a broken
    target (e.g. no ref_ligand.sdf, or unreadable SDF), the pipeline fails
    here — not after 4 hours of GPU compute on stage 03.
    
    Called only when stage 00 will be SKIPPED (receptor already exists).
    If stage 00 will RUN, it creates the files itself, so preflight would
    incorrectly fail.
    """
    d = BASE / 'data' / target
    problems = []
    
    # 1. Receptor files (produced by 00_prepare_target)
    for f in ('receptor.pdbqt', 'receptor_vina.pdbqt', 'receptor.json'):
        if not (d / f).exists():
            problems.append(f"missing {f}")
    
    # 2. Reference ligand: must be READABLE by RDKit
    #    (6vei bug: file existed but was empty; 7kcc bug: file existed but
    #    had non-aromatic rings that crashed RMSD in stage 05)
    ref_ok = False
    last_err = None
    for f in ('ref_ligand.sdf', 'ref_ligand_prot.sdf'):
        p = d / f
        if not p.exists():
            last_err = f"{f} not found"
            continue
        try:
            found_any = False
            for m in Chem.SDMolSupplier(str(p), removeHs=True):
                if m is not None:
                    found_any = True
                    break
            if found_any:
                ref_ok = True
                break
            last_err = f"{f} exists but contains no valid molecules"
        except Exception as e:
            last_err = f"{f} unreadable: {e}"
            continue
    
    if not ref_ok:
        problems.append(f"no readable reference ligand (last error: {last_err})")
    
    if problems:
        msg = (
            f"\n{'='*60}\n"
            f"❌ PREFLIGHT FAILED [{target}]\n"
            f"{'='*60}\n"
            + "\n".join(f"   • {p}" for p in problems) +
            f"\n\n"
            f"   Pipeline остановлен ДО траты GPU-часов на стадии 03.\n"
            f"   Почини подготовку:\n"
            f"     python scripts/00_prepare_target.py --pdb-id {target}\n"
            f"   или вручную:\n"
            f"     cd data/{target} && cp ref_ligand_prot.sdf ref_ligand.sdf\n"
            f"   Затем запусти run_all.py снова.\n"
        )
        log.say(msg)
        raise SystemExit(1)
    
    log.say(f"✅ Preflight OK for {target}: receptor + reference ligand ready")


def redock_done(t):
    """Safe check for stage 05 completion.
    
    Returns True only if redock_report.txt exists AND contains verdict=PASS.
    Tolerates missing/corrupt files (returns False, triggering a rerun).
    """
    rpt = BASE / 'validation' / 'redock' / t / 'redock_report.txt'
    if not rpt.exists():
        return False
    try:
        return 'verdict=PASS' in rpt.read_text()
    except Exception:
        return False


def stages(t, a):
    res = BASE / 'results' / t
    refine = res / 'refine'
    supp = BASE / 'supplementary' / t
    db = BASE / 'data' / 'chembl_37.db'

    cmd01 = [sys.executable, 'scripts/01_build_library.py', '--min-phase', str(a.min_phase)]
    if db.exists(): cmd01 += ['--db', str(db)]

    cmd03 = [sys.executable, 'scripts/03_screen_gnina.py', t]
    if a.gpus: cmd03 += ['--gpus', a.gpus]
    if a.cpu: cmd03 += ['--cpu', str(a.cpu)]

    cmd04 = [sys.executable, 'scripts/04_refine_consensus.py', t]
    if a.gpus: cmd04 += ['--gpus', a.gpus]
    if a.workers: cmd04 += ['--workers', str(a.workers)]

    return [
        ('00 prepare_target',
         [sys.executable, 'scripts/00_prepare_target.py', '--pdb-id', t],
         (BASE / 'data' / t / 'receptor.pdbqt').exists()),
        ('01 build_library', cmd01,
         (BASE / 'data' / 'library' / 'library.tsv').exists()),
        ('02 prepare_ligands',
         [sys.executable, 'scripts/02_prepare_ligands.py'],
         n_lib() > 0 and n_pdbqt() >= COMPLETION * n_lib()),
        ('03 screen_gnina', cmd03,
         n_pdbqt() > 0 and n_rows(res / 'all_docking_results.csv') >= COMPLETION * n_pdbqt()),
        ('04 refine_consensus', cmd04,
         (refine / 'consensus.csv').exists()),
        ('05 redock_validation',
         [sys.executable, 'scripts/05_redock_validation.py', '--target', t],
         redock_done(t)),
        ('06 admet',
         [sys.executable, 'scripts/06_admet.py', '--all',
          '--input', str(refine / 'consensus.csv')],
         (refine / 'admet_all.csv').exists()),
        ('07 final_candidates',
         [sys.executable, 'scripts/07_final_candidates.py', '--target', t],
         (refine / 'final_candidates.csv').exists()),
        # 09 before 08: figures must exist before the bundle/ZIP is packed
        ('09 make_figures',
         [sys.executable, 'scripts/09_make_figures.py', '--target', t],
         (supp / 'figures' / 'FIGURES_DONE.txt').exists()),
        ('08 supplementary',
         [sys.executable, 'scripts/08_supplementary.py', '--target', t, '--zip'],
         (supp / 'S0_manifest.txt').exists()
         and (supp / 'figures' / 'FIGURES_DONE.txt').exists()),
    ]


def main():
    ap = argparse.ArgumentParser(description='Unattended 00->09 pipeline')
    ap.add_argument('--target', default='7gqu', help='PDB ID of the target')
    ap.add_argument('--force', action='store_true', help='recompute everything, ignoring cache')
    ap.add_argument('--keep-going', action='store_true', help='do not stop on errors')
    ap.add_argument('--start', type=int, default=0, help='start from stage N')
    ap.add_argument('--min-phase', type=int, default=3, choices=[0, 1, 2, 3, 4],
                    help='library depth: 3=approved+PhIII (default), 2=+PhII, 0=all')
    ap.add_argument('--gpus', default=None, help='GPU list for 03/04 (e.g. 0,1)')
    ap.add_argument('--cpu', type=int, default=None, help='CPU threads per GNINA worker (03)')
    ap.add_argument('--workers', type=int, default=None, help='CPU pool size for 04 (Vina/LeDock)')
    a = ap.parse_args()

    t = a.target.lower()

    # logs/ is created automatically; the run is mirrored to the log file
    logdir = BASE / 'logs'
    logdir.mkdir(parents=True, exist_ok=True)
    log = Logger(logdir / f'run_all_{t}.log')

    t0 = time.time(); report = []
    log.say(f'=== run_all started {datetime.now().isoformat(timespec="seconds")} ===')
    log.say(f'Parameters: target={t}, min_phase={a.min_phase}, force={a.force}, start={a.start}')

    # ---- PREFLIGHT: check target integrity BEFORE spending GPU-hours ----
    # Only runs when stage 00 will be SKIPPED (receptor already exists).
    # If stage 00 will RUN (first launch / --force), skip preflight because
    # stage 00 itself will create the files.
    receptor_exists = (BASE / 'data' / t / 'receptor.pdbqt').exists()
    if receptor_exists and not a.force:
        preflight(t, log)
    # --------------------------------------------------------------------

    try:
        for name, cmd, done in stages(t, a):
            num = int(name[:2])
            if num < a.start: continue
            if done and not a.force:
                log.say(f'\n=== [{name}] already ready — skipping ===')
                report.append((name, 'SKIP', 0)); continue
            ts = time.time()
            log.say(f'\n{"=" * 20} [{name}] start {"=" * 20}')
            rc = log.pipe(cmd)
            dt = time.time() - ts
            status = 'OK' if rc == 0 else 'FAIL'
            report.append((name, status, dt))
            if rc != 0 and not a.keep_going:
                log.say(f'\n❌ [{name}] exited with code {rc} — stopping.')
                break
    finally:
        log.say(f'\n{"=" * 50}\nSUMMARY (total {time.time() - t0:.0f} s)')
        for name, st, dt in report:
            log.say(f'  {name:<20} {st:<6} {dt:7.0f} s')
        log.say(f'=== run_all finished {datetime.now().isoformat(timespec="seconds")} ===')
        log.close()


if __name__ == '__main__':
    main()