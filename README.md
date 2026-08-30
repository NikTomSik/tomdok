# Universal Drug Repurposing Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22099251.svg)](https://doi.org/10.5281/zenodo.22099251)

A reproducible virtual screening pipeline for any protein target. Runs three
independent docking engines (GNINA, AutoDock Vina, LeDock), applies rule-based
and data-driven ADMET profiling, produces publication figures for all final
candidates, and ships a safety-aware prioritised hit list with a complete
Zenodo/journal supplementary bundle. One command runs every stage; every
stage supports resume; every output column is documented below.

The redocking validation runs as stage 01, immediately after target
preparation, as a PROTOCOL GATE: a broken target or a failed self-dock stops
the pipeline before any GPU-hours are spent on the screen.

Stage numbering matches execution order: 00 prepare, 01 redock gate,
02 library, 03 ligands, 04 screen, 05 consensus, 06 ADMET, 07 prioritisation,
08 figures, 09 supplementary bundle (runs LAST so the ZIP contains figures).

---

## Install

### Automated (recommended)

    bash scripts/setup_env.sh            # installs Miniconda + `screen` env + all deps
    bash scripts/setup_env.sh --check    # verify only (no installs)
    conda activate screen

`setup_env.sh` is idempotent. It handles known quirks: `meeko` and
`dimorphite-dl` are installed via pip, `obabel` via `openbabel`, `vina` via
`autodock-vina`, `gnina` via bioconda; the broken conda-bundled `curl` is
worked around by falling back to the system `curl` or Python `urllib`.

### Manual alternative

    conda create -n screen python=3.11 && conda activate screen
    conda install -c conda-forge rdkit numpy pandas openbabel pdbfixer openmm gnina \
                                 pymol-open-source matplotlib pillow
    conda install -c bioconda autodock-vina
    python -m pip install meeko dimorphite-dl plip
    # LeDock: download ledock/lepro from lephar.com and place in PATH (optional)

---

## Quick start

    conda activate screen
    python scripts/run_all.py --target <your_pdb_id>          # all stages 00 to 09

`run_all.py` flags:

    --force                recompute everything
    --keep-going           do not stop on errors
    --start N              start from stage N
    --min-phase {0..4}     library depth (3 = approved + Phase III, default)
    --target <pdb_id>      PDB identifier of the protein target
    --gpus 0,1             GPU list for GNINA (04/05)
    --cpu N                threads per GNINA worker (04)
    --workers N            CPU pool for Vina/LeDock (05)

Each run is mirrored to `logs/run_all_<target>.log`; the previous log is
auto-archived with a timestamp suffix, so a short resume run never
overwrites a full `--force` log. The bundle (stage 09) picks the largest
(most complete) log for `S9_run_all_log.txt`.

Live monitor (optional, run in a second terminal while stage 04 is going):

    python scripts/monitor_gnina.py --target <your_pdb_id>
    python scripts/monitor_gnina.py --target <id> --interval 15 --top 15 --cnn 0.7 --once

    --interval S           refresh period, seconds (default 30)
    --cnn X                CNN threshold for the table (strictly >, default 0.7)
    --top N                rows in the table (default 15)
    --once                 render once and exit

The monitor prints ONE table: top-N by affinity among molecules with
CNN > X, with drug name and clinical phase (4 = approved, 3 = PhIII).

Manual step-by-step:

    python scripts/00_prepare_target.py --pdb-id <id>
    python scripts/01_redock_validation.py --target <id>        # protocol gate (--basis best|top)
    python scripts/02_build_library.py --min-phase 3
    python scripts/03_prepare_ligands.py
    python scripts/04_screen_gnina.py <id> --gpus 0,1
    python scripts/05_refine_consensus.py <id>
    python scripts/06_admet.py --all
    python scripts/07_final_candidates.py --target <id>
    python scripts/08_make_figures.py --target <id>             # 3D + 2D figures
    python scripts/09_supplementary.py --target <id> --zip      # bundle (includes figures + poses)

---

## Pipeline stages

    00_prepare_target.py    protein cleaning, two receptors (GNINA/AD4/LeDock
                            and Vina), ref-ligand preparation, docking box.
                            Robustness: two-tier ligand selection (drug-like
                            residues preferred; cofactors only as fallback),
                            mmCIF fallback when RCSB has no legacy .pdb,
                            obabel-first bond perception, --reextract flag.
    01_redock_validation.py PROTOCOL GATE: self-docking RMSD right after prep
                            (best-pose or top-pose criterion). Robust ref
                            loading (ref_ligand.sdf -> ref_ligand_prot.sdf),
                            sanitisation + MCS fallback for RMSD. FAIL exits
                            non-zero and stops the pipeline.
    02_build_library.py     ChEMBL library with tunable phase depth + data-driven
                            hERG liability from ChEMBL KCNH2/CHEMBL240 activities
    03_prepare_ligands.py   PDBQT library (pH 7.4, ETKDGv3 + MMFF94, Meeko)
    04_screen_gnina.py      primary GNINA screen (multi-GPU, exh=4)
    05_refine_consensus.py  three-engine consensus (GNINA exh=32, Vina exh=32,
                            LeDock) with filter aff <= -9 and cnn >= 0.7
    06_admet.py             rule-based ADMET + ChEMBL-derived hERG/withdrawn
    07_final_candidates.py  final CLEAN / FLAGGED prioritisation
    08_make_figures.py      publication figures for ALL final candidates
                            (PyMOL 3D per ligand + overlays + redock;
                             RDKit 2D structures + PLIP interaction summary;
                             .pse sessions for overlays)
    09_supplementary.py     build Zenodo/journal supplementary bundle
                            (S0-S9 files, figures, poses, manifest, ZIP)

---

## Repository layout

    scripts/
      setup_env.sh            install / verify environment
      run_all.py              orchestrator 00 to 09 (resume-aware, log rotation,
                              preflight check before GPU stages)
      monitor_gnina.py        live top-N monitor for stage 04 (second terminal)
      00_prepare_target.py    protein + receptors + box + ref-ligand
      01_redock_validation.py RMSD protocol gate (--basis best|top)
      02_build_library.py     ChEMBL library (--min-phase) + hERG liability
      03_prepare_ligands.py   PDBQT library
      04_screen_gnina.py      primary GNINA screen (multi-GPU)
      05_refine_consensus.py  3-engine consensus
      06_admet.py             ADMET profiling
      07_final_candidates.py  final CLEAN / FLAGGED prioritisation
      08_make_figures.py      3D + 2D figures (PyMOL + PLIP)
      09_supplementary.py     supplementary bundle builder (Zenodo/journal)
    data/<target>/            receptor.pdbqt, receptor_vina.pdbqt, receptor.json,
                              ref_ligand.sdf, ref_ligand_prot.sdf, ref_ligand.pdbqt
    data/library/library.tsv  curated library (SMILES + properties + phase + herg)
    pdbqt/                    prepared ligands
    logs/                     run_all execution logs (auto-created, rotated)
    supplementary/<target>/   generated supplementary bundle
      S0_manifest.txt         manifest: target, funnel, parameters, RMSD, figures
      S1_library.tsv          full repurposing library
      S2_primary_screen_full.csv primary GNINA screen
      S3_gnina_refined.csv    refined GNINA top hits
      S4_consensus.csv        3-engine consensus + phase/approval + safety class
      S5_admet.csv            full ADMET profiles
      S6_final_candidates.csv final CLEAN / FLAGGED prioritisation
      S7_redock_validation.txt redocking RMSD report
      S8_parameters.json      machine-readable parameters
      S9_run_all_log.txt      most complete run_all execution log
      figures/                publication figures (stage 08)
        overlay/Fig_overlay_{CLEAN,FLAGGED}.{png,pse}
        redock/Fig_redock.png
        per_ligand/{CLEAN,FLAGGED}/<name>.png
        interactions_2d/{CLEAN,FLAGGED}/<name>.png
      poses/                  universal 3D data for any viewer
        receptor.pdb          cleaned protein
        <chembl_id>.sdf       per-ligand docked pose
    results/<target>/
      all_docking_results.csv primary screen
      refine/consensus.csv    consensus
      refine/admet_all.csv    ADMET
      refine/final_candidates.csv final prioritisation
    validation/redock/<target>/ redock_report.txt (per-target)
    environment.yml  LICENSE  DATA_LICENSE  CITATION.cff  .gitignore

---

## Output columns

### 02 build_library -> library.tsv

| Column|Meaning|
| ---|---|
| chembl_id, pref_name|ChEMBL identifier and name|
| max_phase|4 = approved, 3 = Phase III, 2 = Phase II|
| first_approval|year of first approval|
| canonical_smiles|input SMILES|
| smiles_clean|parent SMILES (desalted, largest fragment)|
| mw, alogp, hba, hbd, psa, rtb|physico-chemistry from ChEMBL|
| ro5_viol, qed|Lipinski violations / drug-likeness|
| herg|hERG liability (RISK if IC50/Ki < 10 uM vs KCNH2 in ChEMBL; ok / '' if no data)|

### 04 screen_gnina -> all_docking_results.csv

| Column|Meaning|
| ---|---|
| affinity|GNINA binding affinity, kcal/mol (more negative = stronger binder)|
| cnn_score|CNN pose confidence, 0 to 1 (threshold >= 0.7)|

### 05 refine_consensus -> consensus.csv

| Column|Meaning|
| ---|---|
| gnina32_aff, cnn|refined GNINA (exhaustiveness 32)|
| vina_aff|AutoDock Vina (exhaustiveness 32)|
| ledock_score|LeDock (systematically softer scale than GNINA/Vina)|
| n_engines_le9|number of engines (out of 3) with score <= -9|

### 01 redock_validation

Two RMSD metrics against the crystallised inhibitor:

- **top-pose RMSD** — RMSD of the best-scored pose (tests *scoring* power)
- **best-pose RMSD** — lowest RMSD among all poses (tests *sampling* power)

Verdict is PASS when `best-pose RMSD < 2.0 A` (community standard for
protocol validation). Use `--basis top` for the stricter scoring criterion.
The reference is loaded robustly (`ref_ligand.sdf`, falling back to
`ref_ligand_prot.sdf`); pose and reference are sanitised with unified
aromaticity/stereo perception, with an MCS-based atom map as final RMSD
fallback. On verdict=FAIL the script exits non-zero so `run_all.py` stops
before the expensive stages. Report is written to
`validation/redock/<target>/redock_report.txt` and contains both RMSD values
plus the rank of the native-like pose.

### 06 admet -> admet_all.csv

| Column|Meaning|
| ---|---|
| MW, LogP, TPSA, RB|mass / lipophilicity / polarity / flexibility|
| HBA, HBD, lipinski|H-bond acceptors/donors, Ro5 violations|
| qed|drug-likeness 0 to 1 ( > 0.5 is good)|
| hia|human intestinal absorption (high / low)|
| bbb|blood-brain barrier penetration (yes / no)|
| herg|hERG-block risk from ChEMBL KCNH2/CHEMBL240 activities (RISK / ok)|
| pains|pan-assay interference (ALERT / ok)|
| withdrawn|withdrawn per ChEMBL withdrawn_flag (YES / no)|

### 07 final_candidates -> final_candidates.csv

| Column|Meaning|
| ---|---|
| phase|max clinical phase (4 = approved, 3 = PhIII, 2 = PhII)|
| approved|YES when phase = 4|
| gnina32, vina, ledock, n_le9|consensus scores|
| qed, bbb|drug-likeness / BBB penetration|
| cls|CLEAN (no safety flags) / FLAGGED|
| flags|semicolon-separated list: hERG; withdrawn; PAINS; Lip >1|

### 08 make_figures -> supplementary/<target>/figures/

| File|Meaning|
| ---|---|
| overlay/Fig_overlay_CLEAN.{png,pse}|all CLEAN candidates overlaid; .pse is an editable PyMOL session|
| overlay/Fig_overlay_FLAGGED.{png,pse}|all FLAGGED candidates overlaid|
| redock/Fig_redock.png|crystal ligand (green) vs top redock pose (magenta)|
| per_ligand/{CLEAN,FLAGGED}/<name>.png|per-ligand 3D: site sticks + polar contacts|
| interactions_2d/{CLEAN,FLAGGED}/<name>.png|2D structure (RDKit) + scores + PLIP interactions|

CLEAN figures are rendered with green ligand sticks, FLAGGED with orange;
the redock figure uses green for the crystallographic pose and magenta for
the top-ranked self-docked pose.

### 09 supplementary -> supplementary/<target>/poses/

| File|Meaning|
| ---|---|
| receptor.pdb|cleaned protein structure (stage 00)|
| <chembl_id>.sdf|docked pose of each final candidate|

The `poses/` folder contains universal 3D coordinates that open in any
viewer (PyMOL, ChimeraX, VMD, Mol*, iCn3D) — independent of PyMOL sessions.

---

## Interpretation

**GNINA and Vina agreement.** The primary robustness signal is that GNINA and
Vina affinities agree to Delta-Delta-G < 0.1 kcal/mol for the top hits. Both
engines use physics-based scoring but differ in implementation, so agreement
across them is strong evidence the signal is real.

**LeDock as orthogonal check.** LeDock scores are on a systematically softer
scale than GNINA/Vina and are therefore used as a rank-level orthogonal
check rather than a threshold voter: for the consensus hits, LeDock scores
consistently fall in the upper part of their score range, independently
supporting the ranking. The binary consensus (`n_engines_le9 <= -9`) is
driven by GNINA and Vina.

**CLEAN vs FLAGGED.** A CLEAN hit has no safety flags and is ready for wet-lab
follow-up. A FLAGGED hit is a strong binder but carries one or more of
hERG / withdrawn / PAINS / Lip>1; this is not a disqualification but must be
acknowledged in any downstream discussion.

**Approved vs Phase III.** Approved compounds (phase = 4) are easier to
repurpose because their pharmacokinetics and safety profile are already
established in humans. Phase-III compounds are promising but require
additional justification.

**Known biases.** ADMET is rule-based for phys-chem / HIA / BBB / PAINS, and
data-driven for hERG (ChEMBL activities vs KCNH2) and withdrawn (ChEMBL
flags). Cationic amphiphilic drugs (antihistamines, antipsychotics) are
enriched in the hits and may be frequent-hitters; consider counter-screens
for non-specific binding.

---

## Example: 7GQU (WRN helicase)

| Funnel step|Count|
| ---|---|
| Library (ChEMBL 37, phase >= 3, deduplicated)|3,369|
| Prepared PDBQT (stage 03, QC attrition −43)|3,326|
| Docked successfully (stage 04)|3,326|
| Primary hits (aff <= −9.0, CNN >= 0.7)|71|
| Consensus hits (GNINA32 + Vina + LeDock)|64|
| CLEAN / FLAGGED|44 / 20|

Self-docking RMSD: top-pose 2.75 A, best-pose 1.28 A → **PASS** (native-like
pose sampled at rank 2; scoring places it second, motivating 3-engine
consensus).

## Example: 6VEI (mutant IDH1 R132H)

| Funnel step|Count|
| ---|---|
| Library (ChEMBL 37, phase >= 3, deduplicated)|3,369|
| Prepared PDBQT (stage 03, QC attrition −43)|3,326|
| Docked successfully (stage 04)|3,135|
| Primary hits (aff <= −9.0, CNN >= 0.7)|33|
| Consensus hits (GNINA32 + Vina + LeDock)|31|
| CLEAN / FLAGGED|20 / 11|

Self-docking RMSD: top-pose 1.14 A, best-pose 0.45 A → **PASS**. The blind
screen recovered two clinically validated mutant-IDH inhibitors: vorasidenib
(the co-crystallised ligand; triple-engine consensus, CNN 0.998) and
enasidenib (approved mutant-IDH2 inhibitor; top affinity −11.84).

---

## Hardware scaling

The pipeline is fully parameterised. `--gpus` lists GPUs (one GNINA worker
per card), `--cpu` is threads per GNINA instance, `--workers` is the CPU pool
for Vina and LeDock in the consensus stage.

| Configuration|Screening (04)|Consensus (05)|
| ---|---|---|
| 4x GPU + 16-core server|--gpus 0,1,2,3 --cpu 8|--workers 24|
| 2x GPU + 8 cores (reference, e.g. 2x RTX 3060)|--gpus 0,1 --cpu 4|--workers 8|
| 1x GPU + 6 cores (workstation)|--gpus 0 --cpu 4|--workers 5|
| 1x GPU with >= 6 GB VRAM (boost mode, two workers share one card)|--gpus 0,0 --cpu 3|--workers 5|

Requirements: about 2 GB VRAM per GNINA instance. Ligand preparation (03)
uses `min(8, n_cores)` processes automatically. Every stage prints a live ETA.

---

## Reproducibility

Random seed is 42 in all engines. Exhaustiveness is 32 for refinement and 4
for primary screening. Primary thresholds are affinity <= -9 kcal/mol and
CNN score >= 0.7. Docking box is 25 A per axis (clip of ptp + 5 A, with a
floor of 25 A and a ceiling of 27 A). Library is built from ChEMBL 37.
GNINA, Vina and LeDock are installed from conda/bioconda. Full recomputation
of all stages:

    python scripts/run_all.py --force --target <id>

---

## Notes

Stage numbering. Numbers match execution order (00, 01, ..., 09). The redock
gate is stage 01 so that a broken target fails in ~1 minute, not after
~4 GPU-hours of screening.

Target preparation robustness. `00_prepare_target.py` uses two-tier ligand
selection (drug-like residues preferred; cofactors/coenzymes such as NDP,
NAD, ATP, FAD, HEM, SAM used only as fallback), obabel-first bond perception
with RDKit fallback, mmCIF fallback when RCSB has no legacy .pdb (HTTP 404),
and a Meeko -> OpenBabel fallback for the reference PDBQT when Gasteiger
charges fail. `--reextract` forces re-extraction of the reference ligand.

Redock robustness. `01_redock_validation.py` loads `ref_ligand.sdf` with
fallback to `ref_ligand_prot.sdf`, sanitises pose and reference with unified
aromaticity/stereo perception, and falls back to an MCS-based atom map for
RMSD when graph matching fails.

AutoDock 4 was attempted as a fourth consensus engine. The conda-provided
binary (4.2.7.x) is incompatible with every publicly available
`AD4_parameters.dat` file (autogrid4 exits with "Van der Waals coefficient cA
is not a number"). For this reason AutoDock 4 is excluded from the consensus.
The remaining three engines give a robust ranking.

Cmax estimation was initially part of the pipeline but was removed. A
simple one-compartment PBPK-lite model does not add meaningful signal for
repurposing of compounds with unknown standard doses, and the question is
better evaluated case-by-case on the final shortlist.

---

## AI disclosure

AI assistants were used for code development and documentation only; all
scientific decisions were made by the author.

---

## License

- **Code** (`scripts/`, `setup_env.sh`): MIT License — see `LICENSE`.
- **Data and results** (docking scores, consensus tables, ADMET profiles,
redock reports, publication figures, 3D poses): Creative Commons
Attribution 4.0 International (CC BY 4.0) — see `DATA_LICENSE`.

CC BY means anyone may reuse the results, including commercially, provided
they give appropriate credit (cite the publication / repository). The MIT
license on the code requires keeping the copyright notice when copying the
software.