[![English](https://img.shields.io/badge/lang-English-blue.svg)](free_energy.md)
[![Español](https://img.shields.io/badge/lang-Español-green.svg)](free_energy.es.md)

# Free Energy — detailed documentation

[← Back to the repository overview](../README.md) ·
[← Basic Analyses](basic_analyses.md)

**Binding energetics of protein–ligand complexes from Autogromacs
trajectories: MM-GBSA, direct LJ/Coulomb interaction energy and binding
validation, collected into one Excel file for all your systems.**

The tool has two scripts:

| Script | Where it runs | What it does |
|---|---|---|
| `free_energy.sh` | Inside each system's `gromacs/` folder | Computes the energetics of **one** system and writes `summary.csv` |
| `collect_free_energy.py` | Anywhere, once | Gathers the `summary.csv` of **all** systems into one Excel workbook |

A Spanish quick guide is also available:
[`LEEME_free_energy.txt`](../LEEME_free_energy.txt).

---

## Contents

1. [The idea](#1-the-idea)
2. [Requirements](#2-requirements)
3. [Quick start](#3-quick-start)
4. [Modes and options](#4-modes-and-options)
5. [Automatic detection of receptor and ligand](#5-automatic-detection-of-receptor-and-ligand)
6. [The conda environment for gmx_MMPBSA](#6-the-conda-environment-for-gmx_mmpbsa)
7. [What is computed, step by step](#7-what-is-computed-step-by-step)
8. [Output files](#8-output-files)
9. [Quality flags](#9-quality-flags)
10. [How to read and report the results](#10-how-to-read-and-report-the-results)
11. [collect_free_energy.py: one Excel for all systems](#11-collect_free_energypy-one-excel-for-all-systems)
12. [Running many systems at once](#12-running-many-systems-at-once)
13. [Troubleshooting](#13-troubleshooting)

---

## 1. The idea

Same philosophy as Autogromacs and Basic Analyses: copy **one file** into the
CHARMM-GUI `gromacs/` folder after production and run it:

```bash
bash free_energy.sh
```

There is no `make_ndx` menu and no group numbers to choose. The script reads
`topol.top` and decides by itself:

- **Receptor** = the protein chain(s) (`PROA`, `PROB`, ...).
- **Ignored** = water (`TIP3`, ...) and ions (`POT`, `CLA`, `SOD`, ...).
- **Ligand** = whatever remains (a fatty acid, a PET hexamer, PIP2, ...).

Before computing anything it prints what it detected (receptor, ligand,
heavy atoms, charges, window, temperature). If the system is ambiguous, it
**stops** and tells you which option to add. It never guesses silently.

All calculations use the **same window of frames**, so MM-GBSA, LJ/Coulomb
and validation are directly comparable. The defaults are the recommended
protocol: if every system is run with the defaults, the results can be
compared with each other.

## 2. Requirements

| Component | Needed for | Notes |
|---|---|---|
| GROMACS (`gmx` or `gmx_mpi`) | Everything | System build; found like in Autogromacs (`PATH` or standard `GMXRC`) |
| Python 3 | Everything | Standard library only for `free_energy.sh` |
| conda + an env with **gmx_MMPBSA** | MM-GBSA | Not needed in `lj` mode |
| `mpirun` | Optional | Parallel MM-GBSA; falls back to serial |
| **openpyxl** | `collect_free_energy.py` | Without it, CSV files are written instead of Excel |
| Finished production | Everything | `dynamic<N>ns.{tpr,xtc}` (Autogromacs) or `step5_production.*` |
| `topol.top`, `index.ndx` | Everything | CHARMM-GUI files in `gromacs/` |

No GPU is needed: the energy re-run uses the CPU (energy groups are not
supported on the GPU) and only evaluates a few hundred frames.

## 3. Quick start

```bash
cd /path/to/my_system/gromacs
cp /path/to/free_energy.sh .

bash free_energy.sh --check      # 1. see what was detected (seconds, computes nothing)
bash free_energy.sh              # 2. compute (MM-GBSA + LJ + validation)

# For long jobs, in the background:
nohup bash free_energy.sh > free_energy.log 2>&1 &
```

Then, once all systems are done, from any folder:

```bash
python3 collect_free_energy.py /path/to/Project_A /path/to/Project_B -o Energies.xlsx
```

The first time, run it in a terminal (not with `nohup`) so it can ask which
conda environment has gmx_MMPBSA. The answer is remembered (see
[section 6](#6-the-conda-environment-for-gmx_mmpbsa)).

## 4. Modes and options

### Modes

| Command | What it does |
|---|---|
| `bash free_energy.sh` | MM-GBSA + LJ/Coulomb + validation (default) |
| `bash free_energy.sh --check` | Only shows what was detected and writes the index; computes nothing |
| `bash free_energy.sh mmgbsa` | Only MM-GBSA (+ validation) |
| `bash free_energy.sh lj` | Only LJ/Coulomb interaction (+ validation); fast, no conda needed |
| `bash free_energy.sh summary` | Rebuilds `summary.csv` from results that already exist |
| `bash free_energy.sh --help` / `--version` | Help / version |

### Options (all optional)

| Option | Default | Use it when |
|---|---|---|
| `--start NS` | 25 % of the run | The RMSD shows the system settles later (e.g. `--start 80`) |
| `--end NS` | End of the run | You want to exclude the end |
| `--frames N` | 300 | More frames (lower SEM, slower) or fewer (faster); minimum 10 |
| `--ligand NAME` | Auto | Detection is ambiguous: two ligands, a cofactor, a protein–protein complex. Comma-separated for several |
| `--receptor NAME` | All proteins | The receptor is not the default |
| `--np N` | min(cores, 16) | Number of MPI processes for gmx_MMPBSA |
| `--env NAME` | Asked once, remembered | Use another conda env; `--env ask` to choose again |
| `--no-decomp` | Decomposition on | Faster run without per-residue decomposition |
| `--base NAME` | Auto | Name the trajectory explicitly, e.g. `dynamic200ns` |
| `--name NAME` | Auto | The system name in the tables is too generic |

> **Golden rule:** to compare systems, run **all of them with the same
> options** (same window and number of frames). With the defaults this is
> automatic.

### Situations and what to type

| Situation | Command |
|---|---|
| First time on a system (recommended) | `bash free_energy.sh --check`, then `bash free_energy.sh` |
| Normal case (protein + one ligand) | `bash free_energy.sh` |
| Only the LJ energy, quickly | `bash free_energy.sh lj` |
| Faster, without key residues | `bash free_energy.sh --no-decomp` |
| RMSD reaches a plateau at 80 ns | `bash free_energy.sh --start 80` (same window in all systems) |
| Stops with "more than one candidate ligand" | `--ligand NAME` (or `NAME1,NAME2`), as the message suggests |
| Protein–protein complex | `--ligand PROB` (receptor = `PROA`) |
| Rebuild the table only | `bash free_energy.sh summary` |
| System name comes out generic | `--name Gelsolin_LPA` |

## 5. Automatic detection of receptor and ligand

The script parses `topol.top` and every `.itp` it includes, and classifies
each molecule in the `[ molecules ]` section:

| Class | Rule |
|---|---|
| **Water** | Known water names (`TIP3`, `SOL`, `SPC`, `WAT`, ...) |
| **Ion** | Known ion names (`POT`, `CLA`, `SOD`, `MG`, `ZN`, ...) or single charged atoms |
| **Protein** | `PROA`, `PROB`, ... or more than half of its residues are amino acids |
| **Other** | Everything else: ligand candidates |

Then:

- **Ligand** = the "other" molecules with **≤ 5 copies**. Molecules with more
  than 5 copies (membrane lipids, co-solvents) are ignored and reported in a
  note.
- If a ligand has several copies (≤ 5), all of them go into the `Ligand`
  group, and a note says so.
- **Receptor** = all protein chains.

The script **stops without computing anything** when:

- there are two different ligand candidates (it lists them and suggests
  `--ligand`);
- there is no small molecule but several protein chains (it suggests
  `--ligand <chain>` for a protein–protein complex);
- there is no ligand at all (only protein, water and ions);
- a molecule in `topol.top` has no definition in the `.itp` files.

It also warns (without stopping) if the **net charge** of the system is not
zero (missing ions?) or the ligand is **highly charged** (|q| ≥ 3), since
then the electrostatic terms will be large and noisy.

Heavy atoms are counted by **mass ≥ 4**, so the count is correct even with
hydrogen mass repartitioning. The index `index_FE.ndx` keeps the CHARMM-GUI
groups and adds `Receptor`, `Ligand`, `Receptor_BB` (N, CA, C),
`Ligand_heavy` and `Complex`.

The **system name** used in tables is the first folder above `gromacs/` that
is not a generic name (`gromacs`, `charmm-gui-*`, `complex`, `system`, `md`,
`run`, ...). Change it with `--name`.

## 6. The conda environment for gmx_MMPBSA

`gmx` / `gmx_mpi` is the system GROMACS (used for `trjconv`, the re-run,
etc.). gmx_MMPBSA lives separately, in a conda environment.

- The **first time** you run the script in a terminal, it looks for the
  environments that contain gmx_MMPBSA and asks which one to use:

  ```
   Which conda environment has gmx_MMPBSA?   (conda: /home/user/miniconda3)
     1) gmxMMPBSA            /home/user/miniconda3/envs/gmxMMPBSA
   Number or env name [Enter = gmxMMPBSA]:
  ```

  Enter accepts the suggestion. The answer is saved in
  `~/.config/free_energy/conda_env` and is not asked again (it is also used
  under `nohup`).
- To choose again: `bash free_energy.sh --env ask`, or directly
  `--env OTHER_ENV`.
- It works even when conda is not in the `PATH` (typical under `nohup`): it
  looks in `~/miniconda3`, `~/anaconda3`, `~/miniforge3`, `~/mambaforge`
  and `/opt/conda`.
- Without a terminal, if there is only one env with gmx_MMPBSA (or one
  called `gmxMMPBSA`), it is used automatically; otherwise the script stops
  and asks you to run it once in a terminal or pass `--env`.
- If the env's own `gmx` is too old to read the `.tpr`, gmx_MMPBSA
  automatically uses the system GROMACS instead. This is expected and needs
  no action.

## 7. What is computed, step by step

### [1] Window trajectory

```
gmx trjconv -b <start> -e <end> -dt <step> -pbc mol -center   (center: Receptor · output: System)
```

From the 25 % mark to the end of the run, about 300 frames (e.g. 200 ns →
50–200 ns, one frame every 500 ps). The frame step is rounded to a multiple
of the saved-frame spacing. The result, `window.xtc`, is shared by all the
following calculations.

### [2] Validation (always computed)

| Quantity | Command | Purpose |
|---|---|---|
| Contacts < 4 Å and minimum distance | `gmx mindist -d 0.4` | Is the ligand still touching the receptor? |
| Ligand RMSD | `gmx rms` (fit `Receptor_BB`, calc `Ligand_heavy`) | Did the ligand stay in its site? |
| Hydrogen bonds | `gmx hbond -r Receptor -t Ligand` (legacy fallback) | Specific polar interactions |

These numbers feed the quality flags (`LIGAND_UNBOUND`, `LIGAND_DISPLACED`).

### [3] LJ / Coulomb interaction energy

The production `.mdp` is copied with `energygrps = Receptor Ligand` and
output turned off, a new `.tpr` is built and the window is re-evaluated with
`gmx mdrun -rerun` on the CPU. `gmx energy` then extracts
`LJ-SR:Receptor-Ligand` and `Coul-SR:Receptor-Ligand`.

Reported (kcal/mol): `E_LJ`, `E_Coul`, `E_int = E_LJ + E_Coul` (mean, SD
and GROMACS error estimate), and `E_LJ` per ligand heavy atom.

### [4] MM-GBSA (gmx_MMPBSA)

Settings:

| Parameter | Value |
|---|---|
| Model | Generalized Born, `igb=5` (OBC II) |
| Radii | `mbondi2` (`PBRadii=3`) |
| Salt | 0.15 M |
| Temperature | Read from the `.tpr` |
| Entropy | C2 (second-order cumulant) |
| Decomposition | Per residue (`idecomp=2`), residues within 6 Å |

It runs with `mpirun -np N` when available, retrying in serial if MPI fails.

Reported (kcal/mol, mean / SD / SEM): `dVDW`, `dEEL`, `dEGB`, `dESURF`,
`dGgas`, `dGsolv` and **`dH`** (the MM-GBSA binding enthalpy). From the
fluctuation of `dGgas` the script computes the C2 entropy term
`-TdS = σ²/(2kT)` and `dG = dH − TdS`, plus ligand efficiency
(`LE = dH / heavy atoms`) and the **10 residues that contribute most**.

### [5] Summary

Everything is gathered in one row (`summary.csv`) and printed:

```
 RESULTS (kcal/mol, mean ± SEM; LJ terms mean ± Err.Est.)   frames: 301
   dVDW     -38.12 ± 0.21    dEEL     -10.45 ± 0.40
   dEGB      17.80 ± 0.35    dESURF    -4.95 ± 0.02
   dH       -35.72 ± 0.25    -TdS(C2) 9.80   dG(C2) -25.92
   E_LJ     -40.10 ± 0.60    E_Coul   -12.30 ± 0.90    E_int -52.40 ± 1.08
   H-bonds 1.20   contacts<4A 310   ligand RMSD 0.18 nm
   Top residues: LEU45(-2.10); PHE17(-1.85); ...
   FLAGS: OK
```

*(Example values for illustration.)*

**Resumable.** Each step is skipped if its output already exists. To redo a
step, delete its output and run again: `MMGBSA/FINAL_RESULTS_MMPBSA.dat`,
`LJ/lj_coul.xvg`, or `window.xtc` (which redoes everything that depends on
it).

## 8. Output files

```
Analyses/FreeEnergy_<system>_<N>ns/
├── detected_system.txt     what was detected (for your lab notebook)
├── index_FE.ndx            CHARMM-GUI groups + Receptor, Ligand, Receptor_BB, Ligand_heavy, Complex
├── window.xtc              the frames used by every calculation
├── validation/             hbonds.xvg, numcont.xvg, mindist.xvg, ligand_rmsd.xvg
├── LJ/                     rerun.mdp, rerun.tpr, lj_coul.xvg, energy_stdout.txt, logs
├── MMGBSA/                 mmgbsa.in, FINAL_RESULTS_MMPBSA.dat/.csv,
│                           FINAL_DECOMP_MMPBSA.dat/.csv, gmx_MMPBSA.out
├── summary.csv             ONE row with every result of this system
└── top_residues.csv        the 10 residues that contribute most
```

Main columns of `summary.csv`:

| Group | Columns |
|---|---|
| Identity | `system`, `receptor`, `ligand`, `ligand_resname`, `ligand_heavy_atoms`, `ligand_charge`, `receptor_charge` |
| Window | `sim_ns`, `window_start_ns`, `window_end_ns`, `frame_step_ps`, `temperature_K`, `n_frames` |
| MM-GBSA | `dVDW_*`, `dEEL_*`, `dEGB_*`, `dESURF_*`, `dGgas_*`, `dGsolv_*`, `dH_*` (`_mean`, `_sd`, `_sem`), `sigma_dGgas`, `mTdS_C2`, `dG_C2`, `LE_dH`, `LE_dG`, `top5_residues` |
| LJ/Coulomb | `E_LJ_*`, `E_Coul_*`, `E_int_*` (`_mean`, `_sd`, `_err`), `E_LJ_per_heavy` |
| Validation | `hbonds_mean`, `hbonds_sd`, `contacts_4A_mean`, `mindist_mean_nm`, `frac_unbound`, `lig_rmsd_mean_nm`, `lig_rmsd_last_nm` |
| Quality | `flags`, `notes` |
| Provenance | `gromacs`, `gmx_MMPBSA`, `date`, `path` |

All energies are in **kcal/mol**.

## 9. Quality flags

| Flag | Condition | What to do |
|---|---|---|
| `OK` | No problem found | — |
| `C2_UNRELIABLE` | σ(dGgas) > 3.6 kcal/mol | Report `dH`, not `dG(C2)` |
| `HIGH_SEM` | SEM of `dH` > 1 kcal/mol | More frames or a longer window |
| `POSITIVE_dH` | `dH` > 0 | Unfavourable binding: check the trajectory |
| `LIGAND_UNBOUND` | Ligand > 4.5 Å from the receptor in > 10 % of frames | The average does **not** represent the complex |
| `LIGAND_DISPLACED` | Mean ligand RMSD > 1 nm | The ligand moved to another site |
| `FEW_FRAMES` | Fewer than 100 frames | Increase `--frames` or the window |

## 10. How to read and report the results

- **`dH` (MM-GBSA)** is the main quantity for comparing affinities between
  ligands.
- **`E_LJ`** is a direct comparison of dispersion contacts, fair between
  charged and neutral ligands.
- **`E_Coul` and `E_int`** are informative only. They are computed in vacuum:
  for charged ligands (oleate, arachidonate, LPA, PIP2...) they are very
  large and must **not** be used alone to claim stronger binding.
- **`-TdS (C2)` and `dG`**: use them only if the `C2_UNRELIABLE` flag is
  absent.
- **`LE` and `E_LJ/HA`** are normalized per heavy atom, useful when ligands
  have very different sizes (e.g. LDPE, PET and PC hexamers).
- Do **not** call `E_LJ` or `E_int` a "binding free energy": they are
  interaction energies.
- Before using a number, check the flags. A system with `LIGAND_UNBOUND` or
  `LIGAND_DISPLACED` needs to be inspected (Basic Analyses RMSD, snapshots)
  before its average is interpreted.

## 11. collect_free_energy.py: one Excel for all systems

Run it **once**, from anywhere, giving the parent folders of your projects:

```bash
python3 collect_free_energy.py ~/Documents/Microplastics \
        ~/Documents/Endogenous_Ligands -o Energies_round2.xlsx
```

- It searches every `.../Analyses/FreeEnergy_*/summary.csv` inside those
  folders, at any depth. There is no need to move or merge folders.
- The **name of each parent folder** becomes the `Group` column, so systems
  from different projects can be compared in the same table. Use
  `--labels A,B` for custom group names (one per folder, in order).
- **Where the file is saved:** in the folder where you run the command. Without
  `-o` it is called `FreeEnergy_summary_<date>.xlsx`; with `-o` you choose
  name and path.

### Workbook sheets

| Sheet | Content |
|---|---|
| **Summary** | Paper-ready table: mean ± SEM (kcal/mol) |
| **Data_kcal** | All numbers (mean, SD, SEM) in kcal/mol |
| **Data_kJ** | The same in kJ/mol |
| **Top_residues** | Key residues per system |
| **Metadata** | Window, frames, temperature, versions, notes, path |
| **Legend** | Meaning of every column and every flag |

### openpyxl

The Excel file needs `openpyxl`. Check with:

```bash
python3 collect_free_energy.py --check-deps
```

If it is missing, the script:

1. looks for a conda env whose Python has it (`gmxMMPBSA` first) and re-runs
   itself with that Python;
2. otherwise prints an alarm with the install commands:
   ```
   python3 -m pip install --user openpyxl          # Ubuntu 22.04
   sudo apt install python3-openpyxl               # Ubuntu 23.04+
   conda activate gmxMMPBSA && pip install openpyxl
   ```
3. in a terminal, offers to install it right away;
4. if it is not installed, writes the same tables as CSV files.

`--no-install` skips the question (useful in scripts or under `nohup`).
`free_energy.sh` also warns at the end if `openpyxl` is missing.

## 12. Running many systems at once

**Step 0.** Once, in a terminal, choose the conda env (it is remembered):

```bash
cd <one_system>/gromacs && bash free_energy.sh --check
```

**Step 1.** Check detection in every system (seconds):

```bash
for d in $(find /path/Project_A /path/Project_B -name "dynamic*ns.xtc" -printf '%h\n' | sort -u); do
    cp ~/scripts/free_energy.sh "$d"/
    (cd "$d" && bash free_energy.sh --check)
done
```

**Step 2.** If everything looks right, compute all of them in the background:

```bash
nohup bash -c 'for d in $(find /path/Project_A /path/Project_B \
    -name "dynamic*ns.xtc" -printf "%h\n" | sort -u); do
    (cd "$d" && bash free_energy.sh); done' > free_energy_all.log 2>&1 &
```

**Step 3.** Collect:

```bash
python3 collect_free_energy.py /path/Project_A /path/Project_B
```

Tip: before launching many systems, run one completely and check
`MMGBSA/gmx_MMPBSA.out`.

## 13. Troubleshooting

| Message / problem | What to do |
|---|---|
| `AUTO-DETECTION STOPPED: ...` | Read the message: it says which option to add (`--ligand`, `--receptor`) |
| `System net charge is not zero` | Check `topol.top` (missing ions?) |
| `conda was not found` | Run `conda init bash` and open a new terminal, or use `bash free_energy.sh lj` |
| `conda env '...' does not contain gmx_MMPBSA` | Choose again: `bash free_energy.sh --env ask` |
| `Cannot ask for the conda env (no terminal...)` | Run once in a terminal, or pass `--env NAME` |
| gmx_MMPBSA fails | See `MMGBSA/gmx_MMPBSA.out` and `gmx_MMPBSA.log` |
| MPI fails | The script retries in serial automatically; or use `--np 1` |
| `Invalid window` | `--start`/`--end` outside the run length |
| Redo a step | Delete its output (`MMGBSA/FINAL_RESULTS_MMPBSA.dat`, `LJ/lj_coul.xvg` or `window.xtc`) and run again |

**Validation note:** the scripts were tested with the real topology of the
gastrotropin–oleate complex and with simulated GROMACS and gmx_MMPBSA
output.

---

**Authors:** Ricardo Aguilar-Garay (ricardoaguilargaray@gmail.com) ·
contributor: Dulce Estefanía Nicolás Álvarez, Departamento de Fisiología,
ENCB-IPN (dnicolas@ipn.mx). MIT License. Please cite the reference article
listed in the [main README](../README.md#citing).
