# Basic Analyses — detailed documentation

[← Back to the repository overview](../README.md) ·
[← Autogromacs](Autogromacs.md) ·
[Free Energy →](free_energy.md)

**Standard post-processing of Autogromacs trajectories, selecting GROMACS
groups by name instead of by number.**

`basic_analyses_with_names.sh` is the companion of
[`Autogromacs.sh`](Autogromacs.md). Copy this single file into the CHARMM-GUI
`gromacs/` folder after production finishes and run it. It produces the
usual structural-stability and interaction analyses (centered trajectory,
RMSD, hydrogen bonds or contacts, RMSF, radius of gyration and PDB
snapshots), with no interactive menus and nothing to edit.

---

## Contents

1. [Why "with names"](#1-why-with-names)
2. [Requirements](#2-requirements)
3. [Quick start](#3-quick-start)
4. [How the trajectory is detected](#4-how-the-trajectory-is-detected)
5. [The index file and the Ligand group](#5-the-index-file-and-the-ligand-group)
6. [The analyses, one by one](#6-the-analyses-one-by-one)
7. [Special cases](#7-special-cases)
8. [Output files](#8-output-files)
9. [Copying to shared project folders](#9-copying-to-shared-project-folders)
10. [Plotting the results](#10-plotting-the-results)
11. [Troubleshooting](#11-troubleshooting)

---

## 1. Why "with names"

The original lab script sent group **numbers** to the GROMACS tools, for
example `echo 4 13 | gmx rms ...`. Those numbers depend on what the system
contains (ligand, ions, lipids...) and change from one system to another, so
the same script could quietly analyse the wrong group.

This version builds its own index file and always passes groups **by name**
(`Protein`, `Backbone`, `Ligand`, `System`), so the same file gives the
correct result on every system.

| Analysis | Old script (numbers) | This script (names) |
|---|---|---|
| Centered trajectory | `1 0` | `Protein`, `System` |
| RMSD | `4 13` | `Backbone`, `Ligand` |
| H-bonds | `1 13` | `Protein`, `Ligand` |
| RMSF | `1` | `Protein` |
| Radius of gyration | `1` | `Protein` |
| Snapshots | `0` | `System` |

## 2. Requirements

- **GROMACS**, either `gmx` (thread-MPI build) or `gmx_mpi` (MPI build). It
  is taken from the `PATH` or loaded from a standard `GMXRC`
  (`/usr/local/gromacs/bin`, `/usr/local/gromacs/*/bin`, `/opt/gromacs/bin`,
  `~/gromacs/bin`). No GPU is needed.
- A **finished production run** in the current folder:
  `dynamic<N>ns.{tpr,xtc}` (Autogromacs output) or
  `step5_production*.{tpr,xtc}` (a CHARMM-GUI run made by hand).
- The CHARMM-GUI **`index.ndx`**, which provides the `SOLU` group.

## 3. Quick start

```bash
cd /path/to/my_system/gromacs
cp /path/to/basic_analyses_with_names.sh .

bash basic_analyses_with_names.sh                 # auto-detects the trajectory
bash basic_analyses_with_names.sh dynamic150ns    # or name the basename
```

For long trajectories you can leave it running in the background:

```bash
nohup bash basic_analyses_with_names.sh > analyses.log 2>&1 &
```

At start-up the script prints what it found, for example:

```
✓ GROMACS: /usr/local/gromacs/bin/gmx
═══ ComplexA · dynamic200ns.xtc · 200 ns → Analyses/ ═══
✓ Ligand group created (SOLU minus Protein).
```

## 4. How the trajectory is detected

In this order:

1. **Basename given as argument** (`dynamic150ns` → `dynamic150ns.tpr` +
   `dynamic150ns.xtc`).
2. The **newest** `dynamic*ns.tpr` that has a matching `.xtc`.
3. `step5_production.xtc` with the newest `step5_production*.tpr`.

If none is found, or `index.ndx` is missing, the script stops with a clear
message.

**Simulation length.** It is taken from the `dynamic<N>ns` name. For other
names it is read from the `.tpr` (`gmx dump`, `nsteps × dt`). The length is
used in the output file names (`_200ns`) and to place the snapshots.

## 5. The index file and the Ligand group

The script never modifies the CHARMM-GUI `index.ndx`. It creates its own
`Analyses/index_analysis.ndx`:

1. `gmx make_ndx` writes the GROMACS default groups (`System`, `Protein`,
   `Backbone`, `C-alpha`, ...).
2. The `SOLU` group is copied from the CHARMM-GUI `index.ndx`. `SOLU` holds
   everything that is not water or ions: protein plus ligand(s).
3. A new group is defined as `"SOLU" & ! "Protein"` and named **`Ligand`**.

So `Ligand` is **everything in `SOLU` that is not protein**, whatever the
residue is called (a small molecule, an oligomer of a microplastic, a
cofactor, several copies of a ligand...). If that group comes out empty, the
system has no ligand and the script switches to the protein-only analyses
(see [Special cases](#7-special-cases)).

## 6. The analyses, one by one

All time axes are in **ns**. Analyses 2–6 use the centered trajectory.

### 6.1 Centered trajectory

```
gmx trjconv -pbc mol -center    (center: Protein · output: System)
```

Makes molecules whole across the periodic boundaries and puts the protein at
the center of the box, so the ligand does not appear to "jump" from one side
of the box to the other. It is the input for every other analysis and is
computed **only once**: if `Analyses/<base>_centered.xtc` exists, it is
reused.

### 6.2 Ligand RMSD

```
gmx rms    (fit: Backbone · calc: Ligand)
```

Each frame is superimposed on the protein backbone of the reference
structure, and then the RMSD of the ligand is measured. This tells you
whether the ligand **stays in its binding pose** relative to the protein: a
low, flat curve means a stable pose; a jump means the ligand moved to
another site or left the pocket.

Without a ligand, the backbone RMSD is computed instead (overall protein
stability).

### 6.3 Protein–ligand hydrogen bonds

```
gmx hbond -r Protein -t Ligand     (new syntax, GROMACS 2024+)
```

Number of hydrogen bonds between protein and ligand in each frame. If the
new syntax is not available, the script falls back to the legacy interactive
selection. A persistent number of H-bonds supports a specific interaction;
zero for most of the run suggests a weak or non-specific binding.

### 6.4 RMSF per residue

```
gmx rmsf -res    (Protein)
```

Average fluctuation of each residue around its mean position. Peaks are
usually loops and termini; changes in the binding-site residues when
comparing complexes (or complex vs. apo protein) point to residues affected
by the ligand.

### 6.5 Radius of gyration

```
gmx gyrate    (Protein)
```

Compactness of the protein over time. A stable value means the fold is
kept; a steady rise means unfolding or a large opening of domains.

### 6.6 PDB snapshots

```
gmx trjconv -dt <total/2>    (System)
```

Writes three frames (start, middle and end of the run) to a single
multi-model PDB file. Open it in PyMOL, VMD or ChimeraX to see how the
complex changed and to prepare figures.

## 7. Special cases

**No ligand (apo protein).** The `Ligand` group is empty, so RMSD is computed
on the `Backbone` and H-bonds within the `Protein`. The script reports:
`✓ No ligand in SOLU: RMSD on Backbone, H-bonds within Protein.`

**Ligands with no N, O, S or F atoms** (polyethylene, polypropylene,
polystyrene and other C/H-only microplastics). They have no H-bond donors or
acceptors, and `gmx hbond` refuses to run. The script reads the ligand atoms
from the structure, sees that none is N/O/S/F, and computes
**protein–ligand contacts** instead:

```
gmx mindist -d 0.4    (Protein · Ligand)
```

- `contacts_<name>_<N>ns.xvg`: number of protein–ligand atom pairs closer
  than 0.4 nm in each frame (how much the ligand touches the protein).
- `mindist_<name>_<N>ns.xvg`: minimum protein–ligand distance in each frame
  (whether the ligand stays bound or drifts away).

**A failing H-bond or contact analysis does not stop the script.** It prints
a warning and carries on with RMSF, radius of gyration and snapshots.

## 8. Output files

Everything goes to `Analyses/` inside the `gromacs/` folder:

```
Analyses/
├── index_analysis.ndx            # default groups + SOLU + Ligand
├── dynamic<N>ns_centered.xtc     # centered trajectory (reused on re-runs)
├── rmsd_<name>_<N>ns.xvg         # ligand RMSD (backbone RMSD if no ligand)
├── hbond_<name>_<N>ns.xvg        # protein–ligand H-bonds
│   or contacts_<name>_<N>ns.xvg  #   contacts < 0.4 nm   (C/H-only ligands)
│      mindist_<name>_<N>ns.xvg   #   minimum distance    (C/H-only ligands)
├── rmsf_<name>_<N>ns.xvg         # RMSF per residue
├── gyrate_<name>_<N>ns.xvg       # radius of gyration
└── trj_<name>_<N>ns.pdb          # snapshots: start / middle / end
```

`<name>` is the current folder name (usually `gromacs`) and `<N>` the
simulation length in ns. Running the script again overwrites the `.xvg` and
`.pdb` files but reuses the centered trajectory.

## 9. Copying to shared project folders

Many projects in the lab keep the results of all systems side by side:

```
Project/
├── rmsd/  hbond/  contacts/  rmsf/  gyrate/  conformation/
├── ComplexA/gromacs/     ← run the script here
└── ComplexB/gromacs/
```

If any of these folders exists **two levels up** from `gromacs/`
(`../../rmsd`, `../../hbond`, `../../contacts`, `../../rmsf`,
`../../gyrate`, `../../conformation`), each result is also copied there,
renamed with the **system name** (the folder that contains `gromacs/`):

```
../../rmsd/rmsd_ComplexA_200ns.xvg
../../conformation/trj_ComplexA_200ns.pdb
```

so results from different systems never overwrite each other. Folders that
do not exist are simply skipped; nothing is created.

## 10. Plotting the results

The `.xvg` files are plain text (lines starting with `#` or `@` are
comments). Some quick options:

```bash
xmgrace Analyses/rmsd_gromacs_200ns.xvg        # Grace, if installed
```

```python
import numpy as np, matplotlib.pyplot as plt
t, y = np.loadtxt("Analyses/rmsd_gromacs_200ns.xvg", comments=("#", "@"), unpack=True)
plt.plot(t, y); plt.xlabel("Time (ns)"); plt.ylabel("RMSD (nm)"); plt.show()
```

They can also be imported directly into Excel, Origin or GraphPad by
skipping the comment lines.

## 11. Troubleshooting

**`neither 'gmx' nor 'gmx_mpi' found`** — load GROMACS first
(`source /path/to/gromacs/bin/GMXRC`).

**`no trajectory found`** — there is no `dynamic<N>ns.{tpr,xtc}` or
`step5_production.*` in the folder, or the `.xtc` does not match the `.tpr`
name. Pass the basename explicitly.

**`index.ndx not found`** — run the script inside the CHARMM-GUI `gromacs/`
folder.

**The Ligand group contains something unexpected** (e.g. a cofactor you did
not want) — the group is "everything in SOLU that is not protein". Check
`Analyses/index_analysis.ndx`, or edit the `Ligand` group there and run the
specific `gmx` command by hand.

**Snapshots look broken across the box** — use the centered trajectory
(`Analyses/*_centered.xtc`) for any extra analysis, not the raw `.xtc`.

---

**Authors:** Ricardo Aguilar-Garay (ricardoaguilargaray@gmail.com) ·
contributor: Dulce Estefanía Nicolás Álvarez, Departamento de Fisiología,
ENCB-IPN (dnicolas@ipn.mx). MIT License. Please cite the reference article
listed in the [main README](../README.md#citing).
