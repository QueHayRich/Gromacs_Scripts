# basic_analyses_with_names.sh

**Standard post-processing of Autogromacs trajectories, selecting GROMACS
groups by name instead of by number.**

Companion of [`Autogromacs.sh`](README.md). Copy this single file into the
CHARMM-GUI `gromacs/` folder after production finishes and run it. It
produces the usual quality and interaction analyses (RMSD, H-bonds, RMSF,
radius of gyration, snapshots) with no interactive menus and no editing.

---

## Why "with names"

The original lab script piped group **numbers** to the GROMACS tools
(`echo 4 13 | gmx rms ...`). Those numbers change from system to system
(depending on the ligand, ions, lipids, etc.), so the same script could
silently analyse the wrong group. This version builds its own index file and
always passes groups **by name** (`Protein`, `Backbone`, `Ligand`, `System`),
so the result is the same on every system.

The `Ligand` group is defined as **everything in the CHARMM-GUI `SOLU` group
that is not `Protein`**, so it works for any ligand name (small molecules,
microplastic oligomers, lipids in SOLU, ...).

## Requirements

- GROMACS (`gmx` or `gmx_mpi`) — found in the `PATH` or loaded automatically
  from a standard `GMXRC` (`/usr/local/gromacs`, `/opt/gromacs`,
  `~/gromacs`).
- A finished production run in the current folder:
  `dynamic<N>ns.{tpr,xtc}` (Autogromacs output) or
  `step5_production*.{tpr,xtc}` (plain CHARMM-GUI).
- The CHARMM-GUI `index.ndx` (it provides the `SOLU` group).

## Usage

```bash
cd /path/to/system/gromacs
bash basic_analyses_with_names.sh                 # auto-detects the trajectory
bash basic_analyses_with_names.sh dynamic150ns    # or give the basename
```

Trajectory detection order: command-line basename → newest
`dynamic<N>ns.tpr` that has a matching `.xtc` → `step5_production.xtc`.

The simulation length (used in file names and snapshot spacing) is taken from
the `dynamic<N>ns` name or, otherwise, read from the `.tpr` (`gmx dump`).

## Analyses

| # | Analysis | Groups | GROMACS tool | Output |
|---|----------|--------|--------------|--------|
| 1 | Centered trajectory | center `Protein`, write `System` | `trjconv -pbc mol -center` | `<base>_centered.xtc` |
| 2 | Ligand RMSD | fit `Backbone`, calc `Ligand` | `rms` | `rmsd_<name>_<N>ns.xvg` |
| 3 | Protein–ligand H-bonds | `Protein` – `Ligand` | `hbond` | `hbond_<name>_<N>ns.xvg` |
| 4 | RMSF per residue | `Protein` | `rmsf -res` | `rmsf_<name>_<N>ns.xvg` |
| 5 | Radius of gyration | `Protein` | `gyrate` | `gyrate_<name>_<N>ns.xvg` |
| 6 | PDB snapshots (start / middle / end) | `System` | `trjconv` | `trj_<name>_<N>ns.pdb` |

All time axes are in **ns**. Steps 2–6 use the centered trajectory, which is
computed only once (it is reused if it already exists).

### Special cases

- **No ligand** (apo protein): RMSD is computed on the `Backbone` and H-bonds
  within the `Protein`.
- **Ligands without N/O/S/F atoms** (e.g. polyethylene, polypropylene,
  polystyrene microplastics): they have no H-bond donors/acceptors and
  `gmx hbond` refuses to run. The script detects this and computes
  **protein–ligand contacts** instead (`gmx mindist`, atom pairs < 0.4 nm):
  `contacts_<name>_<N>ns.xvg` and `mindist_<name>_<N>ns.xvg`.
- **GROMACS versions**: H-bonds use the new `gmx hbond -r/-t` syntax and fall
  back to the legacy interactive selection on older versions. If the H-bond
  or contact analysis fails, the script warns and continues with the rest.

## Output

Everything is written to `Analyses/` inside the `gromacs/` folder:

```
Analyses/
├── index_analysis.ndx            # default groups + SOLU + Ligand
├── dynamic<N>ns_centered.xtc
├── rmsd_<name>_<N>ns.xvg
├── hbond_<name>_<N>ns.xvg        # or contacts_/mindist_ for C/H-only ligands
├── rmsf_<name>_<N>ns.xvg
├── gyrate_<name>_<N>ns.xvg
└── trj_<name>_<N>ns.pdb
```

`<name>` is the current folder name (usually `gromacs`).

### Copy to shared folders

If the project uses the lab's shared-folder layout, i.e. these folders exist
two levels up:

```
../../rmsd  ../../hbond  ../../contacts  ../../rmsf  ../../gyrate  ../../conformation
```

each result is also copied there, renamed with the **system name** (the
folder that contains `gromacs/`), e.g. `../../rmsd/rmsd_ComplexA_200ns.xvg`,
so results from different systems do not overwrite each other. Folders that
do not exist are simply skipped.

## Authors

- **Ricardo Aguilar-Garay** — author (ricardoaguilargaray@gmail.com)
- **Dra. Dulce Estefanía Nicolás Álvarez** — contributor. Departamento de
  Fisiología, Escuela Nacional de Ciencias Biológicas, Instituto Politécnico
  Nacional (ENCB-IPN), Mexico (dnicolas@ipn.mx)

Questions, bug reports and suggestions are welcome. Released under the MIT
License (see `LICENSE`); please cite the reference article listed in the
main [`README.md`](README.md) if this script is useful in your research.
