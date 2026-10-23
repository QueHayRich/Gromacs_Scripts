# Gromacs_Scripts — Autogromacs & Basic Analyses

**Automated GROMACS molecular dynamics for CHARMM-GUI systems: from the
downloaded `gromacs/` folder to ready-to-plot results with two commands.**

Developed at the Laboratorio de Hormonas y Conducta, Escuela Nacional de
Ciencias Biológicas, Instituto Politécnico Nacional (ENCB-IPN), Mexico.

```
CHARMM-GUI gromacs/ ──► Autogromacs.sh ──► dynamic<N>ns.xtc ──► basic_analyses_with_names.sh ──► Analyses/*.xvg, *.pdb
                        (minimization,                          (RMSD, H-bonds/contacts, RMSF,
                         equilibration,                          radius of gyration, snapshots)
                         production)
```

---

## Tools

### 🧬 Autogromacs — run the simulation

`Autogromacs.sh` runs the whole MD workflow (energy minimization →
equilibration → production) on an NVIDIA GPU for any length you ask for.

- **One file, nothing to edit:** copy it into the CHARMM-GUI `gromacs/`
  folder and give the length in ns.
- **Any length:** production runs in 100 ns segments appended to **one
  continuous trajectory** (`dynamic<N>ns.xtc`), no `trjcat` needed.
- **Resumable:** after a crash, power cut or Ctrl-C, run it again and it
  continues from the last checkpoint.
- **GPU failure recovery:** if the full-GPU mode fails (e.g. CUDA error
  #700), it retries automatically in a stable mode.
- **Portable:** same file on RTX 2080 / 3080 / 5090, GROMACS 2024 and 2026.

```bash
nohup bash Autogromacs.sh 200 > autogromacs.log 2>&1 &
```

📖 **[Full Autogromacs documentation →](docs/Autogromacs.md)**

### 📊 Basic Analyses — analyse the trajectory

`basic_analyses_with_names.sh` runs the standard post-production analyses on
the Autogromacs output, selecting GROMACS groups **by name** so it works on
any system without editing.

- **Centered trajectory**, **ligand RMSD**, **protein–ligand hydrogen
  bonds**, **RMSF per residue**, **radius of gyration** and **PDB snapshots**
  (start / middle / end).
- **Any ligand:** the ligand is detected automatically as everything in
  CHARMM-GUI's `SOLU` group that is not protein.
- **Microplastics and other C/H-only ligands:** H-bonds are replaced
  automatically by protein–ligand contacts and minimum distance.
- **Organized output** in `Analyses/`, optionally copied to shared project
  folders with the system name.

```bash
bash basic_analyses_with_names.sh
```

📖 **[Full Basic Analyses documentation →](docs/basic_analyses.md)**

---

## Quick start

```bash
cd /path/to/my_system/gromacs            # CHARMM-GUI Solution Builder output
cp /path/to/Gromacs_Scripts/*.sh .

nohup bash Autogromacs.sh 200 > autogromacs.log 2>&1 &   # 1. simulate 200 ns
bash basic_analyses_with_names.sh                         # 2. analyse
```

**Requirements:** Linux, GROMACS 2024/2026 with CUDA, an NVIDIA GPU (for the
simulation; the analyses run on CPU) and a CHARMM-GUI Solution Builder
`gromacs/` folder. Details in each tool's documentation.

## Roadmap

- **Short-range Lennard-Jones / Coulomb interaction energies**
  (`gmx energy` on energygrps).
- **MM-PBSA / MM-GBSA module**: binding free energy (ΔG_bind ± SEM) with
  per-residue decomposition, following the protocol validated in the
  reference publication below.
- Support for CHARMM-GUI Membrane Builder systems.

## Authors and contributors

- **Ricardo Aguilar-Garay**: author and maintainer
  (ricardoaguilargaray@gmail.com)
- **Dulce Estefanía Nicolás Álvarez**: contributor. Departamento de
  Fisiología, Escuela Nacional de Ciencias Biológicas, Instituto Politécnico
  Nacional (ENCB-IPN), Mexico (dnicolas@ipn.mx)

Questions, bug reports and suggestions are welcome.

## Citing

If this pipeline is useful in your research, **please cite the article in
which the underlying MD + MM-PBSA methodology was developed and validated**:

> Hernandez-Munguia, A.M.; Reyes-Chaparro, A.; Fregoso-Aguilar, T.A.;
> Tenorio-Barajas, A.Y.; Mendoza-Pérez, J.A.; **Aguilar-Garay, R.**;
> Nicolás-Álvarez, D.E. Binding Affinity and Interaction Profiles of
> Erinacines and Erinacerins with iNOS and NF-κB Revealed by Molecular
> Dynamics Simulations. *Int. J. Mol. Sci.* **2026**, *27*(7), 3145.
> https://doi.org/10.3390/ijms27073145

BibTeX:

```bibtex
@Article{ijms27073145,
  AUTHOR         = {Hernandez-Munguia, April Michelle and Reyes-Chaparro, Andrés and Fregoso-Aguilar, Tomas Alejandro and Tenorio-Barajas, Aldo Yair and Mendoza-Pérez, Jorge Alberto and Aguilar-Garay, Ricardo and Nicolás-Álvarez, Dulce E.},
  TITLE          = {Binding Affinity and Interaction Profiles of Erinacines and Erinacerins with iNOS and NF-κB Revealed by Molecular Dynamics Simulations},
  JOURNAL        = {International Journal of Molecular Sciences},
  VOLUME         = {27},
  YEAR           = {2026},
  NUMBER         = {7},
  ARTICLE-NUMBER = {3145},
  URL            = {https://www.mdpi.com/1422-0067/27/7/3145},
  ISSN           = {1422-0067},
  DOI            = {10.3390/ijms27073145}
}
```

GitHub users can also use the **"Cite this repository"** button (powered by
`CITATION.cff` in this repo).

## License

Released under the **MIT License** (see `LICENSE`). You are free to use,
modify, and redistribute this software; attribution via the citation above
is kindly requested in any academic work that benefits from it.

## Acknowledgments

Developed with support from the SECIHTI postdoctoral program. GPU workflow
validated on gelsolin and microplastic–protein complex systems at ENCB-IPN.
