# Autogromacs

**Automated, resumable, chunked GROMACS molecular dynamics pipeline for CHARMM-GUI systems on NVIDIA GPUs.**

A single self-contained Bash script that takes a freshly downloaded CHARMM-GUI
`gromacs/` folder and runs the entire workflow — energy minimization →
equilibration → production MD — for any requested simulation length, on any
lab machine, with automatic crash recovery.

Developed at the Laboratorio de Hormonas y Conducta, Escuela Nacional de
Ciencias Biológicas, Instituto Politécnico Nacional (ENCB-IPN), Mexico.

---

## Key features

- **Single file, zero configuration.** Copy `Autogromacs.sh` into the
  CHARMM-GUI `gromacs/` folder and run it. No editing of `.mdp` files, no
  manual `nsteps` arithmetic.
- **Arbitrary simulation length via chunked production.** Production runs in
  100 ns segments: the script runs a segment, extends the run input file with
  `gmx convert-tpr -extend`, and continues with `gmx mdrun -cpi -append`
  (`-cpi/-append` are only passed once a checkpoint exists, so the very first
  segment starts cleanly on every GROMACS version).
  Because output is appended to the same files, the final trajectory is a
  **single continuous `.xtc`** with no duplicated frames and no need for
  `gmx trjcat`.
- **Crash- and power-failure-proof (resumable).** If the run is interrupted
  for any reason, re-running the script (no arguments needed) resumes from
  the last GROMACS checkpoint (`.cpt`). The requested total length is
  persisted in a hidden state file, and the *authoritative* simulation target
  is always read back from the current `.tpr` itself (`gmx dump`), so no
  desynchronization between script state and simulation state is possible.
- **Clean handling of interruptions.** Ctrl-C, `kill`, or a shutdown make
  `mdrun` stop gracefully and write a checkpoint. The script recognizes this
  (signal trap + the `Received the ... signal` line in the current segment of
  the `mdrun` log) and exits with code 130 instead of mistaking it for a GPU
  failure; re-running the script resumes the run.
- **Maximum-performance GPU offload with automatic fallback.** Each
  production segment first attempts the full-offload configuration
  (`-nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu`). If the run fails
  mid-segment (e.g., CUDA error #700 / GPU instability), the script
  automatically restarts **from the last checkpoint** in a stable
  configuration (`-bonded cpu -update cpu`, set explicitly because recent
  GROMACS versions keep `update` on the GPU under `auto`) without any user
  intervention.
- **Multi-version GROMACS support.** Designed for GROMACS 2024 and 2026 builds
  (tested on 2026.4).
  The detected version is printed at start-up.
- **No experimental features in production.** If `GMX_CUDA_GRAPH` is set in
  the environment (e.g., from `~/.bashrc`), it is unset for the run, since
  GROMACS flags CUDA graphs as experimental.
- **Hardware-portable across the lab.** The script locates the GROMACS
  installation (PATH first, then standard `GMXRC` locations, including
  `~/software/gromacs-*/`), verifies that
  an NVIDIA GPU is present (aborting cleanly otherwise), and sets OpenMP
  threads from `nproc`. The same file runs unmodified on an RTX 2080, 3080,
  or 5090 workstation.
- **Idempotent stage skipping.** Minimization and equilibration are skipped
  automatically when their output already exists, so the script can be
  re-run safely at any point.
- **Safe file handling.** `convert-tpr` writes to a temporary file that is
  atomically moved into place (an interrupted extension can never corrupt
  the run input). An existing final trajectory is never overwritten.
- **Self-documenting output naming.** On completion, all production files are
  renamed to `dynamic<N>ns.*` (e.g., `dynamic200ns.xtc`, `dynamic200ns.tpr`),
  ready for downstream analysis tooling.
- **Background GPU telemetry.** `nvidia-smi dmon` logs power, utilization,
  clocks, and temperature to `gpu_monitor.log` for the duration of the run.

## Requirements

| Component | Requirement |
|---|---|
| OS | Linux (tested on Ubuntu) |
| GROMACS | 2024 or 2026, GPU-enabled (CUDA) build; `gmx` on `PATH` or `GMXRC` in a standard location (`/usr/local/gromacs/bin/GMXRC` preferred) |
| GPU | NVIDIA GPU with working driver (`nvidia-smi`) — **mandatory** |
| Input system | CHARMM-GUI **Solution Builder** output (`gromacs/` folder) |

> Note for RTX 5090 (Blackwell) machines: GROMACS must be compiled against
> CUDA ≥ 12.8. The script itself needs no changes.

## Expected input files (CHARMM-GUI naming)

```
step3_input.gro              # solvated, neutralized system
topol.top                    # topology (with toppar/ folder)
index.ndx                    # index groups
step4.0_minimization.mdp
step4.1_equilibration.mdp
step5_production.mdp
```

The script verifies these before starting and aborts with a clear message if
any are missing. Membrane Builder systems (`step6.x`/`step7` naming) are not
yet supported.

## Usage

```bash
# Interactive (asks for the total production length in ns):
bash Autogromacs.sh

# Non-interactive (recommended for long runs under nohup):
nohup bash Autogromacs.sh 200 > autogromacs.log 2>&1 &

# After a power failure / crash / Ctrl-C — just run it again, it resumes:
bash Autogromacs.sh
```

## How it works

1. **Environment checks** — locate `gmx`, verify NVIDIA GPU, count CPU
   threads, verify CHARMM-GUI input files.
2. **Target registration** — the requested total (argument or prompt) is
   stored in `.autogromacs_total_ns`; on resume it is read back so the user
   is never asked twice (and cannot accidentally answer inconsistently).
3. **Minimization / equilibration** — standard CHARMM-GUI protocol
   (position restraints referenced to `step3_input.gro`), skipped if outputs
   exist. Equilibration runs on GPU in the stable configuration.
4. **Production, chunked**:
   - First segment: `grompp` builds `step5_production.tpr` (reading the exact
     equilibrated state from `step4.1_equilibration.cpt` with `-t`, as
     required by `continuation = yes`), then
     `convert-tpr -nsteps` trims it to min(total, 100) ns.
   - Loop: `mdrun` (with `-cpi -append` once a checkpoint exists; full GPU
     offload; automatic stable-mode retry on failure, clean exit on
     interruption) → read achieved target from the `.tpr` → if the total
     is not yet reached, `convert-tpr -extend` by the next chunk → repeat.
   - The timestep `dt` is parsed from the production `.mdp` (fallback
     0.002 ps), so step arithmetic is always consistent with the actual run.
5. **Finalization** — production files are renamed to `dynamic<N>ns.*` and
   the state file is removed.

## Output

```
dynamic<N>ns.xtc    # single continuous production trajectory
dynamic<N>ns.tpr    # run input matching the full length
dynamic<N>ns.edr    # energies
dynamic<N>ns.log    # mdrun log
dynamic<N>ns.gro    # final coordinates
dynamic<N>ns.cpt    # final checkpoint
gpu_monitor.log     # GPU telemetry during the run
```

## Roadmap (in development)

- **`basic_analyses.sh`** — companion post-processing script that
  auto-detects `dynamic<N>ns.*` output and produces the standard quality and
  interaction analyses used in our published workflow: backbone RMSD, RMSF,
  radius of gyration, hydrogen bonds, ligand–protein minimum distance, and
  **short-range Lennard-Jones / Coulomb interaction energies**
  (`gmx energy` on energygrps).
- **MM-PBSA module** — binding free energy calculation (ΔG_bind ± SEM) with
  per-residue energy decomposition, following the protocol validated in the
  reference publication below. *This module will be updated in a future
  release.*

## Citing this pipeline

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
