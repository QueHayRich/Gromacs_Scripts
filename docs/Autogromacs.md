# Autogromacs — detailed documentation

[← Back to the repository overview](../README.md) ·
[Basic Analyses →](basic_analyses.md)

**Automated, resumable, chunked GROMACS molecular dynamics pipeline for
CHARMM-GUI systems on NVIDIA GPUs.**

`Autogromacs.sh` is a single self-contained Bash script. You copy it into a
freshly downloaded CHARMM-GUI `gromacs/` folder and it runs the whole
workflow (energy minimization → equilibration → production MD) for any
simulation length you ask for. It runs unchanged on any lab machine and
recovers automatically from crashes, power failures and GPU errors.

---

## Contents

1. [Why Autogromacs](#1-why-autogromacs)
2. [Requirements](#2-requirements)
3. [Expected input files](#3-expected-input-files)
4. [Quick start](#4-quick-start)
5. [What the script does, step by step](#5-what-the-script-does-step-by-step)
6. [Chunked production in detail](#6-chunked-production-in-detail)
7. [Resuming after an interruption](#7-resuming-after-an-interruption)
8. [GPU modes and automatic fallback](#8-gpu-modes-and-automatic-fallback)
9. [Output files](#9-output-files)
10. [Messages and exit codes](#10-messages-and-exit-codes)
11. [Troubleshooting](#11-troubleshooting)
12. [FAQ](#12-faq)
13. [Limitations](#13-limitations)

---

## 1. Why Autogromacs

Running a CHARMM-GUI system by hand means typing a chain of `grompp` and
`mdrun` commands, working out `nsteps` for the length you want, choosing GPU
flags for each machine and, for long runs, extending and joining
trajectories. Doing all that by hand invites mistakes, and a power cut in the
middle of a 500 ns run can cost days.

Autogromacs automates all of it:

| Problem | What Autogromacs does |
|---|---|
| Editing `.mdp` files to change the simulation length | You give the length in ns; the script computes `nsteps` from the `dt` in the `.mdp`. |
| Long runs must be split and joined with `trjcat` | Production runs in 100 ns segments appended to **one continuous `.xtc`**. |
| A crash or power cut loses the run | Re-run the script with no arguments; it continues from the last checkpoint. |
| GPU flags differ between machines | The same flags work on RTX 2080 / 3080 / 5090; GROMACS adapts. |
| CUDA error #700 / GPU instability | Automatic retry from the checkpoint in a stable mode. |
| Ctrl-C mistaken for a GPU failure | Interruptions are detected and the script exits cleanly (code 130). |
| Different GROMACS versions | Works with GROMACS 2024 and 2026 (tested on 2026.4). |

## 2. Requirements

| Component | Requirement |
|---|---|
| OS | Linux (tested on Ubuntu) |
| GROMACS | 2024 or 2026, GPU-enabled (CUDA) build |
| GPU | NVIDIA GPU with a working driver (`nvidia-smi`): **mandatory** |
| Input system | CHARMM-GUI **Solution Builder** output (`gromacs/` folder) |

**How GROMACS is located.** The script first looks for `gmx` in the `PATH`.
If it is not there, it sources the first `GMXRC` found in:

```
/usr/local/gromacs/bin/GMXRC
/usr/local/gromacs/*/bin/GMXRC
/opt/gromacs/bin/GMXRC
~/software/gromacs-*/bin/GMXRC
~/gromacs/bin/GMXRC
```

The detected version and path are printed at start-up
(`✓ GROMACS 2026.4: /usr/local/gromacs/bin/gmx`).

> **RTX 5090 (Blackwell) machines:** GROMACS must be compiled against
> CUDA ≥ 12.8. The script itself needs no changes.

## 3. Expected input files

These are the standard file names of a CHARMM-GUI Solution Builder
`gromacs/` folder:

```
step3_input.gro              # solvated, neutralized system
topol.top                    # topology (plus the toppar/ folder)
index.ndx                    # index groups (SOLU, MEMB, SOLV, ...)
step4.0_minimization.mdp
step4.1_equilibration.mdp
step5_production.mdp
```

The script checks for `step5_production.mdp`, `topol.top` and `index.ndx`
before starting and stops with a clear message if one is missing, which
usually means it was not run inside the `gromacs/` folder.

## 4. Quick start

```bash
cd /path/to/my_system/gromacs
cp /path/to/Autogromacs.sh .

# Interactive: asks how many ns of production you want
bash Autogromacs.sh

# Non-interactive (recommended for long runs; survives closing the terminal)
nohup bash Autogromacs.sh 200 > autogromacs.log 2>&1 &

# Follow progress
tail -f autogromacs.log

# After a crash, power cut or Ctrl-C: run it again with no arguments
bash Autogromacs.sh
```

The length must be a whole number of ns ≥ 1 (e.g. `100`, `200`, `250`).
Under `nohup` there is no terminal to ask in, so on the first run the length
**must** be passed as an argument.

## 5. What the script does, step by step

```
 ┌───────────────────────┐
 │ 0. Environment checks │  gmx, GPU, CPU threads, input files
 └──────────┬────────────┘
            ▼
 ┌───────────────────────┐
 │ Target registration   │  argument > saved state > prompt
 └──────────┬────────────┘
            ▼
 ┌───────────────────────┐
 │ 1. Minimization       │  skipped if step4.0_minimization.gro exists
 └──────────┬────────────┘
            ▼
 ┌───────────────────────┐
 │ 2. Equilibration      │  skipped if step4.1_equilibration.gro exists
 └──────────┬────────────┘
            ▼
 ┌───────────────────────┐      ┌──────────────────────────┐
 │ 3. Production segment │ ───► │ total reached?           │
 │    (mdrun, ≤100 ns)   │      │  no → convert-tpr -extend│──┐
 └──────────▲────────────┘      │  yes → finish            │  │
            └───────────────────┴──────────────────────────┘◄─┘
            ▼
 ┌───────────────────────┐
 │ 4. Rename to          │  dynamic<N>ns.{xtc,tpr,edr,log,gro,cpt}
 │    dynamic<N>ns.*     │
 └───────────────────────┘
```

### 0. Environment checks

- Locates GROMACS (see [Requirements](#2-requirements)) and prints its
  version.
- If `GMX_CUDA_GRAPH` is set in the environment (for example from
  `~/.bashrc`), it is **unset** for this run, because GROMACS marks CUDA
  graphs as experimental. Only this script is affected.
- Checks that an NVIDIA GPU is present (`nvidia-smi -L`) and stops if not.
- Sets the number of OpenMP threads from `nproc`.

### Target registration

The total production length comes from, in this order:

1. the command-line argument (`bash Autogromacs.sh 200`);
2. the hidden state file `.autogromacs_total_ns`, written on the first run;
3. an interactive prompt.

The value is saved to `.autogromacs_total_ns`, so on a resume you are never
asked again and cannot type a different number by mistake.

The timestep `dt` is read from `step5_production.mdp` (fallback 0.002 ps)
and used for every ns ↔ steps conversion.

If `dynamic<N>ns.xtc` already exists in the folder, the script stops without
touching anything, so finished results are never overwritten.

### 1. Minimization

```
gmx grompp -f step4.0_minimization.mdp -c step3_input.gro -r step3_input.gro ...
gmx mdrun  -deffnm step4.0_minimization
```

Standard CHARMM-GUI protocol, with position restraints referenced to
`step3_input.gro`. Skipped if `step4.0_minimization.gro` already exists.

### 2. Equilibration

```
gmx grompp -f step4.1_equilibration.mdp -c step4.0_minimization.gro -r step3_input.gro ...
gmx mdrun  -deffnm step4.1_equilibration   (SAFE GPU mode)
```

Runs on the GPU in the stable configuration (see
[GPU modes](#8-gpu-modes-and-automatic-fallback)). Skipped if
`step4.1_equilibration.gro` already exists.

### 3. Production

Described in detail in the [next section](#6-chunked-production-in-detail).

### 4. Finalization

All `step5_production.*` files are renamed to `dynamic<N>ns.*`, the state
file is deleted and the script prints the next step: running
[Basic Analyses](basic_analyses.md).

## 6. Chunked production in detail

Production runs in **segments of 100 ns** (`CHUNK_NS=100`). For a 250 ns
run:

| Segment | `.tpr` target | What happens |
|---|---|---|
| 1 | 100 ns | `grompp` + `convert-tpr -nsteps` → `mdrun` |
| 2 | 200 ns | `convert-tpr -extend 100000` → `mdrun -cpi -append` |
| 3 | 250 ns | `convert-tpr -extend 50000` → `mdrun -cpi -append` |

**First segment.** `grompp` builds `step5_production.tpr` from the
equilibrated structure. If `step4.1_equilibration.cpt` exists it is passed
with `-t`, so the run starts from the exact equilibrated state (velocities
included), as `continuation = yes` in the CHARMM-GUI `.mdp` requires. Then
`gmx convert-tpr -nsteps` sets the length of that first segment to
min(total, 100) ns.

**Next segments.** When a segment ends and the total has not been reached,
`gmx convert-tpr -extend` adds the next segment (100 ns, or what remains) to
the `.tpr`, and `mdrun -cpi -append` continues from the checkpoint, writing
to the **same** `.xtc`, `.edr` and `.log`. The result is one continuous
trajectory with no duplicated frames and no need for `gmx trjcat`.

**Why `-cpi/-append` only when a checkpoint exists.** GROMACS 2026 aborts if
`-append` is given without a checkpoint. The script adds these flags only
when `step5_production.cpt` exists, so the very first segment starts cleanly
on every version.

**The `.tpr` is the source of truth.** After each segment the script reads
the current target back from the `.tpr` itself (`gmx dump`, `nsteps × dt`).
Even if a power cut happens between extending the `.tpr` and anything else,
script and simulation can never get out of step.

**Safe writes.** `convert-tpr` always writes to a temporary file
(`step5_production_tmp.tpr`) that is then moved into place, so an
interrupted extension cannot leave a half-written `.tpr`.

## 7. Resuming after an interruption

Just run the script again in the same folder, with no arguments:

```bash
bash Autogromacs.sh
# or
nohup bash Autogromacs.sh > autogromacs.log 2>&1 &
```

What happens depends on where the run stopped:

| Interrupted during | On resume |
|---|---|
| Minimization | Minimization runs again (it is short). |
| Equilibration | Minimization is skipped; equilibration runs again. |
| Production (any segment) | Minimization and equilibration are skipped; `mdrun` continues from `step5_production.cpt`. |
| Between two segments | The `.tpr` already holds the right target; the run continues. |
| After renaming | `dynamic<N>ns.xtc` exists, so the script reports it and exits. |

**Changing your mind about the length.** You can make a run longer by passing
a bigger number on resume (`bash Autogromacs.sh 300`). If the `.tpr` already
targets more than the number you pass, the script finishes up to the `.tpr`
target and names the output accordingly.

## 8. GPU modes and automatic fallback

Two sets of `mdrun` flags are used:

| Mode | Flags | Used for |
|---|---|---|
| **HARD** (maximum performance) | `-nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu -ntmpi 1 -ntomp <nproc> -pin on` | First attempt of every production segment |
| **SAFE** (stable) | `-nb gpu -pme gpu -pmefft gpu -bonded cpu -update cpu -ntmpi 1 -ntomp <nproc> -pin on` | Equilibration, and automatic retry when HARD fails |

In SAFE mode, `-bonded cpu -update cpu` are set **explicitly**, because
recent GROMACS versions keep `update` on the GPU when it is left on `auto`.

**Automatic fallback.** If a segment fails in HARD mode (for example CUDA
error #700 or another GPU instability), the script prints a warning and
retries the same segment in SAFE mode **from the last checkpoint**, with no
user action. If SAFE mode also fails, the script stops and points you to
`step5_production.log`: the system itself is probably unstable (see
[Troubleshooting](#11-troubleshooting)).

**Interruptions are not failures.** When you press Ctrl-C, `kill` the
process or shut the machine down, `mdrun` stops gracefully, writes a
checkpoint and exits with an error code, just like a real failure. To avoid
wrongly switching to SAFE mode, the script detects interruptions in two
ways:

1. a signal trap (`INT`, `TERM`, `HUP`) on the script itself;
2. the `Received the ... signal` line that `mdrun` writes to its log. Only
   the lines written in the current attempt are checked, so an old message
   from a previous run cannot cause a false alarm.

In either case it prints `⏹ mdrun interrumpido por señal` and exits with
code **130**. Run the script again to continue.

**GPU telemetry.** During production, `nvidia-smi dmon -s puct -d 10`
records power, utilization, clocks and temperature every 10 s to
`gpu_monitor.log`. It stops automatically when the script ends.

## 9. Output files

At the end, in the `gromacs/` folder:

```
dynamic<N>ns.xtc        # single continuous production trajectory
dynamic<N>ns.tpr        # run input for the full length (use it for analysis)
dynamic<N>ns.edr        # energies
dynamic<N>ns.log        # mdrun log (all segments, appended)
dynamic<N>ns.gro        # final coordinates
dynamic<N>ns.cpt        # final checkpoint (lets you extend the run later)
dynamic<N>ns_prev.cpt   # previous checkpoint (backup)
gpu_monitor.log         # GPU telemetry
step4.0_minimization.*  # minimization files
step4.1_equilibration.* # equilibration files
```

While the run is in progress, production files are called
`step5_production.*` and the hidden file `.autogromacs_total_ns` stores the
target.

## 10. Messages and exit codes

| Message | Meaning |
|---|---|
| `🔄 Reanudación detectada: objetivo guardado = N ns` | Resume: the target was read from the state file. |
| `✓ Minimización ya hecha, se omite.` | Stage skipped because its output exists. |
| `⚠ El modo HARD falló ... Reintentando en MODO SAFE` | GPU failure; automatic retry in stable mode. |
| `⏹ mdrun interrumpido por señal` | Ctrl-C / kill / shutdown; checkpoint saved. |
| `⚠ Ya existe dynamic<N>ns.xtc` | Run already finished; nothing is done. |
| `❌ mdrun falló incluso en modo estable` | Both modes failed; check the log. |

| Exit code | Meaning |
|---|---|
| `0` | Finished, or the final trajectory already existed. |
| `1` | Error (missing file, no GPU, `grompp`/`mdrun` failure). |
| `130` | Interrupted by a signal; run the script again to resume. |

## 11. Troubleshooting

**`no encuentro 'gmx'`** — GROMACS is not in the `PATH` or in a standard
location. Run `source /path/to/gromacs/bin/GMXRC` before the script, or add
it to `~/.bashrc`.

**`este pipeline requiere GPU NVIDIA`** — `nvidia-smi` does not work. Check
the driver (`nvidia-smi` in a terminal) and reboot if it was just updated.

**`falta 'step5_production.mdp'` (or `topol.top`, `index.ndx`)** — you are not
inside the CHARMM-GUI `gromacs/` folder.

**`Sin terminal interactiva y sin objetivo previo`** — under `nohup` the first
run needs the length as an argument: `nohup bash Autogromacs.sh 200 ...`.

**Fails in SAFE mode too** — look in `step5_production.log` for
`LINCS warnings`, `Segmentation fault` or a large positive potential energy.
These point to an unstable system (bad ligand parameters, atoms overlapping,
too large a timestep) rather than to the GPU. Check the minimization and
equilibration and the ligand topology.

**HARD mode fails on every segment** — the run still completes in SAFE mode,
only more slowly. Repeated CUDA #700 errors may point to overheating
(see `gpu_monitor.log`), an unstable overclock or a driver problem.

**I want to start over from scratch** — delete the generated files
(`step4.0_minimization.*`, `step4.1_equilibration.*`, `step5_production.*`,
`.autogromacs_total_ns`), keeping the original CHARMM-GUI inputs, and run the
script again.

## 12. FAQ

**Can I extend a finished run?** Yes. With GROMACS directly:
`gmx convert-tpr -s dynamic200ns.tpr -extend 100000 -o dynamic300ns.tpr`,
then `gmx mdrun -s dynamic300ns.tpr -cpi dynamic200ns.cpt -deffnm dynamic200ns -append`.
Rename the files afterwards if you want the `dynamic300ns` name.

**Can I run several systems at once?** Each system runs in its own
`gromacs/` folder, so there is no file conflict. On a single GPU, though,
simultaneous runs share the card and each one will be slower.

**Why 100 ns segments?** It is a balance: segments are short enough that the
GPU-failure retry repeats little work, and few enough that the overhead of
extending the `.tpr` is negligible.

**Does it change the `.mdp` files?** No. All length changes are made on the
`.tpr` with `convert-tpr`; the CHARMM-GUI `.mdp` files are left as they are.

## 13. Limitations

- Only CHARMM-GUI **Solution Builder** naming (`step4.0`, `step4.1`,
  `step5`). Membrane Builder systems (`step6.x` / `step7`) are not yet
  supported.
- An NVIDIA GPU is required; there is no CPU-only mode.
- Single GPU per run (`-ntmpi 1`).

---

**Authors:** Ricardo Aguilar-Garay (ricardoaguilargaray@gmail.com) ·
contributor: Dulce Estefanía Nicolás Álvarez, Departamento de Fisiología,
ENCB-IPN (dnicolas@ipn.mx). MIT License. Please cite the reference article
listed in the [main README](../README.md#citing).
