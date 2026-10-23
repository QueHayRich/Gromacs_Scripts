#!/bin/bash
###############################################################################
# free_energy.sh — binding energetics for Autogromacs trajectories
#
# Author : Ricardo Aguilar-Garay
# Contact: ricardoaguilargaray@gmail.com
#          (questions, bug reports, and suggestions are welcome)
#
# Contributor: Dulce Estefanía Nicolás Álvarez
#              Departamento de Fisiología, Escuela Nacional de Ciencias
#              Biológicas, Instituto Politécnico Nacional (ENCB-IPN), Mexico
#              dnicolas@ipn.mx
#
# Companion of Autogromacs.sh and basic_analyses_with_names.sh. Copy this SINGLE file into
# any CHARMM-GUI gromacs/ folder (after production) and run:
#
#       bash free_energy.sh                 MM-GBSA + LJ + validation (default)
#       bash free_energy.sh --check         only show what was detected; computes nothing
#       bash free_energy.sh mmgbsa          only MM-GBSA (+ validation)
#       bash free_energy.sh lj              only LJ/Coulomb interaction (+ validation)
#       bash free_energy.sh summary         rebuild summary.csv from existing results
#
# Options (all optional; defaults are the recommended protocol):
#       --start NS      first ns of the analysis window   (default: 25% of the run)
#       --end NS        last ns of the analysis window    (default: end of the run)
#       --frames N      approximate number of frames      (default: 300)
#       --ligand NAME   ligand molecule(s) from topol.top, comma-separated
#                       (only needed if auto-detection is ambiguous)
#       --receptor NAME receptor molecule(s), comma-separated (default: all proteins)
#       --np N          MPI processes for gmx_MMPBSA      (default: min(cores, 16))
#       --env NAME      conda env that has gmx_MMPBSA. If not given, the script
#                       finds the envs that contain gmx_MMPBSA and ASKS you the
#                       first time (Enter = suggested); the answer is remembered
#                       in ~/.config/free_energy/conda_env. --env ask = choose again
#       --no-decomp     skip per-residue decomposition    (faster)
#       --base NAME     trajectory basename, e.g. dynamic200ns (default: auto)
#       --name NAME     system name for tables/Excel (default: first folder above
#                       gromacs/ that is not a generic name like complex, charmm-gui-*)
#
# NOTHING has to be chosen by hand: receptor and ligand are detected from the
# [ molecules ] section of topol.top (proteins = receptor; water and ions are
# ignored; whatever remains is the ligand). If the system is ambiguous (two
# different ligands, a protein-protein complex, a membrane...), the script
# STOPS before computing anything and tells you which option to add.
#
# What is computed (all on the SAME window of frames):
#   1. Window trajectory, PBC-fixed and centered on the receptor
#   2. MM-GBSA (gmx_MMPBSA, igb=5, mbondi2, 0.15 M, T from the .mdp):
#      dVDW, dEEL, dEGB, dESURF, dGgas, dGsolv, dH (mean / SD / SEM),
#      -TdS by the C2 (second-order cumulant) method, dG = dH - TdS,
#      ligand efficiency, per-residue decomposition (top residues)
#   3. Direct interaction energy receptor-ligand (mdrun -rerun, energygrps):
#      E_LJ, E_Coul, E_int (mean / SD / Err.Est.), E_LJ per heavy atom
#   4. Validation: receptor-ligand H-bonds, contacts < 4 A, minimum distance,
#      ligand RMSD in the receptor frame (flags unbound / displaced ligands)
#
# Output:  Analyses/FreeEnergy_<system>_<N>ns/
#              window.xtc  index_FE.ndx  detected_system.txt
#              MMGBSA/  LJ/  validation/
#              summary.csv  top_residues.csv
# Then, from any parent folder:  python3 collect_free_energy.py <folder> ...
#
# RESUMABLE: finished steps are skipped when the script is run again.
###############################################################################
set -uo pipefail
FE_VERSION="1.4 (2026-10-14)"

# ── 0. Arguments ────────────────────────────────────────────────────────────
MODE="both"; START_NS=""; END_NS=""; NFRAMES=300
LIGAND_OPT=""; RECEPTOR_OPT=""; NP=""; CONDA_ENV="${GMXMMPBSA_ENV:-}"
DECOMP=1; BASE=""; CHECK_ONLY=0; NAME_OPT=""

usage(){ awk 'NR>=9 && /^# NOTHING/{exit} NR>=9' "$0" | sed 's/^# \{0,1\}//'; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        both|mmgbsa|lj|summary) MODE="$1" ;;
        --check)     CHECK_ONLY=1 ;;
        --start)     START_NS="${2:-}"; shift ;;
        --end)       END_NS="${2:-}"; shift ;;
        --frames)    NFRAMES="${2:-}"; shift ;;
        --ligand)    LIGAND_OPT="${2:-}"; shift ;;
        --receptor)  RECEPTOR_OPT="${2:-}"; shift ;;
        --np)        NP="${2:-}"; shift ;;
        --env)       CONDA_ENV="${2:-}"; shift ;;
        --no-decomp) DECOMP=0 ;;
        --base)      BASE="${2:-}"; shift ;;
        --name)      NAME_OPT="${2:-}"; shift ;;
        -h|--help)   usage; exit 0 ;;
        -v|--version) echo "free_energy.sh v${FE_VERSION}"; exit 0 ;;
        *) echo "❌ Unknown option: $1"; echo; usage; exit 1 ;;
    esac
    shift
done
is_num(){ [[ "$1" =~ ^[0-9]+([.][0-9]+)?$ ]]; }
for v in "$START_NS" "$END_NS"; do
    [[ -z "$v" ]] || is_num "$v" || { echo "❌ --start/--end must be numbers (ns)."; exit 1; }
done
[[ "$NFRAMES" =~ ^[0-9]+$ ]] && (( NFRAMES >= 10 )) || { echo "❌ --frames must be an integer >= 10."; exit 1; }

PYTHON=$(command -v python3 || command -v python || true)
[[ -n "$PYTHON" ]] || { echo "❌ python3 not found."; exit 1; }

# ── 1. Locate GROMACS (system build, same logic as Autogromacs.sh) ──────────
find_gmx(){
    if command -v gmx >/dev/null 2>&1; then echo "gmx"
    elif command -v gmx_mpi >/dev/null 2>&1; then echo "gmx_mpi"
    else echo ""
    fi
}
GMX=$(find_gmx)
if [[ -z "$GMX" ]]; then
    for GMXRC in /usr/local/gromacs/bin/GMXRC /usr/local/gromacs/*/bin/GMXRC \
                 /opt/gromacs/bin/GMXRC "$HOME/gromacs/bin/GMXRC"; do
        if [[ -f "$GMXRC" ]]; then
            # shellcheck disable=SC1090
            set +u; source "$GMXRC"; set -u     # GMXRC uses unset variables
            GMX=$(find_gmx); [[ -n "$GMX" ]] && break
        fi
    done
fi
[[ -n "$GMX" ]] || { echo "❌ ERROR: neither 'gmx' nor 'gmx_mpi' found."; exit 1; }
GMX_PATH=$(command -v "$GMX"); GMX_DIR=$(dirname "$GMX_PATH")
GMX_VER=$("$GMX" --version 2>/dev/null | grep -m1 -oE '[0-9]{4}(\.[0-9]+)?' | head -n1)
NTMPI_FLAG="-ntmpi 1"; [[ "$GMX" == "gmx_mpi" ]] && NTMPI_FLAG=""
NCORES=$(nproc 2>/dev/null || echo 4)

# ── 2. Trajectory (same detection as basic_analyses_with_names.sh) ─────────────────────
TPR=""; XTC=""
if [[ -n "$BASE" ]]; then
    TPR="${BASE}.tpr"; XTC="${BASE}.xtc"
else
    for f in $(ls -t dynamic*ns.tpr 2>/dev/null); do
        b="${f%.tpr}"
        if [[ -f "${b}.xtc" ]]; then TPR="$f"; XTC="${b}.xtc"; BASE="$b"; break; fi
    done
    if [[ -z "$TPR" ]]; then
        for f in $(ls -t step5_production*.tpr 2>/dev/null); do
            if [[ -f step5_production.xtc ]]; then
                TPR="$f"; XTC="step5_production.xtc"; BASE="${f%.tpr}"; break
            fi
        done
    fi
fi
if [[ -z "$TPR" || ! -f "$TPR" || ! -f "$XTC" ]]; then
    echo "❌ ERROR: no trajectory found (dynamic<N>ns.{tpr,xtc} or step5_production.*)."
    echo "   Run this script inside the CHARMM-GUI gromacs/ folder after Autogromacs.sh."
    exit 1
fi
for f in topol.top index.ndx; do
    [[ -f "$f" ]] || { echo "❌ ERROR: '$f' not found. Are you inside a CHARMM-GUI gromacs/ folder?"; exit 1; }
done

TPR_DUMP=$("$GMX" dump -s "$TPR" 2>/dev/null)
dump_val(){ grep -m1 -E "^[[:space:]]*$1[[:space:]]*[=:]" <<<"$TPR_DUMP" | sed -E 's/.*[=:][[:space:]]*//' | awk '{print $1}'; }
DT=$(dump_val "delta[-_]t"); DT=${DT:-0.002}
NSTEPS=$(dump_val "nsteps")
NSTXOUTC=$(dump_val "nstxout[-_]compressed")
REF_T=$(grep -m1 -E "^[[:space:]]*ref[-_]t:" <<<"$TPR_DUMP" | awk '{print $2}')
REF_T=${REF_T:-303.15}

if [[ "$BASE" =~ dynamic([0-9]+)ns ]]; then
    TOTAL_NS="${BASH_REMATCH[1]}"
else
    TOTAL_NS=$(awk -v s="${NSTEPS:-0}" -v dt="$DT" 'BEGIN{printf "%.0f", s*dt/1000}')
fi
(( TOTAL_NS > 0 )) || { echo "❌ ERROR: could not determine the simulation length."; exit 1; }

# Frame spacing of the xtc (ps); fallback 100 ps (CHARMM-GUI default)
FRAME_PS=$(awk -v n="${NSTXOUTC:-0}" -v dt="$DT" 'BEGIN{v=n*dt; if(v<=0)v=100; printf "%g", v}')

# Analysis window
[[ -z "$START_NS" ]] && START_NS=$(awk -v t="$TOTAL_NS" 'BEGIN{printf "%g", t*0.25}')
[[ -z "$END_NS"   ]] && END_NS="$TOTAL_NS"
awk -v s="$START_NS" -v e="$END_NS" -v t="$TOTAL_NS" 'BEGIN{exit !(s>=0 && e>s && e<=t+1e-6)}' \
    || { echo "❌ Invalid window: ${START_NS}-${END_NS} ns (run length ${TOTAL_NS} ns)."; exit 1; }
WIN_DT=$(awk -v s="$START_NS" -v e="$END_NS" -v n="$NFRAMES" -v f="$FRAME_PS" \
    'BEGIN{d=(e-s)*1000/n; k=int(d/f+0.5); if(k<1)k=1; printf "%g", k*f}')
B_PS=$(awk -v s="$START_NS" 'BEGIN{printf "%g", s*1000}')
E_PS=$(awk -v e="$END_NS"   'BEGIN{printf "%g", e*1000}')

# System name: first folder above gromacs/ that is not a generic name
# (gromacs, charmm-gui-*, complex, system, md, run, ...). Override with --name.
GENERIC='^(gromacs|charmm-gui-[0-9]+|charmm-gui|openmm|namd|amber|complex|complejo|complexes|system|sistema|systems|md|mds|run|runs|sim|sims|simulation|simulations|simulacion|production|prod|dinamica|dynamics)$'
sysname(){
    local d="$PWD" n nl
    while [[ "$d" != "/" && "$d" != "$HOME" ]]; do
        n="${d##*/}"; nl=$(echo "$n" | tr '[:upper:]' '[:lower:]')
        if [[ ! "$nl" =~ $GENERIC ]]; then echo "$n"; return; fi
        d=$(dirname "$d")
    done
    echo "${PWD##*/}"
}
SYSTEM="${NAME_OPT:-$(sysname)}"
SYS_TAG=$(echo "$SYSTEM" | tr -c 'A-Za-z0-9._,+-' '_' | sed 's/_*$//')
OUTDIR="Analyses/FreeEnergy_${SYS_TAG}_${TOTAL_NS}ns"
# Adopt results from an older run of this script that used another name
# (e.g. FreeEnergy_complex_200ns): if there is exactly one, rename it.
if [[ ! -d "$OUTDIR" ]]; then
    shopt -s nullglob; OLD=(Analyses/FreeEnergy_*_"${TOTAL_NS}"ns); shopt -u nullglob
    if (( ${#OLD[@]} == 1 )); then
        echo "ℹ Previous results found in ${OLD[0]} → renamed to ${OUTDIR} (nothing is recomputed)"
        mv "${OLD[0]}" "$OUTDIR"
    elif (( ${#OLD[@]} > 1 )); then
        echo "⚠ Several previous result folders found (${OLD[*]}); a new one is created: ${OUTDIR}"
    fi
fi
mkdir -p "$OUTDIR"/{MMGBSA,LJ,validation}
OUTABS=$(cd "$OUTDIR" && pwd)
NDX="${OUTDIR}/index_FE.ndx"

# ── 3. Detect receptor / ligand from topol.top and build index_FE.ndx ───────
# The helper prints KEY='value' lines (sourced below) and writes the index:
# CHARMM-GUI groups (SOLU/SOLV/SYSTEM, needed by grompp) + Receptor, Ligand,
# Receptor_BB, Ligand_heavy, Complex.
DETECT=$("$PYTHON" - "$LIGAND_OPT" "$RECEPTOR_OPT" "$NDX" <<'PYDETECT'
import os, re, sys
lig_opt = [x for x in sys.argv[1].split(",") if x.strip()]
rec_opt = [x for x in sys.argv[2].split(",") if x.strip()]
ndx_out = sys.argv[3]

WATER = {"TIP3","SOL","TIP4","TIP4P","TIP5P","SPC","SPCE","WAT","HOH","TP3","T3P","SWM4","TIP3P"}
IONS  = {"POT","CLA","SOD","MG","CAL","ZN2","ZN","CES","LIT","RUB","BAR","CD2","NA","CL","K",
         "CA","LI","CS","RB","FE2","FE3","MN","CU","CU2","NI","CO","IOD","BR","F","FLU"}
AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","HSD","HSE","HSP","HID","HIE","HIP",
      "ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","CYX","CYM","ASH","GLH",
      "LYN","ARN","SEP","TPO","PTR","NLE","MSE","ACE","NME","NMA","CT3"}

def q(s): return "'" + str(s).replace("'", "'\\''") + "'"
def fail(msg):
    print("FAIL=" + q(msg)); sys.exit(0)

moltypes, molecules, atomtype_mass = {}, [], {}
def parse(path, top_level=False, seen=None):
    seen = seen if seen is not None else set()
    path = os.path.normpath(path)
    if path in seen or not os.path.isfile(path): return
    seen.add(path)
    base, section, cur = os.path.dirname(path), None, None
    for raw in open(path, errors="replace"):
        line = raw.split(";")[0].strip()
        if not line: continue
        if line.startswith("#include"):
            m = re.search(r'["<](.+?)[">]', line)
            if m:
                inc = m.group(1)
                cand = os.path.join(base, inc)
                parse(cand if os.path.isfile(cand) else inc, False, seen)
            continue
        if line.startswith("#"): continue
        m = re.match(r"\[\s*(\w+)\s*\]", line)
        if m: section = m.group(1).lower(); continue
        f = line.split()
        if section == "atomtypes" and len(f) >= 4:
            for tok in f[1:4]:
                try: atomtype_mass[f[0]] = float(tok); break
                except ValueError: pass
        elif section == "moleculetype":
            cur = f[0]; moltypes[cur] = []; section = "moleculetype_done"
        elif section == "atoms" and cur is not None and len(f) >= 5:
            try: charge = float(f[6]) if len(f) > 6 else 0.0
            except ValueError: charge = 0.0
            try: mass = float(f[7]) if len(f) > 7 else None
            except ValueError: mass = None
            if mass is None: mass = atomtype_mass.get(f[1], 1.008 if f[4].upper().startswith("H") else 12.0)
            moltypes[cur].append((f[3], f[4], charge, mass))
        elif section == "molecules" and len(f) >= 2:
            try: molecules.append((f[0], int(f[1])))
            except ValueError: pass

parse("topol.top", True)
if not molecules: fail("could not read the [ molecules ] section of topol.top")
missing = [n for n, _ in molecules if n not in moltypes]
if missing: fail("moleculetype(s) not found in the included .itp files: " + ", ".join(missing))

def is_heavy(a): return a[3] >= 4.0          # robust to hydrogen mass repartitioning
def kind(name):
    atoms = moltypes[name]; resn = [a[0].upper() for a in atoms]
    if name.upper() in WATER or set(resn) <= WATER: return "water"
    if name.upper() in IONS or (len(atoms) == 1 and abs(sum(a[2] for a in atoms)) > 0.5): return "ion"
    if re.match(r"^PRO[A-Z]$", name) or (resn and sum(r in AA for r in resn) / len(resn) > 0.5): return "protein"
    return "other"

# Global atom ranges in topology order (1-based)
blocks, start = [], 1
for name, count in molecules:
    n = len(moltypes[name])
    blocks.append(dict(name=name, count=count, natoms=n, start=start, kind=kind(name)))
    start += n * count

def match(blk, names):
    up = [x.upper() for x in names]
    return blk["name"].upper() in up or any(a[0].upper() in up for a in moltypes[blk["name"]])

notes = []
if lig_opt:
    lig = [b for b in blocks if match(b, lig_opt)]
    if not lig: fail("--ligand %s does not match any molecule in topol.top. Molecules: %s"
                     % (",".join(lig_opt), ", ".join(b["name"] for b in blocks)))
else:
    cands = [b for b in blocks if b["kind"] == "other"]
    bulk = [b for b in cands if b["count"] > 5]
    lig = [b for b in cands if b["count"] <= 5]
    if bulk:
        notes.append("ignored as bulk/membrane/cosolvent (>5 copies): " +
                     ", ".join("%s(x%d)" % (b["name"], b["count"]) for b in bulk))
    if not lig:
        prots = [b for b in blocks if b["kind"] == "protein"]
        if len(prots) > 1:
            fail("no small-molecule ligand found but %d protein chains (%s). For a protein-protein "
                 "complex re-run with --ligand <chain>, e.g. --ligand %s"
                 % (len(prots), ", ".join(b["name"] for b in prots), prots[-1]["name"]))
        fail("no ligand found: topol.top only contains protein, water and ions ("
             + ", ".join(b["name"] for b in blocks) + ")")
    if len({b["name"] for b in lig}) > 1:
        names = ", ".join("%s(x%d, %d atoms)" % (b["name"], b["count"], b["natoms"]) for b in lig)
        fail("more than one candidate ligand: %s. Re-run choosing one (or several), e.g. "
             "--ligand %s   or   --ligand %s" % (names, lig[0]["name"], ",".join(b["name"] for b in lig)))

lig_ids = {id(b) for b in lig}
if rec_opt:
    rec = [b for b in blocks if match(b, rec_opt) and id(b) not in lig_ids]
    if not rec: fail("--receptor %s does not match any molecule (or it is the ligand)." % ",".join(rec_opt))
else:
    rec = [b for b in blocks if b["kind"] == "protein" and id(b) not in lig_ids]
    if not rec: fail("no protein found to use as receptor. Use --receptor NAME.")
for b in lig:
    if b["count"] > 1: notes.append("ligand %s has %d copies: all are included in the Ligand group" % (b["name"], b["count"]))

def indices(blks, pick=lambda a: True):
    out = []
    for b in blks:
        atoms = moltypes[b["name"]]
        for c in range(b["count"]):
            off = b["start"] + c * b["natoms"]
            out += [off + i for i, a in enumerate(atoms) if pick(a)]
    return out

rec_idx = indices(rec); lig_idx = indices(lig)
bb_idx  = indices(rec, lambda a: a[1].upper() in ("N", "CA", "C"))
if not bb_idx: bb_idx = indices(rec, is_heavy)
lh_idx  = indices(lig, is_heavy)
cpx_idx = sorted(rec_idx + lig_idx)

orig = open("index.ndx").read()
ngroups = len(re.findall(r"^\s*\[", orig, re.M))
def group(name, idx):
    s = "[ %s ]\n" % name
    for i in range(0, len(idx), 15): s += " ".join("%6d" % x for x in idx[i:i+15]) + "\n"
    return s
with open(ndx_out, "w") as fh:
    fh.write(orig.rstrip("\n") + "\n")
    for name, idx in (("Receptor", rec_idx), ("Ligand", lig_idx), ("Receptor_BB", bb_idx),
                      ("Ligand_heavy", lh_idx), ("Complex", cpx_idx)):
        fh.write(group(name, idx))

charge = lambda blks: sum(sum(a[2] for a in moltypes[b["name"]]) * b["count"] for b in blks)
lig_resn = sorted({a[0] for b in lig for a in moltypes[b["name"]]})
total_q = charge(blocks)
print("REC_NAMES=" + q(",".join(b["name"] + ("(x%d)" % b["count"] if b["count"] > 1 else "") for b in rec)))
print("LIG_NAMES=" + q(",".join(b["name"] + ("(x%d)" % b["count"] if b["count"] > 1 else "") for b in lig)))
print("LIG_RESN=" + q(",".join(lig_resn)))
print("REC_ATOMS=%d" % len(rec_idx)); print("LIG_ATOMS=%d" % len(lig_idx))
print("LIG_HEAVY=%d" % len(lh_idx))
z = lambda x: round(x, 2) + 0.0
print("REC_Q=%.2f" % z(charge(rec)));    print("LIG_Q=%.2f" % z(charge(lig))); print("SYS_Q=%.2f" % z(total_q))
print("G_REC=%d" % ngroups); print("G_LIG=%d" % (ngroups + 1))
print("G_BB=%d" % (ngroups + 2)); print("G_LH=%d" % (ngroups + 3)); print("G_CPX=%d" % (ngroups + 4))
print("ALL_MOLS=" + q(", ".join("%s x%d [%s]" % (b["name"], b["count"], b["kind"]) for b in blocks)))
print("NOTES=" + q(" | ".join(notes)))
PYDETECT
)
if [[ "$DETECT" != *"LIG_NAMES="* && "$DETECT" != *"FAIL="* ]]; then
    echo "❌ ERROR: could not parse topol.top / the .itp files (python helper crashed)."
    exit 1
fi
eval "$DETECT"
if [[ -n "${FAIL:-}" ]]; then
    echo "❌ AUTO-DETECTION STOPPED: ${FAIL}"
    echo "   Nothing was computed. Fix the option and run again (add --check to only test detection)."
    exit 1
fi

# ── 3b. Conda environment with gmx_MMPBSA (asked once, then remembered) ─────
# GROMACS (gmx / gmx_mpi) is the system build used for trjconv/rerun; gmx_MMPBSA
# lives in a conda env. Non-interactive shells (nohup) usually do NOT have conda
# in PATH, so conda is also searched in the usual install locations.
ENV_SAVE="${XDG_CONFIG_HOME:-$HOME/.config}/free_energy/conda_env"
CONDA_BASE=""; ENV_PATH=""; ENV_SRC=""
find_conda_base(){
    local c
    if command -v conda >/dev/null 2>&1; then conda info --base 2>/dev/null && return; fi
    for c in "${CONDA_EXE:-}" "$HOME"/miniconda3/bin/conda "$HOME"/anaconda3/bin/conda \
             "$HOME"/miniforge3/bin/conda "$HOME"/mambaforge/bin/conda /opt/conda/bin/conda \
             /opt/miniconda3/bin/conda /opt/anaconda3/bin/conda; do
        [[ -n "$c" && -x "$c" ]] && { dirname "$(dirname "$c")"; return; }
    done
}
env_names(){ local e o=""; for e in "${ENVS[@]}"; do o+="${e%%|*} "; done; echo "${o:-none}"; }
mmpbsa_envs(){                     # prints "name|path" of envs that contain gmx_MMPBSA
    local d
    for d in "$CONDA_BASE" "$CONDA_BASE"/envs/*; do
        [[ -x "$d/bin/gmx_MMPBSA" ]] || continue
        if [[ "$d" == "$CONDA_BASE" ]]; then echo "base|$d"; else echo "${d##*/}|$d"; fi
    done
}
env_path(){ [[ "$1" == "base" ]] && echo "$CONDA_BASE" || echo "$CONDA_BASE/envs/$1"; }

if [[ "$MODE" == "both" || "$MODE" == "mmgbsa" ]]; then
    CONDA_BASE=$(find_conda_base | tail -n1)
    if [[ -z "$CONDA_BASE" || ! -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
        echo "❌ conda was not found (PATH, ~/miniconda3, ~/anaconda3, ~/miniforge3, /opt/conda)."
        echo "   MM-GBSA needs the conda env with gmx_MMPBSA. Options:"
        echo "     - run 'conda init bash', open a new terminal and try again"
        echo "     - or compute only the LJ part:  bash free_energy.sh lj"
        exit 1
    fi
    mapfile -t ENVS < <(mmpbsa_envs)
    [[ "$CONDA_ENV" == "ask" ]] && { CONDA_ENV=""; rm -f "$ENV_SAVE"; }
    if [[ -n "$CONDA_ENV" ]]; then
        ENV_SRC="--env"
    elif [[ -f "$ENV_SAVE" ]]; then
        CONDA_ENV=$(head -n1 "$ENV_SAVE"); ENV_SRC="saved in $ENV_SAVE"
    fi
    if [[ -z "$CONDA_ENV" ]]; then
        SUGG=""
        for e in "${ENVS[@]}"; do [[ "${e%%|*}" == "gmxMMPBSA" ]] && SUGG="gmxMMPBSA"; done
        [[ -z "$SUGG" && -n "${CONDA_DEFAULT_ENV:-}" ]] && \
            for e in "${ENVS[@]}"; do [[ "${e%%|*}" == "$CONDA_DEFAULT_ENV" ]] && SUGG="$CONDA_DEFAULT_ENV"; done
        [[ -z "$SUGG" && ${#ENVS[@]} -gt 0 ]] && SUGG="${ENVS[0]%%|*}"
        if [[ -t 0 ]]; then
            echo "───────────────────────────────────────────────────────────────────────"
            echo " Which conda environment has gmx_MMPBSA?   (conda: ${CONDA_BASE})"
            if (( ${#ENVS[@]} )); then
                i=1; for e in "${ENVS[@]}"; do printf "   %d) %-20s %s\n" "$i" "${e%%|*}" "${e#*|}"; i=$((i+1)); done
            else
                echo "   (no env with gmx_MMPBSA was found automatically)"
            fi
            read -rp " Number or env name [Enter = ${SUGG:-none}]: " ans
            ans="${ans// /}"
            if [[ -z "$ans" ]]; then CONDA_ENV="$SUGG"
            elif [[ "$ans" =~ ^[0-9]+$ ]] && (( ans >= 1 && ans <= ${#ENVS[@]} )); then
                CONDA_ENV="${ENVS[$((ans-1))]%%|*}"
            else CONDA_ENV="$ans"; fi
            ENV_SRC="chosen now"
        elif [[ ${#ENVS[@]} -eq 1 ]]; then
            CONDA_ENV="$SUGG"; ENV_SRC="only env with gmx_MMPBSA"
        elif [[ -n "$SUGG" && "$SUGG" == "gmxMMPBSA" ]]; then
            CONDA_ENV="$SUGG"; ENV_SRC="auto (gmxMMPBSA)"
        else
            echo "❌ Cannot ask for the conda env (no terminal, e.g. nohup) and it is ambiguous."
            echo "   Envs with gmx_MMPBSA: $(env_names)"
            echo "   Run once in a terminal (it will ask and remember) or pass --env NAME."
            exit 1
        fi
    fi
    [[ -n "$CONDA_ENV" ]] || { echo "❌ No conda env chosen. Use --env NAME."; exit 1; }
    ENV_PATH=$(env_path "$CONDA_ENV")
    if [[ ! -x "$ENV_PATH/bin/gmx_MMPBSA" ]]; then
        echo "❌ conda env '${CONDA_ENV}' does not contain gmx_MMPBSA (${ENV_PATH}/bin/gmx_MMPBSA)."
        echo "   Envs with gmx_MMPBSA: $(env_names)"
        echo "   Choose again with:  bash free_energy.sh --env ask"
        rm -f "$ENV_SAVE"
        exit 1
    fi
    if [[ "$ENV_SRC" == "chosen now" || "$ENV_SRC" == "--env" || ! -f "$ENV_SAVE" ]]; then
        mkdir -p "$(dirname "$ENV_SAVE")" && echo "$CONDA_ENV" > "$ENV_SAVE"
    fi
fi

# ── 4. Pre-flight report ────────────────────────────────────────────────────
NFR_EST=$(awk -v s="$B_PS" -v e="$E_PS" -v d="$WIN_DT" 'BEGIN{printf "%d", (e-s)/d+1}')
REPORT=$(cat <<EOF
═══════════════════════════════════════════════════════════════════════
 free_energy v${FE_VERSION} — ${SYSTEM}
───────────────────────────────────────────────────────────────────────
 Trajectory : ${XTC}  (${TOTAL_NS} ns, one frame every ${FRAME_PS} ps)
 Molecules  : ${ALL_MOLS}
 RECEPTOR   : ${REC_NAMES}   ${REC_ATOMS} atoms, charge ${REC_Q}
 LIGAND     : ${LIG_NAMES} (resname ${LIG_RESN})   ${LIG_ATOMS} atoms, ${LIG_HEAVY} heavy, charge ${LIG_Q}
 System net charge : ${SYS_Q}
 Window     : ${START_NS}-${END_NS} ns, every ${WIN_DT} ps  (~${NFR_EST} frames)
 Temperature: ${REF_T} K (from the tpr)
 GROMACS    : ${GMX_PATH} (${GMX_VER:-?})
 Mode       : ${MODE}$( ((DECOMP)) || echo "  (no decomposition)")
 gmx_MMPBSA : $( [[ -n "$ENV_PATH" ]] && echo "conda env '${CONDA_ENV}' (${ENV_SRC}; change with --env ask)" || echo "not needed (mode ${MODE})")
 Output     : ${OUTDIR}/
EOF
)
echo "$REPORT"
[[ -n "${NOTES:-}" ]] && echo " ⚠ ${NOTES}"
awk -v q="$SYS_Q" 'BEGIN{exit !(q>0.01||q<-0.01)}' && \
    echo " ⚠ System net charge is not zero: check topol.top (missing ions?)."
awk -v q="$LIG_Q" 'BEGIN{exit !(q>=3||q<=-3)}' && \
    echo " ⚠ Highly charged ligand (${LIG_Q}): dEEL and dEGB will be large and opposite;" \
         "expect a large SD in dH."
echo "═══════════════════════════════════════════════════════════════════════"
echo "$REPORT" > "${OUTDIR}/detected_system.txt"
[[ -n "${NOTES:-}" ]] && echo "NOTES: ${NOTES}" >> "${OUTDIR}/detected_system.txt"

if (( CHECK_ONLY )); then
    echo "✓ --check: detection OK. Index written to ${NDX}. Nothing else was computed."
    exit 0
fi

# ── helpers ─────────────────────────────────────────────────────────────────
WIN_XTC="${OUTDIR}/window.xtc"
DO_MMGBSA=0; DO_LJ=0
[[ "$MODE" == "both" || "$MODE" == "mmgbsa" ]] && DO_MMGBSA=1
[[ "$MODE" == "both" || "$MODE" == "lj"     ]] && DO_LJ=1

# ── 5. Window trajectory (PBC-fixed, centered on the receptor) ──────────────
if [[ "$MODE" != "summary" ]]; then
    if [[ ! -s "$WIN_XTC" ]]; then
        echo "=== [1] Window trajectory ${START_NS}-${END_NS} ns every ${WIN_DT} ps ==="
        printf "Receptor\nSystem\n" | "$GMX" trjconv -s "$TPR" -f "$XTC" -n "$NDX" \
            -o "${WIN_XTC%.xtc}_tmp.xtc" -b "$B_PS" -e "$E_PS" -dt "$WIN_DT" \
            -pbc mol -center >"${OUTDIR}/trjconv.log" 2>&1 \
            || { echo "❌ trjconv failed, see ${OUTDIR}/trjconv.log"; exit 1; }
        mv -f "${WIN_XTC%.xtc}_tmp.xtc" "$WIN_XTC"
    else
        echo "✓ [1] Window trajectory already exists, skipping."
    fi
fi

# ── 6. Validation (cheap; always computed) ──────────────────────────────────
VAL="${OUTDIR}/validation"
if [[ "$MODE" != "summary" ]]; then
    echo "=== [2] Validation: contacts, minimum distance, ligand RMSD, H-bonds ==="
    if [[ ! -s "${VAL}/numcont.xvg" ]]; then
        printf "Receptor\nLigand\n" | "$GMX" mindist -s "$TPR" -f "$WIN_XTC" -n "$NDX" \
            -od "${VAL}/mindist.xvg" -on "${VAL}/numcont.xvg" -d 0.4 -tu ns \
            >"${VAL}/mindist.log" 2>&1 || echo "  ⚠ mindist failed (see ${VAL}/mindist.log)"
    fi
    if [[ ! -s "${VAL}/ligand_rmsd.xvg" ]]; then
        printf "Receptor_BB\nLigand_heavy\n" | "$GMX" rms -s "$TPR" -f "$WIN_XTC" -n "$NDX" \
            -o "${VAL}/ligand_rmsd.xvg" -tu ns >"${VAL}/rms.log" 2>&1 \
            || echo "  ⚠ rms failed (see ${VAL}/rms.log)"
    fi
    if [[ ! -s "${VAL}/hbonds.xvg" ]]; then
        if ! "$GMX" hbond -s "$TPR" -f "$WIN_XTC" -n "$NDX" -r Receptor -t Ligand \
                -num "${VAL}/hbonds.xvg" >"${VAL}/hbond.log" 2>&1; then
            printf "Receptor\nLigand\n" | "$GMX" hbond -s "$TPR" -f "$WIN_XTC" -n "$NDX" \
                -num "${VAL}/hbonds.xvg" >"${VAL}/hbond.log" 2>&1 \
                || echo "  ⚠ hbond failed (see ${VAL}/hbond.log)"
        fi
    fi
    rm -f "${VAL}"/\#*\# 2>/dev/null
fi

# ── 7. LJ / Coulomb interaction energy (mdrun -rerun with energygrps) ───────
# energygrps are not supported on GPU, so the rerun runs on CPU (fast: only
# ~NFRAMES force evaluations). 1-4 terms are always zero between different
# molecules and are not reported.
LJD="${OUTDIR}/LJ"
if (( DO_LJ )); then
    if [[ ! -s "${LJD}/lj_coul.xvg" ]]; then
        echo "=== [3] LJ / Coulomb receptor-ligand interaction (rerun on CPU) ==="
        PROD_MDP=$(ls step5_production.mdp 2>/dev/null || ls step5*.mdp 2>/dev/null | head -n1)
        [[ -n "$PROD_MDP" ]] || { echo "❌ No step5 production .mdp found for the rerun."; exit 1; }
        awk 'BEGIN{IGNORECASE=1}
             { line=$0; key=$0; sub(/;.*/,"",key); split(key,kv,"="); k=kv[1];
               gsub(/[[:space:]]/,"",k); k=tolower(k); gsub(/_/,"-",k)
               if (k ~ /^(energygrps|nstcalcenergy|nstenergy|nstlog|nstxout|nstvout|nstfout|nstxout-compressed|continuation|gen-vel|define)$/) next
               print line }
             END{ print "; --- free_energy.sh rerun overrides ---"
                  print "energygrps          = Receptor Ligand"
                  print "nstcalcenergy       = 1"
                  print "nstenergy           = 1"
                  print "nstlog              = 0"
                  print "nstxout             = 0"
                  print "nstvout             = 0"
                  print "nstfout             = 0"
                  print "nstxout-compressed  = 0"
                  print "continuation        = yes"
                  print "gen-vel             = no" }' "$PROD_MDP" > "${LJD}/rerun.mdp"
        GRO=$(ls step3_input.gro 2>/dev/null || ls "${BASE}.gro" 2>/dev/null)
        "$GMX" grompp -f "${LJD}/rerun.mdp" -c "$GRO" -p topol.top -n "$NDX" \
            -o "${LJD}/rerun.tpr" -po "${LJD}/rerun_out.mdp" -maxwarn 3 >"${LJD}/grompp.log" 2>&1 \
            || { echo "❌ grompp for the rerun failed, see ${LJD}/grompp.log"; exit 1; }
        "$GMX" mdrun -rerun "$WIN_XTC" -s "${LJD}/rerun.tpr" -deffnm "${LJD}/rerun" \
            -nb cpu -pme cpu ${NTMPI_FLAG} -ntomp "$NCORES" >"${LJD}/mdrun.log" 2>&1 \
            || { echo "❌ mdrun -rerun failed, see ${LJD}/mdrun.log"; exit 1; }
        printf "LJ-SR:Receptor-Ligand\nCoul-SR:Receptor-Ligand\n\n" | \
            "$GMX" energy -f "${LJD}/rerun.edr" -o "${LJD}/lj_coul_tmp.xvg" \
            > "${LJD}/energy_stdout.txt" 2>&1 \
            || { echo "❌ gmx energy failed, see ${LJD}/energy_stdout.txt"; exit 1; }
        mv -f "${LJD}/lj_coul_tmp.xvg" "${LJD}/lj_coul.xvg"
        rm -f "${LJD}/rerun.trr" "${LJD}"/\#*\# 2>/dev/null
    else
        echo "✓ [3] LJ/Coulomb already computed, skipping."
    fi
fi

# ── 8. MM-GBSA with gmx_MMPBSA (inside the conda env) ───────────────────────
MMD="${OUTDIR}/MMGBSA"
if (( DO_MMGBSA )); then
    if [[ ! -s "${MMD}/FINAL_RESULTS_MMPBSA.dat" ]]; then
        echo "=== [4] MM-GBSA (gmx_MMPBSA, conda env '${CONDA_ENV}') ==="
        DECOMP_BLOCK=""
        (( DECOMP )) && DECOMP_BLOCK=$'&decomp\n  idecomp=2, dec_verbose=0, print_res="within 6",\n/'
        cat > "${MMD}/mmgbsa.in" <<EOF
# Generated by free_energy.sh — window ${START_NS}-${END_NS} ns every ${WIN_DT} ps
&general
  sys_name="${SYS_TAG}",
  startframe=1, endframe=9999999, interval=1,
  temperature=${REF_T},
  PBRadii=3,
  c2_entropy=1,
  verbose=2, keep_files=2,
/
&gb
  igb=5, saltcon=0.150,
/
${DECOMP_BLOCK}
EOF
        [[ -z "$NP" ]] && NP=$(( NCORES < 16 ? NCORES : 16 ))
        (
            set +u
            # shellcheck disable=SC1091
            source "${CONDA_BASE}/etc/profile.d/conda.sh" && conda activate "$CONDA_ENV" \
                || { echo "❌ could not activate conda env '$CONDA_ENV' (use --env NAME)."; exit 1; }
            command -v gmx_MMPBSA >/dev/null 2>&1 || { echo "❌ gmx_MMPBSA not found in '$CONDA_ENV'."; exit 1; }
            # gmx_MMPBSA reads the tpr with the gmx found in PATH. If the env's gmx is
            # older than the GROMACS used for the run, fall back to the system build.
            # This is expected when the env has an older GROMACS: it is handled
            # automatically and is NOT an error.
            ENV_GMX=$(command -v gmx || command -v gmx_mpi || true)
            if [[ -n "$ENV_GMX" && "$(dirname "$ENV_GMX")" == "$GMX_DIR" ]]; then
                :   # already the system build
            elif [[ -z "$ENV_GMX" ]]; then
                export PATH="${GMX_DIR}:$PATH"
                echo "  ✓ gmx_MMPBSA uses the system GROMACS ${GMX_VER:-} (the env has no gmx)"
            elif ! "$ENV_GMX" dump -s "$TPR" >/dev/null 2>&1; then
                ENV_GMX_VER=$("$ENV_GMX" --version 2>/dev/null | grep -m1 -oE '[0-9]{4}(\.[0-9]+)?' | head -n1)
                export PATH="${GMX_DIR}:$PATH"
                echo "  ✓ gmx_MMPBSA uses the system GROMACS ${GMX_VER:-} (the env has GROMACS ${ENV_GMX_VER:-?},"
                echo "    too old to read ${TPR}). Expected and handled automatically: nothing to fix."
            fi
            GMXMMPBSA_VER=$(gmx_MMPBSA --version 2>&1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -n1)
            echo "${GMXMMPBSA_VER:-unknown}" > "${OUTABS}/MMGBSA/gmx_MMPBSA_version.txt"
            WORKDIR="$PWD"
            cd "${OUTABS}/MMGBSA" || exit 1
            ARGS=(-O -i mmgbsa.in -cs "${WORKDIR}/${TPR}" -ci "${OUTABS}/index_FE.ndx"
                  -cg "$G_REC" "$G_LIG" -ct "${OUTABS}/window.xtc" -cp "${WORKDIR}/topol.top"
                  -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv -nogui)
            (( DECOMP )) && ARGS+=(-do FINAL_DECOMP_MMPBSA.dat -deo FINAL_DECOMP_MMPBSA.csv)
            ok=1
            if (( NP > 1 )) && command -v mpirun >/dev/null 2>&1; then
                echo "  running with mpirun -np ${NP} ..."
                mpirun -np "$NP" gmx_MMPBSA "${ARGS[@]}" > gmx_MMPBSA.out 2>&1 && ok=0
                if (( ok )); then
                    echo "  ⚠ MPI run failed; retrying in serial (slower) ..."
                    rm -f _GMXMMPBSA_* 2>/dev/null
                fi
            fi
            if (( ok )); then
                gmx_MMPBSA "${ARGS[@]}" > gmx_MMPBSA.out 2>&1 && ok=0
            fi
            (( ok )) && { echo "❌ gmx_MMPBSA failed, see ${OUTABS}/MMGBSA/gmx_MMPBSA.out (and gmx_MMPBSA.log)"; exit 1; }
            rm -f _GMXMMPBSA_*.mdcrd* _GMXMMPBSA_*.nc* 2>/dev/null
            exit 0
        ) || exit 1
    else
        echo "✓ [4] MM-GBSA already computed, skipping."
    fi
fi

# ── 9. Summary (summary.csv + top_residues.csv) ─────────────────────────────
echo "=== [5] Summary ==="
export FE_SYSTEM="$SYSTEM" FE_PATH="$PWD" FE_REC="$REC_NAMES" FE_LIG="$LIG_NAMES" \
       FE_LIG_RESN="$LIG_RESN" FE_LIG_HEAVY="$LIG_HEAVY" FE_LIG_ATOMS="$LIG_ATOMS" \
       FE_LIG_Q="$LIG_Q" FE_REC_Q="$REC_Q" FE_SIM_NS="$TOTAL_NS" FE_START="$START_NS" \
       FE_END="$END_NS" FE_STEP="$WIN_DT" FE_T="$REF_T" FE_GMX="${GMX_VER:-}"
"$PYTHON" - "$OUTDIR" <<'PYSUM'
import csv, datetime, math, os, re, statistics as st, sys
out = sys.argv[1]; E = os.environ
KJ = 1 / 4.184                      # kJ/mol -> kcal/mol
KB = 0.0019872041                   # kcal/(mol K)
nan = float("nan")
def f(x):
    try: return float(x)
    except (TypeError, ValueError): return nan
def isnum(x): return isinstance(x, float) and not math.isnan(x)

def xvg(path):
    """Return (legends, rows) from a GROMACS .xvg file."""
    leg, rows = {}, []
    if not os.path.isfile(path): return leg, rows
    for line in open(path, errors="replace"):
        if line.startswith("#"): continue
        if line.startswith("@"):
            m = re.match(r'@\s+s(\d+)\s+legend\s+"(.*)"', line)
            if m: leg[int(m.group(1))] = m.group(2)
            continue
        p = line.split()
        if p:
            try: rows.append([float(x) for x in p])
            except ValueError: pass
    return leg, rows

def stats(v):
    v = [x for x in v if isnum(x)]
    if not v: return nan, nan, nan
    m = st.fmean(v); sd = st.pstdev(v) if len(v) > 1 else 0.0
    return m, sd, sd / math.sqrt(len(v))
def block_err(v, nb=5):
    v = [x for x in v if isnum(x)]
    if len(v) < 2 * nb: return nan
    n = len(v) // nb
    means = [st.fmean(v[i*n:(i+1)*n]) for i in range(nb)]
    return st.stdev(means) / math.sqrt(nb)

row = dict(system=E["FE_SYSTEM"], receptor=E["FE_REC"], ligand=E["FE_LIG"],
           ligand_resname=E["FE_LIG_RESN"], ligand_heavy_atoms=int(E["FE_LIG_HEAVY"]),
           ligand_charge=f(E["FE_LIG_Q"]), receptor_charge=f(E["FE_REC_Q"]),
           sim_ns=f(E["FE_SIM_NS"]), window_start_ns=f(E["FE_START"]), window_end_ns=f(E["FE_END"]),
           frame_step_ps=f(E["FE_STEP"]), temperature_K=f(E["FE_T"]))
flags, notes = [], []
heavy = row["ligand_heavy_atoms"] or nan

# ---------- MM-GBSA ----------
mm = os.path.join(out, "MMGBSA")
terms = ["VDWAALS", "EEL", "EGB", "ESURF", "GGAS", "GSOLV", "TOTAL"]
names = dict(VDWAALS="dVDW", EEL="dEEL", EGB="dEGB", ESURF="dESURF", GGAS="dGgas", GSOLV="dGsolv", TOTAL="dH")
for t in terms:
    for s in ("mean", "sd", "sem"): row["%s_%s" % (names[t], s)] = nan
row.update(mmgbsa_frames=0, sigma_dGgas=nan, mTdS_C2=nan, dG_C2=nan, LE_dH=nan, LE_dG=nan)

def per_frame_delta(path):
    """Delta (complex - receptor - ligand) per-frame table of the -eo CSV."""
    if not os.path.isfile(path): return None
    lines = open(path, errors="replace").read().splitlines()
    for i, line in enumerate(lines):
        if "delta" in line.lower() and "energy terms" in line.lower():
            j = i + 1
            while j < len(lines) and "frame" not in lines[j].lower(): j += 1
            if j >= len(lines): return None
            head = [h.strip().upper() for h in next(csv.reader([lines[j]]))]
            data = []
            for l in lines[j+1:]:
                cells = next(csv.reader([l])) if l.strip() else []
                if not cells or not re.match(r"^\s*\d", cells[0]): break
                data.append([f(c) for c in cells])
            return head, data
    return None

got = per_frame_delta(os.path.join(mm, "FINAL_RESULTS_MMPBSA.csv"))
if got and got[1]:
    head, data = got
    col = lambda name: [r[head.index(name)] for r in data] if name in head else None
    vals = {t: col(t) for t in terms}
    if vals["GGAS"] is None and vals["VDWAALS"] and vals["EEL"]:
        vals["GGAS"] = [a + b for a, b in zip(vals["VDWAALS"], vals["EEL"])]
    if vals["GSOLV"] is None and vals["EGB"] and vals["ESURF"]:
        vals["GSOLV"] = [a + b for a, b in zip(vals["EGB"], vals["ESURF"])]
    for t in terms:
        if vals[t]:
            m, sd, sem = stats(vals[t])
            row[names[t] + "_mean"], row[names[t] + "_sd"], row[names[t] + "_sem"] = m, sd, sem
    row["mmgbsa_frames"] = len(data)
    if vals["GGAS"]:
        sig = st.pstdev(vals["GGAS"]); row["sigma_dGgas"] = sig
        row["mTdS_C2"] = sig ** 2 / (2 * KB * row["temperature_K"])
else:
    # Fallback: parse the text report (old and new gmx_MMPBSA layouts)
    dat = os.path.join(mm, "FINAL_RESULTS_MMPBSA.dat")
    if os.path.isfile(dat):
        txt = open(dat, errors="replace").read()
        k = max(txt.find("Delta (Complex - Receptor - Ligand)"), txt.find("Differences (Complex - Receptor - Ligand)"))
        if k >= 0:
            for line in txt[k:].splitlines()[1:60]:
                m = re.match(r"^\s*(?:Δ|DELTA\s+)?(G\s*gas|G\s*solv|TOTAL|[A-Z0-9\-]+)\s+(-?\d.*)$", line.strip())
                if not m: continue
                key = m.group(1).replace(" ", "").upper()
                nums = [f(x) for x in m.group(2).split()]
                if key in names and nums:
                    i = (0, 2, 4) if len(nums) >= 5 else (0, 1, 2)
                    row[names[key] + "_mean"] = nums[i[0]]
                    row[names[key] + "_sd"] = nums[i[1]] if len(nums) > i[1] else nan
                    row[names[key] + "_sem"] = nums[i[2]] if len(nums) > i[2] else nan
            notes.append("MM-GBSA read from .dat (per-frame CSV not found): C2 entropy not available")
if isnum(row["dH_mean"]):
    if isnum(row["mTdS_C2"]): row["dG_C2"] = row["dH_mean"] + row["mTdS_C2"]
    row["LE_dH"] = row["dH_mean"] / heavy
    if isnum(row["dG_C2"]): row["LE_dG"] = row["dG_C2"] / heavy
    if isnum(row["sigma_dGgas"]) and row["sigma_dGgas"] > 3.6:
        flags.append("C2_UNRELIABLE")
        notes.append("sigma(dGgas)=%.1f kcal/mol > 3.6: C2 entropy is not reliable; report dH" % row["sigma_dGgas"])
    if isnum(row["dH_sem"]) and row["dH_sem"] > 1.0:
        flags.append("HIGH_SEM"); notes.append("SEM(dH) > 1 kcal/mol")
    if row["dH_mean"] > 0:
        flags.append("POSITIVE_dH"); notes.append("dH > 0: unfavourable or unbound complex")
vf = os.path.join(mm, "gmx_MMPBSA_version.txt")
row["gmx_MMPBSA"] = open(vf).read().strip() if os.path.isfile(vf) else ""

# ---------- per-residue decomposition ----------
top = []
dec = os.path.join(mm, "FINAL_DECOMP_MMPBSA.dat")
lig_resn = {x.upper() for x in row["ligand_resname"].split(",") if x}
def short(label):
    m = re.match(r"^[RL]:[A-Z]?:?([A-Z0-9]{2,4}):(-?\d+)", label.strip())
    if m: return m.group(1) + m.group(2)
    m = re.match(r"^\s*([A-Z0-9]{2,4})\s+(-?\d+)", label)
    return (m.group(1) + m.group(2)) if m else label.strip()
if os.path.isfile(dec):
    lines = open(dec, errors="replace").read().splitlines()
    start = 0
    for i, l in enumerate(lines):
        if l.strip().upper().startswith("DELTA"): start = i; break
    i = start
    while i < len(lines) and not lines[i].strip().lower().startswith("residue"): i += 1
    i += 1
    while i < len(lines) and lines[i].strip().startswith(","): i += 1
    for l in lines[i:]:
        c = next(csv.reader([l])) if l.strip() else []
        if len(c) < 17 or not re.search(r"\d", c[0]): break
        lab = c[0].strip()
        if lab.startswith("L:") or any(tok in lig_resn for tok in re.split(r"[:\s]+", lab.upper())):
            continue
        top.append(dict(residue=short(lab), label=lab, total=f(c[16]), total_sd=f(c[17]) if len(c) > 17 else nan,
                        vdw=f(c[4]), elec=f(c[7]), polar_solv=f(c[10]), nonpolar_solv=f(c[13])))
    top.sort(key=lambda d: d["total"] if isnum(d["total"]) else 1e9)
    top = top[:10]
    if not top: notes.append("per-residue decomposition present but could not be parsed")
row["top5_residues"] = "; ".join("%s(%.2f)" % (d["residue"], d["total"]) for d in top[:5])

# ---------- LJ / Coulomb ----------
for k in ("E_LJ", "E_Coul", "E_int"):
    for s in ("mean", "sd", "err"): row["%s_%s" % (k, s)] = nan
row["E_LJ_per_heavy"] = nan
leg, rows = xvg(os.path.join(out, "LJ", "lj_coul.xvg"))
if rows:
    idx = {}
    for s, name in leg.items():
        if name.startswith("LJ-SR"): idx["E_LJ"] = s + 1
        elif name.startswith("Coul-SR"): idx["E_Coul"] = s + 1
    if not idx and len(rows[0]) >= 3: idx = {"E_Coul": 1, "E_LJ": 2}
    errs = {}
    so = os.path.join(out, "LJ", "energy_stdout.txt")
    if os.path.isfile(so):
        for line in open(so, errors="replace"):
            p = line.split()
            if len(p) >= 3 and p[0].startswith("LJ-SR:"): errs["E_LJ"] = f(p[2]) * KJ
            if len(p) >= 3 and p[0].startswith("Coul-SR:"): errs["E_Coul"] = f(p[2]) * KJ
    series = {}
    for k, c in idx.items():
        v = [r[c] * KJ for r in rows if len(r) > c]; series[k] = v
        m, sd, _ = stats(v)
        row[k + "_mean"], row[k + "_sd"] = m, sd
        e = errs.get(k, nan)
        row[k + "_err"] = e if isnum(e) else block_err(v)
    if "E_LJ" in series and "E_Coul" in series:
        tot = [a + b for a, b in zip(series["E_LJ"], series["E_Coul"])]
        row["E_int_mean"], row["E_int_sd"], _ = stats(tot)
        row["E_int_err"] = math.sqrt(row["E_LJ_err"] ** 2 + row["E_Coul_err"] ** 2)
    row["E_LJ_per_heavy"] = row["E_LJ_mean"] / heavy
    row["lj_frames"] = len(rows)

# ---------- validation ----------
val = os.path.join(out, "validation")
_, hb = xvg(os.path.join(val, "hbonds.xvg"))
row["hbonds_mean"], row["hbonds_sd"], _ = stats([r[1] for r in hb if len(r) > 1])
_, nc = xvg(os.path.join(val, "numcont.xvg"))
row["contacts_4A_mean"], _, _ = stats([r[1] for r in nc if len(r) > 1])
_, md = xvg(os.path.join(val, "mindist.xvg"))
d = [r[1] for r in md if len(r) > 1]
row["mindist_mean_nm"] = stats(d)[0]
row["frac_unbound"] = (sum(x > 0.45 for x in d) / len(d)) if d else nan
_, rm = xvg(os.path.join(val, "ligand_rmsd.xvg"))
r = [x[1] for x in rm if len(x) > 1]
row["lig_rmsd_mean_nm"] = stats(r)[0]; row["lig_rmsd_last_nm"] = r[-1] if r else nan
row["n_frames"] = len(rm) or len(md) or row.get("lj_frames", 0) or row["mmgbsa_frames"]
if isnum(row["frac_unbound"]) and row["frac_unbound"] > 0.10:
    flags.append("LIGAND_UNBOUND")
    notes.append("ligand >4.5 A from the receptor in %.0f%% of frames" % (100 * row["frac_unbound"]))
if isnum(row["lig_rmsd_mean_nm"]) and row["lig_rmsd_mean_nm"] > 1.0:
    flags.append("LIGAND_DISPLACED"); notes.append("ligand RMSD > 1 nm vs. start of production (moved site?)")
if row["n_frames"] and row["n_frames"] < 100:
    flags.append("FEW_FRAMES"); notes.append("fewer than 100 frames in the window")
if abs(row["ligand_charge"]) >= 3:
    notes.append("highly charged ligand: dEEL/dEGB large and opposite")

row["flags"] = ";".join(flags) if flags else "OK"
row["notes"] = " | ".join(notes)
row.update(gromacs=E["FE_GMX"], date=datetime.date.today().isoformat(), path=E["FE_PATH"])

def fmt(v):
    if isinstance(v, float): return "" if math.isnan(v) else "%.4f" % v
    return v
with open(os.path.join(out, "summary.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(row.keys()); w.writerow([fmt(v) for v in row.values()])
with open(os.path.join(out, "top_residues.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["system", "rank", "residue", "label", "total", "total_sd", "vdw", "elec", "polar_solv", "nonpolar_solv"])
    for i, dd in enumerate(top, 1):
        w.writerow([row["system"], i, dd["residue"], dd["label"]] +
                   [fmt(dd[k]) for k in ("total", "total_sd", "vdw", "elec", "polar_solv", "nonpolar_solv")])

def pm(k, e="sem"):
    m, s = row.get(k + "_mean", nan), row.get("%s_%s" % (k, e), nan)
    return "   n/a" if not isnum(m) else "%8.2f ± %.2f" % (m, s if isnum(s) else 0)
print("───────────────────────────────────────────────────────────────────────")
print(" RESULTS (kcal/mol, mean ± SEM; LJ terms mean ± Err.Est.)   frames: %s" % row["n_frames"])
print("   dVDW   %s    dEEL   %s" % (pm("dVDW"), pm("dEEL")))
print("   dEGB   %s    dESURF %s" % (pm("dEGB"), pm("dESURF")))
print("   dH     %s    -TdS(C2) %s   dG(C2) %s" % (pm("dH"),
      "%.2f" % row["mTdS_C2"] if isnum(row["mTdS_C2"]) else "n/a",
      "%.2f" % row["dG_C2"] if isnum(row["dG_C2"]) else "n/a"))
print("   E_LJ   %s    E_Coul %s    E_int %s" % (pm("E_LJ", "err"), pm("E_Coul", "err"), pm("E_int", "err")))
print("   H-bonds %.2f   contacts<4A %.0f   ligand RMSD %.2f nm" % (
      row["hbonds_mean"] if isnum(row["hbonds_mean"]) else nan,
      row["contacts_4A_mean"] if isnum(row["contacts_4A_mean"]) else nan,
      row["lig_rmsd_mean_nm"] if isnum(row["lig_rmsd_mean_nm"]) else nan))
if row["top5_residues"]: print("   Top residues: " + row["top5_residues"])
print("   FLAGS: " + row["flags"] + (("  — " + row["notes"]) if row["notes"] else ""))
PYSUM

echo "═══════════════════════════════════════════════════════════════════════"
echo " DONE. Results in ${OUTDIR}/  (summary.csv, top_residues.csv)"
echo " Collect all systems into one Excel from a parent folder:"
echo "     python3 collect_free_energy.py <parent_folder> [<other_parent> ...]"
if ! "$PYTHON" -c "import openpyxl" >/dev/null 2>&1; then
    echo " ⚠ ALARM: 'openpyxl' is NOT installed for ${PYTHON}: collect_free_energy.py"
    echo "   needs it to write the Excel (otherwise it writes CSV). Install with:"
    echo "       ${PYTHON} -m pip install --user openpyxl      (or: sudo apt install python3-openpyxl)"
    echo "   Check:  python3 collect_free_energy.py --check-deps"
fi
echo "═══════════════════════════════════════════════════════════════════════"
