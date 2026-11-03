#!/bin/bash
###############################################################################
# Autogromacs.sh — automated GROMACS pipeline with ADAPTABLE length (v2)
#                  (Spanish version: Autogromacs.es.sh)
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
# Works for any protein–ligand (or apo protein) system prepared with
# CHARMM-GUI: nothing in the script depends on the protein or the ligand.
#
# Copy this SINGLE file into any CHARMM-GUI gromacs/ folder and run:
#       bash Autogromacs.sh            (asks for the ns the first time)
#       bash Autogromacs.sh 200        (or give them as an argument; handy with nohup)
#
# Runs the WHOLE workflow: minimization -> equilibration -> production.
#
# PRODUCTION runs in SEGMENTS of 100 ns: it runs 100 ns, then extends
# +100 ns (convert-tpr -extend + mdrun -cpi -append) until the total is
# reached. With -append everything ends up in ONE continuous file, with no
# duplicated frames.
#
# RESUMABLE: if the run is interrupted, run the script again (with no
# arguments) and it continues from the last checkpoint. The total target is
# remembered in a state file: there is NO need to type it again.
#
# UNIVERSAL: finds gmx (PATH or GMXRC), checks for an NVIDIA GPU and sets the
# mdrun flags; uses all available cores (nproc). If the aggressive GPU mode
# fails (e.g. CUDA error #700), it retries automatically in a stable mode
# (bonded/update on CPU) from the last checkpoint.
# Compatible with GROMACS 2024 and 2026: -append is only used when a .cpt
# exists (2026 aborts with -append and no checkpoint) and experimental
# features (CUDA graphs) are disabled even if set in the environment.
#
# When finished, the trajectory gets a length-based name:
#       dynamic<N>ns.{xtc,tpr,edr,log,gro,cpt}      (e.g. dynamic200ns.xtc)
###############################################################################
set -uo pipefail

# ── 0a. Locate GROMACS ──────────────────────────────────────────────────────
# First the PATH; otherwise, the most common install locations.
if ! command -v gmx >/dev/null 2>&1; then
    for GMXRC in /usr/local/gromacs/bin/GMXRC \
                 /usr/local/gromacs/*/bin/GMXRC \
                 /opt/gromacs/bin/GMXRC \
                 "$HOME"/software/gromacs-*/bin/GMXRC \
                 "$HOME/gromacs/bin/GMXRC"; do
        if [[ -f "$GMXRC" ]]; then
            # shellcheck disable=SC1090
            source "$GMXRC"
            break
        fi
    done
fi
if ! command -v gmx >/dev/null 2>&1; then
    echo "❌ ERROR: 'gmx' not found in the PATH nor a GMXRC in the known locations."
    echo "   Check the GROMACS installation on this machine."
    exit 1
fi
GMX_VER=$(gmx --version 2>/dev/null | awk -F: '/GROMACS version/{gsub(/ /,"",$2); print $2; exit}')
echo "✓ GROMACS ${GMX_VER:-?}: $(command -v gmx)"

# No experimental features in production (CUDA graphs may be enabled from
# ~/.bashrc on some machines). Only affects this script.
if [[ -n "${GMX_CUDA_GRAPH:-}" ]]; then
    echo "ℹ GMX_CUDA_GRAPH was set in the environment: disabled (experimental)."
    unset GMX_CUDA_GRAPH
fi

# ── 0b. Check GPU (REQUIRED) and available cores ───────────────────────────
NTOMP=$(nproc 2>/dev/null || echo 4)
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
    echo "✓ GPU detected: $(nvidia-smi -L | head -n1)"
else
    echo "❌ ERROR: this pipeline requires an NVIDIA GPU and none was detected."
    echo "   Check the drivers (nvidia-smi) before running on this machine."
    exit 1
fi

# mdrun flags (same scheme for 2080 / 3080 / 5090; GROMACS adapts):
# HARD = EVERYTHING on GPU including -update: maximum performance
# SAFE = bonded/update on CPU: automatic rescue for CUDA #700 / instability
#        (explicit -bonded/-update cpu: in 2026 'auto' leaves both on GPU)
FLAGS_HARD="-pme gpu -nb gpu -bonded gpu -update gpu -pmefft gpu -ntmpi 1 -ntomp ${NTOMP} -pin on"
FLAGS_SAFE="-pme gpu -nb gpu -bonded cpu -update cpu -pmefft gpu -ntmpi 1 -ntomp ${NTOMP} -pin on"

# ── Fixed parameters ────────────────────────────────────────────────────────
CHUNK_NS=100                     # size of each segment (ns)
PROD_MDP="step5_production.mdp"
PROD="step5_production"          # deffnm base (renamed to dynamic<N>ns at the end)
STATE=".autogromacs_total_ns"    # the TOTAL target is stored here (for resuming)

# Check the minimum CHARMM-GUI inputs before starting
for f in "$PROD_MDP" topol.top index.ndx; do
    if [[ ! -f "$f" ]]; then
        echo "❌ ERROR: '$f' is missing. Are you inside the CHARMM-GUI gromacs/ folder?"
        exit 1
    fi
done

# ── 0c. Total ns target: argument > saved state > ask ───────────────────────
TOTAL_NS="${1:-}"
if [[ -z "$TOTAL_NS" && -f "$STATE" ]]; then
    TOTAL_NS=$(cat "$STATE")
    echo "🔄 Resume detected: saved target = ${TOTAL_NS} ns (from ${STATE})"
fi
if [[ -z "$TOTAL_NS" ]]; then
    if [[ -t 0 ]]; then
        read -rp "How many ns of production do you want to simulate? " TOTAL_NS
    else
        echo "❌ No interactive terminal and no previous target: give the ns as an argument."
        echo "   Example:  nohup bash Autogromacs.sh 200 > autogromacs.log 2>&1 &"
        exit 1
    fi
fi
if ! [[ "$TOTAL_NS" =~ ^[0-9]+$ ]] || (( TOTAL_NS < 1 )); then
    echo "❌ You must enter an integer number of ns >= 1 (e.g. 100, 200, 250)."; exit 1
fi
echo "$TOTAL_NS" > "$STATE"

# dt from the .mdp (fallback 0.002 ps)
DT=$(awk -F= '/^[[:space:]]*dt[[:space:]]*=/{gsub(/[^0-9.]/,"",$2); print $2; exit}' "$PROD_MDP" 2>/dev/null)
DT=${DT:-0.002}
ns_to_steps(){ awk -v ns="$1" -v dt="$DT" 'BEGIN{printf "%.0f", ns*1000/dt}'; }

# Read the REAL target from the current tpr (source of truth; avoids
# mismatches if a cut happens right between extending the tpr and updating
# any record)
tpr_target_ns(){
    local nsteps
    nsteps=$(gmx dump -s "${PROD}.tpr" 2>/dev/null \
             | grep -m1 -oE 'nsteps[^0-9-]*[0-9]+' | grep -oE '[0-9]+$')
    [[ -n "$nsteps" ]] || { echo 0; return; }
    awk -v s="$nsteps" -v dt="$DT" 'BEGIN{printf "%.0f", s*dt/1000}'
}

FINAL="dynamic${TOTAL_NS}ns"

echo "═══════════════════════════════════════════════════════════════════════"
echo " Autogromacs — target: ${TOTAL_NS} ns   (dt=${DT} ps · ${CHUNK_NS} ns segments)"
echo " Hardware: GPU + ${NTOMP} CPU threads"
echo " Final output: ${FINAL}.{xtc,tpr,...}"
echo "═══════════════════════════════════════════════════════════════════════"

# If the final output already exists, do not redo it (avoids overwriting results)
if [[ -f ${FINAL}.xtc ]]; then
    echo "⚠ ${FINAL}.xtc already exists in this folder. Delete it or choose another length to redo."
    exit 0
fi

# ── 1. Minimization (skipped if already done) ───────────────────────────────
if [[ ! -f step4.0_minimization.gro ]]; then
    echo "=== MINIMIZATION ==="
    gmx grompp -f step4.0_minimization.mdp -o step4.0_minimization.tpr \
        -c step3_input.gro -r step3_input.gro -p topol.top -n index.ndx -maxwarn 1 || exit 1
    gmx mdrun -v -deffnm step4.0_minimization -ntmpi 1 -ntomp "${NTOMP}" || exit 1
else
    echo "✓ Minimization already done, skipped."
fi

# ── 2. Equilibration (skipped if already done) ──────────────────────────────
if [[ ! -f step4.1_equilibration.gro ]]; then
    echo "=== EQUILIBRATION ==="
    gmx grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr \
        -c step4.0_minimization.gro -r step3_input.gro -p topol.top -n index.ndx -maxwarn 1 || exit 1
    gmx mdrun -v -deffnm step4.1_equilibration $FLAGS_SAFE || exit 1
else
    echo "✓ Equilibration already done, skipped."
fi

# ── 3. Production in ${CHUNK_NS} ns segments ────────────────────────────────
# Prepare the tpr of the first segment if it does not exist yet
if [[ ! -f ${PROD}.tpr ]]; then
    FIRST_NS=$(( TOTAL_NS < CHUNK_NS ? TOTAL_NS : CHUNK_NS ))
    echo "=== PRODUCTION: preparing first segment of ${FIRST_NS} ns ==="
    # continuation=yes in the .mdp: take the exact velocities/state from the
    # equilibration .cpt (more precise than the .gro) if it exists
    EQ_CPT=""
    [[ -f step4.1_equilibration.cpt ]] && EQ_CPT="-t step4.1_equilibration.cpt"
    gmx grompp -f "$PROD_MDP" -o ${PROD}.tpr \
        -c step4.1_equilibration.gro $EQ_CPT -p topol.top -n index.ndx -maxwarn 1 || exit 1
    # Set nsteps of the first segment by writing to a temporary file, then mv
    gmx convert-tpr -s ${PROD}.tpr -nsteps "$(ns_to_steps "$FIRST_NS")" -o ${PROD}_tmp.tpr || exit 1
    mv -f ${PROD}_tmp.tpr ${PROD}.tpr
fi

# Current target ALWAYS read from the tpr itself (robust to interruptions)
TARGET_NS=$(tpr_target_ns)
if (( TARGET_NS == 0 )); then
    echo "❌ ERROR: could not read nsteps from ${PROD}.tpr (corrupt file?)."; exit 1
fi
if (( TARGET_NS > TOTAL_NS )); then
    echo "⚠ The current tpr already targets ${TARGET_NS} ns, more than the goal (${TOTAL_NS} ns)."
    echo "  It will run up to ${TARGET_NS} ns and be renamed dynamic${TARGET_NS}ns.*"
    TOTAL_NS=$TARGET_NS
    FINAL="dynamic${TOTAL_NS}ns"
    echo "$TOTAL_NS" > "$STATE"
fi

# GPU monitor in the background (optional, does not affect the simulation)
GPU_PID=""
nvidia-smi dmon -s puct -d 10 >> gpu_monitor.log 2>&1 &
GPU_PID=$!
cleanup(){ [[ -n "$GPU_PID" ]] && kill "$GPU_PID" 2>/dev/null || true; }
trap cleanup EXIT

# Continuation flags: -cpi/-append ONLY if the .cpt exists. GROMACS 2026
# aborts with "-append" and no checkpoint (other versions may tolerate it).
# Evaluated on every attempt so the SAFE retry continues from the .cpt left
# by HARD.
cpt_flags(){
    if [[ -f ${PROD}.cpt ]]; then
        echo "-cpi ${PROD}.cpt -append"
    fi
}

# Interruptions (Ctrl-C, kill, shutdown): mdrun stops cleanly, writes a
# checkpoint and exits with code 1 — just like a real failure. To NOT mistake
# it for a GPU failure (and jump to SAFE) it is detected by:
#   a) the signal reaching the script itself (Ctrl-C in the terminal), or
#   b) the "Received the ... signal" message mdrun writes to its log
#      (when only the mdrun process is killed, e.g. with nohup + kill).
INTERRUPTED=0
trap 'INTERRUPTED=1' INT TERM HUP
# Only the log lines written in THIS attempt (LOG_START) are checked, so a
# signal notice from a previous run is not confused with the current one.
log_lines(){ local n; n=$(wc -l < ${PROD}.log 2>/dev/null); echo "${n:-0}"; }
log_says_signal(){
    local start=$LOG_START
    (( $(log_lines) < start )) && start=0     # mdrun created a new log
    awk -v a="$start" 'NR>a && /Received the .* signal/{f=1} END{exit !f}' \
        ${PROD}.log 2>/dev/null
}
stop_if_interrupted(){
    if (( INTERRUPTED )) || log_says_signal; then
        echo "⏹ mdrun interrupted by a signal; the checkpoint was saved."
        echo "   Run 'bash Autogromacs.sh' again (no arguments) to resume."
        exit 130
    fi
}

# One mdrun segment with automatic retry in stable mode.
run_chunk(){
    echo "── mdrun [HARD: everything on GPU, maximum performance] ──"
    LOG_START=$(log_lines)
    gmx mdrun -v -deffnm ${PROD} -s ${PROD}.tpr $(cpt_flags) $FLAGS_HARD && return 0
    stop_if_interrupted
    echo "⚠ HARD mode failed (CUDA #700 / instability?)."
    echo "🚑 Retrying in SAFE MODE (bonded/update on CPU) from the checkpoint..."
    LOG_START=$(log_lines)
    gmx mdrun -v -deffnm ${PROD} -s ${PROD}.tpr $(cpt_flags) $FLAGS_SAFE
}

# Segment loop: run/continue up to TARGET_NS, then extend, until TOTAL_NS
while :; do
    echo "── Production: running/continuing up to ${TARGET_NS} ns ──"
    if ! run_chunk; then
        stop_if_interrupted
        echo "❌ mdrun failed even in stable mode. Check ${PROD}.log"
        echo "   (look for LINCS warnings or positive energy: the system may have blown up)."
        exit 1
    fi

    if (( TARGET_NS >= TOTAL_NS )); then
        echo "✓ Production complete: ${TARGET_NS} ns reached."
        break
    fi

    REMAIN=$(( TOTAL_NS - TARGET_NS ))
    NEXT_NS=$(( REMAIN < CHUNK_NS ? REMAIN : CHUNK_NS ))
    EXT_PS=$(( NEXT_NS * 1000 ))
    echo "=== Extending +${NEXT_NS} ns (${EXT_PS} ps) ==="
    # Write to a temporary file, then mv: never leave a half-written tpr
    gmx convert-tpr -s ${PROD}.tpr -extend "$EXT_PS" -o ${PROD}_tmp.tpr || exit 1
    mv -f ${PROD}_tmp.tpr ${PROD}.tpr
    TARGET_NS=$(tpr_target_ns)
done

cleanup; trap - EXIT

# ── 4. Rename final output to dynamic<N>ns.* ────────────────────────────────
echo "=== Renaming output to ${FINAL}.* ==="
for ext in xtc trr tpr edr log gro cpt; do
    [[ -f ${PROD}.${ext} ]] && mv -f ${PROD}.${ext} ${FINAL}.${ext} && echo "  ✓ ${PROD}.${ext} → ${FINAL}.${ext}"
done
[[ -f ${PROD}_prev.cpt ]] && mv -f ${PROD}_prev.cpt ${FINAL}_prev.cpt
rm -f "$STATE"

echo "═══════════════════════════════════════════════════════════════════════"
echo " DONE. ${TOTAL_NS} ns trajectory:  ${FINAL}.xtc  /  ${FINAL}.tpr"
echo " Next step: bash basic_analyses_with_names.sh  (auto-detects ${FINAL})"
echo "═══════════════════════════════════════════════════════════════════════"
