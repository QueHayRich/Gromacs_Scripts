#!/bin/bash
###############################################################################
# basic_analyses_with_names.sh — standard post-processing for Autogromacs trajectories
#
# Author : Ricardo Aguilar-Garay
# Contact: ricardoaguilargaray@gmail.com
#          (questions, bug reports, and suggestions are welcome)
#
# Contributor: Dra. Dulce Estefanía Nicolás Álvarez
#              Departamento de Fisiología, Escuela Nacional de Ciencias
#              Biológicas, Instituto Politécnico Nacional (ENCB-IPN), Mexico
#              dnicolas@ipn.mx
#
# Companion of Autogromacs.sh. Run it inside the CHARMM-GUI gromacs/ folder
# after production finishes:
#       bash basic_analyses_with_names.sh               (auto-detects dynamic<N>ns.*)
#       bash basic_analyses_with_names.sh dynamic150ns  (or name the basename explicitly)
#
# Same analyses as the original lab script, but using GROMACS default groups
# BY NAME instead of by number (the numbers shift between systems):
#   1. Centered trajectory     center Protein, write System      (old: 1 0)
#   2. RMSD of the ligand      fit Backbone, calc Ligand          (old: 4 13)
#   3. H-bonds                 Protein – Ligand                   (old: 1 13)
#      (ligands with only C/H, e.g. microplastics: Protein–Ligand contacts)
#   4. RMSF per residue        Protein                            (old: 1)
#   5. Radius of gyration      Protein                            (old: 1)
#   6. PDB snapshots           System at start / middle / end     (old: 0)
#
# "Ligand" = everything in the CHARMM-GUI SOLU group that is not Protein,
# so it works for any ligand name. If there is no ligand, RMSD is computed
# on the Backbone and H-bonds within the Protein.
#
# Results go to Analyses/ and, if the folders ../../rmsd, ../../hbond,
# ../../rmsf, ../../gyrate and ../../conformation exist, are copied there
# too, as the original script did (renamed with the system name).
###############################################################################
set -uo pipefail

# ── Locate GROMACS: gmx (thread-MPI build) or gmx_mpi (real MPI build) ──────
find_gmx(){ command -v gmx >/dev/null 2>&1 && echo gmx && return
            command -v gmx_mpi >/dev/null 2>&1 && echo gmx_mpi; }
GMX=$(find_gmx)
if [[ -z "$GMX" ]]; then
    for GMXRC in /usr/local/gromacs/bin/GMXRC /usr/local/gromacs/*/bin/GMXRC \
                 /opt/gromacs/bin/GMXRC "$HOME/gromacs/bin/GMXRC"; do
        # shellcheck disable=SC1090
        [[ -f "$GMXRC" ]] && source "$GMXRC" && GMX=$(find_gmx) && [[ -n "$GMX" ]] && break
    done
fi
[[ -z "$GMX" ]] && { echo "❌ ERROR: neither 'gmx' nor 'gmx_mpi' found."; exit 1; }
echo "✓ GROMACS: $(command -v "$GMX")"

# ── Names: files in Analyses/ are tagged with the current folder name
#    (e.g. rmsd_gromacs_200ns.xvg); copies to the shared folders use the
#    system name (the folder that contains gromacs/) so systems don't clash.
NAME="${PWD##*/}"
SYS="$NAME"
[[ "$NAME" == "gromacs" ]] && SYS="$(basename "$(dirname "$PWD")")"

# ── Detect trajectory: argument > newest dynamic<N>ns.* > step5_production ──
BASE="${1:-}"; TPR=""; XTC=""
if [[ -n "$BASE" ]]; then
    TPR="${BASE}.tpr"; XTC="${BASE}.xtc"
else
    for f in $(ls -t dynamic*ns.tpr 2>/dev/null); do
        [[ -f "${f%.tpr}.xtc" ]] && { BASE="${f%.tpr}"; TPR="$f"; XTC="${BASE}.xtc"; break; }
    done
    if [[ -z "$TPR" && -f step5_production.xtc ]]; then
        TPR=$(ls -t step5_production*.tpr 2>/dev/null | head -n1)
        XTC="step5_production.xtc"; BASE="${TPR%.tpr}"
    fi
fi
if [[ -z "$TPR" || ! -f "$TPR" || ! -f "$XTC" ]]; then
    echo "❌ ERROR: no trajectory found (dynamic<N>ns.{tpr,xtc} or step5_production.*)."
    exit 1
fi
[[ -f index.ndx ]] || { echo "❌ ERROR: index.ndx not found. Run inside gromacs/."; exit 1; }

# ── Simulation length (for file labels and snapshot spacing) ────────────────
if [[ "$BASE" =~ dynamic([0-9]+)ns ]]; then
    TOTAL_NS="${BASH_REMATCH[1]}"
else
    TOTAL_NS=$("$GMX" dump -s "$TPR" 2>/dev/null | awk '
        $1=="nsteps"{s=$3} $1=="delta_t"||$1=="dt"{d=$3} END{printf "%.0f", s*d/1000}')
fi
LABEL="${TOTAL_NS}ns"

OUT="Analyses"; mkdir -p "$OUT"
NDX="${OUT}/index_analysis.ndx"
CENTERED="${OUT}/${BASE}_centered.xtc"
echo "═══ ${SYS} · ${XTC} · ${TOTAL_NS} ns → ${OUT}/ ═══"

# ── 0. Index: GROMACS default groups + SOLU + Ligand (= SOLU & !Protein) ────
echo q | "$GMX" make_ndx -f "$TPR" -o "$NDX" >/dev/null 2>&1
awk '/^\[/{p = ($0 ~ /^\[ *SOLU *\]/)} p' index.ndx >> "$NDX"
NG=$(grep -c '^\[' "$NDX")              # the new group will get number NG
printf '"SOLU" & ! "Protein"\nname %d Ligand\nq\n' "$NG" \
    | "$GMX" make_ndx -f "$TPR" -n "$NDX" -o "${OUT}/tmp_index.ndx" >/dev/null 2>&1
mv -f "${OUT}/tmp_index.ndx" "$NDX" 2>/dev/null

if grep -q '^\[ Ligand \]' "$NDX"; then
    LIG="Ligand"; echo "✓ Ligand group created (SOLU minus Protein)."
else
    LIG=""; echo "✓ No ligand in SOLU: RMSD on Backbone, H-bonds within Protein."
fi

# ── 1. Centered trajectory: center Protein, write System ────────────────────
if [[ ! -f "$CENTERED" ]]; then
    printf "Protein\nSystem\n" | "$GMX" trjconv -s "$TPR" -f "$XTC" -n "$NDX" \
        -o "$CENTERED" -pbc mol -center || exit 1
fi

# ── 2. RMSD: fit Backbone, calc Ligand (Backbone if no ligand) ──────────────
printf "Backbone\n%s\n" "${LIG:-Backbone}" | "$GMX" rms -s "$TPR" -f "$CENTERED" \
    -n "$NDX" -o "${OUT}/rmsd_${NAME}_${LABEL}.xvg" -tu ns || exit 1

# ── 3. H-bonds Protein–Ligand (new -r/-t syntax, legacy fallback) ───────────
# Ligands made only of C and H (polyethylene, polypropylene, polystyrene...)
# have no H-bond donors/acceptors and gmx hbond refuses to run. For those,
# Protein–Ligand contacts (atom pairs < 0.4 nm) are computed instead.
POLAR=1
if [[ -n "$LIG" ]]; then
    "$GMX" editconf -f "$TPR" -o "${OUT}/.ref.gro" >/dev/null 2>&1
    POLAR=$(awk 'FNR==NR{ if(/^\[/){g=($0 ~ /^\[ *Ligand *\]/); next}
                          if(g) for(i=1;i<=NF;i++) L[$i]=1; next }
                 FNR>2 && ((FNR-2) in L){ a=substr($0,11,5); gsub(/ /,"",a)
                          if(a ~ /^[NOSF]/) n++ }
                 END{print n+0}' "$NDX" "${OUT}/.ref.gro")
    rm -f "${OUT}/.ref.gro"
fi

if (( POLAR > 0 )); then
    HB="${OUT}/hbond_${NAME}_${LABEL}.xvg"
    "$GMX" hbond -f "$CENTERED" -s "$TPR" -n "$NDX" -num "$HB" -tu ns \
        -r Protein -t "${LIG:-Protein}" 2>/dev/null \
    || printf "Protein\n%s\n" "${LIG:-Protein}" | "$GMX" hbond -f "$CENTERED" \
        -s "$TPR" -n "$NDX" -num "$HB" -tu ns \
    || echo "⚠ H-bond analysis failed; continuing with the rest."
else
    echo "✓ Ligand has no N/O/S/F atoms (e.g. microplastic): H-bonds skipped,"
    echo "  computing Protein–Ligand contacts (< 0.4 nm) instead."
    printf "Protein\nLigand\n" | "$GMX" mindist -s "$TPR" -f "$CENTERED" -n "$NDX" \
        -d 0.4 -tu ns -on "${OUT}/contacts_${NAME}_${LABEL}.xvg" \
        -od "${OUT}/mindist_${NAME}_${LABEL}.xvg" \
    || echo "⚠ Contact analysis failed; continuing with the rest."
fi

# ── 4. RMSF per residue (Protein) ───────────────────────────────────────────
printf "Protein\n" | "$GMX" rmsf -s "$TPR" -f "$CENTERED" -n "$NDX" \
    -o "${OUT}/rmsf_${NAME}_${LABEL}.xvg" -res || exit 1

# ── 5. Radius of gyration (Protein) ─────────────────────────────────────────
printf "Protein\n" | "$GMX" gyrate -s "$TPR" -f "$CENTERED" -n "$NDX" \
    -o "${OUT}/gyrate_${NAME}_${LABEL}.xvg" || exit 1

# ── 6. PDB snapshots at t = 0, middle and end (System) ──────────────────────
SNAP_DT=$(( TOTAL_NS * 500 )); (( SNAP_DT < 1 )) && SNAP_DT=1
printf "System\n" | "$GMX" trjconv -s "$TPR" -f "$CENTERED" -n "$NDX" \
    -dt "$SNAP_DT" -o "${OUT}/trj_${NAME}_${LABEL}.pdb" -pbc mol || exit 1

# ── Copy results to the shared folders, as the original script did ─────────
for pair in rmsd:rmsd hbond:hbond contacts:contacts mindist:contacts \
            rmsf:rmsf gyrate:gyrate trj:conformation; do
    src="${pair%%:*}"; dst="../../${pair##*:}"
    [[ -d "$dst" ]] || continue
    for f in "${OUT}/${src}_${NAME}_${LABEL}".*; do
        [[ -f "$f" ]] || continue
        cp "$f" "$dst/${src}_${SYS}_${LABEL}.${f##*.}"
    done
done

echo "═══ DONE: rmsd, hbond/contacts, rmsf, gyrate (.xvg) and trj (.pdb) in ${OUT}/ ═══"
