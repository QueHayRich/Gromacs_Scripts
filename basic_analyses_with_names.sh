#!/bin/bash
###############################################################################
# basic_analyses_with_names.sh — standard post-processing for Autogromacs trajectories
#                                (Spanish version: basic_analyses_with_names.es.sh)
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
# Companion of Autogromacs.sh. Run it inside the CHARMM-GUI gromacs/ folder
# after production finishes:
#       bash basic_analyses_with_names.sh               (auto-detects dynamic<N>ns.*)
#       bash basic_analyses_with_names.sh dynamic150ns  (or name the basename explicitly)
#
# Options (only for special cases; by default everything is detected):
#       --ligand NAME[,NAME]    analyse only these molecules (from [ molecules ]
#                               in topol.top) as ligands, e.g. to leave out a
#                               cofactor, or to treat a peptide chain as ligand
#       --receptor NAME[,NAME]  receptor molecule(s) (default: all protein chains)
#
# Works for any protein and any ligand (microplastic, drug, metabolite...).
# GROMACS groups are selected BY NAME, never by number (the numbers shift
# between systems), so nothing has to be edited:
#   1. Centered trajectory     center Receptor, write System
#   2. RMSD of each ligand     fit receptor backbone, calc ligand
#   3. H-bonds                 Receptor – each ligand
#      (ligands with no N or O atoms, e.g. polyethylene, polystyrene, PVC,
#       PTFE: Receptor–ligand contacts and minimum distance instead)
#   4. RMSF per residue        Receptor
#   5. Radius of gyration      Receptor
#   6. PDB snapshots           System at start / middle / end
#
# Receptor and ligands are read from the [ molecules ] section of topol.top:
# proteins = receptor; water and ions are ignored; molecules with more than 5
# copies (lipids, cosolvents) are ignored; every other molecule type is a
# ligand. Each ligand gets its own group (Ligand_<name>) and its own results,
# so a microplastic, an endogenous ligand and a cofactor are never mixed.
# Copies of the same molecule go into one group. With one ligand the output
# names are rmsd_<folder>_<N>ns.xvg; with several, the ligand name is added:
# rmsd_<folder>_<ligand>_<N>ns.xvg. If there is no ligand, RMSD is computed
# on the receptor backbone and H-bonds within the receptor.
#
# Results go to Analyses/. Optionally, if the folders ../../rmsd, ../../hbond,
# ../../contacts, ../../rmsf, ../../gyrate and ../../conformation exist (a
# project folder that gathers many systems), they are copied there too,
# renamed with the system name.
###############################################################################
set -uo pipefail

# ── Arguments ───────────────────────────────────────────────────────────────
BASE=""; LIGAND_OPT=""; RECEPTOR_OPT=""
usage(){ awk 'NR>=3 && /^####/{exit} NR>=3' "$0" | sed 's/^# \{0,1\}//'; }
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ligand|--receptor)
            [[ -n "${2:-}" ]] || { echo "❌ ERROR: $1 needs a value (see --help)."; exit 1; }
            [[ "$1" == "--ligand" ]] && LIGAND_OPT="$2" || RECEPTOR_OPT="$2"
            shift ;;
        -h|--help) usage; exit 0 ;;
        -*)        echo "❌ ERROR: unknown option '$1' (see --help)."; exit 1 ;;
        *)         BASE="$1" ;;
    esac
    shift
done

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
PYTHON=$(command -v python3 || command -v python || true)
[[ -n "$PYTHON" ]] || { echo "❌ ERROR: python3 not found (needed to read topol.top)."; exit 1; }

# ── Names: files in Analyses/ are tagged with the current folder name
#    (e.g. rmsd_gromacs_200ns.xvg); copies to the shared folders use the
#    system name (the folder that contains gromacs/) so systems don't clash.
NAME="${PWD##*/}"
SYS="$NAME"
[[ "$NAME" == "gromacs" ]] && SYS="$(basename "$(dirname "$PWD")")"

# ── Detect trajectory: argument > newest dynamic<N>ns.* > step5_production ──
TPR=""; XTC=""
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
[[ -f topol.top ]] || { echo "❌ ERROR: topol.top not found. Run inside gromacs/."; exit 1; }

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

# ── 0. Index: GROMACS default groups + Receptor, Receptor_BB, Ligand_<name> ──
# The helper reads [ molecules ] from topol.top (following #include), appends
# the groups to the index and prints KEY='value' lines (evaluated below).
echo q | "$GMX" make_ndx -f "$TPR" -o "$NDX" >/dev/null 2>&1
[[ -s "$NDX" ]] || { echo "❌ ERROR: gmx make_ndx could not read ${TPR}."; exit 1; }
DETECT=$("$PYTHON" - "$LIGAND_OPT" "$RECEPTOR_OPT" "$NDX" <<'PYDETECT'
import os, re, sys
lig_opt = [x.strip() for x in sys.argv[1].split(",") if x.strip()]
rec_opt = [x.strip() for x in sys.argv[2].split(",") if x.strip()]
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

moltypes, molecules = {}, []
def parse(path, seen=None):
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
                cand = os.path.join(base, m.group(1))
                parse(cand if os.path.isfile(cand) else m.group(1), seen)
            continue
        if line.startswith("#"): continue
        m = re.match(r"\[\s*(\w+)\s*\]", line)
        if m: section = m.group(1).lower(); continue
        f = line.split()
        if section == "moleculetype":
            cur = f[0]; moltypes[cur] = []; section = "moleculetype_done"
        elif section == "atoms" and cur is not None and len(f) >= 5:
            try: charge = float(f[6]) if len(f) > 6 else 0.0
            except ValueError: charge = 0.0
            moltypes[cur].append((f[3], f[4], charge, f[2]))    # resname, atom, charge, resnr
        elif section == "molecules" and len(f) >= 2:
            try: molecules.append((f[0], int(f[1])))
            except ValueError: pass

parse("topol.top")
if not molecules: fail("could not read the [ molecules ] section of topol.top")
missing = [n for n, _ in molecules if n not in moltypes]
if missing: fail("moleculetype(s) not found in the included .itp files: " + ", ".join(missing))

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
    lig = [b for b in cands if b["count"] <= 5]
    bulk = [b for b in cands if b["count"] > 5]
    if bulk:
        notes.append("molecules with >5 copies are not analysed as ligands (use --ligand to include them): "
                     + ", ".join("%s(x%d)" % (b["name"], b["count"]) for b in bulk))

lig_ids = {id(b) for b in lig}
if rec_opt:
    rec = [b for b in blocks if match(b, rec_opt) and id(b) not in lig_ids]
    if not rec: fail("--receptor %s does not match any molecule (or it is the ligand)." % ",".join(rec_opt))
else:
    rec = [b for b in blocks if b["kind"] == "protein" and id(b) not in lig_ids]
    if not rec: fail("no protein found to use as receptor. Use --receptor NAME.")
rec_ids = {id(b) for b in rec}
ignored = [b for b in blocks if id(b) not in lig_ids | rec_ids]

def indices(blks, pick=lambda a: True):
    out = []
    for b in blks:
        atoms = moltypes[b["name"]]
        for c in range(b["count"]):
            off = b["start"] + c * b["natoms"]
            out += [off + i for i, a in enumerate(atoms) if pick(a)]
    return out

def group(name, idx):
    s = "[ %s ]\n" % name
    for i in range(0, len(idx), 15): s += " ".join("%6d" % x for x in idx[i:i+15]) + "\n"
    return s

rec_idx = indices(rec)
bb_idx  = indices(rec, lambda a: a[1].upper() in ("N", "CA", "C"))
if not bb_idx: bb_idx = indices(rec, lambda a: not a[1].upper().startswith("H"))
groups = [("Receptor", rec_idx), ("Receptor_BB", bb_idx)]

# One group per ligand molecule type (all its copies together)
ligs, used = [], set()
for b in lig:
    known = next((L for L in ligs if L["name"] == b["name"]), None)
    if known: known["blocks"].append(b); continue
    tag = re.sub(r"[^A-Za-z0-9]+", "_", b["name"]).strip("_") or "LIG"
    while tag.upper() in used: tag += "_"
    used.add(tag.upper())
    ligs.append(dict(name=b["name"], tag=tag, blocks=[b]))
for L in ligs:
    idx = indices(L["blocks"])
    copies = sum(b["count"] for b in L["blocks"])
    polar = any(a[1].upper()[:1] in ("N", "O") for b in L["blocks"] for a in moltypes[b["name"]])
    L.update(copies=copies, natoms=len(idx) // copies, polar=polar)
    groups.append(("Ligand_" + L["tag"], idx))
    if copies > 1: notes.append("ligand %s has %d copies: analysed together" % (L["name"], copies))

with open(ndx_out, "a") as fh:
    for name, idx in groups: fh.write(group(name, idx))

def show(b): return b["name"] + ("(x%d)" % b["count"] if b["count"] > 1 else "")
nres = sum(len({a[3] for a in moltypes[b["name"]]}) * b["count"] for b in rec)
print("REC_SHOW=" + q("%s (%d residues)" % (", ".join(show(b) for b in rec), nres)))
print("LIG_SHOW=" + q(" · ".join("%s%s (%d atoms%s)" % (L["name"], "(x%d)" % L["copies"] if L["copies"] > 1 else "",
                                L["natoms"], "" if L["polar"] else ", no N/O → contacts") for L in ligs)))
print("IGN_SHOW=" + q(", ".join(show(b) for b in ignored)))
print("LIGS=" + q(" ".join("%s:%d" % (L["tag"], L["polar"]) for L in ligs)))
print("NOTES=" + q(" | ".join(notes)))
PYDETECT
)
if [[ "$DETECT" != *"LIGS="* && "$DETECT" != *"FAIL="* ]]; then
    echo "❌ ERROR: could not read topol.top (python helper failed)."; exit 1
fi
FAIL=""; eval "$DETECT"
[[ -n "$FAIL" ]] && { echo "❌ ERROR: ${FAIL}"; exit 1; }

read -ra LIG_LIST <<< "$LIGS"
NLIG=${#LIG_LIST[@]}
echo " Receptor : ${REC_SHOW}"
echo " Ligands  : ${LIG_SHOW:-none (apo): RMSD on the receptor backbone, H-bonds within the receptor}"
[[ -n "$IGN_SHOW" ]] && echo " Ignored  : ${IGN_SHOW}"
[[ -n "$NOTES" ]] && echo " ℹ ${NOTES}"

# Output name of an analysis: the ligand name is added only if there are several
out(){ # $1 = analysis, $2 = ligand tag (may be empty)
    if (( NLIG > 1 )); then echo "${OUT}/${1}_${NAME}_${2}_${LABEL}.xvg"
    else echo "${OUT}/${1}_${NAME}_${LABEL}.xvg"; fi
}

# ── 1. Centered trajectory: center Receptor, write System ───────────────────
if [[ ! -f "$CENTERED" ]]; then
    printf "Receptor\nSystem\n" | "$GMX" trjconv -s "$TPR" -f "$XTC" -n "$NDX" \
        -o "$CENTERED" -pbc mol -center || exit 1
fi

# ── 2-3. Per ligand: RMSD, then H-bonds or contacts ─────────────────────────
# gmx hbond only uses N and O as donors/acceptors. Ligands with no N or O
# atoms (polyethylene, polypropylene, polystyrene, PVC, PTFE...) have none,
# and gmx hbond refuses to run. For those, Receptor–ligand contacts (atom
# pairs < 0.4 nm) and the minimum distance are computed instead.
hbonds(){ # $1 = second group, $2 = output (new -r/-t syntax, legacy fallback)
    "$GMX" hbond -f "$CENTERED" -s "$TPR" -n "$NDX" -num "$2" -tu ns \
        -r Receptor -t "$1" 2>/dev/null \
    || printf "Receptor\n%s\n" "$1" | "$GMX" hbond -f "$CENTERED" \
        -s "$TPR" -n "$NDX" -num "$2" -tu ns \
    || echo "⚠ H-bond analysis failed ($1); continuing with the rest."
}
contacts(){ # $1 = ligand group, $2 = ligand tag
    printf "Receptor\n%s\n" "$1" | "$GMX" mindist -s "$TPR" -f "$CENTERED" -n "$NDX" \
        -d 0.4 -tu ns -on "$(out contacts "$2")" -od "$(out mindist "$2")" \
    || echo "⚠ Contact analysis failed ($1); continuing with the rest."
}

if (( NLIG == 0 )); then
    printf "Receptor_BB\nReceptor_BB\n" | "$GMX" rms -s "$TPR" -f "$CENTERED" \
        -n "$NDX" -o "$(out rmsd "")" -tu ns || exit 1
    hbonds Receptor "$(out hbond "")"
fi
for entry in "${LIG_LIST[@]}"; do
    tag="${entry%:*}"; polar="${entry##*:}"; grp="Ligand_${tag}"
    (( NLIG > 1 )) && echo "── Ligand ${tag} ──"
    printf "Receptor_BB\n%s\n" "$grp" | "$GMX" rms -s "$TPR" -f "$CENTERED" \
        -n "$NDX" -o "$(out rmsd "$tag")" -tu ns || exit 1
    if (( polar )); then
        hbonds "$grp" "$(out hbond "$tag")"
    else
        echo "✓ ${tag} has no N or O atoms (e.g. polyolefin microplastic): H-bonds skipped,"
        echo "  computing Receptor–ligand contacts (< 0.4 nm) instead."
        contacts "$grp" "$tag"
    fi
done

# ── 4. RMSF per residue (Receptor) ──────────────────────────────────────────
printf "Receptor\n" | "$GMX" rmsf -s "$TPR" -f "$CENTERED" -n "$NDX" \
    -o "${OUT}/rmsf_${NAME}_${LABEL}.xvg" -res || exit 1

# ── 5. Radius of gyration (Receptor) ────────────────────────────────────────
printf "Receptor\n" | "$GMX" gyrate -s "$TPR" -f "$CENTERED" -n "$NDX" \
    -o "${OUT}/gyrate_${NAME}_${LABEL}.xvg" || exit 1

# ── 6. PDB snapshots at t = 0, middle and end (System) ──────────────────────
SNAP_DT=$(( TOTAL_NS * 500 )); (( SNAP_DT < 1 )) && SNAP_DT=1
printf "System\n" | "$GMX" trjconv -s "$TPR" -f "$CENTERED" -n "$NDX" \
    -dt "$SNAP_DT" -o "${OUT}/trj_${NAME}_${LABEL}.pdb" -pbc mol || exit 1

# ── Optional: copy results to the shared project folders, if they exist ─────
for pair in rmsd:rmsd hbond:hbond contacts:contacts mindist:contacts \
            rmsf:rmsf gyrate:gyrate trj:conformation; do
    src="${pair%%:*}"; dst="../../${pair##*:}"
    [[ -d "$dst" ]] || continue
    for f in "${OUT}/${src}_${NAME}_"*; do
        [[ -f "$f" ]] || continue
        rest="${f##*/}"; rest="${rest#"${src}_${NAME}_"}"     # [<ligand>_]<N>ns.<ext>
        [[ "$rest" == "${LABEL}."* || "$rest" == *"_${LABEL}."* ]] || continue
        cp "$f" "$dst/${src}_${SYS}_${rest}"
    done
done

echo "═══ DONE: rmsd, hbond/contacts, rmsf, gyrate (.xvg) and trj (.pdb) in ${OUT}/ ═══"
