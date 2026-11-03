#!/bin/bash
###############################################################################
# basic_analyses_with_names.es.sh — análisis estándar de trayectorias de Autogromacs
#                                   (versión en español de basic_analyses_with_names.sh)
#
# Autor  : Ricardo Aguilar-Garay
# Contacto: ricardoaguilargaray@gmail.com
#          (dudas, reportes de errores y sugerencias son bienvenidos)
#
# Contribuidora: Dulce Estefanía Nicolás Álvarez
#                Departamento de Fisiología, Escuela Nacional de Ciencias
#                Biológicas, Instituto Politécnico Nacional (ENCB-IPN), México
#                dnicolas@ipn.mx
#
# Compañero de Autogromacs.es.sh. Ejecútalo dentro de la carpeta gromacs/ de
# CHARMM-GUI cuando termine la producción:
#       bash basic_analyses_with_names.es.sh               (detecta dynamic<N>ns.* solo)
#       bash basic_analyses_with_names.es.sh dynamic150ns  (o indica el nombre base)
#
# Opciones (solo para casos especiales; por defecto todo se detecta solo):
#       --ligand NOMBRE[,NOMBRE]    analizar solo estas moléculas (de [ molecules ]
#                                   en topol.top) como ligandos, p.ej. para dejar
#                                   fuera un cofactor, o para tratar una cadena
#                                   peptídica como ligando
#       --receptor NOMBRE[,NOMBRE]  molécula(s) receptora(s) (por defecto: todas
#                                   las cadenas de proteína)
#
# Sirve para cualquier proteína y cualquier ligando (microplástico, fármaco,
# metabolito...). Los grupos de GROMACS se eligen POR NOMBRE, nunca por
# número (los números cambian entre sistemas), así que no hay que editar nada:
#   1. Trayectoria centrada    centrar Receptor, escribir System
#   2. RMSD de cada ligando    ajustar backbone del receptor, calcular ligando
#   3. Puentes de H            Receptor – cada ligando
#      (ligandos sin átomos de N ni O, p.ej. polietileno, poliestireno, PVC,
#       PTFE: en su lugar, contactos y distancia mínima Receptor–ligando)
#   4. RMSF por residuo        Receptor
#   5. Radio de giro           Receptor
#   6. Snapshots PDB           System al inicio / mitad / final
#
# Receptor y ligandos se leen de la sección [ molecules ] de topol.top:
# proteínas = receptor; agua e iones se ignoran; las moléculas con más de 5
# copias (lípidos, cosolventes) se ignoran; cualquier otro tipo de molécula
# es un ligando. Cada ligando tiene su propio grupo (Ligand_<nombre>) y sus
# propios resultados, así que un microplástico, un ligando endógeno y un
# cofactor nunca se mezclan. Las copias de una misma molécula van en un solo
# grupo. Con un ligando los archivos se llaman rmsd_<carpeta>_<N>ns.xvg; con
# varios se agrega el nombre del ligando: rmsd_<carpeta>_<ligando>_<N>ns.xvg.
# Si no hay ligando, el RMSD se calcula sobre el backbone del receptor y los
# puentes de H dentro del receptor.
#
# Los resultados van a Analyses/. Opcionalmente, si existen las carpetas
# ../../rmsd, ../../hbond, ../../contacts, ../../rmsf, ../../gyrate y
# ../../conformation (una carpeta de proyecto que reúne muchos sistemas),
# también se copian ahí, renombrados con el nombre del sistema.
###############################################################################
set -uo pipefail

# ── Argumentos ──────────────────────────────────────────────────────────────
BASE=""; LIGAND_OPT=""; RECEPTOR_OPT=""
usage(){ awk 'NR>=3 && /^####/{exit} NR>=3' "$0" | sed 's/^# \{0,1\}//'; }
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ligand|--receptor)
            [[ -n "${2:-}" ]] || { echo "❌ ERROR: $1 necesita un valor (ver --help)."; exit 1; }
            [[ "$1" == "--ligand" ]] && LIGAND_OPT="$2" || RECEPTOR_OPT="$2"
            shift ;;
        -h|--help) usage; exit 0 ;;
        -*)        echo "❌ ERROR: opción desconocida '$1' (ver --help)."; exit 1 ;;
        *)         BASE="$1" ;;
    esac
    shift
done

# ── Localizar GROMACS: gmx (compilación thread-MPI) o gmx_mpi (MPI real) ────
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
[[ -z "$GMX" ]] && { echo "❌ ERROR: no encuentro ni 'gmx' ni 'gmx_mpi'."; exit 1; }
echo "✓ GROMACS: $(command -v "$GMX")"
PYTHON=$(command -v python3 || command -v python || true)
[[ -n "$PYTHON" ]] || { echo "❌ ERROR: no encuentro python3 (se necesita para leer topol.top)."; exit 1; }

# ── Nombres: los archivos de Analyses/ llevan el nombre de la carpeta actual
#    (p.ej. rmsd_gromacs_200ns.xvg); las copias a las carpetas compartidas
#    usan el nombre del sistema (la carpeta que contiene gromacs/) para que
#    los sistemas no se pisen.
NAME="${PWD##*/}"
SYS="$NAME"
[[ "$NAME" == "gromacs" ]] && SYS="$(basename "$(dirname "$PWD")")"

# ── Detectar trayectoria: argumento > dynamic<N>ns.* más reciente > step5_production
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
    echo "❌ ERROR: no encuentro trayectoria (dynamic<N>ns.{tpr,xtc} o step5_production.*)."
    exit 1
fi
[[ -f topol.top ]] || { echo "❌ ERROR: no encuentro topol.top. Ejecútalo dentro de gromacs/."; exit 1; }

# ── Duración de la simulación (para nombres de archivo y espaciado de snapshots)
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

# ── 0. Índice: grupos por defecto de GROMACS + Receptor, Receptor_BB, Ligand_<nombre>
# El ayudante lee [ molecules ] de topol.top (siguiendo los #include), agrega
# los grupos al índice e imprime líneas CLAVE='valor' (se evalúan abajo).
echo q | "$GMX" make_ndx -f "$TPR" -o "$NDX" >/dev/null 2>&1
[[ -s "$NDX" ]] || { echo "❌ ERROR: gmx make_ndx no pudo leer ${TPR}."; exit 1; }
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
if not molecules: fail("no pude leer la sección [ molecules ] de topol.top")
missing = [n for n, _ in molecules if n not in moltypes]
if missing: fail("moleculetype(s) que no están en los .itp incluidos: " + ", ".join(missing))

def kind(name):
    atoms = moltypes[name]; resn = [a[0].upper() for a in atoms]
    if name.upper() in WATER or set(resn) <= WATER: return "water"
    if name.upper() in IONS or (len(atoms) == 1 and abs(sum(a[2] for a in atoms)) > 0.5): return "ion"
    if re.match(r"^PRO[A-Z]$", name) or (resn and sum(r in AA for r in resn) / len(resn) > 0.5): return "protein"
    return "other"

# Rangos globales de átomos en el orden de la topología (desde 1)
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
    if not lig: fail("--ligand %s no coincide con ninguna molécula de topol.top. Moléculas: %s"
                     % (",".join(lig_opt), ", ".join(b["name"] for b in blocks)))
else:
    cands = [b for b in blocks if b["kind"] == "other"]
    lig = [b for b in cands if b["count"] <= 5]
    bulk = [b for b in cands if b["count"] > 5]
    if bulk:
        notes.append("las moléculas con >5 copias no se analizan como ligandos (usa --ligand para incluirlas): "
                     + ", ".join("%s(x%d)" % (b["name"], b["count"]) for b in bulk))

lig_ids = {id(b) for b in lig}
if rec_opt:
    rec = [b for b in blocks if match(b, rec_opt) and id(b) not in lig_ids]
    if not rec: fail("--receptor %s no coincide con ninguna molécula (o es el ligando)." % ",".join(rec_opt))
else:
    rec = [b for b in blocks if b["kind"] == "protein" and id(b) not in lig_ids]
    if not rec: fail("no hay proteína para usar como receptor. Usa --receptor NOMBRE.")
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

# Un grupo por tipo de molécula de ligando (todas sus copias juntas)
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
    if copies > 1: notes.append("el ligando %s tiene %d copias: se analizan juntas" % (L["name"], copies))

with open(ndx_out, "a") as fh:
    for name, idx in groups: fh.write(group(name, idx))

def show(b): return b["name"] + ("(x%d)" % b["count"] if b["count"] > 1 else "")
nres = sum(len({a[3] for a in moltypes[b["name"]]}) * b["count"] for b in rec)
print("REC_SHOW=" + q("%s (%d residuos)" % (", ".join(show(b) for b in rec), nres)))
print("LIG_SHOW=" + q(" · ".join("%s%s (%d átomos%s)" % (L["name"], "(x%d)" % L["copies"] if L["copies"] > 1 else "",
                                L["natoms"], "" if L["polar"] else ", sin N/O → contactos") for L in ligs)))
print("IGN_SHOW=" + q(", ".join(show(b) for b in ignored)))
print("LIGS=" + q(" ".join("%s:%d" % (L["tag"], L["polar"]) for L in ligs)))
print("NOTES=" + q(" | ".join(notes)))
PYDETECT
)
if [[ "$DETECT" != *"LIGS="* && "$DETECT" != *"FAIL="* ]]; then
    echo "❌ ERROR: no pude leer topol.top (falló el ayudante de python)."; exit 1
fi
FAIL=""; eval "$DETECT"
[[ -n "$FAIL" ]] && { echo "❌ ERROR: ${FAIL}"; exit 1; }

read -ra LIG_LIST <<< "$LIGS"
NLIG=${#LIG_LIST[@]}
echo " Receptor  : ${REC_SHOW}"
echo " Ligandos  : ${LIG_SHOW:-ninguno (apo): RMSD sobre el backbone del receptor, puentes de H dentro del receptor}"
[[ -n "$IGN_SHOW" ]] && echo " Ignorados : ${IGN_SHOW}"
[[ -n "$NOTES" ]] && echo " ℹ ${NOTES}"

# Nombre de salida de un análisis: el nombre del ligando solo se agrega si hay varios
out(){ # $1 = análisis, $2 = etiqueta del ligando (puede ir vacía)
    if (( NLIG > 1 )); then echo "${OUT}/${1}_${NAME}_${2}_${LABEL}.xvg"
    else echo "${OUT}/${1}_${NAME}_${LABEL}.xvg"; fi
}

# ── 1. Trayectoria centrada: centrar Receptor, escribir System ──────────────
if [[ ! -f "$CENTERED" ]]; then
    printf "Receptor\nSystem\n" | "$GMX" trjconv -s "$TPR" -f "$XTC" -n "$NDX" \
        -o "$CENTERED" -pbc mol -center || exit 1
fi

# ── 2-3. Por ligando: RMSD, luego puentes de H o contactos ──────────────────
# gmx hbond solo usa N y O como donadores/aceptores. Los ligandos sin átomos
# de N ni O (polietileno, polipropileno, poliestireno, PVC, PTFE...) no tienen
# ninguno y gmx hbond se niega a correr. Para ellos se calculan en su lugar
# los contactos Receptor–ligando (pares de átomos < 0.4 nm) y la distancia
# mínima.
hbonds(){ # $1 = segundo grupo, $2 = salida (sintaxis nueva -r/-t, respaldo antiguo)
    "$GMX" hbond -f "$CENTERED" -s "$TPR" -n "$NDX" -num "$2" -tu ns \
        -r Receptor -t "$1" 2>/dev/null \
    || printf "Receptor\n%s\n" "$1" | "$GMX" hbond -f "$CENTERED" \
        -s "$TPR" -n "$NDX" -num "$2" -tu ns \
    || echo "⚠ Falló el análisis de puentes de H ($1); se continúa con lo demás."
}
contacts(){ # $1 = grupo del ligando, $2 = etiqueta del ligando
    printf "Receptor\n%s\n" "$1" | "$GMX" mindist -s "$TPR" -f "$CENTERED" -n "$NDX" \
        -d 0.4 -tu ns -on "$(out contacts "$2")" -od "$(out mindist "$2")" \
    || echo "⚠ Falló el análisis de contactos ($1); se continúa con lo demás."
}

if (( NLIG == 0 )); then
    printf "Receptor_BB\nReceptor_BB\n" | "$GMX" rms -s "$TPR" -f "$CENTERED" \
        -n "$NDX" -o "$(out rmsd "")" -tu ns || exit 1
    hbonds Receptor "$(out hbond "")"
fi
for entry in "${LIG_LIST[@]}"; do
    tag="${entry%:*}"; polar="${entry##*:}"; grp="Ligand_${tag}"
    (( NLIG > 1 )) && echo "── Ligando ${tag} ──"
    printf "Receptor_BB\n%s\n" "$grp" | "$GMX" rms -s "$TPR" -f "$CENTERED" \
        -n "$NDX" -o "$(out rmsd "$tag")" -tu ns || exit 1
    if (( polar )); then
        hbonds "$grp" "$(out hbond "$tag")"
    else
        echo "✓ ${tag} no tiene átomos de N ni O (p.ej. microplástico poliolefínico): se omiten"
        echo "  los puentes de H y se calculan contactos Receptor–ligando (< 0.4 nm)."
        contacts "$grp" "$tag"
    fi
done

# ── 4. RMSF por residuo (Receptor) ──────────────────────────────────────────
printf "Receptor\n" | "$GMX" rmsf -s "$TPR" -f "$CENTERED" -n "$NDX" \
    -o "${OUT}/rmsf_${NAME}_${LABEL}.xvg" -res || exit 1

# ── 5. Radio de giro (Receptor) ─────────────────────────────────────────────
printf "Receptor\n" | "$GMX" gyrate -s "$TPR" -f "$CENTERED" -n "$NDX" \
    -o "${OUT}/gyrate_${NAME}_${LABEL}.xvg" || exit 1

# ── 6. Snapshots PDB en t = 0, mitad y final (System) ───────────────────────
SNAP_DT=$(( TOTAL_NS * 500 )); (( SNAP_DT < 1 )) && SNAP_DT=1
printf "System\n" | "$GMX" trjconv -s "$TPR" -f "$CENTERED" -n "$NDX" \
    -dt "$SNAP_DT" -o "${OUT}/trj_${NAME}_${LABEL}.pdb" -pbc mol || exit 1

# ── Opcional: copiar resultados a las carpetas compartidas del proyecto, si existen
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

echo "═══ LISTO: rmsd, hbond/contacts, rmsf, gyrate (.xvg) y trj (.pdb) en ${OUT}/ ═══"
