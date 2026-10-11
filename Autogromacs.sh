#!/bin/bash
###############################################################################
# Autogromacs.sh — pipeline GROMACS automático con duración ADAPTABLE (v2)
#
# Copia este ÚNICO archivo a cualquier carpeta gromacs/ de CHARMM-GUI y ejecuta:
#       bash Autogromacs.sh            (pregunta los ns la primera vez)
#       bash Autogromacs.sh 200        (o pásalos como argumento; útil con nohup)
#
# Corre TODO el flujo: minimización -> equilibración -> producción.
#
# La PRODUCCIÓN se hace en TRAMOS de 100 ns: corre 100 ns, luego extiende
# +100 ns (convert-tpr -extend + mdrun -cpi -append, método validado en
# Gelsolina) hasta completar el total. Con -append todo queda en UN SOLO
# archivo continuo, sin frames duplicados.
#
# REANUDABLE: si la corrida se corta, vuelve a ejecutar el script (sin
# argumentos) y continúa desde el último checkpoint. El objetivo total se
# recuerda en un archivo de estado: NO hace falta volver a teclearlo.
#
# UNIVERSAL: detecta gmx (PATH o GMXRC), detecta si hay GPU NVIDIA y ajusta
# los flags de mdrun; usa todos los núcleos disponibles (nproc). Si el modo
# GPU agresivo falla (p.ej. CUDA error #700), reintenta automáticamente en
# modo estable (bonded/update en CPU) desde el último checkpoint.
# Compatible con GROMACS 2024 y 2026: -append solo se usa si ya hay .cpt
# (2026 aborta con -append sin checkpoint) y se desactivan funciones
# experimentales (CUDA graphs) aunque estén activas en el entorno.
#
# Al terminar deja la trayectoria con nombre adaptable:
#       dynamic<N>ns.{xtc,tpr,edr,log,gro,cpt}      (p.ej. dynamic200ns.xtc)
#
# Proyecto: Maydelis Viltres-Acosta — MPs <-> Daphnia magna
###############################################################################
set -uo pipefail

# ── 0a. Localizar GROMACS ───────────────────────────────────────────────────
# Primero el PATH; si no, la instalación estándar del laboratorio; si no,
# rutas alternativas comunes (por si el script se distribuye a otros equipos).
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
    echo "❌ ERROR: no encuentro 'gmx' en el PATH ni un GMXRC en las rutas conocidas."
    echo "   Verifica la instalación de GROMACS en esta máquina."
    exit 1
fi
GMX_VER=$(gmx --version 2>/dev/null | awk -F: '/GROMACS version/{gsub(/ /,"",$2); print $2; exit}')
echo "✓ GROMACS ${GMX_VER:-?}: $(command -v gmx)"

# Sin funciones experimentales en producción (CUDA graphs puede venir activado
# desde ~/.bashrc en algunas máquinas). Solo afecta a este script.
if [[ -n "${GMX_CUDA_GRAPH:-}" ]]; then
    echo "ℹ GMX_CUDA_GRAPH estaba activo en el entorno: se desactiva (experimental)."
    unset GMX_CUDA_GRAPH
fi

# ── 0b. Verificar GPU (OBLIGATORIA) y núcleos disponibles ───────────────────
NTOMP=$(nproc 2>/dev/null || echo 4)
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
    echo "✓ GPU detectada: $(nvidia-smi -L | head -n1)"
else
    echo "❌ ERROR: este pipeline requiere GPU NVIDIA y no se detectó ninguna."
    echo "   Verifica drivers (nvidia-smi) antes de correr en esta máquina."
    exit 1
fi

# Flags de mdrun (mismo esquema para 2080 / 3080 / 5090; GROMACS se adapta):
# HARD = TODO en GPU incluido -update: máximo rendimiento (validado en Gelsolina)
# SAFE = bonded/update en CPU: rescate automático para CUDA #700 / inestabilidad
#        (-bonded/-update cpu explícitos: en 2026 'auto' deja ambos en GPU)
FLAGS_HARD="-pme gpu -nb gpu -bonded gpu -update gpu -pmefft gpu -ntmpi 1 -ntomp ${NTOMP} -pin on"
FLAGS_SAFE="-pme gpu -nb gpu -bonded cpu -update cpu -pmefft gpu -ntmpi 1 -ntomp ${NTOMP} -pin on"

# ── Parámetros fijos ────────────────────────────────────────────────────────
CHUNK_NS=100                     # tamaño de cada tramo (ns)
PROD_MDP="step5_production.mdp"
PROD="step5_production"          # deffnm base (se renombra a dynamic<N>ns al final)
STATE=".autogromacs_total_ns"    # aquí se persiste el objetivo TOTAL (para reanudar)

# Verificar insumos mínimos de CHARMM-GUI antes de arrancar
for f in "$PROD_MDP" topol.top index.ndx; do
    if [[ ! -f "$f" ]]; then
        echo "❌ ERROR: falta '$f'. ¿Estás dentro de la carpeta gromacs/ de CHARMM-GUI?"
        exit 1
    fi
done

# ── 0c. Objetivo total de ns: argumento > estado guardado > preguntar ───────
TOTAL_NS="${1:-}"
if [[ -z "$TOTAL_NS" && -f "$STATE" ]]; then
    TOTAL_NS=$(cat "$STATE")
    echo "🔄 Reanudación detectada: objetivo guardado = ${TOTAL_NS} ns (de ${STATE})"
fi
if [[ -z "$TOTAL_NS" ]]; then
    if [[ -t 0 ]]; then
        read -rp "¿Cuántos ns quieres simular en producción? " TOTAL_NS
    else
        echo "❌ Sin terminal interactiva y sin objetivo previo: pasa los ns como argumento."
        echo "   Ejemplo:  nohup bash Autogromacs.sh 200 > autogromacs.log 2>&1 &"
        exit 1
    fi
fi
if ! [[ "$TOTAL_NS" =~ ^[0-9]+$ ]] || (( TOTAL_NS < 1 )); then
    echo "❌ Debes ingresar un entero de ns >= 1 (ej. 100, 200, 250)."; exit 1
fi
echo "$TOTAL_NS" > "$STATE"

# dt desde el .mdp (fallback 0.002 ps)
DT=$(awk -F= '/^[[:space:]]*dt[[:space:]]*=/{gsub(/[^0-9.]/,"",$2); print $2; exit}' "$PROD_MDP" 2>/dev/null)
DT=${DT:-0.002}
ns_to_steps(){ awk -v ns="$1" -v dt="$DT" 'BEGIN{printf "%.0f", ns*1000/dt}'; }

# Leer el objetivo REAL del tpr vigente (fuente de verdad; evita desfases
# si un corte ocurre justo entre extender el tpr y actualizar cualquier registro)
tpr_target_ns(){
    local nsteps
    nsteps=$(gmx dump -s "${PROD}.tpr" 2>/dev/null \
             | grep -m1 -oE 'nsteps[^0-9-]*[0-9]+' | grep -oE '[0-9]+$')
    [[ -n "$nsteps" ]] || { echo 0; return; }
    awk -v s="$nsteps" -v dt="$DT" 'BEGIN{printf "%.0f", s*dt/1000}'
}

FINAL="dynamic${TOTAL_NS}ns"

echo "═══════════════════════════════════════════════════════════════════════"
echo " Autogromacs — objetivo: ${TOTAL_NS} ns   (dt=${DT} ps · tramos de ${CHUNK_NS} ns)"
echo " Hardware: GPU + ${NTOMP} hilos CPU"
echo " Salida final: ${FINAL}.{xtc,tpr,...}"
echo "═══════════════════════════════════════════════════════════════════════"

# Si ya existe la salida final, no rehacer (evita sobrescribir resultados)
if [[ -f ${FINAL}.xtc ]]; then
    echo "⚠ Ya existe ${FINAL}.xtc en esta carpeta. Bórralo o elige otro tiempo para rehacer."
    exit 0
fi

# ── 1. Minimización (omitir si ya está hecha) ───────────────────────────────
if [[ ! -f step4.0_minimization.gro ]]; then
    echo "=== MINIMIZACIÓN ==="
    gmx grompp -f step4.0_minimization.mdp -o step4.0_minimization.tpr \
        -c step3_input.gro -r step3_input.gro -p topol.top -n index.ndx -maxwarn 1 || exit 1
    gmx mdrun -v -deffnm step4.0_minimization -ntmpi 1 -ntomp "${NTOMP}" || exit 1
else
    echo "✓ Minimización ya hecha, se omite."
fi

# ── 2. Equilibración (omitir si ya está hecha) ──────────────────────────────
if [[ ! -f step4.1_equilibration.gro ]]; then
    echo "=== EQUILIBRACIÓN ==="
    gmx grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr \
        -c step4.0_minimization.gro -r step3_input.gro -p topol.top -n index.ndx -maxwarn 1 || exit 1
    gmx mdrun -v -deffnm step4.1_equilibration $FLAGS_SAFE || exit 1
else
    echo "✓ Equilibración ya hecha, se omite."
fi

# ── 3. Producción por tramos de ${CHUNK_NS} ns ──────────────────────────────
# Preparar el tpr del primer tramo si aún no existe
if [[ ! -f ${PROD}.tpr ]]; then
    FIRST_NS=$(( TOTAL_NS < CHUNK_NS ? TOTAL_NS : CHUNK_NS ))
    echo "=== PRODUCCIÓN: preparando primer tramo de ${FIRST_NS} ns ==="
    # continuation=yes en el .mdp: tomar velocidades/estado exactos del .cpt de
    # la equilibración (más preciso que el .gro) si existe
    EQ_CPT=""
    [[ -f step4.1_equilibration.cpt ]] && EQ_CPT="-t step4.1_equilibration.cpt"
    gmx grompp -f "$PROD_MDP" -o ${PROD}.tpr \
        -c step4.1_equilibration.gro $EQ_CPT -p topol.top -n index.ndx -maxwarn 1 || exit 1
    # Fijar nsteps del primer tramo escribiendo a un temporal y luego mv
    gmx convert-tpr -s ${PROD}.tpr -nsteps "$(ns_to_steps "$FIRST_NS")" -o ${PROD}_tmp.tpr || exit 1
    mv -f ${PROD}_tmp.tpr ${PROD}.tpr
fi

# Objetivo vigente leído SIEMPRE del propio tpr (a prueba de cortes)
TARGET_NS=$(tpr_target_ns)
if (( TARGET_NS == 0 )); then
    echo "❌ ERROR: no pude leer nsteps de ${PROD}.tpr (¿archivo corrupto?)."; exit 1
fi
if (( TARGET_NS > TOTAL_NS )); then
    echo "⚠ El tpr actual ya apunta a ${TARGET_NS} ns, más que el objetivo (${TOTAL_NS} ns)."
    echo "  Se completará hasta ${TARGET_NS} ns y se renombrará como dynamic${TARGET_NS}ns.*"
    TOTAL_NS=$TARGET_NS
    FINAL="dynamic${TOTAL_NS}ns"
    echo "$TOTAL_NS" > "$STATE"
fi

# Monitor de GPU en segundo plano (opcional, no altera la simulación)
GPU_PID=""
nvidia-smi dmon -s puct -d 10 >> gpu_monitor.log 2>&1 &
GPU_PID=$!
cleanup(){ [[ -n "$GPU_PID" ]] && kill "$GPU_PID" 2>/dev/null || true; }
trap cleanup EXIT

# Flags de continuación: -cpi/-append SOLO si ya existe el .cpt. GROMACS 2026
# aborta con "-append" sin checkpoint (otras versiones pueden tolerarlo). Se evalúa en cada
# intento para que el reintento SAFE continúe desde el .cpt que dejó el HARD.
cpt_flags(){
    if [[ -f ${PROD}.cpt ]]; then
        echo "-cpi ${PROD}.cpt -append"
    fi
}

# Interrupciones (Ctrl-C, kill, apagado): mdrun se detiene ordenadamente,
# guarda checkpoint y sale con código 1 — igual que un fallo real. Para NO
# confundirlo con un fallo de GPU (y saltar a SAFE) se detecta por:
#   a) la señal llegando al propio script (Ctrl-C en la terminal), o
#   b) el mensaje "Received the ... signal" que mdrun escribe en su log
#      (cuando se mata solo el proceso mdrun, p.ej. con nohup + kill).
INTERRUPTED=0
trap 'INTERRUPTED=1' INT TERM HUP
# Solo se miran las líneas del log escritas en ESTE intento (LOG_START), para
# no confundir un aviso de señal de una corrida anterior con la actual.
log_lines(){ local n; n=$(wc -l < ${PROD}.log 2>/dev/null); echo "${n:-0}"; }
log_says_signal(){
    local start=$LOG_START
    (( $(log_lines) < start )) && start=0     # mdrun creó un log nuevo
    awk -v a="$start" 'NR>a && /Received the .* signal/{f=1} END{exit !f}' \
        ${PROD}.log 2>/dev/null
}
stop_if_interrupted(){
    if (( INTERRUPTED )) || log_says_signal; then
        echo "⏹ mdrun interrumpido por señal; el checkpoint quedó guardado."
        echo "   Re-ejecuta 'bash Autogromacs.sh' (sin argumentos) para reanudar."
        exit 130
    fi
}

# Un tramo de mdrun con reintento automático en modo estable.
run_chunk(){
    echo "── mdrun [HARD: todo en GPU, máximo rendimiento] ──"
    LOG_START=$(log_lines)
    gmx mdrun -v -deffnm ${PROD} -s ${PROD}.tpr $(cpt_flags) $FLAGS_HARD && return 0
    stop_if_interrupted
    echo "⚠ El modo HARD falló (¿CUDA #700 / inestabilidad?)."
    echo "🚑 Reintentando en MODO SAFE (bonded/update en CPU) desde el checkpoint..."
    LOG_START=$(log_lines)
    gmx mdrun -v -deffnm ${PROD} -s ${PROD}.tpr $(cpt_flags) $FLAGS_SAFE
}

# Bucle de tramos: corre/continúa hasta TARGET_NS, luego extiende, hasta TOTAL_NS
while :; do
    echo "── Producción: corriendo/continuando hasta ${TARGET_NS} ns ──"
    if ! run_chunk; then
        stop_if_interrupted
        echo "❌ mdrun falló incluso en modo estable. Revisa ${PROD}.log"
        echo "   (busca LINCS warnings o energía positiva: el sistema pudo explotar)."
        exit 1
    fi

    if (( TARGET_NS >= TOTAL_NS )); then
        echo "✓ Producción completa: ${TARGET_NS} ns alcanzados."
        break
    fi

    REMAIN=$(( TOTAL_NS - TARGET_NS ))
    NEXT_NS=$(( REMAIN < CHUNK_NS ? REMAIN : CHUNK_NS ))
    EXT_PS=$(( NEXT_NS * 1000 ))
    echo "=== Extendiendo +${NEXT_NS} ns (${EXT_PS} ps) ==="
    # Escribir a temporal y luego mv: nunca dejamos un tpr a medio escribir
    gmx convert-tpr -s ${PROD}.tpr -extend "$EXT_PS" -o ${PROD}_tmp.tpr || exit 1
    mv -f ${PROD}_tmp.tpr ${PROD}.tpr
    TARGET_NS=$(tpr_target_ns)
done

cleanup; trap - EXIT

# ── 4. Renombrar salida final a dynamic<N>ns.* ──────────────────────────────
echo "=== Renombrando salida a ${FINAL}.* ==="
for ext in xtc trr tpr edr log gro cpt; do
    [[ -f ${PROD}.${ext} ]] && mv -f ${PROD}.${ext} ${FINAL}.${ext} && echo "  ✓ ${PROD}.${ext} → ${FINAL}.${ext}"
done
[[ -f ${PROD}_prev.cpt ]] && mv -f ${PROD}_prev.cpt ${FINAL}_prev.cpt
rm -f "$STATE"

echo "═══════════════════════════════════════════════════════════════════════"
echo " LISTO. Trayectoria de ${TOTAL_NS} ns:  ${FINAL}.xtc  /  ${FINAL}.tpr"
echo " Siguiente paso: bash basic_analyses.sh   (auto-detecta ${FINAL})"
echo "═══════════════════════════════════════════════════════════════════════"
