[![English](https://img.shields.io/badge/lang-English-blue.svg)](Autogromacs.md)
[![Español](https://img.shields.io/badge/lang-Español-green.svg)](Autogromacs.es.md)

# Autogromacs — documentación detallada

[← Volver a la portada](../README.es.md) ·
[Basic Analyses →](basic_analyses.es.md)

**Pipeline automático, reanudable y por tramos de dinámica molecular en
GROMACS para sistemas de CHARMM-GUI en GPUs NVIDIA.**

`Autogromacs.es.sh` es un solo script de Bash, autosuficiente. Lo copias a una
carpeta `gromacs/` de CHARMM-GUI recién descargada y corre todo el flujo
(minimización de energía → equilibración → producción) por el tiempo que le
pidas. Funciona sin cambios en cualquier máquina y con cualquier sistema
proteína–ligando, y se recupera solo de fallos, cortes de luz y errores de la GPU.

---

## Contenido

1. [Por qué Autogromacs](#1-por-qué-autogromacs)
2. [Requisitos](#2-requisitos)
3. [Archivos de entrada esperados](#3-archivos-de-entrada-esperados)
4. [Inicio rápido](#4-inicio-rápido)
5. [Qué hace el script, paso a paso](#5-qué-hace-el-script-paso-a-paso)
6. [La producción por tramos en detalle](#6-la-producción-por-tramos-en-detalle)
7. [Reanudar después de una interrupción](#7-reanudar-después-de-una-interrupción)
8. [Modos de GPU y rescate automático](#8-modos-de-gpu-y-rescate-automático)
9. [Archivos de salida](#9-archivos-de-salida)
10. [Mensajes y códigos de salida](#10-mensajes-y-códigos-de-salida)
11. [Solución de problemas](#11-solución-de-problemas)
12. [Preguntas frecuentes](#12-preguntas-frecuentes)
13. [Limitaciones](#13-limitaciones)

---

## 1. Por qué Autogromacs

Correr a mano un sistema de CHARMM-GUI implica teclear una cadena de
comandos `grompp` y `mdrun`, calcular `nsteps` para la duración deseada,
elegir los flags de GPU de cada máquina y, en corridas largas, extender y
unir trayectorias. Hacer todo eso a mano invita a errores, y un corte de luz
a media corrida de 500 ns puede costar días.

Autogromacs automatiza todo esto:

| Problema | Qué hace Autogromacs |
|---|---|
| Editar los `.mdp` para cambiar la duración | Le das los ns; el script calcula `nsteps` a partir del `dt` del `.mdp`. |
| Las corridas largas hay que partirlas y unirlas con `trjcat` | La producción corre en tramos de 100 ns agregados a **un solo `.xtc` continuo**. |
| Un fallo o corte de luz pierde la corrida | Vuelve a correr el script sin argumentos; continúa desde el último checkpoint. |
| Los flags de GPU cambian entre máquinas | Los mismos flags sirven en RTX 2080 / 3080 / 5090; GROMACS se adapta. |
| CUDA error #700 / inestabilidad de la GPU | Reintento automático desde el checkpoint en modo estable. |
| Un Ctrl-C confundido con un fallo de GPU | Las interrupciones se detectan y el script sale limpio (código 130). |
| Distintas versiones de GROMACS | Funciona con GROMACS 2024 y 2026 (probado en 2026.4). |

## 2. Requisitos

| Componente | Requisito |
|---|---|
| Sistema operativo | Linux (probado en Ubuntu) |
| GROMACS | 2024 o 2026, compilado con GPU (CUDA) |
| GPU | GPU NVIDIA con el driver funcionando (`nvidia-smi`): **obligatoria** |
| Sistema de entrada | Salida de CHARMM-GUI **Solution Builder** (carpeta `gromacs/`) |

**Cómo se localiza GROMACS.** El script busca primero `gmx` en el `PATH`. Si
no está, carga el primer `GMXRC` que encuentre en:

```
/usr/local/gromacs/bin/GMXRC
/usr/local/gromacs/*/bin/GMXRC
/opt/gromacs/bin/GMXRC
~/software/gromacs-*/bin/GMXRC
~/gromacs/bin/GMXRC
```

Al arrancar imprime la versión y la ruta detectadas
(`✓ GROMACS 2026.4: /usr/local/gromacs/bin/gmx`).

> **Máquinas con RTX 5090 (Blackwell):** GROMACS debe estar compilado con
> CUDA ≥ 12.8. El script no necesita cambios.

## 3. Archivos de entrada esperados

Son los nombres estándar de una carpeta `gromacs/` de CHARMM-GUI Solution
Builder:

```
step3_input.gro              # sistema solvatado y neutralizado
topol.top                    # topología (más la carpeta toppar/)
index.ndx                    # grupos de índice (SOLU, MEMB, SOLV, ...)
step4.0_minimization.mdp
step4.1_equilibration.mdp
step5_production.mdp
```

Antes de empezar, el script verifica que existan `step5_production.mdp`,
`topol.top` e `index.ndx`, y se detiene con un mensaje claro si falta
alguno; normalmente significa que no se corrió dentro de la carpeta
`gromacs/`.

## 4. Inicio rápido

```bash
cd /ruta/a/mi_sistema/gromacs
cp /ruta/a/Autogromacs.es.sh .

# Interactivo: pregunta cuántos ns de producción quieres
bash Autogromacs.es.sh

# No interactivo (recomendado para corridas largas; sobrevive a cerrar la terminal)
nohup bash Autogromacs.es.sh 200 > autogromacs.log 2>&1 &

# Seguir el avance
tail -f autogromacs.log

# Después de un fallo, corte de luz o Ctrl-C: vuelve a correrlo sin argumentos
bash Autogromacs.es.sh
```

La duración debe ser un número entero de ns ≥ 1 (p. ej. `100`, `200`,
`250`). Con `nohup` no hay terminal donde preguntar, así que en la primera
corrida la duración **debe** pasarse como argumento.

## 5. Qué hace el script, paso a paso

```
 ┌──────────────────────────┐
 │ 0. Revisión del entorno  │  gmx, GPU, hilos de CPU, archivos de entrada
 └──────────┬───────────────┘
            ▼
 ┌──────────────────────────┐
 │ Registro del objetivo    │  argumento > estado guardado > pregunta
 └──────────┬───────────────┘
            ▼
 ┌──────────────────────────┐
 │ 1. Minimización          │  se omite si existe step4.0_minimization.gro
 └──────────┬───────────────┘
            ▼
 ┌──────────────────────────┐
 │ 2. Equilibración         │  se omite si existe step4.1_equilibration.gro
 └──────────┬───────────────┘
            ▼
 ┌──────────────────────────┐      ┌────────────────────────────┐
 │ 3. Tramo de producción   │ ───► │ ¿se alcanzó el total?      │
 │    (mdrun, ≤100 ns)      │      │  no → convert-tpr -extend  │──┐
 └──────────▲───────────────┘      │  sí → terminar             │  │
            └──────────────────────┴────────────────────────────┘◄─┘
            ▼
 ┌──────────────────────────┐
 │ 4. Renombrar a           │  dynamic<N>ns.{xtc,tpr,edr,log,gro,cpt}
 │    dynamic<N>ns.*        │
 └──────────────────────────┘
```

### 0. Revisión del entorno

- Localiza GROMACS (ver [Requisitos](#2-requisitos)) e imprime su versión.
- Si `GMX_CUDA_GRAPH` está activo en el entorno (por ejemplo desde
  `~/.bashrc`), se **desactiva** para esta corrida, porque GROMACS marca los
  CUDA graphs como experimentales. Solo afecta a este script.
- Verifica que haya una GPU NVIDIA (`nvidia-smi -L`) y se detiene si no.
- Fija el número de hilos OpenMP con `nproc`.

### Registro del objetivo

La duración total de la producción se toma, en este orden, de:

1. el argumento de la línea de comandos (`bash Autogromacs.es.sh 200`);
2. el archivo oculto de estado `.autogromacs_total_ns`, escrito en la
   primera corrida;
3. una pregunta interactiva.

El valor se guarda en `.autogromacs_total_ns`, así que al reanudar nunca te
lo vuelve a preguntar y no puedes teclear por error un número distinto.

El paso de tiempo `dt` se lee de `step5_production.mdp` (si no, 0.002 ps) y
se usa en todas las conversiones ns ↔ pasos.

Si ya existe `dynamic<N>ns.xtc` en la carpeta, el script se detiene sin
tocar nada: nunca sobrescribe resultados terminados.

### 1. Minimización

```
gmx grompp -f step4.0_minimization.mdp -c step3_input.gro -r step3_input.gro ...
gmx mdrun  -deffnm step4.0_minimization
```

Protocolo estándar de CHARMM-GUI, con restricciones de posición referidas a
`step3_input.gro`. Se omite si ya existe `step4.0_minimization.gro`.

### 2. Equilibración

```
gmx grompp -f step4.1_equilibration.mdp -c step4.0_minimization.gro -r step3_input.gro ...
gmx mdrun  -deffnm step4.1_equilibration   (modo GPU SAFE)
```

Corre en la GPU en la configuración estable (ver
[Modos de GPU](#8-modos-de-gpu-y-rescate-automático)). Se omite si ya existe
`step4.1_equilibration.gro`.

### 3. Producción

Se describe en detalle en la [siguiente sección](#6-la-producción-por-tramos-en-detalle).

### 4. Cierre

Todos los archivos `step5_production.*` se renombran a `dynamic<N>ns.*`, se
borra el archivo de estado y el script indica el siguiente paso: correr
[Basic Analyses](basic_analyses.es.md).

## 6. La producción por tramos en detalle

La producción corre en **tramos de 100 ns** (`CHUNK_NS=100`). Para una
corrida de 250 ns:

| Tramo | Objetivo del `.tpr` | Qué pasa |
|---|---|---|
| 1 | 100 ns | `grompp` + `convert-tpr -nsteps` → `mdrun` |
| 2 | 200 ns | `convert-tpr -extend 100000` → `mdrun -cpi -append` |
| 3 | 250 ns | `convert-tpr -extend 50000` → `mdrun -cpi -append` |

**Primer tramo.** `grompp` construye `step5_production.tpr` a partir de la
estructura equilibrada. Si existe `step4.1_equilibration.cpt`, se pasa con
`-t`, de modo que la corrida arranca del estado equilibrado exacto
(velocidades incluidas), como exige `continuation = yes` en el `.mdp` de
CHARMM-GUI. Luego `gmx convert-tpr -nsteps` fija la duración de ese primer
tramo en min(total, 100) ns.

**Tramos siguientes.** Al terminar un tramo, si aún no se alcanza el total,
`gmx convert-tpr -extend` agrega al `.tpr` el siguiente tramo (100 ns, o lo
que falte) y `mdrun -cpi -append` continúa desde el checkpoint, escribiendo
en los **mismos** `.xtc`, `.edr` y `.log`. El resultado es una sola
trayectoria continua, sin frames duplicados y sin necesidad de `gmx trjcat`.

**Por qué `-cpi/-append` solo cuando existe un checkpoint.** GROMACS 2026
aborta si se da `-append` sin checkpoint. El script agrega estos flags solo
cuando existe `step5_production.cpt`, así que el primer tramo arranca bien
en cualquier versión.

**El `.tpr` es la fuente de verdad.** Después de cada tramo, el script lee
el objetivo vigente del propio `.tpr` (`gmx dump`, `nsteps × dt`). Aunque
haya un corte de luz justo entre extender el `.tpr` y cualquier otra cosa,
el script y la simulación nunca se desfasan.

**Escritura segura.** `convert-tpr` siempre escribe a un archivo temporal
(`step5_production_tmp.tpr`) que después se mueve a su lugar, así que una
extensión interrumpida no puede dejar un `.tpr` a medio escribir.

## 7. Reanudar después de una interrupción

Solo vuelve a correr el script en la misma carpeta, sin argumentos:

```bash
bash Autogromacs.es.sh
# o
nohup bash Autogromacs.es.sh > autogromacs.log 2>&1 &
```

Qué pasa depende de dónde se cortó:

| Se interrumpió durante | Al reanudar |
|---|---|
| Minimización | La minimización se repite (es corta). |
| Equilibración | Se omite la minimización; la equilibración se repite. |
| Producción (cualquier tramo) | Se omiten minimización y equilibración; `mdrun` continúa desde `step5_production.cpt`. |
| Entre dos tramos | El `.tpr` ya tiene el objetivo correcto; la corrida continúa. |
| Después de renombrar | Existe `dynamic<N>ns.xtc`; el script lo avisa y sale. |

**Si cambias de opinión sobre la duración.** Puedes alargar una corrida
dando un número mayor al reanudar (`bash Autogromacs.es.sh 300`). Si el `.tpr`
ya apunta a más de lo que indicas, el script termina hasta el objetivo del
`.tpr` y nombra la salida de acuerdo con él.

## 8. Modos de GPU y rescate automático

Se usan dos juegos de flags de `mdrun`:

| Modo | Flags | Se usa para |
|---|---|---|
| **HARD** (máximo rendimiento) | `-nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu -ntmpi 1 -ntomp <nproc> -pin on` | Primer intento de cada tramo de producción |
| **SAFE** (estable) | `-nb gpu -pme gpu -pmefft gpu -bonded cpu -update cpu -ntmpi 1 -ntomp <nproc> -pin on` | Equilibración, y reintento automático cuando HARD falla |

En modo SAFE, `-bonded cpu -update cpu` se ponen **explícitamente**, porque
las versiones recientes de GROMACS dejan `update` en la GPU cuando se deja en
`auto`.

**Rescate automático.** Si un tramo falla en modo HARD (por ejemplo CUDA
error #700 u otra inestabilidad de la GPU), el script avisa y repite el mismo
tramo en modo SAFE **desde el último checkpoint**, sin que tengas que hacer
nada. Si SAFE también falla, el script se detiene y te remite a
`step5_production.log`: lo más probable es que el sistema mismo sea
inestable (ver [Solución de problemas](#11-solución-de-problemas)).

**Las interrupciones no son fallos.** Cuando presionas Ctrl-C, haces `kill`
al proceso o apagas la máquina, `mdrun` se detiene ordenadamente, guarda un
checkpoint y sale con código de error, igual que en un fallo real. Para no
saltar por error al modo SAFE, el script detecta las interrupciones de dos
formas:

1. una trampa de señales (`INT`, `TERM`, `HUP`) en el propio script;
2. la línea `Received the ... signal` que `mdrun` escribe en su log. Solo se
   revisan las líneas escritas en el intento actual, así que un mensaje
   viejo de una corrida anterior no causa una falsa alarma.

En ambos casos imprime `⏹ mdrun interrumpido por señal` y sale con código
**130**. Vuelve a correr el script para continuar.

**Telemetría de la GPU.** Durante la producción, `nvidia-smi dmon -s puct -d 10`
registra potencia, uso, frecuencias y temperatura cada 10 s en
`gpu_monitor.log`. Se detiene solo al terminar el script.

## 9. Archivos de salida

Al final, en la carpeta `gromacs/`:

```
dynamic<N>ns.xtc        # trayectoria de producción, una sola y continua
dynamic<N>ns.tpr        # archivo de entrada de la corrida completa (úsalo para analizar)
dynamic<N>ns.edr        # energías
dynamic<N>ns.log        # log de mdrun (todos los tramos, agregados)
dynamic<N>ns.gro        # coordenadas finales
dynamic<N>ns.cpt        # checkpoint final (permite extender la corrida después)
dynamic<N>ns_prev.cpt   # checkpoint anterior (respaldo)
gpu_monitor.log         # telemetría de la GPU
step4.0_minimization.*  # archivos de la minimización
step4.1_equilibration.* # archivos de la equilibración
```

Mientras la corrida está en curso, los archivos de producción se llaman
`step5_production.*` y el archivo oculto `.autogromacs_total_ns` guarda el
objetivo.

## 10. Mensajes y códigos de salida

| Mensaje | Significado |
|---|---|
| `🔄 Reanudación detectada: objetivo guardado = N ns` | Reanudación: el objetivo se leyó del archivo de estado. |
| `✓ Minimización ya hecha, se omite.` | Etapa omitida porque su salida ya existe. |
| `⚠ El modo HARD falló ... Reintentando en MODO SAFE` | Fallo de GPU; reintento automático en modo estable. |
| `⏹ mdrun interrumpido por señal` | Ctrl-C / kill / apagado; checkpoint guardado. |
| `⚠ Ya existe dynamic<N>ns.xtc` | La corrida ya había terminado; no se hace nada. |
| `❌ mdrun falló incluso en modo estable` | Fallaron ambos modos; revisa el log. |

| Código de salida | Significado |
|---|---|
| `0` | Terminó, o la trayectoria final ya existía. |
| `1` | Error (falta un archivo, no hay GPU, falla de `grompp`/`mdrun`). |
| `130` | Interrumpido por una señal; vuelve a correr el script para reanudar. |

## 11. Solución de problemas

**`no encuentro 'gmx'`** — GROMACS no está en el `PATH` ni en una ubicación
estándar. Corre `source /ruta/a/gromacs/bin/GMXRC` antes del script, o
agrégalo a `~/.bashrc`.

**`este pipeline requiere GPU NVIDIA`** — `nvidia-smi` no funciona. Revisa el
driver (`nvidia-smi` en una terminal) y reinicia si se acaba de actualizar.

**`falta 'step5_production.mdp'` (o `topol.top`, `index.ndx`)** — no estás
dentro de la carpeta `gromacs/` de CHARMM-GUI.

**`Sin terminal interactiva y sin objetivo previo`** — con `nohup`, la primera
corrida necesita la duración como argumento: `nohup bash Autogromacs.es.sh 200 ...`.

**Falla también en modo SAFE** — busca en `step5_production.log`
`LINCS warnings`, `Segmentation fault` o una energía potencial grande y
positiva. Apuntan a un sistema inestable (parámetros del ligando
incorrectos, átomos encimados, paso de tiempo demasiado grande) más que a la
GPU. Revisa la minimización, la equilibración y la topología del ligando.

**El modo HARD falla en todos los tramos** — la corrida se completa igual en
modo SAFE, solo más lento. Errores CUDA #700 repetidos pueden indicar
sobrecalentamiento (ver `gpu_monitor.log`), un overclock inestable o un
problema del driver.

**Quiero empezar desde cero** — borra los archivos generados
(`step4.0_minimization.*`, `step4.1_equilibration.*`, `step5_production.*`,
`.autogromacs_total_ns`), conservando los archivos originales de CHARMM-GUI,
y vuelve a correr el script.

## 12. Preguntas frecuentes

**¿Puedo extender una corrida terminada?** Sí, directamente con GROMACS:
`gmx convert-tpr -s dynamic200ns.tpr -extend 100000 -o dynamic300ns.tpr`,
y luego `gmx mdrun -s dynamic300ns.tpr -cpi dynamic200ns.cpt -deffnm dynamic200ns -append`.
Después renombra los archivos si quieres el nombre `dynamic300ns`.

**¿Puedo correr varios sistemas a la vez?** Cada sistema corre en su propia
carpeta `gromacs/`, así que no hay conflicto de archivos. Pero en una sola
GPU, las corridas simultáneas comparten la tarjeta y cada una irá más lenta.

**¿Por qué tramos de 100 ns?** Es un equilibrio: los tramos son lo bastante
cortos para que el reintento por fallo de GPU repita poco trabajo, y lo
bastante pocos para que el costo de extender el `.tpr` sea despreciable.

**¿Modifica los archivos `.mdp`?** No. Todos los cambios de duración se hacen
sobre el `.tpr` con `convert-tpr`; los `.mdp` de CHARMM-GUI quedan como
están.

## 13. Limitaciones

- **Solo sistemas solubles.** Solo la nomenclatura de CHARMM-GUI **Solution
  Builder** (`step4.0`, `step4.1`, `step5`). Los sistemas con membrana
  (Membrane Builder, `step6.x` / `step7`) quedan fuera del alcance de este
  repositorio y no se corren.
- Requiere una GPU NVIDIA; no hay modo solo-CPU.
- Una sola GPU por corrida (`-ntmpi 1`).

---

**Autores:** Ricardo Aguilar-Garay (ricardoaguilargaray@gmail.com) ·
contribuidora: Dulce Estefanía Nicolás Álvarez, Departamento de Fisiología,
ENCB-IPN (dnicolas@ipn.mx). Licencia MIT. Por favor cita el artículo de
referencia indicado en la [portada](../README.es.md#cómo-citar).
