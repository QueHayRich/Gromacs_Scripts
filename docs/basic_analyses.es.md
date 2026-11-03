[![English](https://img.shields.io/badge/lang-English-blue.svg)](basic_analyses.md)
[![Español](https://img.shields.io/badge/lang-Español-green.svg)](basic_analyses.es.md)

# Basic Analyses — documentación detallada

[← Volver a la portada](../README.es.md) ·
[← Autogromacs](Autogromacs.es.md) ·
[Free Energy →](free_energy.es.md)

**Análisis estándar de las trayectorias de Autogromacs, eligiendo los grupos
de GROMACS por nombre en lugar de por número.**

`basic_analyses_with_names.es.sh` es el compañero de
[`Autogromacs.es.sh`](Autogromacs.es.md). Copia este único archivo a la carpeta
`gromacs/` de CHARMM-GUI cuando termine la producción y córrelo. Genera los
análisis habituales de estabilidad estructural y de interacción (trayectoria
centrada, RMSD, puentes de hidrógeno o contactos, RMSF, radio de giro y
snapshots PDB), sin menús interactivos y sin nada que editar.

---

## Contenido

1. [Por qué "con nombres"](#1-por-qué-con-nombres)
2. [Requisitos](#2-requisitos)
3. [Inicio rápido](#3-inicio-rápido)
4. [Cómo se detecta la trayectoria](#4-cómo-se-detecta-la-trayectoria)
5. [Receptor, ligandos y el archivo de índice](#5-receptor-ligandos-y-el-archivo-de-índice)
6. [Los análisis, uno por uno](#6-los-análisis-uno-por-uno)
7. [Casos especiales](#7-casos-especiales)
8. [Archivos de salida](#8-archivos-de-salida)
9. [Copia a las carpetas compartidas del proyecto](#9-copia-a-las-carpetas-compartidas-del-proyecto)
10. [Graficar los resultados](#10-graficar-los-resultados)
11. [Solución de problemas](#11-solución-de-problemas)

---

## 1. Por qué "con nombres"

Una forma común de automatizar GROMACS es enviar **números** de grupo a las
herramientas de GROMACS, por ejemplo `echo 4 13 | gmx rms ...`. Esos números
dependen de lo que contiene el sistema (ligando, iones, lípidos...) y cambian
de un sistema a otro, así que el mismo script podía analizar en silencio el
grupo equivocado.

Esta versión construye su propio archivo de índice y siempre pasa los grupos
**por nombre** (`Receptor`, `Receptor_BB`, `Ligand_<nombre>`, `System`), así
que el mismo archivo da el resultado correcto en cualquier sistema, sean cuales
sean la proteína y el ligando.

| Análisis | Con números | Este script (nombres) |
|---|---|---|
| Trayectoria centrada | `1 0` | `Receptor`, `System` |
| RMSD | `4 13` | `Receptor_BB`, `Ligand_<nombre>` |
| Puentes de H | `1 13` | `Receptor`, `Ligand_<nombre>` |
| RMSF | `1` | `Receptor` |
| Radio de giro | `1` | `Receptor` |
| Snapshots | `0` | `System` |

## 2. Requisitos

- **GROMACS**, ya sea `gmx` (compilación thread-MPI) o `gmx_mpi`
  (compilación MPI). Se toma del `PATH` o se carga desde un `GMXRC` estándar
  (`/usr/local/gromacs/bin`, `/usr/local/gromacs/*/bin`, `/opt/gromacs/bin`,
  `~/gromacs/bin`). No se necesita GPU.
- **Python 3** (solo la biblioteca estándar), para leer `topol.top`.
- Una **producción terminada** en la carpeta actual:
  `dynamic<N>ns.{tpr,xtc}` (salida de Autogromacs) o
  `step5_production*.{tpr,xtc}` (corrida de CHARMM-GUI hecha a mano).
- El **`topol.top`** de CHARMM-GUI y los `.itp` que incluye; de ahí se
  detectan el receptor y los ligandos.
- Un **sistema de proteína soluble** (CHARMM-GUI Solution Builder / Ligand
  Reader). Este repositorio no está pensado para sistemas con membrana.

## 3. Inicio rápido

```bash
cd /ruta/a/mi_sistema/gromacs
cp /ruta/a/basic_analyses_with_names.es.sh .

bash basic_analyses_with_names.es.sh                 # detecta la trayectoria solo
bash basic_analyses_with_names.es.sh dynamic150ns    # o indica el nombre base
bash basic_analyses_with_names.es.sh --help          # todas las opciones
```

Para trayectorias largas puedes dejarlo corriendo en segundo plano:

```bash
nohup bash basic_analyses_with_names.es.sh > analyses.log 2>&1 &
```

Al arrancar, el script imprime lo que encontró, por ejemplo:

```
✓ GROMACS: /usr/local/gromacs/bin/gmx
═══ ComplejoA · dynamic200ns.xtc · 200 ns → Analyses/ ═══
 Receptor  : PROA (128 residuos)
 Ligandos  : LIG (104 átomos, sin N/O → contactos)
 Ignorados : POT(x23), CLA(x23), TIP3(x8179)
```

Lee estas tres líneas antes de usar los resultados: dicen exactamente qué se
tomó como receptor, qué como ligando y qué quedó fuera.

## 4. Cómo se detecta la trayectoria

En este orden:

1. **Nombre base dado como argumento** (`dynamic150ns` → `dynamic150ns.tpr`
   + `dynamic150ns.xtc`).
2. El `dynamic*ns.tpr` **más reciente** que tenga su `.xtc` correspondiente.
3. `step5_production.xtc` con el `step5_production*.tpr` más reciente.

Si no encuentra ninguno, o falta `topol.top`, el script se detiene con un
mensaje claro.

**Duración de la simulación.** Se toma del nombre `dynamic<N>ns`. Con otros
nombres se lee del `.tpr` (`gmx dump`, `nsteps × dt`). La duración se usa en
los nombres de los archivos de salida (`_200ns`) y para ubicar los
snapshots.

## 5. Receptor, ligandos y el archivo de índice

El receptor y los ligandos se detectan a partir de la **sección
`[ molecules ]` de `topol.top`** (siguiendo sus `#include`), igual que lo hace
[Free Energy](free_energy.es.md), así que las dos herramientas siempre
coinciden en qué es el ligando. Cada tipo de molécula se clasifica así:

| Tipo | Regla | Se usa como |
|---|---|---|
| Proteína | se llama `PROA`, `PROB`, ... o es mayoritariamente aminoácidos | **Receptor** |
| Agua | `TIP3`, `SOL`, `HOH`, ... | se ignora |
| Iones | `POT`, `CLA`, `SOD`, ... o un solo átomo con carga | se ignoran |
| Cualquier cosa con más de 5 copias | lípidos, cosolventes, crowders | se ignora (con aviso) |
| Todo lo demás | moléculas pequeñas, oligómeros de microplásticos, cofactores... | **un ligando por tipo de molécula** |

**Cada ligando se analiza por separado.** Si un sistema tiene, por ejemplo,
un microplástico y un ligando endógeno, o un ligando y un cofactor, cada uno
recibe su propio grupo (`Ligand_<nombre>`) y sus propios archivos de RMSD y
de puentes de H o contactos, así que sus resultados nunca se mezclan. Las
copias de una misma molécula (p. ej. tres oligómeros idénticos) se analizan
juntas en un solo grupo.

El script nunca modifica el `index.ndx` de CHARMM-GUI. Escribe su propio
`Analyses/index_analysis.ndx` con los grupos por defecto de GROMACS
(`System`, `Protein`, `Backbone`, ...) más `Receptor`, `Receptor_BB` (N, CA, C
del receptor) y un `Ligand_<nombre>` por ligando.

**Elegir a mano.** Solo hace falta en casos especiales; los nombres son los
de `[ molecules ]` (o un nombre de residuo):

```bash
bash basic_analyses_with_names.es.sh --ligand LIG            # dejar fuera un cofactor (solo LIG)
bash basic_analyses_with_names.es.sh --ligand PROB --receptor PROA   # cadena peptídica como ligando
bash basic_analyses_with_names.es.sh --ligand LIG,HEME       # varios, cada uno por separado
```

`--ligand` también permite incluir una molécula con más de 5 copias (p. ej.
muchos oligómeros de microplástico en la misma caja).

## 6. Los análisis, uno por uno

Todos los ejes de tiempo están en **ns**. Los análisis 2–6 usan la
trayectoria centrada.

### 6.1 Trayectoria centrada

```
gmx trjconv -pbc mol -center    (centrar: Receptor · salida: System)
```

Recompone las moléculas que cruzan los límites periódicos y pone la
proteína en el centro de la caja, para que el ligando no parezca "saltar" de
un lado a otro de la caja. Es la entrada de todos los demás análisis y se
calcula **una sola vez**: si `Analyses/<base>_centered.xtc` ya existe, se
reutiliza.

### 6.2 RMSD del ligando

```
gmx rms    (ajuste: Receptor_BB · cálculo: Ligand_<nombre>)    uno por ligando
```

Cada frame se superpone sobre el esqueleto de la proteína de la estructura
de referencia y luego se mide el RMSD del ligando. Indica si el ligando
**se mantiene en su pose de unión** respecto a la proteína: una curva baja y
plana significa una pose estable; un salto significa que el ligando se movió
a otro sitio o salió del bolsillo.

Sin ligando, se calcula en cambio el RMSD del esqueleto (estabilidad global
de la proteína).

### 6.3 Puentes de hidrógeno receptor–ligando

```
gmx hbond -r Receptor -t Ligand_<nombre>     (sintaxis nueva, GROMACS 2024+), uno por ligando
```

Número de puentes de hidrógeno entre proteína y ligando en cada frame. Si la
sintaxis nueva no está disponible, el script usa la selección interactiva
antigua. Un número persistente de puentes de H respalda una interacción
específica; cero durante casi toda la corrida sugiere una unión débil o
inespecífica.

### 6.4 RMSF por residuo

```
gmx rmsf -res    (Receptor)
```

Fluctuación promedio de cada residuo alrededor de su posición media. Los
picos suelen ser asas y extremos; los cambios en los residuos del sitio de
unión al comparar complejos (o complejo contra proteína sola) señalan
residuos afectados por el ligando.

### 6.5 Radio de giro

```
gmx gyrate    (Receptor)
```

Compactación de la proteína a lo largo del tiempo. Un valor estable
significa que se conserva el plegamiento; un aumento sostenido indica
desplegamiento o una apertura grande de dominios.

### 6.6 Snapshots PDB

```
gmx trjconv -dt <total/2>    (System)
```

Escribe tres frames (inicio, mitad y final de la corrida) en un solo archivo
PDB de varios modelos. Ábrelo en PyMOL, VMD o ChimeraX para ver cómo cambió
el complejo y para preparar figuras.

## 7. Casos especiales

**Sin ligando (proteína sola).** El RMSD se calcula sobre el backbone del
receptor y los puentes de H dentro del receptor. El script lo indica:
`Ligandos  : ninguno (apo): RMSD sobre el backbone del receptor, puentes de H dentro del receptor`

**Varios ligandos.** Cada uno se analiza por separado y su nombre se agrega a
los archivos de RMSD y de puentes de H o contactos (ver
[Archivos de salida](#8-archivos-de-salida)). RMSF, radio de giro y snapshots
se calculan una sola vez, para el receptor y el sistema completo.

**Péptido o proteína como ligando.** Para GROMACS un péptido es proteína, así
que por defecto se toma como parte del receptor. Usa
`--ligand PROB --receptor PROA` (con los nombres de tus cadenas).

**Ligandos sin átomos de N ni O** (polietileno, polipropileno, poliestireno,
PVC, PTFE y otros microplásticos similares). `gmx hbond` solo usa N y O como
donadores/aceptores, así que no tiene nada que contar y se niega a correr.
El script lee los átomos del ligando en la estructura, ve que ninguno es N
ni O y calcula
en su lugar **contactos receptor–ligando**:

```
gmx mindist -d 0.4    (Receptor · Ligand_<nombre>)
```

- `contacts_<nombre>_<N>ns.xvg`: número de pares de átomos receptor–ligando
  a menos de 0.4 nm en cada frame (cuánto toca el ligando a la proteína).
- `mindist_<nombre>_<N>ns.xvg`: distancia mínima receptor–ligando en cada
  frame (si el ligando sigue unido o se aleja).

**Si falla el análisis de puentes de H o de contactos, el script no se
detiene.** Muestra un aviso y continúa con RMSF, radio de giro y snapshots.

## 8. Archivos de salida

Todo va a `Analyses/` dentro de la carpeta `gromacs/`:

```
Analyses/
├── index_analysis.ndx            # grupos por defecto + Receptor, Receptor_BB, Ligand_<lig>
├── dynamic<N>ns_centered.xtc     # trayectoria centrada (se reutiliza)
├── rmsd_<nombre>_<N>ns.xvg       # RMSD del ligando (del esqueleto si no hay ligando)
├── hbond_<nombre>_<N>ns.xvg      # puentes de H receptor–ligando
│   o contacts_<nombre>_<N>ns.xvg #   contactos < 0.4 nm   (ligandos sin N/O)
│     mindist_<nombre>_<N>ns.xvg  #   distancia mínima     (ligandos sin N/O)
├── rmsf_<nombre>_<N>ns.xvg       # RMSF por residuo
├── gyrate_<nombre>_<N>ns.xvg     # radio de giro
└── trj_<nombre>_<N>ns.pdb        # snapshots: inicio / mitad / final
```

`<nombre>` es el nombre de la carpeta actual (normalmente `gromacs`) y `<N>`
la duración de la simulación en ns. Con **más de un ligando**, los archivos de
RMSD y de puentes de H o contactos llevan también el nombre del ligando, p.
ej. `rmsd_gromacs_LIG_200ns.xvg` y `rmsd_gromacs_HEME_200ns.xvg`; con un
ligando los nombres son exactamente los de arriba. Si vuelves a correr el script, se
sobrescriben los `.xvg` y el `.pdb`, pero se reutiliza la trayectoria
centrada.

## 9. Copia a las carpetas compartidas del proyecto

Muchos proyectos guardan los resultados de todos los
sistemas lado a lado:

```
Proyecto/
├── rmsd/  hbond/  contacts/  rmsf/  gyrate/  conformation/
├── ComplejoA/gromacs/     ← aquí se corre el script
└── ComplejoB/gromacs/
```

Si alguna de estas carpetas existe **dos niveles arriba** de `gromacs/`
(`../../rmsd`, `../../hbond`, `../../contacts`, `../../rmsf`,
`../../gyrate`, `../../conformation`), cada resultado se copia también ahí,
renombrado con el **nombre del sistema** (la carpeta que contiene a
`gromacs/`):

```
../../rmsd/rmsd_ComplejoA_200ns.xvg
../../conformation/trj_ComplejoA_200ns.pdb
../../rmsd/rmsd_ComplejoA_LIG_200ns.xvg      (con varios ligandos)
```

así los resultados de distintos sistemas nunca se sobrescriben. Las carpetas
que no existen simplemente se omiten; no se crea ninguna.

## 10. Graficar los resultados

Los `.xvg` son archivos de texto (las líneas que empiezan con `#` o `@` son
comentarios). Algunas opciones rápidas:

```bash
xmgrace Analyses/rmsd_gromacs_200ns.xvg        # Grace, si está instalado
```

```python
import numpy as np, matplotlib.pyplot as plt
t, y = np.loadtxt("Analyses/rmsd_gromacs_200ns.xvg", comments=("#", "@"), unpack=True)
plt.plot(t, y); plt.xlabel("Tiempo (ns)"); plt.ylabel("RMSD (nm)"); plt.show()
```

También se pueden importar directo a Excel, Origin o GraphPad saltando las
líneas de comentario.

## 11. Solución de problemas

**`no encuentro ni 'gmx' ni 'gmx_mpi'`** — carga GROMACS primero
(`source /ruta/a/gromacs/bin/GMXRC`).

**`no encuentro trayectoria`** — no hay `dynamic<N>ns.{tpr,xtc}` ni
`step5_production.*` en la carpeta, o el `.xtc` no coincide con el nombre
del `.tpr`. Indica el nombre base explícitamente.

**`no encuentro topol.top`** — corre el script dentro de la carpeta
`gromacs/` de CHARMM-GUI.

**`no encuentro python3`** — instala Python 3 (cualquier versión; solo se
usa la biblioteca estándar).

**Los ligandos no son los que esperabas** (un cofactor que no querías, un
ligando con muchas copias que quedó fuera, un péptido dentro del receptor) —
revisa las líneas `Receptor / Ligandos / Ignorados` del inicio y elige a mano
con `--ligand` y `--receptor` (ver
[sección 5](#5-receptor-ligandos-y-el-archivo-de-índice)).

**`no hay proteína para usar como receptor`** — la carpeta no tiene proteína
(p. ej. una corrida de CHARMM-GUI del ligando solo). Usa `--receptor NOMBRE`
si el receptor no es una proteína.

**Los snapshots se ven rotos a través de la caja** — usa la trayectoria
centrada (`Analyses/*_centered.xtc`) para cualquier análisis adicional, no
el `.xtc` original.

---

**Autores:** Ricardo Aguilar-Garay (ricardoaguilargaray@gmail.com) ·
contribuidora: Dulce Estefanía Nicolás Álvarez, Departamento de Fisiología,
ENCB-IPN (dnicolas@ipn.mx). Licencia MIT. Por favor cita el artículo de
referencia indicado en la [portada](../README.es.md#cómo-citar).
