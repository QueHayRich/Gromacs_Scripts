[![English](https://img.shields.io/badge/lang-English-blue.svg)](basic_analyses.md)
[![Español](https://img.shields.io/badge/lang-Español-green.svg)](basic_analyses.es.md)

# Basic Analyses — documentación detallada

[← Volver a la portada](../README.es.md) ·
[← Autogromacs](Autogromacs.es.md) ·
[Free Energy →](free_energy.es.md)

**Análisis estándar de las trayectorias de Autogromacs, eligiendo los grupos
de GROMACS por nombre en lugar de por número.**

`basic_analyses_with_names.sh` es el compañero de
[`Autogromacs.sh`](Autogromacs.es.md). Copia este único archivo a la carpeta
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
5. [El archivo de índice y el grupo Ligand](#5-el-archivo-de-índice-y-el-grupo-ligand)
6. [Los análisis, uno por uno](#6-los-análisis-uno-por-uno)
7. [Casos especiales](#7-casos-especiales)
8. [Archivos de salida](#8-archivos-de-salida)
9. [Copia a las carpetas compartidas del proyecto](#9-copia-a-las-carpetas-compartidas-del-proyecto)
10. [Graficar los resultados](#10-graficar-los-resultados)
11. [Solución de problemas](#11-solución-de-problemas)

---

## 1. Por qué "con nombres"

El script original del laboratorio enviaba **números** de grupo a las
herramientas de GROMACS, por ejemplo `echo 4 13 | gmx rms ...`. Esos números
dependen de lo que contiene el sistema (ligando, iones, lípidos...) y cambian
de un sistema a otro, así que el mismo script podía analizar en silencio el
grupo equivocado.

Esta versión construye su propio archivo de índice y siempre pasa los grupos
**por nombre** (`Protein`, `Backbone`, `Ligand`, `System`), así que el mismo
archivo da el resultado correcto en cualquier sistema.

| Análisis | Script anterior (números) | Este script (nombres) |
|---|---|---|
| Trayectoria centrada | `1 0` | `Protein`, `System` |
| RMSD | `4 13` | `Backbone`, `Ligand` |
| Puentes de H | `1 13` | `Protein`, `Ligand` |
| RMSF | `1` | `Protein` |
| Radio de giro | `1` | `Protein` |
| Snapshots | `0` | `System` |

## 2. Requisitos

- **GROMACS**, ya sea `gmx` (compilación thread-MPI) o `gmx_mpi`
  (compilación MPI). Se toma del `PATH` o se carga desde un `GMXRC` estándar
  (`/usr/local/gromacs/bin`, `/usr/local/gromacs/*/bin`, `/opt/gromacs/bin`,
  `~/gromacs/bin`). No se necesita GPU.
- Una **producción terminada** en la carpeta actual:
  `dynamic<N>ns.{tpr,xtc}` (salida de Autogromacs) o
  `step5_production*.{tpr,xtc}` (corrida de CHARMM-GUI hecha a mano).
- El **`index.ndx`** de CHARMM-GUI, que aporta el grupo `SOLU`.

## 3. Inicio rápido

```bash
cd /ruta/a/mi_sistema/gromacs
cp /ruta/a/basic_analyses_with_names.sh .

bash basic_analyses_with_names.sh                 # detecta la trayectoria solo
bash basic_analyses_with_names.sh dynamic150ns    # o indica el nombre base
```

Para trayectorias largas puedes dejarlo corriendo en segundo plano:

```bash
nohup bash basic_analyses_with_names.sh > analyses.log 2>&1 &
```

Al arrancar, el script imprime lo que encontró, por ejemplo:

```
✓ GROMACS: /usr/local/gromacs/bin/gmx
═══ ComplexA · dynamic200ns.xtc · 200 ns → Analyses/ ═══
✓ Ligand group created (SOLU minus Protein).
```

## 4. Cómo se detecta la trayectoria

En este orden:

1. **Nombre base dado como argumento** (`dynamic150ns` → `dynamic150ns.tpr`
   + `dynamic150ns.xtc`).
2. El `dynamic*ns.tpr` **más reciente** que tenga su `.xtc` correspondiente.
3. `step5_production.xtc` con el `step5_production*.tpr` más reciente.

Si no encuentra ninguno, o falta `index.ndx`, el script se detiene con un
mensaje claro.

**Duración de la simulación.** Se toma del nombre `dynamic<N>ns`. Con otros
nombres se lee del `.tpr` (`gmx dump`, `nsteps × dt`). La duración se usa en
los nombres de los archivos de salida (`_200ns`) y para ubicar los
snapshots.

## 5. El archivo de índice y el grupo Ligand

El script nunca modifica el `index.ndx` de CHARMM-GUI. Crea el suyo propio,
`Analyses/index_analysis.ndx`:

1. `gmx make_ndx` escribe los grupos por defecto de GROMACS (`System`,
   `Protein`, `Backbone`, `C-alpha`, ...).
2. Se copia el grupo `SOLU` del `index.ndx` de CHARMM-GUI. `SOLU` contiene
   todo lo que no es agua ni iones: proteína más ligando(s).
3. Se define un grupo nuevo, `"SOLU" & ! "Protein"`, con el nombre
   **`Ligand`**.

Así, `Ligand` es **todo lo que hay en `SOLU` que no es proteína**, se llame
como se llame el residuo (una molécula pequeña, un oligómero de un
microplástico, un cofactor, varias copias de un ligando...). Si ese grupo
queda vacío, el sistema no tiene ligando y el script pasa a los análisis
solo de proteína (ver [Casos especiales](#7-casos-especiales)).

## 6. Los análisis, uno por uno

Todos los ejes de tiempo están en **ns**. Los análisis 2–6 usan la
trayectoria centrada.

### 6.1 Trayectoria centrada

```
gmx trjconv -pbc mol -center    (centrar: Protein · salida: System)
```

Recompone las moléculas que cruzan los límites periódicos y pone la
proteína en el centro de la caja, para que el ligando no parezca "saltar" de
un lado a otro de la caja. Es la entrada de todos los demás análisis y se
calcula **una sola vez**: si `Analyses/<base>_centered.xtc` ya existe, se
reutiliza.

### 6.2 RMSD del ligando

```
gmx rms    (ajuste: Backbone · cálculo: Ligand)
```

Cada frame se superpone sobre el esqueleto de la proteína de la estructura
de referencia y luego se mide el RMSD del ligando. Indica si el ligando
**se mantiene en su pose de unión** respecto a la proteína: una curva baja y
plana significa una pose estable; un salto significa que el ligando se movió
a otro sitio o salió del bolsillo.

Sin ligando, se calcula en cambio el RMSD del esqueleto (estabilidad global
de la proteína).

### 6.3 Puentes de hidrógeno proteína–ligando

```
gmx hbond -r Protein -t Ligand     (sintaxis nueva, GROMACS 2024+)
```

Número de puentes de hidrógeno entre proteína y ligando en cada frame. Si la
sintaxis nueva no está disponible, el script usa la selección interactiva
antigua. Un número persistente de puentes de H respalda una interacción
específica; cero durante casi toda la corrida sugiere una unión débil o
inespecífica.

### 6.4 RMSF por residuo

```
gmx rmsf -res    (Protein)
```

Fluctuación promedio de cada residuo alrededor de su posición media. Los
picos suelen ser asas y extremos; los cambios en los residuos del sitio de
unión al comparar complejos (o complejo contra proteína sola) señalan
residuos afectados por el ligando.

### 6.5 Radio de giro

```
gmx gyrate    (Protein)
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

**Sin ligando (proteína sola).** El grupo `Ligand` queda vacío, así que el
RMSD se calcula sobre el `Backbone` y los puentes de H dentro de la
`Protein`. El script lo indica:
`✓ No ligand in SOLU: RMSD on Backbone, H-bonds within Protein.`

**Ligandos sin átomos de N, O, S o F** (polietileno, polipropileno,
poliestireno y otros microplásticos de solo C/H). No tienen donadores ni
aceptores de puentes de H, y `gmx hbond` se niega a correr. El script lee
los átomos del ligando en la estructura, ve que ninguno es N/O/S/F y calcula
en su lugar **contactos proteína–ligando**:

```
gmx mindist -d 0.4    (Protein · Ligand)
```

- `contacts_<nombre>_<N>ns.xvg`: número de pares de átomos proteína–ligando
  a menos de 0.4 nm en cada frame (cuánto toca el ligando a la proteína).
- `mindist_<nombre>_<N>ns.xvg`: distancia mínima proteína–ligando en cada
  frame (si el ligando sigue unido o se aleja).

**Si falla el análisis de puentes de H o de contactos, el script no se
detiene.** Muestra un aviso y continúa con RMSF, radio de giro y snapshots.

## 8. Archivos de salida

Todo va a `Analyses/` dentro de la carpeta `gromacs/`:

```
Analyses/
├── index_analysis.ndx            # grupos por defecto + SOLU + Ligand
├── dynamic<N>ns_centered.xtc     # trayectoria centrada (se reutiliza)
├── rmsd_<nombre>_<N>ns.xvg       # RMSD del ligando (del esqueleto si no hay ligando)
├── hbond_<nombre>_<N>ns.xvg      # puentes de H proteína–ligando
│   o contacts_<nombre>_<N>ns.xvg #   contactos < 0.4 nm   (ligandos de solo C/H)
│     mindist_<nombre>_<N>ns.xvg  #   distancia mínima     (ligandos de solo C/H)
├── rmsf_<nombre>_<N>ns.xvg       # RMSF por residuo
├── gyrate_<nombre>_<N>ns.xvg     # radio de giro
└── trj_<nombre>_<N>ns.pdb        # snapshots: inicio / mitad / final
```

`<nombre>` es el nombre de la carpeta actual (normalmente `gromacs`) y `<N>`
la duración de la simulación en ns. Si vuelves a correr el script, se
sobrescriben los `.xvg` y el `.pdb`, pero se reutiliza la trayectoria
centrada.

## 9. Copia a las carpetas compartidas del proyecto

Muchos proyectos del laboratorio guardan los resultados de todos los
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

**`neither 'gmx' nor 'gmx_mpi' found`** — carga GROMACS primero
(`source /ruta/a/gromacs/bin/GMXRC`).

**`no trajectory found`** — no hay `dynamic<N>ns.{tpr,xtc}` ni
`step5_production.*` en la carpeta, o el `.xtc` no coincide con el nombre
del `.tpr`. Indica el nombre base explícitamente.

**`index.ndx not found`** — corre el script dentro de la carpeta `gromacs/`
de CHARMM-GUI.

**El grupo Ligand contiene algo inesperado** (p. ej. un cofactor que no
querías) — el grupo es "todo lo de SOLU que no es proteína". Revisa
`Analyses/index_analysis.ndx`, o edita ahí el grupo `Ligand` y corre a mano
el comando `gmx` que necesites.

**Los snapshots se ven rotos a través de la caja** — usa la trayectoria
centrada (`Analyses/*_centered.xtc`) para cualquier análisis adicional, no
el `.xtc` original.

---

**Autores:** Ricardo Aguilar-Garay (ricardoaguilargaray@gmail.com) ·
contribuidora: Dulce Estefanía Nicolás Álvarez, Departamento de Fisiología,
ENCB-IPN (dnicolas@ipn.mx). Licencia MIT. Por favor cita el artículo de
referencia indicado en la [portada](../README.es.md#cómo-citar).
