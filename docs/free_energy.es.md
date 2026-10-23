[![English](https://img.shields.io/badge/lang-English-blue.svg)](free_energy.md)
[![Español](https://img.shields.io/badge/lang-Español-green.svg)](free_energy.es.md)

# Free Energy — documentación detallada

[← Volver a la portada](../README.es.md) ·
[← Basic Analyses](basic_analyses.es.md)

**Energética de unión de complejos proteína–ligando a partir de las
trayectorias de Autogromacs: MM-GBSA, energía de interacción LJ/Coulomb
directa y validación de la unión, reunidas en un solo Excel para todos tus
sistemas.**

La herramienta tiene dos scripts:

| Script | Dónde se corre | Qué hace |
|---|---|---|
| `free_energy.sh` | Dentro de la carpeta `gromacs/` de cada sistema | Calcula la energética de **un** sistema y escribe `summary.csv` |
| `collect_free_energy.py` | Donde sea, una vez | Junta los `summary.csv` de **todos** los sistemas en un solo Excel |

También hay una guía rápida en texto plano:
[`LEEME_free_energy.txt`](../LEEME_free_energy.txt).

---

## Contenido

1. [La idea](#1-la-idea)
2. [Requisitos](#2-requisitos)
3. [Inicio rápido](#3-inicio-rápido)
4. [Modos y opciones](#4-modos-y-opciones)
5. [Detección automática de receptor y ligando](#5-detección-automática-de-receptor-y-ligando)
6. [El ambiente de conda de gmx_MMPBSA](#6-el-ambiente-de-conda-de-gmx_mmpbsa)
7. [Qué se calcula, paso a paso](#7-qué-se-calcula-paso-a-paso)
8. [Archivos de salida](#8-archivos-de-salida)
9. [Flags de calidad](#9-flags-de-calidad)
10. [Cómo leer y reportar los resultados](#10-cómo-leer-y-reportar-los-resultados)
11. [collect_free_energy.py: un Excel para todos los sistemas](#11-collect_free_energypy-un-excel-para-todos-los-sistemas)
12. [Correr muchos sistemas de una vez](#12-correr-muchos-sistemas-de-una-vez)
13. [Solución de problemas](#13-solución-de-problemas)

---

## 1. La idea

La misma filosofía que Autogromacs y Basic Analyses: copias **un archivo** a
la carpeta `gromacs/` de CHARMM-GUI después de la producción y lo corres:

```bash
bash free_energy.sh
```

No hay menú de `make_ndx` ni números de grupo que elegir. El script lee
`topol.top` y decide solo:

- **Receptor** = la(s) cadena(s) de proteína (`PROA`, `PROB`, ...).
- **Se ignoran** = agua (`TIP3`, ...) e iones (`POT`, `CLA`, `SOD`, ...).
- **Ligando** = lo que queda (un ácido graso, un hexámero de PET, PIP2, ...).

Antes de calcular nada, imprime lo que detectó (receptor, ligando, átomos
pesados, cargas, ventana, temperatura). Si el sistema es ambiguo, **se
detiene** y te dice qué opción agregar. Nunca adivina en silencio.

Todos los cálculos usan la **misma ventana de frames**, así que MM-GBSA,
LJ/Coulomb y validación son directamente comparables. Los valores por
defecto son el protocolo recomendado: si todos los sistemas se corren con
ellos, los resultados se pueden comparar entre sí.

## 2. Requisitos

| Componente | Se necesita para | Notas |
|---|---|---|
| GROMACS (`gmx` o `gmx_mpi`) | Todo | El del sistema; se localiza como en Autogromacs (`PATH` o `GMXRC` estándar) |
| Python 3 | Todo | Solo la biblioteca estándar para `free_energy.sh` |
| conda + un ambiente con **gmx_MMPBSA** | MM-GBSA | No se necesita en modo `lj` |
| `mpirun` | Opcional | MM-GBSA en paralelo; si falla, corre en serie |
| **openpyxl** | `collect_free_energy.py` | Sin él, se escriben archivos CSV en lugar del Excel |
| Producción terminada | Todo | `dynamic<N>ns.{tpr,xtc}` (Autogromacs) o `step5_production.*` |
| `topol.top`, `index.ndx` | Todo | Archivos de CHARMM-GUI en `gromacs/` |

No se necesita GPU: el recálculo de energías usa la CPU (los grupos de
energía no están soportados en GPU) y solo evalúa unos cientos de frames.

## 3. Inicio rápido

```bash
cd /ruta/a/mi_sistema/gromacs
cp /ruta/a/free_energy.sh .

bash free_energy.sh --check      # 1. ver qué detectó (segundos, no calcula nada)
bash free_energy.sh              # 2. calcular (MM-GBSA + LJ + validación)

# Para trabajos largos, en segundo plano:
nohup bash free_energy.sh > free_energy.log 2>&1 &
```

Después, cuando terminen todos los sistemas, desde cualquier carpeta:

```bash
python3 collect_free_energy.py /ruta/a/Proyecto_A /ruta/a/Proyecto_B -o Energias.xlsx
```

La primera vez, córrelo en una terminal (no con `nohup`) para que pueda
preguntar qué ambiente de conda tiene gmx_MMPBSA. La respuesta se recuerda
(ver la [sección 6](#6-el-ambiente-de-conda-de-gmx_mmpbsa)).

## 4. Modos y opciones

### Modos

| Comando | Qué hace |
|---|---|
| `bash free_energy.sh` | MM-GBSA + LJ/Coulomb + validación (por defecto) |
| `bash free_energy.sh --check` | Solo muestra lo detectado y escribe el índice; no calcula nada |
| `bash free_energy.sh mmgbsa` | Solo MM-GBSA (+ validación) |
| `bash free_energy.sh lj` | Solo la interacción LJ/Coulomb (+ validación); rápido, sin conda |
| `bash free_energy.sh summary` | Rehace `summary.csv` con los resultados que ya existen |
| `bash free_energy.sh --help` / `--version` | Ayuda / versión |

### Opciones (todas opcionales)

| Opción | Por defecto | Úsala cuando |
|---|---|---|
| `--start NS` | 25 % de la corrida | El RMSD muestra que el sistema se estabiliza después (p. ej. `--start 80`) |
| `--end NS` | Final de la corrida | Quieres excluir el final |
| `--frames N` | 300 | Más frames (menor SEM, más lento) o menos (más rápido); mínimo 10 |
| `--ligand NOMBRE` | Automático | La detección es ambigua: dos ligandos, un cofactor, un complejo proteína–proteína. Separados por comas para varios |
| `--receptor NOMBRE` | Todas las proteínas | El receptor no es el de por defecto |
| `--np N` | min(núcleos, 16) | Número de procesos MPI para gmx_MMPBSA |
| `--env NOMBRE` | Se pregunta una vez y se recuerda | Usar otro ambiente de conda; `--env ask` para volver a elegir |
| `--no-decomp` | Con descomposición | Corrida más rápida sin descomposición por residuo |
| `--base NOMBRE` | Automático | Indicar la trayectoria explícitamente, p. ej. `dynamic200ns` |
| `--name NOMBRE` | Automático | El nombre del sistema en las tablas sale demasiado genérico |

> **Regla de oro:** para comparar sistemas, corre **todos con las mismas
> opciones** (misma ventana y mismo número de frames). Con los valores por
> defecto esto se cumple solo.

### Situaciones y qué teclear

| Situación | Comando |
|---|---|
| Primera vez en un sistema (recomendado) | `bash free_energy.sh --check`, luego `bash free_energy.sh` |
| Caso normal (proteína + un ligando) | `bash free_energy.sh` |
| Solo la energía LJ, rápido | `bash free_energy.sh lj` |
| Más rápido, sin residuos clave | `bash free_energy.sh --no-decomp` |
| El RMSD llega a una meseta a los 80 ns | `bash free_energy.sh --start 80` (misma ventana en todos los sistemas) |
| Se detiene con "more than one candidate ligand" | `--ligand NOMBRE` (o `NOMBRE1,NOMBRE2`), como sugiere el mensaje |
| Complejo proteína–proteína | `--ligand PROB` (receptor = `PROA`) |
| Solo rehacer la tabla | `bash free_energy.sh summary` |
| El nombre del sistema sale genérico | `--name Gelsolina_LPA` |

## 5. Detección automática de receptor y ligando

El script lee `topol.top` y todos los `.itp` que incluye, y clasifica cada
molécula de la sección `[ molecules ]`:

| Clase | Regla |
|---|---|
| **Agua** | Nombres de agua conocidos (`TIP3`, `SOL`, `SPC`, `WAT`, ...) |
| **Ion** | Nombres de iones conocidos (`POT`, `CLA`, `SOD`, `MG`, `ZN`, ...) o átomos sueltos con carga |
| **Proteína** | `PROA`, `PROB`, ... o más de la mitad de sus residuos son aminoácidos |
| **Otra** | Todo lo demás: candidatos a ligando |

Luego:

- **Ligando** = las moléculas "otras" con **≤ 5 copias**. Las que tienen más
  de 5 copias (lípidos de membrana, cosolventes) se ignoran y se reportan en
  una nota.
- Si un ligando tiene varias copias (≤ 5), todas entran al grupo `Ligand`, y
  una nota lo indica.
- **Receptor** = todas las cadenas de proteína.

El script **se detiene sin calcular nada** cuando:

- hay dos candidatos a ligando distintos (los enumera y sugiere `--ligand`);
- no hay molécula pequeña pero sí varias cadenas de proteína (sugiere
  `--ligand <cadena>` para un complejo proteína–proteína);
- no hay ligando (solo proteína, agua e iones);
- una molécula de `topol.top` no está definida en los `.itp`.

También avisa (sin detenerse) si la **carga neta** del sistema no es cero
(¿faltan iones?) o si el ligando está **muy cargado** (|q| ≥ 3), porque
entonces los términos electrostáticos serán grandes y ruidosos.

Los átomos pesados se cuentan por **masa ≥ 4**, así que el conteo es
correcto incluso con reparto de masa de hidrógenos (HMR). El índice
`index_FE.ndx` conserva los grupos de CHARMM-GUI y agrega `Receptor`,
`Ligand`, `Receptor_BB` (N, CA, C), `Ligand_heavy` y `Complex`.

El **nombre del sistema** en las tablas es la primera carpeta arriba de
`gromacs/` que no tiene un nombre genérico (`gromacs`, `charmm-gui-*`,
`complex`, `system`, `md`, `run`, ...). Se cambia con `--name`.

## 6. El ambiente de conda de gmx_MMPBSA

`gmx` / `gmx_mpi` es el GROMACS del sistema (se usa para `trjconv`, el
recálculo, etc.). gmx_MMPBSA vive aparte, en un ambiente de conda.

- La **primera vez** que corres el script en una terminal, busca los
  ambientes que contienen gmx_MMPBSA y pregunta cuál usar:

  ```
   Which conda environment has gmx_MMPBSA?   (conda: /home/usuario/miniconda3)
     1) gmxMMPBSA            /home/usuario/miniconda3/envs/gmxMMPBSA
   Number or env name [Enter = gmxMMPBSA]:
  ```

  Enter acepta la sugerencia. La respuesta se guarda en
  `~/.config/free_energy/conda_env` y ya no se vuelve a preguntar (también
  se usa con `nohup`).
- Para volver a elegir: `bash free_energy.sh --env ask`, o directamente
  `--env OTRO_ENV`.
- Funciona aunque conda no esté en el `PATH` (lo típico con `nohup`): lo
  busca en `~/miniconda3`, `~/anaconda3`, `~/miniforge3`, `~/mambaforge` y
  `/opt/conda`.
- Sin terminal, si solo hay un ambiente con gmx_MMPBSA (o uno llamado
  `gmxMMPBSA`), se usa automáticamente; si no, el script se detiene y pide
  correrlo una vez en una terminal o pasar `--env`.
- Si el `gmx` propio del ambiente es demasiado viejo para leer el `.tpr`,
  gmx_MMPBSA usa automáticamente el GROMACS del sistema. Es lo esperado y no
  hay que hacer nada.

## 7. Qué se calcula, paso a paso

### [1] Trayectoria de la ventana

```
gmx trjconv -b <inicio> -e <fin> -dt <paso> -pbc mol -center   (centrar: Receptor · salida: System)
```

Del 25 % de la corrida al final, unos 300 frames (p. ej. 200 ns → 50–200 ns,
un frame cada 500 ps). El paso entre frames se redondea a un múltiplo del
intervalo con que se guardaron los frames. El resultado, `window.xtc`, lo
comparten todos los cálculos siguientes.

### [2] Validación (siempre se calcula)

| Cantidad | Comando | Para qué |
|---|---|---|
| Contactos < 4 Å y distancia mínima | `gmx mindist -d 0.4` | ¿El ligando sigue tocando al receptor? |
| RMSD del ligando | `gmx rms` (ajuste `Receptor_BB`, cálculo `Ligand_heavy`) | ¿El ligando se quedó en su sitio? |
| Puentes de hidrógeno | `gmx hbond -r Receptor -t Ligand` (con respaldo a la sintaxis antigua) | Interacciones polares específicas |

Estos números alimentan los flags de calidad (`LIGAND_UNBOUND`,
`LIGAND_DISPLACED`).

### [3] Energía de interacción LJ / Coulomb

Se copia el `.mdp` de producción con `energygrps = Receptor Ligand` y la
escritura de salida apagada, se construye un `.tpr` nuevo y se reevalúa la
ventana con `gmx mdrun -rerun` en la CPU. Después `gmx energy` extrae
`LJ-SR:Receptor-Ligand` y `Coul-SR:Receptor-Ligand`.

Se reporta (kcal/mol): `E_LJ`, `E_Coul`, `E_int = E_LJ + E_Coul` (promedio,
SD y estimación de error de GROMACS), y `E_LJ` por átomo pesado del ligando.

### [4] MM-GBSA (gmx_MMPBSA)

Parámetros:

| Parámetro | Valor |
|---|---|
| Modelo | Born generalizado, `igb=5` (OBC II) |
| Radios | `mbondi2` (`PBRadii=3`) |
| Sal | 0.15 M |
| Temperatura | Leída del `.tpr` |
| Entropía | C2 (cumulante de segundo orden) |
| Descomposición | Por residuo (`idecomp=2`), residuos a menos de 6 Å |

Corre con `mpirun -np N` cuando está disponible, y reintenta en serie si MPI
falla.

Se reporta (kcal/mol, promedio / SD / SEM): `dVDW`, `dEEL`, `dEGB`,
`dESURF`, `dGgas`, `dGsolv` y **`dH`** (la entalpía de unión MM-GBSA). A
partir de la fluctuación de `dGgas`, el script calcula el término de
entropía C2 `-TdS = σ²/(2kT)` y `dG = dH − TdS`, además de la eficiencia de
ligando (`LE = dH / átomos pesados`) y los **10 residuos que más aportan**.

### [5] Resumen

Todo se reúne en una fila (`summary.csv`) y se imprime:

```
 RESULTS (kcal/mol, mean ± SEM; LJ terms mean ± Err.Est.)   frames: 301
   dVDW     -38.12 ± 0.21    dEEL     -10.45 ± 0.40
   dEGB      17.80 ± 0.35    dESURF    -4.95 ± 0.02
   dH       -35.72 ± 0.25    -TdS(C2) 9.80   dG(C2) -25.92
   E_LJ     -40.10 ± 0.60    E_Coul   -12.30 ± 0.90    E_int -52.40 ± 1.08
   H-bonds 1.20   contacts<4A 310   ligand RMSD 0.18 nm
   Top residues: LEU45(-2.10); PHE17(-1.85); ...
   FLAGS: OK
```

*(Valores de ejemplo, solo para ilustrar.)*

**Reanudable.** Cada paso se omite si su salida ya existe. Para rehacer un
paso, borra su salida y vuelve a correr: `MMGBSA/FINAL_RESULTS_MMPBSA.dat`,
`LJ/lj_coul.xvg` o `window.xtc` (este último rehace todo lo que depende de
él).

## 8. Archivos de salida

```
Analyses/FreeEnergy_<sistema>_<N>ns/
├── detected_system.txt     lo que se detectó (para tu bitácora)
├── index_FE.ndx            grupos de CHARMM-GUI + Receptor, Ligand, Receptor_BB, Ligand_heavy, Complex
├── window.xtc              los frames que usan todos los cálculos
├── validation/             hbonds.xvg, numcont.xvg, mindist.xvg, ligand_rmsd.xvg
├── LJ/                     rerun.mdp, rerun.tpr, lj_coul.xvg, energy_stdout.txt, logs
├── MMGBSA/                 mmgbsa.in, FINAL_RESULTS_MMPBSA.dat/.csv,
│                           FINAL_DECOMP_MMPBSA.dat/.csv, gmx_MMPBSA.out
├── summary.csv             UNA fila con todos los resultados de este sistema
└── top_residues.csv        los 10 residuos que más aportan
```

Columnas principales de `summary.csv`:

| Grupo | Columnas |
|---|---|
| Identidad | `system`, `receptor`, `ligand`, `ligand_resname`, `ligand_heavy_atoms`, `ligand_charge`, `receptor_charge` |
| Ventana | `sim_ns`, `window_start_ns`, `window_end_ns`, `frame_step_ps`, `temperature_K`, `n_frames` |
| MM-GBSA | `dVDW_*`, `dEEL_*`, `dEGB_*`, `dESURF_*`, `dGgas_*`, `dGsolv_*`, `dH_*` (`_mean`, `_sd`, `_sem`), `sigma_dGgas`, `mTdS_C2`, `dG_C2`, `LE_dH`, `LE_dG`, `top5_residues` |
| LJ/Coulomb | `E_LJ_*`, `E_Coul_*`, `E_int_*` (`_mean`, `_sd`, `_err`), `E_LJ_per_heavy` |
| Validación | `hbonds_mean`, `hbonds_sd`, `contacts_4A_mean`, `mindist_mean_nm`, `frac_unbound`, `lig_rmsd_mean_nm`, `lig_rmsd_last_nm` |
| Calidad | `flags`, `notes` |
| Procedencia | `gromacs`, `gmx_MMPBSA`, `date`, `path` |

Todas las energías están en **kcal/mol**.

## 9. Flags de calidad

| Flag | Condición | Qué hacer |
|---|---|---|
| `OK` | No se encontró ningún problema | — |
| `C2_UNRELIABLE` | σ(dGgas) > 3.6 kcal/mol | Reportar `dH`, no `dG(C2)` |
| `HIGH_SEM` | SEM de `dH` > 1 kcal/mol | Más frames o una ventana más larga |
| `POSITIVE_dH` | `dH` > 0 | Unión desfavorable: revisar la trayectoria |
| `LIGAND_UNBOUND` | Ligando a > 4.5 Å del receptor en > 10 % de los frames | El promedio **no** representa al complejo |
| `LIGAND_DISPLACED` | RMSD promedio del ligando > 1 nm | El ligando se movió a otro sitio |
| `FEW_FRAMES` | Menos de 100 frames | Aumentar `--frames` o la ventana |

## 10. Cómo leer y reportar los resultados

- **`dH` (MM-GBSA)** es la cantidad principal para comparar afinidades entre
  ligandos.
- **`E_LJ`** es una comparación directa de los contactos de dispersión,
  justa entre ligandos cargados y neutros.
- **`E_Coul` y `E_int`** son solo informativos. Se calculan en vacío: en
  ligandos cargados (oleato, araquidonato, LPA, PIP2...) salen muy grandes y
  **no** deben usarse solos para afirmar que un ligando se une más fuerte.
- **`-TdS (C2)` y `dG`**: úsalos solo si no aparece el flag `C2_UNRELIABLE`.
- **`LE` y `E_LJ/HA`** están normalizados por átomo pesado; son útiles
  cuando los ligandos tienen tamaños muy distintos (p. ej. hexámeros de
  LDPE, PET y PC).
- **No** llames "energía libre de unión" a `E_LJ` ni a `E_int`: son energías
  de interacción.
- Antes de usar un número, revisa los flags. Un sistema con
  `LIGAND_UNBOUND` o `LIGAND_DISPLACED` debe revisarse (RMSD y snapshots de
  Basic Analyses) antes de interpretar su promedio.

## 11. collect_free_energy.py: un Excel para todos los sistemas

Se corre **una vez**, desde donde quieras, dándole las carpetas padre de tus
proyectos:

```bash
python3 collect_free_energy.py ~/Documents/Microplasticos \
        ~/Documents/Ligandos_Endogenos -o Energias_ronda2.xlsx
```

- Busca todos los `.../Analyses/FreeEnergy_*/summary.csv` dentro de esas
  carpetas, a cualquier profundidad. No hace falta mover ni juntar carpetas.
- El **nombre de cada carpeta padre** se vuelve la columna `Group`, así que
  se pueden comparar en la misma tabla sistemas de distintos proyectos. Usa
  `--labels A,B` para nombres de grupo propios (uno por carpeta, en orden).
- **Dónde se guarda:** en la carpeta desde donde corres el comando. Sin `-o`
  se llama `FreeEnergy_summary_<fecha>.xlsx`; con `-o` eliges nombre y ruta.

### Hojas del Excel

| Hoja | Contenido |
|---|---|
| **Summary** | Tabla lista para artículo: promedio ± SEM (kcal/mol) |
| **Data_kcal** | Todos los números (promedio, SD, SEM) en kcal/mol |
| **Data_kJ** | Lo mismo en kJ/mol |
| **Top_residues** | Residuos clave por sistema |
| **Metadata** | Ventana, frames, temperatura, versiones, notas, ruta |
| **Legend** | Qué significa cada columna y cada flag |

### openpyxl

El Excel necesita `openpyxl`. Para verificarlo:

```bash
python3 collect_free_energy.py --check-deps
```

Si falta, el script:

1. busca un ambiente de conda cuyo Python lo tenga (`gmxMMPBSA` primero) y
   se vuelve a ejecutar con ese Python;
2. si no hay, muestra una alarma con los comandos de instalación:
   ```
   python3 -m pip install --user openpyxl          # Ubuntu 22.04
   sudo apt install python3-openpyxl               # Ubuntu 23.04+
   conda activate gmxMMPBSA && pip install openpyxl
   ```
3. en una terminal, ofrece instalarlo en ese momento;
4. si no se instala, escribe las mismas tablas como archivos CSV.

`--no-install` evita la pregunta (útil en scripts o con `nohup`).
`free_energy.sh` también avisa al final si falta `openpyxl`.

## 12. Correr muchos sistemas de una vez

**Paso 0.** Una sola vez, en una terminal, elige el ambiente de conda (se
recuerda):

```bash
cd <un_sistema>/gromacs && bash free_energy.sh --check
```

**Paso 1.** Revisa la detección en todos los sistemas (segundos):

```bash
for d in $(find /ruta/Proyecto_A /ruta/Proyecto_B -name "dynamic*ns.xtc" -printf '%h\n' | sort -u); do
    cp ~/scripts/free_energy.sh "$d"/
    (cd "$d" && bash free_energy.sh --check)
done
```

**Paso 2.** Si todo se ve bien, calcula todos en segundo plano:

```bash
nohup bash -c 'for d in $(find /ruta/Proyecto_A /ruta/Proyecto_B \
    -name "dynamic*ns.xtc" -printf "%h\n" | sort -u); do
    (cd "$d" && bash free_energy.sh); done' > free_energy_all.log 2>&1 &
```

**Paso 3.** Junta los resultados:

```bash
python3 collect_free_energy.py /ruta/Proyecto_A /ruta/Proyecto_B
```

Consejo: antes de lanzar muchos sistemas, corre uno completo y revisa
`MMGBSA/gmx_MMPBSA.out`.

## 13. Solución de problemas

| Mensaje / problema | Qué hacer |
|---|---|
| `AUTO-DETECTION STOPPED: ...` | Lee el mensaje: dice qué opción agregar (`--ligand`, `--receptor`) |
| `System net charge is not zero` | Revisa `topol.top` (¿faltan iones?) |
| `conda was not found` | Corre `conda init bash` y abre una terminal nueva, o usa `bash free_energy.sh lj` |
| `conda env '...' does not contain gmx_MMPBSA` | Vuelve a elegir: `bash free_energy.sh --env ask` |
| `Cannot ask for the conda env (no terminal...)` | Córrelo una vez en una terminal, o pasa `--env NOMBRE` |
| gmx_MMPBSA falla | Revisa `MMGBSA/gmx_MMPBSA.out` y `gmx_MMPBSA.log` |
| MPI falla | El script reintenta solo en serie; o usa `--np 1` |
| `Invalid window` | `--start`/`--end` fuera de la duración de la corrida |
| Rehacer un paso | Borra su salida (`MMGBSA/FINAL_RESULTS_MMPBSA.dat`, `LJ/lj_coul.xvg` o `window.xtc`) y vuelve a correr |

**Nota de validación:** los scripts se probaron con la topología real del
complejo gastrotropina–oleato y con salidas simuladas de GROMACS y
gmx_MMPBSA.

---

**Autores:** Ricardo Aguilar-Garay (ricardoaguilargaray@gmail.com) ·
contribuidora: Dulce Estefanía Nicolás Álvarez, Departamento de Fisiología,
ENCB-IPN (dnicolas@ipn.mx). Licencia MIT. Por favor cita el artículo de
referencia indicado en la [portada](../README.es.md#cómo-citar).
