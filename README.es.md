[![English](https://img.shields.io/badge/lang-English-blue.svg)](README.md)
[![Español](https://img.shields.io/badge/lang-Español-green.svg)](README.es.md)

# Gromacs_Scripts — Autogromacs, Basic Analyses y Free Energy

**Dinámica molecular automatizada con GROMACS para sistemas de CHARMM-GUI:
de la carpeta `gromacs/` recién descargada a los análisis y las energías de
unión, con un comando por paso.**

Desarrollado en el Laboratorio de Hormonas y Conducta, Escuela Nacional de
Ciencias Biológicas, Instituto Politécnico Nacional (ENCB-IPN), México.

```
gromacs/ de CHARMM-GUI
      │
      ▼
 Autogromacs.sh ──────────────► dynamic<N>ns.xtc / .tpr
 (minimización, equilibración,        │
  producción)                         ├──► basic_analyses_with_names.sh ──► Analyses/*.xvg, *.pdb
                                      │    (RMSD, puentes H/contactos, RMSF, Rg, snapshots)
                                      │
                                      └──► free_energy.sh ──► Analyses/FreeEnergy_*/summary.csv
                                           (MM-GBSA, LJ/Coulomb, validación)
                                                 │  todos los sistemas
                                                 ▼
                                           collect_free_energy.py ──► un solo Excel
```

---

## Herramientas

### 🧬 Autogromacs — correr la simulación

`Autogromacs.sh` corre todo el flujo de dinámica molecular (minimización →
equilibración → producción) en una GPU NVIDIA, por el tiempo que le pidas.

- **Un solo archivo, nada que editar:** cópialo a la carpeta `gromacs/` de
  CHARMM-GUI e indica los ns.
- **Cualquier duración:** la producción corre en tramos de 100 ns que se
  agregan a **una sola trayectoria continua** (`dynamic<N>ns.xtc`), sin
  necesidad de `trjcat`.
- **Reanudable:** tras un fallo, un corte de luz o un Ctrl-C, vuelve a
  correrlo y continúa desde el último checkpoint.
- **Recuperación de fallos de GPU:** si el modo todo-en-GPU falla (p. ej.
  CUDA error #700), reintenta solo en un modo estable.
- **Portátil:** el mismo archivo funciona en RTX 2080 / 3080 / 5090, con
  GROMACS 2024 y 2026.

```bash
bash Autogromacs.sh                                       # pregunta cuántos ns → teclea 200
nohup bash Autogromacs.sh 200 > autogromacs.log 2>&1 &   # o dalos directo, en segundo plano
```

📖 **[Documentación completa de Autogromacs →](docs/Autogromacs.es.md)**

### 📊 Basic Analyses — analizar la trayectoria

`basic_analyses_with_names.sh` hace los análisis estándar posteriores a la
producción, eligiendo los grupos de GROMACS **por nombre**, así que funciona
en cualquier sistema sin editar nada.

- **Trayectoria centrada**, **RMSD del ligando**, **puentes de hidrógeno
  proteína–ligando**, **RMSF por residuo**, **radio de giro** y **snapshots
  PDB** (inicio / mitad / final).
- **Cualquier ligando:** se detecta solo, como todo lo que hay en el grupo
  `SOLU` de CHARMM-GUI que no es proteína.
- **Microplásticos y otros ligandos de solo C/H:** los puentes de H se
  sustituyen automáticamente por contactos y distancia mínima
  proteína–ligando.
- **Salida ordenada** en `Analyses/` y, si existen, copias en las carpetas
  compartidas del proyecto con el nombre del sistema.

```bash
bash basic_analyses_with_names.sh
```

📖 **[Documentación completa de Basic Analyses →](docs/basic_analyses.es.md)**

### ⚡ Free Energy — energética de unión

`free_energy.sh` calcula qué tan fuerte se une el ligando a la proteína, y
`collect_free_energy.py` junta los resultados de todos tus sistemas en un
solo Excel.

- **Nada que elegir:** receptor y ligando se detectan desde `topol.top`
  (proteínas = receptor; agua e iones se ignoran; lo demás = ligando). Si el
  sistema es ambiguo, se detiene y te dice qué opción agregar.
- **MM-GBSA** con gmx_MMPBSA (GB-OBC II, mbondi2, sal 0.15 M): ΔH con
  promedio ± SEM, entropía C2, eficiencia de ligando y los **residuos
  clave** (descomposición por residuo).
- **Energía de interacción LJ / Coulomb directa** receptor–ligando
  (`mdrun -rerun` con grupos de energía).
- **Validación** sobre los mismos frames: puentes de H, contactos, distancia
  mínima y RMSD del ligando, con **flags de calidad** (ligando separado o
  desplazado, entropía poco confiable, error alto...).
- **Un Excel para todos los sistemas:** tabla resumen lista para artículo
  (promedio ± SEM), hojas en kcal/mol y kJ/mol, residuos clave, metadatos y
  leyenda.
- **Reanudable**, igual que las otras herramientas.

```bash
bash free_energy.sh --check      # ver qué detectó (no calcula nada)
bash free_energy.sh              # MM-GBSA + LJ/Coulomb + validación
python3 collect_free_energy.py /ruta/Proyecto_A /ruta/Proyecto_B -o Energias.xlsx
```

📖 **[Documentación completa de Free Energy →](docs/free_energy.es.md)** ·
guía rápida en texto plano: [`LEEME_free_energy.txt`](LEEME_free_energy.txt)

---

## Inicio rápido

```bash
cd /ruta/a/mi_sistema/gromacs            # salida de CHARMM-GUI Solution Builder
cp /ruta/a/Gromacs_Scripts/*.sh .         # las tres herramientas

# 1. Simular — elige UNA de las dos opciones:
bash Autogromacs.sh                                       # a) interactivo: pregunta
                                                          #    "¿Cuántos ns...?" → teclea 200 (o los que quieras)
nohup bash Autogromacs.sh 200 > autogromacs.log 2>&1 &   # b) en segundo plano: los ns como argumento

# 2. Analizar
bash basic_analyses_with_names.sh

# 3. Energética de unión
bash free_energy.sh

# 4. Cuando terminen todos los sistemas, desde cualquier carpeta: un Excel para todos
python3 collect_free_energy.py /ruta/a/carpeta_del_proyecto
```

La opción **a** es la más sencilla: corres el script y tecleas los ns cuando
te los pida. Deja la terminal abierta hasta que termine. La opción **b**
corre en segundo plano, así que puedes cerrar la terminal o la sesión; el
avance se sigue con `tail -f autogromacs.log`. Para reanudar después de una
interrupción, vuelve a correr `bash Autogromacs.sh` sin número.

**Requisitos:** Linux, GROMACS 2024/2026 con CUDA, una GPU NVIDIA (para la
simulación; los análisis corren en CPU) y una carpeta `gromacs/` de
CHARMM-GUI Solution Builder. Free Energy necesita además un ambiente de conda
con [gmx_MMPBSA](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/) y, para el
Excel, `openpyxl`. Los detalles están en la documentación de cada
herramienta.

## Próximamente

- Opción **MM-PBSA (Poisson–Boltzmann)** en Free Energy, además del MM-GBSA
  actual.
- Soporte para sistemas de CHARMM-GUI Membrane Builder.

## Autores y contribuidores

- **Ricardo Aguilar-Garay**: autor y responsable
  (ricardoaguilargaray@gmail.com)
- **Dulce Estefanía Nicolás Álvarez**: contribuidora. Departamento de
  Fisiología, Escuela Nacional de Ciencias Biológicas, Instituto Politécnico
  Nacional (ENCB-IPN), México (dnicolas@ipn.mx)

Preguntas, reportes de errores y sugerencias son bienvenidos.

## Cómo citar

Si estas herramientas te son útiles en tu investigación, **por favor cita el
artículo en el que se desarrolló y validó la metodología de DM + MM-PBSA**:

> Hernandez-Munguia, A.M.; Reyes-Chaparro, A.; Fregoso-Aguilar, T.A.;
> Tenorio-Barajas, A.Y.; Mendoza-Pérez, J.A.; **Aguilar-Garay, R.**;
> Nicolás-Álvarez, D.E. Binding Affinity and Interaction Profiles of
> Erinacines and Erinacerins with iNOS and NF-κB Revealed by Molecular
> Dynamics Simulations. *Int. J. Mol. Sci.* **2026**, *27*(7), 3145.
> https://doi.org/10.3390/ijms27073145

BibTeX:

```bibtex
@Article{ijms27073145,
  AUTHOR         = {Hernandez-Munguia, April Michelle and Reyes-Chaparro, Andrés and Fregoso-Aguilar, Tomas Alejandro and Tenorio-Barajas, Aldo Yair and Mendoza-Pérez, Jorge Alberto and Aguilar-Garay, Ricardo and Nicolás-Álvarez, Dulce E.},
  TITLE          = {Binding Affinity and Interaction Profiles of Erinacines and Erinacerins with iNOS and NF-κB Revealed by Molecular Dynamics Simulations},
  JOURNAL        = {International Journal of Molecular Sciences},
  VOLUME         = {27},
  YEAR           = {2026},
  NUMBER         = {7},
  ARTICLE-NUMBER = {3145},
  URL            = {https://www.mdpi.com/1422-0067/27/7/3145},
  ISSN           = {1422-0067},
  DOI            = {10.3390/ijms27073145}
}
```

En GitHub también puedes usar el botón **"Cite this repository"** (generado
a partir del archivo `CITATION.cff` de este repositorio).

## Licencia

Publicado bajo la **Licencia MIT** (ver `LICENSE`). Puedes usar, modificar y
redistribuir este software libremente; en trabajos académicos que se
beneficien de él, se agradece citar el artículo de arriba.

## Agradecimientos

Desarrollado con apoyo del programa de estancias posdoctorales de SECIHTI.
El flujo en GPU se validó en sistemas de gelsolina y de complejos
microplástico–proteína en la ENCB-IPN.
