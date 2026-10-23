#!/usr/bin/env python3
###############################################################################
# collect_free_energy.py — gather free_energy.sh results into ONE Excel file
#
# Author : Ricardo Aguilar-Garay
# Contact: ricardoaguilargaray@gmail.com
#
# Contributor: Dulce Estefanía Nicolás Álvarez
#              Departamento de Fisiología, Escuela Nacional de Ciencias
#              Biológicas, Instituto Politécnico Nacional (ENCB-IPN), Mexico
#              dnicolas@ipn.mx
#
# Usage:
#   python3 collect_free_energy.py <parent_folder> [<parent_folder> ...] [options]
#
#   Each parent folder is searched recursively for
#       .../Analyses/FreeEnergy_*/summary.csv     (written by free_energy.sh)
#   The NAME of each parent folder becomes the "Group" column (for example
#   "Microplasticos" or "Ligandos_Endogenos"), so systems from different
#   projects can be compared in the same table without moving any folder.
#
# Options:
#   -o FILE.xlsx        output file (default: ./FreeEnergy_summary_<date>.xlsx,
#                       i.e. in the folder where you RUN this command)
#   --labels A,B        custom group labels, one per parent folder, in order
#   --check-deps        only check that openpyxl is available (and how to install it)
#   --no-install        never offer to install openpyxl (write CSV files instead)
#
# Example:
#   python3 collect_free_energy.py ~/Documents/Microplasticos \
#           ~/Documents/Maydelis/Ligandos_Endogenos -o Energias_ronda2.xlsx
#
# Sheets: Summary (mean ± SEM, ready for a paper table), Data_kcal, Data_kJ,
#         Top_residues, Metadata, Legend.
# Requires openpyxl. If it is missing, the script: (1) looks for a conda env
# python that has it (gmxMMPBSA first) and re-runs itself with it, (2) prints an
# ALARM with the install commands, (3) offers to install it (interactive only),
# and (4) if not installed, writes the same tables as CSV files.
###############################################################################
import argparse, csv, datetime, glob, math, os, subprocess, sys

KCAL2KJ = 4.184
TEXT = {"system", "receptor", "ligand", "ligand_resname", "gromacs", "gmx_MMPBSA",
        "date", "flags", "notes", "path", "top5_residues"}

def num(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return x if x not in ("", None) else None

def find_runs(root):
    pat = os.path.join(root, "**", "Analyses", "FreeEnergy_*", "summary.csv")
    return sorted(glob.glob(pat, recursive=True))

def is_energy(col):
    """Columns in kcal/mol that must be converted for the kJ sheet."""
    if col.endswith(("_frames", "_nm")) or col in ("frac_unbound",):
        return False
    return (col.startswith(("dVDW", "dEEL", "dEGB", "dESURF", "dGgas", "dGsolv", "dH_",
                            "mTdS", "dG_C2", "LE_", "E_LJ", "E_Coul", "E_int", "sigma_")))

def pm(row, key, err="sem", nd=2):
    m, e = row.get(key + "_mean"), row.get("%s_%s" % (key, err))
    if not isinstance(m, float):
        return None
    return ("%.*f ± %.*f" % (nd, m, nd, e)) if isinstance(e, float) else "%.*f" % (nd, m)

def one(row, key, nd=2):
    v = row.get(key)
    return round(v, nd) if isinstance(v, float) else None

def has_openpyxl(py):
    try:
        return subprocess.run([py, "-c", "import openpyxl"], stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL, timeout=60).returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False

def conda_pythons():
    """Python executables of conda envs, gmxMMPBSA first."""
    bases = []
    for v in ("CONDA_PREFIX", "CONDA_EXE"):
        p = os.environ.get(v)
        if p:
            p = os.path.dirname(os.path.dirname(p)) if v == "CONDA_EXE" else p
            bases.append(p)
    home = os.path.expanduser("~")
    bases += [os.path.join(home, d) for d in ("miniconda3", "anaconda3", "miniforge3", "mambaforge")]
    bases += ["/opt/conda", "/opt/miniconda3", "/opt/anaconda3"]
    found = []
    for b in bases:
        b = os.path.realpath(b)
        root = b.split(os.sep + "envs" + os.sep)[0]
        cands = [os.path.join(b, "bin", "python")]
        cands += sorted(glob.glob(os.path.join(root, "envs", "*", "bin", "python")))
        cands.append(os.path.join(root, "bin", "python"))
        for c in cands:
            if os.path.isfile(c) and c not in found:
                found.append(c)
    found.sort(key=lambda c: 0 if "gmxMMPBSA" in c else 1)
    return [c for c in found if os.path.abspath(c) != os.path.abspath(sys.executable)]

def alarm_openpyxl():
    py = sys.executable
    bar = "!" * 75
    print(bar)
    print(" ⚠  ALARMA: no se encontró 'openpyxl' (necesario para escribir el Excel)")
    print("    Python usado: %s" % py)
    print(bar)
    print(" Cómo instalarlo (elige UNA opción):")
    print("   1) pip en tu usuario (Ubuntu 22.04 y anteriores):")
    print("        %s -m pip install --user openpyxl" % py)
    print("   2) Ubuntu 23.04+ / error 'externally-managed-environment':")
    print("        sudo apt install python3-openpyxl")
    print("      o bien")
    print("        %s -m pip install --user --break-system-packages openpyxl" % py)
    print("   3) En el ambiente de conda (no toca el python del sistema):")
    print("        conda activate gmxMMPBSA && pip install openpyxl")
    print("        python collect_free_energy.py <carpetas> ...")
    print(" Verificar después:   python3 collect_free_energy.py --check-deps")
    print(bar)

def ensure_openpyxl(allow_install=True, quiet_ok=False):
    """True if openpyxl can be used (possibly after re-running with another python)."""
    try:
        import openpyxl  # noqa: F401
        if not quiet_ok:
            print("✓ openpyxl %s  (%s)" % (openpyxl.__version__, sys.executable))
        return True
    except ImportError:
        pass
    if not os.environ.get("CFE_REEXEC"):
        for py in conda_pythons():
            if has_openpyxl(py):
                print("ℹ openpyxl no está en %s, pero sí en %s" % (sys.executable, py))
                print("  → se vuelve a ejecutar automáticamente con ese python.")
                os.environ["CFE_REEXEC"] = "1"
                sys.stdout.flush()
                os.execv(py, [py] + sys.argv)
    alarm_openpyxl()
    if allow_install and sys.stdin.isatty() and sys.stdout.isatty():
        try:
            ans = input(" ¿Instalarlo ahora con 'pip install --user openpyxl'? [s/N]: ").strip().lower()
        except EOFError:
            ans = ""
        if ans in ("s", "si", "sí", "y", "yes"):
            cmd = [sys.executable, "-m", "pip", "install", "--user", "openpyxl"]
            r = subprocess.run(cmd)
            if r.returncode != 0:
                print(" pip falló; reintentando con --break-system-packages ...")
                r = subprocess.run(cmd + ["--break-system-packages"])
            if r.returncode == 0:
                print(" ✓ openpyxl instalado. Reiniciando el script ...")
                sys.stdout.flush()
                os.execv(sys.executable, [sys.executable] + sys.argv)
            print(" ❌ No se pudo instalar automáticamente. Usa una de las opciones de arriba.")
    return False

def main():
    ap = argparse.ArgumentParser(description="Collect free_energy.sh results into one Excel file.")
    ap.add_argument("roots", nargs="*", default=["."], help="parent folder(s) to search")
    ap.add_argument("-o", "--output", default=None, help="output .xlsx")
    ap.add_argument("--labels", default=None, help="comma-separated group labels (one per root)")
    ap.add_argument("--check-deps", action="store_true", help="only check that openpyxl is available")
    ap.add_argument("--no-install", action="store_true", help="never offer to install openpyxl")
    a = ap.parse_args()

    if a.check_deps:
        sys.exit(0 if ensure_openpyxl(allow_install=not a.no_install) else 1)
    XLSX_OK = ensure_openpyxl(allow_install=not a.no_install)
    if not XLSX_OK:
        print("⚠ Se continúa SIN Excel: se escribirán archivos CSV con las mismas tablas.\n")

    labels = a.labels.split(",") if a.labels else None
    if labels and len(labels) != len(a.roots):
        sys.exit("❌ --labels must have one label per folder (%d folders, %d labels)" % (len(a.roots), len(labels)))

    rows, residues = [], []
    for i, root in enumerate(a.roots):
        root = os.path.abspath(os.path.expanduser(root))
        if not os.path.isdir(root):
            sys.exit("❌ Folder not found: " + root)
        group = labels[i] if labels else os.path.basename(root.rstrip("/"))
        runs = find_runs(root)
        print("📂 %s  →  group '%s': %d system(s)" % (root, group, len(runs)))
        for s in runs:
            with open(s, newline="") as fh:
                for r in csv.DictReader(fh):
                    r = {k: (v if k in TEXT else num(v)) for k, v in r.items()}
                    r["group"] = group
                    r["run_folder"] = os.path.dirname(s)
                    rows.append(r)
                    print("    ✓ %-45s %s" % (r.get("system"), r.get("flags") or ""))
            tr = os.path.join(os.path.dirname(s), "top_residues.csv")
            if os.path.isfile(tr):
                with open(tr, newline="") as fh:
                    for r in csv.DictReader(fh):
                        r = {k: num(v) for k, v in r.items()}
                        r["group"] = group
                        residues.append(r)
    if not rows:
        sys.exit("❌ No results found. Run free_energy.sh first (it writes Analyses/FreeEnergy_*/summary.csv).")
    rows.sort(key=lambda r: (str(r["group"]), str(r.get("system"))))

    # Several runs of the same system (e.g. different windows) are kept and marked
    seen = {}
    for r in rows:
        seen.setdefault((r["group"], r.get("system")), []).append(r)
    for k, lst in seen.items():
        if len(lst) > 1:
            for r in lst:
                r["notes"] = ((r.get("notes") or "") + " | duplicate system: %d runs found" % len(lst)).strip(" |")

    # ---------- sheet contents ----------
    summary_head = ["Group", "System", "Ligand", "Heavy atoms", "Charge",
                    "ΔVDW", "ΔEEL", "ΔEGB", "ΔESURF", "ΔG gas", "ΔG solv", "ΔH (MM-GBSA)",
                    "−TΔS (C2)", "ΔG (C2)", "LE = ΔH/HA",
                    "E_LJ", "E_Coul", "E_int", "E_LJ/HA",
                    "H-bonds", "Contacts <4Å", "Top residues", "Flags"]
    summary = []
    for r in rows:
        summary.append([r["group"], r.get("system"), r.get("ligand"), r.get("ligand_heavy_atoms"),
                        r.get("ligand_charge"),
                        pm(r, "dVDW"), pm(r, "dEEL"), pm(r, "dEGB"), pm(r, "dESURF"),
                        pm(r, "dGgas"), pm(r, "dGsolv"), pm(r, "dH"),
                        one(r, "mTdS_C2"), one(r, "dG_C2"), one(r, "LE_dH", 3),
                        pm(r, "E_LJ", "err"), pm(r, "E_Coul", "err"), pm(r, "E_int", "err"),
                        one(r, "E_LJ_per_heavy", 3),
                        pm({"x_mean": r.get("hbonds_mean"), "x_sd": r.get("hbonds_sd")}, "x", "sd"),
                        one(r, "contacts_4A_mean", 0), r.get("top5_residues"), r.get("flags")])

    meta_cols = ["group", "system", "receptor", "ligand", "ligand_resname", "sim_ns",
                 "window_start_ns", "window_end_ns", "frame_step_ps", "n_frames", "mmgbsa_frames",
                 "temperature_K", "gromacs", "gmx_MMPBSA", "date", "flags", "notes", "path"]
    skip = set(meta_cols) - {"group", "system", "flags"} | {"run_folder", "top5_residues"}
    first = ["group", "system", "ligand", "ligand_heavy_atoms", "ligand_charge", "receptor_charge"]
    data_cols = first + [c for c in rows[0].keys() if c not in first and c not in skip]
    for r in rows[1:]:
        data_cols += [c for c in r.keys() if c not in data_cols and c not in skip]

    def data_rows(conv):
        out = []
        for r in rows:
            line = []
            for c in data_cols:
                v = r.get(c)
                if conv and is_energy(c) and isinstance(v, float):
                    v = v * KCAL2KJ
                line.append(v)
            out.append(line)
        return out

    res_cols = ["group", "system", "rank", "residue", "total", "total_sd", "vdw", "elec",
                "polar_solv", "nonpolar_solv", "label"]
    res_rows = [[r.get(c) for c in res_cols] for r in residues]
    meta_rows = [[r.get(c) for c in meta_cols] for r in rows]

    legend = [
        ("Unidades", "Summary y Data_kcal en kcal/mol; Data_kJ en kJ/mol (×4.184)."),
        ("mean ± SEM", "Términos MM-GBSA: promedio ± error estándar de la media (SD/√N)."),
        ("E_LJ, E_Coul ± Err.Est.", "Energía de interacción directa receptor–ligando (mdrun -rerun, energygrps). "
                                    "Err.Est. = estimación de error por bloques de gmx energy. NO es energía libre."),
        ("ΔVDW, ΔEEL", "Van der Waals y electrostática en fase gas (MM)."),
        ("ΔEGB, ΔESURF", "Solvatación polar (Generalized Born, igb=5) y no polar (SASA)."),
        ("ΔG gas / ΔG solv", "ΔVDW+ΔEEL  /  ΔEGB+ΔESURF."),
        ("ΔH (MM-GBSA)", "ΔG gas + ΔG solv, sin entropía. Es la cantidad recomendada para comparar."),
        ("−TΔS (C2)", "Entropía por cumulante de 2º orden: σ²(ΔG gas)/(2kT). Poco confiable si σ > 3.6 kcal/mol (flag C2_UNRELIABLE)."),
        ("ΔG (C2)", "ΔH − TΔS(C2). Usar solo si no hay flag C2_UNRELIABLE."),
        ("LE = ΔH/HA", "Eficiencia de ligando: ΔH por átomo pesado (normaliza por tamaño)."),
        ("E_LJ/HA", "E_LJ por átomo pesado del ligando."),
        ("H-bonds", "Puentes de H receptor–ligando por frame (promedio ± SD)."),
        ("Contacts <4Å", "Pares de átomos receptor–ligando a menos de 0.4 nm (promedio)."),
        ("Flags", "OK, C2_UNRELIABLE, HIGH_SEM (SEM ΔH > 1), POSITIVE_dH, LIGAND_UNBOUND "
                  "(ligando a >4.5 Å en >10% de frames), LIGAND_DISPLACED (RMSD > 1 nm), FEW_FRAMES."),
        ("Group", "Nombre de la carpeta padre de donde se recolectó el sistema."),
    ]

    out = a.output or "FreeEnergy_summary_%s.xlsx" % datetime.date.today().isoformat()
    try:
        if not XLSX_OK:
            raise ImportError
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Font, PatternFill
        from openpyxl.utils import get_column_letter
    except ImportError:
        base = os.path.splitext(out)[0]
        for name, head, body in (("summary", summary_head, summary), ("data_kcal", data_cols, data_rows(False)),
                                 ("data_kJ", data_cols, data_rows(True)), ("top_residues", res_cols, res_rows),
                                 ("metadata", meta_cols, meta_rows)):
            with open("%s_%s.csv" % (base, name), "w", newline="") as fh:
                w = csv.writer(fh); w.writerow(head); w.writerows(body)
        print("═" * 75)
        print(" ⚠ SIN EXCEL (falta openpyxl). Tablas escritas como CSV:")
        for name in ("summary", "data_kcal", "data_kJ", "top_residues", "metadata"):
            print("     %s_%s.csv" % (os.path.abspath(base), name))
        print(" Instala openpyxl (ver ALARMA arriba) y vuelve a correr para obtener el .xlsx")
        print("═" * 75)
        return

    wb = Workbook()
    head_font, head_fill = Font(bold=True, color="FFFFFF"), PatternFill("solid", fgColor="1F4E78")
    flag_fill = PatternFill("solid", fgColor="FCE4D6")

    def sheet(ws, head, body, numfmt="0.00", flag_col=None):
        ws.append(head)
        for c in ws[1]:
            c.font, c.fill = head_font, head_fill
            c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        for line in body:
            ws.append(line)
        for row in ws.iter_rows(min_row=2):
            for c in row:
                if isinstance(c.value, float):
                    c.number_format = numfmt
            if flag_col is not None:
                v = row[flag_col].value
                if v and v != "OK":
                    row[flag_col].fill = flag_fill
        for j, h in enumerate(head, 1):
            width = max([len(str(h))] + [len(str(l[j-1])) for l in body if l[j-1] is not None])
            ws.column_dimensions[get_column_letter(j)].width = min(max(10, width + 2), 60)
        ws.freeze_panes = "C2"
        ws.row_dimensions[1].height = 32

    ws = wb.active; ws.title = "Summary"
    sheet(ws, summary_head, summary, flag_col=summary_head.index("Flags"))
    sheet(wb.create_sheet("Data_kcal"), data_cols, data_rows(False), "0.000",
          data_cols.index("flags") if "flags" in data_cols else None)
    sheet(wb.create_sheet("Data_kJ"), data_cols, data_rows(True), "0.000")
    sheet(wb.create_sheet("Top_residues"), res_cols, res_rows)
    sheet(wb.create_sheet("Metadata"), meta_cols, meta_rows, flag_col=meta_cols.index("flags"))
    wl = wb.create_sheet("Legend")
    wl.append(["Columna", "Significado"])
    for c in wl[1]:
        c.font, c.fill = head_font, head_fill
    for k, v in legend:
        wl.append([k, v])
    wl.column_dimensions["A"].width = 24; wl.column_dimensions["B"].width = 120
    wb.save(out)
    print("═" * 71)
    print(" Excel written: %s   (%d systems)" % (os.path.abspath(out), len(rows)))
    print("═" * 71)

if __name__ == "__main__":
    main()
