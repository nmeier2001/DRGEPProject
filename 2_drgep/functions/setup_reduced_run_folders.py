import re
import shutil
from pathlib import Path

def _extract_numeric_key(path: Path) -> tuple:
    """
    Extrahiert alle Zahlen (epsilon, delta, ...) aus dem Dateinamen für eine
    korrekte NUMERISCHE Sortierung.

    Wichtig: sorted(list(src_dir.glob("*.dat"))) sortiert Path-Objekte
    lexikographisch als Strings. Das führt bei Dateinamen wie
    "model_0.5_0.0.dat" vs. "model_0.52_0.0.dat" zu falscher
    Reihenfolge, weil '_' (ASCII 95) lexikographisch NACH Ziffern wie
    '2','4','6','8' (ASCII 50-56) einsortiert wird:
        "...0.5_..."  ->  Zeichen nach der '5' ist '_'  (95)
        "...0.52_..." ->  Zeichen nach der '5' ist '2'  (50)
    => "0.5_..." sortiert nach "0.52_...", "0.54_...", "0.56_...", "0.58_..."
    Numerisches Sortieren nach den extrahierten Float-Werten behebt das.
    """
    stem = path.stem.replace("model_", "")
    numbers = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", stem)
    return tuple(float(n) for n in numbers) if numbers else (float("inf"),)


def prepare_reduced_run_folders(source_dat_dir: str | Path, target_base_dir: str | Path):
    src_dir = Path(source_dat_dir)
    target_dir = Path(target_base_dir)

    if not src_dir.exists():
        raise FileNotFoundError(f"Quellordner nicht gefunden: {src_dir}")

    # Alle .dat-Dateien finden - NICHT mit sorted() (lexikographisch!), siehe
    # _extract_numeric_key() für die Begründung. Hier erstmal unsortiert einlesen.
    dat_files = list(src_dir.glob("*.dat"))
    if not dat_files:
        print(f"Keine .dat-Dateien in {src_dir} gefunden.")
        return

    print(f"Gefundene Dateien ({len(dat_files)}): {[f.name for f in dat_files]}")

    run1_file = None
    other_files = []

    # Exakte Prüfung auf die Epsilon=0.0 / Delta=0.0 Datei
    for f in dat_files:
        if "_0.0_0.0.dat" in f.name or f.name == "model_0.0.dat":
            run1_file = f
        else:
            other_files.append(f)

    # NUMERISCH sortieren (statt lexikographisch) anhand der im Dateinamen
    # enthaltenen Zahlen (z.B. epsilon, delta)
    other_files.sort(key=_extract_numeric_key)

    # Zusammenführen: Run 1 zuerst, danach der Rest numerisch aufsteigend sortiert
    ordered_files = []
    if run1_file:
        ordered_files.append(run1_file)
    else:
        print("WARNUNG: Keine explizite 0.0_0.0 Datei für Run 1 gefunden!")

    ordered_files.extend(other_files)

    print(f"Reihenfolge nach numerischer Sortierung: {[f.name for f in ordered_files]}")
    print(f"Starte Erstellung der Run-Ordner in: {target_dir.resolve()}\n")

    # Ordner anlegen und Dateien kopieren
    for idx, dat_path in enumerate(ordered_files, start=1):
        # Parameter-Suffix extrahieren (z.B. "0.52_0.0" aus "model_0.52_0.0.dat")
        param_suffix = dat_path.stem.replace("model_", "").replace("model", "").strip("_")
        if not param_suffix:
            param_suffix = "0.0_0.0"
        # HINWEIS (Fix): das ursprüngliche
        #   param_suffix, _ = param_suffix.split("_", 1)
        # hat bei Mehrfachparametern (z.B. "0.52_0.0") nur den ERSTEN Teil
        # ("0.52") behalten und delta stillschweigend verworfen. Das ist hier
        # nur deshalb nicht aufgefallen, weil delta bisher immer 0.0 war.
        # -> param_suffix bleibt jetzt vollständig erhalten (z.B. "0.52_0.0").

        param_suffix, _ = param_suffix.split("_", 1) if "_" in param_suffix else (param_suffix, "")
        folder_name = f"Run_{idx}_{param_suffix}"
        run_folder = target_dir / folder_name

        run_folder.mkdir(parents=True, exist_ok=True)
        target_file_path = run_folder / "model.dat"

        shutil.copy2(dat_path, target_file_path)
        print(f"  [+] Ordner '{folder_name}' erstellt -> 'model.dat' gespeichert.")

    print(f"\nFertig! {len(ordered_files)} Run-Ordner wurden erfolgreich vorbereitet.")