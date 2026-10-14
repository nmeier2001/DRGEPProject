"""
create_run_report.py
=====================

Erzeugt einen übersichtlichen Markdown-Report für einen einzelnen
Model Parameterstudien-Run (inkl. optionaler DRGEP-Mechanismus-
reduktion), abgelegt im jeweiligen Run-Ordner.

Erwartete Ordnerstruktur (wie bei param_study_analysis.py):

    Model_Parameterstudie/
        create_run_report.py          <- dieses Script
        param_study_analysis.py
        functions/...
        Runs/
            <RUN_NAME>/
                <RUN_NAME>.pkl        <- Study-Pickle (von load_save_study erzeugt)
                dic_reduced_*.pkl     <- optional: DRGEP-Reduktionsmetadaten

Aufruf:
    python create_run_report.py <RUN_NAME>
    python create_run_report.py                     # fragt interaktiv nach dem Run-Namen
    python create_run_report.py <RUN_NAME> --reason "Testen von eps=2e6 Reduktion"

Das Script fragt (falls nicht per --reason übergeben) immer nach dem
Grund für den Run und legt anschließend im Run-Ordner folgende
Unterordner an:

    Report/               -> das eigentliche Markdown-Dokument + CSV-Exporte
    Error_Analysis/        -> Plots zu CO2/CO-Fehler & NRMSE der Elektronendichte
    Reduction_Overview/    -> Plots & Tabelle zur DRGEP-Reduktion (falls vorhanden)

WICHTIGE ANNAHMEN (bitte bei Bedarf anpassen):
  1. Die dic_reduced_*.pkl-Struktur wird als
        dic_reduced_data[powplug_key][epsilon][delta] = {...}
     angenommen (siehe compare_pkl.py). Zusätzliche Felder wie
     'targets', 'species_kept', 'species_eliminated' werden nur
     übernommen, wenn sie unter genau diesem Namen im dict stehen
     -> siehe REDUCTION_EXTRA_KEYS unten, dort ggf. Feldnamen ergänzen.
  2. Die Zuordnung dic_reduced-Eintrag -> Run-Index erfolgt (wie in
     compare_pkl.py) über die Reihenfolge der Einträge, da dic_reduced
     selbst keine Run-Indizes enthält. Zur Kontrolle werden die daraus
     abgeleiteten epsilon/delta-Werte in der Tabelle den direkt aus
     dem Study-Objekt gelesenen Parameterwerten gegenübergestellt.
  3. Als Referenz-/Baseline-Run gilt der Run, dessen param_value
     (bzw. alle Einträge, falls mehrere Parameter variiert wurden)
     numerisch ~0 ist (voller Mechanismus, eps=0, delta=0).
  4. NRMSE der Elektronendichte (via DeviationAnalyzer) wird jetzt nach
     (epsilon, sigma) gestaffelt, nicht mehr nur nach param_value. Die
     Zuordnung DeviationAnalyzer-Ergebnis -> Run-Index läuft über
     da.get_by_run_index(run_idx), nicht mehr über param_value/_hashable.
"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless: nur Dateien schreiben, kein plt.show()
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

from functions.functions_analysis import (
    load_study_cache,
    DeviationAnalyzer,
    add_EN_to_study,
)

BASE_DIR = Path.cwd()
RUNS_DIR = BASE_DIR / "1_create_data_and_preparation" / "Runs"

# Bekannte Felder pro (epsilon, delta)-Eintrag laut create_data_reduced_dic():
#   kept_species (Liste), eliminated_species (Liste), kept_reactions (int),
#   powplug, epsilon
KNOWN_REDUCTION_FIELDS = {"kept_species", "eliminated_species", "kept_reactions", "powplug", "epsilon"}

DEFAULT_CO2_ERROR_THRESHOLD = 2.0    # %
DEFAULT_NRMSE_THRESHOLD = 5.0        # %


# ---------------------------------------------------------------------------
# Kleine Hilfsfunktionen
# ---------------------------------------------------------------------------

def is_number(x) -> bool:
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


def df_to_markdown(df: pd.DataFrame, float_fmt: str = "{:.4g}") -> str:
    """Minimalistischer Markdown-Tabellen-Export ohne 'tabulate'-Abhängigkeit."""
    if df is None or df.empty:
        return "_Keine Daten verfügbar._"

    def fmt_cell(v):
        if isinstance(v, float) and not pd.isna(v):
            return float_fmt.format(v)
        if isinstance(v, (list, tuple, np.ndarray)):
            return ", ".join(str(x) for x in v)
        if pd.isna(v):
            return ""
        return str(v)

    headers = list(df.columns)
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(fmt_cell(row[c]) for c in headers) + " |")
    return "\n".join(lines)


def find_study_pickle(run_dir: Path) -> Path:
    """Findet die Haupt-Study-Pickle-Datei (schließt dic_reduced_*.pkl aus)."""
    candidates = [p for p in run_dir.glob("*.pkl") if not p.name.startswith("dic_reduced")]
    if not candidates:
        raise FileNotFoundError(
            f"Keine Study-Pickle-Datei in {run_dir} gefunden.\n"
            f"-> Erst param_study_analysis.py für diesen Run ausführen."
        )
    if len(candidates) == 1:
        return candidates[0]
    for c in candidates:
        if c.stem == run_dir.name:
            return c
    return max(candidates, key=lambda p: p.stat().st_mtime)


def find_dic_reduced_pickles(run_dir: Path) -> list[Path]:
    return sorted(run_dir.glob("dic_reduced*.pkl"))


def load_dic_reduced(path: Path) -> dict:
    with open(path, "rb") as f:
        return pickle.load(f)


# ---------------------------------------------------------------------------
# DRGEP-Reduktionsdaten (dic_reduced_*.pkl)
# ---------------------------------------------------------------------------

def build_reduction_table(dic_reduced_data: dict) -> pd.DataFrame:
    """
    Reale Struktur (siehe Reduktions-Script):
        dic_reduced[powplug] = {
            'targets': [...],                 <- Geschwister-Key, KEIN epsilon
            epsilon1: {delta1: {kept_species, eliminated_species, kept_reactions, powplug, epsilon}},
            epsilon2: {delta1: {...}},
            ...
        }
    -> flache Tabelle, eine Zeile pro (POWPLUG, epsilon, delta)-Kombination.
    """
    rows = []
    for powplug_key, powplug_data in dic_reduced_data.items():
        if not isinstance(powplug_data, dict):
            continue
        for epsilon, delta_dict in powplug_data.items():
            if not is_number(epsilon) or not isinstance(delta_dict, dict):
                continue  # z.B. der 'targets'-Key auf dieser Ebene - kein epsilon
            for delta, data in delta_dict.items():
                if not isinstance(data, dict):
                    continue
                kept_species = data.get("kept_species") or []
                eliminated_species = data.get("eliminated_species") or []
                row = {
                    "POWPLUG (Reduktion)": data.get("powplug", powplug_key),
                    "epsilon": float(data.get("epsilon", epsilon)),
                    "delta (sigma)": float(delta),
                    "num_species": len(kept_species) if kept_species else None,
                    "kept_reactions": data.get("kept_reactions"),
                    "num_eliminated_species": len(eliminated_species) if eliminated_species else None,
                }
                # Unbekannte/zusätzliche Felder generisch mit übernehmen
                for key, val in data.items():
                    if key not in KNOWN_REDUCTION_FIELDS and key not in row:
                        row[key] = val
                rows.append(row)
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["POWPLUG (Reduktion)", "epsilon", "delta (sigma)"]).reset_index(drop=True)
    return df


def map_dic_reduced_to_runs(dic_reduced_data: dict, study: dict) -> dict:
    """
    Ordnet dic_reduced-Einträge den Study-Run-Indizes zu, basierend auf der
    Reihenfolge (analog compare_pkl.py: Model nummeriert Runs ab 1).
    ACHTUNG: reine Annahme über die Reihenfolge - Kontrollspalten
    ('epsilon (dic_reduced)' etc.) im Report prüfen!
    """
    run_indices = sorted(study["runs"].keys())
    mapping = {}
    counter = 0
    for powplug_key, powplug_data in dic_reduced_data.items():
        if not isinstance(powplug_data, dict):
            continue
        for epsilon, delta_dict in powplug_data.items():
            if not is_number(epsilon) or not isinstance(delta_dict, dict):
                continue
            for delta, data in delta_dict.items():
                if not isinstance(data, dict):
                    continue
                if counter >= len(run_indices):
                    break
                run_idx = run_indices[counter]
                kept_species = data.get("kept_species") or []
                mapping[run_idx] = {
                    "epsilon (dic_reduced)": float(data.get("epsilon", epsilon)),
                    "delta (dic_reduced)": float(delta),
                    "num_species": len(kept_species) if kept_species else None,
                    "kept_reactions": data.get("kept_reactions"),
                    "POWPLUG (Reduktion)": data.get("powplug", powplug_key),
                }
                counter += 1
    return mapping


def extract_reduction_metadata(dic_reduced_data: dict) -> dict:
    """Extrahiert die Targets (Geschwister-Key von epsilon) je POWPLUG."""
    meta = {"POWPLUG-Werte (Reduktion)": list(dic_reduced_data.keys())}
    for powplug_data in dic_reduced_data.values():
        if isinstance(powplug_data, dict) and "targets" in powplug_data and "targets" not in meta:
            meta["targets"] = powplug_data["targets"]
    return meta


def build_eliminated_species_table(dic_reduced_data: dict) -> pd.DataFrame:
    """Zusatz-Übersicht: welche Spezies fallen bei welchem epsilon zusätzlich weg.
    Nützlich um zu sehen, ob bei einer bestimmten epsilon-Stufe eine wichtige
    Spezies (z.B. ein Zielprodukt) plötzlich eliminiert wird."""
    rows = []
    for powplug_key, powplug_data in dic_reduced_data.items():
        if not isinstance(powplug_data, dict):
            continue
        for epsilon, delta_dict in powplug_data.items():
            if not is_number(epsilon) or not isinstance(delta_dict, dict):
                continue
            for delta, data in delta_dict.items():
                if not isinstance(data, dict):
                    continue
                eliminated = data.get("eliminated_species") or []
                eliminated = [s for s in eliminated if s]  # leere Strings raus
                if not eliminated:
                    continue
                rows.append({
                    "POWPLUG": data.get("powplug", powplug_key),
                    "epsilon": float(data.get("epsilon", epsilon)),
                    "delta (sigma)": float(delta),
                    "Anzahl eliminiert": len(eliminated),
                    "Eliminierte Spezies": ", ".join(str(s) for s in eliminated),
                })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["POWPLUG", "epsilon", "delta (sigma)"]).reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# Fehleranalyse: CO2/CO ggü. Referenz, NRMSE der Elektronendichte
# ---------------------------------------------------------------------------

def find_baseline_run(study: dict, tol: float = 1e-9):
    """Sucht den Run mit param_value ~ 0 (voller Mechanismus / kein Offset)."""
    for idx, run in study["runs"].items():
        pv = run.get("param_value")
        if isinstance(pv, (list, tuple)):
            if pv and all(is_number(v) and abs(float(v)) < tol for v in pv):
                return idx
        elif is_number(pv) and abs(float(pv)) < tol:
            return idx
    return None


def compute_species_errors(study: dict, baseline_idx: int, species=("CO2", "CO")) -> dict:
    """Mittlerer relativer Fehler [%] je Run ggü. Baseline (analog compare_pkl.py)."""
    base_tec = study["runs"][baseline_idx]["tec"]
    errors = {sp: {} for sp in species}
    for idx, run in study["runs"].items():
        tec = run.get("tec")
        if tec is None:
            continue
        for sp in species:
            if sp not in tec.columns or sp not in base_tec.columns:
                continue
            base = base_tec[sp]
            diff = (base - tec[sp]) / base
            errors[sp][idx] = float(diff.mean() * 100)
    return errors


def compute_nrmse_E(study: dict, baseline_idx: int):
    """
    Nutzt DeviationAnalyzer. Die Baseline hat per Konvention epsilon=0.0 und
    sigma=0.0 (kein DRGEP-Reduktionsschritt) - das wird direkt so übergeben,
    unabhängig davon, in welchem Format study['runs'][...]['param_value']
    vorliegt (String/Liste/Float). DeviationAnalyzer identifiziert den
    Referenzrun jetzt über (epsilon, sigma) aus dic_reduced_*.pkl statt über
    den rohen param_value - baseline_idx wird hier nur noch für die
    CO2/CO-Fehleranalyse (compute_species_errors) gebraucht, nicht mehr für
    DeviationAnalyzer selbst.
    """
    try:
        da = DeviationAnalyzer(study, ref_param_value=0.0, ref_sigma=0.0, value_col="E")
        da.run()
        return da
    except Exception as e:
        print(f"  ⚠ NRMSE (Elektronendichte) konnte nicht berechnet werden: {e}")
        return None


def save_nrmse_plot(da: DeviationAnalyzer, out_path: Path) -> None:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)

    ax1.plot(da.t_ref, da.E_ref, "k-", lw=2,
             label=f"Referenz ({da.fmt_key(da.ref_param_value, da.ref_sigma)})")
    for (eps, sigma), r in sorted(da.results.items()):
        ax1.plot(r.t, r.values, label=da.fmt_key(eps, sigma))
    ax1.set_ylabel(da.value_col)
    ax1.legend(fontsize=7, ncol=2)
    ax1.set_title(f"{da.study.get('param_name', '')} – {da.value_col}(t) Abweichung")

    for (eps, sigma), r in sorted(da.results.items()):
        ax2.plot(r.t, r.rel_dev * 100, label=da.fmt_key(eps, sigma))
    ax2.axhline(0, color="k", lw=0.5)
    ax2.axvspan(*da.window, color="gray", alpha=0.1, label="Auswertefenster")
    ax2.set_xscale("log")
    ax2.set_ylabel("rel. Abweichung [%]")
    ax2.set_xlabel("Time (s)")

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def save_species_error_plot(errors: dict, species, out_path: Path,
                             epsilon_map: dict | None = None) -> None:
    """
    epsilon_map: {Run-Idx: epsilon}. Wenn übergeben, wird epsilon statt
    Run-Idx auf der x-Achse verwendet (Punkte werden nach epsilon sortiert
    geplottet); Runs ohne bekanntes epsilon werden ausgelassen.
    """
    fig, ax = plt.subplots(figsize=(9, 5))
    for sp in species:
        if sp not in errors or not errors[sp]:
            continue
        idxs = sorted(errors[sp].keys())
        if epsilon_map:
            pairs = [(epsilon_map[i], errors[sp][i]) for i in idxs if i in epsilon_map]
            pairs.sort(key=lambda p: p[0])
            xs = [p[0] for p in pairs]
            vals = [p[1] for p in pairs]
        else:
            xs = idxs
            vals = [errors[sp][i] for i in idxs]
        ax.plot(xs, vals, "o-", label=sp)

    ax.set_xlabel(f"$\\epsilon$" if epsilon_map else "Run-Index")
    ax.set_ylabel("Mittlerer Fehler [%]")
    # ax.set_title("Mittlerer Fehler ggü. Referenz (voller Mechanismus)")
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def save_nrmse_vs_epsilon_plot(run_table: pd.DataFrame, out_path: Path,
                                epsilon_map: dict) -> bool:
    """
    Plottet NRMSE der Elektronendichte (gesamt & ggf. Auswertefenster) gegen
    epsilon - analog zu save_species_error_plot, nur mit den bereits in
    run_table stehenden NRMSE-Spalten (NRMSE_E_gesamt_%, NRMSE_E_fenster_%)
    statt CO2/CO-Fehlern. Gibt True zurück, wenn tatsächlich etwas geplottet
    (und gespeichert) wurde.
    """
    if not epsilon_map:
        print("  ⚠ NRMSE-vs-epsilon-Plot übersprungen: keine epsilon-Zuordnung "
              "('epsilon' bzw. 'epsilon (dic_reduced)') in run_table vorhanden.")
        return False

    series_specs = [
        ("NRMSE_E_gesamt_%", "o-", "NRMSE (gesamt)"),
    ]
    available = [s for s in series_specs if s[0] in run_table.columns]
    if not available:
        print("  ⚠ NRMSE-vs-epsilon-Plot übersprungen: keine NRMSE-Spalten "
              "('NRMSE_E_gesamt_%' / 'NRMSE_E_fenster_%') in run_table vorhanden "
              "(vermutlich konnte die NRMSE-Berechnung nicht durchgeführt werden).")
        return False

    fig, ax = plt.subplots(figsize=(9, 5))
    plotted = False
    for col, marker, label in available:
        pairs = [
            (epsilon_map[run_idx], val)
            for run_idx, val in zip(run_table["Run"], run_table[col])
            # epsilon ≈ 0 ist der Referenzrun selbst - NRMSE ist dort
            # trivial 0 und verzerrt nur die Skala bei kleinen epsilon.
            if run_idx in epsilon_map and pd.notna(val) and abs(epsilon_map[run_idx]) > 1e-9
        ]
        if not pairs:
            print(f"  ⚠ NRMSE-vs-epsilon-Plot: Spalte '{col}' vorhanden, aber keine "
                  f"gültigen (epsilon, Wert)-Paare gefunden - übersprungen.")
            continue
        pairs.sort(key=lambda p: p[0])
        xs = [p[0] for p in pairs]
        vals = [p[1] for p in pairs]
        ax.plot(xs, vals, marker, label=label)
        plotted = True

    if not plotted:
        plt.close(fig)
        print("  ⚠ NRMSE-vs-epsilon-Plot übersprungen: keine plottbaren Daten.")
        return False

    print(f"  ✅ NRMSE-vs-epsilon-Plot gespeichert: {out_path}")

    ax.set_xlabel(r"$\epsilon$")
    ax.set_ylabel("NRMSE [%]")
    ax.set_title("NRMSE (Elektronendichte) vs. epsilon")
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return True


def save_species_error_vs_sigma_plot(errors: dict, species, out_path: Path,
                                      sigma_map: dict | None = None) -> bool:
    """
    Wie save_species_error_plot(), aber mit sigma (DRGEP-Reduktionsschwelle)
    statt epsilon auf der x-Achse. Der Referenzrun (epsilon=sigma=0) wird
    ausgeschlossen, da sein Fehler per Definition 0 ist und nur die Skala
    verzerrt (analog zur epsilon≈0-Filterung in save_nrmse_vs_epsilon_plot).
    """
    if not sigma_map:
        print("  ⚠ Fehler-vs-sigma-Plot übersprungen: keine sigma-Zuordnung vorhanden.")
        return False

    fig, ax = plt.subplots(figsize=(9, 5))
    plotted = False
    for sp in species:
        if sp not in errors or not errors[sp]:
            continue
        pairs = [
            (sigma_map[i], errors[sp][i])
            for i in errors[sp]
            if i in sigma_map and abs(sigma_map[i]) > 1e-9
        ]
        if not pairs:
            continue
        pairs.sort(key=lambda p: p[0])
        xs = [p[0] for p in pairs]
        vals = [p[1] for p in pairs]
        ax.plot(xs, vals, "o-", label=sp)
        plotted = True

    if not plotted:
        plt.close(fig)
        print("  ⚠ Fehler-vs-sigma-Plot übersprungen: keine plottbaren Daten.")
        return False

    print(f"  ✅ Fehler-vs-sigma-Plot gespeichert: {out_path}")

    ax.set_xlabel(r"$\sigma$")
    ax.set_ylabel("Mittlerer Fehler [%]")
    ax.set_title("Mittlerer Fehler ggü. Referenz (voller Mechanismus) vs. sigma")
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return True


def save_nrmse_vs_sigma_plot(run_table: pd.DataFrame, out_path: Path,
                              sigma_map: dict) -> bool:
    """
    Wie save_nrmse_vs_epsilon_plot(), aber mit sigma statt epsilon auf der
    x-Achse.
    """
    if not sigma_map:
        print("  ⚠ NRMSE-vs-sigma-Plot übersprungen: keine sigma-Zuordnung vorhanden.")
        return False

    series_specs = [
        ("NRMSE_E_gesamt_%", "o-", "NRMSE (gesamt)"),
    ]
    available = [s for s in series_specs if s[0] in run_table.columns]
    if not available:
        print("  ⚠ NRMSE-vs-sigma-Plot übersprungen: keine NRMSE-Spalten "
              "('NRMSE_E_gesamt_%' / 'NRMSE_E_fenster_%') in run_table vorhanden "
              "(vermutlich konnte die NRMSE-Berechnung nicht durchgeführt werden).")
        return False

    fig, ax = plt.subplots(figsize=(9, 5))
    plotted = False
    for col, marker, label in available:
        pairs = [
            (sigma_map[run_idx], val)
            for run_idx, val in zip(run_table["Run"], run_table[col])
            if run_idx in sigma_map and pd.notna(val) and abs(sigma_map[run_idx]) > 1e-9
        ]
        if not pairs:
            print(f"  ⚠ NRMSE-vs-sigma-Plot: Spalte '{col}' vorhanden, aber keine "
                  f"gültigen (sigma, Wert)-Paare gefunden - übersprungen.")
            continue
        pairs.sort(key=lambda p: p[0])
        xs = [p[0] for p in pairs]
        vals = [p[1] for p in pairs]
        ax.plot(xs, vals, marker, label=label)
        plotted = True

    if not plotted:
        plt.close(fig)
        print("  ⚠ NRMSE-vs-sigma-Plot übersprungen: keine plottbaren Daten.")
        return False

    print(f"  ✅ NRMSE-vs-sigma-Plot gespeichert: {out_path}")

    ax.set_xlabel(r"$\sigma$")
    ax.set_ylabel("NRMSE [%]")
    ax.set_title("NRMSE (Elektronendichte) vs. sigma")
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return True


def save_reduction_plots(run_table: pd.DataFrame, out_dir: Path) -> None:
    if "num_species" not in run_table.columns:
        return
    df = run_table.dropna(subset=["num_species"])
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    plotted = False
    if "CO2-Fehler [%]" in df.columns:
        ax.plot(df["num_species"], df["CO2-Fehler [%]"].abs(), "o-", label="CO2")
        plotted = True
    if "CO-Fehler [%]" in df.columns:
        ax.plot(df["num_species"], df["CO-Fehler [%]"].abs(), "s-", label="CO")
        plotted = True
    if plotted:
        ax.set_xlabel("Anzahl Spezies (kept)")
        ax.set_ylabel("|Mittlerer Fehler| [%]")
        ax.set_title("Fehler vs. Anzahl Spezies")
        ax.grid(True)
        ax.legend()
        fig.tight_layout()
        fig.savefig(out_dir / "error_vs_species.png", dpi=200)
    plt.close(fig)

    if "kept_reactions" in df.columns:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(df["num_species"], df["kept_reactions"], "o-")
        ax.set_xlabel("Anzahl Spezies")
        ax.set_ylabel("Anzahl Reaktionen")
        ax.set_title("Reaktionen vs. Spezies")
        ax.grid(True)
        fig.tight_layout()
        fig.savefig(out_dir / "species_vs_reactions.png", dpi=200)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 5))
        plotted = False
        if "CO2-Fehler [%]" in df.columns:
            ax.plot(df["kept_reactions"], df["CO2-Fehler [%]"].abs(), "o-", label="CO2")
            plotted = True
        if "CO-Fehler [%]" in df.columns:
            ax.plot(df["kept_reactions"], df["CO-Fehler [%]"].abs(), "s-", label="CO")
            plotted = True
        if plotted:
            ax.set_xlabel("Anzahl Reaktionen")
            ax.set_ylabel("|Mittlerer Fehler| [%]")
            ax.set_title("Fehler vs. Anzahl Reaktionen")
            ax.grid(True)
            ax.legend()
            fig.tight_layout()
            fig.savefig(out_dir / "error_vs_reactions.png", dpi=200)
        plt.close(fig)


# ---------------------------------------------------------------------------
# Run-Tabelle (ein Eintrag pro Run: Parameterwerte, POWPLUG, Konversion, CO-Produktion...)
# ---------------------------------------------------------------------------

def expand_param_values(study: dict, run_table: pd.DataFrame) -> pd.DataFrame:
    """Zieht bei mehreren variierten Parametern (param_name = Liste) jeden
    Parameterwert in eine eigene Spalte, direkt aus dem Study-Objekt."""
    param_name = study.get("param_name")
    if isinstance(param_name, list):
        for i, pname in enumerate(param_name):
            run_table[pname] = run_table["param_value"].apply(
                lambda v, i=i: v[i] if isinstance(v, (list, tuple)) and len(v) > i else None
            )
    return run_table


def build_run_table(study: dict, errors: dict | None = None,
                     reduction_lookup: dict | None = None) -> pd.DataFrame:
    rows = []
    for idx, run in study["runs"].items():
        nam = run.get("nam")
        tec = run.get("tec")
        row = {"Run": idx, "param_value": run.get("param_value")}

        if nam is not None:
            powplug = nam.get("POWPLUG")
            row["POWPLUG_peak [W/cm3]"] = max(powplug) if isinstance(powplug, list) else powplug
            row["XPOW"] = nam.get("XPOW")

        if tec is not None and "CO2" in tec.columns and len(tec["CO2"]) > 0:
            n_CO2_0 = tec["CO2"].iloc[0]
            n_CO2_final = tec["CO2"].min()
            row["CO2-Konversion [%]"] = (1 - n_CO2_final / n_CO2_0) * 100 if n_CO2_0 else None

            if "CO" in tec.columns and len(tec["CO"]) > 0:
                n_CO_0 = tec["CO"].iloc[0]
                n_CO_final = tec["CO"].max()
                # CO-Produktion als Anteil des anfänglichen CO2-Bestands (analog zur
                # CO2-Konversion oben) - zeigt, wie viel vom CO2-Feed als CO endet.
                row["CO-Produktion [%]"] = (n_CO_final - n_CO_0) / n_CO2_0 * 100 if n_CO2_0 else None

        if errors:
            row["CO2-Fehler [%]"] = errors.get("CO2", {}).get(idx)
            row["CO-Fehler [%]"] = errors.get("CO", {}).get(idx)

        if reduction_lookup and idx in reduction_lookup:
            row.update(reduction_lookup[idx])

        rows.append(row)

    df = pd.DataFrame(rows)
    df = expand_param_values(study, df)
    return df


def param_value_sort_key(v):
    """Sortierschlüssel für param_value (Skalar, String oder Liste/Tupel bei
    Mehrparameter-Runs), damit Anzeige-Tabellen numerisch aufsteigend sortiert
    werden können statt nach roher Run-Reihenfolge im Pickle."""
    if isinstance(v, (list, tuple)):
        return tuple(float(x) if is_number(x) else float("inf") for x in v)
    if is_number(v):
        return (float(v),)
    return (float("inf"), str(v))


def find_duplicate_error_groups(run_table: pd.DataFrame, cols=("CO2-Fehler [%]", "CO-Fehler [%]", "NRMSE_E_gesamt_%")) -> list[str]:
    """Findet Runs mit (nahezu) identischen Fehlerwerten - typischerweise ein
    Hinweis, dass DRGEP für benachbarte epsilon-Werte dieselbe Spezies-Menge
    behalten hat und der reduzierte Mechanismus dadurch identisch ist."""
    present = [c for c in cols if c in run_table.columns]
    if not present:
        return []
    df = run_table.dropna(subset=present).copy()
    if df.empty:
        return []
    df["_key"] = df[present].round(4).apply(lambda r: tuple(r), axis=1)
    notes = []
    for key, group in df.groupby("_key"):
        if len(group) > 1:
            runs = ", ".join(str(r) for r in sorted(group["Run"]))
            notes.append(f"Runs {runs}: identische Fehlerwerte ({', '.join(f'{c}={key[i]:.4g}' for i, c in enumerate(present))})")
    return notes


def recommend_best_reduction(run_table: pd.DataFrame, baseline_idx,
                              co2_threshold: float = DEFAULT_CO2_ERROR_THRESHOLD,
                              nrmse_threshold: float = DEFAULT_NRMSE_THRESHOLD) -> str:
    needed = {"CO2-Fehler [%]", "NRMSE_E_gesamt_%", "kept_reactions"}
    if not needed.issubset(run_table.columns):
        return ("_Für eine automatische Empfehlung fehlen Spalten "
                "(CO2-Fehler, NRMSE der Elektronendichte, kept_reactions)._")

    candidates = run_table.dropna(subset=list(needed))
    candidates = candidates[candidates["Run"] != baseline_idx]  # Baseline ist keine "Reduktion"
    if candidates.empty:
        return "_Keine reduzierten Runs (außer der Baseline) mit vollständigen Daten vorhanden._"
    ok = candidates[
        (candidates["CO2-Fehler [%]"].abs() <= co2_threshold)
        & (candidates["NRMSE_E_gesamt_%"].abs() <= nrmse_threshold)
    ]
    if ok.empty:
        return (f"_Kein Run erfüllt die Kriterien "
                f"(CO2-Fehler ≤ {co2_threshold:.1f}% und NRMSE(E) ≤ {nrmse_threshold:.1f}%). "
                f"Ggf. Schwellwerte über --co2-threshold / --nrmse-threshold anpassen._")

    best = ok.loc[ok["kept_reactions"].idxmin()]
    return (
        f"**Empfohlener reduzierter Mechanismus: Run {best['Run']}** "
        f"(epsilon={best.get('epsilon', best.get('epsilon (dic_reduced)', '?'))}, "
        f"delta={best.get('delta', best.get('delta (dic_reduced)', '?'))}) "
        f"mit {best.get('kept_reactions', '?')} Reaktionen / {best.get('num_species', '?')} Spezies — "
        f"CO2-Fehler={best['CO2-Fehler [%]']:.3f}%, NRMSE(E)={best['NRMSE_E_gesamt_%']:.3f}%.\n\n"
        f"Das ist unter den Runs, die beide Fehlerkriterien einhalten, derjenige mit "
        f"der stärksten Reduktion (wenigste Reaktionen)."
    )


# ---------------------------------------------------------------------------
# Markdown-Report
# ---------------------------------------------------------------------------

def generate_markdown(run_name: str, run_dir: Path, study: dict, reason: str,
                       run_table: pd.DataFrame, baseline_idx, reduction_df: pd.DataFrame,
                       reduction_meta: dict, has_reduction: bool,
                       co2_threshold: float, nrmse_threshold: float,
                       reduction_plot_files: set | None = None,
                       nrmse_plot_written: bool = False,
                       species_plot_written: bool = False,
                       nrmse_vs_epsilon_plot_written: bool = False,
                       species_vs_sigma_plot_written: bool = False,
                       nrmse_vs_sigma_plot_written: bool = False,
                       eliminated_df: pd.DataFrame | None = None) -> str:
    reduction_plot_files = reduction_plot_files or set()

    # Für die Anzeige nach param_value sortieren (nicht nach roher Run-Reihenfolge
    # im Pickle) - die 'Run'-Spalte bleibt als Rückverweis auf den Model-
    # Laufindex/Ordner erhalten.
    run_table_disp = run_table.copy()
    run_table_disp["_sort_key"] = run_table_disp["param_value"].apply(param_value_sort_key)
    run_table_disp = run_table_disp.sort_values("_sort_key").drop(columns="_sort_key")

    param_name = study.get("param_name")
    param_name_str = ", ".join(param_name) if isinstance(param_name, list) else str(param_name)
    n_runs = len(study.get("runs", {}))

    md = []
    md.append(f"# Run-Report: {run_name}")
    md.append("")
    md.append(f"*Erstellt am {datetime.now():%Y-%m-%d %H:%M}*")
    md.append("")
    md.append("## Inhalt")
    md.append("1. [Allgemeine Informationen](#1-allgemeine-informationen)")
    md.append("2. [Parameterwerte je Run](#2-parameterwerte-je-run)")
    if has_reduction:
        md.append("3. [Mechanismus-Reduktion (DRGEP)](#3-mechanismus-reduktion-drgep)")
    md.append("4. [Fehleranalyse ggü. Referenz](#4-fehleranalyse-ggü-referenz)")
    md.append("5. [Empfehlung](#5-empfehlung)")
    md.append("6. [Hinweise / Annahmen](#6-hinweise--annahmen)")
    md.append("")

    # 1. Allgemeine Informationen
    md.append("## 1. Allgemeine Informationen")
    md.append("")
    md.append(f"- **Run-Name:** {run_name}")
    md.append(f"- **Speicherpfad:** `{run_dir}`")
    md.append(f"- **Study-Datum / Uhrzeit:** {study.get('date')} / {study.get('time')}")
    md.append(f"- **Variierter Parameter:** {param_name_str}")
    md.append(f"- **Anzahl Runs:** {n_runs}")
    if baseline_idx is not None:
        ref_val = study["runs"][baseline_idx].get("param_value")
        md.append(f"- **Referenz-/Baseline-Run:** Run {baseline_idx} (param_value={ref_val})")
    else:
        md.append("- **Referenz-/Baseline-Run:** ⚠ nicht gefunden (kein Run mit param_value ≈ 0)")
    md.append(f"- **Grund für den Run:** {reason}")
    md.append("")

    # 2. Parameterwerte je Run
    md.append("## 2. Parameterwerte je Run")
    md.append("")
    display_cols = [c for c in [
        "Run", "param_value", "epsilon", "delta", "epsilon (dic_reduced)", "delta (dic_reduced)",
        "POWPLUG_peak [W/cm3]", "XPOW", "CO2-Konversion [%]", "CO-Produktion [%]",
    ] if c in run_table_disp.columns]
    md.append(df_to_markdown(run_table_disp[display_cols]))
    md.append("")
    md.append("_(sortiert nach param_value; `Run` = originaler Model-Laufindex. "
               "Vollständige Tabelle inkl. Fehlerwerten: siehe `Report/run_table.csv`)_")
    md.append("")

    # 3. Mechanismus-Reduktion
    if has_reduction:
        md.append("## 3. Mechanismus-Reduktion (DRGEP)")
        md.append("")
        powplug_vals = reduction_meta.get("POWPLUG-Werte (Reduktion)", [])
        md.append(f"- **POWPLUG der Daten, die zur Reduktion führten:** {', '.join(str(p) for p in powplug_vals)}")
        if "targets" in reduction_meta:
            md.append(f"- **Target-Spezies:** {', '.join(str(s) for s in reduction_meta['targets'])}")
        md.append("")
        md.append("### Übersicht epsilon / delta / Spezies- & Reaktionsanzahl")
        md.append("")
        core_cols = [c for c in [
            "POWPLUG (Reduktion)", "epsilon", "delta (sigma)",
            "num_species", "kept_reactions", "num_eliminated_species",
        ] if c in reduction_df.columns]
        md.append(df_to_markdown(reduction_df[core_cols]))
        md.append("")
        for fname, caption in (
            ("species_vs_reactions.png", "Reaktionen vs. Spezies"),
            ("error_vs_species.png", "Fehler vs. Anzahl Spezies"),
            ("error_vs_reactions.png", "Fehler vs. Anzahl Reaktionen"),
        ):
            if fname in reduction_plot_files:
                md.append(f"![{caption}](../Reduction_Overview/{fname})")
                md.append("")
        if not reduction_plot_files:
            md.append("_Keine Reduktions-Plots erzeugt (fehlende num_species/kept_reactions-Daten)._")
            md.append("")

        if eliminated_df is not None and not eliminated_df.empty:
            md.append("### Eliminierte Spezies je Epsilon-Stufe")
            md.append("")
            md.append(
                "_Zur Kontrolle, ob bei einer bestimmten epsilon-Stufe eine wichtige "
                "Spezies (z.B. ein Target) unerwartet wegfällt._"
            )
            md.append("")
            md.append(df_to_markdown(eliminated_df))
            md.append("")

    # 4. Fehleranalyse
    md.append("## 4. Fehleranalyse ggü. Referenz")
    md.append("")
    if baseline_idx is None:
        md.append("⚠ Keine Fehleranalyse möglich, da kein Referenz-Run (param_value ≈ 0) gefunden wurde.")
    else:
        error_cols = [c for c in [
            "Run", "param_value", "epsilon", "delta", "CO2-Fehler [%]", "CO-Fehler [%]",
            "NRMSE_E_gesamt_%", "NRMSE_E_fenster_%",
        ] if c in run_table_disp.columns]
        md.append("### CO2 / CO Fehler & NRMSE der Elektronendichte")
        md.append("")
        md.append(df_to_markdown(run_table_disp[error_cols]))
        md.append("")
        md.append("_(sortiert nach param_value; `Run` = originaler Model-Laufindex.)_")
        md.append("")

        dup_notes = find_duplicate_error_groups(run_table)
        if dup_notes:
            md.append(
                "**Hinweis:** Folgende Runs liefern (auf 4 Nachkommastellen) identische "
                "Fehlerwerte — meist ein Zeichen, dass DRGEP für die betroffenen "
                "epsilon/delta-Werte dieselbe Spezies-/Reaktionsmenge behalten hat, "
                "der reduzierte Mechanismus also identisch ist:"
            )
            for note in dup_notes:
                md.append(f"- {note}")
            md.append("")
        md.append("")
        if species_plot_written:
            md.append("![CO2/CO Fehler](../Error_Analysis/species_error.png)")
            md.append("")
        if species_vs_sigma_plot_written:
            md.append("![CO2/CO Fehler vs. sigma](../Error_Analysis/species_error_vs_sigma.png)")
            md.append("")
        if nrmse_plot_written:
            md.append("![NRMSE Elektronendichte](../Error_Analysis/nrmse_electron_density.png)")
            md.append("")
        else:
            md.append(
                "_NRMSE-Plot der Elektronendichte konnte nicht erzeugt werden "
                "(siehe Konsolenausgabe beim Ausführen des Scripts)._"
            )
            md.append("")
        if nrmse_vs_epsilon_plot_written:
            md.append("![NRMSE vs. epsilon](../Error_Analysis/nrmse_vs_epsilon.png)")
            md.append("")
        if nrmse_vs_sigma_plot_written:
            md.append("![NRMSE vs. sigma](../Error_Analysis/nrmse_vs_sigma.png)")
            md.append("")

    # 5. Empfehlung
    md.append("## 5. Empfehlung")
    md.append("")
    if has_reduction and baseline_idx is not None:
        md.append(recommend_best_reduction(run_table, baseline_idx, co2_threshold, nrmse_threshold))
    else:
        md.append("_Keine automatische Empfehlung möglich (keine Reduktionsdaten und/oder kein Referenz-Run)._")
    md.append("")

    # 6. Hinweise / Annahmen
    md.append("## 6. Hinweise / Annahmen")
    md.append("")
    md.append("- Baseline-Run wird automatisch als der Run mit param_value ≈ 0 identifiziert.")
    if has_reduction:
        md.append(
            "- Die Zuordnung von dic_reduced-Einträgen zu Run-Indizes erfolgt über die "
            "Reihenfolge im dic_reduced-dict (wie in `compare_pkl.py`). Zur Kontrolle stehen "
            "in der Tabelle sowohl die direkt aus dem Study-Objekt gelesenen Parameterwerte "
            "(`epsilon`, `delta`) als auch die aus dic_reduced abgeleiteten "
            "(`epsilon (dic_reduced)`, `delta (dic_reduced)`) — bei Abweichung bitte die "
            "Zuordnungslogik in `map_dic_reduced_to_runs()` prüfen."
        )
        md.append(
            "- Felder pro (epsilon, delta)-Eintrag werden aus `kept_species`, "
            "`eliminated_species`, `kept_reactions`, `powplug`, `epsilon` gelesen "
            "(wie in `create_data_reduced_dic()` erzeugt). Zusätzliche, unbekannte "
            "Felder werden automatisch mit in die Tabelle übernommen."
        )
    md.append(
        f"- Schwellwerte für die Empfehlung: CO2-Fehler ≤ {co2_threshold:.1f}%, "
        f"NRMSE(E) ≤ {nrmse_threshold:.1f}% (änderbar über --co2-threshold / --nrmse-threshold)."
    )
    md.append("")

    return "\n".join(md)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Erzeugt einen Markdown-Report für einen Model Parameterstudien-Run."
    )
    parser.add_argument("run", nargs="?", help="Name des Run-Ordners unter 'Runs/' (oder absoluter Pfad)")
    parser.add_argument("--reason", default=None, help="Grund für den Run (überspringt die interaktive Abfrage)")
    parser.add_argument("--co2-threshold", type=float, default=DEFAULT_CO2_ERROR_THRESHOLD)
    parser.add_argument("--nrmse-threshold", type=float, default=DEFAULT_NRMSE_THRESHOLD)
    args = parser.parse_args()

    run_name_arg = args.run or input("Name des Run-Ordners (unter 'Runs/'): ").strip()
    run_dir = Path(run_name_arg)
    if not run_dir.is_absolute():
        run_dir = RUNS_DIR / run_name_arg
    if not run_dir.exists():
        raise SystemExit(f"❌ Run-Ordner nicht gefunden: {run_dir}")

    reason = args.reason or input("Warum wurde dieser Run durchgeführt? ").strip()

    report_dir = run_dir / "Report"
    error_dir = run_dir / "Error_Analysis"
    reduction_dir = run_dir / "Reduction_Overview"
    for d in (report_dir, error_dir, reduction_dir):
        d.mkdir(exist_ok=True)

    print(f"\n📂 Lade Study aus: {run_dir}")
    study_path = find_study_pickle(run_dir)
    study = load_study_cache(study_path)

    first_tec = next((r["tec"] for r in study["runs"].values() if r.get("tec") is not None), None)
    if first_tec is not None and "EN_Td" not in first_tec.columns:
        try:
            add_EN_to_study(study)
        except Exception as e:
            print(f"  ⚠ E/N konnte nicht ergänzt werden: {e}")

    # --- DRGEP-Reduktionsdaten ---
    dic_paths = find_dic_reduced_pickles(run_dir)
    reduction_df = pd.DataFrame()
    reduction_lookup = {}
    reduction_meta = {}
    eliminated_df = pd.DataFrame()
    has_reduction = bool(dic_paths)
    if has_reduction:
        print(f"🔎 Gefundene dic_reduced-Dateien: {[p.name for p in dic_paths]}")
        combined_dic = {}
        for p in dic_paths:
            combined_dic.update(load_dic_reduced(p))
        reduction_df = build_reduction_table(combined_dic)
        reduction_lookup = map_dic_reduced_to_runs(combined_dic, study)
        reduction_meta = extract_reduction_metadata(combined_dic)
        eliminated_df = build_eliminated_species_table(combined_dic)

    # --- Baseline & Fehleranalyse ---
    baseline_idx = find_baseline_run(study)
    errors = {}
    da = None
    if baseline_idx is not None:
        print(f"📏 Baseline-Run: {baseline_idx}")
        errors = compute_species_errors(study, baseline_idx)
        da = compute_nrmse_E(study, baseline_idx)
    else:
        print("⚠ Kein Referenz-Run (param_value ≈ 0) gefunden — Fehleranalyse wird übersprungen.")

    # --- Run-Tabelle zusammenbauen ---
    run_table = build_run_table(study, errors=errors, reduction_lookup=reduction_lookup)

    nrmse_plot_written = False
    if da is not None:
        # da.results ist jetzt ueber (epsilon, sigma) geschluesselt, nicht
        # mehr ueber param_value - Zuordnung zur run_table laeuft daher ueber
        # den originalen Run-Index (get_by_run_index), nicht mehr ueber
        # param_key_str(param_value).
        nrmse_full_by_run = {}
        nrmse_window_by_run = {}
        for run_idx in study["runs"]:
            try:
                r = da.get_by_run_index(run_idx)
            except KeyError:
                continue
            nrmse_full_by_run[run_idx] = r.nrmse_full * 100
            nrmse_window_by_run[run_idx] = r.nrmse_window * 100
        run_table["NRMSE_E_gesamt_%"] = run_table["Run"].map(nrmse_full_by_run)
        run_table["NRMSE_E_fenster_%"] = run_table["Run"].map(nrmse_window_by_run)
        save_nrmse_plot(da, error_dir / "nrmse_electron_density.png")
        nrmse_plot_written = True

    # epsilon-Spalte/-Map einmal bestimmen (aus dem Study-Objekt oder, falls
    # nicht vorhanden, aus den über dic_reduced zugeordneten Werten) - wird für
    # den Fehler-vs-epsilon- und den NRMSE-vs-epsilon-Plot gemeinsam genutzt.
    eps_col = "epsilon" if "epsilon" in run_table.columns else (
        "epsilon (dic_reduced)" if "epsilon (dic_reduced)" in run_table.columns else None
    )
    if eps_col:
        epsilon_map = dict(zip(run_table["Run"], run_table[eps_col]))
    else:
        # Fallback: kein dic_reduced vorhanden (has_reduction=False) - aber
        # epsilon kann auch direkt der variierte Studienparameter sein, dann
        # steckt er als 'param_value' im Study-Objekt (bzw. als eigene Spalte,
        # falls param_name eine Liste mehrerer Parameter ist, siehe
        # expand_param_values()).
        param_name = study.get("param_name")
        if isinstance(param_name, str) and param_name.strip().lower() == "epsilon":
            eps_col = "param_value"
        elif isinstance(param_name, list):
            for pname in param_name:
                if isinstance(pname, str) and pname.strip().lower() == "epsilon" and pname in run_table.columns:
                    eps_col = pname
                    break
        epsilon_map = (
            {
                row["Run"]: float(row[eps_col])
                for _, row in run_table.iterrows()
                if eps_col is not None and is_number(row[eps_col])
            }
            if eps_col is not None else None
        ) or None
    print(f"  ℹ epsilon-Quelle für Plots: {eps_col!r} "
          f"({'gefunden' if epsilon_map else 'NICHT gefunden - epsilon-Plots werden übersprungen'})")

    # sigma-Zuordnung fuer die Sigma-Plots: da.run_meta wird bereits von
    # DeviationAnalyzer aus dic_reduced_*.pkl aufgebaut (siehe
    # compute_nrmse_E) und ist damit die verlaesslichste Quelle fuer
    # run_idx -> sigma, ohne die Zuordnungslogik hier noch einmal zu
    # duplizieren.
    sigma_map = (
        {idx: meta["sigma"] for idx, meta in da.run_meta.items()}
        if da is not None and da.run_meta else None
    ) or None
    print(f"  ℹ sigma-Quelle für Plots: "
          f"{'da.run_meta' if sigma_map else 'NICHT gefunden - sigma-Plots werden übersprungen'}")

    species_plot_written = False
    if errors and any(errors.get(sp) for sp in errors):
        save_species_error_plot(errors, ["CO2", "CO"], error_dir / "species_error.png",
                                 epsilon_map=epsilon_map)
        species_plot_written = True

    species_vs_sigma_plot_written = False
    if errors and any(errors.get(sp) for sp in errors):
        species_vs_sigma_plot_written = save_species_error_vs_sigma_plot(
            errors, ["CO2", "CO"], error_dir / "species_error_vs_sigma.png", sigma_map=sigma_map
        )

    nrmse_vs_epsilon_plot_written = save_nrmse_vs_epsilon_plot(
        run_table, error_dir / "nrmse_vs_epsilon.png", epsilon_map
    )

    nrmse_vs_sigma_plot_written = save_nrmse_vs_sigma_plot(
        run_table, error_dir / "nrmse_vs_sigma.png", sigma_map
    )

    if has_reduction:
        save_reduction_plots(run_table, reduction_dir)
    reduction_plot_files = {p.name for p in reduction_dir.glob("*.png")}

    # --- CSV-Exporte ---
    run_table.to_csv(report_dir / "run_table.csv", index=False)
    if not reduction_df.empty:
        reduction_df.to_csv(reduction_dir / "reduction_table.csv", index=False)
    if not eliminated_df.empty:
        eliminated_df.to_csv(reduction_dir / "eliminated_species_by_epsilon.csv", index=False)

    # --- Markdown schreiben ---
    md_text = generate_markdown(
        run_name=run_dir.name,
        run_dir=run_dir,
        study=study,
        reason=reason,
        run_table=run_table,
        baseline_idx=baseline_idx,
        reduction_df=reduction_df,
        reduction_meta=reduction_meta,
        has_reduction=has_reduction,
        co2_threshold=args.co2_threshold,
        nrmse_threshold=args.nrmse_threshold,
        reduction_plot_files=reduction_plot_files,
        nrmse_plot_written=nrmse_plot_written,
        species_plot_written=species_plot_written,
        nrmse_vs_epsilon_plot_written=nrmse_vs_epsilon_plot_written,
        species_vs_sigma_plot_written=species_vs_sigma_plot_written,
        nrmse_vs_sigma_plot_written=nrmse_vs_sigma_plot_written,
        eliminated_df=eliminated_df,
    )
    report_path = report_dir / f"Run_Report_{run_dir.name}.md"
    report_path.write_text(md_text, encoding="utf-8")

    # --- Verlauf protokollieren (Grund je Report-Lauf) ---
    history_path = report_dir / "run_history.md"
    with open(history_path, "a", encoding="utf-8") as f:
        f.write(f"- {datetime.now():%Y-%m-%d %H:%M} — {reason}\n")

    print(f"\n✅ Report gespeichert: {report_path}")
    print(f"   Plots:  {error_dir}" + (f"  &  {reduction_dir}" if has_reduction else ""))
    print(f"   Tabellen: {report_dir / 'run_table.csv'}"
          + (f", {reduction_dir / 'reduction_table.csv'}" if not reduction_df.empty else ""))


if __name__ == "__main__":
    main()