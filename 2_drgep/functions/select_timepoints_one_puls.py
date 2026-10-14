import math

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

def sliding_window(EN, window_size=5):
    """
    Berechnet den gleitenden Mittelwert von E/N über ein Fenster.
    """
    EN_log = []
    for x in EN:
        EN_log.append(math.log(x))
    score = np.zeros_like(EN)
    for idx, value in enumerate(EN):
        mean_before = np.mean(EN_log[max(0, idx - window_size):idx])
        mean_after = np.mean(EN_log[idx:min(len(EN_log), idx + window_size)])
        score[idx] = mean_before - mean_after

    # Find all peaks over the threshold of 0.1
    peaks = np.zeros_like(score)
    for idx, value in enumerate(score):
        if value > 0.01:
            peaks[idx] = value

    peaks_isolated, _ = find_peaks(peaks, height=0.03, distance=40)
    #print(f"Peaks found at indices: {peaks_isolated}, corresponding times: {[EN[i] for i in peaks_isolated]}")

    return score, peaks_isolated

def select_timepoints_plateau(time: np.array, EN: np.array, n_points=0):
    """
    Zeitpunktselektion für E/N-Verläufe mit Plateau-Struktur.
    Punkte:
    1. Absolutes Maximum (Spitze beim Anstieg)
    2. Plateau 1: stabiler Bereich nach Peak
    3. Übergang zwischen Plateaus: stärkster relativer Abfall
    4. Quasi-Gleichgewicht: letzter Zeitpunkt
    """

    _, peaks_idx = sliding_window(EN, window_size=5)
    #print(f"Peaks found at indices: {peaks_idx}")

    timepoints = {}
    # --- 1. Absolutes Maximum (Spitze beim Anstieg) ---
    idx_peak = np.argmax(EN)
    timepoints['EN_peak'] = {
        'idx':    idx_peak,
        'time':   time[idx_peak],
        'EN':     EN[idx_peak],
        'reason': 'Absolutes E/N Maximum'
    }

    # 2. Peak Plateau: schwacher Abstieg nach dem Peak
    idx_plateau = idx_peak + 16  # Beispiel: 73 Schritte nach dem Peak
    timepoints['plateau_1'] = {
        'idx':    idx_plateau,
        'time':   time[idx_plateau],
        'EN':     EN[idx_plateau],
        'reason': 'Plateau 1: stabiler Bereich nach Peak'
    }

    # 3. zweiter Peak: anstieg nahe null kurz vor abfall
    timepoints['2ed peak'] = {
        'idx':    peaks_idx[1],
        'time':   time[peaks_idx[1]],
        'EN':     EN[peaks_idx[1]],
        'reason': 'Übergang zwischen Plateaus'
    }

    # 4. Abfall vor letztem Plateau: stärkster Abfall
    timepoints['fall'] = {
        'idx':    peaks_idx[2],
        'time':   time[peaks_idx[2]],
        'EN':     EN[peaks_idx[2]],
        'reason': 'Stärkster Abfall vor letztem Plateau'
    }

    # 5. Afterglow: stabiler Bereich nach Abfall
    timepoints['afterglow'] = {
        'idx':    peaks_idx[2] + 3000,
        'time':   time[peaks_idx[2] + 3000],
        'EN':     EN[peaks_idx[2] + 3000],
        'reason': 'Afterglow: stabiler Bereich nach Abfall'
    }

    # 6. Quasi-Gleichgewicht: letzter Zeitpunkt
    timepoints['quasi_steady'] = {
        'idx':    len(time) - 1,
        'time':   time[len(time) - 1],
        'EN':     EN[len(EN) - 1],
        'reason': 'Quasi-Gleichgewicht, thermische Chemie dominiert'
    }

    # 7. Optional: weitere Zeitpunkte, falls n_points > 0
    if n_points > 0:
        existing_keys = list(timepoints.keys())
        
        for i in range(0, len(existing_keys) - 1):
            start_idx = timepoints[existing_keys[i]]['idx']
            end_idx   = timepoints[existing_keys[i + 1]]['idx']
            
            t_start = time[start_idx]
            t_end   = time[end_idx]
            
            # Gleichmäßig im Log-Raum → dann nächsten Index im Array finden
            t_between = np.logspace(
                np.log10(t_start), 
                np.log10(t_end), 
                n_points + 2
            )[1:-1]  # start und end abschneiden
            
            for j, t_target in enumerate(t_between):
                # Nächsten tatsächlichen Zeitpunkt im Array finden
                idx = np.argmin(np.abs(time - t_target))
                
                timepoints[f'between_{i}_{i+1}_pt{j}'] = {
                    'idx':    idx,
                    'time':   time[idx],
                    'EN':     EN[idx],
                    'reason': f'Zwischenzeitpunkt {j+1}/{n_points} zwischen '
                            f'{existing_keys[i]} und {existing_keys[i+1]}'
                }

    print(len(timepoints))
    return timepoints