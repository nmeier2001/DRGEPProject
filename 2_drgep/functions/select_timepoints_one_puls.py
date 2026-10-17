import math

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

def sliding_window(EN, window_size=5):
    """
    Calculates the moving average of E/N over a window.
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
    Selects time points for E/N profiles with a plateau structure.
    Points:
    1. Absolute maximum (peak during the rise)
    2. Plateau 1: stable region after the peak
    3. Transition between plateaus: strongest relative decrease
    4. Quasi-equilibrium: final time point
    """

    _, peaks_idx = sliding_window(EN, window_size=5)
    #print(f"Peaks found at indices: {peaks_idx}")

    timepoints = {}
    # --- 1. Absolute maximum (peak during the rise) ---
    idx_peak = np.argmax(EN)
    timepoints['EN_peak'] = {
        'idx':    idx_peak,
        'time':   time[idx_peak],
        'EN':     EN[idx_peak],
        'reason': 'Absolutes E/N Maximum'
    }

    # 2. Peak plateau: slight decrease after the peak
    idx_plateau = idx_peak + 16  # Example: 73 steps after the peak
    timepoints['plateau_1'] = {
        'idx':    idx_plateau,
        'time':   time[idx_plateau],
        'EN':     EN[idx_plateau],
        'reason': 'Plateau 1: stabiler Bereich nach Peak'
    }

    # 3. Second peak: rise close to zero shortly before the decrease
    timepoints['2ed peak'] = {
        'idx':    peaks_idx[1],
        'time':   time[peaks_idx[1]],
        'EN':     EN[peaks_idx[1]],
        'reason': 'Übergang zwischen Plateaus'
    }

    # 4. Decrease before the final plateau: strongest decrease
    timepoints['fall'] = {
        'idx':    peaks_idx[2],
        'time':   time[peaks_idx[2]],
        'EN':     EN[peaks_idx[2]],
        'reason': 'Stärkster Abfall vor letztem Plateau'
    }

    # 5. Afterglow: stable region after the decrease
    timepoints['afterglow'] = {
        'idx':    peaks_idx[2] + 3000,
        'time':   time[peaks_idx[2] + 3000],
        'EN':     EN[peaks_idx[2] + 3000],
        'reason': 'Afterglow: stabiler Bereich nach Abfall'
    }

    # 6. Quasi-equilibrium: final time point
    timepoints['quasi_steady'] = {
        'idx':    len(time) - 1,
        'time':   time[len(time) - 1],
        'EN':     EN[len(EN) - 1],
        'reason': 'Quasi-Gleichgewicht, thermische Chemie dominiert'
    }

    # 7. Optional: additional time points if n_points > 0
    if n_points > 0:
        existing_keys = list(timepoints.keys())
        
        for i in range(0, len(existing_keys) - 1):
            start_idx = timepoints[existing_keys[i]]['idx']
            end_idx   = timepoints[existing_keys[i + 1]]['idx']
            
            t_start = time[start_idx]
            t_end   = time[end_idx]
            
            # Distribute points evenly in log space, then find the nearest array index
            t_between = np.logspace(
                np.log10(t_start), 
                np.log10(t_end), 
                n_points + 2
            )[1:-1]  # Exclude the start and end points
            
            for j, t_target in enumerate(t_between):
                # Find the nearest actual time point in the array
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