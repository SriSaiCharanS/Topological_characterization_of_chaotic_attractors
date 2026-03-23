from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import numpy as np

# FRM points with coordinates (yn, yn1) 
# Sort
sort_idx = np.argsort(yn)
yn_sorted = yn[sort_idx]
yn1_sorted = yn1[sort_idx]

# Smooth
yn1_smooth = gaussian_filter1d(yn1_sorted, sigma=5)

# Find extrema
peaks, _ = find_peaks(yn1_smooth, prominence=0.5)
troughs, _ = find_peaks(-yn1_smooth, prominence=0.5)
critical_indices = np.sort(np.concatenate([peaks, troughs]))
critical_values = yn_sorted[critical_indices]

def symbolic_dynamics_general_numeric(y, crit_vals):
    symbols = []
    crit_vals = np.sort(crit_vals)
    if len(crit_vals) == 0:
        symbols = ['0'] * len(y)
    else:
        for val in y:
            index = np.searchsorted(crit_vals, val, side='right')
            symbols.append(str(index))
    return symbols

symbols_sorted = symbolic_dynamics_general_numeric(yn_sorted, critical_values)

# ----------- Map Symbols to Time Order -----------

symbols_time_ordered = [''] * len(sort_idx)
for sorted_idx, original_idx in enumerate(sort_idx):
    symbols_time_ordered[original_idx] = symbols_sorted[sorted_idx]