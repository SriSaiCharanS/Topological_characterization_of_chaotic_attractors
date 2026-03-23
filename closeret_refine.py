#System that is integrated should be deined and maxima 
#should be computed from the timeseries before running this code.

import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import argrelextrema
from scipy.optimize import leastsq


P = np.column_stack((a0_n, a0_n1)) #FRM points
#close return method of extracting UPOs from FRM
def close_return(P, max_period=6, threshold=1e-2):
    N = len(P)
    upos = []
    for tau in range(1, max_period + 1):
        for i in range(N - tau):
            dist = np.linalg.norm(P[i] - P[i + tau])
            if dist < threshold:
                upos.append((i, tau))
    return upos

upos = close_return(P)
print(f"Found {len(upos)} candidate UPOs")

# ----------- Shooting Residual Function -----------
#max_indices corresponds to index of maxima from the timeseies of variable a0 in the original trajectory
peak_states = np.column_stack((a0[max_indices], a1[max_indices], psi[max_indices]))

def shooting_residual(v, tau, gamma, delta):
    T = v[0]
    state = np.array(v[1:])
    init_state = state.copy()

    for _ in range(tau):
        sol = solve_ivp(plasma_system, [0, T], state, t_eval=np.linspace(0, T, 1000),
                        args=(gamma, delta), rtol=1e-9, atol=1e-9)
        a0_seg, a1_seg = sol.y[0], sol.y[1]
        peak_idx = argrelextrema(a0_seg, np.greater)[0]
        if len(peak_idx) == 0:
            return [1e6] * 4
        idx = peak_idx[0]
        state = [a0_seg[idx], a1_seg[idx], sol.y[2][idx]]

    return np.concatenate((np.array(state) - init_state, [0]))

def refine_shooting(init_state, tau, T_guess, gamma, delta):
    v0 = np.array([T_guess, *init_state])
    result, cov_x, infodict, mesg, ier = leastsq(
        shooting_residual, v0, args=(tau, gamma, delta),
        full_output=True, ftol=1e-6, maxfev=2000
    )
    residual = shooting_residual(result, tau, gamma, delta)
    error = np.linalg.norm(residual)
    print(f"   🔁 Final residual: {error:.3e}")
    if error > 1e-3:
        return None
    return result[1:], tau, result[0]

def integrate_upo(full_state, tau, T, gamma, delta):
    state = list(full_state)
    a0_list, a1_list = [], []
    for _ in range(tau):
        sol = solve_ivp(plasma_system, [0, T], state, t_eval=np.linspace(0, T, 1000),
                        args=(gamma, delta), rtol=1e-9, atol=1e-9)
        a0_seg, a1_seg, psi_seg = sol.y
        a0_list.append(a0_seg)
        a1_list.append(a1_seg)
        peak_idx = argrelextrema(a0_seg, np.greater)[0]
        if len(peak_idx) == 0:
            break
        idx = peak_idx[0]
        state = [a0_seg[idx], a1_seg[idx], psi_seg[idx]]
    return a0_list, a1_list
