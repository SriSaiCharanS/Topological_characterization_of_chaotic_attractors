import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import os

# --- Parameters (bounded-chaos / EE-prone) ---
gamma = 1.978e5
tau   = 3.5e-9
N0    = 0.175
k0    = 0.17

# modulation: use Hz (paper lists kHz; convert)
a = 0.07                # modulation depth (example for bounded-chaos / EEs)
f = 185640.0            # Hz  (185.650 kHz in paper -> 185650 Hz)

# integration controls
T_total = 0.2          # total time (s). Increase if you want longer time-series.
dt_sample = 2e-8        # sampling interval for output (s)
t_eval = np.arange(0.0, T_total, dt_sample)

# initial condition near steady-state (I_S, N_S) with small perturbation
N_S = k0
I_S = gamma * (N0 / k0 - 1.0)   # analytic steady-state when no modulation
I0 = I_S * (1.0 + 1e-4)
N0_ic = N_S * (1.0 + 1e-4)
z0 = 0.0

OUTDIR = "results_laser_185640"
os.makedirs(OUTDIR, exist_ok=True)

# RHS where z acts as time-phase (dz/dt = 1)
def laser_system(t, y, gamma, tau, N0, k0, f, a):
    I, N, z = y
    dI_dt = ((N - k0 * (1 + a * np.cos(2.0 * np.pi * f * z))) * I) / tau
    dN_dt = (N0 - N) * gamma - I * N
    dz_dt = 1.0
    return [dI_dt, dN_dt, dz_dt]

print("Integrating laser system...")
sol = solve_ivp(
    fun=lambda tt, ss: laser_system(tt, ss, gamma, tau, N0, k0, f, a),
    t_span=(0.0, T_total),
    y0=[I0, N0_ic, z0],
    t_eval=t_eval,
    rtol=1e-8,
    atol=1e-9,
    method="RK45",      # adaptive; change to 'RK23' or custom RK4 if desired
    max_step=1e-4       # safety: limit internal step size
)

t_all = sol.t
I_all = sol.y[0]
N_all = sol.y[1]
z_all = sol.y[2]
print("Integration done. Samples:", len(t_all))

# discard a small transient (if present)
discard_until = 0.1
mask = t_all >= discard_until
t = t_all[mask]
I = I_all[mask]
N = N_all[mask]
z = z_all[mask]
print("Post-transient samples:", len(t))

# Phase-plane plot I vs N
plt.figure(figsize=(6,4))
plt.plot(I, N)
plt.xlabel("Intensity I")
plt.ylabel("Carrier N")
plt.title("Laser attractor (I vs N)")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "phase_I_N.png"))
plt.close()

# Time series and peaks on I(t)
plt.figure(figsize=(8,3))
plt.plot(t, I, linewidth=0.7, label="I(t)")
peaks, _ = find_peaks(I, height=None, distance=1)  # tweak distance/height if needed
if peaks.size > 0:
    plt.plot(t[peaks], I[peaks], "ro", markersize=3, label="peaks")
plt.xlabel("Time (s)")
plt.ylabel("Intensity I")
plt.title("I(t) with detected peaks")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "I_time_peaks.png"))
plt.close()

# First return map (I_n vs I_{n+1}) using detected peaks
if peaks.size >= 2:
    I_pk = I[peaks]
    plt.figure(figsize=(5,5))
    plt.plot(I_pk[:-1], I_pk[1:], "k.", markersize=2)
    plt.xlabel("$I_n$")
    plt.ylabel("$I_{n+1}$")
    plt.title("First return map (peak I)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "frm_Ipeaks.png"))
    plt.close()
else:
    print("Too few peaks for a FRM. Adjust T_total or peak-detection parameters.")

print("Saved plots to", OUTDIR)
