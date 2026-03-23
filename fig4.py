import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.optimize import leastsq

from scipy.optimize import minimize
import csv
import os
import matplotlib as mpl

mpl.rcParams['agg.path.chunksize'] = 10000

# -------------------------
# Parameters (modify if needed)
# -------------------------

r= 1/400
alpha = 1/180
mu=0.0026
Tr= 29.5
Tr0= 16
H = 100
L= 15e6
zeta = 1.3
z0= 75
h = 62
eps = 0.0984
blmu_beta = 22

T_total = 20e5      # total integration time
dt = 0.01
t_eval = np.arange(0, T_total, dt)

discard_until = 10e5  # transient to discard before analysis
cutoff_frac = 0.10   # discard first 10% of FRM peaks as FRM transient

# UPO detection settings
max_period = 16
close_tol = 1e-2

OUTDIR = "results_upo_enso"
os.makedirs(OUTDIR, exist_ok=True)

# -------------------------
# RHS
# -------------------------
def enso_system(t, y, r, alpha, mu, Tr, Tr0, H, L, zeta, z0, h, eps, blmu_beta):
    h1, T1, T2 = y
    dh1_dt = r*(-h1-(blmu_beta*(T2-T1))/2)
    dT1_dt = -alpha*(T1-Tr)-eps*mu*(T2-T1)**2
    dT2_dt = -alpha*(T2-Tr)+ zeta*mu*(T2-T1)*(T2-Tr+0.5*(Tr-Tr0)*(1-np.tanh((H+h1+blmu_beta*(T2-T1)-z0)/h)))
    return [dh1_dt, dT1_dt, dT2_dt]

print("Integrating forced Liénard system (this can be slow)...")
sol = solve_ivp(lambda tt, ss: enso_system(tt, ss,r, alpha, mu, Tr, Tr0, H, L, zeta, z0, h, eps, blmu_beta ),
                (0.0, T_total), [76.406200,   27.275902,   20.329884], t_eval=t_eval, rtol=1e-8, atol=1e-9)
t_all = sol.t
h1_all = sol.y[0]
T1_all = sol.y[1]
T2_all = sol.y[2]
print("Integration done. Total samples:", len(t_all))

# -------------------------
# Discard transient
# -------------------------
mask = t_all >= discard_until
t = t_all[mask].copy()
h1 = h1_all[mask].copy()
T1 = T1_all[mask].copy()
T2 = T2_all[mask].copy()
print("Using post-transient samples:", len(t))

plt.plot(T1, T2, color='blue', linewidth='0.8')
plt.xlabel(r'$T_1$', fontsize = "15")
plt.ylabel(r'$T_2$', fontsize = "15")
plt.xticks(fontsize = "15")
plt.yticks(fontsize = "15")
plt.title("ENSO attractor- Bounded Chaos")
plt.savefig("ENSO_bc.png")

# Detect local maxima of a0
max_indices = argrelextrema(T2, np.greater)[0]
T2_max = T2[max_indices]

# Optional: plot FRM
T2_n = T2_max[:-1]
T2_n1 = T2_max[1:]

plt.figure(figsize=(8,6))
plt.plot(T2_n, T2_n1, 'k.', markersize=1)
#plt.axvline(a0_critical, color='r', linestyle='--', label='Critical point')
plt.xlabel(r'$T_2^{(n)}$', fontsize = "15")
plt.ylabel(r'$T_2^{(n+1)}$', fontsize = "15")
plt.title("First Return Map (FRM)")
plt.xticks(fontsize = "15")
plt.yticks(fontsize = "15")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("frm_enso.png")
