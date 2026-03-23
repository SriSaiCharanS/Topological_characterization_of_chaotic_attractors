#This code can plot both phase space and UPOs in I vs dI/dt phase space


import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt("upo_005_tau4.dat")

t = data[:, 0]
y = data[:, 1]

# Numerical derivative dy/dt
dy_dt = np.gradient(y, t)

# Plot y vs dy/dt
plt.figure(figsize=(6, 6))
plt.plot(y, dy_dt, lw=1)
plt.xlabel(r"$y$")
plt.ylabel(r"$\dot y$")
plt.title(r"Phase Plot: $y$ vs $\dot y$")
plt.grid(True)

plt.tight_layout()
plt.show()

out = np.column_stack((t, y, dy_dt))

# Save WITHOUT header
np.savetxt(
    "time_y_ydot_ee4.dat",
    out,
    fmt="%.8e"
)

