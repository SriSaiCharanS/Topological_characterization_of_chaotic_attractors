import numpy as np
import glob
from topoly import gln

def load_orbit(file, usecols=(-3, -2, -1)):
    """
    Load an orbit file and return its 3D trajectory (N x 3).
    By default, uses the last 3 columns: (a0, a1, psi).
    """
    data = np.loadtxt(file, comments=["#", "@"], skiprows = 2)
    return data[:, list(usecols)]

def compute_pairwise_linking(orbits):
    """
    Compute pairwise linking numbers between all orbits.
    """
    results = {}
    names = list(orbits.keys())
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            name1, name2 = names[i], names[j]
            xyz1 = orbits[name1].tolist()
            xyz2 = orbits[name2].tolist()
            lk = gln(xyz1,xyz2)
            results[(name1, name2)] = lk
    return results

if __name__ == "__main__":
    # Load all UPOs in current folder (change pattern if needed)
    files = sorted(glob.glob("upo_*.dat"))
    orbits = {f: load_orbit(f) for f in files}

        
    # --- Pairwise linking ---
    pairwise_linking = compute_pairwise_linking(orbits)
    print("\nPairwise linking numbers:")
    for (f1, f2), v in pairwise_linking.items():
        print(f"{f1} <-> {f2}: {v:.6f}")

