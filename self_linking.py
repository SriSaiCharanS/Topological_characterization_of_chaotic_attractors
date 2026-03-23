import numpy as np
from topoly import gln

def compute_self_linking(traj, epsilon=1e-2):
    """
    Compute the self-linking number of a 3D curve (traj) using a small push along
    the Frenet frame's normal vector.

    Parameters:
    traj : ndarray, shape (N,3)
        3D trajectory points of the curve
    epsilon : float
        Small displacement for the pushed curve

    Returns:
    float
        Self-linking number
    """
    traj = np.array(traj)
    N = len(traj)
    
    # Compute tangent vectors
    tangents = np.gradient(traj, axis=0)
    tangents /= np.linalg.norm(tangents, axis=1)[:, None]
    
    # Approximate normal using derivative of tangent
    dt = np.gradient(tangents, axis=0)
    normals = dt / np.linalg.norm(dt, axis=1)[:, None]
    
    # If dt is zero (straight segment), fallback to an arbitrary perpendicular
    for i in range(N):
        if not np.isfinite(normals[i]).all() or np.linalg.norm(normals[i]) < 1e-8:
            # Find a vector perpendicular to tangent
            t = tangents[i]
            if abs(t[0]) < 0.9:
                n = np.cross(t, [1,0,0])
            else:
                n = np.cross(t, [0,1,0])
            normals[i] = n / np.linalg.norm(n)
    
    # Create pushed curve
    pushed = traj + epsilon * normals
    
    # Compute linking number
    SL = gln(traj.tolist(), pushed.tolist())
    return SL

# Example usage
file = "upo_480_tau5_word11111.dat"
data = np.loadtxt(file, comments=["#", "@"], skiprows=2) 
xyz = data[:, -3:]  # shape (N, 3)

self_link = compute_self_linking(xyz)
print("Self-linking number:", self_link)


