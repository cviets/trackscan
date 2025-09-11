from typing import List, Dict, Union
import numpy as np
from numpy.typing import NDArray
from tqdm import tqdm

def turning_angle_analysis(tracks: List[Dict[str, Union[int, float]]], dt: int) -> NDArray[np.float64]:

    out = []
    for cell in tqdm(tracks, desc="Computing turning angles"):

        t = np.array(cell['t'])
        idx_t0, idx_t1 = get_timelag_indices(t, dt)
        vx = np.gradient(cell['x'], t)
        vy = np.gradient(cell['y'], t)

        velocities_final = np.array([vx[idx_t1], vy[idx_t1]]).T
        velocities_initial = np.array([vx[idx_t0], vy[idx_t0]]).T

        for i in range(len(velocities_final)):
            u = velocities_initial[i]
            v = velocities_final[i]

            angle = np.arctan2(u[0]*v[1] - u[1]*v[0], u[0]*v[0] + u[1]*v[1])
            out.append(angle)
    return np.array(out)*180/np.pi
    
def get_timelag_indices(t: NDArray[np.int_], tau: int):
    
    t1_vals = t[t>=tau]
    indices_t1 = np.searchsorted(t,t1_vals)
    indices_t0 = np.searchsorted(t,t1_vals-tau)
    
    valid_indices = t[indices_t0] == t1_vals - tau
    indices_t0 = indices_t0[valid_indices]
    indices_t1 = indices_t1[valid_indices]

    return indices_t0, indices_t1

def get_ensemble_mean_speed(tracks: List[Dict[str, Union[int, float]]]) -> float:
    
    means = np.zeros(shape=len(tracks))
    for i, cell in tqdm(enumerate(tracks), desc="Computing mean speed"):
        t = np.array(cell['t'])
        vx = np.gradient(cell['x'], t)
        vy = np.gradient(cell['y'], t)
        speeds = np.sqrt(vx**2+ vy**2)
        means[i] = np.mean(speeds)

    return np.mean(means)

        