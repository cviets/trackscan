import matplotlib.pyplot as plt
import os
import numpy as np
from typing import List
from .msd_analysis import msd_fit
from numpy.typing import NDArray

def plot_msd(filename: str, t: List[int], y: List[float], dy: List[float], alpha: float, D: float, a_err: float, D_err: float, r2: float) -> None:

    fig, axs = plt.subplots(1, 1, figsize=(10,10))
    
    axs.errorbar(t[1:], y[1:], yerr=dy[1:])
    axs.set_xscale('log')
    axs.set_yscale('log')
        
    t_all = np.logspace(np.log10(t[1]), np.log10(t[-1]), 100)
    y_fit = [msd_fit(tau, alpha, D) for tau in t_all]
    
    axs.plot(t_all, y_fit, ls="--")
    axs.axvline(len(t)//2, ls="--", color='black')
    axs.text(0.05, 0.9, f"alpha = {alpha:.2f} +/- {a_err:.2f}\nD = {D:.2f} +/- {D_err:.2f}\nr^2 = {r2:.3f}", transform=axs.transAxes)

    fig.savefig(os.path.expanduser(filename), dpi=400, bbox_inches='tight')
        
    return None

def plot_turning_angles(filename: str, angles: NDArray[np.float64]):

    fig, axs = plt.subplots(1,1,figsize=(10,10))
    counts, bins = np.histogram(angles, bins='auto', density=True)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    axs.scatter(bin_centers, counts)
    axs.set_xlim(-185, 185)

    fig.savefig(os.path.expanduser(filename), dpi=400, bbox_inches='tight')

    return None