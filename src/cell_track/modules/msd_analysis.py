from scipy.stats import linregress
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Sequence
from tqdm import tqdm
import os
from .analysis_utils import get_timelag_indices

def msd_fit(x,alpha,D):
    
    return 4*D*x**alpha

def msd(track_dict):
    
    # get all possible tau values
    tau = set()
    for cell in track_dict:
        tau.update({t - cell['t'][0] for t in cell['t']})
    tau = sorted(tau)
    raw_data = {T: [] for T in tau}
    
    for cell in tqdm(track_dict, desc="Computing MSDs"):
        
        for tau_ in tau:
        
            t = np.array(cell['t'])
            indices_t0, indices_t1 = get_timelag_indices(t, tau_)
            
            x = np.array(cell['x'])
            y = np.array(cell['y'])
            dx = x[indices_t1] - x[indices_t0]
            dy = y[indices_t1] - y[indices_t0]
            squared_disp = dx**2 + dy**2
            
            raw_data[tau_].extend(squared_disp)
    
    MSDs = [np.mean(raw_data[tau_]) for tau_ in tau]
    dMSDs = [np.std(raw_data[tau_])/np.sqrt(len(raw_data[tau_])) for tau_ in tau]
    # hack to avoid infinite error at zeroth timepoint
    dMSDs[0] = dMSDs[1]
    
    return tau, MSDs, dMSDs

    
def fit_msd(t: List[int], y: List[float], fit_window: Sequence[int]):

    for i in range(1, len(t)):
        assert t[i] > t[i-1], "Time data must be sorted."
        
    [fit_from, fit_to] = fit_window
    
    xfit = np.log(t[fit_from:fit_to])
    yfit = np.log(y[fit_from:fit_to])
    
    LinReg = linregress(xfit, yfit)
    alpha = LinReg.slope
    diffusivity = np.exp(LinReg.intercept)/4
    
    alpha_err = LinReg.stderr
    diffusivity_err = diffusivity*LinReg.intercept_stderr
    
    return (alpha, diffusivity), (alpha_err, diffusivity_err), LinReg.rvalue**2


    