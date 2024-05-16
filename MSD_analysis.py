#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 11:24:15 2024

@author: chrisviets
"""

from scipy.stats import linregress
import numpy as np
import matplotlib.pyplot as plt

def MSD_fit(x,M,alpha):
    
    return 4*M*x**alpha

def MSD(track_dict):
    """
    

    Parameters
    ----------
    track_dict : TYPE
        DESCRIPTION.

    Returns
    -------
    dict
        DESCRIPTION.

    """
    
    # get all possible tau values
    tau = set()
    for cell in track_dict:
        tau.update({t - cell['t'][0] for t in cell['t']})
    tau = sorted(tau)
    raw_data = {T: [] for T in tau}
    
    for cell in track_dict:
        
        for T in tau:
        
            t = np.array(cell['t'])
            t1_vals = t[t>=T]
            # gives indices of t for which t >= T:
            indices_t1 = np.searchsorted(t,t1_vals)
            # gives indices of t where t1_vals - T would be inserted:
            indices_t0 = np.searchsorted(t,t1_vals-T)
            
            valid_indices = t[indices_t0] == t1_vals - T
            
            indices_t0 = indices_t0[valid_indices]
            indices_t1 = indices_t1[valid_indices]
            
            x = np.array(cell['x'])
            y = np.array(cell['y'])
            dx = x[indices_t1] - x[indices_t0]
            dy = y[indices_t1] - y[indices_t0]
            squared_disp = dx**2 + dy**2
            
            raw_data[T].extend(squared_disp)
    
    MSDs = [np.mean(raw_data[T]) for T in tau]
    dMSDs = [np.std(raw_data[T])/np.sqrt(len(raw_data[T])) for T in tau]
    dMSDs[0] = dMSDs[1]
    
    return tau, MSDs, dMSDs

    
def fit_msd(t, y, fit_window):
    """
    

    """
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

def plot_msd(t, y, dy=None, fit_params=None):
    
    plt.errorbar(t[1:], y[1:], yerr=dy[1:])
    plt.xscale('log')
    plt.yscale('log')
    
    if fit_params is not None:
        
        t_all = np.logspace(t[1], t[-1], 100)
        y_fit = [MSD_fit(tau, fit_params[0], fit_params[1]) for tau in t_all]
        
        plt.plot(t_all, y_fit)
        
    return None
    