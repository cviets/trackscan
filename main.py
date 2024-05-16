#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 12:38:47 2024

@author: chrisviets
"""

from Tracks import Tracks
import MSD_analysis as msd
import matplotlib.pyplot as plt
import numpy as np
import os

def get_csvfiles(folder):
    
    if not os.path.isdir(folder):
        raise Exception("Not a valid directory.")
    if folder[-1] == '/':
        folder = folder[:-1]
    
    ret = np.array([])
    
    subdirs = np.array([x for x in os.listdir(folder) if x != '.DS_Store'],dtype=object)
    subdirs = folder +'/'+subdirs
    
    for elt in subdirs:
        
        if os.path.isdir(elt):
            newcsv = get_csvfiles(elt)
            ret = np.concatenate([ret,newcsv])
            
        else:
            if os.path.isfile(elt):
                if elt[-4:] == '.csv':
                    ret = np.append(ret,elt)
    
    return ret

if __name__ == '__main__':
    
    filename = '/Users/chrisviets/Documents/ribbeck/E. coli/231026_LB/wt/z1/Run1.csv'
    cells = Tracks(filename, 4, 5, 8, 2)
    cells.dedrift()
    
    t, y, dy = msd.MSD(cells.tracks)
    msd.plot_msd(t, y, dy)
    
    # fit_params, fit_errors, r2 = msd.fit_msd(t, y, [1, 20])
    