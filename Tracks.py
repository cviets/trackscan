#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 12:58:21 2023

@author: chrisviets

TODO: 
    implement turning angle histogram maker
    
    write something that checks if something is a true splitting artifact by
     verifying that it's not just two tracks that started and ended nearby
     
    Do more safety testing on splitting and linking algorithms. Eg, for linking 
     algorithm, generate a hyperbolic-ish track with a split and see if linking
     algorithm can fix it. For splitting algorithm, see if hyperbolic-ish tracks
     are safe from breakpoints. 
"""

import curve_fitting as cf
import linking_functions as lf
import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress


class Tracks:
    
    def __init__(self, filename, x_column, y_column, t_column, cell_column):
        """
        Tracks object has filename and tracks attributes. 
        
        filename: file containing track data (must be .csv)
        tracks: list of dicts with keys 'x', 'y', 't' and values list of x, y, 
                and t data. Each dict represents one cell's track.
        x,y,t_column: column in filename containing x,y,t data. Must be input 
                      as an int (e.g., column A is 0)
        cell_column: column in filename identifying cell or track
        """
        
        def sort_track_by_time(track_dict):
            """
            Track data is often not ordered chronologically. This function

            """
            i_sort = np.argsort(track_dict['t'])
            track_dict['t'] = list(np.array(track_dict['t'])[i_sort])
            track_dict['x'] = list(np.array(track_dict['x'])[i_sort])
            track_dict['y'] = list(np.array(track_dict['y'])[i_sort])
        
        with open(filename, newline='') as csvfile:
            
            csv_read = csv.reader(csvfile,delimiter=',')
            
            prev_cell = None
            cur_dict = {'x': [], 'y': [], 't': []}
            tracks = []
                
            for row in csv_read:
                
                if row[x_column] == '' or not all([elt in '0123456789,. ' for elt in row[x_column]]):
                    continue
                
                cur_cell = row[cell_column]
                x = row[x_column]
                y = row[y_column]
                t = row[t_column]
                
                if cur_cell != prev_cell and prev_cell is not None:
                    sort_track_by_time(cur_dict)
                    tracks.append(cur_dict)
                    cur_dict = {'x': [], 'y': [], 't': []}
                
                cur_dict['x'].append(float(x))
                cur_dict['y'].append(float(y))
                cur_dict['t'].append(float(t))

                prev_cell = cur_cell
        
        sort_track_by_time(cur_dict)
        tracks.append(cur_dict)        
        self.tracks = tracks
        
        # save x,y,t,cell columns from original csvfile in case of new save
        self.filename = filename
        self.x_column = x_column
        self.y_column = y_column
        self.t_column = t_column
        self.cell_column = cell_column           
            
    def scale_space(self, factor):
        """
        Scale all x and y values by factor
        """
        for cell in self.tracks:
            cell['x'] = [x*factor for x in cell['x']]
            cell['y'] = [y*factor for y in cell['y']]
    
    def save_data(self, filename=None):
        
        if filename is None:
            filename = self.filename
        
        with open(filename, 'w', newline = '') as csvfile:
            csv_write = csv.writer(csvfile)
            header = [''] * (max(self.x_column, self.y_column, 
                                self.t_column, self.cell_column)+1)
            header[self.x_column] = 'POSITION_X'
            header[self.y_column] = 'POSITION_Y'
            header[self.t_column] = 'FRAME'
            header[self.cell_column] = 'TRACK_ID'
            
            csv_write.writerow(header)
            cur_cell = 0
            for cell in self.tracks:
                cur_cell += 1
                for i, x in enumerate(cell['x']):
                    newrow = header = [''] * (max(self.x_column, self.y_column, 
                                                 self.t_column, self.cell_column)+1)
                    newrow[self.x_column] = x
                    newrow[self.y_column] = cell['y'][i]
                    newrow[self.t_column] = cell['t'][i]
                    newrow[self.cell_column] = cur_cell
                    csv_write.writerow(newrow)
                    
            
    def dedrift(self):
        """
        Mutates Tracks object to dedrift all tracks using mean velocity at every 
        time point. Returns mean velocity at every time point for reference.
        """
        velocities = {t: [] for cell in self.tracks for t in cell['t']}
        del velocities[min(velocities.keys())]
        
        for cell in self.tracks:
            for i in range(1,len(cell['t'])):
                
                dx = cell['x'][i] - cell['x'][i-1]
                dy = cell['y'][i] - cell['y'][i-1]
                
                velocities[cell['t'][i]].append((dx, dy))
                
        vbar = {t: (sum(elt[0] for elt in velocities[t])/len(velocities[t]), 
                    sum(elt[1] for elt in velocities[t])/len(velocities[t])) 
                for t in velocities}
        
        for cell in self.tracks:
            for i,t in enumerate(cell['t']):
                if i == 0:
                    continue
                cell['x'][i] -= sum(vbar[T][0] for T in cell['t'][1:i+1])
                cell['y'][i] -= sum(vbar[T][1] for T in cell['t'][1:i+1])
        return vbar
            
    def link(self, max_time_gap, max_distance):
        """
        Links tracks together in three steps. 
        First, generates candidates of tracks to link to each track using criteria 
        specified by parameters. Then chooses best track to link using parabolic 
        least squares fits to both tracks. Last, interpolates the space between 
        the tracks, links them, and deletes the old one.
        
        Parameters
        ----------
        max_time_gap : int
            max time between end of one track and start of another to consider 
            linking the two.
        max_distance : float
            max distance between two tracks to consider linking them.

        Returns
        -------
        None.
        """
        
        linking_candidates = lf.get_linking_candidates(self.tracks, max_time_gap, max_distance)
        link_dict = lf.choose_linking_partners(linking_candidates, self.tracks)
        lf.link_partners(link_dict, self.tracks)
        
    def split_tracks(self):
        """
        Finds where a track abruptly changes direction and splits track into 
        two at that location. Abrupt changes in track direction are indicative 
        of erroneous track switching (ie, track corresponds to two different 
                                      cells that passed near each other).
        
        Acts in place.
        """
        
        new_tracks = []
        
        for i, cell in enumerate(self.tracks):
            
            X = cell['x']
            if len(X) <= 4:
                continue
            
            Y = cell['y']
            T = cell['t']
            
            x_brks = cf.find_breakpoints(X, T)
            y_brks = cf.find_breakpoints(Y, T)
            
            brks = x_brks.union(y_brks)
            
            removals = set()
            for idx in brks:
                if idx - 1 in brks or idx - 2 in brks:
                    removals.add(idx)
            brks -= removals
            
            brks = sorted(brks, reverse=True)
            
            for idx in brks:
                # idx = cell['t'].index(brk)
                assert len(cell['x'][idx + 2:]) > 1, f"{i=}"
                new_tracks.append({'x': cell['x'][idx + 2:], 
                                   'y': cell['y'][idx + 2:], 
                                   't': cell['t'][idx + 2:]})
                cell['t'] = cell['t'][:idx+1]
                cell['x'] = cell['x'][:idx+1]
                cell['y'] = cell['y'][:idx+1]
                
        self.tracks.extend(new_tracks)        
    
    def get_cell_msd(self, cell):
        """
        Returns MSD of particular cell in Tracks.
        
        Parameters
        ----------
        cell : dict
            Cell whose MSD is returned.

        Returns
        -------
        t : list of ints
            List of time points over which MSD is measured.
        MSD : list of floats
            List of MSD values corresponding to times in t.
        dMSD : list of floats
            List of errors in MSD values (1 sigma).

        """
        tau = {t - cell['t'][0] for t in cell['t']}
        tau = sorted(tau)
        raw_data = {T: [] for T in tau}
        
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
                
        MSD_data = {T: (np.mean(raw_data[T]), 
                    np.std(raw_data[T])/np.sqrt(len(raw_data[T]))) 
                for T in raw_data}
        
        t = []
        MSD = []
        dMSD = []
        
        for T in MSD_data:
            t.append(T)
            MSD.append(MSD_data[T][0])
            dMSD.append(MSD_data[T][1])
        dMSD[0] = dMSD[1]
        dMSD[-1] = dMSD[-2]
        
        return t, MSD, dMSD
    
    # def sort_cells_by_alpha(self, tf, linreg=False):
    #     alphas = []
    #     flag = False
    #     for i, cell in enumerate(self.tracks):
    #         t, MSD, dMSD = self.get_cell_msd(cell)
    #         if not linreg:
    #             try:
    #                 if 0 not in dMSD:
    #                     popt, pcov = curve_fit(MSD_fit, t[:tf], MSD[:tf], sigma=dMSD[:tf])
    #                 else:
    #                     popt, pcov = curve_fit(MSD_fit, t[:tf], MSD[:tf], sigma=None)
    #             except:
    #                 continue
    #             alphas.append(popt[1])
    #         else:
                
    #             x = np.log(t[1:tf])
    #             y = np.log(MSD[1:tf])
    #             alpha = linregress(x, y)[0]
    #             if alpha<0 and not flag:
    #                 plt.plot(t[:tf], MSD[:tf])
    #                 plt.show()
    #                 flag = True
    #             alphas.append(alpha)
            
    #     return alphas
        
    def get_cell_mean_stepsize(self, cell):
        dists = []
        for i in range(1, len(cell['t'])):
            cur_x = cell['x'][i]
            prev_x = cell['x'][i-1]
            
            cur_y = cell['y'][i]
            prev_y = cell['y'][i-1]
            dists.append(((cur_x - prev_x)**2 + (cur_y - prev_y)**2)**(1/2))
            
        return sum(dists)/len(dists)
    
    
    # def get_diffusing_cells(self, tumbling_freq, tolerance=0.25):
    #     """
    #     Returns indices of cells that do not exhibit ballistic motion at start 
    #     of their trajectories (MSD !~ t^2). The idea is to find the indices of
    #     cells that are "stuck", i.e., don't exhibit run and tumble behavior

    #     Parameters
    #     ----------
    #     tumbling_freq : int
    #         cell tumbling frequency (typically a few seconds at most)
    #     tolerance : float, optional
    #         Return cells that have alpha +/- tolerance away from 2. 
    #         The default is 0.25.

    #     Returns
    #     -------
    #     out : list
    #         indices of cells with no ballistic motion.

    #     """
    #     out = []
    #     tf = tumbling_freq
    #     for i, cell in enumerate(self.tracks):
    #         t, MSD, dMSD = self.get_cell_msd(cell)
    #         try:
    #             if 0 not in dMSD:
    #                 popt, pcov = curve_fit(MSD_fit, t[:tf], MSD[:tf], sigma=dMSD[:tf])
    #             else:
    #                 popt, pcov = curve_fit(MSD_fit, t[:tf], MSD[:tf], sigma=None)
    #         except:
    #             continue
    #         if abs(popt[1] - 2) > tolerance:
    #             out.append(i)
    #     return out
    
    # def delete_cells(self, indices):
    #     """
    #     Delete tracks from Tracks object with specified indices
    #     """
        
    #     indices = sorted(indices, reverse=True)
    #     for i in indices:
    #         del self.tracks[i]
            
if __name__ == '__main__':
    
    filename = '/Users/chrisviets/Documents/ribbeck/Neutrophils/tracks/MUC2/01_imdm/230627_imdm_01muc2_D1.csv'
    M = Tracks(filename, 4, 5, 8, 2)
    
    
    # # Divide tracks in two to solve track switching artifacts
    # M.split_tracks()
    # # Link tracks more appropriately to solve track splitting artifacts
    # M.link(3, 30)
    
    # v_drift = M.dedrift()
    # # ind = M.get_diffusing_cells(4, 0.3)
    # # M.delete_cells(ind)
    
    # # t, MSD, dMSD = M.get_msd()
    
    # # # plot MSD
    # # plt.figure(figsize=(10, 7))
    # # plt.errorbar(t, MSD, lw=2, yerr=dMSD)
    # # plt.xscale('log')
    # # plt.yscale('log')
    # # plt.xticks(fontsize=18)
    # # plt.yticks(fontsize=18)
    # # plt.xlabel('Time [frames]', fontsize=20)
    # # plt.ylabel('Mean squared displacement [$\mu$m$^2$]', fontsize=20)
    
    # # # get fitted MSD
    # # max_frame = 11
    # # popt, pcov = curve_fit(MSD_fit, t[:max_frame], MSD[:max_frame], sigma=dMSD[:max_frame])
    # # perr = np.sqrt(np.diag(pcov))
    
    # # # if got a bad fit the first time, try again
    # # if perr[0] > 50 or popt[1] < 0.05:
    # #     popt, pcov = curve_fit(MSD_fit, t[:max_frame], MSD[:max_frame], sigma=None)
    # #     perr = np.sqrt(np.diag(pcov))
    
    # # plt.plot(t, MSD_fit(t, *popt), lw=2, ls='--')
    # # plt.show()
    
    
    # # diffs1 = []
    # # kappa = []
    
    # for tf in [20]:
        
    #     alphas = []
    #     for i, cell in enumerate(M.tracks):
    #         if len(cell['t']) < 4:
    #             continue
    #         t, MSD, dMSD = M.get_cell_msd(cell)
    #         try:
    #             if 0 not in dMSD:
    #                 popt, pcov = curve_fit(MSD_fit, t[:tf], MSD[:tf], sigma=dMSD[:tf])
    #             else:
    #                 popt, pcov = curve_fit(MSD_fit, t[:tf], MSD[:tf], sigma=None)
    #             alphas.append(popt[1])
    #         except:
    #             continue
        
    #     # alphas_ = M.sort_cells_by_alpha(i, False)
        
    #     plt.figure(figsize=(10, 7))
    #     plt.title(r'Distribution of $\alpha$ values', fontsize=22)
    #     plt.xticks(fontsize=18)
    #     plt.yticks(fontsize=18)
    #     plt.xlabel(r'Cell MSD $\alpha$ values', fontsize=20)
    #     plt.ylabel('Frequency', fontsize=20)
    #     alphas = [elt for elt in alphas if 0 <= elt <= 2]
    #     cnts, bin_edges, _ = plt.hist(alphas, bins='auto', density=True)
    #     bins = (bin_edges[1:] + bin_edges[:-1])/2
    #     cnts_maxima = sorted(cnts, reverse=True)[:3]
    #     starred_ix = sorted([i for i, elt in enumerate(cnts) if elt in cnts_maxima])
    #     starred_bins = [bins[ix] for ix in starred_ix]
        
    #     midline = bins[starred_ix[2]] 
    #     plt.axvline(midline, color='red', ls='--')
        
    #     # plt.xlim(0, 2)
    #     plt.savefig('a_new.png', bbox_inches='tight', dpi=400)
    #     plt.show()
        
    #     plt.figure(figsize=(7, 7))
    #     plt.xticks(fontsize=18)
    #     plt.yticks(fontsize=18)
    #     for i, cell in enumerate(M.tracks):
    #         if len(cell['t']) < 4:
    #             continue
    #         t, MSD, dMSD = M.get_cell_msd(cell)
    #         try:
    #             if 0 not in dMSD:
    #                 popt, pcov = curve_fit(MSD_fit, t[:tf], MSD[:tf], sigma=dMSD[:tf])
    #             else:
    #                 popt, pcov = curve_fit(MSD_fit, t[:tf], MSD[:tf], sigma=None)
    #             if popt[1] >= midline:
    #                 plt.plot(cell['x'], [800 - elt for elt in cell['y']], color='blue')
    #             else:
    #                 plt.plot(cell['x'], [800 - elt for elt in cell['y']], color='red', alpha=0.5)
    #         except:
    #             continue
    #     plt.plot([1, 1], [2, 2], color='blue', label='Crawling')
    #     plt.plot([1, 1], [2, 2], color='red', alpha=0.5, label='Stationary')
    #     plt.xlabel('x [um]', fontsize=20)
    #     plt.ylabel('y [um]', fontsize=20)
    #     plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=18)
    #     plt.title('Separated trajectories', fontsize=22)
    #     plt.xlim([0, 800])
    #     plt.ylim([0, 800])
        
    #     plt.savefig('b_new.png', bbox_inches='tight', dpi=400)
    #     plt.show()
        
    # # out = []
    # # for cell in M.tracks:
    # #     if len(cell['t']) < 4:
    # #         continue
    # #     out.append(M.get_cell_max_displacement(cell))
        
    # # cnts, bin_edges, _ = plt.hist(out, bins='auto')
    # # bins = (bin_edges[1:] + bin_edges[:-1])/2
    # # popt, pcov = curve_fit(gauss_sum, bins, cnts)
    # # xdata = np.linspace(bins[0], bins[-1], 100)
    # # plt.plot(xdata, gauss_sum(xdata, *popt), color='red', ls='--')
    # # plt.show()
    
    # # cutoff = popt[4]
    # # for cell in M.tracks:
    # #     if len(cell['t']) < 4:
    # #         continue
    # #     stat = M.get_cell_max_displacement(cell)
    # #     if stat > cutoff:
    # #         plt.plot(cell['x'], cell['y'], color='blue')
    # #     else:
    # #         plt.plot(cell['x'], cell['y'], color='red', alpha=0.5)
            
    # # plt.show()
        
