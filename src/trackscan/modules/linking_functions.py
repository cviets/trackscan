import numpy as np
from .curve_fitting import least_squares
from typing import Dict, List
"""
Start by finding candidates (newborn cells) to link to each cell when it dies. 
`linking_candidates` is a dict: {cellid: {set of cell ids}}. 
Keys indicate a cell that has died and value contains all cells that are candidates to link to key, 
depending on whether they meet the max_time_gap and max_distance requirements. 
No values in dict are empty sets.
"""
        
def get_linking_candidates(track_dict: List[Dict], max_time_gap, max_distance: float):
    """
    Generates candidates of tracks to link to each track using criteria specified by parameters
    """
    linking_candidates = {i: set() for i in range(len(track_dict))}
    
    t_all = sorted({T for cell in track_dict for T in cell['t']})
    queue = {time_gap: set() for time_gap in t_all 
             if time_gap <= max_time_gap}
    
    remaining_cells = {i for i in range(len(track_dict))}
    for t in t_all:
        
        # add newly dead cells to queue. 
        # Remove checked cells from remaining_cells
        queue[0] = {cell_id for cell_id in remaining_cells 
                    if max(track_dict[cell_id]['t']) == t}
        remaining_cells -= queue[0]
        
        # for cell in queue, check for newborns at distance <= max_distance
        newborns = {cell_id for cell_id in remaining_cells 
                    if min(track_dict[cell_id]['t']) == t}
        cells_currently_in_queue = {c for cell_set in queue.values() 
                                    for c in cell_set}
        for dead_cell_id in cells_currently_in_queue:
            dead_cell = track_dict[dead_cell_id]
            death_point = (dead_cell['x'][-1], dead_cell['y'][-1])
            for newborn_cell_id in newborns:
                newborn_cell = track_dict[newborn_cell_id]
                spawn_point = (newborn_cell['x'][0], newborn_cell['y'][0])
                dist = np.sqrt((spawn_point[0] - death_point[0])**2 + 
                               (spawn_point[1] - death_point[1])**2)
                if dist <= max_distance:
                    linking_candidates[dead_cell_id].add(newborn_cell_id)
        # update queue such that queue[n+1] = queue[n], queue[0] = set()
        time_gaps= sorted(queue.keys(), reverse=True)
        for idx, time_gap in enumerate(time_gaps[:-1]):
            queue[time_gap] = queue[time_gaps[idx+1]]
        queue[0] = set()         
    linking_candidates = {key: val 
                          for key, val in linking_candidates.items() 
                          if val != set()}
    
    return linking_candidates

def choose_linking_partners(linking_candidates, track_dict):
    """
    Chooses best track to link using parabolic least squares fits to both tracks.
    """
    link_dict = {}
    for dead_cell_id in linking_candidates:
        min_dist = -1
        best_candidate = -1
        cell_1 = track_dict[dead_cell_id]
        
        # arbitrarily fit last 3 points of track to a quadratic
        if len(cell_1['t']) >= 3:
            t1_fit = cell_1['t'][-4:]
            x1 = least_squares(2, cell_1['x'][-4:], t1_fit)
            y1 = least_squares(2, cell_1['y'][-4:], t1_fit)
        else:
            x1 = least_squares(1, cell_1['x'], cell_1['t'])
            y1 = least_squares(1, cell_1['y'], cell_1['t'])
        
        for spawn_cell_id in linking_candidates[dead_cell_id]: 
            
            cell_2 = track_dict[spawn_cell_id]
                
            if len(cell_2['t']) >= 3:
                t2_fit = cell_2['t'][:4]
                x2 = least_squares(2, cell_2['x'][:4], t2_fit)
                y2 = least_squares(2, cell_2['y'][:4], t2_fit)
            else:
                x2 = least_squares(1, cell_2['x'], cell_2['t'])
                y2 = least_squares(1, cell_2['y'], cell_2['t'])
            
            t_max = min(cell_2['t'])
            t_min = max(cell_1['t'])
            
            t_all =  {T for cell in track_dict for T in cell['t'] 
                      if T <= t_max and T >= t_min}
            
            
            dist = sum(np.sqrt((x1(t) - x2(t))**2 + (y1(t) - y2(t))**2) 
                       for t in t_all)/len(t_all)
        
            if dist < min_dist or min_dist == -1:
                best_candidate = spawn_cell_id
                min_dist = dist 
                interpolation_domain = sorted(t_all)
                x2_best = x2
                y2_best = y2
                
        link_dict[dead_cell_id] = (best_candidate, interpolation_domain, 
                                   (x1, y1), (x2_best, y2_best))
        
    return link_dict

def link_partners(link_dict, track_dict):
    """
    Interpolates the space between tracks, links them, and deletes the old one.
    """
    
    removals = set()
    
    while link_dict:
        key_removals = set()
        for key, (val, t_domain, (x1, y1), (x2, y2)) in link_dict.items():
            if val not in link_dict:
                
                dead_cell = track_dict[key]
                new_cell = track_dict[val]
                
                t_max = t_domain[-1]
                t_min = t_domain[0]
                
                if t_max != t_min:
                    weight = lambda t: (t - t_min)/(t_max - t_min)
                else:
                    weight = lambda t: 0.5
                X = lambda t: (1 - weight(t)) * x1(t) + weight(t) * x2(t)
                Y = lambda t: (1 - weight(t)) * y1(t) + weight(t) * y2(t)
                
                dead_cell['t'] = dead_cell['t'][:-1]
                dead_cell['x'] = dead_cell['x'][:-1]
                dead_cell['y'] = dead_cell['y'][:-1]
                dead_cell['t'].extend(t_domain)
                dead_cell['x'].extend([X(t) for t in t_domain])
                dead_cell['y'].extend([Y(t) for t in t_domain])
                dead_cell['t'].extend(new_cell['t'][1:])
                dead_cell['x'].extend(new_cell['x'][1:])
                dead_cell['y'].extend(new_cell['y'][1:])


                key_removals.add(key)
                removals.add(val)
        for key in key_removals:
            del link_dict[key]
                
    for idx in sorted(removals, reverse=True):
        del track_dict[idx]
    
