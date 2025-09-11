from typing import Dict, List, Union, Tuple
import numpy as np
from .curve_fitting import find_breakpoints

def scale_space(track_data: List[Dict[str, Union[int, float]]], factor: float) -> None:
    """
    Scale all x and y values by factor.
    """
    for cell in track_data:
        cell['x'] = [x*factor for x in cell['x']]
        cell['y'] = [y*factor for y in cell['y']]

def print_first_n(track_data: List[Dict[str, Union[int, float]]], n: int=5) -> None:
    count = 0
    cell_num = 0
    idx = 0

    header = ["Track ID", "Position_X", "Position_Y", "Frame"]

    print(f"{header[0]:^12}", "|", f"{header[1]:^12}", "|", f"{header[2]:^12}", "|", f"{header[3]:^12}")
    while count < n:
        cur_cell = track_data[cell_num]
        if idx >= len(cur_cell['x']):
            idx = 0
            cell_num += 1
            cur_cell = track_data[cell_num]

        print(f"{cell_num:<12}", "|", f"{cur_cell['x'][idx]:<12.2f}", "|", f"{cur_cell['y'][idx]:<12.2f}", "|", f"{cur_cell['t'][idx]:<12}")
        idx += 1
        count += 1

def sort_track_by_time(track_dict: Dict[str, Union[int, float]]) -> None:
    """
    Track data is often not ordered chronologically. 
    This function ensures that position and time data are chronologoically ordered.
    Acts in place.
    """
    i_sort = np.argsort(track_dict['t'])
    track_dict['t'] = list(np.array(track_dict['t'])[i_sort])
    track_dict['x'] = list(np.array(track_dict['x'])[i_sort])
    track_dict['y'] = list(np.array(track_dict['y'])[i_sort])

def dedrift(track_data: List[Dict[str, Union[int, float]]]) -> Dict[int, Tuple[float, float]]:
        """
        Mutates track_data  to dedrift all tracks using mean velocity at every 
        time point. Returns mean velocity at every time point.
        """
        velocities = {t: [] for cell in track_data for t in cell['t']}
        del velocities[min(velocities.keys())]
        
        for cell in track_data:
            for i in range(1,len(cell['t'])):
                
                dx = cell['x'][i] - cell['x'][i-1]
                dy = cell['y'][i] - cell['y'][i-1]
                
                velocities[cell['t'][i]].append((dx, dy))

        # get mean velocity at each available time point        
        vbar = {t: (sum(elt[0] for elt in velocities[t])/len(velocities[t]), 
                    sum(elt[1] for elt in velocities[t])/len(velocities[t])) 
                for t in velocities}
        
        # subtract summed drift velocities from each available timepoint
        for cell in track_data:
            for i,t in enumerate(cell['t']):
                if i == 0:
                    continue
                cell['x'][i] -= sum(vbar[T][0] for T in cell['t'][1:i+1])
                cell['y'][i] -= sum(vbar[T][1] for T in cell['t'][1:i+1])

        out = {'x': [], 'y': [], 't': []}
        for key, val in vbar.items():
            out['t'].append(int(key))
            out['x'].append(val[0])
            out['y'].append(val[1])
        sort_track_by_time(out)
        return out

def split_tracks(self) -> None:
    """
    Splits tracks where abrupt changes in direction occur.
    """
    
    new_tracks = []
    
    for i, cell in enumerate(self.tracks):
        
        X = cell['x']
        if len(X) <= 4:
            continue
        
        Y = cell['y']
        T = cell['t']
        
        x_brks = find_breakpoints(X, T)
        y_brks = find_breakpoints(Y, T)
        
        brks = x_brks.union(y_brks)
        
        removals = set()
        for idx in brks:
            if idx - 1 in brks or idx - 2 in brks:
                removals.add(idx)
        brks -= removals
        
        brks = sorted(brks, reverse=True)
        
        for idx in brks:
            assert len(cell['x'][idx + 2:]) > 1, f"{i=}"
            new_tracks.append({'x': cell['x'][idx + 2:], 
                                'y': cell['y'][idx + 2:], 
                                't': cell['t'][idx + 2:]})
            cell['t'] = cell['t'][:idx+1]
            cell['x'] = cell['x'][:idx+1]
            cell['y'] = cell['y'][:idx+1]
            
    self.tracks.extend(new_tracks) 