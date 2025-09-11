import csv
from typing import Dict, Union, List
from tqdm import tqdm
import os
from .track_manipulations import sort_track_by_time
from numpy.typing import NDArray
import numpy as np

def load_data_from_csv(filename: str) -> List[Dict[str, Union[float, int]]]:

    with open(os.path.expanduser(filename), newline='') as csvfile:
        
        csv_read = csv.reader(csvfile,delimiter=',')
        
        prev_cell = None
        cur_dict = {'x': [], 'y': [], 't': []}
        tracks = []

        header = next(csv_read)
        header = [elt.upper() for elt in header]
        assert "POSITION_X" in header, "File header must contain `POSITION_X` column"
        assert "POSITION_Y" in header, "File header must contain `POSITION_Y` column"
        assert "FRAME" in header, "File must contain `FRAME` column"
        assert "TRACK_ID" in header, "File must contain `TRACK_ID` column"

        counter = [0,0,0,0]
        for i, label in enumerate(header):
            if label == "POSITION_X":
                x_column = i
                counter[0] += 1
            elif label == "POSITION_Y":
                y_column = i
                counter[1] += 1
            elif label == "FRAME":
                t_column = i
                counter[2] += 1
            elif label == "TRACK_ID":
                cell_column = i
                counter[3] += 1

        assert counter == [1,1,1,1], "Position_x, Position_y, Frame, and Track_ID columns must each appear exactly once"
            
        for i, row in tqdm(enumerate(csv_read)):
            
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
            cur_dict['t'].append(int(float(t)))

            prev_cell = cur_cell
    
    sort_track_by_time(cur_dict)
    tracks.append(cur_dict) 

    return tracks

def save_track_data(track_data, filename: str, label=None) -> None:
    with open(os.path.expanduser(filename), 'w', newline = '') as csvfile:
        csv_write = csv.writer(csvfile)
        header = [''] * 4
        header[0] = 'TRACK_ID'
        header[1] = 'POSITION_X'
        header[2] = 'POSITION_Y'
        header[3] = 'FRAME'
        
        csv_write.writerow(header)
        cur_cell = 0
        for cell in tqdm(track_data):
            cur_cell += 1
            for i, x in enumerate(cell['x']):
                newrow = [''] * 4
                if label is None:
                    newrow[0] = str(cur_cell)
                else:
                    newrow[0] = label
                newrow[1] = x
                newrow[2] = cell['y'][i]
                newrow[3] = cell['t'][i]
                
                csv_write.writerow(newrow)

def save_MSD_data(filename: str, t: List[int], y: List[float], dy: List[float]) -> None:
    
    assert len(t) == len(y), "Input lists must have same length"
    assert len(t) == len(dy), "Input lists must have same length"

    with open(os.path.expanduser(filename), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Time Lag", "MSD", "Error_MSD"])
        for i, elt in tqdm(enumerate(t)):
            writer.writerow([str(elt), str(y[i]), str(dy[i])])

def save_turning_angles_data(filename: str, data: NDArray[np.float64], dt: int) -> None:
    with open(os.path.expanduser(filename), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([f"Turning angles [deg], delta_t={dt} frames"])
        for elt in tqdm(data):
            writer.writerow([str(elt)])