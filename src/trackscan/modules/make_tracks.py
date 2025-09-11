from .linking_functions import get_linking_candidates, choose_linking_partners, link_partners
from .io_utils import load_data_from_csv, save_track_data, save_MSD_data, save_turning_angles_data
from .track_manipulations import scale_space, print_first_n, dedrift, split_tracks
from .msd_analysis import msd, fit_msd
from .analysis_utils import turning_angle_analysis, get_ensemble_mean_speed
from .visualization_utils import plot_msd, plot_turning_angles
import argparse
import os
import cmd
import shlex

class Tracks(cmd.Cmd):
    intro = "Welcome to the trackscan command line interface for cell track post-processing and analysis.\n" \
    "Type help or ? to list commands.\n"
    prompt = "(trackscan) "

    def __init__(self, filename: str):
        """
        Tracks object has filename and tracks attributes. 
        
        filename: file containing track data (must be .csv)
        tracks: list of dicts with keys 'x', 'y', 't' and values list of x, y, 
                and t data. Each dict represents one cell's track.
        """
        super().__init__()
        self.tracks = load_data_from_csv(filename)
        self.filename = filename 
            
    def do_scale_space(self, arg: str) -> None:
        """Scale x and y data by a constant factor. Type scale_space -h to see available options."""

        parser = argparse.ArgumentParser(prog="scale_space", add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("-c", "--constant", default=1, type=float, help="Scale x and y data by a constant factor")
        parser.add_argument("-h", "--help", action="help", help="Show this message")

        try:
            args = parser.parse_args(shlex.split(arg))
        except SystemExit:
            return
        
        scale_space(self.tracks, args.constant)

    def do_show(self, arg: str) -> None:
        """Print the first n lines of track data. Type show -h to see available options."""
        parser = argparse.ArgumentParser(prog="show", add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("-n", "--n", default=5, type=int, help="Print the first n lines of track data.")
        parser.add_argument("-h", "--help", action="help", help="Show this message")

        try:
            args = parser.parse_args(shlex.split(arg))
        except SystemExit:
            return
        
        print_first_n(self.tracks, args.n)
    
    def do_save(self, arg: str) -> None:
        """Save cell track data in CSV format. Type save -h to see available options."""

        root_name, ext = os.path.splitext(self.filename)
        default_filename = root_name + "PROCESSED" + ext

        parser = argparse.ArgumentParser(prog="save", add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("-o", "--out", default=default_filename, type=str, help="Path to CSV file to save data")
        parser.add_argument("-h", "--help", action="help", help="Show this message")
        
        try:
            args = parser.parse_args(shlex.split(arg))
        except SystemExit:
            return

        save_track_data(self.tracks, args.output)
                    
    def do_dedrift(self, arg: str) -> None:
        """Dedrift track data by cell-averaged velocity at each timepoint. Type dedrift -h to see available options."""

        parser = argparse.ArgumentParser(prog="dedrift", add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("-s", "--save", type=bool, default=False, help="Whether to save the drift velocities")
        parser.add_argument("-o", "--out", type=str, default="", help="Path to CSV file to save drift velocities if --save is True")
        parser.add_argument("-h", "--help", action="help", help="Show this message")

        try:
            args = parser.parse_args(shlex.split(arg))
        except SystemExit:
            return
        
        drift_velocities = dedrift(self.tracks)
        if args.save:
            save_track_data([drift_velocities], args.out, label="DRIFT")
            
    def do_link(self, arg:str) -> None:
        """
        Links tracks together by identifying candidates and choosing best one.
        """

        parser = argparse.ArgumentParser(prog="link", add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("-t", "--time", type=int, default=3, help="Max time gap to consider linking two tracks")
        parser.add_argument("d", "--dist", type=float, default=10, help="Max distance to consider linking two tracks")
        parser.add_argument("-h", "--help", action="help", help="Show this message")

        try:
            args = parser.parse_args(shlex.split(arg))
        except SystemExit:
            return
        
        linking_candidates = get_linking_candidates(self.tracks, args.time, args.dist)
        link_dict = choose_linking_partners(linking_candidates, self.tracks)
        link_partners(link_dict, self.tracks)
    
    def do_split(self, arg: str) -> None:
        """
        Splits tracks where abrupt changes in direction occur.
        """
        parser = argparse.ArgumentParser(prog="split", add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("-h", "--help", action="help", help="Show this message")

        try:
            args = parser.parse_args(shlex.split(arg))
        except SystemExit:
            return
        
        split_tracks(self.tracks)

    def do_correct_artifacts(self, arg:str) -> None:
        """
        Fix tracking artifacts by first splitting tracks at locations of artifacts, then linking tracks correctly
        """ 

        parser = argparse.ArgumentParser(prog="correct_artifacts", add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("-t", "--time", type=int, default=3, help="Max time gap to consider linking two tracks")
        parser.add_argument("d", "--dist", type=float, default=10, help="Max distance to consider linking two tracks")
        parser.add_argument("-h", "--help", action="help", help="Show this message")

        try:
            args = parser.parse_args(shlex.split(arg))
        except SystemExit:
            return
        
        print("Splitting tracks...")
        split_tracks(self.tracks)
        print("Linking tracks...")
        linking_candidates = get_linking_candidates(self.tracks, args.time, args.dist)
        link_dict = choose_linking_partners(linking_candidates, self.tracks)
        link_partners(link_dict, self.tracks)

    def do_MSD_analysis(self, arg:str) -> None:
        """Compute mean squared displacements"""

        root_name, ext = os.path.splitext(self.filename)
        default_filename = root_name + "_MSD" + ext

        parser = argparse.ArgumentParser(prog="MSD_analysis", add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("-o", "--out", type=str, default=default_filename, help="Path to CSV file to save MSD data")
        parser.add_argument("-p", "--plot", type=bool, default=False, help="Whether to output a plot of MSD data")
        parser.add_argument("-h", "--help", action="help", help="Show this message")

        try:
            args = parser.parse_args(shlex.split(arg))
        except SystemExit:
            return
        
        t, y, dy = msd(self.tracks)
        save_MSD_data(args.out, t, y, dy)

        if args.plot:
            root_name, ext = os.path.splitext(args.out)
            save_to = root_name + "_MSD_plot" + ".png"

            (a, d), (a_err, D_err), r2 = fit_msd(t, y, (1, len(t)//2))
            plot_msd(save_to, t, y, dy, a, d, a_err, D_err, r2)

    def do_turning_angle(self, arg:str) -> None:
        """
        Compute turning angles across a specified time
        """

        root_name, ext = os.path.splitext(self.filename)
        default_filename = root_name + "_turning_angles" + ext

        parser = argparse.ArgumentParser(prog="turning_angles", add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("-t", "--time", type=int, default=10, help="Time lag to calculate turning angles")
        parser.add_argument("-o", "--out", type=str, default=default_filename, help="Path to CSV file to save turning angle data")
        parser.add_argument("-p", "--plot", type=bool, default=False, help="Whether to output a plot of turning angle data")
        parser.add_argument("-h", "--help", action="help", help="Show this message")

        try:
            args = parser.parse_args(shlex.split(arg))
        except SystemExit:
            return
        
        if args.time <= 0:
            print("turning_angle: time must be positive")
            return
        
        turning_angles = turning_angle_analysis(self.tracks, args.time)
        save_turning_angles_data(args.out, turning_angles, args.time)

        if args.plot:
            root_name, ext = os.path.splitext(args.out)
            save_to = root_name + "_turning_angles_plot" + ".png"
            plot_turning_angles(save_to, turning_angles)    
        

    
    def do_mean_speed(self, arg:str) -> None:
        """Measure and print ensemble mean cell speed"""

        parser = argparse.ArgumentParser(prog="mean_speed", add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("-h", "--help", action="help", help="Show this message")

        try:
            args = parser.parse_args(shlex.split(arg))
        except SystemExit:
            return
        print(get_ensemble_mean_speed(self.tracks))
        

    def do_exit(self, arg:str) -> bool:
        """Exit the program"""
        return True
    
    
    def get_cell_mean_stepsize(self, cell):
        dists = []
        for i in range(1, len(cell['t'])):
            cur_x = cell['x'][i]
            prev_x = cell['x'][i-1]
            
            cur_y = cell['y'][i]
            prev_y = cell['y'][i-1]
            dists.append(((cur_x - prev_x)**2 + (cur_y - prev_y)**2)**(1/2))
            
        return sum(dists)/len(dists)