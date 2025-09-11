import argparse
from .modules.make_tracks import Tracks

def main():
    parser = argparse.ArgumentParser(prog="trackscan", description="Cell track post-processing and analysis")
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to CSV file containing tracking data.")
    args = parser.parse_args()

    tracks = Tracks(args.input)
    tracks.cmdloop()

if __name__ == '__main__':
    main()