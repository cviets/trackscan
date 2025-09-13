<h2 align="center">TrackScan</h2>
A Python package for post-processing and analyzing cell track data. Includes functionality for correcting automated tracking artifacts, measuring mean squared displacements, and measuring cell turning angles.

## Introduction
Trackscan is a Python package that performs post-processing and analysis on 2D cell tracking data by first addressing two types of automated tracking artifacts: splitting and switching artifacts. Splitting artifacts arise from momentary failure in cell detection, e.g., if the cell briefly exits the focal plane, while splitting artifacts occur when two cells' tracks are confused with each other. Both of these artifact types obscure the long-term behavior of individual cells and inject large amounts of inaccuracy into common motility measurements, and must therefore be removed prior to analysis. Trackscan removes automated tracking artifacts by purposefully splitting tracks, then re-linking them in the correct configuration. Specifically, trackscan performs polynomial regression on moving windows of each track's $x(t)$ and $y(t)$ trajectories to determine where to place track splits and how to optimally mend the splits. Trackscan also contains functionality for track de-drifting and measurement of cell motility parameters such as mean squared displacement, turning angle distribution, and cell speed.

## Cite
```bibtex
@inbook{viets_measuring_2025,
author = {Viets, Chris and Stevens, Corey A.},
editor = {Brockhausen, Inka},
title = {Measuring and Analyzing Bacterial Movement in Mucus},
bookTitle = {Dynamics of Bacteria-Mucus Interactions},
year = {2025},
publisher = {Springer US},
address = {New York, NY},
pages = {187--197},
isbn = {978-1-0716-4627-4},
doi = {10.1007/978-1-0716-4627-4_16},
url = {https://doi.org/10.1007/978-1-0716-4627-4_16}
}

```

## Installation
```
pip install trackscan
```

## Usage
The first step to using `trackscan` is to read in a CSV file containing (note that the CSV file must contain columns labeled "Position_X", "Position_Y", "Frame", and "Track_ID -- not case-sensitive).
```
trackscan -i /path/to/track_data.csv
```
This command launches an interactive shell where the track data can be manipulated and measured. Once the interactive shell has appeared, simply type `?` to view the available commands. For example, the interactive shell contains commands to de-drift track data, correct artifacts arising from automated tracking, and measure mean squared displacement, turning angles, or mean cell speed.
