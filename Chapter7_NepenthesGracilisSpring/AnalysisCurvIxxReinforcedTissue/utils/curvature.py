""" Calculate curvature

Python module to load in a .csv-file with x and y coordinates, caluclate the
curvature of the graph at every coordinate pair and save the curvature into a
.csv-file.

This module requires the following packages:
    * pathlib
    * numpy
    * logging

This module contains the following functions:
    * load_file - save index and x, y coordinates (cols 0,2,3) into an array
    * calc_curvature - curvature = (x'y" - x"y') / (x'**2 + y'**2)**1.5
    * export_curvature - export index, x, y coordinates and curvature to file
"""

from pathlib import Path
import numpy as np
import logging

def load_file(filename, path=None):
    """Load in file and return array with index, x and y coordinate"""

    logging.info("Loading curvature file: {}".format(filename))

    if path:
        filename = path / Path(filename)

    coords = np.loadtxt(filename, delimiter=',', skiprows=3,
                        usecols=(0, 2, 3))
    return coords


def calc_curvature(coords):
    """Calculate the curvature of a set of x,y coordinates and return 1D array

    The general parametrisation equation:
    curvature = (x"y' - x'y") / (x'**2 + y'**2)**1.5

    The function uses numpy.gradient() for solving the equation.
    """

    # Calculate all necessary gradients from coordinates
    dx = np.gradient(coords[:, 1])  # x'
    dy = np.gradient(coords[:, 2])  # y'
    d2x = np.gradient(dx)           # x"
    d2y = np.gradient(dy)           # y"

    # Calculate curvature = (x'y" - x"y') / (x'**2 + y'**2)**1.5
    curv = (dx*d2y - d2x*dy) / (dx*dx + dy*dy)**1.5

    # Reverse the order of the coordinates and redo curvature calculation
    # Doesn't change anything!!!
    #revcoords = coords[::-1]

    #revdx = np.gradient(revcoords[:, 1])  # x'
    #revdy = np.gradient(revcoords[:, 2])  # y'
    #revd2x = np.gradient(revdx)           # x"
    #revd2y = np.gradient(revdy)           # y"

    #revcurv = -(revdx*revd2y - revd2x*revdy) / (revdx*revdx + revdy*revdy)**1.5

    # Turn order of revcurv around to match original coordinates
    #revcurv = revcurv[::-1]

    # Non lagged curvature is mean of curvature and reverse calculated curvature
    #finalcurv = (curv + revcurv) / 2

    #print(curv, revcurv, finalcurv)

    logging.info("Curvature calculated: {}...".format(curv[0:2]))

    return curv


def export_curvature(filename, coords, curv, path=None):
    """Export index, x, y coordinates and curvature into .csv-file"""

    # Set export path and filename
    if path:
        filepath = path / Path(filename + "_curvature.csv")
    else:
        filepath = Path(filename + "_curvature.csv")

    # Merge data and save into file
    header = "Index, x, y, curvature"
    curv = curv.reshape((len(curv), 1))
    data = np.hstack((coords, curv))
    np.savetxt(filepath, data, delimiter=',', header=header)
    logging.info("Curvature exported to: {}".format(filepath))

    return None
