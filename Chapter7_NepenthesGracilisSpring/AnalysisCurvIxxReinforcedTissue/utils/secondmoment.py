""" Calculate the second moment of area for a tiff stack

Python module to load in every .tiff-file form a CT Scan stack, calculate the
second moment of area and area for each file and export as a single .csv-file.

This module requires the following packages:
    * pathlib
    * numpy
    * PIL
    * matplotlib
    * scipy
    * logging

This module contains the following functions:
    * load_tif - load single image into np.array
    * load_res - load the resolution and unit for image stack
    * plot_histogram - plot the histogram of a single tiff image with threshold
    * appl_thresh - set all pixels below threshold to zero and above to one
    * appl_thresh_grey - only set pixels below threshold to zero and above to
                         fraction of one
    * find_xaxis - find the mean x axis of all non black pixels in a image
    * calc_I - calculate the second moment of area for a single image
    * calc_I_stack - loop through whole image stack to calculate I with all
                     other functions above

Exceptions are raised in the following functions:
    * load_res - IOError: if .xtekct file cannot be found
    * calc_I_stack - IOError: if folder with image stack does not exist
                   - IOError: if folder is empty
"""

import numpy as np
from pathlib import Path
from PIL import Image
import matplotlib.pyplot as plt
from scipy import ndimage
from skimage.segmentation import flood_fill
import logging


def load_tif(filename, path=None):
    """Load a single tif file and return image as np.array"""

    logging.debug("Loading tif-file: {}".format(filename))

    if path:
        filename = path / Path(filename)

    im = Image.open(filename)
    image = np.array(im)
    im.close()

    return image


def load_res(pathImport):
    """Load .xtekct file and read the voxel size and unit from it
    
    IOError: File not found
    """

    # Check if file can be found
    if not pathImport.is_file():
        raise IOError("File not found: {}!".format(pathImport))
        return

    # Extract resolution and the unit form file
    with open(pathImport) as f:
        for line in f:
            if line[:10] == "VoxelSizeX":
                resolution = float(line[12:])
            if line[:5] == "Units":
                unit = line[6:-1]

    logging.info("Resolution and unit extracted from: {}".format(pathImport))

    return resolution, unit


def plot_histogram(filename, image, thresh=None, path=None):
    """
    Plot the histogram for a single image array and mark the threshold line in
    the plot, use path and filename to export histogram
    """

    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1)

    # Histogram cut at 500 so high peak from black pixels is cut off and more
    # interesting area scaled up
    ax1.hist(image.reshape(image.size, 1), bins=256)
    ax1.set_ylim(0, 500)
    ax1.set_ylabel('No. of pixels', size=14)
    # Line for threshold with annotation
    if thresh:
        ax1.axvline(x=thresh, color='orange')
        ax1.annotate('Threshold: {}'.format(thresh), (thresh, 432),
                     xytext=(0.52, 0.9), textcoords='axes fraction',
                     arrowprops=dict(width=0.001, headwidth=5,
                     headlength=5, facecolor='black', shrink=0.05),
                     fontsize=10, horizontalalignment='right',
                     verticalalignment='top')

    # Histogram with log-scale and no cut off
    ax2.hist(image.reshape(image.size, 1), bins=256, log=True)
    ax2.set_ylabel('Log-scale: No. of pixels', size=14)
    ax2.set_xlabel('black - - - - - -  grey - - - - -  white', size=14)
    if thresh:
        ax2.axvline(x=thresh, color='orange')

    fig.suptitle('Histogram image {}'.format(filename), size=16)

    # Save figure
    if not path:
        path = Path.cwd()
    figname = path / (filename[:-4] + '_histogram.pdf')
    fig.savefig(figname)
    plt.close()

    logging.info("Saved histogram plot: {}".format(figname))

    return None


def appl_tresh(image, thresh):
    """Apply a threshold to the image to make it black white"""

    bwimage = np.zeros(image.shape, dtype=int)
    np.copyto(bwimage, image)

    bwimage[bwimage < thresh] = 0
    bwimage[bwimage >= thresh] = 1

    logging.debug("Black white threshold {} applied to image.".format(thresh))

    return bwimage


def fill_holes(image):
    """Fill all holes within structure to have real area and I"""

    # Fill everything but the holes
    invertedholes = flood_fill(image, (0, 0), 1, connectivity=1)

    # For ring cross sections the inside of the ring needs to be filled
    if np.sum(invertedholes) < 0.95 * np.size(invertedholes):
        centrex = int(image.shape[0] / 2)
        centrey = int(image.shape[1] / 2)
        invertedholes = flood_fill(invertedholes, (centrex, centrey), 1,
                                   connectivity=1)

        # If there is still more than 5% black pixels, 
        if np.sum(invertedholes) < 0.95 * np.size(invertedholes):
            centrex = int(image.shape[0] * 0.8)
            centrey = int(image.shape[1] / 2)
            invertedholes = flood_fill(invertedholes, (centrex, centrey), 1,
                                   connectivity=1)

            if np.sum(invertedholes) < 0.95 * np.size(invertedholes):
                centrex = int(image.shape[0] * 0.3)
                centrey = int(image.shape[1] / 2)
                invertedholes = flood_fill(invertedholes, (centrex, centrey), 1,
                                   connectivity=1)

                if np.sum(invertedholes) < 0.95 * np.size(invertedholes):
                    logging.warning("Holes in cross section could not be filled!")
                    holes = (-1 * invertedholes) + 1
                    filledimage = image + holes
                    return filledimage

    holes = (-1 * invertedholes) + 1
    filledimage = image + holes

    logging.debug("Holes filled successfully.")

    return filledimage


def save_bwimage(filename, image, path=None):
    """Export the black white image when holes are filled"""

    filename = str(filename)[:-4] + "_bwexport.png"

    if path:
        filename = path / filename

    Image.fromarray((image* 255).astype(np.uint8)).save(filename)

    logging.info("Black white cross section saved as: {}".format(filename))

    return


def appl_tresh_grey(image, thresh, maxval=None):
    """Apply a threshold to the image to make it black below threshold and
    greyvalues in percentage of how white they are
    """

    if not maxval:
        maxval = 65535

    greyimage = np.zeros(image.shape, dtype=int)
    np.copyto(greyimage, image)

    greyimage[greyimage < thresh] = 0
    greyimage = greyimage / maxval

    logging.debug("Grey threshold {} applied to image.".format(thresh))

    return greyimage


def find_xaxis(image):
    """Find the axis of an mage by choosing the centre of gravity"""

    axis = ndimage.measurements.center_of_mass(image)
    logging.debug("Axis of gravity found.")

    return round(axis[0], 0)


def calc_I(image, axis, resolution):
    """
    Calculate the second moment of area for an image (np.array)
    around a defined axis
    I_x1 = b*h**3 / 12
    I_x1' = I_x1 + A_1 * d_1**2
    I_x' = n * I_x + sum(n * A * d**2)
    """
    # Variables per pixel - Ix, A, n
    i_x_pixel = resolution**4 / 12
    area_pixel = resolution**2
    n_pixel = image.sum()

    # Variables based on distance to axis - d, n
    distances = (np.arange(0, len(image), 1) - axis) * resolution
    n_row = image.sum(axis=1)

    # n * A * d**2
    area_dsq_row = n_row * area_pixel * distances**2

    # I_x' of whole array
    secondmoment = n_pixel * i_x_pixel + area_dsq_row.sum()

    logging.debug("Second moment calculated.")

    return secondmoment

def load_conversion(basename):
    filepath = Path.cwd() / "ConversionScans.csv"
    
    if not filepath.is_file():
        raise IOError("Files not found: {}!".format(filepath))
        return np.linspace(0,-15,16), np.zeros(16)

    nrows = 16
    if int(basename[-1]) == 3:
        nrows = 14
    if int(basename[-1]) == 4:
        nrows = 15
    if int(basename[-1]) == 5:
        nrows = 15
    if int(basename[-1]) == 6:
        nrows = 15
    indices, slices = np.loadtxt(filepath, delimiter=',', skiprows=1, 
                                 usecols=(0,int(basename[-1])), unpack=True,
                                 dtype='i', max_rows=nrows)

    return indices, slices


def calc_I_stack(pathImport, thresh, thresh_venation, resolution, unit, pathExport=None):
    """
    Calculate the second moment of area for an image stack from pathImport and
    export as a .csv file

    IO Error: Folder not found or no files in folder
    """

    # Check that folder exists, if not throw exception and return
    if not pathImport.is_dir():
        raise IOError("Folder not found: {}!".format(pathImport))
        return None

    # If no export path is given, use import path and add export folder
    if not pathExport:
        pathExport = pathImport / "Export"
        logging.info("Export path created: {}".format(pathExport))
        if not pathExport.is_dir():
            pathExport.mkdir()

    pathExportNew = pathExport / (pathImport.name[:-11] + "Export")
    if not pathExportNew.is_dir():
            pathExportNew.mkdir()

    # Make list of all files in folder
    files = sorted([x.name for x in pathImport.glob('*.tif')])
    logging.info("{} files found in folder for I calc.".format(len(files)))

    # If no files are found, throw exception and return
    if not files:
        raise IOError("No tif-files found in folder: {}!".format(pathImport))
        return None

    # Create variable for saving I
    i_total_bw = []
    i_total_grey = []
    area_total = []
    area_veins = []

    if isinstance(thresh_venation, list):
        thresh_venation = thresh_venation[int(pathImport.name[12])-1]

    # indices and slices to save out images only at index numbers
    indices, slices = load_conversion(pathImport.name[12])
    maxlid = int(files[-1][-7:-4])
    lidindex = np.array([0, 1, 2, 3, 4, 5])
    lidslice = len(files) - np.flip(np.round(lidindex * maxlid / 5).astype(int))

    # Loop through all files in list
    for i, file in enumerate(files):

        # Load image
        image = load_tif(file, pathImport)

        # Plot histogram for every 100th file
        #if i % 100 == 0:
        #    plot_histogram(file, image, thresh, pathExportNew)

        # Apply the threshold for every file and fill holes in bw image
        bwimage = appl_tresh(image, thresh)
        #logging.info("Filling holes for file {}.".format(file))
        bwimage = fill_holes(bwimage)
        # Plot png with filled holes as black white image
        # Uncomment if loop for exporting all frames
        #if i % 100 == 0:
        #save_bwimage(file, bwimage, pathExportNew)

        greyimage = appl_tresh_grey(image, thresh)

        veinimage = appl_tresh(image, thresh_venation)
        
        if -(i+1) in slices:
            filename = "Veins_Index{}xxxxx".format(indices[np.where(slices == -(i+1))])
            save_bwimage(filename, veinimage, pathExportNew)
            filename = "BW_Index{}xxxxx".format(indices[np.where(slices == -(i+1))])
            save_bwimage(filename, bwimage, pathExportNew)
        
        if i in lidslice:
            filename = "Veins_Index{}xxxx".format(lidindex[np.where(lidslice == i)])
            save_bwimage(filename, veinimage, pathExportNew)
            filename = "BW_Index{}xxxx".format(lidindex[np.where(lidslice == i)])
            save_bwimage(filename, bwimage, pathExportNew)
        # Find axis in first file
        #if i == 0:
        axis = find_xaxis(bwimage)
        axisgrey = find_xaxis(greyimage)

        # Calculate the area of all white pixels and add to list of areas
        area_total.append(bwimage.sum() * resolution**2)
        area_veins.append(veinimage.sum() * resolution**2)

        # Calculate I and add to list of I
        i_total_bw.append(calc_I(bwimage, axis, resolution))
        i_total_grey.append(calc_I(greyimage, axisgrey, resolution))

    # Create index row and stack data into array
    index = [int(x.split("_")[-1][:-4]) for x in files]
    data = np.column_stack((index, i_total_bw, i_total_grey, area_total, area_veins))
    # Create header and filename
    header = "Image, Ix bw in {0}^4, Ix grey in {0}^4, area in {0}^2, vein area in {0}^2".format(
                                                                      unit)
    filename = pathExport / (pathImport.name[:-11] + "secondmoment.csv")
    # Save list as .csv file in export folder
    np.savetxt(filename, data, delimiter=',', header=header)

    logging.info("Saved I calculation in file: {}".format(filename))

    return None
