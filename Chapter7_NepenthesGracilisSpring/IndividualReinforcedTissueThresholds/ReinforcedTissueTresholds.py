import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy import ndimage


def load_tif(filename, path=None):
    """Load a single tif file and return image as np.array"""

    logging.debug("Loading tif-file: {}".format(filename))

    if path:
        filename = path / Path(filename)

    im = Image.open(filename)
    image = np.array(im)
    im.close()

    return image


def appl_tresh(image, thresh):
    """Apply a threshold to the image to make separate the object"""

    ctobject = np.zeros(image.shape, dtype=int)
    np.copyto(ctobject, image)

    ctobject[ctobject < thresh] = 0

    logging.debug("Background threshold {} applied to image.".format(thresh))

    return ctobject


def plot_histogram(filename, image, thresh=None, path=None):
    """
    Plot the histogram for a single image array and mark the threshold line in
    the plot, use path and filename to export histogram
    """

    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1)

    # Histogram cut at 500 so high peak from black pixels is cut off and more
    # interesting area scaled up
    ax1.hist(image.reshape(image.size, 1), bins=65535-thresh, range=(thresh,65535))
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
    data = ax2.hist(image.reshape(image.size, 1), bins=65535-thresh, range=(thresh,65535), density=True, cumulative=True)
    ax2.set_ylabel('Cumulative No. of pixels', size=14)
    ax2.set_xlabel('black - - - - - -  grey - - - - -  white', size=14)
    if thresh:
        ax2.axvline(x=thresh, color='orange')

    fig.suptitle('Histogram image {}'.format(filename[:-10]), size=16)

    # Save figure
    if not path:
        path = Path.cwd()
    figname = path / (filename[:-9] + '_histogram.pdf')
    fig.savefig(figname)
    plt.close()

    logging.info("Saved histogram plot: {}".format(figname))

    percentage = data[0]
    greyvals = data[1]

    return percentage, greyvals

def save_distribution(filename, percentage, greyvals, path):

    # Set export path and filename
    if path:
        filepath = path / Path(filename[:-9] + "_distribution.csv")
    else:
        filepath = Path(filename + "_distribution.csv")


    percentage_selection = np.linspace(0,1,101)
    indices = np.searchsorted(percentage, percentage_selection)
    greyval_selection = greyvals[indices]
    
    # Merge data and save into file
    header = "Percentage, Greyvalue"
    data = np.vstack((percentage_selection*100 , greyval_selection)).T.astype(int)
    np.savetxt(filepath, data, delimiter=',', fmt='%5i' ,header=header)
    logging.info("Distribution exported to: {}".format(filepath))

    return


def find_distribution(pathImport, thresh, pathExport=None):


    # Check that folder exists, if not throw exception and return
    if not pathImport.is_dir():
        raise IOError("Folder not found: {}!".format(pathImport))
        return None

    # If no export path is given, use import path and add export folder
    if not pathExport:
        pathExport = Path.cwd()
        logging.info("Export path created: {}".format(pathExport))
        if not pathExport.is_dir():
            pathExport.mkdir()

    # Make list of all files in folder
    files = sorted([x.name for x in pathImport.glob('*.tif')])
    logging.info("{} files found in folder for histogram and grey value distribution.".format(len(files)))

    # If no files are found, throw exception and return
    if not files:
        raise IOError("No tif-files found in folder: {}!".format(pathImport))
        return None

    # Loop through all files in list
    for i,file in enumerate(files):

        # Load image
        image = load_tif(file, pathImport)

        # If it is first file create a variable to hold all data in array
        if i == 0:
            imsize = image.size
            imshape = image.shape
            stack = np.zeros((len(files), imshape[0], imshape[1]))

        # Lid images need to be adjusted so they have same shape as body images
        if image.size != imsize:
            image = image.reshape(image.size)
            image = np.pad(image, (0, imsize-image.size))
            image = image.reshape(imshape)

        # Add grey values of object to image stack, background is np.nan
        stack[i] = appl_tresh(image, thresh)

    # Plot histogram of ct object stack
    percentage, greyvals = plot_histogram(file, stack, thresh, pathExport)

    save_distribution(file, percentage, greyvals, pathExport)

    return


def main():
    # Setup logger
    logging.basicConfig(filename="DistributionDensity.log", level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info('Start...')

    # Folder
    pathImport = Path.cwd() / "N_gracilis_p6_middle[tif]"

    # Threshold for separating background and object
    thresh = 10000

    # Find distribution in whole image stack
    find_distribution(pathImport, thresh)

    logging.info('...script finished!')


if __name__ == '__main__':
    main()
