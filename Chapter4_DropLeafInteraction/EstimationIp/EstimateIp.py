#Ip calc
from pathlib import Path
import csv
import glob
from PIL import Image
import numpy as np
from scipy import ndimage
from skimage.segmentation import flood_fill


def load_tif(filename, path=None):
    """Load a single tif file and return image as np.array"""

    if path:
        filename = path / Path(filename)

    im = Image.open(filename)
    image = np.array(im)
    im.close()

    return image

def appl_tresh(image, thresh=100):
    """Apply a threshold to the image to make it black white"""

    bwimage = np.zeros(image.shape, dtype=int)
    np.copyto(bwimage, image)

    bwimage[bwimage < thresh] = 1
    bwimage[bwimage >= thresh] = 0

    return bwimage

def fill_holes(image):
    """Fill all holes within structure to have real area and I"""
    # Fill everything but the holes
    invertedholes = flood_fill(image, (0, 0), 1, connectivity=1)

    holes = (-1 * invertedholes) + 1
    filledimage = image + holes

    return filledimage

def find_xaxis(image):
    """Find the axis of an mage by choosing the centre of gravity"""

    axis = ndimage.center_of_mass(image)

    return round(axis[0], 0)

def calc_I(image, axis, resolution, mass):
    """
    Calculate the polar moment of area for an image (np.array)
    around a defined axis based on Bhosale 2020
    I_p = mass * distance**2
    """
    # Variables per pixel - M, n
    n_pixel = image.sum()
    mass_pixel = mass/ n_pixel

    # Variables based on distance to axis - d, n
    distances = (np.arange(0, len(image), 1) - axis) * resolution
    n_row = image.sum(axis=1)

    # n * M * d**2
    mass_dsq_row = n_row * mass_pixel * distances**2

    # I_p of whole array
    polarmoment = mass_dsq_row.sum()

    return polarmoment

def save_bwimage(filename, image, axis, path=None):
    """Export the black white image when holes are filled"""

    filename = str(filename)[:-4] + "_bwexport.png"
    #print(axis)
    image[int(axis), :] = 0

    if path:
        filename = path / filename

    Image.fromarray((image* 255).astype(np.uint8)).save(filename)

    return


def main():

    pathData = Path.cwd() / "AreaOutlines"
    thresh = 100

    # Check if data path exists
    if not pathData.is_dir():
        print("No data found to analyse!")
        return

    # Load resolution
    species = np.genfromtxt("Resolution.csv", delimiter=',', skip_header = 1,
        usecols=0,  dtype=str)
    resolution, masses = np.genfromtxt("Resolution.csv", delimiter=',', skip_header = 1,
        usecols=(1,2), unpack=True)

    # Convert resolution from pixel/mm to mm/pixel
    resolution = 1/resolution

    files = sorted([x.name for x in pathData.glob('*.tif')])

    area_total = []
    i_total_bw = []

    for i,file in enumerate(files):

        # Find correct resolution
        if file[:4] == species[i][1:5]:
            res = resolution[i]
            mass = masses[i]

        # Load image
        image = load_tif(file, pathData)
        
        # Apply the threshold for every file and fill holes in bw image
        bwimage = appl_tresh(image, thresh)
        #logging.info("Filling holes for file {}.".format(file))
        bwimage = fill_holes(bwimage)

        axis = find_xaxis(bwimage)
         # Calculate the area of all white pixels and add to list of areas
        area_total.append(bwimage.sum() * res**2)

        # Calculate Ip and add to list of I
        i_total_bw.append(calc_I(bwimage, axis, res, mass))

        save_bwimage(file, bwimage, axis, (Path.cwd() / "IpOutput"))


    data = zip(species,i_total_bw,area_total)

    #header = "Species, Ip in g*mm^2, area in mm^2"
    filename = Path.cwd() / "IpOutput" / "IpVariables.csv"

    # Write data to the CSV file
    with open(filename, 'w', newline='') as newfile:
        # Create a CSV writer object
        csv_writer = csv.writer(newfile)

        # Write the header (if needed)
        csv_writer.writerow(['Species', 'Ip in g*mm^2', 'area in mm^2'])

        # Write the data
        csv_writer.writerows(data)


if __name__ == '__main__':
    main()
