"""
TODO: Write comments/ Logging and check pep8 guidelines with linter
TODO: Ask for threshold
TODO: Replace get explanation if with dictionary
TODO: Add no fit for correction, but smoothed real values
"""


from pathlib import Path
import numpy as np
from PIL import Image
import logging
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
import random
import re

# Setup logger
logger = logging.getLogger('PrepareCTScan')
logger.setLevel(logging.INFO)

logfile = logging.FileHandler('PrepareCTScan.log')
logfile.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
logfile.setFormatter(formatter)

logger.addHandler(logfile)


# Help functions
def get_explanation(length):
    """Explanation text for using script interactively."""

    if length == "Long":
        text = """\n********************************************************\n
Interactive script to prepare CT scans for further analysis.\n
The standard steps in this script are:
    1) Separate the calibration object and object from the original scan images
    2) Find the shift in grey values over calibration object height and fit correction
    3) Apply correction to object images
    4) Clean chosen object images (from other objects in scan)
    5) Crop object images to minimal size\n
Usually the steps are performed in numerical order. The script stops after
each step and asks for the next step to be performed leaving the user time
to check if the previous step has been performed correctly or has to be
redone.\n
The script can be stopped at each break of the script by typing
"quit".\n
Which step of the analysis should be performed first? (Options: 1-5, quit)\n
\033[33m>>> \033[00m"""
    elif length == "Short":
        text = """\nSteps:
    1) Separate the calibration object and object from the original scan images
    2) Find the shift in grey values over calibration object height and fit correction
    3) Apply correction to object images
    4) Clean chosen object images (from other objects in scan)
    5) Crop object images to minimal size\n
Which step of the analysis should be performed? (Options: 1-5, quit)\n
\033[33m>>> \033[00m"""
    elif length == "Fit":
        text = """Which fit should be performed? 
    1) Linear fit or
    2) Exponential fit\n
\033[33m>>> \033[00m"""
    elif length == "Lower":
        text = """\nSet the lower boundary for the fit. 
    (In most scans the lowest 200 images are not as important and would distort
    the correction for the major part of the scan. This is due to the CT scanner
    not working well in the lower few cm of the scanning volume.)\n
Type in the lowest image number to be included.\n
\033[33m>>> \033[00m"""
    elif length == "Fit Error":
        text = """
    \033[36m'Don't Panic'   -- Douglas Adams\033[00m\n
Input not one of the options. Please try again or quit.\n
\033[33m>>> \033[00m"""
    elif length == "Scaling":
        text = """\nShould the correction include a general grey value scaling?
    (Usually a value in the range of the mean value of the calibration object
    should be used to avoid clipping, unless the grey values of the object and
    calibration object differ a lot. The mean of the calibration object can be
    found in the Correction_Fit_Plot.png)\n
Type in the grey value that should be used for scaling (int < 65535).\n
\033[33m>>> \033[00m"""
    elif length == "Cleaning":
        text = """\nWhich mask should be used for the cleaning process? It has
to be a black and white .tif- image that is in the main folder.\n
\033[33m>>> \033[00m"""
    elif length == "Range":
        text = """\nWhich images should be cleaned? Name lowest and highest
image to be cleaned with previously named mask, separated by '-'.\n
\033[33m>>> \033[00m"""
    elif length == "Cropping":
        text = """\nEither choose a threshold value to automatically separate
the background from the object and crop to size of identified object
    or
manually set x and y boundaries for cropping images to size.\n
For threshold enter a number under 65535.
For manual cropping enter xmin-xmax,ymin-ymax (x = vertical, y  = horizontal).\n
\033[33m>>> \033[00m"""
    elif length == "End":
        text = """\nScript will be terminated!\n
    \033[36m'So long, and thanks for all the fish'   -- Douglas Adams\033[00m\n
"""
    else:
        text = """
    \033[36m'For a moment, nothing happened. 
    Then, after a second or so, nothing continued to happen.'   -- Douglas Adams\033[00m\n
Input not one of the options, please choose either the numbers 1-5 or "quit"\n
\033[33m>>> \033[00m"""

    return text

def random_quote():
    """Quote list from The Hitchhiker's Guide to the Galaxy"""

    quotelist = [
    "I'm 50,000 times more intelligent than you, and even I don't know the answer.",
    "Life! Don't talk to me about life!",
    "Freeze? I'm robot, not a refrigerator.",
    "I have a million ideas. They all point to certain death.",
    "I could calculate your chance of survival, but you won't like it!",
    "Funny, how just when you think life can't possibly get worse it suddenly does.",
    "This'll all ends in tears... I just know it.",
    "I think you ought to know I'm feeling very depressed."]
    random_int = random.randint(0, len(quotelist)-1)

    quote = """\033[36m
    '{}'  -- Marvin \033[00m
    """.format(quotelist[random_int])

    return quote


def load_mask(filename, path=None):
    """Load mask and set values to zero and one."""

    if not path:
        filename = Path.cwd() / Path(filename)

    try:
        image = Image.open(filename)
        mask = np.array(image)
        mask[mask > 0.] = 1
        logger.info("Mask ({}) loaded successfully.".format(filename))
    except FileNotFoundError:
        mask = np.zeros((1, 1))

    return mask


def load_tif(filename, path=None):
    """Import .tif file."""

    logger.debug("Loading file: {}".format(filename))

    if path:
        filename = path / Path(filename)

    image = Image.open(filename)
    im = np.array(image)

    return im


def apply_mask(image, mask):
    """Multiply the mask array and image array."""
    if len(mask.shape) == 3:
        mask = mask[ : , : , 0]
    newimage = image * mask

    logger.debug("Mask applied.")
    # For debugging
    # im = Image.fromarray(newimage)
    # im.show()

    return newimage


def invert_mask(mask):
    """Invert mask to erase calibration object from images."""

    inverted_mask = -1 * mask + 1

    logger.debug("Mask inverted.")

    return inverted_mask.astype('uint16')


def save_image(image, filename, path=None):
    """Save image under export path or generate folder."""

    if not path:
        path = Path.cwd() / "Exported"

    if not path.is_dir():
        path.mkdir()

    filename = path / filename
    Image.fromarray(image).save(filename)

    logger.debug("File saved successfully!")

    return


def find_mean_value(image, threshold):
    """Find mean value of pixels in image above given threshold."""

    image = image.astype(float)
    image[image < threshold] = np.nan
    mean = np.nanmean(image)

    return mean


def exp_func(x, a, b, c):
    """Exponential function for fitting exponential curve."""

    return a * np.exp(b * x) + c

def plain_correction_plot(xs, ys):
    """Plot grey values over image number for mean calibration object values."""
    # Get overall mean and max
    mean_val = np.around(np.mean(ys), 0)
    max_val = np.around(np.max(ys), 0)

    # Plot grey values
    plt.plot(xs, ys, c='orange', label="""Grey value in calibration object
Overall mean: {0}\n Max value: {1}""".format(mean_val, max_val))
    plt.xlabel('Image number in stack')
    plt.ylabel('Mean grey value in calibration object')
    plt.legend()
    plt.show()

    return

def correction_fit_plot(xs, ys, fit, lower_boundary, savepath):
    
    # Get overall mean and max
    mean_val = np.around(np.mean(ys), 0)
    max_val = np.around(np.max(ys), 0)

     # Plot grey values
    fig, ax1 = plt.subplots()
    ax1.plot(xs, ys, c='orange', label="""Grey value in calibration object
Overall mean: {0}\n Max value: {1}""".format(mean_val, max_val))
    ax1.axvline(x=np.where(xs==lower_boundary), color='blue', linestyle='--')
    ax1.set_xlabel('Image number in stack')
    ax1.set_ylabel('Mean grey value in calibration object')

    # Calculate fit and plot
    if re.search(r"(1|lin|Lin)", fit):
        # Perform linear regression
        gradient, intercept, r_value, p_value, std_err = linregress(
            xs[lower_boundary:], ys[lower_boundary:])
        # Add linear regression line to plot
        ax1.plot(xs, intercept + gradient*xs, c='blue',
                 label='Linear fit: {0}x + {1}\nfrom image {2} - {3}'.format(
                 round(gradient, 2), round(intercept, 2), lower_boundary,
                 len(xs)))

        correction_factor = np.array(1 / (gradient * xs + intercept))
        logger.info("Correction factor from linear fit: 1/({0}*x + {1})".format(
                    gradient, intercept))
        print("\nCorrection factor from linear fit: 1/({0}*x + {1})".format(
                    gradient, intercept))

        ax2 = ax1.twinx()
        ax2.plot(xs, ys * correction_factor * np.mean(ys), c='darkred', 
                 label='Corrected grey values for calibration object')

    elif re.search(r"(2|exp|Exp)", fit):
        # Perform exponential curve fit
        params, cov = curve_fit(exp_func, xs[lower_boundary:],
                                ys[lower_boundary:])
        # Add exponential fit to plot
        ax1.plot(xs, exp_func(xs, *params), c='blue',
             label='Exponential fit: {0} e^({1} x) + {2}\nfrom image {3} - {4}'.format(
             round(params[0], 2), round(params[1], 2), round(params[2], 2),
             lower_boundary, len(xs)))

        correction_factor = np.array(1/ exp_func(xs, *params))
        logger.info("Correction factor from exponential fit: 1/({0}* e^({1}x) + {2})".format(
                    params[0], params[1], params[2]))
        print("\nCorrection factor from exponential fit: 1/({0}* e^({1}x) + {2})".format(
                    params[0], params[1], params[2]))

        ax2 = ax1.twinx()
        ax2.plot(xs, ys * correction_factor * np.mean(ys), c='darkred', 
                 label='Corrected grey values for calibration object')
    
    # Show and save plot
    ax2.set_ylim(np.round(np.min(ys * correction_factor * np.mean(ys)), 0) - 10,
                 np.round(np.max(ys * correction_factor * np.mean(ys)), 0) + 10)
    ax2.set_ylabel('Corrected mean grey value in calibration object', c='darkred')
    ax2.tick_params(axis='y', color='darkred', labelcolor='darkred')
    ax1.legend()
    plt.show()
    figname = savepath / "Correction_Fit_Plot.png"
    fig.savefig(figname)

    logger.info("Saved figure as Correction_Fit_Plot.png in main folder.")
    print("Saved figure as Correction_Fit_Plot.png in main folder.")

    return correction_factor

def apply_factor(image, factor, scale=1):
    """
    Multiply the image array with a factor and add the offset
    """
    newimage = np.around(image * factor * scale, 0)
    newimage[newimage > 65535] = 65535
    newimage = np.asarray(newimage, dtype = 'uint16')

    logger.debug("Factor applied.")

    return newimage


def find_max_value(image):
    """Find max value of pixels in image."""

    image = image.astype(float)
    maxval = np.max(image)

    return maxval


def corrected_plot(max_before, max_after, mean_before, mean_after, scale, savepath):
    
     # Plot max grey values
    fig, ax1 = plt.subplots()
    ax1.plot(max_after[:, 0], max_after[:, 1], c='darkorange', 
             label="Max grey value in object after correction")
    ax1.plot(max_before[:, 0], max_before[:, 1], '--', c='orange', 
             label="Max grey value in object before correction")
    ax1.plot(mean_after[:, 0], mean_after[:, 1], c='darkred', 
             label="Mean grey value in object after correction")
    ax1.plot(mean_before[:, 0], mean_before[:, 1], '--', c='red', 
             label="Mean grey value in object before correction")
    fig.suptitle("""Mean and max values before and after applying correction
factor with scaling {}""". format(scale))
    ax1.set_xlabel('Image number in stack')
    ax1.set_ylabel('Mean/ Max grey value in object')
    ax1.legend()
    plt.show()
    figname = savepath / "Corrected_Plot.png"
    fig.savefig(figname)

    logger.info("Saved figure as Corrected_Plot.png in main folder.")
    print("Saved figure as Corrected_Plot.png in main folder.")

    return

def find_object_bounds(image, threshold):

    bounds = np.array([0,image.shape[0],0,image.shape[1]])

    image[image < threshold] = 0
    cols = image.sum(axis = 0)
    rows = image.sum(axis = 1)

    bounds[0] = np.min(np.nonzero(rows))
    bounds[1] = np.max(np.nonzero(rows))
    bounds[2] = np.min(np.nonzero(cols))
    bounds[3] = np.max(np.nonzero(cols))

    return bounds

# Functions for each step of analysis
def separate_object(inputpath, maskpath, outputpaths=None):
    """
    Step 1: Separate calibration object and object from each other using a mask
    (white for rough area of calibration object, black for area of object) and
    save both in different folders for further analysis.
    """
    logger.info("Starting with step 1 - Separating calibration object and object ...")
    print("\nStarting with step 1...")
    print(random_quote())

    # Create output paths and folders if none are given
    if not outputpaths or len(outputpaths) != 2:
        outputpaths = [(Path.cwd() / "01_CalibrationObject"), (Path.cwd() / "01_Object")]
        for path in outputpaths:
            if not path.is_dir():
               path.mkdir()

    # Load mask
    mask = load_mask(maskpath)
    if not mask.any():
        logger.error("Mask not found. No transformation possible!")
        print("Mask not found. No transformation of images possible!")
        return
    # Invert mask
    inverted_mask = invert_mask(mask)

    # Load filenames
    files = [x.name for x in inputpath.glob('*.tif')]
    if not files:
        logger.error("No files detected. No transform performed!")
        print("No original files found in folder {}.".format(inputpath))
        return
    logger.info("{} files detected.".format(len(files)))

    # Apply masks and save files out once only calibration object and once object
    for i, file in enumerate(sorted(files)):
        # Load file
        image = load_tif(file, inputpath)
        # Separate calibration object and only save this
        calib_image = apply_mask(image, mask)
        save_image(calib_image, file, outputpaths[0])
        # Repeat procedure with inverted mask and object
        object_image = apply_mask(image, inverted_mask)
        save_image(object_image, file, outputpaths[1])

    logger.info("Step 1 complete.")
    print("Step 1 complete.")

    return


def shift_correction(inputpath):
    """
    Step 2: Take calibration object images, find shift in grey values and plot graph,
    ask which type of fit should be chosen (exponential or linear), find fit
    and plot again, ask if fit is good or should be redone.
    """
    logger.info("Starting with step 2 - Finding correction through calibration object ...")
    print("\nStarting with step 2...")
    print(random_quote())

    # Load filenames
    files = [x.name for x in inputpath.glob('*.tif')]
    if not files:
        logger.error("No files detected. No shift correction performed!")
        print("No calibration object images found in folder {}.".format(inputpath))
        return
    logger.info("{} files detected.".format(len(files)))

    # Threshold for separating background from object
    threshold = 10000
    # TODO: Threshold could be different for other scans...Ask Liz

    # Calculate mean for each calibration object image with only values from
    # calibration object
    meanarray = np.zeros((len(files), 2))
    for i, file in enumerate(sorted(files)):
        image = load_tif(file, inputpath)
        meanarray[i, :] = [i, find_mean_value(image, threshold)]

    # Plot calibration object grey value data
    plain_correction_plot(meanarray[:, 0], meanarray[:, 1])

    # Get user input for correct fit
    fit = input(get_explanation("Fit"))
    while not re.search(r"(1|lin|Lin|2|exp|Exp)", fit):
        if fit == "quit":
            return fit
        fit = input(get_explanation("Fit Error"))

    # Get input for lower boundary for fit
    lower = input(get_explanation("Lower"))
    while not re.match('[0-9]+', lower) or not int(lower) < len(files):
        if lower == "quit":
            return lower
        lower = input(get_explanation("Fit Error"))

    # Perform fit and plot
    correction_factor = correction_fit_plot(meanarray[:, 0], meanarray[:, 1],
                                            fit, int(lower), Path.cwd())

    logger.info("Step 2 complete.")
    print("\nStep 2 complete.")

    return correction_factor


def apply_correction(inputpath, correction_factor, outputpath=None):
    """
    Step 3: 
    """
    logger.info("Starting with step 3 - Applying correction to object ...")
    print("\nStarting with step 3...")

    # Create output path and folder if none is given
    if not outputpath:
        outputpath = Path.cwd() / "03_CorrectedObject"
        if not outputpath.is_dir():
               outputpath.mkdir()

    # Load filenames
    files = [x.name for x in inputpath.glob('*.tif')]
    if not files:
        logger.error("No files detected. No correction applied!")
        print("No original files found in folder {}.".format(inputpath))
        return
    logger.info("{} files detected.".format(len(files)))

    # Get user input for scaling
    scale = input(get_explanation("Scaling"))
    while not re.match('[0-9]+', scale) or not int(scale) < 65535:
        if scale == "quit":
            return scale
        scale = input(get_explanation("Fit Error"))

    logger.info("Scaling factor used: {}".format(scale))

    print("\nApplying correction...")
    print(random_quote())

    # Arrays for before and after correction plot
    maxarray_before = np.zeros((len(files), 2))
    maxarray_after = np.zeros((len(files), 2))
    threshold = 10000
    meanarray_before = np.zeros((len(files), 2))
    meanarray_after = np.zeros((len(files), 2))
    for i, file in enumerate(sorted(files)):
        # Load image
        image = load_tif(file, inputpath)
        # Get mean value of image before applying correction
        maxarray_before[i, :] = [i, find_max_value(image)]
        meanarray_before[i, :] = [i, find_mean_value(image, threshold)]
        # Apply the correction
        newimage = apply_factor(image, correction_factor[i], int(scale))
        # Get mean value of image after applying correction
        maxarray_after[i, :] = [i, find_max_value(newimage)]
        meanarray_after[i, :] = [i, find_mean_value(newimage, threshold)]
        # Save image
        save_image(newimage, file, outputpath)

    # Perform fit and plot
    corrected_plot(maxarray_before, maxarray_after, meanarray_before, 
                   meanarray_after, scale, Path.cwd())

    logger.info("Step 3 complete.")
    print("Step 3 complete.")

    return None


def clean_images(inputpath, outputpath=None):
    """
    Step 4: TODO: Continue here with coding interactive part of step 4, add logger infos
    """
    logger.info("Starting with step 4 - Cleaning images from other objects ...")
    print("\nStarting with step 4...")

    # Create output path and folder if none are given
    if not outputpath:
        outputpath = (Path.cwd() / "04_ObjectCleaned")
        if not outputpath.is_dir():
               outputpath.mkdir()

    # Ask for mask name
    maskpath = input(get_explanation("Cleaning"))
    while not re.match('\w*.tif', maskpath):
        if maskpath == "quit":
            return maskpath
        maskpath = input(get_explanation("Fit Error"))

    # Load mask
    mask = load_mask(maskpath)
    if not mask.any():
        logger.error("Mask not found. No cleaning possible!")
        print("\nMask not found. No cleaning of images possible!")
        return None

    # Load filenames
    files = [x.name for x in inputpath.glob('*.tif')]
    if not files:
        logger.error("No files detected. No cleaning performed!")
        print("No original files found in folder {}.".format(inputpath))
        return None
    logger.info("{} files detected.".format(len(files)))

    # Ask for image range to clean
    cleaning_range = input(get_explanation("Range"))
    success = False
    while not success == True:
        while not re.match('[0-9]{1,}-[0-9]{1,}', cleaning_range):
            if cleaning_range == "quit":
                return cleaning_range
            cleaning_range = input(get_explanation("Fit Error"))

        lower_range, upper_range = list(map(int,cleaning_range.split('-')))[:2]
        upper_range += 1
        if lower_range < upper_range and lower_range >= 0 and upper_range <= len(files):
            success = True
        else:
            cleaning_range = input(get_explanation("Fit Error"))

    # List all files to clean
    files_tochange = sorted(files)[lower_range:upper_range]
    logger.info("{} files to change.".format(len(files_tochange)))

    print("\nCleaning images...")
    print(random_quote())

    # Copy files from input to output that don't need cleaning, if they are not yet in folder
    files_copied = [x.name for x in outputpath.glob('*.tif')]
    files_to_copy = [x for x in files if x not in files_copied and x not in files_tochange]
    for i, file in enumerate(files_to_copy):
        copy_image = load_tif(file, inputpath)
        save_image(copy_image, file, outputpath)

    # Apply mask and save files out if they need to be cleaned
    for i, file in enumerate(files_tochange):
        image = load_tif(file, inputpath)
        cleaned_image = apply_mask(image, mask)
        save_image(cleaned_image, file, outputpath)

    logger.info("Step 4 complete.")
    print("Step 4 complete.")

    return None


def crop_images(inputpath, outputpath=None):
    """
    Step 5: 
    """
    logger.info("Starting with step 5 - Cropping images ...")
    print("\nStarting with step 5...")

    # Create output path and folder if none are given
    if not outputpath:
        outputpath = (Path.cwd() / "05_ObjectCropped")
        if not outputpath.is_dir():
               outputpath.mkdir()

    # Ask for threshold to separate background and object to use automated cropping
    # Alternatively choose max and min x and y values to crop to wish
    crop_input = input(get_explanation("Cropping"))
    success = False
    xbounds = [0,0]
    ybounds = [0,0]
    while not success == True:
        while not re.match('[0-9]{1,5}', crop_input) or re.match('[0-9]{1,}-[0-9]{1,},[0-9]{1,}-[0-9]{1,}', crop_input):
            if crop_input == "quit":
                return crop_input
            crop_input = input(get_explanation("Fit Error"))

        if re.match('[0-9]{1,5}', crop_input) and int(crop_input) < 65535:
            threshold = int(crop_input)
            success = True
        elif re.match('[0-9]{1,}-[0-9]{1,},[0-9]{1,}-[0-9]{1,}', crop_input):
            xbounds, ybounds = crop_input.split(',')[:2]
            xbounds = list(map(int,xbounds.split('-')))[:2]
            ybounds = list(map(int,ybounds.split('-')))[:2]
            if xbounds[1] > xbounds[0] and ybounds[1] > ybounds[0]:
                success = True
            else: 
                crop_input = input(get_explanation("Fit Error"))
        else: 
            crop_input = input(get_explanation("Fit Error"))

    # Load filenames
    files = [x.name for x in inputpath.glob('*.tif')]
    if not files:
        logger.error("No files detected. No cleaning performed!")
        print("No original files found in folder {}.".format(inputpath))
        return None
    logger.info("{} files detected.".format(len(files)))

    # If threshold is given, load all files and find all cropping parameters
    if threshold:
        print("\nFinding cropping parameters...")
        boundsarray = np.zeros([len(files), 4])
        for i, file in enumerate(files):
            image = load_tif(file, inputpath)
            boundsarray[i,:] = find_object_bounds(image, threshold)

        xbounds[0] = int(np.min(boundsarray[:, 0]) - 5)
        xbounds[1] = int(np.max(boundsarray[:, 1]) + 5) 
        ybounds[0] = int(np.min(boundsarray[:, 2]) - 5)
        ybounds[1] = int(np.max(boundsarray[:, 3]) + 5)

    # Check that all cropping parameters are smaller than the images
    image = load_tif(files[0], inputpath)
    if xbounds[0] < 0:
        xbounds[0] = 0 
    if xbounds[1] > image.shape[0]:
        xbounds[1] = image.shape[0]
    if ybounds[0] < 0:
        ybounds[0] = 0
    if ybounds[1] > image.shape[1]:
        ybounds[1] = image.shape[1]

    print("\nCropping images with parameters {0}-{1}, {2}-{3}".format(
        xbounds[0], xbounds[1], ybounds[0], ybounds[1]))
    print(random_quote())

    # Go through all images to crop and save in output path
    for i, file in enumerate(files):
        image = load_tif(file, inputpath)
        cropped_image = image[xbounds[0]:xbounds[1], ybounds[0]:ybounds[1]]
        save_image(cropped_image, file, outputpath)

    logger.info("Step 5 complete.")
    print("Step 5 complete.")

    return



def main():
    """
    Step by step analysis to prepare CT scan data for curvature/I/venation
    script.
    """

    logger.info("********************************************************************")
    logger.info("Interactive script to prepare CT scans for further analysis started.")

    path_original_images = Path.cwd() / "OriginalScan"
    path_separation_mask = Path.cwd() / "Mask.tif"

    correction_factor = "Empty"

    step = input(get_explanation("Long"))

    while step != "quit":
        if step == "1":
            separate_object(path_original_images, path_separation_mask)
        elif step == "2":
            path_calib = Path.cwd() / "01_CalibrationObject"
            if path_calib.is_dir():
                correction_factor = shift_correction(path_calib)
                if isinstance(correction_factor, str) and correction_factor == "quit":
                    logger.warning("Step 2 terminated early. Correction factor not calculated.")
                    break
            else:
                print("Calibration object images not found in folder '01_CalibrationObject'")
        elif step == "3":
            if isinstance(correction_factor, str):
                print("\nStep 2 has to be executed before step 3 to pass on correction factor.")
            else:
                path_object = Path.cwd() / "01_Object"
                if path_object.is_dir():
                    value = apply_correction(path_object, correction_factor)
                    if isinstance(value, str) and value == "quit":
                        logger.warning("Step 3 terminated early. Correction not applied.")
                        break
                else:
                    print("Object images not found in folder '01_Object'")

        elif step == "4":
            path_object = Path.cwd() / "03_CorrectedObject"
            if path_object.is_dir():
                cleaning = clean_images(path_object)
                if isinstance(cleaning, str) and cleaning == "quit":
                    logger.warning("Step 4 terminated early. Correction not applied.")
                    break
            else:
                print("Images not found in folder '03_CorrectedObject'")
        elif step == "5":
            inputpath = Path.cwd() / "04_ObjectCleaned"
            if inputpath.is_dir():
                cropping = crop_images(inputpath)
                if isinstance(cropping, str) and cropping == "quit":
                    logger.warning("Step 5 terminated early. Images not cropped.")
                    break
            else:
                print("Images not found in folder '04_ObjectCleaned'")
            
        else:
            step = input(get_explanation("Wrong Input Error"))
            continue

        step = input(get_explanation("Short"))

    print(get_explanation("End"))
    logger.info("Script terminated.")
    return


if __name__ == '__main__':
    main()
