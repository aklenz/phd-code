# Typochecker and Open CV camera calibration script for automated analysis of
# the drop impact on leaf surfaces experiments conducted in the Botanical
# Garden 2020 by Anne-Kristin Lenz

# ------------------------- IMPORT ------------------------------------------ #
from pathlib import Path
import glob
import re
import numpy as np
import cv2
# --------------------------------------------------------------------------- #

# ------------------------- VARIABLES --------------------------------------- #
# Path to calibration images that work, if calib images are not working
alternativeCalibPath = Path('') / 'CalibrationImages'
# Used calibration pattern size (internal corners of chessboard)
ROWS = 6
COLS = 9
# Cropping parameter for undistort function
alpha = 0  # 1 for retaining all source image pixels in the corrected frame
		   # (black half moons at edges)
           # 0 for cutting the image, so no black pixels are added
# --------------------------------------------------------------------------- #


# ------------------------- HELP FUNCTIONS ---------------------------------- #
# This function checks if there are typos and other formatting issues in the
# folders to be used for further analysis
def typocheck(cwd):
    ## Typo check for main folder and taking species name for further typo
    ## checks
    # Define reg expression for main folder
    mainreg = re.compile(r'2020_0[7-9]_[0-9][0-9]_([A-Z][a-z]*[A-Z][a-z]*)')
    # Obtain name of main folder - is accessed in main function
    # cwd = Path('').absolute()
    mainName = list(cwd.parent.glob('2020_*'))[0].stem
    # print(mainName)
    # Check for Typos, throw exception if folder doesn't match reg expression
    if mainreg.match(mainName):
        print('No format typos in main folder name detected')
    else:
        raise IOError('Folder did not match reg expression')
    # Take species name and save in variable
    species = mainreg.match(mainName).group(1)
    print('Species: ' + species + '\n')

    ## List all folders in new path, rename species name automatically and
    ## check for wrong endings
    # New path to subfolders of HighSpeedCameras
    rawDataDir = cwd.parent / mainName / 'HighSpeedCameras'
    # Check that path was found
    if not rawDataDir.exists():
        raise IOError(
'\n \nCouldn\'t find folder HighSpeedCameras. Check for typos or existence\n')
    # Reg ex for checking the endings of the folder names
    subregs = re.compile(r'{}\.(Calib|3d|Central|Lateral|Tip)\d'.format(
                         species))
    # Save found directories in list
    foundDirectories = [f for f in rawDataDir.iterdir() if f.is_dir()]
    # Only check for typos, if no folder is missing
    if len(foundDirectories) % 2 == 0 and len(foundDirectories) >= 10:
        # Obtain all folder names in new directory
        for file in rawDataDir.iterdir():
            # Check for typos in spcies name
            if file.stem != species:
                print('Folder has typo in species name: ' + str(file).split(
                      '/')[-1])
                # Automatically removes typos
                file.rename(rawDataDir.absolute() / (species + '.' + str(
                            file.name).split('.')[-1]))
                # Rename file variable so next check just checks ending
                file = rawDataDir.absolute() / (species + '.' + str(
                        file.name).split('.')[-1])
                print('Name corrected to: ' + file.name)
            # Check for typos in ending of folder name
            if not subregs.match(file.name):
                print('Ending has typo:' + file.name)
    else:
        print('One or more folders might be missing!')
    print('Rename folders manually until no typos are found and rerun script.')
    return
# --------------------------------------------------------------------------- #
# This function obtains a calibration matrix using OpenCV
# This script is based on the Camera Calibration Tutorial by Alexander
# Mordvintsev and Abid K. 2013 (https://opencv-python-tutroals.readthedocs.io/
# en/latest/py_tutorials/py_calib3d/py_calibration/py_calibration.html)
# written by George Hale and adjusted by Anne-Kristin Lenz in 2020
def obtainCalibMatrix(path_calibImages, alternativePath, path_saveSuccessful,
                      rows = 6, cols =9, alpha=0):
    # print('Start obtaining calibration matrix...')
    # Set termination criteria for finding corners to high precision in
    # pattern: Terminate after a maximum of 20 iterations or when the accuracy
    # is below 0.001
    term_criteria = (cv2.TERM_CRITERIA_MAX_ITER + cv2.TERM_CRITERIA_EPS ,
                     30, 0.001)

    # Arrays to store object points and image points from all the frames
    objectpoints_array = [] #3D points in the real world space
    imagepoints_array = [] #2D points in the image plane.
    # Count for number of calibration images
    noOfCalibImg = 0
    # Make directory for writing processed data
    if not path_saveSuccessful.exists():
        path_saveSuccessful.mkdir()

    # Make list of coordinates for every point in chessboard grid and fills
    # with zeros
    objectpoints = np.zeros((rows*cols, 3), np.float32)
    # Replace first and second coordinate to create a grid [0,0,0], [1,0,0],
    # [2,0,0]
    objectpoints[:,:2] = np.mgrid[0:rows,0:cols].T.reshape(-1,2)
    # Load all images in current working directory and sort
    images = glob.glob('{}/*.tiff'.format(path_calibImages))
    images.sort(reverse =True)
    # Take 50 images spread over all images in list to find points in
    # chessboard-pattern
    for tiff_name in images[1::int(round(len(images)/50,0))]:
        # Read image
        img = cv2.imread(tiff_name)
        # Convert to greyscale format of openCV - cornerSubPix needs different
        # format
        img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        # Find the chessboard internal corners (only approximate positions)
        # Return value not 0 for successful found corners
        returnVal, corners = cv2.findChessboardCorners(img, (rows,cols))

        # Refine successful found corners, add image and object points to
        # arrays
        if returnVal:
            # Refine the corners within a search window of 23*23 pixels,
            # without leaving any area out (-1,-1) and termination criteria as
            # specified above
            ref_corners = cv2.cornerSubPix(img, corners, (11,11), (-1,-1),
                                           term_criteria)
            # Save the refined corners (in pixels) to imagepoints_array array
            imagepoints_array.append(ref_corners)

            # Draw the chessboard corners on the image
            img = cv2.drawChessboardCorners(img, (rows,cols), ref_corners,
                                            returnVal)

            # New name for saving in folder Successful
            name = path_saveSuccessful / ("Chess_" + str(noOfCalibImg) +
                                          ".tiff")
            # Save image with chessboard drawn to it
            cv2.imwrite(str(name), img)

            # Save the object points to the objectpoints_array array
            objectpoints_array.append(objectpoints)
            # Count add 1
            noOfCalibImg += 1
            # Take max 30 calibration images to reduce calculation time of
            # camera matrix
            if noOfCalibImg == 30:
                break

    print("Number of Calibration Images successfully analysed: {}".format(
          noOfCalibImg))

    # If more than  20 calibration images were found, obtain camera matrix and
    # calculate error, if calib images were to blurry, use external folder to
    # calibrate
    if noOfCalibImg > 20:
        # Use object and image point arrays to calibrate the camera matrix
        returnVal, camMtx, distCoeff, rotVecs, transVecs = cv2.calibrateCamera(
         objectpoints_array, imagepoints_array, img.shape[::-1], None, None)
        # Debugging
        print("Camera Matrix: \n", camMtx)

        # Comment: Don't need ROI when alpha = 0, it automatically crops
        newCamMtx, =cv2.getOptimalNewCameraMatrix(camMtx,distCoeff,
                            img.shape[::-1],alpha)
        print('Finished obtaining calibration matrix.')

    # ------------------------- RE-PROJECTION ERROR ------------------------- #
        print('Start obtaining reprojection error...')
        # Variable for summing up error
        tot_error = 0
        # Calculate error for each calibration image
        for i in range(len(objectpoints_array)):
            # Transform the object points to image points
            imagepoints_array2, _ = cv2.projectPoints(objectpoints_array[i],
             rotVecs[i], transVecs[i], camMtx, distCoeff)
            #print(imagepoints_array[i].shape, imagepoints_array2.shape)
            # Calculate the norm between transformation and corner finding
            # algorithm
            error = cv2.norm(imagepoints_array[i], imagepoints_array2,
                             cv2.NORM_L2)/len(imagepoints_array2)
            error = 0
            #print(error)
            # Sum up the single error terms
            tot_error += error
            meanError = tot_error/len(objectpoints_array)
        # Print the mean error
        print("Mean Error in pixels: {} (Should be under 1)\n".format(
              meanError))
        # The matrix calculated is more accurate the closer it is to 0.
        # print('Finished obtaining reprojection error.\n')
    elif path_calibImages == alternativePath:
        camMtx, newCamMtx, distCoeff, img, = False,False,False,False
        meanError = False
        print('Calib not possible')
    else:
        print('Using BackUp folder for calibration')
        # If I make new calib images, I have to adjust the rows and cols
        # to 6 and 9
        camMtx, newCamMtx, distCoeff, img, meanError = obtainCalibMatrix(
         alternativeCalibPath, alternativeCalibPath, path_saveSuccessful, 5, 8)

    return camMtx, newCamMtx, distCoeff, img, meanError
# --------------------------------------------------------------------------- #
# This function undistorts images from a given directory using OpenCV
# This script is based on the Camera Calibration Tutorial by Alexander
# Mordvintsev and Abid K. 2013 (https://opencv-python-tutroals.readthedocs.io/
# en/latest/py_tutorials/ py_calib3d/py_calibration/py_calibration.html)
# written by George Hale and adjusted by Anne-Kristin Lenz in 2020
def undistortImages(camMtx, newCamMtx, distCoeff, img, path_origImages,
                    path_corImages):
    # Load all images in current working directory
    images = glob.glob('{}/*.tiff'.format(path_origImages))
    #print('Start undistorting images...')
    # check if undistorted images already exist
    if len(images) == len(glob.glob('{}/*.tiff'.format(path_corImages))):
        returnVal = True
        print('Folder already filled with undistorted images.')
        return returnVal
    # Debugging
    # print("New Camera Matrix: \n", newCamMtx)
    # Execute undistorsion for every image in images list
    for image in sorted(images):
        # Read in image
        originalIm = cv2.imread(image)
        # Convert to greyscale format of openCV
        originalIm = cv2.cvtColor(originalIm, cv2.COLOR_BGR2GRAY)
        # Undistort the image using previously generated parameters
        correctedIm = cv2.undistort(originalIm, camMtx, distCoeff, None,
                                    newCamMtx)

        # New name for saving
        name = path_corImages / ('Frame_' + image[-11:])
        # Save the corrected image with original file name in specified folder
        cv2.imwrite(str(name), correctedIm)

    if len(images) == len(glob.glob('{}/*.tiff'.format(path_corImages))):
        returnVal = True
    else:
        returnVal = False
    #print('Finished undistorting images.\n')
    return returnVal
# --------------------------------------------------------------------------- #
# This function uses obtainCalibMatrix() and undistortImages() to automatically
# go through one species folder of the drop impact on leaf measurements and
# obtains calibration matrices for both used cameras and undistorts all images.
# All important parameters are saved in the SummaryUndistortionOpenCV.txt
def automatedOpenCV(cwd):
    mainName = list(cwd.parent.glob('2020_*'))[0].stem

    # Define reg expression for main folder
    mainreg = re.compile(r'2020_0[7-9]_[0-9][0-9]_([A-Z][a-z]*[A-Z][a-z]*)')
    species = mainreg.match(mainName).group(1)
    # Path for all folders that need to be accessed (only read, no write)
    rawDataDir = cwd.parent / mainName / 'HighSpeedCameras'
    # Path for writing processed data
    analysedPath = cwd.parent / mainName / 'OpenCVCorrectedImages'

    # Distinguish into camera 1 and camera 2 folders in read only folder
    folderlist = [[],[]]
    #calibUnevenReg = re.compile(r'[a-zA-Z]*\.(Calib1)')
    #calibEvenReg = re.compile(r'[a-zA-Z]*\.(Calib2)')
    unevenReg = re.compile(r'[a-zA-Z]*\.([a-zA-Z3]*[dlp][13579])')
    evenReg = re.compile(r'[a-zA-Z]*\.([a-zA-Z3]*[dlp][2468])')

    # Find all folders in high speed cameras that are not Calib folders and
    # save in separate lists
    for file in rawDataDir.iterdir():
        if unevenReg.match(file.name):
            folderlist[0].append(file.name)
        elif evenReg.match(file.name):
            folderlist[1].append(file.name)

    # Make directory for writing processed data
    if not analysedPath.exists():
        analysedPath.mkdir()

    # Path to find both Calib folders
    pathCalib = [rawDataDir / (str(species) + '.Calib1'), rawDataDir / (str(
                 species) + '.Calib2')]
    # Path to save Calib Images
    pathSuccessful = [analysedPath / (str(species) + '.Calib1'),
                      analysedPath / (str(species) + '.Calib2')]
    # Loop to run over both camera calibrations - first Calib1 then Calib2
    for i in range(2):
        #if i == 1:
        print('\n Camera {}:'.format(i+1))
        # Check that Calib folder exists and obtain the calibration matrix
        if pathCalib[i].exists() and pathCalib[i].is_dir():
            # Obtain the calibration matrix and error
            calibMat, newCalibMat, distCoeff, img, error = obtainCalibMatrix(
             pathCalib[i], alternativeCalibPath, pathSuccessful[i], ROWS, COLS)
            # Count successfully undistorted folders
            undistFolder = 0
            # Go through all folders from folderlist[0] for uneven and
            # folderlist[1] for even
            for j in folderlist[i]:
                print('Undistorting {} ...'.format(j))
                # create new folder for undistorted images
                joinedPath = analysedPath / j
                if not joinedPath.exists():
                    joinedPath.mkdir()
                # undistort images by rewriting function!
                returnVal = undistortImages(calibMat, newCalibMat, distCoeff,
                                            img, rawDataDir / j , joinedPath)

                if returnVal:
                    undistFolder += 1
                else:
                    print('Undistortion of folder {} not complete'.format(j))
            # Open summary file to save all important parameters from camera
            # calib and undistortion
            summary_file = open('SummaryUndistortionOpenCV.txt', 'a')
            summary_file.write('''\n{0}, Camera{1}, {2}, {3}, {4}, {5}, {6},
 {7}, {8}, {9}, {10}, {11}'''.format(species, i+1, error, calibMat[0,0],
                             calibMat[1,1], calibMat[0,2], calibMat[1,2],
                             newCalibMat[0,0], newCalibMat[1,1],
                             newCalibMat[0,2], newCalibMat[1,2], undistFolder))
            summary_file.close()
        else:
            raise IOError('Couldn\'t find calibration folder.')
# --------------------------------------------------------------------------- #


# ------------------------- MAIN CODE --------------------------------------- #
def main():
    print('Execute main function')
    # Obtain name of main folder
    cwd = Path('').absolute()

    # typo check and check folder structure
    typocheck(cwd)

    # automated camera calibration and undistortion of hierarchical folder
    # structure
    automatedOpenCV(cwd)


if __name__=="__main__":
    main()
