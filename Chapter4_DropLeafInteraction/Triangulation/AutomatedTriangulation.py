# This code takes calibration points from two camera perspectives, obtains
# a calibration matrix and performs a triangulation on tracked markers 
# The code was developed first with George Hale during a MSci project.
# This is the second more automated version and was written by Anne-Kristin
# Lenz in November 2020 for the analysis of drop impact on leaf surfaces 
# experiments conducted in the Botanical Garden 2020.

# ------------------------- IMPORT ------------------------------------------------ #
from pathlib import Path
import glob
import re
import cv2
import numpy as np 
import statistics as stats
#-----
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pdb 	# for debugging 	Stop code at pdb.set_trace()
#----
import triangulation_OpenCV as triCV
import triangulation_Simple as triSim
# --------------------------------------------------------------------------------- #

# ------------------------- VARIABLES --------------------------------------------- #
# var a = []
# pairFactor for pairing Calibration Videos to Drop Impact Videos
# a, when just 2 calib folders,
#  , 4 calib folders: 3d1, 3d2 + Central & Lateral 			  | 3d3 + 3d4 + Tip
#  , 6 calib folders: 3d1, 3d2 + Central | 3d3, 3d4 + Lateral | 3d5, 3d6 + Tip
# b, 4 calib folders: 3d1, 3d2 + Central | 3d3 + 3d4 + Lateral & Tip
#  , 6 calib folders: 3d1, 3d2 + Central | 3d3, 3d4 + Tip 	  | 3d5, 3d6 + Lateral
# c, when none of the options match the pairing: you have to code new pairing! 
pairFactor = 'a'

# Sync factor for Central, Lateral and Tip impact
sync = [1,-1,0,0,0,0]

# --------------------------------------------------------------------------------- #

# ------------------------- HELP FUNCTIONS ---------------------------------------- #
# Function to combine all functions in hierarchical order. 
# Takes calibration files from two camera perspectives, does the camera calibration,
# pairs folders, so calib files and corresponding data files are used for triangulation
# and triangulates all files from all species folders and saves them
def cameraCalibAndTriangulation_Simple():
	# Current path
	cwd = Path('').absolute()
	# Get to Tracking folder with subfolders in it
	mainName = list(cwd.parent.glob('2020_*'))[0].stem
	rawDataDir = cwd.parent / mainName / 'Tracking'

	# Save found directories in folder 'Tracking' into list (just names)
	foundDirs = [f.name for f in sorted(rawDataDir.iterdir()) if f.is_dir()]
	
	# Make two lists with calibration folders and data folders
	# Reg expression for calib folders
	calibreg = re.compile(r'{}.3d[1-6]'.format(mainName.split('_')[-1]))
	# Make list of all folders that match calib reg expression
	calibfolders = [folder for folder in foundDirs if calibreg.match(folder)]
	# Reg expression for data folders
	datareg = re.compile(r'{}.(Central|Lateral|Tip)[12]'.format(mainName.split('_')[-1]))
	# Make list of all folders that match data reg expression	
	datafolders = [folder for folder in foundDirs if datareg.match(folder)]

	# Dictionary for saving all Calibration Matrices with their folder names
	calibMatrices = {}

	# For single camera calibration approach
	# Go through all calibfolders for obtaining calibration matrices
	for folder in calibfolders:
		# Get all files in folder and save paths in files
		folderpath = rawDataDir / folder
		# Obtain the Camera Calibration Matrix for the given folder
		matP = triSim.cameraCalib(folder, folderpath)

		calibMatrices[folder] = matP

	# Pair folders so matching calibration and datafolders names are in same line
	# of the matrix pairedFolders
	pairedFolders = triSim.pairFolders(calibfolders, datafolders, pairFactor)
	
	# Path for writing processd data
	analysedPath = cwd.parent / mainName / '3DCoordinates'
	# Make new directory for writing processed data, if it does not yet exists
	if not analysedPath.exists():
		analysedPath.mkdir()
	
	#--- TRIANGULATION IN EACH FOLDERPAIR ---------------------------------
	# Iterate over the 3 drop impacts and 3D calib videos and do 
	# triangulation for all markers using the pairedFolders
	for i in range(len(pairedFolders)):
		
		triSim.triangulateFolder(pairedFolders[i], rawDataDir, analysedPath, calibMatrices, i)


# ------------------------- MAIN CODE --------------------------------------------- #
def main():
	print('Execute main function')

	# Perform camera calibration and triangulation for current species sample
	cameraCalibAndTriangulation_Simple()


if __name__=="__main__":
	main()