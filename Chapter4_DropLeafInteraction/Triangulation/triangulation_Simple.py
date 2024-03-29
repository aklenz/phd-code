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
# --------------------------------------------------------------------------------- #

sync = [-1,-2,-1,0,0,0]

# Function to obtain a Camera Calibration Matrix P (3*4) from a given folder with
# a minimum of 6 files with 3D coordinates in the name and 2D coordinates tracked
# (frame, x, y). Returns the Calibration Matrix P (3*4)
def cameraCalib(foldername, folderpath):
	# Get all files in folder and save paths in files
	files = sorted(list(folderpath.glob('{}_*.csv'.format(foldername.split('.')[-1]))))

	# Matrix A for each camera and perspective (each folder)
	matA = np.zeros((len(files)*2,12))
	# Get 2D and 3D coordinates from all files in files
	for i,file in enumerate(files):
		# Get 3D coordinates from name of file and convert to mm
		X, Y, Z = file.name.split('_')[1][0:3]
		X, Y, Z = float(X)*10, float(Y)*10, float(Z)*10
		# print(X, Y, Z)
		# Load data from file
		data = np.loadtxt(file, delimiter=';', skiprows=1)
		# Take mode (most often occuring value) for x and y coordinate (2D)
		x = stats.mode(data[:,1])
		y = stats.mode(data[:,2])
		# print(x, y)
		# Write lines into matA
		matA[(i*2)] = [X, Y, Z, 1., 0., 0., 0., 0., -X*x, -Y*x, -Z*x, -x]
		matA[(i*2)+1] = [0., 0., 0., 0., X, Y, Z, 1., -X*y, -Y*y, -Z*y, -y]

	# print(matA)

	# Mathematical part to obtain Camera Calibration Matrix P from A
	eigenvalues, eigenvectors = np.linalg.eig(np.dot(matA.transpose(),matA))
	# p is Calibration matrix in vector form
	vecp = eigenvectors[:,eigenvalues.argmin()]
	# P is Calibration Matrix in a 3*4 Matrix
	matP = vecp.reshape((3,4))

	return matP;

# Function to obtain a folder pairing depending on number of calibration videos
# and selected pairFactor. Returns a 3*4 Matrix with 2 calib folder names and 2
# data folder names in each line - further analysis steps for each line of matrix
def pairFolders(calibfolders, datafolders, pairFactor):
	# Reshape the names of the datafolders so Central1 and 2, 
	# Lateral1 and 2 and Tip1 and 2 are in same row
	pairdatafolders = np.array(datafolders)
	pairdatafolders = pairdatafolders.reshape(3,2)

	# Reshape the names of the calibration videos depending on how they should be paired
	paircalibfolders = np.array(calibfolders)

	# Save all calib folder pairs in extra variable to add lines to paired folders
	# So also calib images are triangulated (for checking result is triangulated correctly)
	addcalibfolders = paircalibfolders.reshape(int(len(calibfolders)/2),2)
	addcalibfolders = np.hstack((addcalibfolders,addcalibfolders))

	# 3d1 and 2 for each Central, Lateral and Tip videos
	if len(paircalibfolders) == 2:
		paircalibfolders = [paircalibfolders, paircalibfolders, paircalibfolders]
	# 3d1 and 2 for Central and Lateral, 3d3 and 4 for Tip
	elif len(paircalibfolders) == 4 and pairFactor == 'a':
		paircalibfolders = [paircalibfolders[:2], paircalibfolders[:2], paircalibfolders[2:]]
	# 3d1 and 2 for Central, 3d3 and 4 for Lateral and Tip
	elif len(paircalibfolders) == 4 and pairFactor == 'b':
		paircalibfolders = [paircalibfolders[:2], paircalibfolders[2:], paircalibfolders[2:]]
	# 3d1 and 2 for Central, 3d3 and 4 for Lateral and 3d5 and 6 for Tip
	elif len(paircalibfolders) == 6 and pairFactor == 'a':
		pass
	# 3d1 and 2 for Central, 3d3 and 4 for Tip and 3d5 and 6 for Lateral
	elif len(paircalibfolders) == 6 and pairFactor == 'b':
		paircalibfolders = [paircalibfolders[:2], paircalibfolders[4:], paircalibfolders[2:4]]
	# When no option matches, new pairing needs to be implemented
	else:
		raise IOError('\n\nFolder pairing not possible, choose variable "pairFactor" a or b or code option c!\n\n')
	# Reshape calibration folder array
	paircalibfolders = np.reshape(paircalibfolders, (3,2))

	# Stack both folder arrays so 2 calib videos and 
	# the corresponding 2 impact videos are in same line of matrix
	pairedFolders = np.hstack((paircalibfolders, pairdatafolders))
    # Stack calibration folderpairs and calibration/ datafolder pairs
	pairedFolders = np.vstack((pairedFolders, addcalibfolders))

	return pairedFolders;

# Function to perform the triangulation of a single coordinate with a given x1,x2
# and y1,y2 and the camera matrices matPs. Returns a X, Y, Z coordinate (x_world)
def triangulateCoords(xs,ys, matPs):
	# Generate R and b
	matR = np.zeros((4,3))
	vecb = np.zeros((4,1))

	# Sort values for calib matrices and 2D coords into new lines of triangulation equation
	# (x1, y1, 1) == P1 * x_world and (x2, y2, 1) == P2 * x_world
	# --> matR * x_world = vecb)
	for k in range(2):
		# Calculate lines for R and b and add to matR and vecb
		matR[k*2] = [matPs[k,2,0]*xs[k] - matPs[k,0,0], matPs[k,2,1]*xs[k] - matPs[k,0,1], matPs[k,2,2]*xs[k] - matPs[k,0,2]]
		vecb[k*2] = -(matPs[k,2,3]*xs[k] - matPs[k,0,3])
		matR[k*2+1] = [matPs[k,2,0]*ys[k] - matPs[k,1,0], matPs[k,2,1]*ys[k] - matPs[k,1,1], matPs[k,2,2]*ys[k] - matPs[k,1,2]]
		vecb[k*2+1] = -(matPs[k,2,3]*ys[k] - matPs[k,1,3])

	# Solving lin eq. - Real triangulation
	x_world = np.linalg.solve(np.dot(matR.transpose(),matR), np.dot(matR.transpose(),vecb))
	x_world = [x_world[0,0], x_world[1,0], x_world[2,0]]
	#print("3D coordinate for 0/0/0: in mm \n{}".format(x_world))
	return x_world;

# Function to calculate the reprojection error for one triangulated pair of co-
# ordinates in a single camera perspective
def reprojectionError(xs, ys, x_world, matPs):

	error = np.zeros(2)
	# Calculate reprojection error for each camera
	for i in range(2):
		# 
		X = np.array([x_world[0], x_world[1], x_world[2], 1], dtype = float)
		# First build dot product of matP and x_world
		dotPX = np.dot(matPs[i],X)
		# Calculate reprojection error (coord in picture - reconstructed coord)²
		error[i] = np.sqrt((xs[i] - dotPX[0]/dotPX[2])**2 + (ys[i] - dotPX[1]/dotPX[2])**2)

	return error;

# Function to load in 2 data arrays of x and y coordinates from 2 camera perspectives
# checks that length and frame number are the same, performs triangulation and saves
# triangulated coordinates in new file according to newDataPath
def triangulateFile(file, paths, newDataPath , matPs, videono):
	# Open data from both camera perspectives and save in two arrays
	camera1 = np.loadtxt(paths[0], delimiter=';', skiprows = 1,
	 usecols=(0,1,2), unpack=True).transpose()
	camera2 = np.loadtxt(paths[1], delimiter=';', skiprows = 1,
	 usecols=(0,1,2), unpack=True).transpose()

	# Adjust framenumber according to sync factor (difference between camera 1 and
	# camera 2 for drop impact frame, so coordinates are in true sync for triangulation)
	for frame in range(len(camera2[:,0])):
		camera2[frame, 0] = camera2[frame, 0] + sync[videono]

	# Select only rows that have the same frame number in both camera perspectives
	common_frames = np.intersect1d(camera1[:, 0], camera2[:, 0])
	camera1 = camera1[np.isin(camera1[:, 0], common_frames)]
	camera2 = camera2[np.isin(camera2[:, 0], common_frames)]

	# # Crop arrays so both start with same frame number
	# # Check if they already start with same number
	# if not camera1[0,0] == camera2[0,0]:
	# 	# If data set 2 has lower frame numbers than data set 1
	# 	if camera1[0,0] > camera2[0,0]:
	# 		camera2 = camera2[int(np.where(camera2[:,0] == camera1[0,0])[0]):]
	# 		# Sometimes Drop export makes issues (delete first lines if there is big gap between first few images)
	# 	# If data set 1 has lower frame numbers than data set 2
	# 	else:
	# 		camera1 = camera1[int(np.where(camera1[:,0] == camera2[0,0])[0]):]

	# # Crop array so they have same length
	# # Check if they already have same length
	# if not len(camera1) == len(camera2):
	# 	# If data set 1 is longer than data set 2
	# 	if len(camera1) > len(camera2):
	# 		camera1 = camera1[:len(camera2)]
	# 	# If data set 2 is longer than data set 1
	# 	else: 
	# 		camera2 = camera2[:len(camera1)]

	# Make variable to save in frame no, X, Y, Z and one to save reprojection error
	txyz = np.zeros((len(camera1),4))
	errors = np.zeros((len(camera1),2))
	
	# ----- TRIANGULATION IN EACH FRAME ------------------
	# Access every frame and do triangulation
	for frame in range(len(camera1)):
		# Read out x and y coordinates for current frame
		xs = [camera1[frame,1], camera2[frame,1]]
		ys = [camera1[frame,2], camera2[frame,2]]

		# Triangulation of single coordinates
		x_world = triangulateCoords(xs,ys, matPs)

		# Write triangulated coordinate into txyz
		txyz[frame, :] = camera1[frame,0], x_world[0], x_world[1], x_world[2]

		# Calculate reprojection error 
		# = error in 2D between original coordinates and back projected 3D coordinates
		errors[frame] = reprojectionError(xs,ys,x_world, matPs)

	# File name and path
	newFileName = newDataPath / ("{}.csv".format(file))
	# Stack txyz and errors into one variable
	tosave = np.hstack((txyz, errors))
	# Save file as .csv with header
	np.savetxt(newFileName, tosave, delimiter = ',', header = 'frameno,x,y,z,error1,error2')
	
# Function to find all files in 2 paired folders and compare common files from folders
# to then triangulate all files using the calibMatrices
def triangulateFolder(pairedFolders, rawDataDir, analysedPath, calibMatrices, videono):
	print("Starting Triangulation in folders {0} and {1} ...".format(pairedFolders[2],
	 pairedFolders[3]))
	
	# Generate new folder for saving data, get name from pairedFolder name
	# For .3d4 and .3d6 folders are named 3d4 and 3d6 
	if pairedFolders[3][-1] != '2':
		folderName = pairedFolders[3].split('.')[-1]
	# For .Central2 .Lateral2 .Tip2 and .3d2 folders are named Central, Lateral, Tip and 3d
	else:
		folderName = pairedFolders[3].split('.')[-1][:-1]
	# Assemble path for new folder
	newDataPath = analysedPath / folderName
	# Make new subdirectory, if it doesn't exist yet
	if not newDataPath.exists():
		newDataPath.mkdir()

	# Write calibration matrices into a 3D array
	matPs = np.array([calibMatrices[pairedFolders[0]], calibMatrices[pairedFolders[1]]])
	
	# Access folder paths with data
	currentDataDir = [rawDataDir / pairedFolders[2], rawDataDir / pairedFolders[3]]
			
	# List all file names from folder 1 and get also the base of the filename
	foundFiles1 = [f.name for f in sorted(currentDataDir[0].iterdir())]
	fileNames1 = [f.split('_')[-1].split('.')[0] for f in foundFiles1]
	baseName1 = foundFiles1[0].split('_')[0]
	# List all file names from folder 2 and get also the base of the filename
	foundFiles2 = [f.name for f in sorted(currentDataDir[1].iterdir())]
	fileNames2 = [f.split('_')[-1].split('.')[0] for f in foundFiles2]
	baseName2 = foundFiles2[1].split('_')[0]

	# Save all files that appear in both lists into new list
	commonFiles = [name for name in fileNames1 if name in fileNames2]
	# Check if all files have match (appear in both folders), if not print which files havent got a match
	if len(commonFiles)*2 != len(fileNames1) + len(fileNames2):
		notMatching1 = [name for name in fileNames1 if name not in commonFiles]
		notMatching2 = [name for name in fileNames2 if name not in commonFiles]
		notMatching = notMatching1 + notMatching2
		print("Some markers only appear in one camera: {}\nNo triangulation possible for those files.".format(notMatching))

	#---- TRIANGULATION IN EACH FILE ------------------
	# Go through all files that appear in both folders and do triangluation
	for file in commonFiles:
		# Access the corresponding files from each folder and save paths
		paths = [currentDataDir[0] / ("{0}_{1}.csv".format(baseName1,file)),
		 currentDataDir[1] / ("{0}_{1}.csv".format(baseName2,file))]

		# Triangulates coordinates for each frame and save in new file 
		triangulateFile(file, paths, newDataPath , matPs, videono)

	print('Finished triangulation of folder {}\n'.format(folderName))