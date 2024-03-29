# This code performs all analysis for 3D coordinates of a drop impacting a leaf
# and the leaf response (tracked marker points on leaf). It was written by 
# Anne-Kristin Lenz in January 2021 for the analysis of drop impact on leaf surfaces 
# experiments conducted in the Botanical Garden 2020.

# ------------------------- IMPORT ------------------------------------------------ #
from pathlib import Path
import csv
import glob
import re
import numpy as np 
import statistics as stats
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.transform import Rotation
from csaps import csaps
from scipy.signal import find_peaks, butter, filtfilt, argrelextrema
from scipy.optimize import minimize
from scipy.optimize import curve_fit
# --------------------------------------------------------------------------------- #

# ------------------------- VARIABLES --------------------------------------------- #
framesPerSecond = 1502 # in 1/s

# --------------------------------------------------------------------------------- #

# ------------------------- HELP FUNCTIONS ---------------------------------------- #
# Find folders (spilt into calib folders and data folders) and paths
def findFoldersAndPaths(speciesFolder):

	rawDataDir = speciesFolder / '3DCoordinates'

	# Save found directories in folder '3DCoordinates' into list (just names)
	foundDirs = [f.name for f in sorted(rawDataDir.iterdir()) if f.is_dir()]
	
	# Make two lists with calibration folders and data folders
	# Reg expression for calib folders
	calibreg = re.compile(r'3d|3d4|3d6')
	# Make list of all folders that match calib reg expression
	calibfolders = [folder for folder in foundDirs if calibreg.match(folder)]
	# Reg expression for data folders
	datareg = re.compile(r'Central|Lateral|Tip')
	# Make list of all folders that match data reg expression	
	datafolders = [folder for folder in foundDirs if datareg.match(folder)]

	# Path for writing processd data
	analysedPath = speciesFolder / 'Analysis'
	# Make new directory for writing processed data, if it does not yet exists
	if not analysedPath.exists():
		analysedPath.mkdir()

	return rawDataDir, datafolders, calibfolders, analysedPath

# Read out the coordinates of the markers for frame 2 to 40 and build mean
def loadMarkerPositions(folderpath):
	# Make filenamelist of files within folder
	files = sorted(list(folderpath.glob('*.csv')))

	# Variables for saving mean position of markers before impact
	meanXs, meanYs, meanZs, marker = np.zeros((len(files),1)), np.zeros((len(files),1)), np.zeros((len(files),1)), []
	# Loop over all files in folder to obtain mean position before impact
	for i,file in enumerate(files):
		if file.stem == "Drop":
			# Load whole drop data
			drop = np.loadtxt(file, delimiter=',', skiprows=1, usecols= (1,2,3))
		else:
			pass
		# Write name in marker
		marker.append(file.stem)
		# Read out lines 2 to 40 and coordinates x,y,z from file
		xs,ys,zs = np.loadtxt(file, delimiter=',', skiprows=1, usecols= (1,2,3),
	 		unpack= True, max_rows = 40)
		# Build mean for each coordinate and save in array for all files
		meanXs[i], meanYs[i], meanZs[i] = np.mean(xs), np.mean(ys), np.mean(zs)

	return marker, meanXs, meanYs, meanZs, drop

# Find best fit line to noisy data points and return unit vector
# https://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
# Accessed on 4th of Feb 2021
def findLinearUnitVec(data):
	# Build mean to normalise points
	meanPoint = data.mean(axis=0)
	# SVD on normalised data
	uu, dd, vv = np.linalg.svd(data-meanPoint)
	# vv[0] is first principal component = unit vector of direction of line with best fit 
	bestFitUnitVector = vv[0] / np.linalg.norm(vv[0])

	return bestFitUnitVector

# Get rotation matrix from old and new coordinate axis
def getRotMatFromMarkers(marker, meanXs, meanYs, meanZs, drop):
	# Coordinate rotation if axis align with leaf surface
	#line2reg = re.compile(r'l2|l1|c0|r1|r2')
	#line2Indices= [marker.index(m) for m in marker if line2reg.match(m)]
	#line2 = np.concatenate(meanXs[line2Indices], meanYs[line2Indices], 
	#				meanZs[line2Indices]), axis=1)

	# Get unit vectors of new coordinate system
	# unit vector of midrib of leaf lamina (lying in direction of length)
	# Will be new x-axis
	#new_x = findLinearUnitVec(line1Xs, line1Ys, line1Zs)
	# unit vector of points across leaf lamina (lying in direction of width)
	#approx_new_y = findLinearUnitVec(line2Xs, line2Ys, line2Zs)

	# New z axis is norm to leaf surface
	#new_z = np.cross(new_x, approx_new_y)/ np.linalg.norm(np.cross(new_x, approx_new_y)) 
	# New y axis is norm to midvein and new z axis 
	# (roughly in direction of width but othorgonal to both x and z)
	#new_y = np.cross(new_x, new_z)/ np.linalg.norm(np.cross(new_x,new_z))
	#print(new_x, new_y, new_z)

	# Coordinate rotation if axis align with gravity and mid vein

	# New z is against drop falling direction = gravity
	new_z = -findLinearUnitVec(drop)
	# Get new x direction from mid vein (line 1)
	# Regular expressions to find coordinates along mid vein
	line1reg = re.compile(r'b2|b1|c0|t1|t2')
	# Find indices matching the marker names from the reg expressions
	line1Indices= [marker.index(m) for m in marker if line1reg.match(m)]
	# Slice coordinates arrays so only coordinates from reg expressions are left
	line1 = np.concatenate((meanXs[line1Indices], meanYs[line1Indices], 
							meanZs[line1Indices]), axis=1)
	# Is not yet perfectly horizontally aligned --> approximate
	approx_new_x = findLinearUnitVec(line1)
	# New y axis is norm to plane going vertically through mid vein
	# Will be across leaf and horizontal
	new_y = np.cross(new_z, approx_new_x)/ np.linalg.norm(np.cross(new_z, 
													approx_new_x))
	# New x is norm to y,z and will be in direction of mid vein, but horizontal
	new_x = np.cross(new_y, new_z)/ np.linalg.norm(np.cross(new_y,new_z))
	
	# Original axes
	axes = np.array([[1,0,0],[0,1,0],[0,0,1]])
	# New axes
	target_axes = np.array([new_x,new_y,new_z])

	# Find rotation matrix between unit vectors and coordinate system
	rotMat = Rotation.from_matrix(np.dot(axes.T, target_axes)).as_matrix()

	return rotMat;

# Get a vector to flip axis to right handed coordinate system
# X from base to leaf tip, Y from middle to left side, Z from surface to sky
def getFlipVecFromRotMarkers(marker, rotPoints):
	# Create vector that does not change coordinates when multiplied linewise
	flipVec = np.array([1,1,1])

	# Check if x-Axsis needs flipping (End state: Leaf-Tip > Leaf-Base)
	if rotPoints[marker.index('t2'),0] < rotPoints[marker.index('b2'),0]:
		flipVec[0] = -1
	# Check if y-Axsis needs flipping (End state: Left-Lamina > Right-Lamina)
	if rotPoints[marker.index('l2'),1] < rotPoints[marker.index('r2'),1]:
		flipVec[1] = -1
	# Check if z-Axsis needs flipping (End state: Drop > Leaf)
	if rotPoints[marker.index('Drop'),2] < rotPoints[marker.index('t2'),2]:
		flipVec[2] = -1

	return flipVec;

# Get all transformation parameters to align leaf and coordinate system
# X along mid vein, Y perpendicular to vein, Z against gravity, b2 at origin
# Optional scaling factor with leaf length = 1
def transformParams(folderpath):
	# Get mean marker positions from before impact (frames 2 to 40)
	marker, meanXs, meanYs, meanZs, drop = loadMarkerPositions(folderpath)

	# Calculate the rotation matrix for reorientation of coordinate system
	rotMat = getRotMatFromMarkers(marker, meanXs, meanYs, meanZs, drop)
	# Rotate points
	rotPoints = np.dot(np.hstack((meanXs, meanYs, meanZs)), rotMat.T)
		
	# Create a vector for flipping axsis
	flipVec = getFlipVecFromRotMarkers(marker, rotPoints)
	# Flip axis with flip vector
	flipPoints = rotPoints * flipVec

	# Get new origin (leaf base) and use as translation vector
	translVec = np.copy(flipPoints[marker.index('b2')])
	# Translate points
	translPoints = flipPoints - translVec

	# Scale coordinates
	# Uniform scaling, so lamina length is 1
	scaleFac = 1/translPoints[marker.index('t2')][0]

	return rotMat, flipVec, translVec, scaleFac

# Figure for controlling bins of FFT
def setupBinFig():
	# Make figure
	binfig = plt.figure(figsize = (10,8))
	binfig.suptitle("FFT and bins")

	binax = []

	# Subplot 0: FFT central impact
	binax.append(binfig.add_subplot(3,1,1))
	binax[0].set_ylabel('FFT in Z')
	binax[0].set_title('Central impact', loc='left')
	binax[0].set_xlim(0, 50)
	plt.grid(True)
	
	# Subplot 1: FFT lateral impact
	binax.append(binfig.add_subplot(312, sharex=binax[0], sharey=binax[0]))
	binax[1].set_ylabel('FFT in Z')
	binax[1].set_title('Lateral impact', loc='left')
	binax[1].set_xlim(0, 50)
	plt.grid(True)

	# Subplot 2: FFT tip impact
	binax.append(binfig.add_subplot(313, sharex=binax[0], sharey=binax[0]))
	binax[2].set_xlabel('Frequency in Hz')
	binax[2].set_ylabel('FFT in Z')
	binax[2].set_title('Tip impact', loc='left')
	binax[2].set_xlim(0, 50)
	plt.grid(True)

	return binfig, binax

# Plot with all marker points after they have been aligned
def overviewFig(folder, files, rotMat, flipVec, translVec, scaleFac, analysedPath, binax):
	# Colours for figure
	colour =  ['#666666', '#ffcf33', '#dba400', '#000000', '#96ed1c', '#578e0b', 
	'#f471c8', '#71094f', '#29a6ff', '#0062a8']

	# Make figure
	fig1 = plt.figure(figsize = (16,8))
	fig1.suptitle("Drop impact: {}".format(folder))

	# Subplot 1: 3D marker track
	ax1 = fig1.add_subplot(3,4,2, projection = '3d')
	ax1.plot([0,1/scaleFac], [0,0], [0,0], '-', color='grey', alpha=0.5)
	font = {'color':  'grey', 'weight': 'normal', 'size': 12,}
	ax1.text(0, 0, 0, 'base', fontdict=font)
	ax1.text(1/scaleFac, 0, 0, 'tip', fontdict=font)
	ax1.set_xlim(-0.05/scaleFac, 1.05/scaleFac)
	ax1.set_ylim(-0.55/scaleFac, 0.55/scaleFac)
	ax1.set_zlim(-0.55/scaleFac, 0.55/scaleFac)
	plt.grid(False)

	# Subplot 4: z trajectory over time
	ax4 = fig1.add_subplot(323)
	ax4.set_xlabel('Time')
	ax4.set_ylabel('Z position in mm')
	ax4.set_title('Up-down oscillation')
	#ax4.set_ylim(-0.6/scaleFac, 0.6/scaleFac)
	plt.grid(True)

	# Subplot 2: x trajectory over time
	ax2 = fig1.add_subplot(348, sharey=ax4)
	ax2.set_xlabel('Time')
	ax2.set_ylabel('X position in mm')
	ax2.set_title('Back-forth oscillation')
	#ax2.set_ylim(-0.1/scaleFac, 1.1/scaleFac)
	plt.grid(True)

	# Subplot 3: y trajectory over time
	ax3 = fig1.add_subplot(347, sharey=ax4)
	ax3.set_xlabel('Time')
	ax3.set_ylabel('Y position in mm')
	ax3.set_title('Sideway oscillation')
	#ax3.set_ylim(-0.6/scaleFac, 0.6/scaleFac)
	plt.grid(True)

	# Subplot 5: error cam 1 over time
	ax5 = fig1.add_subplot(343)
	ax5.set_xlabel('Time')
	ax5.set_ylabel('Reprojection error in cam 1')
	plt.grid(True)

	# Subplot 6: error cam 2 over time
	ax6 = fig1.add_subplot(344)
	ax6.set_xlabel('Time')
	ax6.set_ylabel('Reprojection error in cam 2')
	plt.grid(True)

	# Subplot 7: Position of markers in top down view
	ax7 = fig1.add_subplot(341)
	ax7.set_xlabel('Marker position x')
	ax7.set_ylabel('Marker position y')
	ax7.plot([0,1/scaleFac], [0,0], '-', color='grey', alpha=0.5)
	ax7.set_ylim(-0.55/scaleFac, 0.55/scaleFac)
	ax7.set_xlim(-0.05/scaleFac, 1.05/scaleFac)

	# Subplot 8: FFT in z
	ax8 = fig1.add_subplot(325)
	ax8.set_xlabel('Frequency in Hz')
	ax8.set_ylabel('FFT in Z')
	ax8.set_xlim(0, 150)
	plt.grid(True)

	# Subplot 9: FFT in x
	ax9 = fig1.add_subplot(3,4,12, sharey=ax8)
	ax9.set_xlabel('Frequency in Hz')
	ax9.set_ylabel('FFT in X')
	ax9.set_xlim(0, 150)
	plt.grid(True)

	# Subplot 10: FFT in y
	ax10 = fig1.add_subplot(3,4,11, sharey=ax8)
	ax10.set_xlabel('Frequency in Hz')
	ax10.set_ylabel('FFT in Y')
	ax10.set_xlim(0, 150)
	plt.grid(True)

	# Define subplot for FFT bin plot for each folder
	if folder == 'Central':
		j = 0
	elif folder == 'Lateral':
		j = 1
	else:
		j = 2

	# Add each marker data to plot
	for i, file in enumerate(files):
		# Ignore drop coordinate
		if not file.stem == 'Drop':				
			# Read out data
			time = np.loadtxt(file, delimiter=',', skiprows=1, usecols= 0,
 				unpack= False)
			data = np.loadtxt(file, delimiter=',', skiprows=1, usecols= (1,2,3),
 				unpack= False)
			errors = np.loadtxt(file, delimiter=',', skiprows=1, usecols= (4,5),
 				unpack= False)

			# Check if time data has holes
			if int(time[-1]-time[0]) > len(time)-1:
				print('Warning: holes in data for marker: {}'.format(file))

			# Transform time
			time = (time-time[0]) / framesPerSecond
			# Transform coordinates
			data = (flipVec* np.dot(data, rotMat.T) - translVec) #*scaleFac
			# Build mean and std
			means = np.mean(data[0:40,], axis=0)
			sds = np.std(data, axis=0)
			# For oscillations set all data relative to start position
			oscilldata = data - means

			# Time start value for all filter actions in analysis
			if file.stem == 'c0':
				start = np.where(oscilldata[:80,2]<-0.1)[-1][0]

			# FFT data in x, y and z
			x_fft = np.fft.fft(data[:, 0])
			y_fft = np.fft.fft(data[:, 1])
			z_fft = np.fft.fft(data[:, 2])
			# Frequencies for fft
			freq = np.fft.fftfreq(len(time), d=1/framesPerSecond)
			# Negative part will be cut off (second half of fft data)
			# as well as taking the abs of the fft
			# https://pythonnumericalmethods.berkeley.edu/notebooks/chapter24.04-FFT-in-Python.html
			n_oneside = len(time)//2
			

			ax1.plot(data[:, 0], data[:, 1], data[:, 2], '-', color=colour[i])
			ax2.plot(time, oscilldata[:,0], '-', color=colour[i])
			ax3.plot(time, oscilldata[:,1], '-', color=colour[i])
			ax4.plot(time, oscilldata[:,2], '-', color=colour[i])
			ax5.plot(time, errors[:,0], '-', color=colour[i])
			ax6.plot(time, errors[:,1], '-', color=colour[i])
			ax7.plot(means[0], means[1],'X', color=colour[i])
			font = {'color':  colour[i], 'weight': 'normal', 'size': 12,}
			ax7.text(means[0], means[1], file.stem, fontdict=font)
			ax8.plot(freq[1:n_oneside], np.abs(z_fft[1:n_oneside]), '-', color=colour[i])
			binax[j].plot(freq[1:n_oneside], np.abs(z_fft[1:n_oneside]), '-', color=colour[i])
			ax9.plot(freq[1:n_oneside], np.abs(x_fft[1:n_oneside]), '-', color=colour[i])
			ax10.plot(freq[1:n_oneside], np.abs(y_fft[1:n_oneside]), '-', color=colour[i])
		# Add drop to 3D plot
		else:
			drop = np.loadtxt(file, delimiter=',', skiprows=1, usecols= (1,2,3),
 				unpack= False)
			drop = (flipVec* np.dot(rotMat, drop.T).T - translVec) #*scaleFac
			ax1.plot(drop[:, 0], drop[:, 1], drop[:, 2], '-', color=colour[i])
	
	# Adjust layout and save plot
	plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.90, wspace=0.25, hspace=0.35)
	fig1name = analysedPath / "Overview_{}.pdf".format(folder)
	fig1.savefig(fig1name)
	plt.close(fig1)

	return start

# Cutoff reprojection error is defined in here!
def findValidMarker(files):
	# Colours for figure
	colour =  ['#666666', '#ffcf33', '#dba400', '#000000', '#96ed1c', '#578e0b', 
	'#f471c8', '#71094f', '#29a6ff', '#0062a8']

	validMarkers = []
	validMarkerNames = []
	validColours = []
	# Find only valid data
	for i,file in enumerate(files):
		# Load data
		errors = np.loadtxt(file, delimiter=',', skiprows=1, usecols= (4,5),
				unpack= False)

		# If mean error are higher than 10 pixel do not analyse image
		meanerror = np.mean(errors, axis = 0)
		if np.any(meanerror > 10):
			continue

		# Add to list of valid data
		validMarkers.append(file)
		validMarkerNames.append(file.stem)
		validColours.append(colour[i])

	return validMarkers, validMarkerNames, validColours

def findFftBins(validMarkers, validMarkerNames, flipVec, rotMat, translVec):

	# Find bins for fft analysis for all data
	#fftBinMarkerOrder = ['r2', 'l2', 'l1', 'r1', 't2', 't1', 'c0', 'b1', 'b2', 'Drop']
	#commonMarkers = [name for name in fftBinMarkerOrder if name in validMarkerNames]

	# Initialise arrays for saving frequencies with FFT peaks and minima between the first 3
	peak1 = np.full(10, np.nan)
	peak2 = np.full(10, np.nan)
	peak3 = np.full(10, np.nan)
	peak4 = np.full(10, np.nan)
	peak5 = np.full(10, np.nan)
	bin1 = np.full(10, np.nan)
	bin2 = np.full(10, np.nan)

	for i, marker in enumerate(validMarkers):
		# Ignore drop marker
		if marker.stem == 'Drop':
			continue
		# Load data
		data = np.loadtxt(marker, delimiter=',', skiprows=1, usecols= (1,2,3), unpack= False)
		# Transform coordinates
		data = (flipVec* np.dot(data, rotMat.T) - translVec)
		means = np.mean(data[0:40,2], axis=0)

		# Only take Z oscillation data
		zdata = data[:, 2] - means
		# Filter data with cubic spline
		x = np.linspace(0, len(zdata)-1, len(zdata))
		zsmoothed = csaps(x, zdata, x, smooth=0.1)

		# FFT data in z
		z_fft = np.fft.fft(zsmoothed)
		# Frequencies for fft
		freq = np.fft.fftfreq(len(zsmoothed), d=1/framesPerSecond)
		# Find peaks in fft
		zpeaks_fft, _ = find_peaks(np.abs(z_fft), prominence=3, height=25)
		
		# Find only relevant peaks for bins, not negative half
		zpeaks_fft = zpeaks_fft[:(len(zpeaks_fft)//2)]
		peakfreq = freq[zpeaks_fft]
		# Also dont take any main frequencies below 3 Hz
		peakfreq = peakfreq[peakfreq>=3]
		zpeakfreq = np.abs(z_fft[zpeaks_fft[-len(peakfreq):]])
		zpeaks_fft = zpeaks_fft[-len(peakfreq):]

		acrossmarker = ['r2', 'l2', 'l1', 'r1']
		# Add peak frequencies and bin frequencies to arrays
		# for some only across markers
		if zpeaks_fft.any():
			peak1[i] = freq[zpeaks_fft[0]]
			indexMax1 = zpeaks_fft[0]
			#indexMax1 = zpeaks_fft[np.where(zpeakfreq == np.max(zpeakfreq))][0]
			if len(zpeaks_fft) == 1:
				if marker.stem in acrossmarker:
					bin1[i] = freq[indexMax1 + argrelextrema(np.abs(z_fft[indexMax1:]), np.less)[0][0]]
		if len(zpeaks_fft) > 1:
			indexMax2 = zpeaks_fft[1]
			if marker.stem in acrossmarker:
				peak2[i] = freq[zpeaks_fft[1]]
				#indexMax2 = zpeaks_fft[np.where(zpeakfreq == np.sort(zpeakfreq)[-2])][0]
				bin1[i] = freq[indexMax1 + np.argmin(np.abs(z_fft[indexMax1:indexMax2]))]
			if len(zpeaks_fft) == 2:
				bin2[i] = freq[indexMax2 + argrelextrema(np.abs(z_fft[indexMax2:]), np.less)[0][0]]
		if len(zpeaks_fft) > 2:
			peak3[i] = freq[zpeaks_fft[2]]
			indexMax3 = zpeaks_fft[2]
			#indexMax3 = zpeaks_fft[np.where(zpeakfreq == np.sort(zpeakfreq)[-3])][0]
			#bin3 = (indexMax2 + indexMax3) // 2
			bin2[i] = freq[indexMax2 + np.argmin(np.abs(z_fft[indexMax2:indexMax3]))]
		if len(zpeaks_fft) > 3:
			peak4[i] = freq[zpeaks_fft[3]]
		if len(zpeaks_fft) > 4:
			peak5[i] = freq[zpeaks_fft[3]]

	return bin1, bin2, peak1, peak2, peak3, peak4, peak5

def find_individual_peaks(peaks, threshold=2):
	peaks = np.sort(peaks[~np.isnan(peaks)].flatten())
	individualpeaks = []
	start_idx = 0

	for i in range(1, len(peaks)):
		if abs(peaks[i] - peaks[i-1]) > threshold:
			individualpeaks.append(np.nanmedian(peaks[start_idx:i]))
			start_idx = i

	individualpeaks.append(np.nanmedian(peaks[start_idx:]))  # Append the last subarray

	return np.array(individualpeaks)

def setupFigAnalysis(folder, scaleFac):

	# Make figure
	fig = plt.figure(figsize = (16,23)) # was 16,23
	fig.suptitle("Drop impact: {}".format(folder))

	axs = []
	# Subplot 0: Position of markers in top down view
	axs.append(fig.add_subplot(741))
	axs[0].set_xlabel('Marker position x')
	axs[0].set_ylabel('Marker position y')
	axs[0].plot([0,1/scaleFac], [0,0], '-', color='grey', alpha=0.5)
	axs[0].set_ylim(-0.55/scaleFac, 0.55/scaleFac)
	axs[0].set_xlim(-0.05/scaleFac, 1.05/scaleFac)
	# Subplot 1: 3D marker track
	axs.append(fig.add_subplot(7,4,2, projection = '3d'))
	axs[1].plot([0,1/scaleFac], [0,0], [0,0], '-', color='grey', alpha=0.5)
	font = {'color':  'grey', 'weight': 'normal', 'size': 12,}
	axs[1].text(0, 0, 0, 'base', fontdict=font)
	axs[1].text(1/scaleFac, 0, 0, 'tip', fontdict=font)
	axs[1].set_xlim(-0.05/scaleFac, 1.05/scaleFac)
	axs[1].set_ylim(-0.55/scaleFac, 0.55/scaleFac)
	axs[1].set_zlim(-0.55/scaleFac, 0.55/scaleFac)
	plt.grid(False)
	# Subplot 2: error cam 1 over time
	axs.append(fig.add_subplot(743))
	axs[2].set_xlabel('Time in s')
	axs[2].set_ylabel('Reprojection error in cam 1')
	plt.grid(True)
	# Subplot 3: error cam 2 over time
	axs.append(fig.add_subplot(744))
	axs[3].set_xlabel('Time in s')
	axs[3].set_ylabel('Reprojection error in cam 2')
	plt.grid(True)

	# Subplot 4: z trajectory over time
	axs.append(fig.add_subplot(734))
	axs[4].set_xlabel('Time in s')
	axs[4].set_ylabel('Z position in mm')
	axs[4].set_title('Up-down oscillation', loc='left')
	#axs[4].set_xlim(0, 2)
	#axs[4].set_ylim(-17, 10.5)
	#ax4.set_ylim(-0.6/scaleFac, 0.6/scaleFac)
	plt.grid(True)
	# Subplot 5: y trajectory over time
	axs.append(fig.add_subplot(735, sharey=axs[4]))
	axs[5].set_xlabel('Time in s')
	axs[5].set_ylabel('Y position in mm')
	axs[5].set_title('Sideway oscillation')
	#ax3.set_ylim(-0.6/scaleFac, 0.6/scaleFac)
	plt.grid(True)
	# Subplot 2: x trajectory over time
	axs.append(fig.add_subplot(736, sharey=axs[4]))
	axs[6].set_xlabel('Time in s')
	axs[6].set_ylabel('X position in mm')
	axs[6].set_title('Back-forth oscillation', loc='left')
	#ax2.set_ylim(-0.1/scaleFac, 1.1/scaleFac)
	plt.grid(True)

	# Subplot 7: FFT in z
	axs.append(fig.add_subplot(737))
	axs[7].set_xlabel('Frequency in Hz')
	axs[7].set_ylabel('FFT in Z')
	axs[7].set_xlim(0, 150)
	plt.grid(True)
	# Subplot 8: FFT in y
	axs.append(fig.add_subplot(738, sharey=axs[7]))
	axs[8].set_xlabel('Frequency in Hz')
	axs[8].set_ylabel('FFT in Y')
	axs[8].set_xlim(0, 150)
	plt.grid(True)
	# Subplot 9: FFT in x
	axs.append(fig.add_subplot(739, sharey=axs[7]))
	axs[9].set_xlabel('Frequency in Hz')
	axs[9].set_ylabel('FFT in X')
	axs[9].set_xlim(0, 150)
	plt.grid(True)

	# Filtering just permanent displacement (under 1 Hz)
	# Subplot 10: z trajectory over time
	axs.append(fig.add_subplot(7,3,10))
	axs[10].set_xlabel('Time in s')
	axs[10].set_ylabel('Z position in mm')
	axs[10].set_title('Permanent displacement and bending', loc='left')
	#axs[10].set_xlim(0, 2)
	#axs[10].set_ylim(-17, 10.5)
	plt.grid(True)
	# Subplot 11: y trajectory over time
	axs.append(fig.add_subplot(7,3,11, sharey=axs[10]))
	axs[11].set_xlabel('Time in s')
	axs[11].set_ylabel('Y position in mm')
	plt.grid(True)
	# Subplot 12: x trajectory over time
	axs.append(fig.add_subplot(7,3,12, sharey=axs[10]))
	axs[12].set_xlabel('Time in s')
	axs[12].set_ylabel('X position in mm')
	plt.grid(True)

	# Filtering Bending peak
	# Subplot 13: z trajectory over time
	axs.append(fig.add_subplot(7,3,13))
	axs[13].set_xlabel('Time in s')
	axs[13].set_ylabel('Z position in mm')
	axs[13].set_title('Pure Twisting', loc='left')
	#axs[13].set_xlim(0, 0.5)
	#axs[13].set_ylim(-4.7, 4.7)
	plt.grid(True)
	# Subplot 14: y trajectory over time
	axs.append(fig.add_subplot(7,3,14, sharey=axs[13]))
	axs[14].set_xlabel('Time in s')
	axs[14].set_ylabel('Y position in mm')
	plt.grid(True)
	# Subplot 15: x trajectory over time
	axs.append(fig.add_subplot(7,3,15, sharey=axs[13]))
	axs[15].set_xlabel('Time in s')
	axs[15].set_ylabel('X position in mm')
	plt.grid(True)

	# Filtering Twisting peak
	# Subplot 16: z trajectory over time
	axs.append(fig.add_subplot(7,3,16))
	axs[16].set_xlabel('Time in s')
	axs[16].set_ylabel('Z position in mm')
	axs[16].set_title('Higher frequencies along midrib', loc='left')
	#axs[16].set_xlim(0.025, 0.225)
	#axs[16].set_ylim(-1.7, 1.5)
	plt.grid(True)
	# Subplot 17: y trajectory over time
	axs.append(fig.add_subplot(7,3,17, sharey=axs[16]))
	axs[17].set_xlabel('Time in s')
	axs[17].set_ylabel('Y position in mm')
	plt.grid(True)
	# Subplot 18: x trajectory over time
	axs.append(fig.add_subplot(7,3,18, sharey=axs[16]))
	axs[18].set_xlabel('Time in s')
	axs[18].set_ylabel('X position in mm')
	plt.grid(True)

	# Filtering Higher frequencies
	# Subplot 19: z trajectory over time
	axs.append(fig.add_subplot(7,3,19))
	axs[19].set_xlabel('Time in s')
	axs[19].set_ylabel('Z position in mm')
	axs[19].set_title('Higher Frequencies across midrib', loc='left')
	plt.grid(True)
	# Subplot 20: y trajectory over time
	axs.append(fig.add_subplot(7,3,20, sharey=axs[19]))
	axs[20].set_xlabel('Time in s')
	axs[20].set_ylabel('Y position in mm')
	plt.grid(True)
	# Subplot 21: x trajectory over time
	axs.append(fig.add_subplot(7,3,21, sharey=axs[19]))
	axs[21].set_xlabel('Time in s')
	axs[21].set_ylabel('X position in mm')
	plt.grid(True)

	return fig, axs

def setupDampingFit(folder, scaleFac):
	# Make figure
	fig2 = plt.figure(figsize = (8,8))
	fig2.suptitle("Drop impact: {}".format(folder))

	axs2 = []

	# Filtering Bending peak
	# Subplot 0: z trajectory over time
	axs2.append(fig2.add_subplot(2,1,1))
	axs2[0].set_xlabel('Time in s')
	axs2[0].set_ylabel('Z position in mm')
	axs2[0].set_title('Pure Bending, First peak of FFT in Z')
	plt.grid(True)
	
	# Filtering Twisting peak
	# Subplot 1: z trajectory over time
	axs2.append(fig2.add_subplot(2,1,2))
	axs2[1].set_xlabel('Time in s')
	axs2[1].set_ylabel('Z position in mm')
	axs2[1].set_title('Pure Twisting, Second peak of FFT in Z')
	plt.grid(True)

	return fig2, axs2

def loadAndTransformData(file, rotMat, flipVec, translVec):

	# Read out data
	time = np.loadtxt(file, delimiter=',', skiprows=1, usecols= 0,
			unpack= False)
	data = np.loadtxt(file, delimiter=',', skiprows=1, usecols= (1,2,3),
			unpack= False)
	errors = np.loadtxt(file, delimiter=',', skiprows=1, usecols= (4,5),
			unpack= False)

	# Transform time
	time = (time-time[0]) / framesPerSecond
	# Transform coordinates
	data = (flipVec* np.dot(data, rotMat.T) - translVec) #*scaleFac
	# Build mean and std
	means = np.mean(data[0:40,], axis=0)
	sds = np.std(data, axis=0)
	# For oscillations set all data relative to start position
	oscilldata = data - means

	return time, oscilldata, means, errors

def plotFirstRow(axs, colour, file, time, data, means, errors):

	axs[0].plot(means[0], means[1],'X', color=colour)
	font = {'color':  colour, 'weight': 'normal', 'size': 12,}
	axs[0].text(means[0], means[1], file.stem, fontdict=font)
	data3D = data + means
	axs[1].plot(data3D[:, 0], data3D[:, 1], data3D[:, 2], '-', color=colour)
	axs[2].plot(time, errors[:,0], '-', color=colour)
	axs[3].plot(time, errors[:,1], '-', color=colour)

	return

# Define the damped oscillator function
def damped_oscillator(time, maxA, zeta, omega, displ):	#, phi, C
    # https://en.wikipedia.org/wiki/Damping
    # https://iopscience.iop.org/article/10.1088/1748-3190/ab68a8/meta#bbab68a8as2
    return maxA * np.exp(-zeta * time) * np.cos(omega * time) + displ# - phi)

# Define the cost function (to be minimized)
def cost_function(params, time, z_amp):
    maxA, zeta, omega, displ = params # , phi
    z_pred = damped_oscillator(time, maxA, zeta, omega, displ)#, phi)
    return np.sum((z_amp - z_pred)**2)

def oscillation_fit(time, bending, init_params):
	# Minimize the cost function
	result = minimize(cost_function, init_params, args=(time, bending), method='Nelder-Mead')
	# Extract the fitted parameters: maxA_fit, zeta_fit, omega_fit, phi_fit
	params = result.x
	#params, covariance = curve_fit(damped_oscillator, time, bending, p0=init_params)

	# Generate fitted curve
	#z_pred = damped_oscillator(time, params[0], params[1], params[2], params[3], params[4])
	z_pred = damped_oscillator(time, params[0], params[1], params[2], params[3])

	return z_pred, params

# Plot with all marker points after they have been aligned
def validDataAnalysis(folder, files, rotMat, flipVec, translVec, scaleFac, analysedPath, start, binfreq, peaks):
	# Find data to analyse
	validMarkers, validMarkerNames, colour = findValidMarker(files)

	# Is there any data to analyse
	if not validMarkers:
		print('No data exported for folder: {}.'.format(folder))
		c0_twisting = np.array([False])
		displacement_export = np.full(4,np.nan)
		amplitude_export = np.full(4,np.nan)
		damping_export =  np.full(4,np.nan)
		return displacement_export, c0_twisting.any(), amplitude_export, damping_export

	# Set up figure for all analysis plots
	fig, axs = setupFigAnalysis(folder, scaleFac)

	# Set up second figure for damping coeff fit
	fig2, axs2 = setupDampingFit(folder, scaleFac)

	# Plot bins as vertical lines
	for i,bfreq in enumerate(binfreq):
		if i < 1:
			col = 'navy'
		else:
			col = 'deeppink'
		axs[7].axvline(bfreq, color = col, alpha = 0.5)
		axs[8].axvline(bfreq, color = col, alpha = 0.5)
		axs[9].axvline(bfreq, color = col, alpha = 0.5)

	# Plot peaks in graph
	for peak in peaks:
		axs[7].plot(peak, -10, '2', color= 'black', markersize=10)
		axs[8].plot(peak, -10, '2', color= 'black', markersize=10)
		axs[9].plot(peak, -10, '2', color= 'black', markersize=10)
		font = {'color':  'black', 'weight': 'normal', 'size': 4,}
		axs[7].text(peak, -120, np.around(peak,1), ha='center', fontdict=font)

	# Initialise variables to save information temporarily or permanently
	c0_twisting_present = False
	displacement_export = np.full(4,np.nan)
	amplitude_export = np.full(4,np.nan)
	damping_export =  np.full(4,np.nan)
	
	# Add each marker data to plot
	for i, file in enumerate(validMarkers):

		# Add drop to 3D plot
		if file.stem == 'Drop':
			drop = np.loadtxt(file, delimiter=',', skiprows=1, usecols= (1,2,3),
 				unpack= False)
			drop = (flipVec* np.dot(rotMat, drop.T).T - translVec) #*scaleFac
			axs[1].plot(drop[:, 0], drop[:, 1], drop[:, 2], '-', color=colour[i])
			continue

		# Load data and transform coordinates for analysis
		time, data, means, errors = loadAndTransformData(file, rotMat, flipVec, translVec)
		
		if len(data) < 200:
			continue
		# Plot first row: marker position and errors
		plotFirstRow(axs, colour[i], file, time, data, means, errors)

		if file.stem == 'c0':
			c0_twisting = np.zeros((len(time)-start,3))

		# All other rows plot column by column, first z, then y, then x
		for column in range(3):
			# Choose data
			oscilldata = data[:, 2-column] 

			# Filter data with cubic spline # Try between 0.01 and 0.5
			x = np.linspace(0, len(oscilldata)-1, len(oscilldata))
			smoothedOscill = csaps(x, oscilldata, x, smooth=0.1)

			# FFT data
			fftData = np.fft.fft(smoothedOscill)
			# Frequencies for fft
			freq = np.fft.fftfreq(len(smoothedOscill), d=1/framesPerSecond)

			# Butterworth filter
			# https://stackoverflow.com/questions/12093594/how-to-implement-band-pass-butterworth-filter-with-scipy-signal-butter
			# https://scipy-cookbook.readthedocs.io/items/ButterworthBandpass.html

			# First bin is for permanent displacement
			# Last bin for higher fs without detected peaks
			if binfreq.any():
				b,a = butter(3, binfreq[0]/(framesPerSecond/2), btype = 'low', 
							analog=False)
				bending = filtfilt(b, a, oscilldata[start:])

				b, a = butter(3, binfreq[0]/(framesPerSecond/2), btype='high',
							analog=False)
				higherbending = filtfilt(b, a, oscilldata[start:])

			# Second bin is for bending frequencies
			if len(binfreq) > 1:
				#b, a = butter(3, [binfreq[0]/(framesPerSecond/2), 
				#				  binfreq[1]/(framesPerSecond/2)], 
				#				  btype='band',analog=False)
				#bending = filtfilt(b, a, oscilldata[start:])
			# Third bin is for twisting
			#if len(binfreq) > 3:
				b, a = butter(3, [binfreq[1]/(framesPerSecond/2),
								  binfreq[2]/(framesPerSecond/2)],
								  btype='band',analog=False)
				twisting = filtfilt(b, a, oscilldata[start:])

				b, a = butter(3, binfreq[2]/(framesPerSecond/2), btype='high',
							analog=False)
				highertwisting = filtfilt(b, a, oscilldata[start:])
				if file.stem == 'c0':
					c0_twisting[:, column] = twisting
					c0_twisting_present = True

			# Plot second row: oscillation data
			axs[4+column].plot(time, oscilldata, '-', color=colour[i], markersize=1)
			#axs[4+column].plot(time, smoothedOscill, '-', color=colour[i])

			# Plot third row: FFT data
			# Negative part will be cut off (second half of fft data)
			# as well as taking the abs of the fft
			# https://pythonnumericalmethods.berkeley.edu/notebooks/chapter24.04-FFT-in-Python.html
			n_oneside = len(oscilldata)//2
			axs[7+column].plot(freq[1:n_oneside], np.abs(fftData[1:n_oneside]),
							  '-', color=colour[i])

			# Plot filtered data
			if binfreq.any():
				# Plot forth row: Permanent displacement and bending
				lastcycle = int(framesPerSecond/peaks[0])
				axs[10+column].plot(time[-lastcycle:], bending[-lastcycle:], '-',
								   color=colour[i], markersize=1)
				# Add z displacement data to save out for central impact
				if folder == 'Central' and column == 0:
					if file.stem == 'c0':
						displacement_export[0] = np.mean(bending[-lastcycle:])
						axs[10+column].plot(time[-1], displacement_export[0], 'x',
								   color=colour[i], markersize=10)
					if file.stem == 'l2':
						displacement_export[1] = np.mean(bending[-lastcycle:])
						axs[10+column].plot(time[-1], displacement_export[1], 'x',
								   color=colour[i], markersize=10)
					if file.stem == 'r2':
						displacement_export[2] = np.mean(bending[-lastcycle:])
						axs[10+column].plot(time[-1], displacement_export[2], 'x',
								   color=colour[i], markersize=10)
					if file.stem == 't2':
						displacement_export[3] = np.mean(bending[-lastcycle:])
						axs[10+column].plot(time[-1], displacement_export[3], 'x',
								   color=colour[i], markersize=10)

				# Plot fifth row if bin detected: Pure bending along mid vein
				if file.stem in ['b2', 'b1', 'c0', 't1', 't2']:
					if len(binfreq) > 1:
						# Plot bending in overview plot
						axs[10+column].plot(time[start:], bending, '-',
									   color=colour[i], markersize=1)
			
						if column == 0 and (file.stem == 'c0' or file.stem == 't2'):

							timeshift = start + np.argmin(bending)
							cyclelength = int(framesPerSecond/peaks[0]) * 10
							if (timeshift+cyclelength) > len(time):
								cyclelength = len(time) - timeshift

							# Fit for bending data with initial params: maxA, zeta, omega=2pi*f, displ
							init_params = [np.min(bending), 2, 2*np.pi*peaks[0], np.mean(bending[-300:])]
							bend_fit, params = oscillation_fit(time[:cyclelength],
								bending[np.argmin(bending):(np.argmin(bending)+cyclelength)], init_params)

							# Find peaks in data for exponential fit
							#timesegment = time[:cyclelength+3]
							#zsegment = np.abs(bending[np.argmin(bending)-3:(np.argmin(bending)+cyclelength)])

							#zpeaks_oscillation, _ = find_peaks(zsegment, height=0.1)
							#timesegment = timesegment[zpeaks_oscillation]
							#zsegment = zsegment[zpeaks_oscillation]

							#exp_fit, params2 = damping_fit(timesegment,
							#	zsegment, init_params[:-1])

							#axs2[0].plot(timesegment, exp_fit, 'x-',
							#		   color=colour[i], alpha=0.6, markersize=1)

							#print(params[1], params2[1])


							# Plot data for fit and fit for damping coeff
							axs2[0].plot(time[:cyclelength], bending[np.argmin(
							   bending):(np.argmin(bending)+cyclelength)], '.',
							   color=colour[i], markersize=1)
							axs2[0].plot(time[:cyclelength], bend_fit, '-',
									   color=colour[i], alpha=0.6, markersize=1)

							# Add max A and damping coeff data to save out
							if file.stem == 'c0':
								amplitude_export[0] = np.min(bending)
								damping_export[0] = params[1]
								axs2[0].plot(time[1:cyclelength], params[0]*
									np.exp(-params[1]*time[1:cyclelength])+params[3], '--',
									   color=colour[i], alpha=0.6, markersize=0.5)
								axs2[0].plot(time[1:cyclelength], -params[0]*
									np.exp(-params[1]*time[1:cyclelength])+params[3], '--',
									   color=colour[i], alpha=0.6, markersize=0.5)

								axs[10+column].plot(time[start+np.argmin(bending)], 
									amplitude_export[0], 'x', color=colour[i], markersize=10)
							if file.stem == 't2':
								amplitude_export[3] = np.min(bending)
								damping_export[3] = params[1]
								axs2[0].plot(time[1:cyclelength], params[0]*
									np.exp(-params[1]*time[1:cyclelength])+params[3], '--',
									   color=colour[i], alpha=0.6, markersize=0.5)
								axs2[0].plot(time[1:cyclelength], -params[0]*
									np.exp(-params[1]*time[1:cyclelength])+params[3], '--',
									   color=colour[i], alpha=0.6, markersize=0.5)
								axs[10+column].plot(time[start+np.argmin(bending)], 
									amplitude_export[3], 'x', color=colour[i], markersize=10)

				# Plot fifth row if bin detected: Pure twisting
				if file.stem in ['r2', 'r1', 'l1', 'l2', 'c0'] and c0_twisting_present:
					if len(binfreq) > 1:
						twisting = twisting - c0_twisting[:, column]
						axs[13+column].plot(time[start:], twisting, '-',
									   color=colour[i], markersize=1)

						if column == 0 and (file.stem == 'l2' or file.stem == 'r2'):

							if np.abs(np.min(twisting)) > np.abs(np.max(twisting)):
								init_maxA = np.min(twisting)
								startindex = np.argmin(twisting)							
								timeshift = start + np.argmin(twisting)
							else:
								init_maxA = np.max(twisting)
								startindex = np.argmax(twisting)							
								timeshift = start + np.argmax(twisting)
							
							init_params = [init_maxA, 4, 2*np.pi*peaks[1],0]
							cyclelength = int(framesPerSecond/peaks[1]) * 10

							if (timeshift+cyclelength) > len(time):
								cyclelength = len(time) - timeshift

							twist_fit, params = oscillation_fit(time[:cyclelength],
								twisting[startindex:(startindex+cyclelength)], init_params)

							# Find peaks in data for exponential fit
							#timesegment = time[:cyclelength+3]
							#zsegment = np.abs(twisting[startindex-3:(startindex+cyclelength)])

							#zpeaks_oscillation, _ = find_peaks(zsegment, height=0.01)
							#timesegment = timesegment[zpeaks_oscillation]
							#zsegment = zsegment[zpeaks_oscillation]

							#exp_fit, params2 = damping_fit(timesegment,
							#	zsegment, init_params[:-1])

							#axs2[1].plot(timesegment, exp_fit, 'x-',
							#		   color=colour[i], alpha=0.6, markersize=1)

							#print(params[1], params2[1])

							# Plot data for fit and fit for damping coeff
							axs2[1].plot(time[:cyclelength],
								twisting[startindex:(startindex+cyclelength)], '.',
							   color=colour[i], markersize=1)
							axs2[1].plot(time[:cyclelength], twist_fit, '-',
									   color=colour[i], alpha=0.6, markersize=1)

							# Add max A and damping coeff data to save out
							if file.stem == 'l2':
								amplitude_export[1] = init_maxA
								damping_export[1] = params[1]
								axs2[1].plot(time[1:cyclelength], params[0]*
									np.exp(-params[1]*time[1:cyclelength])+params[3], '--',
									   color=colour[i], alpha=0.6, markersize=0.5)
								axs2[1].plot(time[1:cyclelength], -params[0]*
									np.exp(-params[1]*time[1:cyclelength])+params[3], '--',
									   color=colour[i], alpha=0.6, markersize=0.5)

								axs[13+column].plot(time[timeshift], 
									amplitude_export[1], 'x', color=colour[i], markersize=10)
							if file.stem == 'r2':
								amplitude_export[2] = init_maxA
								damping_export[2] = params[1]
								axs2[1].plot(time[1:cyclelength], params[0]*
									np.exp(-params[1]*time[1:cyclelength])+params[3], '--',
									   color=colour[i], alpha=0.6, markersize=0.5)
								axs2[1].plot(time[1:cyclelength], -params[0]*
									np.exp(-params[1]*time[1:cyclelength])+params[3], '--',
									   color=colour[i], alpha=0.6, markersize=0.5)
								axs[13+column].plot(time[timeshift], 
									amplitude_export[2], 'x', color=colour[i], markersize=10)

				#if file.stem in ['b2', 'b1', 'c0', 't1', 't2']:
				axs[16+column].plot(time[start:start+int(len(time)/10)], 
						highertwisting[:int(len(time)/10)], '-', color=colour[i], markersize=1)
				
				if file.stem == 't2' and column == 0:
					amplitude_export[3] = np.min(bending)
					damping_export[3] = params[1]
					newstart = start+np.argmin(bending)
					#print(start, time[start:start+cyclelength].shape, bending[np.argmin(
					#		   bending):(np.argmin(bending)+cyclelength)].shape)
					axs[19+column].plot(time[newstart:newstart+cyclelength], bending[np.argmin(
							   bending):(np.argmin(bending)+cyclelength)], '.',
							   color=colour[i], markersize=1)
					axs[19+column].plot(time[newstart:newstart+cyclelength], bend_fit, '-',
									   color=colour[i], alpha=0.6, markersize=1)
					axs[19+column].plot(time[newstart+1:newstart+cyclelength], params[0]*
						np.exp(-params[1]*time[1:cyclelength])+params[3], '--',
					   color=colour[i], alpha=0.6, markersize=0.5)
					axs[19+column].plot(time[newstart+1:newstart+cyclelength], -params[0]*
						np.exp(-params[1]*time[1:cyclelength])+params[3], '--',
					   color=colour[i], alpha=0.6, markersize=0.5)
					axs[19+column].plot(time[start+np.argmin(bending)], 
						amplitude_export[3], 'x', color=colour[i], markersize=10)
					
				#	axs[19+column].plot(time[start:start+int(len(time)/10)], 
				#		highertwisting[:int(len(time)/10)], '-', color=colour[i], markersize=1)

		#break

	# Adjust layout and save plot
	#plt.tight_layout()
	plt.subplots_adjust(left=0.02, bottom=0.01, right=0.98, top=0.99, wspace=0.25, hspace=0.35)
	#fig2.show()
	#plt.show()
	figname = analysedPath / "ValidData_{}.pdf".format(folder)
	fig.savefig(figname)
	plt.close(fig)
	plt.subplots_adjust(left=0.1, bottom=0.05, right=0.95, top=0.90, hspace=0.3)
	figname2 = analysedPath / "ValidDampingFit_{}.pdf".format(folder)
	fig2.savefig(figname2)
	plt.close(fig2)

	return displacement_export, c0_twisting_present, amplitude_export, damping_export

# Open 3D Coordinate files and perform analysis on all files
def analysisComplete(speciesFolder):

	# Prep folders and paths for analysis
	rawDataDir, datafolders, calibfolders, analysedPath = findFoldersAndPaths(speciesFolder)

	# Choose folder from calibration folders to open files
	for folder in calibfolders:
		# Get all files in folder and save paths in files
		folderpath = rawDataDir / folder
		# Make filenamelist of files within folder
		files = sorted(list(folderpath.glob('*.csv')))

		# Loop over all files in folder
		for file in files:
			# Tasks for calibration markers
			pass

	# Set default start values for oscillation, will be replaced if found in overviewFig
	start = np.zeros(3, int) + 50
	# Initialise variables for saving out for each impact location
	bin1 = np.full((10,3), np.nan)
	bin2 = np.full((10,3), np.nan)
	peak1 = np.full((10,3), np.nan)
	peak2 = np.full((10,3), np.nan)
	peak3 = np.full((10,3), np.nan)
	peak4 = np.full((10,3), np.nan)
	peak5 = np.full((10,3), np.nan)
	binfig, binax = setupBinFig()
	# Choose folder from data folders to open files
	for i, folder in enumerate(datafolders):

		# Get all files in folder and save paths in files
		folderpath = rawDataDir / folder
		files = sorted(list(folderpath.glob('*.csv')))

		# Find transformation parameter to align leaf and coordinate system
		rotMat, flipVec, translVec, scaleFac = transformParams(folderpath)

		# Plot and save overview over all data per impact location
		start[i] = overviewFig(folder, files, rotMat, flipVec, translVec,
								scaleFac, analysedPath, binax)


		# Find data to analyse
		validMarkers, validMarkerNames, colour = findValidMarker(files)

		# Is there any data to analyse
		if not validMarkers:
			print('No data with good enough precision for analysis in folder {}.'.format(folder))
			#data_output = '''{0},,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n'''.format(
				#speciesFolder.stem)
			#return data_output
			continue
		# Find common bins for filtering oscillation and peak frequencies
		# Bins can be any length between 0 and 3 indices
		bin1[:,i], bin2[:,i], peak1[:,i], peak2[:,i], peak3[:,i], peak4[:,i], peak5[:,i] = findFftBins(validMarkers, validMarkerNames, flipVec, rotMat, translVec)
	
	# Combine all found bins and peaks into one unique set
	if not np.isnan(bin1).all():
		meanbin1 = np.nanmedian(bin1)
	else:
		meanbin1 = np.nan
	if not np.isnan(bin2).all():
		meanbin2 = np.nanmedian(bin2)
	else:
		meanbin2 = np.nan
	binfreq = np.array([1, meanbin1, meanbin1, meanbin2]) # Hz
	#if not np.isnan(peak2).all():
	#	meanpeak2 = np.nanmedian(peak2)
	#else:
	#	meanpeak2 = np.nan
	#if not np.isnan(peak3).all():
	#	meanpeak3 = np.nanmedian(peak3)
	#else:
	#	meanpeak3 = np.nan

	#peakfreq = [np.nanmedian(peak1), meanpeak2, meanpeak3]
	individualpeaks = find_individual_peaks(np.stack([peak1, peak2, peak3, peak4, peak5]), threshold=2)
	
	# Plot bins as vertical lines
	for i,bfreq in enumerate(binfreq):
		if i < 2:
			col = 'navy'
		else:
			col = 'deeppink'
		binax[0].axvline(bfreq, color = col, alpha = 0.5)
		binax[1].axvline(bfreq, color = col, alpha = 0.5)
		binax[2].axvline(bfreq, color = col, alpha = 0.5)
		font = {'color':  'black', 'weight': 'normal', 'size': 10,}
		binax[2].text(bfreq+0.2, 1000, np.around(bfreq,4), ha='left', fontdict=font)


	# Plot peaks in graph
	for peak in individualpeaks:
		binax[0].plot(peak, -10, '2', color= 'black', markersize=10)
		binax[1].plot(peak, -10, '2', color= 'black', markersize=10)
		binax[2].plot(peak, -10, '2', color= 'black', markersize=10)
		font = {'color':  'black', 'weight': 'normal', 'size': 10,}
		binax[2].text(peak, -400, np.around(peak,2), ha='center', fontdict=font)

	plt.show()
	binfigname = analysedPath / "FFTbins.pdf"
	binfig.savefig(binfigname)
	plt.close(binfig)
		
	newbins = input('Input new bins in format: b1,b2,b3\nDetected bins:{0}\n'.format(
		np.round(binfreq[1:],2)))
	# Take input if any and convert to np array to replace old binfreq
	if newbins:
		newbins = [float(i) for i in newbins.split(',')]
		if len(newbins) == 3:
			binfreq = np.array(newbins)

	newpeaks = input('Input new peaks in format: p1,p2\nDetected peaks:{0}\n'.format(
		np.round(individualpeaks,2)))
	# Take input if any and convert to np array to replace old binfreq
	if newpeaks:
		newpeaks = [float(i) for i in newpeaks.split(',')]
		if len(newpeaks) == 2 and len(individualpeaks)>1:
			individualpeaks[0] = newpeaks[0]
			individualpeaks[1] = newpeaks[1]
		elif len(newpeaks) == 2:
			individualpeaks = np.array((newpeaks[0],newpeaks[1]))


	c0twist = np.full((3,1), np.nan)
	displ_data = np.full((4), np.nan)
	amplitude_export = np.full((3,4), np.nan)
	damping_export = np.full((3,4), np.nan)

	for i, folder in enumerate(datafolders):
		print("Analysing drop impact: {}".format(folder))

		# Get all files in folder and save paths in files
		folderpath = rawDataDir / folder
		files = sorted(list(folderpath.glob('*.csv')))

		# Find transformation parameter to align leaf and coordinate system
		rotMat, flipVec, translVec, scaleFac = transformParams(folderpath)

		# Do analysis for certain information and write out data
		displacement_export, c0twist[i, 0], amplitude_export[i, :], damping_export[i, :] = validDataAnalysis(folder, files, rotMat, flipVec, translVec, scaleFac,
						  analysedPath, start[i], binfreq, individualpeaks)

		# Just save data for certain drop impact locations
		if folder == 'Central':
			displ_data = displacement_export
		

		#break

	# Save out variables
	# Species, bending f, twisting f, higher f1, higher f2, higher f3, higher f4, 
	# displacement, maxA and damping Coefficient data
	if len(individualpeaks) < 6:
		individualpeaks = np.append(individualpeaks, np.full(6-len(individualpeaks), np.nan))

	data_output = '''{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37}\n'''.format(
		speciesFolder.stem, individualpeaks[0], individualpeaks[1],
		individualpeaks[2], individualpeaks[3], individualpeaks[4],
		individualpeaks[5],
		displ_data[0], displ_data[1], displ_data[2], displ_data[3], 
		c0twist[0,0], amplitude_export[0,0], amplitude_export[0,1], 
		amplitude_export[0,2], amplitude_export[0,3],
		damping_export[0,0], damping_export[0,1], 
		damping_export[0,2], damping_export[0,3],
		c0twist[1,0], amplitude_export[1,0], amplitude_export[1,1], 
		amplitude_export[1,2], amplitude_export[1,3],
		damping_export[1,0], damping_export[1,1], 
		damping_export[1,2], damping_export[1,3],
		c0twist[2,0], amplitude_export[2,0], amplitude_export[2,1], 
		amplitude_export[2,2], amplitude_export[2,3],
		damping_export[2,0], damping_export[2,1], 
		damping_export[2,2], damping_export[2,3], 
		binfreq[0], binfreq[1], binfreq[2])

	return data_output

# --------------------------------------------------------------------------------- #

# ------------------------- MAIN CODE --------------------------------------------- #
def main():
	# Current path
	cwd = Path('').absolute()
	# Get to Tracking folder with subfolders in it
	mainPaths = list(cwd.parent.glob('2020_*'))

	# Make file to save out data
	file = cwd / 'SummaryImpactData.csv'
	if not file.exists():
		with open(file, 'w', newline='') as csvfile:
			header = csv.writer(csvfile)
			header.writerow(['Species', 'f_bend', 'f_twist', 'f3', 'f4', 'f5', 'f6', 
'center_c0_displ', 'center_l2_displ', 'center_r2_displ', 'center_t2_displ',
'center_c0_valid', 'center_c0_maxA', 'center_l2_maxA', 'center_r2_maxA', 'center_t2_maxA',
'center_c0_damping', 'center_l2_damping', 'center_r2_damping', 'center_t2_damping',
'lateral_c0_valid', 'lateral_c0_maxA', 'lateral_l2_maxA', 'lateral_r2_maxA', 'lateral_t2_maxA',
'lateral_c0_damping', 'lateral_l2_damping', 'lateral_r2_damping', 'lateral_t2_damping',
'tip_c0_valid', 'tip_c0_maxA', 'tip_l2_maxA', 'tip_r2_maxA', 'tip_t2_maxA',
'tip_c0_damping', 'tip_l2_damping', 'tip_r2_damping', 'tip_t2_damping', 'cutoff_bending',
'cutoff_twist1', 'cutoff_twist2'])


	# Analyse species by species
	for folder in sorted(mainPaths):
		print('Start analysis of: {}'.format(folder.stem))
		speciesFolder = cwd.parent / folder.stem
		dataline = analysisComplete(speciesFolder)

		# Add row of data output to summary document
		summary_file = open('SummaryImpactData.csv', 'a')
		summary_file.write(dataline)
		summary_file.close()

if __name__=="__main__":
	main()
