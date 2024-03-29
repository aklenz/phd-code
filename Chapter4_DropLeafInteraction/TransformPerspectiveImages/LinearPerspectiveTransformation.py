# Linear Perspective Transformation Script
#
# This Script performs a linear perspective transformation on images (.jpg) and can be used
# to deskew images that have not been taken with the camera being parallel to the 
# object that ahould be measured. The user needs to be able to identify 4 corners of
# a rectangle or square that is skewed because of the perspective of the image. 
#
# The script was written by Anne-Kristin Lenz on 21st of January 2021
#
# Instructions when image opens:  
# - Left mouse button click: select coordinate 
#				(on four corners that should lie on a rectangle or better square)
# - Right mouse button click: delete the last selected coordinate
# - Any key: finish selecting coordinates
# - ESC key: stop the script immediately - current image will not be transformed

# ------------------------- IMPORT ------------------------------------------------ #
import cv2
import time
import glob
import numpy as np 
from pathlib import Path
# --------------------------------------------------------------------------------- #


# ------------------------- VARIABLES --------------------------------------------- #
# Possible path leading to images to be transformed 
# or None, if folder name should be used to locate images in same directory
imagefolder_path = None
# Alternatively folder name in which images can be found 
# (has to be in same directory as this script)
imgfolder_name = 'Images'
# Variable to decide if square or rectangle tranform should be chosen
# square performs transform with x and y to same scale 
# 		(important, if area or diagonal lengths should be measured)
# rectangle only performs perspective transform x and y might have different scale
transForm = 'square' # or 'rectangle' # 

# Global variables for function mouseclick to save selected points and 
# only choose a maximum of 4
coords = np.zeros([4,2])
i = 0
# --------------------------------------------------------------------------------- #


# ------------------------- HELP FUNCTIONS ---------------------------------------- #
# Create list of images that need to be transformed 
# (Images in Images folder - images in TransformedImages folder)
def listImages(imagefolder_path, imgfolder_name):
	# Get current working directory if path is not given
	if not imagefolder_path:
		imagefolder_path = Path('').absolute() / imgfolder_name

	# Load all images in current working directory
	images = glob.glob('{}/*.jpg'.format(imagefolder_path))

	# Check if images were found
	if not images:
		print('No images to transform. Check if correct path to images is given')
		quit()

	# If TransformedImages folder exists, check which Images have already been transformed
	write_dir = Path('').absolute() / 'TransformedImages'

	# Delete already transformed images from images list
	if write_dir.exists():
		# Get all image paths from Transformed folder
		transfImages = glob.glob('{}/*.jpg'.format(write_dir))

		# Create filename lists for the images in folders 'Images' and 'TransformedImages'
		images_names = [Path(f).name for f in sorted(images)]
		transfImages_names = [Path(f).name for f in sorted(transfImages)]

		# Compare lists and write only the ones that don't exist in written_images
		notYetTransformed_names = [name for name in images_names if name not in transfImages_names]

		# Create new path for each image that has not yet been transformed
		images = [ imagefolder_path / name for name in notYetTransformed_names]

	if not images:
		print('No images to transform. All images already in TransformedImages folder.')
		quit()

	return images;

# Events that happen when a mousebutton is clicked
# Adjusted code from (accessed on 20th Jan 2021):
# https://www.geeksforgeeks.org/displaying-the-coordinates-of-the-points-clicked-on-the-image-using-python-opencv/
def mouseclick(event, x, y, flags, params):

	# Grab global variables
	global coords, i

	# Checking for left mouse clicks
	if event == cv2.EVENT_LBUTTONDOWN:

		# If less than 4 coordinates are chosen, save new coordinates
		if i < 4:
			# Save new coordinates into variable
			coords[i,:] = [x,y]

			# Display the coordinates on the image window
			font = cv2.FONT_HERSHEY_SIMPLEX 
			text = '. (' + str(x) + ',' + str(y) + ')'
			cv2.putText(img,text,(x,y),font,2,(200,0,0),7)
			cv2.imshow('Opened Image', img)
			# Set counter +1
			i += 1

			# If 4 coordinates are chosen, display message how to close window
			if i == 4:
				# Display the coordinates on the image window
				font = cv2.FONT_HERSHEY_SIMPLEX 
				text = 'If you are happy with all positions,' 
				text2 = 'press any key to close window.'
				cv2.putText(img,text,(1100, 1200),font,4,(13,150,200),8)
				cv2.putText(img,text2,(1100, 1400),font,4,(13,150,200),8)
				cv2.imshow('Opened Image', img)

		# If already 4 coordinates are chosen
		else:
        	# Display a message on the screen to instruct how to delete points
			font = cv2.FONT_HERSHEY_SIMPLEX
			text = 'Already 4 points chosen.'
			text2 = 'Please first delete earlier selected points by right clicks'
			cv2.putText(img,text,(100,100),font,4,(50,50,50),8)
			cv2.putText(img,text2,(100,280),font,4,(50,50,50),8)
			cv2.imshow('Opened Image', img)

	# Checking for right mouse clicks
	if event == cv2.EVENT_RBUTTONDOWN:

		# Delete the last pair of coordinates if there are any saved
		if i > 0:
			# Save new coordinates into variable
			x,y = int(coords[i-1,0]), int(coords[i-1,1])
			coords[i-1,:] = 0

			# Display image with grey coordinates
			font = cv2.FONT_HERSHEY_SIMPLEX 
			text = '. (' + str(x) + ',' + str(y) + ')'
			cv2.putText(img,text,(int(x),int(y)),font,2,(200,200,200),7)
			cv2.imshow('Opened Image', img)
			# Set counter -1
			i -= 1

		# If no coordinates exist (anymore)
		else:
        	# Display a message on the screen
			font = cv2.FONT_HERSHEY_SIMPLEX 
			text = 'No points to delete. First select points by left clicks.'
			cv2.putText(img,text,(100,2800),font,4,(50,50,50),10)
			cv2.imshow('Opened Image', img)
			
	return;

# Select coordinates on image
def selectCoords(img):

	# Set window properties to full screen
	cv2.namedWindow('Opened Image', cv2.WND_PROP_FULLSCREEN)
	cv2.setWindowProperty('Opened Image' ,cv2.WND_PROP_FULLSCREEN,cv2.WINDOW_FULLSCREEN)
	
	# Open image as full screen
	cv2.imshow('Opened Image', img)

	# Wait for mouse clicks in the opened window and perform actions according to 
	# function mouseclick
	cv2.setMouseCallback('Opened Image', mouseclick)

	# Wait for a key to be pressed to close window
	k = cv2.waitKey(0)

	# Only if key pressed is ESC the script will be stopped immediately
	if k == 27:
		quit()

	# Close the window
	cv2.destroyAllWindows()

	return coords;

# Set counter of loop one back, so the selection process for the current image is 
# started again
def redoCurrentImage():

	# Print message in shell to explain how to select correct coordinates
	print('Not enough coordinates chosen to perform image transformation!')
	print('Please choose 4 corners of a rectangle.')
	print('Try again in 3... ', end='', flush=True)
	time.sleep(1)
	print('2... ', end='', flush=True)
	time.sleep(1)
	print('1... \n')
	time.sleep(1)

	# Reset global coordinates and i
	global coords
	coords = np.zeros([4,2])
	global i
	i = 0

	return;

# Sort coordinates so they always match the order in points2
# 0 = upper left, 1 = lower left, 2 = upper right, 3 = lower right
def sortCoords(coords):

	# Make new variable to write sorted coordinates
	sorted_coords = np.zeros([4,2])
	
	# Make list with sum of xs and ys
	sumXY = [coords[i,0]+coords[i,1] for i in range(len(coords))]
	# Make list with subracted ys from xs
	subXY = [coords[i,0]-coords[i,1] for i in range(len(coords))]

	# Find all corners of the rectangle and sort them
	for i in range(len(coords)):
		# upper left: sum of x and y is smallest
		if sumXY[i] == np.min(sumXY):
			sorted_coords[0,:] = coords[i,0], coords[i,1]
		# lower right: sum of x and y is biggest
		elif sumXY[i] == np.max(sumXY):
			sorted_coords[3,:] = coords[i,0], coords[i,1]
		# lower left: x smaller than y
		elif subXY[i] < 0:
			sorted_coords[1,:] = coords[i,0], coords[i,1]
		# upper right: y smaller than x
		else:
			sorted_coords[2,:] = coords[i,0], coords[i,1]

	return sorted_coords;

# Perform transform on image
def performTransform(img, coords, transForm):

	# Sort coordinates into position, so coordinates in points1 and points 2 are 
	# the same corners
	points1 = np.float32(sortCoords(coords))

	if transForm == 'rectangle':
		# Make new coordinates out of points1 that actually lie in the corners of 
		# a rectangle
		points2 = np.float32([	[np.min(points1[:,0]),	np.min(points1[:,1])],
								[np.min(points1[:,0]),	np.max(points1[:,1])],
								[np.max(points1[:,0]),	np.min(points1[:,1])],
								[np.max(points1[:,0]),	np.max(points1[:,1])]	])
	elif transForm == 'square':
		# calculate centre point of the new square
		centrePoint = [round(np.mean(points1[:,0]),0), round(np.mean(points1[:,1]),0)]
		# calculate side length of the new square
		halfSideLength = round(np.max([np.max(points1[:,0])-np.min(points1[:,0]),
						np.max(points1[:,1])-np.min(points1[:,1])])/2,0)
		# Make new coordinates out of points1 that lie on the corners of a square
		points2 = np.float32([	[centrePoint[0] - halfSideLength,	centrePoint[1] - halfSideLength],
								[centrePoint[0] - halfSideLength,	centrePoint[1] + halfSideLength],
								[centrePoint[0] + halfSideLength,	centrePoint[1] - halfSideLength],
								[centrePoint[0] + halfSideLength,	centrePoint[1] + halfSideLength]	])
	else:
		transf_img = None
		print('Choose a transform method by specifiying variable transForm!')
		return transf_img


	# Get size of original image
	y_size,x_size = img.shape[:2]

	# Generate a perspective transform matrix
	matrix = cv2.getPerspectiveTransform(points1,points2)

	# Perfom transformation
	transf_img = cv2.warpPerspective(img, matrix,(x_size,y_size))

	return transf_img;

# Save transformed image in directory 'TransformedImages'
def saveTransfImg(transf_img, filename):

	# Create path for tranformed images folder
	write_dir = Path('').absolute() / 'TransformedImages' 

	# Make directory for writing processed data
	if not write_dir.exists():
		write_dir.mkdir()

	# Check if transform was successfull and only write if successful
	if transf_img.any():
		# Make path for tranformed image
		write_path = write_dir / filename

		if not write_path.exists():
			# Save image
			cv2.imwrite(str(write_path), transf_img)
		else:
			# Display message that image already exsists and was not saved
			print('Image {} already exsists in folder TransformedImages.'.format(filename))
			print('New transformed image was not saved.\n')

	return;


# ------------------------- MAIN CODE --------------------------------------------- #

# Create list of images that need to be transformed
images = listImages(imagefolder_path, imgfolder_name)

# Perform transformation for each image
for file in images:
	# Counter for while loop
	n = 0

	# Redo until 4 coordinates are choosen correctly and transformation can be performed
	while n < 1:
		# Open image to transform
		org_img = cv2.imread(str(file))

		# Make copy of original image to draw on when selecting coordinates
		img = org_img.copy()

		# Select 4 coordinates that should be the corners of a rectangle
		coords = selectCoords(img)

		# Check if four coordinates have been selected
		if not np.all(coords):
			# Less than 4 coordinates, jump back to the beginning of the while loop
			redoCurrentImage()

		# 4 Coordinates selected: continue with transformation
		else: 
			# Transform original image
			transf_img = performTransform(org_img, coords, transForm)

			# Save transformed image
			saveTransfImg(transf_img, Path(file).name)

			# Reset global coordinates and i
			coords = np.zeros([4,2])
			i = 0

			# Set counter +1 so loop ends and new image is loaded
			n +=1