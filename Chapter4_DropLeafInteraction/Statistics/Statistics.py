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
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import linregress
from scipy.stats import t
from scipy.optimize import curve_fit
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.anova import AnovaRM
from scipy.stats import chi2, shapiro
# --------------------------------------------------------------------------------- #

# ------------------------- DATA HEADER ------------------------------------------- #
# Column header of summary impact data
#header = ['Species', 'f_bend', 'f_twist', 'f3', 'f4', 'f5', 'f6', 									#  0- 6  species and freq in Hz
#'center_c0_displ', 'center_l2_displ', 'center_r2_displ', 'center_t2_displ',						#  7-10  center displacement in mm
#'center_c0_valid', 'center_c0_maxA', 'center_l2_maxA', 'center_r2_maxA', 'center_t2_maxA',			# 11-15  center amplitude in mm
#'center_c0_damping', 'center_l2_damping', 'center_r2_damping', 'center_t2_damping',				# 16-19  center damping coefficient
#'lateral_c0_valid', 'lateral_c0_maxA', 'lateral_l2_maxA', 'lateral_r2_maxA', 'lateral_t2_maxA', 	# 20-24  lateral amplitude in mm
#'lateral_c0_damping', 'lateral_l2_damping', 'lateral_r2_damping', 'lateral_t2_damping',			# 25-28  lateral damping coefficient
#'tip_c0_valid', 'tip_c0_maxA', 'tip_l2_maxA', 'tip_r2_maxA', 'tip_t2_maxA',						# 29-33  tip amplitude in mm
#'tip_c0_damping', 'tip_l2_damping', 'tip_r2_damping', 'tip_t2_damping']							# 34-37  tip damping coefficient

# Column header of summary variables
# Species, Width in mm,	Lamina Length in mm, Thickness in mm, Petiole Length in mm,					# 0-4 species and dimensions
# Area in mm2, Contact Angle in degree, Drop Retention Angle in degree,								# 5-7 area and surface properties
# Fresh mass in g, Rehydrated mass in g, Dry mass in g, Ip in gmm2									# 8-10 mass, 11 Ip

# Column header of lever arms
# Species, c0_bend_dist in mm, l2_twist_dist in mm, r2_twist_dist in mm, t2_bend_dist in mm 		# 0-4 species and lever arms
# --------------------------------------------------------------------------------- #

# ------------------------- HELP FUNCTIONS ---------------------------------------- #
# Open impact file and save columns into variables
def openImpactData(filepath):
	filepath = filepath / 'SummaryImpactData_CleanUp.csv'
	# Load data
	species = np.genfromtxt(filepath, delimiter=',', skip_header = 1,
		usecols=0, dtype=str)
	data = np.genfromtxt(filepath, delimiter=',', skip_header = 1,
		unpack= False)
	# Reformat species data
	for i in range(len(species)):
		species[i] = species[i][11:]
	# Sort whole data based on species names
	index = species.argsort()
	species = species[index]
	data = data[index]
	# Sort data into single variables
	freq = np.array(data[:,1:7])
	displ = np.array(data[:, 7:11])
	c0valid = np.array([data[:, 11:12], data[:, 20:21], data[:, 29:30]])
	maxAmp = np.array([data[:, 12:16], data[:, 21:25], data[:, 30:34]])
	damping = np.array([data[:, 16:20], data[:, 25:29], data[:, 34:]])

	return species, freq, displ, c0valid, maxAmp, damping

# Open variable file and save columns into variables
def openVariableData(filepath):
	filepath = filepath / 'SummaryVariables.csv'
	# Load data
	species = np.genfromtxt(filepath, delimiter=',', skip_header = 2,
		usecols=0, dtype=str)
	data = np.genfromtxt(filepath, delimiter=',', skip_header = 2,
		unpack= False)
	# Sort whole data based on species names
	index = species.argsort()
	species = species[index]
	data = data[index]
	# Sort data into single variables
	dims = np.array(data[:, 1:6])
	surface = np.array(data[:, 6:8])
	masses = np.array(data[:, 8:])

	return species, dims, surface, masses

# Open lever file and save columns into variables
def openLeverData(filepath):
	filepath = filepath / 'LeverArms.csv'
	# Load data
	species = np.genfromtxt(filepath, delimiter=',', skip_header = 1,
		usecols=0, dtype=str)
	data = np.genfromtxt(filepath, delimiter=',', skip_header = 1,
		unpack= False)
	# Sort whole data based on species names
	index = species.argsort()
	species = species[index]
	data = data[index]
	# Sort data into single variables
	lever = np.array(data[:, 1:])

	return species, lever

# Open water shed file and save columns into variables
def openWatershedData(filepath):
	filepath = filepath / 'WaterShed.csv'
	# Load data
	species = np.genfromtxt(filepath, delimiter=',', skip_header = 1,
		usecols=0, dtype=str)
	data = np.genfromtxt(filepath, delimiter=',', skip_header = 1,
		unpack= False)
	# Sort whole data based on species names
	index = species.argsort()
	species = species[index]
	data = data[index]
	# Sort data into single variables
	watershed = np.array(data[:, 1:])

	return species, watershed

# Check sorting of all data by comparing order of species
def checkAlignment(species, species_imp, species_lev, species_water):
	# Cut strings to the first 4 letters, so spelling errors won't affect check
	for i in range(len(species)):
		species[i] = species[i][1:5]
		species_imp[i] = species_imp[i][:4]
		species_lev[i] = species_lev[i][:4]
		species_water[i] = species_water[i][1:5]
	# If lists don't match, return error
	if list(species) != list(species_imp) or list(species) != list(species_lev) or list(species) != list(species_water):
		print('Error: Species alignment is off!')
		return False

	return True

# Calculate angles for bending and twisting
def calculateAngles(maxAmp, lever):
	angles = np.arcsin(maxAmp/lever) * 180 / np.pi

	return angles

def calculateCIs(data):
	# Prep for interval calculation
	mean_data = np.mean(data)
	std_data = np.std(data, ddof=1)
	n = len(data)
	ci_level = 0.95
	standard_error = std_data / np.sqrt(n)
	df = n - 1

	# Calculate confidence interval
	ci95 = t.interval(ci_level, df, loc=mean_data, scale=standard_error)

	return ci95

def fit_EI(length, coeff):
	# EI = coeff * L ^ 3
	return coeff * length ** 3

def fit_mass(length, coeff):
	# Mass = coeff * L ^ 2
	return coeff * length ** 2

def fit_GJ(length, coeff):
	# GJ = coeff * L
	return coeff * length

def fit_Ip(length, coeff):
	# GJ = coeff * L
	return coeff * length**4

def fit_exp(length, coeff, exponent):
	# variable = coeff * L ^ exponent
	return coeff * length ** exponent

def exponential_fit(length, variable, init_params, exponent):

	if exponent == 0:
		params, covariance = curve_fit(fit_exp, length, variable, p0=init_params)
		fitVar = fit_exp(length, params[0], params[1])
	elif exponent == 1:
		params, covariance = curve_fit(fit_GJ, length, variable, p0=init_params)
		fitVar = fit_GJ(length, params)
	elif exponent == 2:
		params, covariance = curve_fit(fit_mass, length, variable, p0=init_params)
		fitVar = fit_mass(length, params)
	elif exponent == 3:
		params, covariance = curve_fit(fit_EI, length, variable, p0=init_params)
		fitVar = fit_EI(length, params)
	elif exponent == 4:
		params, covariance = curve_fit(fit_Ip, length, variable, p0=init_params)
		fitVar = fit_Ip(length, params) 
	else:
		print('No fit possible')

	r_squared = calculate_rsq(variable, fitVar)
	s_value = calculate_S(variable, fitVar)

	return params, r_squared, s_value

def calculate_rsq(variable, predictedVar):
	# Calculate R-squared
	residuals = variable - predictedVar
	ss_total = np.sum((variable - np.mean(variable))**2)
	ss_residual = np.sum(residuals**2)
	r_squared = 1 - (ss_residual / ss_total)

	return r_squared

def calculate_S(variable, predictedVar, p=2):

	residuals = variable - predictedVar
	ss_residual = np.sum(residuals**2)
	mse = ss_residual / (len(variable) - p - 1)
	standardErrorOfRegression = np.sqrt(mse)

	return standardErrorOfRegression

def setupFigure(rows, cols, title):
	# Colours for figure
	colour =  ['#666666', '#ffcf33', '#dba400', '#000000', '#96ed1c', '#578e0b', 
	'#f471c8', '#71094f', '#29a6ff', '#0062a8']

	# Make figure
	fig = plt.figure(figsize = (12,rows*3))
	fig.suptitle("{}".format(title))

	axs = []

	for i in range(rows):
		for j in range(cols):
			axs.append(fig.add_subplot(rows, cols, 2*i+j+1))
			#axs[i].set_xlabel('X var')
			#axs[i].set_ylabel('Y var')
			plt.grid(True)

	return fig, axs, colour


def frequencyResponse(species, freq, mass, length, width, polarI):
	# Calculate total leaf length
	leaflength = length[:,0] + length[:,1]
	leaflength = leaflength/1000 		# convert to m
	mass = mass/1000  					# convert to kg
	polarI = polarI/10**9  	    		# convert to kg m2
	# Mask for petiole where petiole measurement was bigger than 5 mm
	maskpeti = length[:,1] > 5
	masknopeti = length[:,1] <= 5

	# Bending and twisiting frequencies
	f_bend = freq[:, 0]
	f_twist = freq[:, 1]

	# Estimated EI and GJ/I_p
	estm_EI = (f_bend/3.5160133)**2 * mass*(leaflength**3) 	# in N m2  
	estm_GJperIp = f_twist**2 * leaflength 					# in 1/s2*m
	estm_GJ = f_twist**2 * leaflength * polarI 				# in N m2
	# Ip was caluclated with g*mm2 so f twist formula is adjusted to this Ip and 
	# does not mean that Ip = Ixx + Iyy!!! It (Bhosale 2020) deviates from 
	# Niklas 1992 twisting formula, which has m in there

	# Least square curve fit on EI = coeff L**3
	init_params = [0.001]
	params, rsq, se = exponential_fit(leaflength, estm_EI, init_params, 3)
	print('Fit EI', params, rsq, se)

	# Least square curve fit on GJ = coeff L**1
	init_params2 = [0.00001]
	params2, rsq2, se2 = exponential_fit(leaflength, estm_GJ, init_params2, 3)
	print('Fit GJ', params2, rsq2, se2)

	# Least square curve fit on Ip = coeff L
	init_params3 = [0.00001, 4]
	params3, rsq3, se3 = exponential_fit(leaflength, polarI, init_params3, 0)
	print('Fit Ip', params3, rsq3, se3)

	# Least square curve fit on M = coeff L**exponent
	init_params4 = [0.001, 2]
	params4, rsq4, se4 = exponential_fit(leaflength, mass, init_params4, 0)
	print('Fit M', params4, rsq4, se4)

	# Theory for estiamted f according to mean EI/(M/L)
	length = np.linspace(30,180,151)/1000
	meanEI = np.nanmean(estm_EI)
	ei_L = params[0]*length**3
	meanM = np.nanmean(mass)
	m_L = params4[0]*length**params4[1]
	theoryBendF1 = 3.5160133 * np.sqrt(meanEI/meanM) * np.sqrt(1/length**3)
	theoryBendF2 = 3.5160133 * np.sqrt(ei_L/meanM) * np.sqrt(1/length**3)
	theoryBendF3 = 3.5160133 * np.sqrt(meanEI/m_L) * np.sqrt(1/length**3)
	theoryBendF4 = 3.5160133 * np.sqrt(ei_L/m_L) * np.sqrt(1/length**3)
	print("Classic bending:", 3.5160133 * np.sqrt(meanEI/meanM), "L^-1.5")
	
	se1 = calculate_S(f_bend, 3.5160133 * np.sqrt(meanEI/meanM) * np.sqrt(1/leaflength**3))
	se2 = calculate_S(f_bend, 3.5160133 * np.sqrt((params[0]*leaflength**3)/(meanM)) * np.sqrt(1/leaflength**3))
	se3 = calculate_S(f_bend, 3.5160133 * np.sqrt((meanEI)/(params4[0]*leaflength**params4[1])) * np.sqrt(1/leaflength**3))
	se4 = calculate_S(f_bend, 3.5160133 * np.sqrt((params[0]*leaflength**3)/(params4[0]*leaflength**params4[1])) * np.sqrt(1/leaflength**3))

	#print('R^2 f bend models', rsq1, rsq2, rsq3, rsq4)
	print('S f bend models', se1, se2, se3, se4)
	#print("3.5160133 * np.sqrt(({0}*L**3)/({1}*L**{2})) * np.sqrt(1/L**3)".format(params[0],params4[0], params4[1]))


	meanGJ = np.nanmean(estm_GJ)
	gj_L = params2[0]*length**3 #**params2[1]
	meanIp = np.nanmean(polarI)
	ip_L = params3[0]*length**params3[1]
	theoryTwistF1 = np.sqrt(meanGJ/meanIp) * np.sqrt(1/length)
	theoryTwistF2 = np.sqrt(gj_L/meanIp) * np.sqrt(1/length)
	theoryTwistF3 = np.sqrt(meanGJ/ip_L) * np.sqrt(1/length)
	theoryTwistF4 = np.sqrt(gj_L/ip_L) * np.sqrt(1/length)
	print("Classic twisting:", np.sqrt(meanGJ/meanIp), "L^-0.5")
	print("Tweaked twisting:", np.sqrt(params2[0]/params3[0]), "L^", 1/2 * (3-1-params3[1]))
	
	se5 = calculate_S(f_twist, np.sqrt(meanGJ/meanIp) * np.sqrt(1/leaflength))
	se6 = calculate_S(f_twist, np.sqrt((params2[0]*leaflength**3)/meanIp) * np.sqrt(1/leaflength))
	se7 = calculate_S(f_twist, np.sqrt(meanGJ/(params3[0]*leaflength**params3[1])) * np.sqrt(1/leaflength))
	se8 = calculate_S(f_twist, np.sqrt((params2[0]*leaflength**3)/(params3[0]*leaflength**params3[1])) * np.sqrt(1/leaflength))
	print('S f twist models',se5, se6, se7, se8)
	
	#meanBendF = 3.5160133 * np.sqrt((params[0]*length**params[1])/ np.nanmean(mass/leaflength)) / (length**2)

	#bendCoeff = 3.5160133 * np.sqrt(estm_EI/ (mass/leaflength))
	#minBendF = np.nanmin(bendCoeff) / (length**2)
	#meanBendF = np.nanmean(bendCoeff) / (length**2)
	#maxBendF = np.nanmax(bendCoeff) / (length**2)
	#intervalsBendCoeff = calculateCIs(bendCoeff)
	#lowerCIBendF = intervalsBendCoeff[0] / (length**2)
	#upperCIBendF = intervalsBendCoeff[1] / (length**2)

	# Theory for estimated f according to mean GJ/Ip
	twistCoeff = np.sqrt(estm_GJperIp)
	#minTwistF = np.nanmin(twistCoeff) * np.sqrt(1/length)
	meanTwistF = np.nanmean(twistCoeff) * np.sqrt(1/length)
	#maxTwistF = np.nanmax(twistCoeff) * np.sqrt(1/length)
	#intervalsTwistCoeff = calculateCIs(twistCoeff)
	#lowerCITwistF = intervalsTwistCoeff[0] * np.sqrt(1/length)
	#upperCITwistF = intervalsTwistCoeff[1] * np.sqrt(1/length)

	# Linregress for twist to bend
	slope, intercept, r_value, p_value, std_err = linregress(f_twist, f_bend)
	#slope2, intercept2, r_value2, p_value2, std_err2 = linregress(estm_GJ, estm_EI)

	slope2, intercept2, r_value2, p_value2, std_err2 = linregress(np.log10(estm_GJ), np.log10(estm_EI))
	
	print('Lin reg f bend to f twist', slope, intercept, r_value**2, p_value)
	print('Lin reg EI to GJ', slope2, intercept2, r_value2**2, p_value2)
	
	# Setup figure
	fig, axs, colour = setupFigure(3, 2, 'Bending and twisting frequencies')

	# Plot theory and extras
	#axs[0].fill_between(length, maxBendF, minBendF, color=colour[8], alpha=0.15)
	#axs[0].fill_between(length, upperCIBendF, lowerCIBendF, color=colour[8], alpha=0.15)
	axs[0].plot(length, theoryBendF1, '--', color=colour[0], alpha=0.7)
	#axs[0].plot(length, theoryBendF2, color=colour[8], alpha=0.7)
	#axs[0].plot(length, theoryBendF3, color=colour[2], alpha=0.7)
	axs[0].plot(length, theoryBendF4, color=colour[9])
	
	axs[2].plot(length, theoryTwistF1, '--', color=colour[0], alpha=0.7)
	#axs[2].plot(length, theoryTwistF2, '.',  color=colour[6], alpha=0.4)
	#axs[2].plot(length, theoryTwistF3, '.', color=colour[5], alpha=0.4)
	axs[2].plot(length, theoryTwistF4, '--', color=colour[7], alpha=0.7)
	
	

	#axs[2].fill_between(length, maxTwistF, minTwistF, color=colour[6], alpha=0.15)
	#axs[2].fill_between(length, upperCITwistF, lowerCITwistF, color=colour[6], alpha=0.15)
	#axs[2].plot(length, meanTwistF, color=colour[6], alpha=0.7)
	#axs[1].axhline(np.nanmean(estm_EI), color=colour[8], alpha=0.7)
	axs[1].plot(length, params[0]*length**3, color=colour[8])
	ax1 = axs[1].twinx()
	ax1.plot(length, params4[0]*length**params4[1], color=colour[2])
	
	#axs[3].axhline(np.nanmean(estm_GJ), color=colour[6], alpha=0.7)
	axs[3].plot(length, params2[0]*length**3, color=colour[6])
	ax3 = axs[3].twinx()
	ax3.plot(length, params3[0]*length**params3[1], color=colour[5])
	axs[4].axline(xy1=(0,0), slope=1, ls='--', color=colour[0], alpha=0.7)
	#axs[4].axline(xy1=(0,0), slope=0.8,ls='--', color=colour[0], alpha=0.7)
	axs[4].axline(xy1=(0,intercept), slope=slope, color=colour[3])

	axs[5].axline(xy1=(1,10**intercept2), slope=slope2, color=colour[0])


	# Plot line of mean bend to twist ratios
	#axs[4].plot(f_twist, f_bend, 'x', color=colour[0])
	#axs[5].plot(estm_GJperIp, estm_EI, 'x', color=colour[0])

	# Plot data
	# Different symbols for petiole and non petiole cutoff at 5 mm
	axs[0].plot(leaflength[maskpeti], f_bend[maskpeti], 'x', color=colour[9], alpha=0.5)
	axs[2].plot(leaflength[maskpeti], f_twist[maskpeti], 'x', color=colour[7], alpha=0.5)
	axs[1].plot(leaflength[maskpeti], estm_EI[maskpeti], 'x', color=colour[8], alpha=0.5)
	axs[3].plot(leaflength[maskpeti], estm_GJ[maskpeti], 'x', color=colour[6], alpha=0.5)
	#axs[4].plot(width[maskpeti]/leaflength[maskpeti], f_bend[maskpeti]/f_twist[maskpeti], 'x', color='black', alpha=0.5)
	#axs[4].plot(leaflength[maskpeti], f_bend[maskpeti]/f_twist[maskpeti], 'x', color='black', alpha=0.5)
	#axs[5].plot(leaflength[maskpeti], estm_EI[maskpeti]/estm_GJ[maskpeti], 'x', color='black', alpha=0.5)
	axs[4].plot(f_twist[maskpeti], f_bend[maskpeti], 'x', color=colour[3], alpha=0.5)
	axs[5].plot(estm_GJ[maskpeti], estm_EI[maskpeti], 'x', color=colour[0], alpha=0.5)

	print("twist to bend", np.mean(estm_EI/estm_GJ))

	axs[0].plot(leaflength[masknopeti], f_bend[masknopeti], '.', color=colour[9], alpha=0.5)
	axs[2].plot(leaflength[masknopeti], f_twist[masknopeti], '.', color=colour[7], alpha=0.5)
	axs[1].plot(leaflength[masknopeti], estm_EI[masknopeti], '.', color=colour[8], alpha=0.5)
	axs[3].plot(leaflength[masknopeti], estm_GJ[masknopeti], '.', color=colour[6], alpha=0.5)
	#axs[4].plot(width[masknopeti]/leaflength[masknopeti], f_bend[masknopeti]/f_twist[masknopeti], '.', color='black', alpha=0.5)
	#axs[4].plot(leaflength[masknopeti], f_bend[masknopeti]/f_twist[masknopeti], '.', color='black', alpha=0.5)
	#axs[5].plot(leaflength[masknopeti], estm_EI[masknopeti]/estm_GJ[masknopeti], '.', color='black', alpha=0.5)
	axs[4].plot(f_twist[masknopeti], f_bend[masknopeti], '.', color=colour[3], alpha=0.5)
	axs[5].plot(estm_GJ[masknopeti], estm_EI[masknopeti], '.', color=colour[0], alpha=0.5)
	
	print(np.max(f_bend), species[np.argmax(f_bend)])
	# Add m/l and Ip to plot 2 and 3
	
	ax1.plot(leaflength[maskpeti], (mass[maskpeti]), 'x', color=colour[2], alpha=0.5)
	ax1.plot(leaflength[masknopeti], (mass[masknopeti]), '.', color=colour[2], alpha=0.5)
	ax1.set_ylabel('Mass in kg')
	ax1.set_yscale('log')

	
	ax3.plot(leaflength[maskpeti], polarI[maskpeti], 'x', color=colour[5], alpha=0.5)
	ax3.plot(leaflength[masknopeti], polarI[masknopeti], '.', color=colour[5], alpha=0.5)
	ax3.set_ylabel('Estimated Ip in kg/m^2')
	ax3.set_yscale('log')

	# Axis labels
	axs[0].set_xlabel('Leaf length in m')
	axs[0].set_ylabel('f bending in Hz')
	axs[2].set_xlabel('Leaf length in m')
	axs[2].set_ylabel('f twisting in Hz')
	axs[1].set_xlabel('Leaf length in m')
	axs[1].set_ylabel('Estimated EI in N*m^2')
	axs[3].set_xlabel('Leaf length in m')
	axs[3].set_ylabel('Estimated GJ in N*m^2')
	axs[4].set_xlabel('F twisting in Hz')
	#axs[4].set_xlabel('Leaf width/length in mm/mm')
	axs[4].set_ylabel('F bending in Hz')
	axs[5].set_xlabel('Estimated GJ in N*m^2')
	axs[5].set_ylabel('Estimated EI in N*m^2')

	# Axis limits
	axs[0].set_xlim([0.03, 0.175])
	axs[0].set_ylim([0, 20])
	axs[2].set_xlim([0.03, 0.175])
	axs[2].set_ylim([0, 80])
	axs[1].set_xlim([0.03, 0.175])
	axs[1].set_yscale('log')
	axs[1].set_ylim([3*10**(-8), 4*10**(-3)])
	ax1.set_ylim([3*10**(-8), 4*10**(-3)])
	#axs[1].set_xscale('log')
	#axs[2].set_ylim([0, ])
	axs[3].set_xlim([0.03, 0.175])
	axs[3].set_yscale('log')
	axs[3].set_ylim([1.5*10**(-9), 2*10**(-5)])
	ax3.set_ylim([1.5*10**(-9), 2*10**(-5)])
	axs[4].set_xlim([0, 80])
	axs[4].set_ylim([0, 20])
	#axs[4].set_xlim([0.03, 0.175])
	#axs[5].set_xlim([0.03, 0.175])
	#axs[5].set_ylim([0, ])
	axs[5].set_xlim([2.5*10**(-8), 4*10**(-5)])
	axs[5].set_ylim([2.5*10**(-8), 4*10**(-5)])
	axs[5].set_xscale('log')
	axs[5].set_yscale('log')

	# Additional info
	axs[4].set_title('f_bend = {0}f_twist + {1}'.format(np.round(slope, 3), np.round(intercept,3)))
	axs[5].set_title('log10(EI) = {0}log10(GJ) + {1}'.format(np.round(slope2, 3), np.round(intercept2,3)))

	plt.subplots_adjust(left=0.1, bottom=0.05, right=0.95, top=0.90, wspace=0.2, hspace=0.3)
	fig.savefig("FrequencyPlot.pdf")

	#plt.show()
	plt.close(fig)

	return


def anglesDamping(species, length, width, angles, damping, mass, freq, polarI):
	# Calculate total leaf length
	leaflength = length[:,0] + length[:,1]
	leaflength = leaflength/1000 # in m
	mass = mass/1000 # in kg
	polarI = polarI/10**9  # in kg m2
	f_bend = freq[:, 0]
	f_twist = freq[:, 1]
	corMLf = 1/(mass*f_bend*leaflength**2)
	corIpf = 1/ (polarI*f_twist)

	# Setup figure
	fig, axs, colour = setupFigure(3, 2, 'Angles and damping coefficient')

	# Masking for only tip markers that are available for all impact locations
	tipangles = np.stack([np.abs(angles[0, :, 3]), np.abs(angles[1, :, 3]),
						  np.abs(angles[2, :, 3])])
	masktip = ~np.isnan(np.mean(tipangles, axis=0))
	tipdamping = np.stack([np.abs(damping[0, :, 3]), np.abs(damping[1, :, 3]),
						   np.abs(damping[2, :, 3])])
	maskdamptip = ~np.isnan(np.mean(tipdamping, axis=0))

	# Masking for only bigger lateral marker that are available for all impact locations
	lateralangles = np.stack([
		np.nanmax([np.abs(angles[0, :, 1]), np.abs(angles[0, :, 2])], axis=0), 
		np.nanmax([np.abs(angles[1, :, 1]), np.abs(angles[1, :, 2])], axis=0), 
		np.nanmax([np.abs(angles[2, :, 1]), np.abs(angles[2, :, 2])], axis=0)])
	masklat = ~np.isnan(np.mean(lateralangles, axis=0))
	lateraldamping = np.stack([
		np.nanmax([np.abs(damping[0, :, 1]), np.abs(damping[0, :, 2])], axis=0), 
		np.nanmax([np.abs(damping[1, :, 1]), np.abs(damping[1, :, 2])], axis=0), 
		np.nanmax([np.abs(damping[2, :, 1]), np.abs(damping[2, :, 2])], axis=0)])
	maskdamplat = ~np.isnan(np.mean(lateraldamping, axis=0))


	#print(np.stack((species, masktip, masklat, maskdamptip, maskdamplat), axis =0).T)
	#print("n for each plot", sum(masktip), sum(masklat), sum(maskdamptip), sum(maskdamplat))

	# Impact location has influence on:
	# log(angle_bend) in interaction with 1/(M L2 fbend)
	# log(angle_twist) in interaction with 1/(Ip ftwist)
	# log(angle_twist) in interaction with 1/(Ip ftwist)
	# --> Make linear regression for each location
	# Differences in bending damping coefficients
	# For Lateral and Central log(damping) in bending
	# But difference to Tip log(damping) in bending
	# No difference in twisting damping coefficients 
	# --> Linear regression for grouped bending(Central+Lateral), 
	# 		bending(Tip), twisting(all locations) on log(angle)

	print('Data for lin regression: slope, intercept, R^2, p value, std error')
	# Bending damping coefficient grouped central and lateral
	groupedDampingBend = np.mean(np.log(tipdamping[:2, maskdamptip]), axis=0)
	slope3, intercept3, r_value3, p_value3, std_err3 = linregress(leaflength[maskdamptip], groupedDampingBend)
	axs[4].axline(xy1=(0,np.exp(intercept3)), slope=slope3/np.exp(1), ls='-', color=colour[3], alpha=0.7)
	print('Damping bend, grouped central and lateral', slope3, intercept3, r_value3**2, p_value3, std_err3)
	# Bending damping coefficient tip
	slope3, intercept3, r_value3, p_value3, std_err3 = linregress(leaflength[maskdamptip], np.log(tipdamping[2, maskdamptip]))
	axs[4].axline(xy1=(0,np.exp(intercept3)), slope=slope3/np.exp(1), ls='-', color=colour[9], alpha=0.7)
	print('Damping bend, tip impact', slope3, intercept3, r_value3**2, p_value3, std_err3)

	# Twisting damping coefficient not grouped
	groupedDampingTwist = np.mean(np.log(lateraldamping[:, maskdamplat]), axis=0)
	slope4, intercept4, r_value4, p_value4, std_err4 = linregress(leaflength[maskdamplat], groupedDampingTwist)
	axs[5].axline(xy1=(0,np.exp(intercept4)), slope=slope4/np.exp(1), ls='-', color=colour[3], alpha=0.7)
	print('Damping twist, grouped all locations', slope4, intercept4, r_value4**2, p_value4, std_err4)

	# Bending angle grouped central and lateral
	#groupedBend = np.mean(np.log(tipdamping[:2, masktip]), axis=0)
	#slope, intercept, r_value, p_value, std_err = linregress(np.log(corMLf[masktip]), groupedBend)
	#axs[2].axline(xy1=(1,np.exp(intercept)), slope=slope, ls='-', color=colour[3], alpha=0.7)
	#print('Angle bend, grouped central and lateral', slope, intercept, r_value**2, p_value, std_err)

	# Plot data
	marker = ['.', '|', '4']
	alphas = [0.7, 0.7, 0.7]
	#bendcols = [colour[3], colour[8], colour[9]]
	#twistcols = [colour[3], colour[6], colour[7]]
	bendcols = [colour[0], colour[2], colour[9]]
	twistcols = [colour[0], "#73bb11", colour[7]]
	#linemarker = ['-', '--', ':']
	linemarker = ['-', '-', '-']

	for loc in range(3): 

		# Lin regression for angle plots
		# Bending angle
		#if loc == 2:
		slope, intercept, r_value, p_value, std_err = linregress(np.log(corMLf[masktip]), np.log(tipangles[loc,masktip]))
		axs[2].axline(xy1=(1,np.exp(intercept)), slope=slope, ls=linemarker[loc], color=bendcols[loc], alpha=0.7)
		print('Bending angle, impact location', loc, slope, intercept, r_value**2, p_value, std_err)

		# Twisting angle
		slope2, intercept2, r_value2, p_value2, std_err2 = linregress(np.log(corIpf[masklat]), np.log(lateralangles[loc,masklat]))
		# Tip impac is not significant different from mean and has too much variation, therefore no regression line in plot
		if loc < 2:
			axs[3].axline(xy1=(1,np.exp(intercept2)), slope=slope2, ls=linemarker[loc], color=twistcols[loc], alpha=0.7)
		if loc == 2:
			axs[3].axhline(np.mean(lateralangles[loc,masklat]), ls="--", color=twistcols[loc], alpha=0.3)
		print('Twisting angle, impact location', loc, slope2, intercept2, r_value2**2, p_value2, std_err2)

		# Plot data for all leaves with 3 impact locations measured for same marker
		axs[0].plot(leaflength[masktip], tipangles[loc, masktip], marker[loc], color=bendcols[loc], alpha=alphas[loc])
		axs[1].plot(leaflength[masklat], lateralangles[loc, masklat], marker[loc], color=twistcols[loc], alpha=alphas[loc])
		axs[2].plot(corMLf[masktip], tipangles[loc, masktip], marker[loc], color=bendcols[loc], alpha=alphas[loc])
		axs[3].plot(corIpf[masklat], lateralangles[loc, masklat], marker[loc], color=twistcols[loc], alpha=alphas[loc])
		axs[4].plot(leaflength[maskdamptip], tipdamping[loc, maskdamptip], marker[loc], color=bendcols[loc], alpha=alphas[loc])
		axs[5].plot(leaflength[maskdamplat], lateraldamping[loc, maskdamplat], marker[loc], color=twistcols[loc], alpha=alphas[loc])

		# All data no masking
		#axs[0].plot(leaflength, angles[loc, :, 0], marker[loc], color=colour[3], alpha=0.1)
		#axs[0].plot(leaflength, angles[loc, :, 3], marker[loc], color=colour[9], alpha=0.1)
		#axs[1].plot(width, -np.abs(angles[loc, :, 1]), marker[loc], color=colour[7], alpha=0.1)
		#axs[1].plot(width, -np.abs(angles[loc, :, 2]), marker[loc], color=colour[5], alpha=0.1)
		#axs[4].plot(leaflength, damping[loc, :, 0], marker[loc], color=colour[3], alpha=0.1)
		#axs[4].plot(leaflength, damping[loc, :, 3], marker[loc], color=colour[9], alpha=0.1)
		#axs[5].plot(width, damping[loc, :, 1], marker[loc], color=colour[7], alpha=0.1)
		#axs[5].plot(width, damping[loc, :, 2], marker[loc], color=colour[5], alpha=0.1)

	# Axis labels
	axs[0].set_xlabel('Leaf length in m')
	axs[0].set_ylabel('Bending angle in degree')
	axs[1].set_xlabel('Leaf length in m')
	axs[1].set_ylabel('Twisting angle in degree')
	axs[2].set_xlabel('1/(M L^2 f_bend) in s/kgm^2')
	axs[2].set_ylabel('Bending angle in degree')
	axs[3].set_xlabel('1/(I_p f_twist) in s/kgm^2')
	axs[3].set_ylabel('Twisting angle in degree')
	axs[4].set_xlabel('Leaf length in m')
	axs[4].set_ylabel('Damping coefficient bending in -')
	axs[5].set_xlabel('Leaf length in m')
	axs[5].set_ylabel('Damping coefficient twisting in -')

	# Axis limits
	axs[0].set_xlim([0.03, 0.175])
	axs[0].set_ylim([0, 40])
	axs[1].set_xlim([0.03, 0.175])
	axs[1].set_ylim([0, 40])
	axs[2].set_xlim([2*10**3, 10**6])
	axs[2].set_ylim([5*10**(-1), 5*10**1])
	axs[2].set_xscale('log')
	axs[2].set_yscale('log')
	axs[3].set_xlim([10**5, 2*10**7])
	axs[3].set_ylim([2*10**(-1), 6*10**1])
	axs[3].set_xscale('log')
	axs[3].set_yscale('log')
	axs[4].set_xlim([0.03, 0.175])
	#axs[4].set_ylim([0, 30])
	axs[4].set_yscale('log')
	axs[5].set_xlim([0.03, 0.175])
	#axs[5].set_ylim([0, 50])
	axs[5].set_yscale('log')

	#plt.show()
	fig.savefig("AngleDampingPlot.pdf")

	plt.close(fig)

	return


def surfaceProps(species, displ, lever, surface, mass, fbend, length):
	displ = np.abs(displ)
	ca = surface[:,0]
	dra = surface[:,1]
	mass = mass/1000
	leaflength = length[:,0] + length[:,1]
	leaflength = leaflength/1000 # in m

	lever[:,1] = lever[:,0]
	lever[:,2] = lever[:,0]
	displangle = np.arcsin(displ/lever) * 180 / np.pi

	# Load of water for tip and central
	loadWater = displ[:,3]/1000 * 8/3.5160133 * mass * fbend**2 / 9.81 * 10**(6) # in N = kg m s-2
	loadWater0 = displ[:,0]/1000 * 8/3.5160133 * mass * fbend**2 / 9.81 * 10**(6)# in N = kg m s-2
	#loadWater0 = displ[:,0] * 8 * mass * fbend**2 * 10**(6) / (1000 * 3.5160133 * 9.81) # in N = kg m s-2

	loadWater1 = displ[:,1]/1000 * 8/3.5160133 * mass * fbend**2 / 9.81 * 10**(6)# in N = kg m s-2
	loadWater2 = displ[:,2]/1000 * 8/3.5160133 * mass * fbend**2 / 9.81 * 10**(6)# in N = kg m s-2

	#loadWater = displangle[:,3]/ mass 
	#loadWater0 = displangle[:,0]/ mass
	#loadWater1 = displangle[:,1]/ mass
	#loadWater2 = displangle[:,2]/ mass


	# Masking
	masktip = ~np.isnan(displ[:,3])
	masktip[np.nanargmax(displ[:,3])] = False
	maskcenter = ~np.isnan(displ[:,0])
	maskleft = ~np.isnan(displ[:,1])
	maskright = ~np.isnan(displ[:,2])

	maskall = masktip & maskcenter & maskright & maskleft
	print(sum(masktip), sum(maskcenter), sum(maskleft), sum(maskright), sum(maskall))
	print(np.stack((species, maskcenter), axis =0).T)

	# Setup figure
	fig, axs, colour = setupFigure(3, 2, 'Water residue and surface props')

	# Plot tip marker
	#axs[0].plot(ca[masktip], loadWater[masktip], '.', color=colour[9], alpha=0.7)
	#axs[1].plot(dra[masktip], loadWater[masktip], '.', color=colour[9], alpha=0.7)
	#axs[2].plot(ca[masktip], displangle[masktip,3], '.', color=colour[9], alpha=0.7)
	#axs[3].plot(dra[masktip], displangle[masktip,3], '.', color=colour[9], alpha=0.7)
	
	
	#axs[4].plot(mass[masktip], displangle[masktip,3], '.', color=colour[9], alpha=0.7)
	
	#axs[5].plot(leaflength[masktip], dra[masktip], '.', color=colour[9], alpha=0.7)
	#axs[5].plot(leaflength[masktip], ca[masktip], '.', color=colour[8], alpha=0.7)
	
	# Plot center marker
	axs[0].plot(ca[maskcenter], loadWater0[maskcenter], '.', color=colour[3], alpha=0.7)
	axs[1].plot(dra[maskcenter], loadWater0[maskcenter], '.', color=colour[3], alpha=0.7)
	axs[2].plot(ca[maskcenter], displangle[maskcenter,0], '.', color=colour[3], alpha=0.7)
	axs[3].plot(dra[maskcenter], displangle[maskcenter,0], '.', color=colour[3], alpha=0.7)
	axs[4].plot(leaflength[maskcenter], displangle[maskcenter,0], '.', color=colour[3], alpha=0.7)
	axs[5].plot(ca, dra, '.', color=colour[3], alpha=0.7)
	# Plot left marker
	#axs[0].plot(ca[maskleft], loadWater1[maskleft], '.', color=colour[5], alpha=0.7)
	#axs[1].plot(dra[maskleft], loadWater1[maskleft], '.', color=colour[5], alpha=0.7)
	#axs[2].plot(ca[maskleft], displangle[maskleft,1], '.', color=colour[5], alpha=0.7)
	#axs[3].plot(dra[maskleft], displangle[maskleft,1], '.', color=colour[5], alpha=0.7)
	# Plot right marker
	#axs[0].plot(ca[maskright], loadWater2[maskright], '.', color=colour[7], alpha=0.7)
	#axs[1].plot(dra[maskright], loadWater2[maskright], '.', color=colour[7], alpha=0.7)
	#axs[2].plot(ca[maskright], displangle[maskright,2], '.', color=colour[7], alpha=0.7)
	#axs[3].plot(dra[maskright], displangle[maskright,2], '.', color=colour[7], alpha=0.7)

	#axs[0].set_yscale('log')
	#axs[1].set_yscale('log')

	axs[0].set_ylabel('Water residue in mg')
	axs[0].set_xlabel('Contact angle in degree')
	axs[1].set_ylabel('Water residue in mg')
	axs[1].set_xlabel('Drop retention angle in degree')
	axs[2].set_ylabel('Permanent displacement angle in degree')
	axs[2].set_xlabel('Contact angle in degree')
	axs[3].set_ylabel('Permanent displacement angle in degree')
	axs[3].set_xlabel('Drop retention angle in degree')
	axs[4].set_ylabel('Permanent displacement angle in degree')
	axs[4].set_xlabel('Leaf length in mm')
	axs[5].set_ylabel('Drop retention angle in degree')
	axs[5].set_xlabel('Contact angle in degree')

	#plt.show()
	fig.savefig("WaterResidue.pdf")

	plt.close(fig)

	return


def waterShedCor(species, watershed, fbend, mass, length, width, surface):
	spilledge = watershed[:,0]
	ejectdrop = watershed[:,1]
	retractspill = watershed[:,2]
	runoff = watershed[:,3]
	roughness = watershed[:,4]

	ca = surface[:,0]
	dra = surface[:,1]
	mass = mass/1000
	leaflength = length[:,0] + length[:,1]
	leaflength = leaflength/1000 # in m
	#width = width/1000 # in m
	estm_EI = (fbend/3.5160133)**2 * mass*(leaflength**3)

	# Define categories and count occurrences
	categories_eject = [0, 1]
	categories_roughness = [0, 0.5, 1, 2]
	counts_roughness_0 = [np.sum((ejectdrop == 0) & (roughness == cat)) for cat in categories_roughness]
	counts_roughness_1 = [np.sum((ejectdrop == 1) & (roughness == cat)) for cat in categories_roughness]

	fig, axs, colour = setupFigure(1, 5, 'Water shedding correlations')

	axs[0].plot(spilledge, width, '.', color=colour[3], alpha=0.5)
	axs[0].violinplot(width[spilledge==1], [1], showmedians=True)
	axs[0].violinplot(width[spilledge==0], [0], showmedians=True)

	bar_width = 0.35
	bottom = [0,0]
	for i in range(len(categories_roughness)):
		counts = [counts_roughness_0[i], counts_roughness_1[i]]
		axs[1].bar(categories_eject, counts, bar_width, bottom = bottom, color=colour[i*2])
		bottom = [bottom[0] + counts_roughness_0[i], bottom[1] + counts_roughness_1[i]]
		

	#axs[1].plot(ejectdrop, ca, '.', color=colour[3], alpha=0.5)


	axs[2].plot(retractspill, ca, '.', color=colour[3], alpha=0.5)
	axs[2].violinplot(ca[retractspill==1], [1], showmedians=True)
	axs[2].violinplot(ca[retractspill==0], [0], showmedians=True)

	#axs[2].plot(retractspill, dra, '.', color=colour[5], alpha=0.7)
	axs[3].plot(runoff, estm_EI, '.', color=colour[3], alpha=0.5)
	axs[3].violinplot(estm_EI[runoff==1], [1], showmedians=True)
	axs[3].violinplot(estm_EI[runoff==0], [0], showmedians=True)

	axs[4].plot(runoff, dra, '.', color=colour[3], alpha=0.5)
	axs[4].violinplot(dra[runoff==1], [1], showmedians=True)
	axs[4].violinplot(dra[runoff==0], [0], showmedians=True)



	#axs[3].plot(runoff, leaflength, '.', color=colour[3], alpha=0.7)

	#axs[3].set_yscale('log')

	axs[0].set_ylabel('Leaf width in mm')
	axs[0].set_xlabel('Spill edges')
	axs[0].set_xticks(np.arange(0, 2), labels=["no","yes"])

	axs[1].set_ylabel('Surface roughness')
	axs[1].set_xlabel('Droplet ejection')
	axs[1].set_xticks(np.arange(0, 2), labels=["no","yes"])

	axs[2].set_ylabel('Static contact angle in degree')
	axs[2].set_xlabel('Retract splitting')
	axs[2].set_xticks(np.arange(0, 2), labels=["no","yes"])

	axs[3].set_ylabel('Estimated EI in N*m^2')
	axs[3].set_xlabel('Run off')
	axs[3].set_xticks(np.arange(0, 2), labels=["no","yes"])

	axs[4].set_ylabel('Drop retention angle in degree')
	axs[4].set_xlabel('Run off')
	axs[4].set_xticks(np.arange(0, 2), labels=["no","yes"])

	
	
	#plt.show()
	plt.tight_layout()
	fig.savefig("WaterShed.pdf")
	plt.close(fig)
	return

# --------------------------------------------------------------------------------- #

# ------------------------- MAIN CODE --------------------------------------------- #
def main():
	# Current path
	cwd = Path('').absolute()

	# Load in all files and save into variables
	species_imp, freq, displ, c0valid, maxAmp, damping = openImpactData(cwd)
	species, dims, surface, masses = openVariableData(cwd)
	species_lev, lever = openLeverData(cwd)
	species_water, watershed = openWatershedData(cwd)

	# Shapes of variables: 
	# species 			50,
	# freq 				50, 6
	# dims 				50, 5
	# displ, lever, watershed 50, 4
	# masses			50, 3
	# surfaces			50, 2
	# c0valid 		 	3, 50, 1
	# maxAmp, damping	3, 50, 4

	# If sorting of data doesn't match, stop analysis
	check = checkAlignment(np.copy(species), np.copy(species_imp), np.copy(species_lev), np.copy(species_water))
	# Stop analysis here
	if not check:
		return

	# Calculate bend and twist angles from max amplitude and lever arms
	angles = calculateAngles(maxAmp, lever)

	# No 1: F bending and F twisting
	frequencyResponse(species, freq[:,:2], masses[:,0], dims[:, (1,3)], dims[:, 0], masses[:,-1])

	# No 2: Angles and damping
	anglesDamping(species, dims[:, (1,3)], dims[:, 0], angles, damping, masses[:,0], freq[:,:2], masses[:,-1])

	# No 3: Surface properties
	surfaceProps(species, displ, lever, surface, masses[:,0], freq[:,0], dims[:, (1,3)])

	# No 4: Water shed
	waterShedCor(species, watershed, freq[:,0], masses[:,0], dims[:, (1,3)], dims[:, 0], surface)

if __name__=="__main__":
	main()
