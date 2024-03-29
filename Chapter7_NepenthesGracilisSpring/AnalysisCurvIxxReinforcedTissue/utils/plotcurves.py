""" Plot curvature data

Python module to load in data from various files (curvatures and second moment
of area) and plot it into 4 subplots.

This module requires the following packages:
    * pathlib
    * numpy
    * matplotlib
    * logging

This module contains the following functions:
    * load_files - load up, middle and down curvature files into 3*n*4 array
    * calc_differences - calculate the curvature differences up-middle &
                        middle-down
    * load_moment - load second moment of area and area from file
    * plot_curves - plot figure with 4 subfigures from all data

Exceptions are raised in the following functions:
    * load_files - IOError: if curvature.csv files cannot be found
    * calc_differences - TypeError: if dataarray has not shape 3*n*4
    * load_moment - IOError: if secondmoment.csv file cannot be found
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import scipy.stats as stat
import logging
from scipy.stats import norm


dark=0

# Colour scheme dark
if dark == True:
    plt.style.use('dark_background')
    col_curve_down = '#f5b700'      #'yellow'
    col_curve_middle = '#ababab'    #'darkgray'
    col_curve_up = '#29a6ff'        #'deepskyblue'
    col_annotation = '#ababab'      #'darkgray'
    col_lines = '#ebebeb'           #'lightgray'
    col_area = '#ababab'            #'darkgray'
    col_I = '#b80f80'               #'darkpink'
    col_I_quater = '#f9b8e3'        #'lightpink'
    col_density = '#578e0b'         #'darkgreen'
    col_density_quater = '#cbf68d'  #'lightgreen'
    colors_up = ['#33acff', '#4cb6ff', '#66c0ff', '#80cbff', '#99d5ff', '#b3e0ff', '#cceaff']
    colors_down = [ '#ffcf33','#ffd54c', '#ffdb66','#ffe180', '#ffe799','#ffedb3', '#fff3cc']
else:
    # Colour scheme light
    col_curve_down = '#dba400'      #'yellow'
    col_curve_middle = '#999999'    #'gray'
    col_curve_up = '#0062a8'        #'deepskyblue'
    col_annotation = '#000000'      #'black'
    col_lines = '#141414'           #'darkgray'
    col_area = '#000000'            #'black'
    col_I = '#71094f'               #'darkpink'
    col_I_quater = '#f471c8'        #'lightpink'
    col_density = '#2c4706'         #'darkgreen'
    col_density_quater = '#96ed1c'  #'lightgreen'
    
    colors_up = ['#33acff', '#4cb6ff', '#66c0ff', '#80cbff', '#99d5ff', '#b3e0ff', '#cceaff']
    colors_down = [ '#ffcf33','#ffd54c', '#ffdb66','#ffe180', '#ffe799','#ffedb3', '#fff3cc']
    
    colors_pitcher = ['#999999', '#dba400', '#0062a8', '#71094f', '#2c4706', '#f471c8', '#96ed1c', '#000000', '#578e0b']
    colors_p7 = ['#00355c', '#757575', '#c29300', '#0081db', '#bdbdbd', '#ffcf3d']

markersel = ['^', 's', 'v', 'o', '+'] #Or 6 and 7 for up and down
alphaline = 0.6
lwidth = 2
alphamarker = 1
mksize = 5
mksize_quater = 7#10


# def single_correlation_old(basename, differencearray, momentarray, path):

#     # end index 
#     #end = 4
#     #end = [-4, -4, -4, -6, -3, -6, -4][int(basename[-1])-1]
#     #end = end-5 dk_down = differencearray[-21:, 0]#[:end] dk_up
#      = -1*differencearray[-21:, 1]#[:end]
#     # Take mean of values at index 0 and stack back together index = np.flip
#       (np.hstack((momentarray[0, :5 , 0], np.mean(momentarray[0, 5:7, 0]),
#       momentarray[0, 7: , 0])))#[:end] I_full = np.flip(np.hstack(
#       (momentarray[0, :5 , 1], np.mean(momentarray[0, 5:7, 1]), momentarray
#       [0, 7: , 1])))#[:end] I_full = 1/ I_full I_quater = np.flip(np.hstack(
#       (momentarray[1, :5 , 1], np.mean(momentarray[1, 5:7, 1]), momentarray
#       [1, 7: , 1])))#[:end] I_quater = 1/ I_quater

#     dens_full = 100*momentarray[0, :, 4]/momentarray[0, :, 3] dens_full =
#     np.flip(np.hstack((dens_full[ :5], np.mean(dens_full[5:7]), dens_full
#     [7: ])))#[:end] dens_quater = 100*momentarray[1, :, 4]/momentarray
#     [1, :, 3] dens_quater = np.flip(np.hstack((dens_quater[ :5], np.mean
#     (dens_quater[5:7]), dens_quater[7: ])))#[:end]
    
#     #index = np.flip(momentarray[0, :, 0])[:end]
#     #I_full = np.flip(momentarray[0, :, 1])[:end]
#     #I_quater = np.flip(momentarray[1, :, 1])[:end]
#     #dens_full = np.flip(100*momentarray[0, :, 4]/momentarray[0, :, 3])
#      [:end]
#     #dens_quater = np.flip(100*momentarray[1, :, 4]/momentarray[1, :, 3])
#      [:end]

#     fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots
#     (3, 3)

#     ax1.plot(dk_down, index, markersel[0], color=col_curve_down) ax1.plot
#     (dk_up, index, markersel[2], color=col_curve_up) ax1.set_xlabel
#     ('Change of curvature\nin 1/mm') ax1.set_ylabel('Spine Index')
#     ax1.yaxis.set_major_locator(MultipleLocator
#     (2)) ax1.yaxis.set_major_formatter(FormatStrFormatter('%d'))

#     ax2.plot(I_full, index, markersel[3], color=col_I) ax2.plot
#     (I_quater, index, markersel[4], color=col_I_quater) ax2.set_xlabel
#     ('I\nin mm4') ax2.set_ylabel('Spine Index') ax2.set_xscale
#     ('log') ax2.yaxis.set_major_locator(MultipleLocator
#     (2)) ax2.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    
#     ax3.plot(dens_full, index, markersel[3], color=col_density) ax3.plot
#     (dens_quater, index, markersel[4], color=col_density_quater)
#     ax3.set_xlabel('Proportion of\nreinforced tissue\nin %') ax3.set_ylabel
#     ('Spine Index') ax3.yaxis.set_major_locator(MultipleLocator
#     (2)) ax3.yaxis.set_major_formatter(FormatStrFormatter('%d'))

#     ax4.set_axis_off()

#     ax5.plot(dk_up, I_full, markersel[3], color=col_curve_up) ax5.plot
#     (dk_down, I_full, markersel[3], color=col_curve_down) ax5.set_xlabel
#     ('Change of curvature\nin 1/mm') ax5.set_ylabel('I\nin mm4')

#     ax6.plot(dk_up, dens_full, markersel[3], color=col_curve_up) ax6.plot
#     (dk_down, dens_full, markersel[3], color=col_curve_down) ax6.set_xlabel
#     ('Change of curvature\nin 1/mm') ax6.set_ylabel
#     ('Proportion of\nreinforced tissue\nin %')

#     ax7.set_axis_off()

#     ax8.plot(dk_up, I_quater, markersel[4], color=col_curve_up) ax8.plot
#     (dk_down, I_quater, markersel[4], color=col_curve_down) ax8.set_xlabel
#     ('Change of curvature\nin 1/mm') ax8.set_ylabel('Local I\nin mm4')

#     ax9.plot(dk_up, dens_quater, markersel[4], color=col_curve_up) ax9.plot
#     (dk_down, dens_quater, markersel[4], color=col_curve_down)
#     ax9.set_xlabel('Change of curvature\nin 1/mm') ax9.set_ylabel
#     ('Local proportion of\nreinforced tissue\nin %')

#     #plt.show()

#     # make lin regress - what are requirements for test
#     #stat.pearsonr(x,y)
#     #slope, intercept, rval, pval, stderr = stat.linregress(x,y)
#     #stat.spearmanr(x,y)
#     #stat.kendalltau(x,y) reg1 = stat.linregress(dk_up, I_full) reg2 =
#      stat.linregress(dk_down, I_full)

#     reg3 = stat.linregress(dk_up, dens_full) reg4 = stat.linregress
#     (dk_down, dens_full)

#     reg5 = stat.linregress(dk_up, I_quater) reg6 = stat.linregress
#     (dk_down, I_quater)

#     reg7 = stat.linregress(dk_up, dens_quater) reg8 = stat.linregress
#     (dk_down, dens_quater)
    

#     # plot lin regress ax5.plot(np.sort(dk_up), reg1[0]*np.sort(dk_up)+reg1
#       [1], color=col_curve_up, label='Rsq={0}, p={1}'.format(round(reg1
#       [2]**2, 3), round(reg1[3], 3))) ax5.plot(np.sort(dk_down), reg2
#       [0]*np.sort(dk_down)+reg2[1], color=col_curve_down, label='Rsq={0}, p=
#       {1}'.format(round(reg2[2]**2, 3), round(reg2[3], 3))) ax5.legend()

#     ax6.plot(np.sort(dk_up), reg3[0]*np.sort(dk_up)+reg3
#     [1], color=col_curve_up, label='Rsq={0}, p={1}'.format(round(reg3
#     [2]**2, 3), round(reg3[3], 3))) ax6.plot(np.sort(dk_down), reg4
#     [0]*np.sort(dk_down)+reg4[1], color=col_curve_down, label='Rsq={0}, p=
#     {1}'.format(round(reg4[2]**2, 3), round(reg4[3], 3))) ax6.legend()

#     ax8.plot(np.sort(dk_up), reg5[0]*np.sort(dk_up)+reg5
#     [1], color=col_curve_up, label='Rsq={0}, p={1}'.format(round(reg5
#     [2]**2, 3), round(reg5[3], 3))) ax8.plot(np.sort(dk_down), reg6
#     [0]*np.sort(dk_down)+reg6[1], color=col_curve_down, label='Rsq={0}, p=
#     {1}'.format(round(reg6[2]**2, 3), round(reg6[3], 3))) ax8.legend()

#     ax9.plot(np.sort(dk_up), reg7[0]*np.sort(dk_up)+reg7
#     [1], color=col_curve_up, label='Rsq={0}, p={1}'.format(round(reg7
#     [2]**2, 3), round(reg7[3], 3))) ax9.plot(np.sort(dk_down), reg8
#     [0]*np.sort(dk_down)+reg8[1], color=col_curve_down, label='Rsq={0}, p=
#     {1}'.format(round(reg8[2]**2, 3), round(reg8[3], 3))) ax9.legend()
    
#     # Save figure fig.tight_layout(pad=0.5, w_pad=-2, h_pad=-4)
#       fig.set_size_inches(20, 20) figname = path /
#       (basename + '_res2.png') fig.savefig(figname) plt.close() logging.info
#       ("Plot saved under: {}".format(figname))

#     data = np.vstack((index, dk_up, dk_down, I_full, I_quater, dens_full,
#     dens_quater)).T

#     if len(data) == 20: data = np.vstack((np.full((1,7), np.nan), data))

#     return data


def load_files(basename, path):
    """
    Load in data from three files with index, x and y coordinates and
    curvature and arrange them into 3*n*4 array

    IO Error: Files not found
    """
    # Create all three paths for file to load in
    filepaths = [(path / Path(basename + "_Down_curvature.csv")),
                 (path / Path(basename + "_Middle_curvature.csv")),
                 (path / Path(basename + "_Up_curvature.csv"))]

    # Try if all 3 files are available
    for path in filepaths:
        if not path.is_file():
            raise IOError("File not found: {}!".format(path))
            return np.zeros(0)

    # Create array and load in all 3 files
    dataarray = np.array([np.loadtxt(file, delimiter=',', skiprows=1)
                         for file in filepaths])

    logging.info("Curvature files loaded successfully.")

    return dataarray

def calc_differences(dataarray):
    """
    Calculate the curvature differences and save in new variable

    Type error: Input has wrong shape
    """
    # Test if dataarray has correct input shape
    if not dataarray.shape[0] == 3 and not dataarray.shape[2] == 4:
        raise TypeError("""Data array for calculating differences has wrong
                        shape!""")
    # Make new array and calculate differences
    differencearray = np.array([dataarray[1, :, 3] - dataarray[0, :, 3],
                                dataarray[1, :, 3] - dataarray[2, :, 3]]).T
    # Differences in percentage of curvature
    #differencearray = np.array([(dataarray[1, :, 3] - dataarray[0, :, 3])/dataarray[1, :, 3],
    #                           (dataarray[1, :, 3] - dataarray[2, :, 3])/dataarray[1, :, 3]]).T

    logging.info("Curvature differences calculated successfully.")

    return differencearray

def load_moment(basename, path):
    """Load in .csv-file with I of bw & greyscale image, cross section area"""

    filepaths = [(path / Path("N_gracilis_" + basename[-2:]
                           + "_secondmoment.csv")),
                (path / Path("N_gracilis_" + basename[-2:]
                           + "_s_secondmoment.csv"))]


    # Check that file exists
    if not filepaths[0].is_file() or not filepaths[1].is_file():
        raise IOError("Files not found: {}!".format(filepaths))
        return np.zeros(0)

    momentarray = np.array([np.loadtxt(file, delimiter=',', skiprows=1)
                         for file in filepaths])

    logging.info("Second moment files loaded successfully.")

    return momentarray

def load_conv_moment(basename, path):
    """Load in .csv-file with I of bw & greyscale image, cross section area"""

    filepaths = [(path / Path("N_gracilis_" + basename[-2:]
                           + "_secondmoment.csv")),
                (path / Path("N_gracilis_" + basename[-2:]
                           + "_030_secondmoment.csv")),
                (path / Path("N_gracilis_" + basename[-2:]
                           + "_045_secondmoment.csv")),
                (path / Path("N_gracilis_" + basename[-2:]
                           + "_060_secondmoment.csv")),
                (path / Path("N_gracilis_" + basename[-2:]
                           + "_075_secondmoment.csv")),
                (path / Path("N_gracilis_" + basename[-2:]
                           + "_090_secondmoment.csv")),
                (path / Path("N_gracilis_" + basename[-2:]
                           + "_105_secondmoment.csv")),
                (path / Path("N_gracilis_" + basename[-2:]
                           + "_120_secondmoment.csv")),
                (path / Path("N_gracilis_" + basename[-2:]
                           + "_180_secondmoment.csv"))
                ]


    # Check that file exists
    if not filepaths[0].is_file() or not filepaths[1].is_file():
        raise IOError("Files not found: {}!".format(filepaths))
        return np.zeros(0)

    momentarray = np.array([np.loadtxt(file, delimiter=',', skiprows=1)
                         for file in filepaths])

    logging.info("Second moment files loaded successfully.")

    return momentarray


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

def adjust_momentarray(basename, momentarray):
    #Adjust index for this (tiff stack to curvature index)

    cutoff = np.where(momentarray[0, :, 0] == 0)

    # Adjust lid array indices from slice number to 0 to 5 and flip array
    lidarray = momentarray[:, cutoff[0][0]: , :]
    lidindices = np.array([0, 1, 2, 3, 4, 5])
    lidslices = np.round(lidindices * np.max(momentarray[0, :, 0]) / 5).astype(int)
    newLidIndex = np.zeros(len(lidarray[0, :, 0]))
    for i in range(5):
        newLidIndex[lidslices[i]:lidslices[i+1]+1] = np.linspace(lidindices[i],
            lidindices[i+1], lidslices[i+1]-lidslices[i]+1)
    lidarray[: ,: , 0] = newLidIndex
    lidarray = np.flip(lidarray, axis=1)

    # Adjust body array indices from slice number to 0 to -15
    indices, slices = load_conversion(basename)
    bodyarray = momentarray[:, :-np.min(slices), :]
    newIndex = np.zeros(len(bodyarray[0, :, 0]))
    for i in range(len(indices)-1):
        newIndex[-slices[i]-1:-slices[i+1]] = np.linspace(indices[i], indices[i+1], 
                                                        slices[i]-slices[i+1]+1)
    bodyarray[: ,: , 0] = newIndex

    # Merge arrays
    momentarray = np.concatenate((lidarray,bodyarray), axis=1)
    # Slice array so only rows with integer indices stay
    rowlist = np.where(momentarray[0, :, 0] == np.round(momentarray[0, :, 0]))[0].tolist()
    shortarray = momentarray[:, rowlist, :]

    return momentarray, shortarray

def prep_cor_data(differencearray, momentarray, inverse):

    dk_down = differencearray[-21:, 0]
    dk_up = -1*differencearray[-21:, 1]
    # Take mean of values at index 0 and stack back together
    index = np.flip(np.hstack((momentarray[0, :5 , 0], 
            np.mean(momentarray[0, 5:7, 0]), momentarray[0, 7: , 0])))
    I_full = np.flip(np.hstack((momentarray[0, :5 , 1], 
            np.mean(momentarray[0, 5:7, 1]), momentarray[0, 7: , 1])))
    I_quater = np.flip(np.hstack((momentarray[1, :5 , 1], 
            np.mean(momentarray[1, 5:7, 1]), momentarray[1, 7: , 1])))
    
    #dens_full = 100*momentarray[0, :, 4]/momentarray[0, :, 3]
    # No proportional
    dens_full = momentarray[0, :, 4]
    dens_full = np.flip(np.hstack((dens_full[ :5], np.mean(dens_full[5:7]),
                dens_full[7: ])))
    #dens_quater = 100*momentarray[1, :, 4]/momentarray[1, :, 3]
    # No proportional
    dens_quater = momentarray[1, :, 4]
    dens_quater = np.flip(np.hstack((dens_quater[ :5], np.mean(dens_quater[5:7]),
                dens_quater[7: ])))

    # I normalisation (to scale all pitchers to same I range)
    #I_full = I_full / np.nanmax(I_full)
    #I_quater = I_quater / np.nanmax(I_quater)
    if inverse:
        I_full = 1/ I_full
        I_quater = 1/ I_quater
        dens_full = 1/dens_full
        dens_quater = 1/dens_quater
    
    data = np.vstack((index, dk_up, dk_down, I_full, I_quater, dens_full, dens_quater)).T

    if len(data) == 20:
        data = np.vstack((np.full((1,7), np.nan), data))

    return data

def single_cor(pitcherno, data, path, method, end):

    # Switch to slice data at different positions    
    if end:
        # Indices of lower max dk
        end = np.array([-9, -9, -9, -11, -8, -11, -9])
    else:
        # Index 0
        end = np.array([-5, -5, -5, -5, -5, -5, -5])
        # Index 4
        #end = np.array([-1, -1, -1, -1, -1, -1, -1])

    # Dataarray to variables
    index = data[:, 0][:end[int(pitcherno)-1]]
    dk_down = data[:, 1][:end[int(pitcherno)-1]]
    dk_up = data[:, 2][:end[int(pitcherno)-1]]
    I_full = data[:, 3][:end[int(pitcherno)-1]]
    I_quater = data[:, 4][:end[int(pitcherno)-1]]
    dens_full = data[:, 5][:end[int(pitcherno)-1]]
    dens_quater = data[:, 6][:end[int(pitcherno)-1]]

    # Plot figure with 7 subfigures
    #fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3)
    fig, ((ax5, ax6), (ax8, ax9)) = plt.subplots(2, 2)

    # ax1.plot(dk_down, index, markersel[0], color=col_curve_down)
    # ax1.plot(dk_up, index, markersel[2], color=col_curve_up)
    # ax1.set_xlabel('Change of curvature\nin 1/mm')
    # ax1.set_ylabel('Spine Index')
    # ax1.yaxis.set_major_locator(MultipleLocator(2))
    # ax1.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    # ax2.plot(I_full, index, markersel[3], color=col_I)
    # ax2.plot(I_quater, index, markersel[4], color=col_I_quater)
    # ax2.set_xlabel('I in mm4 \nor 1/I in mm-4')
    # ax2.set_ylabel('Spine Index')
    # ax2.set_xscale('log')
    # ax2.yaxis.set_major_locator(MultipleLocator(2))
    # ax2.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    
    # ax3.plot(dens_full, index, markersel[3], color=col_density)
    # ax3.plot(dens_quater, index, markersel[4], color=col_density_quater)
    # #ax3.set_xlabel('Proportion of reinforced tissue in \% \nor 1/proportion')
    # ax3.set_xlabel('Reinforced tissue in mm2 \nor 1/reinforcement')
    # ax3.set_ylabel('Spine Index')
    # ax3.yaxis.set_major_locator(MultipleLocator(2))
    # ax3.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    # ax4.set_axis_off()

    ax5.plot(dk_up, I_full, markersel[3], color=col_curve_up)
    ax5.plot(dk_down, I_full, markersel[3], color=col_curve_down)
    ax5.set_xlabel('Change of curvature\nin 1/mm')
    ax5.set_ylabel('I in mm4 \nor 1/I in mm-4')

    ax6.plot(dk_up, dens_full, markersel[3], color=col_curve_up)
    ax6.plot(dk_down, dens_full, markersel[3], color=col_curve_down)
    ax6.set_xlabel('Change of curvature\nin 1/mm')
    #ax6.set_ylabel('Proportion of reinforced tissue in \% \nor 1/proportion')
    ax6.set_ylabel('Reinforced tissue in mm2 \nor 1/reinforcement')

    #ax7.set_axis_off()

    ax8.plot(dk_up, I_quater, markersel[4], color=col_curve_up)
    ax8.plot(dk_down, I_quater, markersel[4], color=col_curve_down)
    ax8.set_xlabel('Change of curvature\nin 1/mm')
    ax8.set_ylabel('Local I in mm4 \nor 1/Local I in mm-4')

    ax9.plot(dk_up, dens_quater, markersel[4], color=col_curve_up)
    ax9.plot(dk_down, dens_quater, markersel[4], color=col_curve_down)
    ax9.set_xlabel('Change of curvature\nin 1/mm')
    #ax9.set_ylabel('Proportion of reinforced tissue in \% \nor 1/Local proportion')
    ax9.set_ylabel('Local reinforced tissue in mm2 \nor 1/local reinforcement')

    # Spearman correlation
    if method == "spear":
        # Stats spearman R
        reg1 = stat.spearmanr(dk_up, I_full)
        reg2 = stat.spearmanr(dk_down, I_full)
        ax5.title.set_text('up: rho = {0}, p = {1}\ndown: rho = {2}, p = {3}'.format(
            round(reg1[0],3), round(reg1[1],3), round(reg2[0],3), round(reg2[1],3)))

        reg3 = stat.spearmanr(dk_up, dens_full)
        reg4 = stat.spearmanr(dk_down, dens_full)
        ax6.title.set_text('up: rho = {0}, p = {1}\ndown: rho = {2}, p = {3}'.format(
            round(reg3[0],3), round(reg3[1],3), round(reg4[0],3), round(reg4[1],3)))

        reg5 = stat.spearmanr(dk_up, I_quater)
        reg6 = stat.spearmanr(dk_down, I_quater)
        ax8.title.set_text('up: rho = {0}, p = {1}\ndown: rho = {2}, p = {3}'.format(
            round(reg5[0],3), round(reg5[1],3), round(reg6[0],3), round(reg6[1],3)))

        reg7 = stat.spearmanr(dk_up, dens_quater)
        reg8 = stat.spearmanr(dk_down, dens_quater)
        ax9.title.set_text('up: rho = {0}, p = {1}\ndown: rho = {2}, p = {3}'.format(
            round(reg7[0],3), round(reg7[1],3), round(reg8[0],3), round(reg8[1],3)))

    # Lin regression with pearson r
    # Todo: Are residuals normally distributed
    else:
        reg1 = stat.linregress(dk_up, I_full)
        reg2 = stat.linregress(dk_down, I_full)

        reg3 = stat.linregress(dk_up, dens_full)
        reg4 = stat.linregress(dk_down, dens_full)

        reg5 = stat.linregress(dk_up, I_quater)
        reg6 = stat.linregress(dk_down, I_quater)

        reg7 = stat.linregress(dk_up, dens_quater)
        reg8 = stat.linregress(dk_down, dens_quater)
        
        # plot lin regress
        ax5.plot(np.sort(dk_up), reg1[0]*np.sort(dk_up)+reg1[1], color=col_curve_up,
            label='Rsq={0}, p={1}'.format(round(reg1[2]**2, 3), round(reg1[3], 3)))
        ax5.plot(np.sort(dk_down), reg2[0]*np.sort(dk_down)+reg2[1], color=col_curve_down,
            label='Rsq={0}, p={1}'.format(round(reg2[2]**2, 3), round(reg2[3], 3)))
        ax5.legend()

        ax6.plot(np.sort(dk_up), reg3[0]*np.sort(dk_up)+reg3[1], color=col_curve_up,
            label='Rsq={0}, p={1}'.format(round(reg3[2]**2, 3), round(reg3[3], 3)))
        ax6.plot(np.sort(dk_down), reg4[0]*np.sort(dk_down)+reg4[1], color=col_curve_down,
            label='Rsq={0}, p={1}'.format(round(reg4[2]**2, 3), round(reg4[3], 3)))
        ax6.legend()

        ax8.plot(np.sort(dk_up), reg5[0]*np.sort(dk_up)+reg5[1], color=col_curve_up,
            label='Rsq={0}, p={1}'.format(round(reg5[2]**2, 3), round(reg5[3], 3)))
        ax8.plot(np.sort(dk_down), reg6[0]*np.sort(dk_down)+reg6[1], color=col_curve_down,
            label='Rsq={0}, p={1}'.format(round(reg6[2]**2, 3), round(reg6[3], 3)))
        ax8.legend()

        ax9.plot(np.sort(dk_up), reg7[0]*np.sort(dk_up)+reg7[1], color=col_curve_up,
            label='Rsq={0}, p={1}'.format(round(reg7[2]**2, 3), round(reg7[3], 3)))
        ax9.plot(np.sort(dk_down), reg8[0]*np.sort(dk_down)+reg8[1], color=col_curve_down,
            label='Rsq={0}, p={1}'.format(round(reg8[2]**2, 3), round(reg8[3], 3)))
        ax9.legend()
    
    # Save figure
    fig.tight_layout(pad=0.5, w_pad=-1, h_pad=-1)
    fig.set_size_inches(10, 10)
    figname = path / ('Curve_p' + pitcherno + '_res2.png')
    fig.savefig(figname)
    plt.close()
    logging.info("Plot saved under: {}".format(figname))

    return

def correlation(data, path, method, end):
    # Switch to slice data at different positions
    if end:
        # Indices of lower max dk
        end = np.array([-9, -9, -9, -11, -8, -11, -9])
    else:
        # Index 0
        end = np.array([-5, -5, -5, -5, -5, -5, -5])
        # Index 4
        #end = np.array([-1, -1, -1, -1, -1, -1, -1])

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    # Array to copy data to for pooled regression
    newdata = np.zeros((21*7+np.sum(end)-1, 7))
    start = 0
    # Plot raw data and regression for each pitcher
    for pitcher in range(len(data)):
        # Slice data into pitcher chunks and add sliced data to new array
        subdata = data[pitcher, :, :]
        if pitcher == 2:
           subdata = data[pitcher, 1: , :]
        subdata = subdata[:end[pitcher], :]

        # Add data to one array for overall regression
        newdata[start : start + len(subdata), :] = subdata
        start = start + len(subdata)

        # Plot raw data for each pitcher
        # I_full
        ax1.plot(subdata[:, 1], subdata[:, 3], markersel[3], color=colors_up[pitcher], label=None)
        ax1.plot(subdata[:, 2], subdata[:, 3], markersel[3], color=colors_down[pitcher], label=None)
        # dens_full
        ax2.plot(subdata[:, 1], subdata[:, 5], markersel[3], color=colors_up[pitcher], label=None)
        ax2.plot(subdata[:, 2], subdata[:, 5], markersel[3], color=colors_down[pitcher], label=None)
        # I_quater
        ax3.plot(subdata[:, 1], subdata[:, 4], markersel[4], color=colors_up[pitcher], label=None)
        ax3.plot(subdata[:, 2], subdata[:, 4], markersel[4], color=colors_down[pitcher], label=None)
        # dens_quater
        ax4.plot(subdata[:, 1], subdata[:, 6], markersel[4], color=colors_up[pitcher], label=None)
        ax4.plot(subdata[:, 2], subdata[:, 6], markersel[4], color=colors_down[pitcher], label=None)

        if method != "spear":
            # Lin regression for each pitcher
            reg1 = stat.linregress(subdata[:, 1], subdata[:, 3])
            reg2 = stat.linregress(subdata[:, 2], subdata[:, 3])

            reg3 = stat.linregress(subdata[:, 1], subdata[:, 5])
            reg4 = stat.linregress(subdata[:, 2], subdata[:, 5])

            reg5 = stat.linregress(subdata[:, 1], subdata[:, 4])
            reg6 = stat.linregress(subdata[:, 2], subdata[:, 4])

            reg7 = stat.linregress(subdata[:, 1], subdata[:, 6])
            reg8 = stat.linregress(subdata[:, 2], subdata[:, 6])

            # Add regression line to plot for each pitcher
            ax1.plot(np.sort(subdata[:, 1]), reg1[0]*np.sort(subdata[:, 1])+reg1[1],
                    '--', color=colors_up[pitcher], label=None)
            ax1.plot(np.sort(subdata[:, 2]), reg2[0]*np.sort(subdata[:, 2])+reg2[1],
                    '--', color=colors_down[pitcher], label=None)

            ax2.plot(np.sort(subdata[:, 1]), reg3[0]*np.sort(subdata[:, 1])+reg3[1],
                    '--', color=colors_up[pitcher], label=None)
            ax2.plot(np.sort(subdata[:, 2]), reg4[0]*np.sort(subdata[:, 2])+reg4[1], 
                    '--', color=colors_down[pitcher], label=None)

            ax3.plot(np.sort(subdata[:, 1]), reg5[0]*np.sort(subdata[:, 1])+reg5[1],
                    '--', color=colors_up[pitcher], label=None)
            ax3.plot(np.sort(subdata[:, 2]), reg6[0]*np.sort(subdata[:, 2])+reg6[1],
                    '--',  color=colors_down[pitcher], label=None)

            ax4.plot(np.sort(subdata[:, 1]), reg7[0]*np.sort(subdata[:, 1])+reg7[1],
                    '--', color=colors_up[pitcher], label=None)
            ax4.plot(np.sort(subdata[:, 2]), reg8[0]*np.sort(subdata[:, 2])+reg8[1],
                    '--', color=colors_down[pitcher], label=None)

    data = newdata

    # Spearman correlation for whole data
    if method == "spear":
        # Stats spearman R
        reg01 = stat.spearmanr(data[:, 1], data[:, 3])
        reg02 = stat.spearmanr(data[:, 2], data[:, 3])
        reg03 = stat.spearmanr(data[:, 1], data[:, 5])
        reg04 = stat.spearmanr(data[:, 2], data[:, 5])
        reg05 = stat.spearmanr(data[:, 1], data[:, 4])
        reg06 = stat.spearmanr(data[:, 2], data[:, 4])
        reg07 = stat.spearmanr(data[:, 1], data[:, 6])
        reg08 = stat.spearmanr(data[:, 2], data[:, 6])

        ax1.title.set_text('up: rho = {0}, p = {1}\ndown: rho = {2}, p = {3}'.format(
            round(reg01[0],3), round(reg01[1],3), round(reg02[0],3), round(reg02[1],3)))
        ax2.title.set_text('up: rho = {0}, p = {1}\ndown: rho = {2}, p = {3}'.format(
            round(reg03[0],3), round(reg03[1],3), round(reg04[0],3), round(reg04[1],3)))
        ax3.title.set_text('up: rho = {0}, p = {1}\ndown: rho = {2}, p = {3}'.format(
            round(reg05[0],3), round(reg05[1],3), round(reg06[0],3), round(reg06[1],3)))
        ax4.title.set_text('up: rho = {0}, p = {1}\ndown: rho = {2}, p = {3}'.format(
            round(reg07[0],3), round(reg07[1],3), round(reg08[0],3), round(reg08[1],3)))

    # Linear Regression for whole data
    else:
        # Regression for whole data
        reg01 = stat.linregress(data[:, 1], data[:, 3])
        reg02 = stat.linregress(data[:, 2], data[:, 3])
        reg03 = stat.linregress(data[:, 1], data[:, 5])
        reg04 = stat.linregress(data[:, 2], data[:, 5])
        reg05 = stat.linregress(data[:, 1], data[:, 4])
        reg06 = stat.linregress(data[:, 2], data[:, 4])
        reg07 = stat.linregress(data[:, 1], data[:, 6])
        reg08 = stat.linregress(data[:, 2], data[:, 6])

        # Regression plots for whole data
        ax1.plot(np.sort(data[:, 1]), reg01[0]*np.sort(data[:, 1])+reg01[1],
                '-', color=col_curve_up, label='Rsq={0}, p={1}'.format(
                round(reg01[2]**2, 3), round(reg01[3], 3)))
        ax1.plot(np.sort(data[:, 2]), reg02[0]*np.sort(data[:, 2])+reg02[1],
                '-', color=col_curve_down, label='Rsq={0}, p={1}'.format(
                round(reg02[2]**2, 3), round(reg02[3], 3)))
        ax1.legend()

        ax2.plot(np.sort(data[:, 1]), reg03[0]*np.sort(data[:, 1])+reg03[1],
                '-', color=col_curve_up, label='Rsq={0}, p={1}'.format(
                round(reg03[2]**2, 3), round(reg03[3], 3)))
        ax2.plot(np.sort(data[:, 2]), reg04[0]*np.sort(data[:, 2])+reg04[1], 
                '-', color=col_curve_down, label='Rsq={0}, p={1}'.format(
                round(reg04[2]**2, 3), round(reg04[3], 3)))
        ax2.legend()

        ax3.plot(np.sort(data[:, 1]), reg05[0]*np.sort(data[:, 1])+reg05[1],
                '-', color=col_curve_up, label='Rsq={0}, p={1}'.format(
                round(reg05[2]**2, 3), round(reg05[3], 3)))
        ax3.plot(np.sort(data[:, 2]), reg06[0]*np.sort(data[:, 2])+reg06[1],
                '-', color=col_curve_down, label='Rsq={0}, p={1}'.format(
                round(reg06[2]**2, 3), round(reg06[3], 3)))
        ax3.legend()

        ax4.plot(np.sort(data[:, 1]), reg07[0]*np.sort(data[:, 1])+reg07[1],
                '-', color=col_curve_up, label='Rsq={0}, p={1}'.format(
                round(reg07[2]**2, 3), round(reg07[3], 3)))
        ax4.plot(np.sort(data[:, 2]), reg08[0]*np.sort(data[:, 2])+reg08[1],
                '-', color=col_curve_down, label='Rsq={0}, p={1}'.format(
                round(reg08[2]**2, 3), round(reg08[3], 3)))
        ax4.legend()

    plt.suptitle('Range from spine Index {0} to {1}'.format(np.min(
        data[ :, 0]), np.max(data[ :, 0])))

    # Axis labels
    ax1.set_xlabel('Change of curvature\nin 1/mm')
    ax1.set_ylabel('I in mm4 \nor 1/I in mm-4')
    #ax1.set_ylabel('Scaled I\nin mm4')
    ax2.set_xlabel('Change of curvature\nin 1/mm')
    #ax2.set_ylabel('Proportion of reinforced tissue in \% \nor 1/proportion')
    ax2.set_ylabel('Reinforced tissue in mm2 \nor 1/reinforcement')
    ax3.set_xlabel('Change of curvature\nin 1/mm')
    ax3.set_ylabel('Local I in mm4 \nor 1/Local I in mm-4')
    #ax3.set_ylabel('Scaled Local I\nin mm4')
    ax4.set_xlabel('Change of curvature\nin 1/mm')
    #ax4.set_ylabel('Local proportion of reinforced tissue in \% \nor 1/Local proportion')
    ax4.set_ylabel('Local reinforced tissue in mm2 \nor 1/Local reinforcement')

    plt.show()
    

def overview_plot(basename, dataarray, differencearray, momentarray, shortmoment, path):
    """
    Create plot with 4 subplots.
    1. Pitcher outline (upper, middle, lower lid position)
    2. Curvature for each pitcher outline
    3. Differences of curvature
    4. Second moments
    5. Reinforced tissue
    6. Areas
    """

    # Create main figure   
    fig, (ax1, ax2, ax3, ax5, ax6, ax4) = plt.subplots(1, 6, gridspec_kw={
        'width_ratios': [4, 2, 2, 2, 2, 2]})
    fig.suptitle('Deformation along midrib for pitcher {0}'.format(
                 basename[-1]))

    # Create subplot 1 - Curves along pitcher midrib
    ax1.axis('equal')
    ax1.axhline(dataarray[1, -6, 2], color=col_lines, lw=0.8)
    ax1.plot(dataarray[2, :, 1], dataarray[2, :, 2], color=col_curve_up)
    ax1.plot(dataarray[1, :, 1], dataarray[1, :, 2], color=col_curve_middle)
    ax1.plot(dataarray[0, :, 1], dataarray[0, :, 2], color=col_curve_down)
    ax1.set_xlabel('X Coordinates\n in mm')
    ax1.set_ylabel('Z Coordinates in mm')
    # Annotation: indices to middle curve for orientation    
    indices = dataarray[1, :, :]
    indices = np.flip(indices)[1:-1:2, :]
    indices = np.flip(indices)
    ax1.plot(indices[:, 1], indices[:, 2], 'x', color=col_annotation, markersize=5)
    for i in range(len(indices)):
        ax1.annotate(int(indices[i, 0]), (indices[i, 1] + 1, indices[i, 2]),
                     va="top", ha="left", size=10, color=col_annotation)

    # Create subplot 2 - Curvature for each point on curves
    ax2.axvline(0, color=col_lines, lw=0.8)
    ax2.axhline(0, color=col_lines, lw=0.8)
    ax2.plot(dataarray[2, :, 3], dataarray[2, :, 0], color=col_curve_up, 
        alpha=alphaline, lw=lwidth)
    ax2.plot(dataarray[1, :, 3], dataarray[1, :, 0], color=col_curve_middle, 
        alpha=alphaline, lw=lwidth)
    ax2.plot(dataarray[0, :, 3], dataarray[0, :, 0], color=col_curve_down, 
        alpha=alphaline)
    ax2.plot(dataarray[2, :, 3], dataarray[2, :, 0], markersel[0], 
        color=col_curve_up, alpha=alphamarker, markersize=mksize)
    ax2.plot(dataarray[1, :, 3], dataarray[1, :, 0], markersel[1], 
        color=col_curve_middle, alpha=alphamarker, markersize=mksize)
    ax2.plot(dataarray[0, :, 3], dataarray[0, :, 0], markersel[2], 
        color=col_curve_down, alpha=alphamarker, markersize=mksize)
    ax2.set_xlabel('Curvature\nin 1/mm')
    ax2.set_ylabel('Spine Index')
    ax2.set_ylim(-15.5,5.5)
    ax2.yaxis.set_major_locator(MultipleLocator(2))
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    # Create subplot 3 - Difference of curvatures
    ax3.axvline(0, color=col_lines, lw=0.8)
    ax3.axhline(0, color=col_lines, lw=0.8)
    ax3.plot(differencearray[:, 1], dataarray[2, :, 0], color=col_curve_up, 
        alpha=alphaline, lw=lwidth)
    ax3.plot(differencearray[:, 0], dataarray[0, :, 0], color=col_curve_down, 
        alpha=alphaline, lw=lwidth)
    ax3.plot(differencearray[:, 1], dataarray[2, :, 0], markersel[0], 
        color=col_curve_up, alpha=alphamarker, markersize=mksize)
    ax3.plot(differencearray[:, 0], dataarray[0, :, 0], markersel[2], 
        color=col_curve_down, alpha=alphamarker, markersize=mksize)
    ax3.set_xlabel('Change of curvature\nin 1/mm')
    ax3.set_xlim(-0.05,0.05)
    ax3.set_ylim(-15.5,5.5)
    ax3.yaxis.set_major_locator(MultipleLocator(2))
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    #Create subplot 4 - Area
    ax4.axhline(0, color=col_lines, lw=0.8)
    ax4.plot(momentarray[0, :, 3], momentarray[0, :, 0], color=col_area, 
        alpha=alphaline, lw=lwidth)
    ax4.plot(momentarray[0, :, 4], momentarray[0, :, 0], color=col_density, 
        alpha=alphaline, lw=lwidth)
    ax4.plot(shortmoment[0, :, 3], shortmoment[0, :, 0], markersel[3], 
        color=col_area, alpha=alphamarker, markersize=mksize)
    ax4.plot(shortmoment[0, :, 4], shortmoment[0, :, 0], markersel[4], 
        color=col_density, alpha=alphamarker, markersize=mksize)
    # Quater section
    ax4.plot(momentarray[1, :, 3], momentarray[1, :, 0], color=col_area, 
        alpha=alphaline, lw=lwidth)
    ax4.plot(momentarray[1, :, 4], momentarray[1, :, 0], color=col_density_quater, 
        alpha=alphaline, lw=lwidth)
    ax4.plot(shortmoment[1, :, 3], shortmoment[1, :, 0], markersel[3], 
        color=col_area, alpha=alphamarker, markersize=mksize)
    ax4.plot(shortmoment[1, :, 4], shortmoment[1, :, 0], markersel[4], 
        color=col_density, alpha=alphamarker, markersize=mksize)
    ax4.set_xlabel('Crossectional &\nreinforced tissue area\nin mm2')
    ax4.set_xlim(0,55)
    ax4.set_ylim(-15.5,5.5)
    ax4.yaxis.set_major_locator(MultipleLocator(2))
    ax4.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    # Create subplot 5 - Second moment of area and area
    ax5.axhline(0, color=col_lines, lw=0.8)
    ax5.plot(momentarray[0, :, 1], momentarray[0, :, 0], color=col_I, 
        alpha=alphaline, lw=lwidth)
    ax5.plot(shortmoment[0, :, 1], shortmoment[0, :, 0], markersel[3], 
        color=col_I, alpha=alphamarker, markersize=mksize)
    # Quater section
    ax5.plot(momentarray[1, :, 1], momentarray[1, :, 0], color=col_I_quater, 
        alpha=alphaline, lw=lwidth)
    ax5.plot(shortmoment[1, :, 1], shortmoment[1, :, 0], markersel[4], 
        color=col_I_quater, alpha=alphamarker, markersize=mksize)
    ax5.set_xlabel('I\nin mm4')
    ax5.set_xscale('log')
    ax5.set_xlim(0.01,5000)
    ax5.set_ylim(-15.5,5.5)
    ax5.yaxis.set_major_locator(MultipleLocator(2))
    ax5.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    # Create subplot - Reinforced tissue
    ax6.axhline(0, color=col_lines, lw=0.8)
    ax6.plot(100*momentarray[0, :, 4]/momentarray[0, :, 3], momentarray[0, :, 0],
            color=col_density, alpha=alphaline, lw=lwidth)
    ax6.plot(100*shortmoment[0, :, 4]/shortmoment[0, :, 3], shortmoment[0, :, 0],
            markersel[3], color=col_density, alpha=alphamarker, markersize=mksize)
    # Quater section
    ax6.plot(100*momentarray[1, :, 4]/momentarray[1, :, 3], momentarray[0, :, 0],
            color=col_density_quater, alpha=alphaline, lw=lwidth)
    ax6.plot(100*shortmoment[1, :, 4]/shortmoment[1, :, 3], shortmoment[0, :, 0],
            markersel[4], color=col_density_quater, alpha=alphamarker, 
            markersize=mksize)
    ax6.set_xlabel('Proportion of\nreinforced tissue\nin %')
    ax6.set_xlim(0,90)
    ax6.set_ylim(-15.5,5.5)
    ax6.yaxis.set_major_locator(MultipleLocator(2))
    ax6.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    # Save figure
    fig.tight_layout(pad=0.3, w_pad=-2, h_pad=0)
    fig.set_size_inches(14, 7)
    figname = path / (basename + '_overview_plot.pdf')
    fig.savefig(figname)
    logging.info("Plot saved under: {}".format(figname))

    plt.close()

    return None

def fig_material0(basename, dataarray, differencearray, path):

    # Create main figure   
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, gridspec_kw={
        'width_ratios': [2,1,1]})
    # Create subplot 1 - Curves along pitcher midrib
    ax1.axis('equal')
    ax1.axhline(dataarray[1, -6, 2], color=col_lines, lw=0.8)
    ax1.plot(dataarray[1, :, 1], dataarray[1, :, 2], color=col_curve_middle, lw=lwidth)
    ax1.plot(dataarray[2, :, 1], dataarray[2, :, 2], color=col_curve_up, lw=lwidth)
    ax1.plot(dataarray[0, :, 1], dataarray[0, :, 2], color=col_curve_down, lw=lwidth)
    #ax1.set_xlabel('X Coordinates\n in mm')
    #ax1.set_ylabel('Z Coordinates in mm')
    # Annotation: indices to middle curve for orientation    
    indices = dataarray[1, :, :]
    indices = np.flip(indices)[1:-1:2, :]
    indices = np.flip(indices)
    ax1.plot(indices[:, 1], indices[:, 2], 'x', color=col_annotation, markersize=7)
    for i in range(len(indices)):
        ax1.annotate(int(indices[i, 0]), (indices[i, 1] + 1, indices[i, 2]),
                     va="center", ha="left", size=10, color=col_annotation)

    # Create subplot 2 - Curvature for each point on curves
    ax2.axvline(0, color=col_lines, lw=0.8)
    ax2.axhline(0, color=col_lines, lw=0.8)
    ax2.plot(dataarray[1, :, 3], dataarray[1, :, 0], color=col_curve_middle, 
        alpha=alphaline, lw=lwidth)
    ax2.plot(dataarray[1, :, 3], dataarray[1, :, 0], markersel[1], 
        color=col_curve_middle, alpha=alphamarker, markersize=mksize)
    ax2.plot(dataarray[2, :, 3], dataarray[2, :, 0], color=col_curve_up, 
        alpha=alphaline, lw=lwidth)
    ax2.plot(dataarray[2, :, 3], dataarray[2, :, 0], markersel[0], 
        color=col_curve_up, alpha=alphamarker, markersize=mksize)
    ax2.plot(dataarray[0, :, 3], dataarray[0, :, 0], color=col_curve_down, 
        alpha=alphaline, lw=lwidth)
    ax2.plot(dataarray[0, :, 3], dataarray[0, :, 0], markersel[2], 
        color=col_curve_down, alpha=alphamarker, markersize=mksize)
    #ax2.set_xlabel('Curvature\nin 1/mm')
    #ax2.set_ylabel('Spine Index')
    ax2.set_ylim(-15.5,5.5)
    ax2.yaxis.set_major_locator(MultipleLocator(2))
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    # Create subplot 3 - Difference of curvatures
    ax3.axvline(0, color=col_lines, lw=0.8)
    ax3.axhline(0, color=col_lines, lw=0.8)
    ax3.plot(differencearray[:, 1], dataarray[2, :, 0], color=col_curve_up, 
        alpha=alphaline, lw=lwidth)
    ax3.plot(differencearray[:, 1], dataarray[2, :, 0], markersel[0], 
        color=col_curve_up, alpha=alphamarker, markersize=mksize)
    ax3.plot(differencearray[:, 0], dataarray[0, :, 0], color=col_curve_down, 
        alpha=alphaline, lw=lwidth)
    ax3.plot(differencearray[:, 0], dataarray[0, :, 0], markersel[2], 
        color=col_curve_down, alpha=alphamarker, markersize=mksize)
    #ax3.set_xlabel('Change of curvature\nin 1/mm')
    ax3.set_xlim(-0.05,0.05)
    ax3.set_ylim(-15.5,5.5)
    ax3.yaxis.set_major_locator(MultipleLocator(2))
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    # Save figure
    fig.tight_layout(pad=0.3, w_pad=-0.5, h_pad=0)
    fig.set_size_inches(9, 7)
    figname = path / ('Mat0_' + basename + '.pdf')
    fig.savefig(figname)
    logging.info("Plot saved under: {}".format(figname))

    plt.close()

    return None

def fig_material1(basename, dataarray, differencearray, momentarray, shortmoment, path):

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, gridspec_kw={
        'width_ratios': [3,1,3]})

    # Create subplot 5 - Second moment of area and area
    ax1.axhline(0, color=col_lines, lw=0.8)
    ax1.plot(momentarray[0, :, 1], momentarray[0, :, 0], color=col_I, 
        alpha=alphaline, lw=lwidth)
    ax1.plot(shortmoment[0, :, 1], shortmoment[0, :, 0], markersel[3], 
        color=col_I, alpha=alphamarker, markersize=mksize)
    # Quater section
    ax1.plot(momentarray[1, :, 1], momentarray[1, :, 0], color=colors_pitcher[0], 
        alpha=alphaline, lw=lwidth)
    ax1.plot(shortmoment[1, :, 1], shortmoment[1, :, 0], markersel[4], 
        color=colors_pitcher[0], alpha=alphamarker, markersize=mksize_quater)
    ax1.plot(np.min(momentarray[1, 540:, 1]), momentarray[1, (int(np.where(
              momentarray[1, 540:, 1] == np.min(momentarray[1, 540:, 1]))[0])+ 540), 0],
              markersize=10, color=colors_pitcher[0], marker='x')

    ax1.plot(momentarray[2, :, 1], momentarray[2, :, 0], color=colors_pitcher[1], 
        alpha=alphaline, lw=lwidth)
    ax1.plot(shortmoment[2, :, 1], shortmoment[2, :, 0], markersel[4], 
        color=colors_pitcher[1], alpha=alphamarker, markersize=mksize_quater)
    ax1.plot(np.min(momentarray[2, 540:, 1]), momentarray[2, (int(np.where(
              momentarray[2, 540:, 1] == np.min(momentarray[2, 540:, 1]))[0])+ 540), 0],
              markersize=10, color=colors_pitcher[1], marker='x')

    ax1.plot(momentarray[3, :, 1], momentarray[3, :, 0], color=colors_pitcher[2], 
        alpha=alphaline, lw=lwidth)
    ax1.plot(shortmoment[3, :, 1], shortmoment[3, :, 0], markersel[4], 
        color=colors_pitcher[2], alpha=alphamarker, markersize=mksize_quater)
    ax1.plot(np.min(momentarray[3, 540:, 1]), momentarray[3, (int(np.where(
          momentarray[3, 540:, 1] == np.min(momentarray[3, 540:, 1]))[0])+ 540), 0],
          markersize=10, color=colors_pitcher[2], marker='x')

    ax1.plot(momentarray[4, :, 1], momentarray[4, :, 0], color=colors_pitcher[4], 
        alpha=alphaline, lw=lwidth)
    ax1.plot(shortmoment[4, :, 1], shortmoment[4, :, 0], markersel[4], 
        color=colors_pitcher[4], alpha=alphamarker, markersize=mksize_quater)
    ax1.plot(np.min(momentarray[4, 540:, 1]), momentarray[4, (int(np.where(
              momentarray[4, 540:, 1] == np.min(momentarray[4, 540:, 1]))[0])+ 540), 0],
              markersize=10, color=colors_pitcher[4], marker='x')

    ax1.plot(momentarray[5, :, 1], momentarray[5, :, 0], color=colors_pitcher[5], 
        alpha=alphaline, lw=lwidth)
    ax1.plot(shortmoment[5, :, 1], shortmoment[5, :, 0], markersel[4], 
        color=colors_pitcher[5], alpha=alphamarker, markersize=mksize_quater)
    ax1.plot(np.min(momentarray[5, 540:, 1]), momentarray[5, (int(np.where(
              momentarray[5, 540:, 1] == np.min(momentarray[5, 540:, 1]))[0])+ 540), 0],
              markersize=10, color=colors_pitcher[5], marker='x')

    ax1.plot(momentarray[6, :, 1], momentarray[6, :, 0], color=colors_pitcher[6], 
        alpha=alphaline, lw=lwidth)
    ax1.plot(shortmoment[6, :, 1], shortmoment[6, :, 0], markersel[4], 
        color=colors_pitcher[6], alpha=alphamarker, markersize=mksize_quater)
    ax1.plot(np.min(momentarray[6, 540:, 1]), momentarray[6, (int(np.where(
              momentarray[6, 540:, 1] == np.min(momentarray[6, 540:, 1]))[0])+ 540), 0],
              markersize=10, color=colors_pitcher[6], marker='x')

    ax1.plot(momentarray[7, :, 1], momentarray[7, :, 0], color=colors_pitcher[7], 
        alpha=alphaline, lw=lwidth)
    ax1.plot(shortmoment[7, :, 1], shortmoment[7, :, 0], markersel[4], 
        color=colors_pitcher[7], alpha=alphamarker, markersize=mksize_quater)
    ax1.plot(np.min(momentarray[7, 540:, 1]), momentarray[7, (int(np.where(
              momentarray[7, 540:, 1] == np.min(momentarray[7, 540:, 1]))[0])+ 540), 0],
              markersize=10, color=colors_pitcher[7], marker='x')

    ax1.plot(momentarray[8, :, 1], momentarray[8, :, 0], color=colors_pitcher[8], 
        alpha=alphaline, lw=lwidth)
    ax1.plot(shortmoment[8, :, 1], shortmoment[8, :, 0], markersel[4], 
        color=colors_pitcher[8], alpha=alphamarker, markersize=mksize_quater)
    ax1.plot(np.min(momentarray[8, 540:, 1]), momentarray[8, (int(np.where(
              momentarray[8, 540:, 1] == np.min(momentarray[8, 540:, 1]))[0])+ 540), 0],
              markersize=10, color=colors_pitcher[8], marker='x')
    #ax1.set_xlabel('I\nin mm4')
    #ax1.set_ylabel('Spine Index')
    ax1.set_xscale('log')
    ax1.set_xlim(0.01,5000)
    ax1.set_ylim(-15.5,5.5)
    ax1.yaxis.set_major_locator(MultipleLocator(2))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    # Create subplot 3 - Difference of curvatures
    ax2.axvline(0, color=col_lines, lw=0.8)
    ax2.axhline(0, color=col_lines, lw=0.8)
    ax2.plot(differencearray[:, 1], dataarray[2, :, 0], color=col_curve_up, 
        alpha=alphaline, lw=lwidth)
    ax2.plot(differencearray[:, 1], dataarray[2, :, 0], markersel[0], 
        color=col_curve_up, alpha=alphamarker, markersize=mksize)
    ax2.plot(differencearray[:, 0], dataarray[0, :, 0], color=col_curve_down, 
        alpha=alphaline, lw=lwidth)
    ax2.plot(differencearray[:, 0], dataarray[0, :, 0], markersel[2], 
        color=col_curve_down, alpha=alphamarker, markersize=mksize)
    #ax2.set_xlabel('Change of curvature\nin 1/mm')
    ax2.set_xlim(-0.04,0.04)
    ax2.set_ylim(-15.5,5.5)
    ax2.yaxis.set_major_locator(MultipleLocator(2))
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    # Create subplot - Reinforced tissue
    ax3.axhline(0, color=col_lines, lw=0.8)
    #ax3.plot(100*momentarray[0, :, 4]/momentarray[0, :, 3], momentarray[0, :, 0],
    #        color=col_density, alpha=alphaline, lw=lwidth)
    #ax3.plot(100*shortmoment[0, :, 4]/shortmoment[0, :, 3], shortmoment[0, :, 0],
    #        markersel[3], color=col_density, alpha=alphamarker, markersize=mksize)
    ax3.plot(momentarray[0, :, 4], momentarray[0, :, 0],
            color=col_I, alpha=alphaline, lw=lwidth)
    ax3.plot(shortmoment[0, :, 4], shortmoment[0, :, 0],
            markersel[3], color=col_I, alpha=alphamarker, markersize=mksize)
    # Quater section
    #ax3.plot(100*momentarray[1, :, 4]/momentarray[1, :, 3], momentarray[0, :, 0],
    #        color=col_density_quater, alpha=alphaline, lw=lwidth)
    #ax3.plot(100*shortmoment[1, :, 4]/shortmoment[1, :, 3], shortmoment[0, :, 0],
    #        markersel[4], color=col_density_quater, alpha=alphamarker, 
    #        markersize=mksize_quater)
    ax3.plot(momentarray[1, :, 4], momentarray[0, :, 0],
            color=colors_pitcher[0], alpha=alphaline, lw=lwidth)
    ax3.plot(shortmoment[1, :, 4], shortmoment[0, :, 0],
            markersel[4], color=colors_pitcher[0], alpha=alphamarker, 
            markersize=mksize_quater)

    ax3.plot(momentarray[2, :, 4], momentarray[0, :, 0],
            color=colors_pitcher[1], alpha=alphaline, lw=lwidth)
    ax3.plot(shortmoment[2, :, 4], shortmoment[0, :, 0],
            markersel[4], color=colors_pitcher[1], alpha=alphamarker, 
            markersize=mksize_quater)

    ax3.plot(momentarray[3, :, 4], momentarray[0, :, 0],
            color=colors_pitcher[2], alpha=alphaline, lw=lwidth)
    ax3.plot(shortmoment[3, :, 4], shortmoment[0, :, 0],
            markersel[4], color=colors_pitcher[2], alpha=alphamarker, 
            markersize=mksize_quater)

    ax3.plot(momentarray[4, :, 4], momentarray[0, :, 0],
            color=colors_pitcher[4], alpha=alphaline, lw=lwidth)
    ax3.plot(shortmoment[4, :, 4], shortmoment[0, :, 0],
            markersel[4], color=colors_pitcher[4], alpha=alphamarker, 
            markersize=mksize_quater)

    ax3.plot(momentarray[5, :, 4], momentarray[0, :, 0],
            color=colors_pitcher[5], alpha=alphaline, lw=lwidth)
    ax3.plot(shortmoment[5, :, 4], shortmoment[0, :, 0],
            markersel[4], color=colors_pitcher[5], alpha=alphamarker, 
            markersize=mksize_quater)

    ax3.plot(momentarray[6, :, 4], momentarray[0, :, 0],
            color=colors_pitcher[6], alpha=alphaline, lw=lwidth)
    ax3.plot(shortmoment[6, :, 4], shortmoment[0, :, 0],
            markersel[4], color=colors_pitcher[6], alpha=alphamarker, 
            markersize=mksize_quater)

    ax3.plot(momentarray[7, :, 4], momentarray[0, :, 0],
            color=colors_pitcher[7], alpha=alphaline, lw=lwidth)
    ax3.plot(shortmoment[7, :, 4], shortmoment[0, :, 0],
            markersel[4], color=colors_pitcher[7], alpha=alphamarker, 
            markersize=mksize_quater)

    ax3.plot(momentarray[8, :, 4], momentarray[0, :, 0],
            color=colors_pitcher[8], alpha=alphaline, lw=lwidth)
    ax3.plot(shortmoment[8, :, 4], shortmoment[0, :, 0],
            markersel[4], color=colors_pitcher[8], alpha=alphamarker, 
            markersize=mksize_quater)
    #ax3.set_xlabel('Proportion of\nreinforced tissue\nin %')
    ax3.set_xlim(0,25)
    #ax3.set_xlim(0,90)
    ax3.set_ylim(-15.5,5.5)
    ax3.yaxis.set_major_locator(MultipleLocator(2))
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    # Save figure
    fig.tight_layout(pad=0.5, w_pad=-2, h_pad=0)
    fig.set_size_inches(10, 7)
    figname = path / ('Mat1_' + basename + '.pdf')
    fig.savefig(figname)
    logging.info("Plot saved under: {}".format(figname))

    plt.close()


    mins = np.array([[30,momentarray[1, (int(np.where(momentarray[1, 540:, 1] == 
        np.min(momentarray[1, 540:, 1]))[0])+ 540), 0]],

        [45,momentarray[2, (int(np.where(momentarray[2, 540:, 1] == 
        np.min(momentarray[2, 540:, 1]))[0])+ 540), 0]],

        [60,momentarray[3, (int(np.where(momentarray[3, 540:, 1] == 
        np.min(momentarray[3, 540:, 1]))[0])+ 540), 0]],

        [75,momentarray[4, (int(np.where(momentarray[4, 540:, 1] == 
        np.min(momentarray[4, 540:, 1]))[0])+ 540), 0]],

        [90,momentarray[5, (int(np.where(momentarray[5, 540:, 1] == 
        np.min(momentarray[5, 540:, 1]))[0])+ 540), 0]],

        [105,momentarray[6, (int(np.where(momentarray[6, 540:, 1] == 
        np.min(momentarray[6, 540:, 1]))[0])+ 540), 0]],

        [120,momentarray[7, (int(np.where(momentarray[7, 540:, 1] == 
        np.min(momentarray[7, 540:, 1]))[0])+ 540), 0]],

        [180,momentarray[8, (int(np.where(momentarray[8, 540:, 1] == 
        np.min(momentarray[8, 540:, 1]))[0])+ 540), 0]],

        #[360,momentarray[0, (int(np.where(momentarray[0, 540:, 1] == 
        #np.min(momentarray[0, 540:, 1]))[0])+ 540), 0]]
        ])
    print(mins)

    fig1, (ax) = plt.subplots(1, 1)

    ax.plot(mins[:,0], mins[:, 1], 'bx')
    ax.set_xlim(190, 20)
    ax.set_ylim(-15, -1)
    #plt.axis([190, 20, -15, -1])
    
    ax.set_xticks(np.arange(30, 181, 15))
    ax.set_xlabel ('Part of cross-section in degree')
    ax.set_ylabel ('Index of I_min')
    ax.grid(axis = 'y')
    #plt.show()

    # Save figure
    fig1.tight_layout(pad=0.5, w_pad=-2, h_pad=0)
    fig1.set_size_inches(10, 7)
    fig1name = path / ('Mat1_Converg_' + basename + '.pdf')
    fig1.savefig(fig1name)
    logging.info("Plot saved under: {}".format(fig1name))
        
    

    return None


def fig_result0(datalist, path):

    fig, (ax1, ax2) = plt.subplots(1, 2)
    # Add subplots and formatting
    # Create subplot 1 - Second moment of area
    #ax1.title.set_text("Second moment of area\nI in mm4")
    ax1.axhline(0, color=col_lines, lw=0.8)
    ax1.set_xscale('log')
    ax1.set_xlim(0.2,7000)
    ax1.set_ylim(-15.5,5.5)
    ax1.yaxis.set_major_locator(MultipleLocator(2))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    # Create subplot 4 - Reinforced tissue
    #ax2.title.set_text("Area of\nreinforced tissue in mm2")
    ax2.axhline(0, color=col_lines, lw=0.8)
    ax2.set_xlim(-2,33)
    ax2.set_ylim(-15.5,5.5)
    ax2.yaxis.set_major_locator(MultipleLocator(2))
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    # Add data curves to each subplot
    for i, pitcher in enumerate(datalist):
        if i != 10:
            momentarray = pitcher[3]
            shortmoment = pitcher[4]

            # Subplot 1 - Second moment of area
            ax1.plot(momentarray[0, :, 1], momentarray[0, :, 0], color=colors_pitcher[i], 
                alpha=alphaline, lw=lwidth)
            ax1.plot(shortmoment[0, :, 1], shortmoment[0, :, 0], markersel[3], 
                color=colors_pitcher[i], alpha=alphamarker, markersize=mksize)
            # Subplot 4 - Reinforced tissue
            ax2.plot(momentarray[0, :, 4], momentarray[0, :, 0],
                    color=colors_pitcher[i], alpha=alphaline, lw=lwidth)
            ax2.plot(shortmoment[0, :, 4], shortmoment[0, :, 0],
                    markersel[3], color=colors_pitcher[i], alpha=alphamarker, markersize=mksize)
    
    # Save figure
    fig.tight_layout(pad=0.5, w_pad=0, h_pad=0)
    fig.set_size_inches(5, 7)
    figname = path / ('Res0_all_pitchers.pdf')
    fig.savefig(figname)
    logging.info("Plot saved under: {}".format(figname))

    plt.close()

    return None

def fig_result1(basename, dataarray, differencearray, momentarray, shortmoment, path):


    cut = [185, 201, 353, 430, 383, 435, 409]
    cut = cut[int(basename[-1])-1]

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

    # Create subplot 5 - Second moment of area and area
    ax1.axhline(0, color=col_lines, lw=0.8)
    ax1.plot(momentarray[0, :, 1], momentarray[0, :, 0], color=col_I, 
        alpha=alphaline, lw=lwidth)
    ax1.plot(shortmoment[0, :, 1], shortmoment[0, :, 0], markersel[3], 
        color=col_I, alpha=alphamarker, markersize=mksize)
    ax1.plot(np.min(momentarray[0, cut:, 1]), momentarray[0, (int(np.where(
        momentarray[0, cut:, 1] == np.min(momentarray[0, cut:, 1]))[0])+ cut), 0],
        markersize=10, color=col_I, marker='x')
    ax1.set_xscale('log')
    ax1.set_xlim(0.2,7000)
    ax1.set_ylim(-15.1,5.1)
    #ax1.yaxis.set_major_locator(MultipleLocator(2))
    #ax1.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax1.set_yticks([-14, -12, -10, -8, -6, -4, -2, 0, 2, 4])
    ax1.set_xticks([1000, 10, 0.1])
    ax1.set_yticklabels([])
    ax1.set_xticklabels([])

    # Quater section
    ax4 = ax1.twiny()
    ax4.plot(momentarray[1, :, 1], momentarray[1, :, 0], color=col_I_quater, 
        alpha=alphaline, lw=lwidth)
    ax4.plot(shortmoment[1, :, 1], shortmoment[1, :, 0], markersel[4], 
        color=col_I_quater, alpha=alphamarker, markersize=mksize_quater)
    ax4.plot(np.min(momentarray[1, cut:, 1]), momentarray[1, (int(np.where(
        momentarray[1, cut:, 1] == np.min(momentarray[1, cut:, 1]))[0])+ cut), 0],
        markersize=10, color=col_I_quater, marker='x')
    #ax1.set_xlabel('I\nin mm4')
    #ax1.set_ylabel('Spine Index')
    ax4.set_xscale('log')
    ax4.set_xlim(0.01,20)
    ax4.set_xticks([10, 1, 0.1, 0.01])
    ax4.set_xticklabels([])
    


    # Create subplot 3 - Difference of curvatures
    ax2.axvline(0, color=col_lines, lw=0.8)
    ax2.axhline(0, color=col_lines, lw=0.8)
    ax2.plot(differencearray[:, 1], dataarray[2, :, 0], color=col_curve_up, 
        alpha=alphaline, lw=lwidth)
    ax2.plot(differencearray[:, 1], dataarray[2, :, 0], markersel[0], 
        color=col_curve_up, alpha=alphamarker, markersize=mksize)
    ax2.plot(differencearray[:, 0], dataarray[0, :, 0], color=col_curve_down, 
        alpha=alphaline, lw=lwidth)
    ax2.plot(differencearray[:, 0], dataarray[0, :, 0], markersel[2], 
        color=col_curve_down, alpha=alphamarker, markersize=mksize)
    #ax2.set_xlabel('Change of curvature\nin 1/mm')
    ax2.set_xlim(-0.05,0.05)
    ax2.set_ylim(-15.1,5.1)
    #ax2.yaxis.set_major_locator(MultipleLocator(2))
    #ax2.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax2.set_yticks([-14, -12, -10, -8, -6, -4, -2, 0, 2, 4])
    ax2.set_xticks([-0.03, 0, 0.03])
    ax2.set_yticklabels([])
    ax2.set_xticklabels([])

    # Create subplot - Reinforced tissue
    ax3.axhline(0, color=col_lines, lw=0.8)
    ax3.plot(100*momentarray[0, :, 4]/momentarray[0, :, 3], momentarray[0, :, 0],
            color=col_density, alpha=alphaline, lw=lwidth)
    ax3.plot(100*shortmoment[0, :, 4]/shortmoment[0, :, 3], shortmoment[0, :, 0],
            markersel[3], color=col_density, alpha=alphamarker, markersize=mksize)
    #ax3.plot(momentarray[0, :, 4], momentarray[0, :, 0],
    #        color=col_density, alpha=alphaline, lw=lwidth)
    #ax3.plot(shortmoment[0, :, 4], shortmoment[0, :, 0],
    #        markersel[3], color=col_density, alpha=alphamarker, markersize=mksize)
    # Quater section
    ax3.plot(100*momentarray[1, :, 4]/momentarray[1, :, 3], momentarray[0, :, 0],
            color=col_density_quater, alpha=alphaline, lw=lwidth)
    ax3.plot(100*shortmoment[1, :, 4]/shortmoment[1, :, 3], shortmoment[0, :, 0],
            markersel[4], color=col_density_quater, alpha=alphamarker, 
            markersize=mksize_quater)
    #ax3.plot(momentarray[1, :, 4], momentarray[0, :, 0],
     #       color=col_density_quater, alpha=alphaline, lw=lwidth)
    #ax3.plot(shortmoment[1, :, 4], shortmoment[0, :, 0],
     #       markersel[4], color=col_density_quater, alpha=alphamarker, 
      #      markersize=mksize_quater)
    #ax3.set_xlabel('Proportion of\nreinforced tissue\nin %')
    #ax3.set_xlim(-2,33)
    ax3.set_xlim(0,100)
    ax3.set_ylim(-15.1,5.1)
    #ax3.yaxis.set_major_locator(MultipleLocator(2))
    #ax3.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax3.set_yticks([-14, -12, -10, -8, -6, -4, -2, 0, 2, 4])
    ax3.set_xticks([0, 25, 50, 75, 100])
    ax3.set_yticklabels([])
    ax3.set_xticklabels([])

    # Save figure
    fig.tight_layout(pad=0, w_pad=0.5, h_pad=0)
    fig.set_size_inches(4, 5)
    figname = path / ('Res1_' + basename + '.pdf')
    fig.savefig(figname)
    logging.info("Plot saved under: {}".format(figname))

    plt.close()

    return None

def plot_mean_cis(datalist, path):
    # Variable to save data in
    sorteddata = np.zeros((7,22,10))
    # Read out data for each pitcher
    for i,pitcher in enumerate(datalist):
        # Extract different arrays
        basename = pitcher[0]
        differencearray = np.vstack((np.flipud(pitcher[2]), np.zeros((2,2))))
        differencearray = np.vstack((differencearray[:6,:], differencearray[5:,:]))
        shortmoment = pitcher[4]

        # Accumulate data per pitcher in single array
        pitcherarray = np.vstack((
        # Pitcher no
        np.ones((len(shortmoment[0])))*int(basename[-1]),
        # Spine Index
        shortmoment[0, :, 0],
        # Second moment of area
        shortmoment[0, :, 1],
        # Local second moment of area
        shortmoment[1, :, 1],
        # Curvature change down
        differencearray[:len(shortmoment[0]), 0],
        # Curvature change up
        differencearray[:len(shortmoment[0]), 1],
        # Area of reinforced tissue
        shortmoment[0, :, 4],
        # Local area of reinforced tissue
        shortmoment[1, :, 4],
        # Reinforced tissue, proportional
        100*shortmoment[0, :, 4]/shortmoment[0, :, 3],
        # Local reinforced tissue, proportional
        100*shortmoment[1, :, 4]/shortmoment[1, :, 3])).T

        # Add to overall data
        nans = np.zeros((22-len(pitcherarray),10))
        nans[:] = np.nan
        sorteddata[i,:,:] = np.vstack((pitcherarray, nans))

    means = np.zeros((22,10))
    means = np.nanmean(sorteddata, axis=0)
    cislow, cishigh = np.zeros((22,10)),np.zeros((22,10))
    cislow[:,:2], cishigh[:,:2] = means[:, :2], means[:, :2]
    cislow[:,2:], cishigh[:,2:] = norm.interval(0.95, loc=np.nanmean(sorteddata[:,:,2:], axis=0), 
                   scale=stat.sem(sorteddata[:,:,2:], axis=0, nan_policy='omit'))
    
    #Todo: Plot all!
    alphacis = 0.4
    fig, (ax1, ax2, ax3) = plt.subplots(3,1)
    index = means[:,1]

    ax4 = ax1.twinx()
    ax1.axvline(0, color=col_lines, lw=0.8)
    ax1.fill_between(index, cislow[:,2], cishigh[:,2], color=col_I, alpha=alphacis, lw=0)
    ax4.fill_between(index, cislow[:,3], cishigh[:,3], color=col_I_quater, alpha=alphacis, lw=0)
    ax1.plot(index, means[:,2], color=col_I, lw=lwidth)
    ax4.plot(index, means[:,3], color=col_I_quater, lw=lwidth)
    ax1.set_yscale('log')
    ax1.set_ylim(7000, 0.2)
    ax1.set_xlim(-15.1,5.1)
    ax1.xaxis.tick_top()
    ax1.set_xticks([-14, -12, -10, -8, -6, -4, -2, 0, 2, 4])
    ax1.set_yticks([1000, 10, 0.1])
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    #ax1.set_xticklabels(ax1.get_xticks(), rotation = -90)
    #ax1.set_yticklabels(ax1.get_yticks(), rotation = -90)
    #ax1.xaxis.set_major_locator(MultipleLocator(2))
    #ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))

    ax4.set_yscale('log')
    ax4.set_ylim(20, 0.01)
    ax4.set_yticks([10, 1, 0.1, 0.01])
    ax4.set_yticklabels([])
    #ax4.set_yticklabels(ax4.get_yticks(), rotation = -90)

    ax2.axvline(0, color=col_lines, lw=0.8)
    ax2.axhline(0, color=col_lines, lw=0.8)
    ax2.fill_between(index, cislow[:,5], cishigh[:,5], color=col_curve_up, alpha=alphacis, lw=0)
    ax2.fill_between(index, cislow[:,4], cishigh[:,4], color=col_curve_down, alpha=alphacis, lw=0)
    ax2.plot(index, means[:,5], color=col_curve_up, lw=lwidth)
    ax2.plot(index, means[:,4], color=col_curve_down, lw=lwidth)
    ax2.set_ylim(0.05,-0.05)
    ax2.set_xlim(-15.1,5.1)
    #ax2.set_xticklabels(ax1.get_xticks(), rotation = -90)
    #ax2.xaxis.set_major_locator(MultipleLocator(2))
    #ax2.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    #ax2.set_xticks([])
    ax2.xaxis.tick_top()
    ax2.set_xticks([-14, -12, -10, -8, -6, -4, -2, 0, 2, 4])
    ax2.set_yticks([-0.03, 0, 0.03])
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    #ax2.yaxis.set_major_formatter(FormatStrFormatter('%3d'))
    #ax2.set_yticklabels(ax2.get_yticks(), rotation = -90)




    ax3.axvline(0, color=col_lines, lw=0.8)
    ax3.fill_between(index, cislow[:,8], cishigh[:,8], color=col_density, alpha=alphacis, lw=0)
    ax3.fill_between(index, cislow[:,9], cishigh[:,9], color=col_density_quater, alpha=alphacis, lw=0)
    ax3.plot(index, means[:,8], color=col_density, lw=lwidth)
    ax3.plot(index, means[:,9], color=col_density_quater, lw=lwidth)
    ax3.set_ylim(100,0)
    ax3.set_xlim(-15.1,5.1)
    #ax3.xaxis.set_major_locator(MultipleLocator(2))
    #ax3.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax3.xaxis.tick_top()
    #ax3.set_xticklabels(ax1.get_xticks(), rotation = -90)
    ax3.set_xticks([-14, -12, -10, -8, -6, -4, -2, 0, 2, 4])
    ax3.set_yticks([0, 25, 50, 75, 100])
    ax3.set_xticklabels([])
    ax3.set_yticklabels([])

    # Save figure
    fig.tight_layout(pad=0, h_pad=0.5, w_pad=0)
    fig.set_size_inches(5,4)
    figname = path / ('Res2_Mean_CIs.pdf')
    fig.savefig(figname)
    logging.info("Plot saved under: {}".format(figname))

    plt.close()
    return

def fig_result2(datalist, path):

    #fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(1, 7)
    #fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5)
    fig, (ax1, ax2, ax3, ax6, ax7) = plt.subplots(1, 5)
    
    # Add subplots and formatting
    # Create subplot 1 - Second moment of area
    ax1.title.set_text("Second moment of area\nI in mm4")
    ax1.axhline(0, color=col_lines, lw=0.8)
    ax1.set_xscale('log')
    ax1.set_xlim(0.2,7000)
    ax1.set_ylim(-15.5,5.5)
    ax1.yaxis.set_major_locator(MultipleLocator(2))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    # Create subplot 2 - Local second moment of area
    ax2.title.set_text("Local second moment of area\nI in mm4")
    ax2.axhline(0, color=col_lines, lw=0.8)
    ax2.set_xscale('log')
    ax2.set_xlim(0.01,20)
    ax2.set_ylim(-15.5,5.5)
    ax2.yaxis.set_major_locator(MultipleLocator(2))
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    # Create subplot 3 - Difference of curvatures
    ax3.title.set_text("Curvature change down\ndK in 1/mm")
    ax3.axvline(0, color=col_lines, lw=0.8)
    ax3.axhline(0, color=col_lines, lw=0.8)
    ax3.set_xlim(-0.06,0.05)
    ax3.set_ylim(-15.5,5.5)
    ax3.yaxis.set_major_locator(MultipleLocator(2))
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    # Create subplot 4 - Reinforced tissue
    # ax4.title.set_text("Area of\nreinforced tissue in mm2")
    # ax4.axhline(0, color=col_lines, lw=0.8)
    # ax4.set_xlim(-2,33)
    # ax4.set_ylim(-15.5,5.5)
    # ax4.yaxis.set_major_locator(MultipleLocator(2))
    # ax4.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    # # Create subplot 5 - Local reinforced tissue
    # ax5.title.set_text("Local are of\nreinforced tissue in mm2")
    # ax5.axhline(0, color=col_lines, lw=0.8)
    # ax5.set_xlim(-1,8)
    # ax5.set_ylim(-15.5,5.5)
    # ax5.yaxis.set_major_locator(MultipleLocator(2))
    # ax5.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    # Create subplot 6 - Reinforced tissue
    ax6.title.set_text("Proportion of\nreinforced tissue in %")
    ax6.axhline(0, color=col_lines, lw=0.8)
    ax6.set_xlim(0,90)
    ax6.set_ylim(-15.5,5.5)
    ax6.yaxis.set_major_locator(MultipleLocator(2))
    ax6.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    # Create subplot 7 - Local reinforced tissue
    ax7.title.set_text("Local proportion of\nreinforced tissue in %")
    ax7.axhline(0, color=col_lines, lw=0.8)
    ax7.set_xlim(0,90)
    ax7.set_ylim(-15.5,5.5)
    ax7.yaxis.set_major_locator(MultipleLocator(2))
    ax7.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    # Add data curves to each subplot
    for i, pitcher in enumerate(datalist):
        if i != 10:
            basename = pitcher[0]
            dataarray = pitcher[1]
            differencearray = pitcher[2]
            momentarray = pitcher[3]
            shortmoment = pitcher[4]

            cut = [185, 201, 353, 430, 383, 435, 409]
            cut = cut[i]

            # Subplot 1 - Second moment of area
            ax1.plot(momentarray[0, :, 1], momentarray[0, :, 0], color=colors_pitcher[i], 
                alpha=alphaline, lw=lwidth)
            ax1.plot(shortmoment[0, :, 1], shortmoment[0, :, 0], markersel[3], 
                color=colors_pitcher[i], alpha=alphamarker, markersize=mksize)
            ax1.plot(np.min(momentarray[0, cut:, 1]), momentarray[0, (int(np.where(
                momentarray[0, cut:, 1] == np.min(momentarray[0, cut:, 1]))[0])+ cut), 0],
                markersize=10, color=colors_pitcher[i], marker='x')
            # Subplot 2 - Local second moment of area
            ax2.plot(momentarray[1, :, 1], momentarray[1, :, 0], color=colors_pitcher[i], 
                alpha=alphaline, lw=lwidth)
            ax2.plot(shortmoment[1, :, 1], shortmoment[1, :, 0], markersel[4], 
                color=colors_pitcher[i], alpha=alphamarker, markersize=mksize_quater)
            ax2.plot(np.min(momentarray[1, cut:, 1]), momentarray[1, (int(np.where(
                momentarray[1, cut:, 1] == np.min(momentarray[1, cut:, 1]))[0])+ cut), 0],
                markersize=10, color=colors_pitcher[i], marker='x')
            # Subplot 3 - Difference of curvatures
            #ax3.plot(differencearray[:, 1], dataarray[2, :, 0], color=colors_pitcher[i], 
            #    alpha=alphaline, lw=lwidth)
            ax3.plot(differencearray[:, 0], dataarray[0, :, 0], color=colors_pitcher[i], 
                alpha=alphaline, lw=lwidth)
            #ax3.plot(differencearray[:, 1], dataarray[2, :, 0], markersel[0], 
            #    color=colors_pitcher[i], alpha=alphamarker, markersize=mksize)
            ax3.plot(differencearray[:, 0], dataarray[0, :, 0], markersel[2], 
                color=colors_pitcher[i], alpha=alphamarker, markersize=mksize)
            # Subplot 4 - Reinforced tissue
            # ax4.plot(momentarray[0, :, 4], momentarray[0, :, 0],
            #         color=colors_pitcher[i], alpha=alphaline, lw=lwidth)
            # ax4.plot(shortmoment[0, :, 4], shortmoment[0, :, 0],
            #         markersel[3], color=colors_pitcher[i], alpha=alphamarker, markersize=mksize)
            # # Subplot 5 - Local reinforced tissue
            # ax5.plot(momentarray[1, :, 4], momentarray[0, :, 0],
            #         color=colors_pitcher[i], alpha=alphaline, lw=lwidth)
            # ax5.plot(shortmoment[1, :, 4], shortmoment[0, :, 0],
            #         markersel[4], color=colors_pitcher[i], alpha=alphamarker, 
            #         markersize=mksize_quater)
            # Subplot 6 - Reinforced tissue, proportional
            ax6.plot(100*momentarray[0, :, 4]/momentarray[0, :, 3], momentarray[0, :, 0],
                    color=colors_pitcher[i], alpha=alphaline, lw=lwidth)
            ax6.plot(100*shortmoment[0, :, 4]/shortmoment[0, :, 3], shortmoment[0, :, 0],
                    markersel[3], color=colors_pitcher[i], alpha=alphamarker, markersize=mksize)
            # Subplot 7 - Local reinforced tissue, proportional
            ax7.plot(100*momentarray[1, :, 4]/momentarray[1, :, 3], momentarray[0, :, 0],
                    color=colors_pitcher[i], alpha=alphaline, lw=lwidth)
            ax7.plot(100*shortmoment[1, :, 4]/shortmoment[1, :, 3], shortmoment[0, :, 0],
                    markersel[4], color=colors_pitcher[i], alpha=alphamarker, 
                    markersize=mksize_quater)

    # Save figure
    fig.tight_layout(pad=0.5, w_pad=-2, h_pad=0)
    #fig.set_size_inches(18, 7)
    fig.set_size_inches(13, 7)
    figname = path / ('Res2_all_pitchers.pdf')
    fig.savefig(figname)
    logging.info("Plot saved under: {}".format(figname))

    plt.close()

    return None

def statistic(x, y):  # permute only `x`
    return stat.spearmanr(x, y)[0]

def do_permutation(x, y):
    x = x.tolist()
    y = y.tolist()

    statistic_function = lambda xlist: statistic(xlist, y)
    
    res_exact = stat.permutation_test((x,), statistic_function, permutation_type='pairings')
    return res_exact.pvalue

def fig_result3(basename, dataarray, dataarray2, differencearray, differencearray2, path):

    # Create main figure   
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, gridspec_kw={
        'width_ratios': [2.2, 1, 1, 1]})

    # Create subplot 1 - Curves along pitcher midrib
    ax1.axis('equal')
    ax1.axhline(dataarray[1, -6, 2], 0.5, color=col_lines, lw=0.8)
    ax1.axhline(dataarray2[1, -6, 2], 0, 0.65, color=col_lines, lw=0.8)
    ax1.plot(dataarray[2, :, 1], dataarray[2, :, 2], color=colors_p7[3],
        lw=lwidth)
    ax1.plot(dataarray[1, :, 1], dataarray[1, :, 2], color=colors_p7[4],
        lw=lwidth)
    ax1.plot(dataarray[0, :, 1], dataarray[0, :, 2], color=colors_p7[5],
        lw=lwidth)
    ax1.plot(dataarray2[2, :, 1], dataarray2[2, :, 2], color=colors_p7[0],
        lw=lwidth)
    ax1.plot(dataarray2[1, :, 1], dataarray2[1, :, 2], color=colors_p7[1],
        lw=lwidth)
    ax1.plot(dataarray2[0, :, 1], dataarray2[0, :, 2], color=colors_p7[2],
        lw=lwidth)
    #ax1.set_xlabel('X Coordinates\n in mm')
    #ax1.set_ylabel('Z Coordinates in mm')
    ax1.set_ylim(0,80)
    # Annotation: indices to middle curve for orientation    
    indices = dataarray[1, :, :]
    indices = np.flip(indices)[1:-1:2, :]
    indices = np.flip(indices)
    ax1.plot(indices[:, 1], indices[:, 2], 'x', color=col_annotation, markersize=5)
    for i in range(len(indices)):
        ax1.annotate(int(indices[i, 0]), (indices[i, 1] + 1, indices[i, 2]),
                     va="center", ha="left", size=10, color=col_annotation[4])
    indices2 = dataarray2[1, :, :]
    indices2 = np.flip(indices2)[1:-1:2, :]
    indices2 = np.flip(indices2)
    ax1.plot(indices2[:, 1], indices2[:, 2], 'x', color=col_annotation, markersize=5)
    for i in range(len(indices2)):
        ax1.annotate(int(indices2[i, 0]), (indices2[i, 1] + 1, indices2[i, 2]),
                     va="center", ha="left", size=10, color=col_annotation[1])

    # Create subplot 2 - Curvature for each point on curves
    ax2.axvline(0, color=col_lines, lw=0.8)
    ax2.axhline(0, color=col_lines, lw=0.8)
    ax2.plot(dataarray2[1, :, 3], dataarray2[1, :, 0], color=colors_p7[1], 
        alpha=alphaline, lw=lwidth)
    ax2.plot(dataarray2[1, :, 3], dataarray2[1, :, 0], markersel[1], 
        color=colors_p7[1], alpha=alphamarker, markersize=mksize)
    ax2.plot(dataarray2[2, :, 3], dataarray2[2, :, 0], color=colors_p7[0], 
        alpha=alphaline, lw=lwidth)
    ax2.plot(dataarray2[2, :, 3], dataarray2[2, :, 0], markersel[0], 
        color=colors_p7[0], alpha=alphamarker, markersize=mksize)
    ax2.plot(dataarray2[0, :, 3], dataarray2[0, :, 0], color=colors_p7[2], 
        alpha=alphaline, lw=lwidth)
    ax2.plot(dataarray2[0, :, 3], dataarray2[0, :, 0], markersel[2], 
        color=colors_p7[2], alpha=alphamarker, markersize=mksize)
    #ax2.set_xlabel('Curvature\nin 1/mm')
    #ax2.set_ylabel('Spine Index')
    ax2.set_ylim(-15.5,5.5)
    ax2.set_xlim(-0.2,0.05)
    ax2.yaxis.set_major_locator(MultipleLocator(2))
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    ax3.axvline(0, color=col_lines, lw=0.8)
    ax3.axhline(0, color=col_lines, lw=0.8)
    ax3.plot(dataarray[1, :, 3], dataarray[1, :, 0], color=colors_p7[4], 
        alpha=alphaline, lw=lwidth)
    ax3.plot(dataarray[1, :, 3], dataarray[1, :, 0], markersel[1], 
        color=colors_p7[4], alpha=alphamarker, markersize=mksize)
    ax3.plot(dataarray[2, :, 3], dataarray[2, :, 0], color=colors_p7[3], 
        alpha=alphaline, lw=lwidth)
    ax3.plot(dataarray[2, :, 3], dataarray[2, :, 0], markersel[0], 
        color=colors_p7[3], alpha=alphamarker, markersize=mksize)
    ax3.plot(dataarray[0, :, 3], dataarray[0, :, 0], color=colors_p7[5], 
        alpha=alphaline, lw=lwidth)
    ax3.plot(dataarray[0, :, 3], dataarray[0, :, 0], markersel[2], 
        color=colors_p7[5], alpha=alphamarker, markersize=mksize)
    #ax3.set_xlabel('Curvature\nin 1/mm')
    ax3.set_ylim(-15.5,5.5)
    ax3.set_xlim(-0.2,0.05)
    ax3.yaxis.set_major_locator(MultipleLocator(2))
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    # Create subplot 3 - Difference of curvatures
    ax4.axvline(0, color=col_lines, lw=0.8)
    ax4.axhline(0, color=col_lines, lw=0.8)
    ax4.plot(differencearray[:, 1], dataarray[2, :, 0], color=colors_p7[3], 
        alpha=alphaline, lw=lwidth)
    ax4.plot(differencearray[:, 1], dataarray[2, :, 0], markersel[0], 
        color=colors_p7[3], alpha=alphamarker, markersize=mksize)
    ax4.plot(differencearray[:, 0], dataarray[0, :, 0], color=colors_p7[5], 
        alpha=alphaline, lw=lwidth)
    ax4.plot(differencearray[:, 0], dataarray[0, :, 0], markersel[2], 
        color=colors_p7[5], alpha=alphamarker, markersize=mksize)

    ax4.plot(differencearray2[:, 1], dataarray2[2, :, 0], color=colors_p7[0], 
        alpha=alphaline, lw=lwidth)
    ax4.plot(differencearray2[:, 1], dataarray2[2, :, 0], markersel[0], 
        color=colors_p7[0], alpha=alphamarker, markersize=mksize)
    ax4.plot(differencearray2[:, 0], dataarray2[0, :, 0], color=colors_p7[2], 
        alpha=alphaline, lw=lwidth)
    ax4.plot(differencearray2[:, 0], dataarray2[0, :, 0], markersel[2], 
        color=colors_p7[2], alpha=alphamarker, markersize=mksize)
    #ax4.set_xlabel('Change of curvature\nin 1/mm')
    ax4.set_xlim(-0.03,0.03)
    ax4.set_ylim(-15.5,5.5)
    ax4.yaxis.set_major_locator(MultipleLocator(2))
    ax4.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    # Save figure
    fig.tight_layout(pad=0.5, w_pad=-2, h_pad=0)
    fig.set_size_inches(11, 7)
    figname = path / ('Res3_' + basename + '.pdf')
    fig.savefig(figname)
    logging.info("Plot saved under: {}".format(figname))

    plt.close()

    # Create correlation figure   
    fig2, (ax) = plt.subplots(1, 1)
    ax.plot(differencearray[:, 1], differencearray2[:, 1], markersel[0],
        color=col_curve_up)
    ax.plot(differencearray[:, 0], differencearray2[:, 0], markersel[2],
        color=col_curve_down)
    # Spearman correlation test for correlating 3d print and plant curvature changes
    # https://docs.scipy.org/doc/scipy-1.8.0/reference/generated/scipy.stats.spearmanr.html
    cor1 = stat.spearmanr(differencearray[:, 1], differencearray2[:, 1])
    cor2 = stat.spearmanr(differencearray[:, 0], differencearray2[:, 0])
    # Permutation test for reliable p val
    # https://docs.scipy.org/doc/scipy-1.8.0/reference/generated/scipy.stats.permutation_test.html
    pval_up = do_permutation(differencearray[:, 1], differencearray2[:, 1])
    pval_down = do_permutation(differencearray[:, 0], differencearray2[:, 0])
    ax.title.set_text('up: rho = {0}, p = {1}, p_permutation = {2}\ndown: rho = {3}, p = {4}, p_permutation = {5}'.format(
            round(cor1[0],6), round(cor1[1],6), pval_up, round(cor2[0],6), round(cor2[1],6), pval_down))


    fig2.tight_layout(pad=0.5, w_pad=-2, h_pad=0)
    fig2.set_size_inches(7, 7)
    figname = path / ('Res3_' + basename + '_correlation.pdf')
    fig2.savefig(figname)
    logging.info("Plot saved under: {}".format(figname))

    plt.close()

    return None

def export_data(datalist, path):
    # Variable to save data in
    sorteddata = np.zeros((1,10))
    # Read out data for each pitcher
    for pitcher in datalist:
        # Extract different arrays
        basename = pitcher[0]
        differencearray = np.vstack((np.flipud(pitcher[2]), np.zeros((2,2))))
        differencearray = np.vstack((differencearray[:6,:], differencearray[5:,:]))
        shortmoment = pitcher[4]

        # Accumulate data per pitcher in single array
        pitcherarray = np.vstack((
        # Pitcher no
        np.ones((len(shortmoment[0])))*int(basename[-1]),
        # Spine Index
        shortmoment[0, :, 0],
        # Second moment of area
        shortmoment[0, :, 1],
        # Local second moment of area
        shortmoment[1, :, 1],
        # Curvature change down
        differencearray[:len(shortmoment[0]), 0],
        # Area of reinforced tissue
        shortmoment[0, :, 4],
        # Local area of reinforced tissue
        shortmoment[1, :, 4],
        # Reinforced tissue, proportional
        100*shortmoment[0, :, 4]/shortmoment[0, :, 3],
        # Local reinforced tissue, proportional
        100*shortmoment[1, :, 4]/shortmoment[1, :, 3],
        # Curvature change up
        differencearray[:len(shortmoment[0]), 1])).T

        # Add to overall data
        sorteddata = np.vstack((sorteddata, pitcherarray))

    # Save data with header as csv
    header = "pitcher number, spine index, second moment of area I in mm4, local second moment of area I in mm4, curvature change down delta K in mm-1, area of reinforced tissue in mm2, local area of reinforced tissue in mm2, proportion of reinforced tissue in %, local proportion of reinforced tissue in %, curvature change up delta K in mm-1"
    filepath = path / "Data_Summary.csv"
    np.savetxt(filepath, sorteddata[1:, :], delimiter=',', header=header)

    return
