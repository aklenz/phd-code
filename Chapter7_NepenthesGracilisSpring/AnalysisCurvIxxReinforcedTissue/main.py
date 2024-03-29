""" Data analysis for Biology Letters Paper and expanded for Chapter 7

Calculate curvature for all coordinates
Calculate second moment of area for tiff stacks
Plot curves and curvatures for each sample

Logging into pyNepenthes.py

This script requires the following packages:
    * pathlib
    * numpy
    * logging
    * all modules within utils
"""

from utils import curvature, plotcurves, secondmoment
from pathlib import Path
import numpy as np
import logging


def all_files_calc_curvature(pathData, pathExport):
    """
    Load all files in pathData, calculate curvatures and export them to
    pathExport

    Comment: Only works if all files are in same folder and no subfolders
    """

    # Check, if there is raw data to analyse
    if pathData.is_dir():
        filesRawData = [x.name[:-5] for x in pathData.rglob('*.fcsv')]
        logging.info("Raw data for curvature calculation: {} files.".format(
                     len(filesRawData)))
    else:
        logging.warning("No files for calculating curvature found!")
        return

    # Check, if there is data already exported
    if pathExport.is_dir():
        filesExported = [x.name[:-14] for x in pathExport.rglob('*.csv') if not
                         x == (pathExport / 'Summary_curve_differences.csv')]
    else:
        pathExport.mkdir()
        filesExported = []
    logging.info("Exported data from curvature calculation: {} files.".format(
                 len(filesExported)))

    # List files that need to be analysed
    filesToAnalyse = [file for file in filesRawData
                      if file not in filesExported]
    logging.info("Files to be analysed: {}".format(len(filesToAnalyse)))

    # Go through the list and calculate curvature for data
    for file in filesToAnalyse:
        filename = file + ".fcsv"
        logging.info("Calculate curvature for file {}.".format(filename))
        coords = curvature.load_file(filename, pathData)
        curv = curvature.calc_curvature(coords)
        curvature.export_curvature(file, coords, curv, pathExport)
        logging.info("...Finished!")

    return None


def all_samples_calc_I(pathData, thresh, thresh_venation, pathExport=None):
    """
    Load all folder from pathData and calculate the second moment of area and
    area for all tiffs within each folder using the given threshold and save
    the calculated variables in .csv file, also save every 100th image as
    bw image and the histogram of the original image with the threshold.
    """

    # Check if data path exists
    if not pathData.is_dir():
        logging.warning("No data found to analyse!")
        return

    # Load all sample folder names but the exported folder
    # (e.g. N_gracilis_p1 from N_gracilis_p1_middle[tif])
    folder = [x.name[:-12] for x in pathData.iterdir() if x.is_dir()
              and not x == (pathData / 'Exported')]
    logging.info("Raw data for second moment: {} folder.".format(len(folder)))

    # List all exported names from csv files from export folder
    # (e.g. N_gracilis_p1 from N_gracilis_p1_middle_secondmoment.csv)
    if not pathExport:
        pathExport = pathData / "Exported"
        if not pathExport.is_dir():
            pathExport.mkdir()
    exported = [x.name[:-17] for x in pathExport.rglob('*_secondmoment.csv')]
    logging.info("Exported data for second moment: {} files.".format(
                 len(exported)))

    # List folders that need to be analysed
    folder = [x for x in folder if x not in exported]
    logging.info("Folders to be analysed: {} folder.".format(len(folder)))

    # Iterate over each sample to calculate and save the second moment of area
    for f in sorted(folder):
        logging.info("Calc I for folder {}".format(f))

        try:
            pathRes = pathData / (f + "_middle.xtekct")
            resolution, unit = secondmoment.load_res(pathRes)
        except IOError as e:
            resolution = 0.030
            unit = "mm"
            logging.warning("Resolution file not found, default values taken.")
            logging.warning(e)

        try:
            pathImport = pathData / (f + "_middle[tif]")
            logging.info("Import folder name: {}".format(pathImport))
            secondmoment.calc_I_stack(pathImport, thresh, thresh_venation,
                                      resolution, unit, pathExport)
        except IOError as e:
            logging.error("Folder or files not found, no analysis done.")
            logging.error(e)
            continue

        logging.info("...Finished!")

    return None


def all_files_plot_curves_export(pathExportC, pathExportI):
    """
    Load curvature and second moment of area data and plot all together.
    """

    # If export directory doesn't exists log and return
    if not pathExportC.is_dir() and not pathExportI.is_dir():
        logging.warning("No data found to plot!")
        return

    # Save all basenames in a list and take only unique names
    files = [x.name[:8] for x in pathExportC.rglob('*.csv') 
             if not x == (pathExportC / "Summary_curve_differences.csv")]
    basenames = sorted(set(files))
    logging.info("Data for plotting from: {}".format(basenames))

    # Make variable for all difference arrays and the regarding headers
    summary_dcurvatures = []
    header = ""
    cordata = np.zeros((len(basenames), 21, 7))
    datalist = []

    # Loop through all basenames to plot single samples
    for name in basenames:
        logging.info("Plot data for pitcher {}".format(name))
        # Load data of all 3 files with same basename in array
        try:
            data = plotcurves.load_files(name, pathExportC)
        except IOError as e:
            logging.error("Files could not be found.")
            logging.error(e)
            continue

        # Calculate the curvature differences from the 3 loaded files
        try:
            differences = plotcurves.calc_differences(data)
            # Add information to each variable
            summary_dcurvatures.extend(differences.T.tolist())
            header = header + name + " Middle-Down, " + name + " Up-Middle,"
        except TypeError as e:
            logging.error("Curvature difference could not be found.")
            logging.error(e)
            continue

        # Load second moment
        try:
            moment = plotcurves.load_moment(name, pathExportI)
            moment, shortmoment = plotcurves.adjust_momentarray(name, moment)
        except IOError as e:
            logging.warning("Second moment file could not be found.")
            logging.error(e)
            moment, shortmoment = np.zeros((1,4)), np.zeros((1,4))

        # Figure Material 0 for each pitcher
        plotcurves.fig_material0(name, data, differences, pathExportC)

        # Convergence study for p3
        if name == "Curve_p3":
            convergmoment = plotcurves.load_conv_moment(name, pathExportI)
            convergmoment, shortconvmoment = plotcurves.adjust_momentarray(name, 
                convergmoment)
            plotcurves.fig_material1(name, data, differences, convergmoment, 
                shortconvmoment, pathExportC)

        # Figure Result 1 for each pitcher
        plotcurves.fig_result1(name, data, differences, moment, shortmoment, 
            pathExportC)

        # Data preparation for correlations
        # True for 1/I and 1/dens,              False for I and dens
        #cordata[int(name[-1])-1] = plotcurves.prep_cor_data(differences,
        #                                                    shortmoment, True)

        # Correlation data for single pitchers
        # "spear" for Spearman correlation,     "lin" for linear regression
        # True for indices max dK to -15,       False for 0 to -15
        #plotcurves.single_cor(name[-1], cordata[int(name[-1])-1], pathExportC,
        #                      "spear", False)

        # Add data to list for comparison plot
        datalist.append([name, data, differences, moment, shortmoment])

        #Figure Result 3 for pitcher 7 and 3D print
        if name == "Curve_p7":
            name2 = "Curve_p7_3DPrint"
            # Load data of all 3 files with same basename in array
            try:
                data2 = plotcurves.load_files(name2, pathExportC)
            except IOError as e:
                logging.error("3D print files could not be found.")
                logging.error(e)
                continue

            # Calculate the curvature differences from the 3 loaded files
            try:
                differences2 = plotcurves.calc_differences(data2)
                # Add information to each variable
                summary_dcurvatures.extend(differences2.T.tolist())
                header = header + name2 + " Middle-Down, " + name2 + " Up-Middle,"
            except TypeError as e:
                logging.error("Curvature difference could not be found.")
                logging.error(e)
                continue

            plotcurves.fig_result3(name, data, data2, differences, differences2, pathExportC)

        logging.info("...Finished!")

    # Figure Result 0 - all pitchers in same figure - Comparison stabilising structures
    plotcurves.fig_result0(datalist, pathExportC)

    # Correlation for whole dataset
    # "spear" for Spearman correlation,     "lin" for linear regression
    # True for indices max dK to -15,       False for 0 to -15
    #plotcurves.correlation(cordata, pathExportC, "spear", False)

    # Figure Result 2 - all pitchers in same figure - Comparison plot for overview
    plotcurves.fig_result2(datalist, pathExportC)

    plotcurves.plot_mean_cis(datalist, pathExportC)


    # Write out comprised data
    plotcurves.export_data(datalist, pathExportC)

    # Reshape list array and set filename
    longestList = len(max(summary_dcurvatures, key=len))
    mylist = np.array([i + [0]*(longestList-len(i)) 
                      for i in summary_dcurvatures]).T
    summaryfile = pathExportC / "Summary_curve_differences.csv"
    # Export data to file
    np.savetxt(summaryfile, mylist, delimiter=',', header=header[:-1])

    logging.info("Summary file with curvature differences saved: {}".format(
                 summaryfile))

    return None


def main():
    # Setup logger
    logging.basicConfig(filename="pyNepenthes.log", level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info('Start...')

    # File structure
    pathCurvature = Path.cwd() / "CurvatureData"
    pathExportC = pathCurvature / "Exported"
    pathSecond = Path.cwd() / "SecondMomentData"
    pathExportI = pathSecond / "Exported"

    # Threshold for second moment of area
    thresh = 10000
    # thresh_venation = 52500
    # Individual threshholds for p1-p7 based on density above 10000 in scans
    # thresh85 = [52718, 49234, 49446, 50531, 48468, 47984, 48719]
    thresh90 = [55700, 52391, 51780, 52955, 51559, 50097, 51523]
    # thresh95 = [59375, 56299, 54538, 55918, 55011, 52936, 55403]
    # Threshold 90 seems to be best one


    # COMMENT: CUrrent version plots with black background and has venation
    # density in plot. For SEB Poster 2022

    # Calculate curvature for all files in pathData and export to pathExport
    #all_files_calc_curvature(pathCurvature, pathExportC)

    # Calculate second moment of area for all data
    #all_samples_calc_I(pathSecond, thresh, thresh90, pathExportI)

    # Combine data for each sample (up, middle, down), build differences, plot
    all_files_plot_curves_export(pathExportC, pathExportI)

    logging.info('...script finished!')


if __name__ == '__main__':
    main()
