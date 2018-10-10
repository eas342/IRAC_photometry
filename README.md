# IRAC_photometry

This repository provides a pipeline that performs high-precision photometric reduction on data from the IRAC instrument of Spitzer Space Telescope. The pipeline is designed for Level 1 IRAC products (i.e. BCD/CBCD files) and only takes input in AORs. An AOR (Astronomical Obeservation Request) is a fundamental unit of Spitzer observing with a unique identification number known as AORKEY. Each AOR consists of multiple BCD/CBCD files which can be thought of as single frame exposures. The end product of this pipeline are two tables that provide photometric analysis for individual BCDs and individual AORs (combination of several BCDs). 

## Contents
- [Dependencies](#dependencies)
- [Important Files In This Repository](#important-files-in-this-repository)
- [How To Use](#how-to-use)
  - [From Script](#from-script)
  - [From Terminal](#from-terminal)
- [Systematics Being Applied Through The Pipeline](#systematics-being-applied-through-the-pipeline)
  - [Array Location Dependent Correction](#array-location-dependent-correction)
  - [Pixel Phase Correction](#pixel-phase-correction)
  - [Photometric Calibration](#photometric-calibration)
  - [Aperture Correction](#aperture-correction)
  - [Linearity Correction](#linearity-correction)
  - [Outlier Rejection](#outlier-rejection)
- [Summary Of Current Tables And Plots](#summary-of-current-tables-and-plots)
  - [Run Summary](#run-summary)
  - [Plots Summary](#plots-summary)
  - [Radius Test Table](#radius-test-table)

## Dependencies

- ### Python Packages
  astropy, photutils, mirpyidl, tqdm

- ### IDL programs
  - [irac_aphot_corr.pro](http://irsa.ipac.caltech.edu/data/SPITZER/docs/dataanalysistools/tools/contributed/irac/iracaphotcorr/)
  - [irac_aphot_corr_cryo.pro](http://irsa.ipac.caltech.edu/data/SPITZER/docs/dataanalysistools/tools/contributed/irac/iracaphotcorrcryo/)


## Important Files In This Repository

- ### AOR_Pipeline.py
  This is the actual pipeline. You can run it directly from terminal or from another script. Its uses, inputs and outputs are discussed in a later section.
  
- ### functions.py
  This is a general library of functions. The pipeline is dependent on this file so, it has to be in the same folder as the pipeline.
  
- ### Reduction_Data_&_Logs (Folder)
  Data files generated from the pipeline are stored in this folder. Every time the pipeline is run, it generates 3 files:
  - **run?_aor_data.csv**: This file provides the average flux per AOR and other important information for individual AORs such as aorkey, date of obsservation, dither information etc.
  - **run?_img_data.csv**: This file provides flux from every single image from every AOR. It includes image specific information like source flux, gaussian fitted center position, FWHMs of fitted gaussian, problem in an image (if any) etc.
  - **run?_log.txt**: This file provides information about the run such as: date of run, target info, mission, instrument, radius used, terminal calling sequence etc.
  
- ### General_Plots (Folder)
  The pipeline itself does not generate plots. However, it holds all the plots I have made using data genarated from the pipeline. 
  
- ### Correction_files (Folder)
  This file contains the idl codes that I used to apply systematics and some correction images as fits files. 


## How To Use

- ### From Script
  - #### Calling sequence
    The pipeline can be called from a script or an IDE like a python module (provided that AOR_pipeline.py is in the same folder as the script or its path is specified). Here is an example:
    ```python
    import AOR_Pipeline as pipeline
    aor_data, img_data, prob_aor = pipeline.run(AORs, sky, r, rIn, rOut, ap_corr, pixArea, channel, N)
    ```
    The input parameters and outputs are discussed below.
    
  - #### Input Parameters
    - **AORs**: A 2D _list/array_ that contains the filenames of desired files. For examples, if you want to work with 2 AORs that have 2 and 3 BCD files respectively, the list would look something like:
      ```python
      AORs = [['aor1_file1_bcd.fits', 'aor1_file2_bcd.fits'],
             ['aor2_file1_bcd.fits', 'aor2_file2_bcd.fits', 'aor2_file3_bcd.fits']]
      ```
      Where the strings in the _list_ are the proper filepath/filename.
    - **sky**: An _astropy SkyCoord object_ that holds the sky coordinates (RA & Dec) of target in degrees. So, you could define _sky_ in any of the following ways:
      ```python
      from astropy import units as u
      from astropy.coordinates import SkyCoord
      sky = SkyCoord(10.625, 41.2, frame='icrs', unit='deg')
      sky = SkyCoord('00 42 30 +41 12 00', unit=(u.hourangle, u.deg))
      sky = SkyCoord('00:42.5 +41:12', unit=(u.hourangle, u.deg))
      ```
    - **r, rIn, rOut**: The pipeline performs photometry using _photutils_ aperture photomotery where a circular aperture is used for source and circular annulus aperture is used for background. r, rIn and rOut are source aperture radius, inner background aperture radius and outer background aperture radius respectively. They can be _int_ or _float_.
    - **ap_corr**: This is the aperture correction factor (_float_). It depends on the radius combination (r, rIn, rOut) and IRAC channel. You can find the proper correction factor from the table provided in [Aperture Correction](#aperture-correction) section.
    - **pixArea**: Pixel size (_float_) in arcseconds("). Find the appropriate pixel size from this [table](https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/5/):
    
    | Channel | Pixel Size (") |
    | ------- | -------------- |
    | 1       | 1.221          |
    | 2       | 1.213          |
    | 3       | 1.222          |
    | 4       | 1.220          |
    
    - **ROMode**:
    - **channel**: Four possible values (_int_) - 1,2,3 or 4. This is the irac channel that was used to collect data. 
    - **N**: Outlier rejection factor. If the value is 'n/a', no outlier rejection will happen. Outlier rejection will happen for any other value you provide.
    
  - #### Output Quantities
    - **aor_data**: This is the table that provides average flux per AOR and other important AOR-specific information like aorkey, date of obsservation, dither information etc. If the pipeline is run from a script, it will not get saved as `run?_aor_data.csv` automatically. You can choose to save the table as a csv file from your script with:
      ```python
      from astropy.io import ascii
      ascii.write(aor_data, 'Reduction_Data_&_Logs/run?_aor_data.csv', delimiter = ',', overwrite = False)
      ```
    - **img_data**: This is the table that would be saved as `run?_img_data.csv` if you ran the pipeline from terminal. Since you're running it from a script, you will get this astropy table which you can save with `ascci.write()` (as shown above) if you wish to.
    - **prob_aor**: This is a numpy array of AORKEYs of those AORs that caused some kind of problem and were not run by the pipeline. This is there so that you can check out these AORs and figure out what the problem is with these AORs. Common problems are usually that the AOR was empty, the target was not in FOV in any of the files in that AOR etc.
  
- ### From Terminal
  - #### Calling sequence
    The pipeline can be called from the terminal like any other python file: `python AOR_Pipeline.py`. Once it's called, it will print the list of input parameters with examples. Parameters need to be separated by a semi-colon (;) and no additional spaces or characters can be used. Here is what it will look like when the pipeline is called from terminal:
    ![screenshot](https://github.com/rafia37/IRAC_photometry/blob/master/pp_from_terminal.png "Running AOR_pipeline.py from terminal.")
  - #### Input Parameters
    Since the input parameters have to be specified in the terminal, you don't have to worry about using quotes when providing a string. For example, instead of writing *'Cryogenic'*, write *Cryogenic* without the quotes. Just list the parameters separated by a semi-colon(;) in the format specified below:
    - **File Path**: File path up to aor names. Could be a file path with wildcard (preferred) or a comma separated list. So if you have 3 AORs in a folder you could specify it in the following ways:  
    `/your/file/path/aor*/` (preferred method) or,  
    `/your/file/path/aor1/,/your/file/path/aor2/,/your/file/path/aor3/`
    - **File Type**: bcd or cbcd. Must be all lowercase.
    - **Target Name**: Name of target. e.g. HD 165459
    - **Target Coordinates**: Sky coordinates of target in hms for RA and dms for Dec. You can write it as `18 02 30.74 +58 37 38.16`, `18:02:30.74 +58:37:38.16` or `18h02m30.74s +58d37m38.16s`.
    - **Readout Mode**: full/sub
    - **Source Aperture Radius**: Source aperture in native pixel. (e.g. 10)
    - **Inner Background Radius**: Inner background aperture in native pixel. (eg. 12)
    - **Outer Background Radius**: Outer Background aperture in native pixel. (eg. 20)
    - **Channel#**: 1,2,3 or 4. IRAC channel that was used to collect data.
    - **Aperture Correction Factor**: A constant factor that depends on the radius combination and IRAC channel. You can find it from the table provided in [Aperture Correction](#aperture-correction).
    - **Pixel Size**: Size of a pixel in arcseconds("). Find the appropriate pixel size from this [table](https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/5/):
    
    | Channel | Pixel Size (") |
    | ------- | -------------- |
    | 1       | 1.221          |
    | 2       | 1.213          |
    | 3       | 1.222          |
    | 4       | 1.220          |
    
    - **Run#**: For output file name. Should be an integer. This number is what goes into "run?_aor_data.csv". So if you don't want the data to be overwritten, use a new run number.
    - **Outlier Rejection**: If the value is *n/a*, no outlier rejection will happen. Outlier rejection will happen for any other value you provide.
    - **Comments**: Write anything that you would want to include in the log.  
    Don't jump to any argument. e.g. You can't skip *Outlier Rejection* to get to *Comments*. *Comments* must be the 14th argument. You could however, use *n/a*, if you don't want outlier rejection. Here is an example parameter list that could be written on the terminal: `/data1/phot_cal/spitzer/hd165459/cryo/r*/;bcd;HD 165459;18 02 30.74 +58 37 38.16;Cryogenic;10;12;20;1;1.0;1.221;5;10;Your comment goes here.`
  - #### Output Quantities
    Once the pipeline completes running, 2 data tables as csv files and a text file will be genarated and saved in the *Reduction_Data_&_Logs* folder. These are the `run?_aor_data.csv`, `run?_img_data.csv` and `run?_log.txt` files described in *Important Files In This Repository* section.
  

## Systematics Being Applied Through The Pipeline

Three of the corrections (array location dependent correction, pixel phase correction and photometric calibration) are applied using [irac_aphot_corr_cryo.pro](http://irsa.ipac.caltech.edu/data/SPITZER/docs/dataanalysistools/tools/contributed/irac/iracaphotcorrcryo/) for cryogenic mission and [irac_aphot_corr.pro](http://irsa.ipac.caltech.edu/data/SPITZER/docs/dataanalysistools/tools/contributed/irac/iracaphotcorr/) for warm mission. You can find these files in the *Correction_files* folder. The BCD/CBCD files are flat fielded and linearity corrected through Spitzer Science Center's BCD pipeline.

- ### Array Location Dependent Correction
  This is a correction that is required for point source photometry or compact sources. IRAC flatfield is made by imaging the high surface brightness zodical background (ZB). Firstly, since the ZB is extended and effectively uniform over the IRAC FOV, its effective gain is slightly different from a point source effective gain. Secondly, spectrum of the ZB peaks redward of the IRAC filters while spectrum of high color temerature sources like stars peak blueward of IRAC filters and are well on the Rayleigh-Jeans side of blackbody spectrum. Thus, IRAC photometry of warm sources require a correction for the spectral slope change between zodical light and Rayleigh-Jeans spectra. Lastly, the effective filter bandpass of IRAC varies as a function of angle of incidence. Thus, all objects in the IRAC FOV needs to be corrected based on their location on the array.  
  All these effects can be corrected by using correction images from [IRAC instrument webpage](http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/calibrationfiles/locationcolor/). You can also find these correction images in the *Correction_files* folder. The correction images are different for the cryogenic and warm missions. These correction images are for compact, point-like sources Rayleigh-Jeans (stellar, Vega-like) spectrum. To apply these correction images, we need to multiply the value from correction image at the central pixel of target to the observed flux of target from BCD file. This gives the corrected flux value. The variation in measured fluxes can reach up to 10% peak-to-peak. This is larger than any other source of uncertainty in IRAC calibration.
  
- ### Pixel Phase Correction
  Due to the variation of quantum efficiency across a pixel, the measured flux density of a point source depends on the exact loaction where the peak of the point spread function (PSF) falls on a pixel. This effect is most severe in channel 1 and the correction of this effect can be as much as 4% peak-to-peak. To correct for this effect, we can use [correction images](http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/20/#_Toc410728308) that are a two dimensional function of pixel phase. You can also find these images in the *Correction_files* folder.
  
- ### Photometric Calibration
  To obtain an absolute flux calibration, about 11 standard stars in the continuous viewing zone were observed in each instrument campaign. IRAC provides flux calibration constants for each channel that needs to be multipled to the measured flux to obtain calibrated flux. The proper constant can be found from a [table](http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/17/#_Toc410728305) in the IRAC instrument handbook but I am providing it here as well:
  
  | Channel | Constant (MJy/sr)/(DN/sec) |
  | ------- | -------------------------- |
  | 1       | 0.1088 (0.1253)            |
  | 2       | 0.1388 (0.1469)            |
  | 3       | 0.5952                     |
  | 4       | 0.2021                     |
  
  The values in parentheses are for warm mission.
  
- ### Aperture Correction
  Aperture correction is a correction that compensates for flux that is lost from point source observation. IRAC instrument handbook provides an estimate of aperture corrections based on channel and radius combination used. It can be found from this [page](http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/27/#_Toc410728317) but I'm providing it here as well:
  
  | Radius on source (px) | Background annulus (px) | Channel 1 | Channel 2 | Channel 3 | channel 4 |
  | --------------------- | ----------------------- | --------- | --------- | --------- | --------- |
  | infinite              | N/A                     | 0.944     | 0.937     | 0.772     | 0.737     |
  | 10                    | 12-20                   | 1.000     | 1.000     | 1.000     | 1.000     |
  | 8                     | 12-20                   | 1.011     | 1.013     | 1.011     | 1.017     |
  | 6                     | 12-20                   | 1.032     | 1.036     | 1.030     | 1.051     |
  | 5                     | 12-20                   | 1.047     | 1.048     | 1.054     | 1.064     |
  | 5                     | 5-10                    | 1.060     | 1.063     | 1.063     | 1.084     |
  | 4                     | 12-20                   | 1.070     | 1.080     | 1.076     | 1.087     |
  | 3                     | 12-20                   | 1.112     | 1.112     | 1.118     | 1.213     |
  | 3                     | 3-7                     | 1.125     | 1.120     | 1.135     | 1.221     |
  | 2                     | 12-20                   | 1.208     | 1.220     | 1.349     | 1.554     |
  | 2                     | 2-6                     | 1.215     | 1.233     | 1.366     | 1.568     |
  
  
- ### Linearity Correction
  TBD
  
- ### Outlier Rejection
  Outliers are rejected per AOR by rejecting the maximum and minimum flux values if the spread of distribution is greater than 2% and the AOR has at least 10 files to work with. If you wish to do sigma clipping, you can do so using the `functions.py` script:
  ```python
  import functions as func
  new_array, mask = func.sigma_clipping(array, N)  
  ```
  Where,
  - **array** : Set of data you want to perform sigma clipping on. Type = Array/list
  - **N** : Number of sigmas to be used for clipping. Data will be clipped by N*sigma. Type = int/float
  - **new_array** : Outlier rejected array. Type = Array
  - **mask** : Mask that was used to clip the original *array*. It contains boolean values and has the same length as *array*. Type = Array
  

## Summary Of Current Tables And Plots

- ### Run Summary
  As of now, I have run the pipeline 12 times and so, I have 2 csv files and a txt file for 12 runs in the *Reduction_Data_&_Logs* folder. I have run the pipeline for 4 groups of data: cryogenic mission data for HD 165459, warm mission data for HD 165459, cryogenic mission data for BD +60 1753, warm mission data for BD +60 1753. So the runs can be grouped into sets of 4. 
  - Run 1-4 were test runs. No outliers were rejected and warm mission data did not have mission-specific systematics applied.
  - Run 5-8 are proper runs for channel 1 data. Which means mission-specific systematics were applied and outliers were rejected in these runs.
  - Run 9-12 are proper runs for channel 2 data.
  
  From all these runs, I consistently found that *HD 165459* was about 16 times brighter than *BD +60 1753*, both for channel 1 and 2.
  
- ### Plots Summary
  I made some plots using the data tables generated by the pipeline and stored them in *General_Plots* folder.
  - *AORtest.png* : This is a time series made using two sample AORs for both calibrator stars (HD 165459 & BD +60 1753). It shows the distribution within an AOR.
  - *calStars_ch1.png* : This is a time series for both calibrator stars using channel 1 data. The cryogenic mission data and warm mission data are shown on the same plot. It seems that there is a shift between these two missions.
  - *calStars_ch1_pDith.png* : This is the time series for both stars using channel 1 data that met certain dither conditions. The AORs used for this plot has cycling dither pattern, large dither scale, 5 dither positions and a frame time of 0.4s for HD 165459/ 2.0s for BD +60 1753.
  - *calStars_ch1_zm.pdf* : Same as *calStars_ch1.png* but zoomed in to show the shift and slope in warm mission data.
  - *calStars_ch1_pDith_zm.pdf* : Same as *calStars_ch1_pDith.png* but zoomed in to show the shift and slope in warm mission data.
  
- ### Radius Test Table
  In the *Reduction_Data_&_Logs* folder, I have 4 csv files that starts with *"rad_test"*. These are radius tests for the 4 groups of data. I have selected a sample of about 30 AORs from both missions, for both stars and 10 different radius combinations. The combinations were chosen from the [aperture correction table](http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/27/#_Toc410728317) so that proper aperture corrections could be aplied.
