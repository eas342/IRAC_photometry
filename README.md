# IRAC_photometry

This repository provides a pipeline that performs high-precision photometric reduction on data from the IRAC instrument of Spitzer Space Telescope. The pipeline is designed for Level 1 IRAC products (i.e. BCD/CBCD files) and only takes input in AORs. An AOR (Astronomical Obeservation Request) is a fundamental unit of Spitzer observing with a unique identification number known as AORKEY. Each AOR consists of multiple BCD/CBCD files which can be thought of as single frame exposures. The end product of this pipeline are two tables that provide photometric analysis for individual BCDs and individual AORs (combination of several BCDs). 

## Contents
- [Dependencies](#dependencies)
- [Contents Of This Repository](#contents-of-this-repository)
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
  - [AOR Table](#aor-table)
  - [IMG Table](#img-table)
  - [Radius Test Table](#radius-test-table)

## Dependencies

- ### Python Packages
  - astropy (Version 2.0 or higher)
  - photutils
  - mirpyidl
  - tqdm

- ### IDL
  - IDL (Version 8.0 or higher)
  - [irac_aphot_corr.pro](http://irsa.ipac.caltech.edu/data/SPITZER/docs/dataanalysistools/tools/contributed/irac/iracaphotcorr/)


## Contents Of This Repository

- ### Pipeline.py
  This is the actual pipeline. You can run it directly from terminal or from another script. Its uses, inputs and outputs are discussed in a later section.
  
- ### functions.py
  This is a general library of functions that is used by the pipeline. The pipeline is dependent on this file so, it has to be in the same folder as the pipeline.
  
- ### Reduction_Data_&_Logs (Folder)
  Data files generated from the pipeline are stored in this folder. Every time the pipeline is run, it generates 3 files:
  - **\*_aor_data.csv**: These files provide the average flux per AOR and other important information for individual AORs such as aorkey, date of observation, dither information etc.
  - **\*_img_data.csv**: These files provide flux from every single image from every AOR. It includes image specific information like source flux, gaussian fitted center position, FWHMs of fitted gaussian, problem in an image (if any) etc.
  - **\*_log.txt**: These files provide information about the run such as: date of run, target info, mission, instrument, radius used, terminal calling sequence etc.
  
- ### General_Plots (Folder)
  The pipeline itself does not generate plots. However, it holds all the plots I made using data genarated from the pipeline. 
  
- ### Correction_files (Folder)
  This folder contains the idl codes that I used to apply systematics and some correction images as fits files. 
  
- ### aorcrd_files (Folder)
  When reducing data for multiple targets using single or multiple AORs, the pipeline requires a csv file that lists AORs with corresponding target coordinates. In this folder, I keep all those csv files that the pipeline requires as input. More on the format and usage of these files can be found in the [How To Use](#how-to-use) section.
  
- ### Experiments (Folder)
  Please ignore this folder. It contains my scrtach work and tests that I did before incorporating various pieces of code into the pipeline.


## How To Use

- ### From Script
  - #### Calling sequence
    The pipeline can be called from a script or an IDE like a python module (provided that Pipeline.py and functions.py are in the same folder as the script or its path is specified through `sys`). Here is an example:
    ```python
    import Pipeline as pipe
    aor_data, img_data, prob_aor = pipe.run(crdFormat, aor_crd, filetype, r, rIn, rOut, ap_corr, channel, pixLen, rejection)
    ```
    The input parameters and outputs are discussed below.
    
  - #### Input Parameters
    - **crdFormat**: A string that specifies the format of target coordinates and the plurality of targets. It can take on one of these 4 possible values: 'single_hms', 'single_deg', 'multiple_hms' or, 'multiple_deg'. The value of *crdFormat* will affect the next imput parameter, *aor_crd*. 
       
    - **aor_crd**: This parameter specifies the AORs and targets you want to work with. The value of this parameter depends on *crdFormat* and so I need to break it down into 4 parts. If:
      - crdFormat = 'single_hms', aor_crd is going to be a *list* with two items: `['path/to/aor', 'hh mm ss dd mm ss']`. Using this option, you can either specify a single target for a single aor or a single target for multiple AORs. If you were to specify the same target for multiple AORs, your calling sequence would look something like:
      ```python
      import Pipeline as pipe
      crdFormat = 'single_hms'
      aor_crd = ['/home/user/data/hd165459/r*', '18 02 30.7410086899 +58 37 38.157415821'] #wildcard provides every AOR in that directory
      aor_data, img_data, prob_aor = pipe.run(crdFormat, aor_crd, 'bcd', 10, 12, 20, 1.0, 1, 1.221)
      ```
      - crdFormat = 'single_deg', aor_crd is going to be a *list* with two items: `['path/to/aor', 'ra_in_deg dec_in_deg']`. This is the same as 'single_hms', the only difference being the coordinates have to be provided in degrees.
      - crdFormat = 'multiple_hms', aor_crd is going to be the filename of a csv file that contains two columns named "AOR Filepath" and "Target Coordinates". Both columns must have string datatype. Target coordinates must be in this format: 'hh mm ss dd mm ss'. With this option, you can specify multiple targets for multiple AORs and multiple targets for single AOR. You could also do single target for multiple AOR and single target for single AOR but that would be redundant because you could do them with the 'single_hms' option which saves you the trouble of making a csv file. You can check the "aorcrd_files" folder in this repository, to find examples of these csv files. If you wanted to test with one of the csv files in that folder, you calling sequence would look like:
      ```python
      import Pipeline as pipe
      crdFormat = 'multiple_hms'
      aor_crd = 'aorcrd_files/taurus_skyCrd.csv'
      aor_data, img_data, prob_aor = pipe.run(crdFormat, aor_crd, 'bcd', 10, 12, 20, 1.0, 1, 1.221)
      ```
      - crdFormat = 'multiple_deg', aor_crd would also be a filename of a csv file as described for "multiple_hms". The only difference is that the coordinates have to be in degrees. An example of *aor_crd*, when *crdFormat* is 'multiple_deg' is `aorcrd_files/upperSco_skyCrd.csv`.
      
    - **filetype**: A *string* that specifies the type of file being used. e.g. 'bcd', 'cbcd', 'sub2d' etc.
    - **r, rIn, rOut**: The pipeline performs photometry using _photutils_ aperture photomotery where a circular aperture is used for source and circular annulus aperture is used for background. r, rIn and rOut are source aperture radius, inner background aperture radius and outer background aperture radius respectively. They can be _int_ or _float_.
    - **ap_corr**: This is the aperture correction factor (_float_). It depends on the radius combination (r, rIn, rOut) and IRAC channel. You can find the proper correction factor from the table provided in [Aperture Correction](#aperture-correction) section.
    - **channel**: IRAC channel that was used to collect data. Must be an integer between 1 and 4 (inclusive).
    - **pixLen**: Length of a pixel (_float_) in arcseconds("). Find the appropriate pixel size from this [table](https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/5/):
    
    | Channel | Pixel Size (") |
    | ------- | -------------- |
    | 1       | 1.221          |
    | 2       | 1.213          |
    | 3       | 1.222          |
    | 4       | 1.220          |
    
    - **rejection**: A boolean, "True" by default. If set to "True", outliers will be rejected by rejecting the maximum and minimum values in an AOR provided that the AOR has more than 10 workable files and the standard deviation of the set is higher than 2%. If "False", no outlier rejection will be done before averaging for an AOR.
    
  - #### Output Quantities
    - **aor_data**: This is the table that provides average flux per AOR and other important AOR-specific information like aorkey, date of obsservation, dither information etc. If the pipeline is run from a script, it will not get saved as `runName_aor_data.csv` automatically. You can choose to save the table as a csv file from your script with:
      ```python
      from astropy.io import ascii
      ascii.write(aor_data, 'Reduction_Data_&_Logs/runName_aor_data.csv', delimiter = ',', overwrite = False)
      ```
    - **img_data**: This is the table that would be saved as `runName_img_data.csv` if you ran the pipeline from terminal. Since you're running it from a script, you will get this astropy table which you can save with `ascci.write()` (as shown above) if you wish to.
    - **prob_aor**: This is a numpy array of AORKEYs of those AORs that caused some kind of problem and were not run by the pipeline. This is there so that you can check out these AORs and figure out what the problem is with these AORs. Common problems are usually that the AOR was empty, the target was not in FOV in any of the files in that AOR etc.
  
- ### From Terminal
  - #### Calling sequence
    The pipeline can be called from the terminal like any other python file: `python Pipeline.py -h`. The '-h' flag at the end is the help menu which will show you the usage and arguments that this Pipeline takes. The pipeline has two positional arguments and 8 optional arguments (excluding the help argument). Which means the calling sequence can be as simple as: 
    ```shell
    python Pipeline.py 'multiple_hms' 'taurus_skyCrd.csv'
    ```
    or as big as,
    ```shell
    python Pipeline.py 'multiple_deg' 'upperSco_skyCrd.csv' -f 'sub2d' -r 8 12 20 -a 1.013 -cp 2 1.213 -o -rn "upperSco_test" -t "Upper Sco stars" -c "Comment describing your run"
    ```
    All arguments must be separated by "space" or "tab". String arguments must be used with quotes. It's safer to use `"` instead of `'` for strings in case your string contains a `'`. 
  - #### Positional Input Parameters
    Positional Arguments are the ones that must be provided for the pipeline to run and they must be in the correct sequence. There are two positional arguments for the pipeline:
    - **crdFormat**: A string that specifies the format of target coordinates and the plurality of targets. It can take on one of these 4 possible values: 'single_hms', 'single_deg', 'multiple_hms' or, 'multiple_deg'. The value of *crdFormat* will affect the next imput parameter, *aor_crd*. 
       
    - **aor_crd**: This parameter specifies the AORs and targets you want to work with. The value of this parameter depends on *crdFormat* and so I need to break it down into 4 parts. If:
      - crdFormat = 'single_hms', aor_crd must be 2 string arguments `'path/to/aor' 'hh mm ss dd mm ss'`. Using this option, you can either specify a single target for a single aor or a single target for multiple AORs. If you were to specify the same target for multiple AORs, you would have to use a wildcard in the filepath and your calling sequence would look something like:
      ```shell
      python Pipeline.py 'single_hms' '/home/user/data/hd165459/r*' '18 02 30.7410086899 +58 37 38.157415821'
      ```
      - crdFormat = 'single_deg', This is the same as 'single_hms', the only difference being the coordinates have to be provided in degrees.
      - crdFormat = 'multiple_hms', aor_crd is going to be the filename of a csv file that contains two columns named "AOR Filepath" and "Target Coordinates". Both columns must have string datatype. Target coordinates must be in this format: 'hh mm ss dd mm ss'. With this option, you can specify multiple targets for multiple AORs and multiple targets for single AOR. You could also do single target for multiple AOR and single target for single AOR but that would be redundant because you could do them with the 'single_hms' option which saves you the trouble of making a csv file. You can check the "aorcrd_files" folder in this repository, to find examples of these csv files. If you wanted to test with one of the csv files in that folder, you calling sequence would look like:
      ```shell
      python Pipeline.py 'multiple_hms' 'taurus_skyCrd.csv'
      ```
      - crdFormat = 'multiple_deg', aor_crd would also be a filename of a csv file as described for "multiple_hms". The only difference is that the coordinates have to be in degrees. An example of *aor_crd*, when *crdFormat* is 'multiple_deg' is `aorcrd_files/upperSco_skyCrd.csv`.
      
  - #### Optional Input Parameters
    Optional arguments are optional as the name suggests. If these optional arguments are not specified, the Pipeline will run with default values. So you only have to specify these arguments if you're not satisfied with the default values. They are not required to be in a prticular sequence. So you can use any optional argument in any order but they have to start after the positional arguments. Here is a list of the optional arguments:
    - **Help**: Specified with the '-h' or '--help' flag. Shows the help menu and exits. This is the only optional argument that can be specified without any positional argument. Calling sequence looks like `python Pipeline.py -h`
    - **radius**: Specified with the '-r' or '--radius' flag. It is followed by 3 integer arguments: r, rIn, rOut. The pipeline performs photometry using _photutils_ aperture photomotery where a circular aperture is used for source and circular annulus aperture is used for background. r, rIn and rOut are source aperture radius, inner background aperture radius and outer background aperture radius respectively. **Default value is 10 12 20**.
    - **ap_corr**: Specified with the '-a' or '--aperture-correction' flag. It is followed by 1 *float* argument. This argument is the aperture correction factor. It depends on the radius combination (r, rIn, rOut) and IRAC channel. You can find the proper correction factor from the table provided in [Aperture Correction](#aperture-correction) section. **Default value is 1.000**
    - **channel & pixel length**: Specified with the '-cp' or '--ch-px' flag. It is followed by 2 *int/float* arguments. The first argument is IRAC channel number, which can be an integer between 1 and 4 (inclusive). The second argument is the length of a pixel for that channel, which you can find from this [table](https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/5/): (**Default value is 1 1.221**)
    
    | Channel | Pixel Size (") |
    | ------- | -------------- |
    | 1       | 1.221          |
    | 2       | 1.213          |
    | 3       | 1.222          |
    | 4       | 1.220          |
    
    - **rejection**: Specified with the '-o' or '--outlier-rejection' flag. It is not followed by any argument. If the flag is specified, outliers will be rejected by rejecting the maximum and minimum values in an AOR provided that the AOR has more than 10 workable files and the standard deviation of the set is higher than 2%. Otherwise, no outlier rejection will be done before averaging for an AOR. **Default value is False (as in, not specified unless you specify it)**
    - **run name**: Specified with the '-rn' or '--run-name' flag. Followed by 1 string argument. This argument will be used as the run name for the output files. **Default value is "Test_run"**. 
    - **Target Name**: Specified with the '-t' or '--target-name' flag. Followed by 1 string argument. This is the name of your target or targets that will be included in the output log file. **Default value is "Many"**
    - **Comments**: Specified with the '-c' or '--comments' flag. Followed by 1 string argument. This is the comment describing your run that will be included in the output log file. **Default value is "No comment was specified"**.
    
  - #### Output Quantities
    Once the pipeline completes running, 2 data tables as csv files and a text file will be genarated and saved in the *Reduction_Data_&_Logs* folder. These are the `runName_aor_data.csv`, `runName_img_data.csv` and `runName_log.txt` files described in [Contents Of This Repository](#contents-of-this-repository) section.
    
    On the terminal, as the pipeline runs, there are two progress bars. The first one shows the progress of AORs, while the second one shows the progress of images in an AOR. 
  

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
  Run summary for each individual run can be easily obtained by running the *run_summar.sh* file in the "Reductions_Data_&_Logs" folder in this is repository. The *comments* parameter gives a smmary of each individual run.
  Here, I still provide the summary of the first 12 runs since they can be easily grouped. For the first 12 runs, I have run the pipeline for 4 groups of data: cryogenic mission data for HD 165459, warm mission data for HD 165459, cryogenic mission data for BD +60 1753, warm mission data for BD +60 1753. So the runs can be grouped into sets of 4. 
  - Run 1-4 were test runs. No outliers were rejected and warm mission data did not have mission-specific systematics applied.
  - Run 5-8 are proper runs for channel 1 data. Which means mission-specific systematics were applied and outliers were rejected in these runs.
  - Run 9-12 are proper runs for channel 2 data.
  - Run 13+ were runs on various sets of target with a new version of code. From this point, it's difficult to group runs based on data sets so, please use `run_summary.sh` for getting a summary of individual runs.
      
- ### Plots Summary
  I made some plots using the data tables generated by the pipeline and stored them in *General_Plots* folder.
  - *AORtest.png* : This is a time series made using two sample AORs for both calibrator stars (HD 165459 & BD +60 1753). It shows the distribution within an AOR.
  - *calStars_ch1.png* : This is a time series for both calibrator stars using channel 1 data. The cryogenic mission data and warm mission data are shown on the same plot. It seems that there is a shift between these two missions.
  - *calStars_ch1_pDith.png* : This is the time series for both stars using channel 1 data that met certain dither conditions. The AORs used for this plot has cycling dither pattern, large dither scale, 5 dither positions and a frame time of 0.4s for HD 165459/ 2.0s for BD +60 1753.
  - *calStars_ch1_zm.pdf* : Same as *calStars_ch1.png* but zoomed in to show the shift and slope in warm mission data.
  - *calStars_ch1_pDith_zm.pdf* : Same as *calStars_ch1_pDith.png* but zoomed in to show the shift and slope in warm mission data.

- ### AOR Table
  The AOR table provides the average flux per AOR and other important information for individual AORs such as aorkey, date of observation, dither information etc. These tables are saved as *\*_aor_data.csv* in the *Reduction_Data_&_Logs* folder and is one of the outputs of the pipeline. It has 16 columns and here is what each of the column names mean:
  - **AORKEY**: The 6-digit AOR identifier. 
  - **Target RA**: RA of the target for thatspecific AOR.
  - **Target Dec**: Dec of the target for thatspecific AOR.
  - **DateObs**: Date of observation as obtained from image header.
  - **Mission**: Spitzer mission: Cryogenic/Warm.
  - **Read Mode**: How the image was read from detector: FULL/SUB array mode.
  - **Workable/Total Files in AOR**: Number of files in AOR that was useful in calculating average/Total number of files in the AOR
  - **Cycling DPattern**: Whether the Dither pattern was cyclic or not.
  - **DScale**: Dither scale. 
  - **DPosition**: Dither Position.
  - **FTime (sec)**: Integration time in Seconds
  - **Time (MJD)**: Date of Observation in Modified Julian Date
  - **Flux (mJy)**: Average Flux (in mJy) of target from this AOR as calculated by pipeline.
  - **Error (mJy)**: Uncertainty in average flux.
  - **Spread (%)**: Spread of distribution in percentage for this AOR.
  - **Outliers Rejected**: Number of outliers that was rejected. No otlier is rejected if number of workable files in AOR is less than 10.

- ### IMG Table
  The IMG table provides flux from every single image from every AOR. It includes image specific information like source flux, gaussian fitted center position, FWHMs of fitted gaussian, problem in an image (if any) etc. These tables are saved as *\*_img_data.csv* in the *Reduction_Data_&_Logs* folder and is one of the outputs of the pipeline. It has 14 columns and here is what each of the column names mean:
  - **File#**: File/Image number in an AOR. The count starts over at 1 when a new AOR is encountered. 
  - **Coord Conversion Issue**: Flags 'Y' if there was an issue with WCS and the given coordinate could not be converted to pixel coordinates. 
  - **Centroiding Issue**: Flags 'Y' if gaussian fitting fails and the PSF center could not be determined.
  - **Bad Center Guess**: Flags 'Y' if gaussian fitted center deviates by 2 px from provided target coordinates.
  - **Not In FOV**: Flags 'Y' if the PSF center lies beyond the detector array limits.
  - **Ap Out Of Bound**: Flags 'Y' if any of the photometry apertures go beyond detector array limits.
  - **Xc**: Gaussian fitted x-center coordinate.
  - **Yc**: Faussian fitted y-center coordinate.
  - **Fx**: FWHM in x direction of fitted gaussian.
  - **Fy**: FWHM in y direction of fitted gaussian.
  - **Time[MJD]**: Date of observation of the image in Modified Julian Date.
  - **Raw_Flux**: Flux extracted from source aperture
  - **Bkg_Flux**: Flux extracted form bakcground aperture
  - **Res_Flux**: Residual flux obtained by subtracting *Bkg_Flux* from *Raw_Flux*.
  
- ### LOG File:
  This is the last output of the pipeline: a txt file containing informations about the run. These files are saved as *\*_log.txt* in the *Reduction_Data_&_Logs* folder and is one of the outputs of the pipeline. It has 8 rows here is what each of these rows mean:
  - **Date reduced**: Date the data was reduced using the pipeline.
  - **Input Parameters**: Parameters as specified through the terminal. If not specified, default values are used.
  - **Instrument**: Instrument used in collecting the data. 
  - **File Type**: Type of file that was used.
  - **Target**: Name or description of target(s).
  - **Radius Used**: Source, inner background and outer background aperture radius that was used to perform photometry.
  - **Problem AORs**: AORKEYs of the AORs that could not be used.
  - **Comments**: Brief description of the run.
  
- ### Radius Test Table
  In the *Reduction_Data_&_Logs* folder, I have 4 csv files that starts with *"rad_test"*. These are radius tests for the 4 groups of data. I have selected a sample of about 30 AORs from both missions, for both stars and 10 different radius combinations. The combinations were chosen from the [aperture correction table](http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/27/#_Toc410728317) so that proper aperture corrections could be aplied.
