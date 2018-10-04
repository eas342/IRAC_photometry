# IRAC_photometry

This repository provides a pipeline that performs high-precision photometric reduction on data from the IRAC instrument of Spitzer Space Telescope. The pipeline is designed for Level 1 IRAC products (i.e. BCD/CBCD files) and only takes input in AORs. An AOR (Astronomical Obeservation Request) is a fundamental unit of Spitzer observing with a unique identification number known as AORKEY. Each AOR consists of multiple BCD/CBCD files which can be thought of as single frame exposures. The end product of this pipeline are two tables that provide photometric analysis for individual BCDs and individual AORs (combination of several BCDs). 

## Dependencies

- ### Python Packages
  astropy, photutils, mirpyidl, tqdm

- ### IDL programs
  - [irac_aphot_corr.pro](http://irsa.ipac.caltech.edu/data/SPITZER/docs/dataanalysistools/tools/contributed/irac/iracaphotcorr/)
  - [irac_aphot_corr_cryo.pro](http://irsa.ipac.caltech.edu/data/SPITZER/docs/dataanalysistools/tools/contributed/irac/iracaphotcorrcryo/)


## Important Files/Folders In This Repository

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
  - The pipeline itself does not generate plots. However, it holds all the plots I have made using data genarated from the pipeline.  


## How To Use

- ### From Script/IDE
  - #### Calling sequence
    The pipeline can be called from a script or an IDL like a python module (provided that AOR_pipeline.py is in the same folder as the script or its path is specified). Here is an example:
    ```python
    import AOR_Pipeline as pipeline
    aor_data, img_data, prob_aor = pipeline.run(AORs, sky, r, rIn, rOut, ap_corr, pixArea, mission, channel, N)
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
    - **ap_corr**: This is the aperture correction factor (_float_). It depends on the radius combination (r, rIn, rOut) and IRAC channel. You can find the proper correction factor from table 4.7 in this [instrument handbook](https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/27/).
    - **pixArea**: Pixel size (_float_) in arcseconds("). Find the appropriate pixel size from table 2.1 in the [instrument handbook](https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/5/). 
    - **mission**: Two possible values (_string_) -'Cryogenic' or 'Warm' (That exact spelling and case). The AORs you provide must all be from the cryogenic mission or the warm mission; it can't take mixed AORs. The spitzer cryogenic mission ended on May 15th, 2009. So, sort input AORs based on that date. 
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
  
- ### From Terminal/Command line
  - #### Calling sequence
    The pipeline can be called from the terminal like any other python file: `python AOR_Pipeline.py`. Once it's called, it will print the list of input parameters with examples. Parameters need to be separated by a semi-colon (;) and no additional spaces or characters can be used. Here is what it will look like when the pipeline is called from terminal:
    ![screenshot](https://github.com/rafia37/IRAC_photometry/blob/master/pp_from_terminal.png "Running AOR_pipeline.py from terminal.")
  - #### Input Parameters
  - #### Output Parameters
  

## Systematics Being Applied Through The Pipeline