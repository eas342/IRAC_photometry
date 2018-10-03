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
  - *run?_aor_data.csv*: This file provides the average flux per AOR and other important information for individual AORs such as aorkey, date of obsservation, dither information etc.
  - *run?_img_data.csv*: This file provides flux from every single image from every AOR. It includes image specific information like source flux, gaussian fitted center position, FWHMs of fitted gaussian, problem in an image (if any) etc.
  - *run?_log.txt*: This file provides information about the run such as: date of run, target info, mission, instrument, radius used, terminal calling sequence etc.
  
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
    - **AORs**: A 2D list/array that contains the filenames of desired files. For examples, if you want to work with 2 AORs that have 2 and 3 BCD files respectively, the list would look something like:
      ```python
      AORs = [['aor1_file1_bcd.fits', 'aor1_file2_bcd.fits'],
             ['aor2_file1_bcd.fits', 'aor2_file2_bcd.fits', 'aor2_file3_bcd.fits']]
      ```
    - **sky**:
    - **r, rIn, rOut**:
    - **ap_corr**:
    - **pixArea**:
    - **mission**:
    - **channel**:
    - **N**:
    
  - #### Output Quantities
    - aor_data:
    - img_data:
    - prob_aor:
  
- ### From Terminal/Command line
  - #### Calling sequence
    The pipeline can be called from the terminal like any other python file: `python AOR_Pipeline.py`. Once it's called, it will print the list of input parameters with examples. Parameters need to be separated by a semi-colon (;) and no additional spaces or characters can be used.
  - #### Input Parameters
  - #### Output Parameters
  

## Systematics Being Applied Through The Pipeline