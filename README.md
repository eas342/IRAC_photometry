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
  This is the actual pipeline.
- ### functions.py
- ### Reduction_Data_&_Logs (Folder)
- ### General_Plots (Folder)


## How To Use

- ### From Script/IDE
- ### From Terminal/Command line