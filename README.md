# IRAC_photometry

This repository provides a pipeline that performs high-precision photometric reduction on data from the IRAC instrument of Spitzer Space Telescope. The pipeline is designed for Level 1 IRAC products (i.e. BCD/CBCD files) and only takes input in AORs. An AOR (Astronomical Obeservation Request) is a fundamental unit of Spitzer observing with a unique identification number known as AORKEY. Each AOR consists of multiple BCD/CBCD files which can be though of as single frame exposures. The end product of this pipeline are two tables that provide photometric analysis for individual BCDs and individual AORs (combination of several BCDs). 

### Dependencies

- #### Python Packages
- #### IDL programs


### Important Files/Folders In This Repository

- #### AOR_Pipeline.py
- #### functions.py
- #### Reduction_Data_&_Logs (Folder)
- #### General_Plots (Folder)


### How To Use

- #### From Script/IDE
- #### From Terminal/Command line