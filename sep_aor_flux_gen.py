import numpy as np
import functions as func
import pdb, glob
import mirpyidl as idl
from tqdm import tqdm
from astropy.table import Table, vstack
from astropy.io import fits, ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
from datetime import datetime


AORnames = np.sort(glob.glob('/data1/phot_cal/spitzer/hd165459/cryo/r*/'))
AORs = []

for aor in AORnames:
    fnames = np.sort(glob.glob(aor + 'ch1/bcd/*_bcd.fits'))
    AORs.append(fnames)

#Defining general constants
#---------------------------

#Provide proper sky coordinates in hms
# '18 02 30.7410086899 +58 37 38.157415821' for HD 165459
# '17 24 52.2772360943 +60 25 50.780790994' for BD +60 1753
sky = SkyCoord('18 02 30.7410086899 +58 37 38.157415821', unit=(u.hourangle, u.deg))

r, rIn, rOut = 10, 12, 20

# Find the aperture correction factor from the following link:
# https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/27/
#from the table at the end of that link, select your value accoring to your aperture size and channel
ap_corr = 1.000

pixLen  = 1.221 #arcsec
pixArea = pixLen**2 #arcsec^2
pixArea = pixArea/(206265**2) #Steradian

#----------------------------


#Generate data for each aor 
for i, aor in tqdm(enumerate(AORs)):
    
    result = Table(names = ('AORKEY', 'DateObs', 'Cycling DPattern', 'DScale', 'DPosition', 'FTime (sec)', 'Time (MJD)', 'Flux (mJy)', 'Error (mJy)', 'Spread (%)'), dtype = ('i4', 'S25', 'S5', 'S10', 'S5', 'f8', 'f8', 'f8', 'f8', 'f8'))
    
    
    #Performing photometry on every image in an AOR
    data, header = func.single_target_phot(aor, sky, r, rIn, rOut)
    
    #Cumulating all the data tables from every aor
    if i==0:
        cum_data = data
    else:
        cum_data = vstack([cum_data, data])
    
    #Defining every table element except for flux terms
    aKey = header['AORKEY']
    dObs = header['DATE_OBS']
    DPat = 'YES' if 'cycl'in header['AORLABEL'] else 'NO'
    fTim = header['FRAMTIME']
    time = header['MJD_OBS']
    try:
        DScl = header['DITHSCAL']
    except:
        DScl = '--'
    try:
        DPos = str(header['DITHPOS'])
    except:
        DPos = '--'
    
    #using idl to remove systematics
    try:
        Xc       = data['Xc'][np.isnan(np.array(data['Xc']).astype('Float64')) == False]
        Yc       = data['Yc'][np.isnan(np.array(data['Yc']).astype('Float64')) == False]
        Res_Flux = data['Res_Flux'][np.isnan(np.array(data['Res_Flux']).astype('Float64')) == False]
        idl.setVariable('cenX', Xc)
        idl.setVariable('cenY', Yc)
        idl.setVariable('obsFlux', Res_Flux)
        idl.execute('corFlux = IRAC_APHOT_CORR_CRYO(obsFlux, cenX, cenY, 1)')
        corFlux = idl.getVariable('corFlux')
    except:
        print "Idl caused some trouble for AOR# %i which has %i files to work with" % ((i+1), len(data['Xc']))
        continue
    
    #applying aperture correction
    flux_bUnit = np.array(corFlux).astype('Float64')*ap_corr  #Flux in bad units (MJy/Sr)!
    
    #Fluxes in proper units (mJy)
    flux_arr = flux_bUnit*(pixArea*(10**9))
    
    #Final analysis
    flux   = np.mean(flux_arr)
    error  = np.std(flux_arr)
    spread = (error/flux)*100
    
    result.add_row([aKey, dObs, DPat, DScl, DPos, fTim, time, flux, error, spread])
    

#Writing generated data tables to csv files
ascii.write(result, 'Reduction_Data_&_Logs/run1_aor_data.csv', overwrite = False)
ascii.write(cum_data, 'Reduction_Data_&_Logs/run1_img_data.csv', overwrite = False) 



#Writing a txt file that keeps log about this reduction run
#----------------------------------------------------------
log = open('Reduction_Data_&_Logs/run1_log.txt', 'x') #The 'x' instead of 'w' prevents overwriting
log.write("Date Reduced : %s" % datetime.now().isoformat())
log.write("Instrument   : IRAC Channel 1")
log.write("File Type    : BCD")
log.write("Mission      : Cryogenic")
log.write("Target       : HD165459")
log.write("Comments     : Test run. No outliers rejected yet. Used a radius combination of (10, 12, 20)px.")
log.close()
#----------------------------------------------------------
