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


#Defining general constants
#---------------------------

print '\n', 'Please input the following arguments and seperate them with a comma-- \n', '1. File Path up to aor names (e.g. /data1/phot_cal/spitzer/hd165459/cryo/r*/) \n', '2. File type (eg. bcd) \n', '3. Target Coordinates (e.g. 18 02 30.7410086899 +58 37 38.157415821) \n', '4. Source Aperture Radius in px (eg. 10) \n', '5. Inner Background Radius in px (eg. 12) \n', '6. Outer Background Radius in px (eg. 20) \n', '7. Aperture Correction Factor (collect from iracinstrumenthandbook/27, e.g. 1.000) \n', '8. Length of a pixel in arcsec (e.g. 1.221 for ch1) \n', '9. Sigma Clipping Number (e.g. 10) [optional argument]', '\n'

constants = raw_input('Input the parameters listed above: ').split(',')

                      
AORnames = np.sort(glob.glob(constants[0]))
AORs = []

for aor in AORnames:
    fnames = np.sort(glob.glob(aor + 'ch1/' + constants[1] + '/*_' + constants[1] + '.fits'))
    AORs.append(fnames)
                      
                      
#Provide proper sky coordinates in hms
# '18 02 30.7410086899 +58 37 38.157415821' for HD 165459
# '17 24 52.2772360943 +60 25 50.780790994' for BD +60 1753
sky = SkyCoord(constants[2], unit=(u.hourangle, u.deg))

r, rIn, rOut = int(constants[3]), int(constants[4]), int(constants[5])

# Find the aperture correction factor from the following link:
# https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/27/
# From the table at the end of that link, select your value accoring to your aperture size and channel
ap_corr = float(constants[6])

pixLen  = float(constants[7]) #arcsec
pixArea = pixLen**2 #arcsec^2
pixArea = pixArea/(206265**2) #Steradian

#----------------------------


#initializing table to hold results
result = Table(names = ('AORKEY', 'DateObs', 'Cycling DPattern', 'DScale', 'DPosition', 'FTime (sec)', 'Time (MJD)', 'Flux (mJy)', 'Error (mJy)', 'Spread (%)'), dtype = ('i4', 'S25', 'S5', 'S10', 'S5', 'f8', 'f8', 'f8', 'f8', 'f8'))

#Generate data for each aor 
for i, aor in tqdm(enumerate(AORs)):
    
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
    Xc       = data['Xc'][np.isnan(np.array(data['Xc']).astype('Float64')) == False]
    Yc       = data['Yc'][np.isnan(np.array(data['Yc']).astype('Float64')) == False]
    Res_Flux = data['Res_Flux'][np.isnan(np.array(data['Res_Flux']).astype('Float64')) == False]
    if len(Res_Flux) > 0:
        idl.setVariable('cenX', Xc)
        idl.setVariable('cenY', Yc)
        idl.setVariable('obsFlux', Res_Flux)
        idl.execute('corFlux = IRAC_APHOT_CORR_CRYO(obsFlux, cenX, cenY, 1)')
        corFlux = idl.getVariable('corFlux')
    else:
        print "Idl caused some trouble for AOR# %i, %i which has %i centers and %i fluxes" % ((i+1), aKey, len(Xc), len(Res_Flux))
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
ascii.write(result, 'Reduction_Data_&_Logs/run1_aor_data.csv', delimiter = ',')
ascii.write(cum_data, 'Reduction_Data_&_Logs/run1_img_data.csv', delimiter = ',') 



#Writing a txt file that keeps log about this reduction run
#----------------------------------------------------------
log = open('Reduction_Data_&_Logs/run1_log.txt', 'w') 
log.write("Date Reduced : %s \n" % datetime.now().isoformat())
log.write("Instrument   : IRAC Channel 1 \n")
log.write("File Type    : %s \n" % constants[1].upper())
log.write("Mission      : Cryogenic \n")
log.write("Target       : HD165459 \n")
log.write("Average Flux : %.2f +- %.2f \n" % (np.mean(result['Flux (mJy)']), np.std(result['Flux (mJy)'])))
log.write("Comments     : Test run. No outliers rejected yet. Used a radius combination of (10, 12, 20)px.")
log.close()
#----------------------------------------------------------
