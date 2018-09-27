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




#Creating a table generating function that can be called from scripts
#---------------------------------------------------------------------
def run(AORs, sky, r, rIn, rOut, ap_corr, pixArea, N):
    #initializing table to hold results
    result = Table(names = ('AORKEY', 'DateObs', 'Cycling DPattern', 'DScale', 'DPosition', 'FTime (sec)', 'Time (MJD)', 'Flux (mJy)', 'Error (mJy)', 'Spread (%)', 'Outliers Rejected'), dtype = ('i4', 'S25', 'S5', 'S10', 'S5', 'f8', 'f8', 'f8', 'f8', 'f8', 'i4'))


    problem = []

    #Generate data for each aor
    for i, aor in tqdm(enumerate(AORs)):


        #Performing photometry on every image in an AOR
        if len(aor)>0:
            data, header = func.single_target_phot(aor, sky, r, rIn, rOut)
        else:
            problem.append(header['AORKEY'])


        #Cumulating all the data tables from every aor
        if i==0:
            cum_data = data
        else:
            cum_data = vstack([cum_data, data])


        #Defining every table element except for flux terms
        #..................................................
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
        #...................................................


        #Using IDL to remove systematics
        #................................
        Xc       = data['Xc'][np.isnan(np.array(data['Xc']).astype('Float64')) == False]
        Yc       = data['Yc'][np.isnan(np.array(data['Yc']).astype('Float64')) == False]
        Res_Flux = data['Res_Flux'][np.isnan(np.array(data['Res_Flux']).astype('Float64')) == False]
        if len(Res_Flux) > 0:
            idl.setVariable('cenX', Xc)
            idl.setVariable('cenY', Yc)
            idl.setVariable('obsFlux', Res_Flux)
            idl.setVariable('ch', constants[8])
            if constants[4] == 'Warm':
                idl.execute('corFlux = IRAC_APHOT_CORR(obsFlux, cenX, cenY, ch)')
            else:
                idl.execute('corFlux = IRAC_APHOT_CORR_CRYO(obsFlux, cenX, cenY, ch)')
            corFlux = idl.getVariable('corFlux')
            corFlux = np.array(corFlux).astype('Float64')
        else:
            problem.append(aKey)
            continue
        #.................................


        #Rejecting Outliers
        #...................
        sig = (np.std(corFlux)/np.std(corFlux))*100
        if (N == 'n/a') | (len(corFlux)<=10) | (sig<2):
            refFlux = corFlux
        else:
            refFlux = np.delete(corFlux, [np.argmin(corFlux), np.argmax(corFlux)])
        clipped = len(corFlux) - len(refFlux)
        #...................


        #Final Flux calculations
        #........................

        #applying aperture correction
        flux_bUnit = refFlux*ap_corr  #Flux in bad units (MJy/Sr)!

        #Fluxes in proper units (mJy)
        flux_arr = flux_bUnit*(pixArea*(10**9))

        #Analysis
        flux   = np.mean(flux_arr)
        error  = np.std(flux_arr)
        spread = (error/flux)*100
        #........................

        
        result.add_row([aKey, dObs, DPat, DScl, DPos, fTim, time, flux, error, spread, clipped])
    
    return result, cum_data, problem

#--------------------------------------------------------------------

    

if __name__=='__main__':
    
    #Defining general constants
    #---------------------------

    print '\n', 'Please input the following arguments and seperate them with a semi-colon(;)-- \n', '1. File Path up to aor names (e.g. /data1/phot_cal/spitzer/hd165459/cryo/r*/) \n', '2. File type (eg. bcd) \n', '3. Target Name (e.g. HD 165459)\n', '4. Target Coordinates (e.g. 18 02 30.7410086899 +58 37 38.157415821) \n', '5. Mission (e.g. Cryogenic or Warm)\n','6. Source Aperture Radius in px (eg. 10) \n', '7. Inner Background Radius in px (eg. 12) \n', '8. Outer Background Radius in px (eg. 20) \n', '9. Channel# (1,2,3 or 4) \n', '10. Aperture Correction Factor (collect from iracinstrumenthandbook/27, e.g. 1.000) \n', '11. Length of a pixel in arcsec (e.g. 1.221 for ch1) \n', '12. Run# (For output file name. Should be an integer)\n', '13. Sigma Clipping Number (e.g. 10) \n', '14. Comments (to be included in the log) \n', 'Don\'t jump to any argument. e.g. You can\'t skip sigma to get to comments. comments must be the 14th argument. \n', 'You could however, use \'n/a\' for sigma if you don\'t want sigma clipping \n', 'Example command: /data1/phot_cal/spitzer/hd165459/cryo/r*/;bcd;HD 165459;18 02 30.7410086899 +58 37 38.157415821;Cryogenic;10;12;20;1;1.0;1.221;5;10;Your comment goes here. \n', '\n'

    const_str = raw_input('Input the parameters listed above: ') 
    constants = const_str.split(';')
    print constants


    AORnames = np.sort(glob.glob(constants[0]))
    AORs = []

    for aor in AORnames:
        fnames = np.sort(glob.glob(aor + 'ch1/' + constants[1] + '/*_' + constants[1] + '.fits'))
        AORs.append(fnames)


    #Provide proper sky coordinates in hms
    # '18 02 30.7410086899 +58 37 38.157415821' for HD 165459
    # '17 24 52.2772360943 +60 25 50.780790994' for BD +60 1753
    sky = SkyCoord(constants[3], unit=(u.hourangle, u.deg))

    r, rIn, rOut = int(constants[5]), int(constants[6]), int(constants[7])

    # Find the aperture correction factor from the following link:
    # https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/27/
    # From the table at the end of that link, select your value accoring to your aperture size and channel
    ap_corr = float(constants[9])

    pixLen  = float(constants[10]) #arcsec
    pixArea = pixLen**2 #arcsec^2
    pixArea = pixArea/(206265**2) #Steradian

    N = constants[12] #For outlier rejection

    #----------------------------
    
    
    #Generating & writing data tables to csv files
    res, data, prob = run(AORs, sky, r, rIn, rOut, ap_corr, pixArea, N)
    ascii.write(res, 'Reduction_Data_&_Logs/run%s_aor_data.csv' % constants[11], delimiter = ',', overwrite = True)
    ascii.write(data, 'Reduction_Data_&_Logs/run%s_img_data.csv' % constants[11], delimiter = ',', overwrite = True) 



    #Writing a txt file that keeps log about this reduction run
    #----------------------------------------------------------
    log = open('Reduction_Data_&_Logs/run%s_log.txt' % constants[11], 'w')
    log.write("Date Reduced : %s \n" % datetime.now().isoformat())
    log.write("Command Used : %s \n" % const_str)
    log.write("Instrument   : IRAC Channel %s \n" % constants[8])
    log.write("File Type    : %s \n" % constants[1].upper())
    log.write("Mission      : %s \n" % constants[4])
    log.write("Target       : %s \n" % constants[2])
    log.write("Radius Used  : src_r: %i, bkg_in: %i, bkg_out: %i \n" % (r, rIn, rOut))
    log.write("Average Flux : %.2f +- %.2f \n" % (np.mean(res['Flux (mJy)']), np.std(res['Flux (mJy)'])))
    log.write("Problem AORs : " + ("%i "*len(prob) % tuple(prob)) + "\n")
    log.write("Comments     : %s" % constants[13])
    log.close()
    #----------------------------------------------------------
