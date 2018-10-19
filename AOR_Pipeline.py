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
def run(crdFormat, aor_crd, channel, filetype, r, rIn, rOut, ap_corr, pixLen, N):
    #initializing table to hold results
    result = Table(names = ('AORKEY', 'Target Coordinates', 'DateObs', 'Mission', 'Read Mode', 'Cycling DPattern', 'DScale', 'DPosition', 'FTime (sec)', 'Time (MJD)', 'Flux (mJy)', 'Error (mJy)', 'Spread (%)', 'Outliers Rejected'), dtype = ('i4', 'S25', 'S25', 'S5', 'S5', 'S5', 'S10', 'S5', 'f8', 'f8', 'f8', 'f8', 'f8', 'i4'))


    problem = []
    
    #Putting coordinates and filepaths in the correct format
    #........................................................
    if crdFormat.lower() == 'single hms':
        AORs    = glob.glob(aor_crd[0])
        skyList = [SkyCoord(aor_crd[1], unit=(u.hourangle, u.deg))]*len(AORs)
        aorList = [glob.glob(aor + '/ch%i/bcd/*_%s.fits' % (channel, filetype)) for aor in AORs]
        
    elif crdFormat.lower() == 'single deg':
        AORs    = glob.glob(aor_crd[0])
        skyList = [SkyCoord(aor_crd[1], unit=u.deg)]*len(AORs)
        aorList = [glob.glob(aor + '/ch%i/bcd/*_%s.fits' % (channel, filetype)) for aor in AORs]
        
    elif crdFormat.lower() == 'multiple hms':
        skyTable = ascii.read(aor_crd)
        AORs     = skyTable['AOR Filepath']
        skyList  = [SkyCoord(st, unit=(u.hourangle, u.deg)) for st in skyTable['Target Coordinates']]
        aorList  = [glob.glob(aor + '/ch%i/bcd/*_%s.fits' % (channel, filetype)) for aor in AORs]
        
    elif crdFormat.lower() == 'multiple deg':
        skyTable = ascii.read(aor_crd)
        AORs     = skyTable['AOR Filepath']
        skyList  = [SkyCoord(st, unit=u.deg) for st in skyTable['Target Coordinates']]
        aorList  = [glob.glob(aor + '/ch%i/bcd/*_%s.fits' % (channel, filetype)) for aor in AORs]
        
    else:
        print('Please input one of these 4 values for crdFormat: single hms, single deg, multiple hms, multiple deg')
        return
    #..........................................................
    
    counter = 0
    #Generate data for each aor
    for i, (aor, sky) in tqdm(enumerate(zip(aorList, skyList))):
        
        #Performing photometry on every image in an AOR=
        if len(aor)>0:
            data, header = func.single_target_phot(aor, sky, r, rIn, rOut)
            counter += 1
        else:
            problem.append(AORs[i].split('/r')[1])
            print('Empty AOR encountered: %s' % AORs[i].split('/r')[1])
            continue


        #Cumulating all the data tables from every aor
        if counter==1:
            cum_data = data
        else:
            cum_data = vstack([cum_data, data])


        #Defining every table element except for flux terms
        #..................................................
        aKey = header['AORKEY']
        dObs = header['DATE_OBS']
        dPat = 'YES' if 'cycl'in header['AORLABEL'] else 'NO'
        fTim = header['FRAMTIME']
        time = header['MJD_OBS']
        mssn = 'CRYO' if (time<=54968) else 'WARM'
        mode = header['READMODE']
        tCrd = str(sky).replace(')','(').split('(')[-2] 
        try:
            dScl = header['DITHSCAL']
        except:
            dScl = '--'
        try:
            dPos = str(header['DITHPOS'])
        except:
            dPos = '--'
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
            idl.setVariable('ch', channel)
            
            if mode == 'SUB':
                if mssn == 'CRYO':
                    idl.execute('corFlux = IRAC_APHOT_CORR(obsFlux, cenX, cenY, ch, /CRYO, /SUBARRAY)')
                else:
                    idl.execute('corFlux = IRAC_APHOT_CORR(obsFlux, cenX, cenY, ch, /SUBARRAY)')
            else:
                if mssn =='CRYO':
                    idl.execute('corFlux = IRAC_APHOT_CORR(obsFlux, cenX, cenY, ch, /CRYO)')
                else:
                    idl.execute('corFlux = IRAC_APHOT_CORR(obsFlux, cenX, cenY, ch)')
            corFlux = idl.getVariable('corFlux')
            corFlux = np.array([corFlux]) if type(corFlux)==float else np.array(corFlux).astype('Float64')
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
        pixArea = (pixLen**2)/(206265**2) #Steradian
        flux_arr = flux_bUnit*(pixArea*(10**9))

        #Analysis
        flux   = np.mean(flux_arr)
        error  = np.std(flux_arr)
        spread = (error/flux)*100
        #........................

        
        result.add_row([aKey, tCrd, dObs, mssn, mode, dPat, dScl, dPos, fTim, time, flux, error, spread, clipped])
    
    return result, cum_data, problem

#--------------------------------------------------------------------

    

if __name__=='__main__':
    
    #Defining general constants
    #---------------------------

    print '\n', 'Please input the following arguments and seperate them with a semi-colon(;)-- \n', '1. Coordinate format (e.g. single hms, single deg, multiple hms, multiple deg) \n', '2. AOR file path and coordinates (e.g. /data1/phot_cal/spitzer/hd165459/cryo/r*,18 02 30.7410086899 +58 37 38.157415821 \n', '3. File type (eg. bcd, cbcd, sub2d) \n', '4. Target Name (e.g. HD 165459, BD + 60 1753) \n','5. Source Aperture Radius in px (eg. 10) \n', '6. Inner Background Radius in px (eg. 12) \n', '7. Outer Background Radius in px (eg. 20) \n', '8. Channel# (1,2,3 or 4) \n', '9. Aperture Correction Factor (collect from iracinstrumenthandbook/27, e.g. 1.000) \n', '10. Length of a pixel in arcsec (e.g. 1.221 for ch1) \n', '11. Run# (For output file name. Should be an integer)\n', '12. Sigma Clipping Number (e.g. 10) \n', '13. Comments (to be included in the log) \n', 'Don\'t jump to any argument. e.g. You can\'t skip sigma to get to comments. comments must be the 14th argument. \n', 'You could however, use \'n/a\' for sigma if you don\'t want sigma clipping \n', 'Example command: single hms;/data1/phot_cal/spitzer/hd165459/cryo/r11638016,18 02 30.7410086899 +58 37 38.157415821;bcd;HD 165459;10;12;20;1;1.0;1.221;14;10;Your comment goes here. \n', '\n'

    const_str = raw_input('Input the parameters listed above: ') 
    constants = const_str.split(';')
    print constants

    crdFormat    = constants[0] 
    aor_crd      = constants[1] if 'csv' in constants[1] else constants[1].split(',') 
    filetype     = constants[2].lower()
    r, rIn, rOut = int(constants[4]), int(constants[5]), int(constants[6])
    channel      = int(constants[7])
    ap_corr      = float(constants[8])
    pixLen       = float(constants[9]) #arcsec
    N            = constants[11] #For outlier rejection

    #----------------------------
    
    #Generating & writing data tables to csv files
    res, data, prob = run(crdFormat, aor_crd, channel, filetype, r, rIn, rOut, ap_corr, pixLen, N)
    ascii.write(res, 'Reduction_Data_&_Logs/run%s_aor_data.csv' % constants[10], delimiter = ',', overwrite = True)
    ascii.write(data, 'Reduction_Data_&_Logs/run%s_img_data.csv' % constants[10], delimiter = ',', overwrite = True) 



    #Writing a txt file that keeps log about this reduction run
    #----------------------------------------------------------
    log = open('Reduction_Data_&_Logs/run%s_log.txt' % constants[10], 'w')
    log.write("Date Reduced : %s \n" % datetime.now().isoformat())
    log.write("Command Used : %s \n" % const_str)
    log.write("Instrument   : IRAC Channel %i \n" % channel)
    log.write("File Type    : %s \n" % filetype.upper())
    log.write("Target       : %s \n" % constants[3])
    log.write("Radius Used  : src_r-%i, bkg_in-%i, bkg_out-%i \n" % (r, rIn, rOut))
    log.write("Problem AORs : " + ("%s "*len(prob) % tuple(prob)) + "\n")
    log.write("Comments     : %s" % constants[12])
    log.close()
    #----------------------------------------------------------
