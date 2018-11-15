import numpy as np
import functions_copy2 as func
import pdb, glob, argparse
import mirpyidl as idl
from tqdm import tqdm
from astropy.table import Table, vstack
from astropy.io import fits, ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
from datetime import datetime




#Creating a table generating function that can be called from scripts
#---------------------------------------------------------------------
def run(crdFormat, aor_crd, filetype, r, rIn, rOut, ap_corr, channel, pixLen, rejection = True):
    #initializing table to hold results
    result = Table(names = ('AORKEY', 'Target Coordinates', 'DateObs', 'Mission', 'Read Mode', 'Workable/Total Files in AOR', 'Cycling DPattern', 'DScale', 'DPosition', 'FTime (sec)', 'Time (MJD)', 'Flux (mJy)', 'Error (mJy)', 'Spread (%)', 'Outliers Rejected'), dtype = ('i4', 'S25', 'S25', 'S5', 'S5', 'S10', 'S5', 'S10', 'S5', 'f8', 'f8', 'f8', 'f8', 'f8', 'i4'))


    problem = []
    
    #Putting coordinates and filepaths in the correct format
    #........................................................
    if crdFormat.lower() == 'single_hms':
        AORs    = glob.glob(aor_crd[0])
        skyList = [SkyCoord(aor_crd[1], unit=(u.hourangle, u.deg))]*len(AORs)
        aorList = [glob.glob(aor + '/ch%i/bcd/*_%s.fits' % (channel, filetype)) for aor in AORs]
        
    elif crdFormat.lower() == 'single_deg':
        AORs    = glob.glob(aor_crd[0])
        skyList = [SkyCoord(aor_crd[1], unit=u.deg)]*len(AORs)
        aorList = [glob.glob(aor + '/ch%i/bcd/*_%s.fits' % (channel, filetype)) for aor in AORs]
        
    elif crdFormat.lower() == 'multiple_hms':
        skyTable = ascii.read(aor_crd)
        AORs     = skyTable['AOR Filepath']
        skyList  = [SkyCoord(st, unit=(u.hourangle, u.deg)) for st in skyTable['Target Coordinates']]
        aorList  = [glob.glob(aor + '/ch%i/bcd/*_%s.fits' % (channel, filetype)) for aor in AORs]
        
    elif crdFormat.lower() == 'multiple_deg':
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
            
        #Getting relevant information from data table 
        Xc       = data['Xc'][np.isnan(np.array(data['Xc']).astype('Float64')) == False]
        Yc       = data['Yc'][np.isnan(np.array(data['Yc']).astype('Float64')) == False]
        Res_Flux = data['Res_Flux'][np.isnan(np.array(data['Res_Flux']).astype('Float64')) == False]


        #Defining every table element except for flux terms
        #..................................................
        aKey = header['AORKEY']                                      #AORKEY
        dObs = header['DATE_OBS']                                    #Date of observation in UT
        dPat = 'YES' if 'cycl'in header['AORLABEL'] else 'NO'        #Dither pattern
        fTim = header['FRAMTIME']                                    #frame time in sec
        time = header['MJD_OBS']                                     #Date of observation in MJD
        mssn = 'CRYO' if (time<=54968) else 'WARM'                   #Spitzer mission
        mode = header['READMODE']                                    #Readout mode: full or sub
        fLen = '%i/%i' % (len(Res_Flux), len(aor))                   #workable file/total file
        tCrd = str(sky).replace(')','(').split('(')[-2]              #Target coordinate
        try:
            dScl = header['DITHSCAL']                                #Dither scale: small, medium or large
        except KeyError:
            dScl = '--'
        try:
            dPos = str(header['DITHPOS'])                            #Dither position
        except KeyError:
            dPos = '--'
        #...................................................
        

        #Using IDL to remove systematics
        #................................
        if len(Res_Flux) > 0:
            idl.setVariable('cenX', Xc)
            idl.setVariable('cenY', Yc)
            idl.setVariable('obsFlux', Res_Flux)
            idl.setVariable('ch', channel)
            
            if mode == 'SUB':
                idl.execute('corFlux = PIXEL_PHASE_CORRECT_GAUSS(obsFlux, cenX, cenY, ch)')
            else:
                if mssn =='CRYO':
                    idl.execute('corFlux = IRAC_APHOT_CORR(obsFlux, cenX, cenY, ch, /CRYO)')
                else:
                    idl.execute('corFlux = IRAC_APHOT_CORR(obsFlux, cenX, cenY, ch)')
            corFlux = idl.getVariable('corFlux')
            corFlux = np.array([corFlux]) if type(corFlux)==float else np.array(corFlux).astype('Float64')
            pdb.set_trace()
        else:
            problem.append(aKey)
            result.add_row([aKey, tCrd, dObs, mssn, mode, fLen, dPat, dScl, dPos, fTim, time, np.nan, np.nan, np.nan, 0])
            continue
        #.................................
        


        #Rejecting Outliers
        #...................
        sig = (np.std(corFlux)/np.std(corFlux))*100
        if (rejection == False) | (len(corFlux)<=10) | (sig<2):
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

        
        result.add_row([aKey, tCrd, dObs, mssn, mode, fLen, dPat, dScl, dPos, fTim, time, flux, error, spread, clipped])
    
    return result, cum_data, problem

#--------------------------------------------------------------------





# Defining a function that parses arguments given in terminal
#-------------------------------------------------------------

def parse_arguments():
    parser = argparse.ArgumentParser()
    
    #crdFormat
    parser.add_argument('crdFormat', nargs = '?', type = str, default = 'single_hms', help = 'A string that specifies format that specifies the format of target coordinate. It tells the unit of your coordinate and whether you want to work with one or more than one corrdinate. It can take in single_hms, single_deg, multiple_hms or multiple_hms as values. The values are not case sensitive.')
    
    #aor_crd
    parser.add_argument('aor_crd', nargs = '+', type = str, help = 'AOR filepaths and coordinates. It can be individual file path with a single coordinate or a text file containing those information. Values must be strings. See github documentation for more on formatting and examples')
    
    #filetype
    parser.add_argument('-f', '--file-type', dest = 'filetype', nargs = '?', const = 'bcd', default = 'bcd', help = 'Which file type to use - bcd, cbcd, sub2d etc.')
    
    #radius combination
    parser.add_argument('-r', '--radius', nargs = 3, default = [10, 12, 20], type = int, help = 'Radius combination. Source aperture radius, inner and outer background aperture radius.')
    
    #ap_corr
    parser.add_argument('-a', '--aperture-correction', dest = 'ap_corr', nargs = '?', const = 1.000, default = 1.000, type = float, help = 'Aperture correction for a particular radius combination. Proper values can be found from github documenation. Defaults to 1.0 when not specified')
    
    #channel + pixLen
    parser.add_argument('-cp', '--ch-px', dest = 'ch_px', nargs = 2, default = [1, 1.221], type = float, help = 'Two arguments. Channel number and length of a pixel for that channel. defaults to channel 1 values.')
    
    #Outlier Rejection
    parser.add_argument('-o', '--outlier-rejection', dest = 'outlier', action = 'store_true', help = 'Whether you want outlier rejection or not. Just the flag, no argument required.')
    
    #Run Name
    parser.add_argument('-rn', '--run-name', dest = 'run', nargs = '?', const = '100', default = '100', help = 'Run number for naming output files')
    
    #Target Name
    parser.add_argument('-t', '--target-name', dest = 'target', nargs = '?', const = 'Many', default = 'Many', help = 'Name of target or targets. This will be included in the log file.')
    
    #Comments
    parser.add_argument('-c', '--comments', nargs = '?', type = str, default = 'No comment was specified.', help = 'Comments about the run that will be included in the log file')
    
    args = parser.parse_args()
    
    return args

#-------------------------------------------------------------





if __name__=='__main__':
    
    
    #Defining general constants
    #---------------------------
    args = parse_arguments()
    crdFormat = args.crdFormat
    aor_crd = args.aor_crd[0] if len(args.aor_crd)==1 else args.aor_crd
    filetype = args.filetype
    r, rIn, rOut = args.radius
    ap_corr = args.ap_corr
    ch, pix = args.ch_px
    rejection = args.outlier
    nRun = args.run
    tName = args.target
    comments = args.comments
    

    #----------------------------
    
    #Generating & writing data tables to csv files
    res, data, prob = run(crdFormat, aor_crd, filetype, r, rIn, rOut, ap_corr, int(ch), pix, rejection)
    ascii.write(res, '../Reduction_Data_&_Logs/%s_aor_data.csv' % nRun, delimiter = ',', overwrite = True)
    ascii.write(data, '../Reduction_Data_&_Logs/%s_img_data.csv' % nRun, delimiter = ',', overwrite = True) 
    



    #Writing a txt file that keeps log about this reduction run
    #----------------------------------------------------------
    log = open('../Reduction_Data_&_Logs/%s_log.txt' % nRun, 'w')
    log.write("Date Reduced     : %s \n" % datetime.now().isoformat())
    log.write("Input Parameters : %s \n" % str(vars(args)))
    log.write("Instrument       : IRAC Channel %i \n" % int(ch))
    log.write("File Type        : %s \n" % filetype.upper())
    log.write("Target           : %s \n" % tName)
    log.write("Radius Used      : r %i, rIn %i, rOut %i \n" % (r, rIn, rOut))
    log.write("Problem AORs     : " + ("%s "*len(prob) % tuple(prob)) + "\n")
    log.write("Comments         : %s" % comments)
    log.close()
    #----------------------------------------------------------
