import numpy as np
import functions as func
import pdb, glob
import mirpyidl as idl
from tqdm import tqdm
from astropy.table import Table
from astropy.io import fits, ascii
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord


AORnames = np.sort(glob.glob('/data1/phot_cal/spitzer/hd165459/cryo/r*/'))
AORs = []

for aor in AORnames:
    fnames = np.sort(glob.glob(aor + 'ch1/bcd/*_bcd.fits'))
    AORs.append(fnames)

#Provide proper sky coordinates in hms
# '18 02 30.7410086899 +58 37 38.157415821' for HD 165459
# '17 24 52.2772360943 +60 25 50.780790994' for BD +60 1753

sky = SkyCoord('18 02 30.7410086899 +58 37 38.157415821', unit=(u.hourangle, u.deg))

#Generate data for each aor 
for i, aor in enumerate(AORs):
    
    result = Table(names = ('AORKEY', 'DateObs', 'Cycling DPattern', 'DScale', 'DPosition', 'FTime (sec)', 'Time (MJD)', 'Flux (mJy)', 'Error (mJy)', 'Spread (%)'), dtype = ('i4', 'S25', 'S5', 'S10', 'S5', 'f8', 'f8', 'f8', 'f8', 'f8'))
    
    
    data, header = func.single_target_phot(aor, sky, 10, 12, 20)
    
    aKey = header['AORKEY']
    dObs = header['DATE_OBS']
    DPat = 'YES' if 'cycl'in header['AORLABEL'] else 'NO'
    fTim = header['FRAMTIME']
    Time = header['MJD_OBS']
    try:
        DScl = header['DITHSCAL']
    except:
        DScl = '--'
    try:
        DPos = str(header['DITHPOS'])
    except:
        DPos = '--'
    
    Xc       = data['Xc'][np.isnan(np.array(data['Xc']).astype('Float64')) == False]
    Yc       = data['Yc'][np.isnan(np.array(data['Yc']).astype('Float64')) == False]
    Res_Flux = data['Res_Flux'][np.isnan(np.array(data['Res_Flux']).astype('Float64')) == False]
    idl.setVariable('cenX', Xc)
    idl.setVariable('cenY', Yc)
    idl.setVariable('obsFlux', Res_Flux)
    idl.execute('corFlux = IRAC_APHOT_CORR_CRYO(obsFlux, cenX, cenY, 1)')
    corFlux = idl.getVariable('corFlux')
    pdb.set_trace()
    
    if i>0:
        data = vstack([data, data])
    