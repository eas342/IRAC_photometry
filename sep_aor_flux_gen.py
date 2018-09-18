import numpy as np
import functions as func
import pdb, glob
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

print len(AORs), len(AORs[0])    

#Provide proper sky coordinates in hms
# '18 02 30.7410086899 +58 37 38.157415821' for HD 165459
# '17 24 52.2772360943 +60 25 50.780790994' for BD +60 1753

sky = SkyCoord('18 02 30.7410086899 +58 37 38.157415821', unit=(u.hourangle, u.deg))

data, hdr = func.single_target_phot(AORs[0], sky, 10, 12, 20)

print data, hdr