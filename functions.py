import numpy as np
import matplotlib.pyplot as plt
import sys, pdb, glob
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.io import fits, ascii
from astropy.table import Table, Column


def photometry(image2d, cen_x, cen_y, index = 0, shape = 'Circ', rad = None, r_in = None, r_out = None, ht = None, wid = None, w_in = None, w_out = None, h_out = None, ang = 0.0):
    """
    PARAMETERS
    ----------
    image2d = 2 dimensional image array. Type = ndarray
    cen_x, cen_y = x & y center position. Type = ndarray/list
    index = if cen_x and cen_Y is a list of more than 1 element, specify the desired index. Type = Integer
    shape = 'Circ':CircularAperture, 'Rect':RectangularAperture, 'CircAnn':CircularAnnulus, 'RectAnn':RectangularAnnulus
    rad, r_in, r_out, ht, wid, w_in, w_out, h_out, ang = Astropy's aperture parameters 
    
    RETURNS
    -------
    flux = flux of the image extracted by the aperture described in the "shape" parameter. Type = Float
    aperture = aperture object created by astropy
    """
    
    mask = np.isnan(image2d) == True
    
    if shape == 'Circ':
        aperture = CircularAperture((cen_x[index], cen_y[index]), r = rad)
    elif shape == 'Rect':
        aperture = RectangularAperture((cen_x[index], cen_y[index]), w = wid, h = ht, theta = ang)
    elif shape == 'CircAnn':
        aperture = CircularAnnulus((cen_x[index], cen_y[index]), r_in = r_in, r_out = r_out)
    elif shape == 'RectAnn':
        aperture = RectangularAnnulus((cen_x[index], cen_y[index]), w_in = w_in, w_out = w_out, h_out = h_out, theta = ang)
            
    phot_table = aperture_photometry(image2d, aperture, mask = mask)
    flux = phot_table[0][0]
    return flux, aperture
