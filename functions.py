import numpy as np
import matplotlib.pyplot as plt
import sys, pdb, glob
from tqdm import tqdm
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.io import fits, ascii
from astropy.table import Table, Column
from astropy.modeling import models, fitting
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord


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
    flux = phot_table['aperture_sum']
    return flux, aperture



def gen_center_g2d(image, center_x, center_y, box_width, amp, x_std, y_std, Theta, model_plotting = False):
    
    """
    PARAMETERS:
        center_x = x coordinate of the circular aperture; Type = float
        center_y = y coordinate of the circular aperture; Type = float
        amp = amplitude of the gaussian.  Find from the projection curve along the center; Type = float
        x_std = Standard deviation of the Gaussian in x before rotating by theta; Type = float
        y_std = Standard deviation of the Gaussian in y before rotating by theta; Type = float
        Theta = Rotation angle in radians. The rotation angle increases counterclockwise; Type = float
    
    RETURNS:
        seperate_centers =  Center of each image; Type = Array [of tuples]
        x_values = x_value of center of each image; Type = Array
        y_values = y_value of center of each image; Type = Array
    """
        
    #Creating a mesh grid with the shape of image to create model  
    y_pos, x_pos = np.mgrid[:image.shape[0],:image.shape[1]]
    
    #defining starting and stopping points for drawing a box to fit the gaussian to
    xA, yA = int(center_x-box_width), int(center_y-box_width)
    xB, yB = int(center_x+box_width), int(center_y+box_width) 
    
    # fitting the gaussian model
    fit_g = fitting.LevMarLSQFitter()
    gauss2D = models.Gaussian2D(amplitude = amp, x_mean = center_x, y_mean = center_y, x_stddev = x_std, y_stddev = y_std, theta = Theta)
    g = fit_g(gauss2D,x_pos[yA:yB,xA:xB],y_pos[yA:yB,xA:xB],image[yA:yB,xA:xB])
    g1 = fit_g(g,x_pos[yA:yB,xA:xB],y_pos[yA:yB,xA:xB],image[yA:yB,xA:xB])
    #pdb.set_trace()
    new_xCen = g1.x_mean[0]
    new_yCen = g1.y_mean[0]
    fwhm_x   = g1.x_fwhm
    fwhm_y   = g1.y_fwhm
    
    if model_plotting == True:
        
        plt.subplot(131)
        plt.imshow(image[yA:yB,xA:xB])
        plt.title('Data')
        
        plt.subplot(132)
        plt.imshow(g1(x_pos[yA:yB,xA:xB],y_pos[yA:yB,xA:xB]))
        plt.title('Model')
        
        plt.subplot(133)
        plt.imshow(image[yA:yB,xA:xB] - g1(x_pos[yA:yB,xA:xB],y_pos[yA:yB,xA:xB]))
        plt.title('Residual')
    
    #Results
    return new_xCen, new_yCen, fwhm_x, fwhm_y



def single_target_phot(fnames, targetCrd, src_r, bkg_rIn, bkg_rOut):
    """
    For a set of images.
    """
    
    data = Table(names = ('File#', 'Coord Conversion Issue', 'Centroiding Issue', 'Bad Center Guess', 'Not In FOV', 'Xc', 'Yc', 'Fx', 'Fy', 'Time[MJD]', 'Raw_Flux', 'Bkg_Flux', 'Res_Flux'), 
                 dtype = ('S25', 'S5', 'S5', 'S5', 'S5', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))


    for i, fn in tqdm(enumerate(fnames)):
        
        #Issues list
        #Initializing values to False
        (crd_conversion, centroiding, bad_cen_guess, not_in_fov) = ('X', 'X', 'X', 'X')
        
        #setting default value to NaN
        (raw_flux, bkg_flux, res_flux, cenX, cenY, fx, fy) = (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)
        
        hdu    = fits.open(fn)
        header = hdu[0].header
        image  = hdu[0].data
        hdu.close()

        Time = header['MJD_OBS']
        
        try:
            w = WCS(header)
            pix = targetCrd.to_pixel(w)
        except:
            crd_conversion = 'O'
            data.add_row([i+1, crd_conversion, centroiding, bad_cen_guess, not_in_fov, cenX, cenY, fx, fy, Time, raw_flux, bkg_flux, res_flux])
            continue

        if (pix[0]>0) & (pix[0]<256) & (pix[1]>0) & (pix[1]<256):
            
            try:
                cenX, cenY, fx, fy = gen_center_g2d(image, pix[0], pix[1], 7, 5, 4, 4, 0)
            except:
                centroiding = 'O'
                data.add_row([i+1, crd_conversion, centroiding, bad_cen_guess, not_in_fov, cenX, cenY, fx, fy, Time, raw_flux, bkg_flux, res_flux])
                continue

            if (np.abs(cenX - pix[0]) <= 2) & (np.abs(cenY-pix[1]) <= 2):
                try:
                    # Extracting raw flux
                    raw_flux, src_ap = photometry(image, [cenX], [cenY], rad = src_r)

                    # Extrating a mean background flux
                    bkg, bkg_ap = photometry(image, [cenX], [cenY], shape = 'CircAnn', r_in = bkg_rIn, r_out = bkg_rOut)
                    bkg_mean = bkg/bkg_ap.area()
                    bkg_flux = bkg_mean*src_ap.area()

                    # Subtracting background
                    res_flux  = raw_flux - bkg_flux
                    

                except:
                    continue

            else:
                bad_cen_guess = 'O'

        else:
            not_in_fov = 'O'
            
        data.add_row([i+1, crd_conversion, centroiding, bad_cen_guess, not_in_fov, cenX, cenY, fx, fy, Time, raw_flux, bkg_flux, res_flux])
        
    return data, header


def sigma_clipping(array, N):
    stdev = np.std(array)
    median  = np.median(array)
    
    upr_lim = median + (N*stdev)
    lwr_lim = median - (N*stdev)
    
    mask = ((array>lwr_lim) & (array<upr_lim))
    new_array = array[mask]
    
    return new_array, mask
    