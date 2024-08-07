#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.skyview import SkyView
from astropy.io import fits, ascii
from astropy.wcs import WCS
import astropy.units as u

import sys

Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"

EXAMPLES = """

search_gaia_cm.py -i input_table.fits -r 1.0 -o output_table.fits

           """

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='''

	Get a list with coordinates and name and run a query on gaia catalog.

        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('-i', '--input', required=True, type=str,
                        help='File containing\
                            the target coordinates and name\
		            in columns called\
                            ra, dec, and name\
			    it is a fits file')

    parser.add_argument('-r', '--radius', required=True, type=float,
                        help='Radius for the search of counterpart\
                            it is a float in arcsec')

    parser.add_argument('-o', '--output', required=True, type=str,
                        help='File output name\
			    it is a fits file')


    return parser.parse_args()

def find_gaia(coord, rad):   
    """
    Find the closest Gaia entry to the position.

    Parameters:
        - coord: astropy.SkyCoord
            Coordinate of the target. The coordinates are automatically
            tranformed to the ICRS system to match the Gaia catalogue
         - radius: astropy.units
            radius for the search. 
    """
    
    coord = coord.transform_to('icrs')

    radius = rad*u.arcsec
    
    # search sources in GAIA
    j = Gaia.cone_search(coord, radius)
    cat = j.get_results()
    #print(len(cat))
    if len(cat) < 0.5:
       print("No counterpart in GAIA DR3")
       return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    if len(cat) > 0:
       selected = cat[0]
       return selected['ra'], selected['dec'], selected['source_id'], selected['phot_g_mean_flux'], selected['phot_g_mean_flux_error'], selected['phot_g_mean_mag']     

if __name__ == '__main__':

    args = parse_arguments()

    table = fits.getdata(args.input)
    ra_s = np.array(table['ra'])
    dec_s = np.array(table['dec'])
    name_s = table['quick_object_id']

    rad = args.radius

    coord_deg = SkyCoord(ra=ra_s*u.degree, dec=dec_s*u.degree)
    
    ra_g_all = []         ;  dec_g_all = []
    source_id_g_all = []  ;  g_mag_g_all = []
    g_flux_g_all = []     ;  g_flux_err_g_all = []
    
    for i in np.arange(len(table)):
    #for i in np.arange(2):
        print("Name source:",name_s[i])
        coord = SkyCoord(coord_deg[i].ra.hms, coord_deg[i].dec.dms, unit=(u.hourangle, u.deg),
                        frame='icrs')
        #stars = find_gaia(coord, 13, 19, 3, 3,plot=False)
        ra_g, dec_g, source_id_g, g_flux_g, g_flux_err_g, g_mag_g = find_gaia(coord, rad)
        ra_g_all.append(ra_g)
        dec_g_all.append(dec_g)
        source_id_g_all.append(source_id_g)
        g_flux_g_all.append(g_flux_g)
        g_flux_err_g_all.append(g_flux_err_g)
        g_mag_g_all.append(g_mag_g)
    
    #---Create output table
    namec   = fits.Column(name='source_id', array=source_id_g_all, format='A12')
    rac     = fits.Column(name='ra_gaia_dr3', array=ra_g_all, format='F')
    decc    = fits.Column(name='dec_gaia_dr3', array=dec_g_all, format='F')
    fluxc   = fits.Column(name='phot_g_mean_flux_gaia', array=g_flux_g_all, format='F')
    fluxec  = fits.Column(name='phot_g_mean_flux_error_gaia', array=g_flux_err_g_all, format='F')
    magc    = fits.Column(name='phot_g_mean_mag_gaia', array=g_mag_g_all, format='F')
    t_out = fits.BinTableHDU.from_columns([namec,rac,decc,fluxc,fluxec,magc])
    name_t_out = args.output
    t_out.writeto(name_t_out,overwrite=True)
    
    print("==============================================")    
    #print(ra_g_all,dec_g_all,source_id_g_all,g_flux_g_all,g_flux_err_g_all,g_mag_g_all)
    print("Output file:",name_t_out)       
    
    
         
