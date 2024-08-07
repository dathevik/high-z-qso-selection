
#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd

from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord

from dl import queryClient as qc

EXAMPLES = """

query_table_delve.py -i input_table.fits

           """


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='''

	Get a list with coordinates and name and run a query on delve_dr2.objects catalog.

        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('-i', '--input', required=True, type=str,
                        help='File containing\
                            the target coordinates and name\
		            in columns called\
                            ra, dec, and name\
			    it is a fits file')

    return parser.parse_args()


def search_survey(ra_s_min, ra_s_max, dec_s_min, dec_s_max, survey_name, name_s,  decals_dr10_important):
    print("Checking the " + str(survey_name) + "table")

    all_results = []

    for ii in range(len(ra_s_min)):
        print("======= Searching", name_s[ii])

        qq = 'SELECT ' + decals_dr10_important
        qq += '  FROM '
        qq += survey_name
        qq += ' WHERE  ra > ('
        qq += str(ra_s_min[ii])
        qq += ') AND ra < ('
        qq += str(ra_s_max[ii])
        qq += ') AND dec > ('
        qq += str(dec_s_min[ii])
        qq += ') AND dec <('
        qq += str(dec_s_max[ii])
        qq += ')'
        print("Query:", qq)
        response_delve_dr2 = qc.query(sql=qq, fmt='pandas')
        response_delve_dr2['query_name'] = name_s[ii]  # Add a column to identify the query
        all_results.append(response_delve_dr2)

    combined_results = pd.concat(all_results, ignore_index=True)
    combined_results.to_csv('racs_sample_decals.csv', index=False)


if __name__ == '__main__':
    args = parse_arguments()

    table = fits.getdata(args.input)
    ra_s = np.array(table['ra'])
    dec_s = np.array(table['dec'])
    name_s = table['ls_id']

    c_s = SkyCoord(ra=ra_s, dec=dec_s, unit=(u.deg, u.deg))
    ra_s_deg = c_s.ra.deg
    dec_s_deg = c_s.dec.deg
    print(ra_s_deg, dec_s_deg)

    ra_s_min = ra_s_deg - 1. / 3600.
    ra_s_max = ra_s_deg + 1. / 3600.
    dec_s_min = dec_s_deg - 1. / 3600.
    dec_s_max = dec_s_deg + 1. / 3600.

    survey_name_decals = 'ls_dr10.tractor'

    # -----------------------------------------DELVE Survey table----------------------------
    # ---coordinates
    # des_dr2_important = survey_name_des + '.ra as ra_des_dr2,' + survey_name_des + '.dec as dec_des_dr2, '
    decals_dr10_important = survey_name_decals + '.ra as ra_decals,' + survey_name_decals + '.dec as dec_decals, '
    # delve_dr2_important = survey_name_delve + '.ra as ra_delve,' + survey_name_delve + '.dec as dec_delve, '
    # gaia_dr3_important = survey_name_gaia+'.ra as ra_gaia_dr3,'+survey_name_gaia+'.dec as dec_gaia_dr3, '
    # vhs_dr5_important = survey_name_vhs+'.ra2000 as ra2000,'+survey_name_vhs+'.dec2000 as dec2000, '
    # allwise_important = survey_name_allwise +'.ra as ra,'+survey_name_allwise+'.dec as dec, '
    # ----magnitudes
    # --mag_auto
    # delve_dr2_important += survey_name_delve + '.mag_auto_g as mag_auto_g, '
    # delve_dr2_important += survey_name_delve + '.magerr_auto_g as magerr_auto_g, '
    # delve_dr2_important += survey_name_delve + '.extended_class_z as extended_class_z '
    # des_dr2_important += survey_name_des + '.mag_auto_g as mag_auto_g, '
    # des_dr2_important += survey_name_des + '.magerr_auto_g as magerr_auto_g, '
    # des_dr2_important += survey_name_des + '.mag_auto_z as mag_auto_z, '
    # des_dr2_important += survey_name_des + '.magerr_auto_z as magerr_auto_z '
    # gaia_dr3_important += survey_name_gaia + '.source_id as source_id, '
    decals_dr10_important += survey_name_decals + '.ls_id as ls_id, '
    decals_dr10_important += survey_name_decals + '.mag_g as mag_g, '
    decals_dr10_important += survey_name_decals + '.snr_g as snr_g, '
    decals_dr10_important += survey_name_decals + '.mag_r as mag_r, '
    decals_dr10_important += survey_name_decals + '.snr_r as snr_r, '
    decals_dr10_important += survey_name_decals + '.mag_i as mag_i, '
    decals_dr10_important += survey_name_decals + '.snr_i as snr_i, '
    decals_dr10_important += survey_name_decals + '.mag_z as mag_z, '
    decals_dr10_important += survey_name_decals + '.snr_z as snr_z '
    # vhs_dr5_important += survey_name_vhs + '.ypetromag as ypetromag, '
    # vhs_dr5_important += survey_name_vhs + '.ypetromagerr as ypetromagerr, '
    # vhs_dr5_important += survey_name_vhs + '.jpetromag as jpetromag, '
    # vhs_dr5_important += survey_name_vhs + '.jpetromagerr as jpetromagerr, '
    # vhs_dr5_important += survey_name_vhs + '.hpetromag as hpetromag, '
    # vhs_dr5_important += survey_name_vhs + '.hpetromagerr as hpetromagerr, '
    # vhs_dr5_important += survey_name_vhs + '.kspetromag as kspetromag, '
    # vhs_dr5_important += survey_name_vhs + '.kspetromagerr as kspetromagerr '
    # allwise_important += survey_name_allwise + '.w1mpro as w1mpro, '
    # allwise_important += survey_name_allwise + '.w1sigmpro as w1sigmpro, '
    # allwise_important += survey_name_allwise + '.w2mpro as w2mpro, '
    # allwise_important += survey_name_allwise + '.w2sigmpro as w2sigmpro '
    # gaia_dr3_important += survey_name_gaia + '.pm as pm, '
    # gaia_dr3_important += survey_name_gaia + '.pmdec as pm_dec, '
    # gaia_dr3_important += survey_name_gaia + '.pmdec_error as pm_dec_err, '
    # gaia_dr3_important += survey_name_gaia + '.pmra as pm_ra, '
    # gaia_dr3_important += survey_name_gaia + '.pmra_error as pm_ra_err, '
    # gaia_dr3_important += survey_name_gaia + '.radial_velocity as vr, '
    # gaia_dr3_important += survey_name_gaia + '.radial_velocity_error vr_err '
    #    gaia_dr3_important += survey_name_gaia + '.classprob_dsc_combmod_galaxy as classprob_dsc_combmod_galaxy, '
    #    gaia_dr3_important += survey_name_gaia + '.classprob_dsc_combmod_quasar as classprob_dsc_combmod_quasar, '
    #    gaia_dr3_important += survey_name_gaia + '.classprob_dsc_combmod_star as classprob_dsc_combmod_star, '
    #    gaia_dr3_important += survey_name_gaia + '.grvs_mag as g_integ_gaia, '
    #    gaia_dr3_important += survey_name_gaia + '.grvs_mag_error as g_integ_gaia_err '
    # gaia_dr3_important += survey_name_gaia + '.phot_g_mean_flux as phot_g_mean_flux_gaia, '
    # gaia_dr3_important += survey_name_gaia + '.phot_g_mean_flux_error as phot_g_mean_flux_error_gaia, '
    # gaia_dr3_important += survey_name_gaia + '.phot_g_mean_mag as phot_g_mean_mag_gaia '
    #    gaia_dr3_important += survey_name_gaia + '.phot_rp_mean_flux as phot_rp_mean_flux_gaia, '
    #    gaia_dr3_important += survey_name_gaia + '.phot_rp_mean_flux_error as phot_rp_mean_flux_error_gaia, '
    #    gaia_dr3_important += survey_name_gaia + '.phot_rp_mean_mag as phot_rp_mean_mag_gaia '
    search_survey(ra_s_min, ra_s_max, dec_s_min, dec_s_max, survey_name_decals, name_s,  decals_dr10_important)

