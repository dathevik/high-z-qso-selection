import shutil

import numpy as np
import os
from astropy.io import ascii, fits
from astropy.table import Table



# ------------------------------- Define functions ---------------------------------

# function for a_scale which is a scaling factor and takes as an input of observed flux, observed flux error of the source and the template model
def a_scale(vec_flux_obs, vec_fluxe_obs, vec_flux_model):
    # ---- Obtain scaling factor for Chi2
    a = np.sum((vec_flux_obs * vec_flux_model) / (vec_fluxe_obs) ** 2) / np.sum(
        (vec_flux_model) ** 2 / (vec_fluxe_obs) ** 2)
    return a


# function for calculation of chi2 statistical parameter and takes as an input of observed flux, observed flux error of the source and the template model and a scaling factor
def chi2_calc(vec_flux_obs, vec_fluxe_obs, vec_flux_model, a):
    # ---- Obtain scaling factor for Chi2
    chi2 = np.sum((vec_flux_obs - a * vec_flux_model) ** 2 / (vec_fluxe_obs) ** 2)
    return chi2


# function for converting cgs parameters to mJy as the whole script works with mJy units of fluxes
def flux_from_cgs_to_mJy(flux_cgs, ll):
    # ---- Obtain flux in mJy from flux in erg/s/cm2/A and considering effective wavelength of the filter in AA
    flux_mJy = (3.33564095E+04 * flux_cgs * ll ** 2) * 1E-3
    return flux_mJy


# ------------------------------- Main code ---------------------------------
# I: Define the path and print the data in the created directory
obj_path = os.path.abspath('input_test/ARAA_decals_vhs_wise_fluxes.dat')
BD_temp_path = os.path.abspath('input_test/BDRA_fluxes_mJy_202301.dat')
QSO_temp_path = os.path.abspath('input_test/Selsing+Matt_temp.dat')
QSO_spec_path = os.path.abspath('input_test/Selsing2015.dat')
output_folder = os.path.join(os.path.abspath(''), 'output_test')
#if os.path.exists(output_folder):
 #   shutil.rmtree(output_folder)
#os.mkdir(output_folder)

# ----------------- Printing and reading files -----------------
print("READING BD:", BD_temp_path, "FILE")
data_bd_temp = ascii.read(BD_temp_path)
# print(data_bd_temp)


print("READING QSOs:", QSO_temp_path, "FILE")
data_qso_temp = ascii.read(QSO_temp_path)
# print(data_qso_temp)

print("READING INPUT CATALOG:", obj_path, "FILE")
data_obj = ascii.read(obj_path)
# print(data_obj[0])

print("READING SELSING 2015:", QSO_spec_path, "FILE")
data_qso_spec = ascii.read(QSO_spec_path)
# print(data_qso_spec)

# IIa : Make the vector arrays for Brown Dwarf templates
BD_g_decam = data_bd_temp.columns[1]
BD_r_decam = data_bd_temp.columns[2]
BD_i_decam = data_bd_temp.columns[3]
BD_z_decam = data_bd_temp.columns[4]
BD_Y_vhs = data_bd_temp.columns[5]
BD_J_vhs = data_bd_temp.columns[6]
BD_H_vhs = data_bd_temp.columns[7]
BD_K_vhs = data_bd_temp.columns[8]
BD_W1 = data_bd_temp.columns[9]
BD_W2 = data_bd_temp.columns[10]

# ----------------- Read BD template rows in a loop -----------------
BD_vec_flux = []
for i in range(len(data_bd_temp)):
    vec_flux = [BD_g_decam[i], BD_r_decam[i], BD_i_decam[i], BD_z_decam[i], BD_Y_vhs[i], BD_J_vhs[i], BD_H_vhs[i],
                BD_K_vhs[i], BD_W1[i], BD_W2[i]]
    BD_vec_flux.append(vec_flux)

# ----------------- Make list with BDs type and print first template as an example -----------------
BD_all_vec_flux = np.array(BD_vec_flux)
BD_type_vec = ["BD1", "BD2", "BD3", "BD4", "MD1", "MD2", "MD3"]

print("------------------------------------------ BD1 Template----------------------------------------------------")
print(BD_all_vec_flux[0])


# IIb : Make vector arrays for Quasars templates
QSO_z_vec = data_qso_temp.columns[1]
QSO_g_decam = data_qso_temp.columns[2]
QSO_r_decam = data_qso_temp.columns[3]
QSO_i_decam = data_qso_temp.columns[4]
QSO_z_decam = data_qso_temp.columns[5]
QSO_Y_vhs = data_qso_temp.columns[6]
QSO_J_vhs = data_qso_temp.columns[7]
QSO_H_vhs = data_qso_temp.columns[8]
QSO_K_vhs = data_qso_temp.columns[9]
QSO_W1 = data_qso_temp.columns[10]
QSO_W2 = data_qso_temp.columns[11]
QSO_temp_type = data_qso_temp.columns[0]
QSO_temp_emline = data_qso_temp.columns[12]
QSO_temp_ebv = data_qso_temp.columns[13]

# ----------------- Read Column of Redshift in a loop and print the first template as an example -----------------
QSO_vec_flux = []
for i in range(len(data_qso_temp)):
    vec_flux = [QSO_g_decam[i], QSO_r_decam[i], QSO_i_decam[i], QSO_z_decam[i], QSO_Y_vhs[i], QSO_J_vhs[i],
                QSO_H_vhs[i], QSO_K_vhs[i], QSO_W1[i], QSO_W2[i]]
    QSO_vec_flux.append(vec_flux)

QSO_all_vec_flux = np.array(QSO_vec_flux)
print("------------------------------------------ QSO1 Template----------------------------------------------------")
print(QSO_all_vec_flux[0])

# IV : Read wavebands, fluxes and errors from source file
header = list(data_obj.columns)
select_err = "err"
del header[:4]
errors = [i for i in header if select_err in i]
wavebands = [i for i in header if i not in errors]
print("-----------------------------------------WAVEBANDS ------------------------------------------------------")
print(wavebands)
obj_name = data_obj.columns[0]
obj_redshift = data_obj["redshift"]

# IVa : Make array for central fluxes of the bands(each of them)
vec_flux = data_obj[
    ["g_prime_delve", "r_prime_delve", "i_prime_delve", "z_prime_delve", "vista.vircam.Y_vhs", "vista.vircam.J_vhs",
     "vista.vircam.H_vhs", "vista.vircam.Ks_vhs", "WISE1", "WISE2"]].copy()
# print("-----------------------------------------FLUXES ------------------------------------------------------")
# print(vec_flux)

# IVb : Make array for errors of the bands(each of them)
vec_fluxe = data_obj[
    ["g_prime_err_delve", "r_prime_err_delve", "i_prime_err_delve", "z_prime_err_delve", "vista.vircam.Y_err_vhs",
     "vista.vircam.J_err_vhs", "vista.vircam.H_err_vhs", "vista.vircam.Ks_err_vhs", "WISE1_err", "WISE2_err"]].copy()
# print("------------------------------------------ WAVEBAND ERRORS----------------------------------------------------")
# print(errors)

# ----------------- Changing flux errors if they are smaller 10 times than the flux -----------------
for i in range(len(vec_flux.columns)):
    x = vec_flux.columns[i]
    y = vec_fluxe.columns[i]
    for j in range(len(data_obj)):
        if y[j] < 0.1 * x[j]:
            y[j] = 0.1 * x[j]
        if x[j] <= 10 ** -20:
            x[j] = np.nan
            y[j] = np.nan

# c. ----------------- Make all arrays with the same type in this case float -----------------
vec_flux_model_BD = BD_all_vec_flux.astype(float)
vec_flux_model_QSO = QSO_all_vec_flux.astype(float)

# ----------------- Making empty lists to write new calculations in it -----------------
ls_id_l = []
BD_min_vec = []
QSO_min_vec =[]
ra = data_obj.columns[1]
dec = data_obj.columns[2]
BD_Chi2_temp_l = []
QSO_Chi2_temp_l = []
QSO_emline_l = []
QSO_ebv_l =[]
QSO_temp_t = []
ratio = []

# V : Calculate Chi2 for BD and QSO templates
for i in range(len(data_obj)):
# for i in range(3):
    print("----------------------------------------- RUN FOR OBJECT NUMBER", i, "----------------------------------------")

    BD_Chi2_array = []
    QSO_Chi2_array = []
    R_Chi = []
    a_BD_array = []
    a_QSO_array = []

    vec_flux_row = vec_flux[i]
    vec_fluxe_row = vec_fluxe[i]
    vec_flux_obs0 = np.array(
        [vec_flux_row[0], vec_flux_row[1], vec_flux_row[2], vec_flux_row[3], vec_flux_row[4], vec_flux_row[5],
         vec_flux_row[6], vec_flux_row[7], vec_flux_row[8], vec_flux_row[9]])
    vec_fluxe_obs0 = np.array(
        [vec_fluxe_row[0], vec_fluxe_row[1], vec_fluxe_row[2], vec_fluxe_row[3], vec_fluxe_row[4], vec_fluxe_row[5],
         vec_fluxe_row[6], vec_fluxe_row[7], vec_fluxe_row[8], vec_fluxe_row[9]])

    # ----------------- mask the nan value, which means it will not be used for the calculation -----------------
    mask_nan = ~np.isnan(vec_flux_obs0)
    vec_flux_obs = vec_flux_obs0[mask_nan]
    vec_fluxe_obs = vec_fluxe_obs0[mask_nan]
    # print(vec_flux_obs)

    # Va : Calculate scaling factor and Chi2 for BDs
    for j in range(len(data_bd_temp)):
        vec_flux_model_BD_j0 = vec_flux_model_BD[j]
        vec_flux_model_BD_j = vec_flux_model_BD_j0[mask_nan]
        a_BD = a_scale(vec_flux_obs, vec_fluxe_obs, vec_flux_model_BD_j)
        a_BD = a_BD.astype(float)
        Chi2_BD = chi2_calc(vec_flux_obs, vec_fluxe_obs, vec_flux_model_BD_j, a_BD)
        BD_Chi2_array.append(Chi2_BD)
        a_BD_array.append(a_BD)
    # print(BD_Chi2_array)

    # Vb : Calculate scaling factor and Chi2 for QSOs
    for k in range(len(data_qso_temp)):
        vec_flux_model_QSO_k0 = vec_flux_model_QSO[k]
        vec_flux_model_QSO_k = vec_flux_model_QSO_k0[mask_nan]
        a_QSO = a_scale(vec_flux_obs, vec_fluxe_obs, vec_flux_model_QSO_k)
        a_QSO = a_QSO.astype(float)
        Chi2_QSO = chi2_calc(vec_flux_obs, vec_fluxe_obs, vec_flux_model_QSO_k, a_QSO)
        QSO_Chi2_array.append(Chi2_QSO)
        a_QSO_array.append(a_QSO)
    # print(QSO_Chi2_array)

    # VI : Calculate and print out the template with the best (lowest) Chi2 for QSOs and BDs
    # VIa : Calculate the lowest value for the Chi2 for the BDs
    BD_Chi2_min = np.min(BD_Chi2_array)  # -- minimum value
    BD_Chi2_min_ind = np.argmin(BD_Chi2_array)  # -- position in array
    a_BD_best = a_BD_array[BD_Chi2_min_ind]  # -- corresponding scaling factor of the best chi2
    BD_Chi2_min_temp = BD_type_vec[BD_Chi2_min_ind]  # -- corresponding best template
    vec_flux_model_BD_best = BD_all_vec_flux[BD_Chi2_min_ind][mask_nan] * a_BD_best
    print("-------------------------------------------------")
    print("Best Chi2 BD:", BD_Chi2_min)
    print("Which BD templates:", BD_Chi2_min_temp)

    # VIb. Calculate the lowest value for the Chi2 for the QSOs
    QSO_Chi2_min = np.min(QSO_Chi2_array)  # -- minimum value
    QSO_Chi2_min_ind = np.argmin(QSO_Chi2_array)  # -- position in array
    a_QSO_best = a_QSO_array[QSO_Chi2_min_ind]  # -- corresponding scaling factor of the best chi2
    QSO_Chi2_min_z = QSO_z_vec[QSO_Chi2_min_ind]  # -- corresponding best template -> best redshift
    vec_flux_model_QSO_best = QSO_all_vec_flux[QSO_Chi2_min_ind][mask_nan] * a_QSO_best
    temp_ebv_best = QSO_temp_ebv[QSO_Chi2_min_ind]
    temp_type_best = QSO_temp_type[QSO_Chi2_min_ind]
    temp_emline_best = QSO_temp_emline[QSO_Chi2_min_ind]
    print("-------------------------------------------------")
    print("Best Chi2 QSO:", QSO_Chi2_min)
    print("Which QSO redshift:", QSO_Chi2_min_z)
    print("Which QSO template:", temp_type_best)
    print("Which template emline:", temp_emline_best)
    print("Which template EBV:", temp_ebv_best)
    # ----Calculate the ratio of chi2 of each QSO and BD templates
    R = QSO_Chi2_min / BD_Chi2_min
    R_Chi.append(R)
    print("-------------------------------------------------")
    print("The Chi2 ratio is", R)
    print(R_Chi)

    # ---- Collect results of the calculation in one array and sort it in descending order of χ2_Ratio
    ls_id_l.append([obj_name[i]])
    BD_min_vec.append(BD_Chi2_min)
    QSO_min_vec.append(QSO_Chi2_min)
    BD_Chi2_temp_l.append(BD_Chi2_min_temp)
    QSO_Chi2_temp_l.append(QSO_Chi2_min_z)
    QSO_ebv_l.append(temp_ebv_best)
    QSO_emline_l.append(temp_emline_best)
    QSO_temp_t.append(QSO_temp_type)
    ratio.append(R)


ls_id_c = fits.Column(name='Object Name', array=ls_id_l, format='A10')
ra_c = fits.Column(name='ra', array=ra, format='F')
dec_c = fits.Column(name='dec', array=dec, format='F')
BD_min_vec_c = fits.Column(name='BD_chi2_min', array=BD_min_vec, format='F')
QSO_min_vec_c = fits.Column(name='QSO_chi2_min', array=QSO_min_vec, format='F')
BD_Chi2_temp_c = fits.Column(name='BD chi2 template', array=BD_Chi2_temp_l, format='A10')
QSO_Chi2_temp_c = fits.Column(name='QSO z', array=QSO_Chi2_temp_l, format='F')
QSO_temp_type_c = fits.Column(name='QSO temptype', array=QSO_temp_t, format='A10')
QSO_ebv_c = fits.Column(name='QSO EBV', array=QSO_ebv_l, format='F')
QSO_emline_c = fits.Column(name='QSO  EMline', array=QSO_emline_l, format='F')
R_chi_c = fits.Column(name='R_chi2_best', array=ratio, format='F')
t_out = fits.BinTableHDU.from_columns([ls_id_c, ra_c, dec_c, BD_min_vec_c, BD_Chi2_temp_c, QSO_min_vec_c, QSO_Chi2_temp_c, R_chi_c, QSO_ebv_c, QSO_emline_c])
name_t_out= os.path.abspath('output_test/results.fits')
t_out.writeto(name_t_out, overwrite=True)


