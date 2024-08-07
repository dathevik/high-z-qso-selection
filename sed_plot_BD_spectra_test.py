import shutil

import numpy as np
from decimal import Decimal
import matplotlib.pyplot as plt
import os
from astropy.io import ascii
from astropy.table import Table, QTable, hstack

from synphot import SourceSpectrum, etau_madau
from synphot.models import Empirical1D

# ------------------------------- Define functions ---------------------------------

# function for a_scale which is a scaling factor and takes as an input of observed flux, observed flux error of the object and the template model
def a_scale(vec_flux_obs, vec_fluxe_obs, vec_flux_model):
    # ---- Obtain scaling factor for Chi2
    a = np.nansum((vec_flux_obs * vec_flux_model) / (vec_fluxe_obs) ** 2) / np.nansum(
        (vec_flux_model) ** 2 / (vec_fluxe_obs) ** 2)
    return a


# function for calculation of chi2 statistical parameter and takes as an input of observed flux, observed flux error of the object and the template model and a scaling factor
def chi2_calc(vec_flux_obs, vec_fluxe_obs, vec_flux_model, a):
    # ---- Obtain scaling factor for Chi2
    chi2 = np.nansum((vec_flux_obs - a * vec_flux_model) ** 2 / (vec_fluxe_obs) ** 2)
    return chi2


# function for converting cgs parameters to mJy as the whole script works with mJy units of fluxes
def flux_from_cgs_to_mJy(flux_cgs, ll):
    # ---- Obtain flux in mJy from flux in erg/s/cm2/micron and considering effective wavelength of the filter in AA
    flux_mJy = (3.33564095E+04 * flux_cgs * ll ** 2) * 1E-3
    return flux_mJy

def plot_files_from_folder(folder_path, color_map, label, alpha, legend_color, empirical_fluxes):
    # Get the list of all files in the folder
    file_list = os.listdir(folder_path)

    # Filter the list to include only .txt files
    txt_files = [f for f in file_list if f.endswith('.txt')]

    # Create a color palette using a color map
    num_files = len(txt_files)
    colors = color_map(np.linspace(0, 1, num_files))

    subtype_fluxes_and_colors = []

    # Plot all files with transparency and dashed lines using the color palette
    for i, file_name in enumerate(txt_files):
        file_path = os.path.join(folder_path, file_name)
        # Extract spectral type from filename
        spectral_type = file_name.split('.')[0]

        # Read the file using astropy.io.ascii
        data = ascii.read(file_path)
        # Assuming the first and second columns are the ones to plot
        x_data = data.columns[0] * 10000
        y_data_cgs = data.columns[1] *0.013616591668795114*1E-8

        # Plot the data with transparency, dashed lines, and color from the palette
        plt.plot(x_data, y_data_cgs, linestyle='--', alpha=alpha, color=colors[i])

    # Plot an invisible point to create a legend entry
    plt.plot([], [], linestyle='--', color=legend_color, alpha=alpha, label=label)

    return subtype_fluxes_and_colors

def get_template_flux(template_name, return_wavelength=False):
    # Define the base folder paths for different types
    base_folders = {
        'L': 'input_test/BD_templates_empirical_CM/BDspectra/LDwarf/specs',
        'T': 'input_test/BD_templates_empirical_CM/BDspectra/TDwarf/specs',
        'M': 'input_test/BD_templates_empirical_CM/BDspectra/MDwarf/specs'
        }

    # Determine the appropriate folder based on the template name
    type_prefix = template_name[0]  # e.g., 'L' or 'T' or 'M'
    if type_prefix in base_folders:
        folder_path = base_folders[type_prefix]
    else:
        raise ValueError(f"Unknown template type for {template_name}")

    # Construct the file path
    file_path = os.path.join(folder_path, f"{template_name}.txt")

    # Check if the file exists
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found")

    # Read the file using astropy.io.ascii
    data = ascii.read(file_path)

     # Assuming the first and second columns are the ones to plot
    x_data = data.columns[0] * 10000  # Convert wavelength from microns to Ångströms
    y_data_cgs = data.columns[1] * 0.013616591668795114 * 1E-8  # Convert flux

    if return_wavelength:
        return x_data, y_data_cgs
    else:
        return y_data_cgs

    # Example usage:
    # template_name = 'L1'
    # try:
    #     y_data_mJy = get_template_flux(template_name)
    #     print(f"Flux values for {template_name} (in mJy): {y_data_mJy}")
    #
    #     x_data, y_data_mJy_with_x = get_template_flux(template_name, return_wavelength=True)
    #     print(f"Wavelength and flux values for {template_name}: {x_data}, {y_data_mJy_with_x}")
    # except Exception as e:
    #     print(e)

# ------------------------------- Main code ---------------------------------
# I: Define the path and print the data in the created directory
BD_temp_path = os.path.abspath('input_test/BDRACM_fluxes_mJy.dat')
data_bd_temp = ascii.read(BD_temp_path)
mdRA_path = os.path.abspath('input_test/lrt_kc04_MS.dat')
bdRA_path = os.path.abspath('input_test/lrt_a07_BDs.dat')
QSO_temp_path = os.path.abspath('input_test/Selsing+Matt_temp.dat')
data_qso_temp = ascii.read(QSO_temp_path)
QSO_spec_path = os.path.abspath('input_test/Selsing2015.dat')
data_qso_spec = ascii.read(QSO_spec_path)
output_folder = os.path.join(os.path.abspath(''), 'output_test')
if os.path.exists(output_folder):
    shutil.rmtree(output_folder)
os.mkdir(output_folder)
output_file = os.path.abspath('output_test/results.csv')


# Input the information
ls_id = input("Input the name of object").strip()
ls_id = str(ls_id)
input_file_name = input("Input the name of catalog ")
input_objects_path = os.path.abspath('input_test/'+input_file_name)

print("READING INPUT CATALOG:", input_objects_path, "FILE")
data_fil = ascii.read(input_objects_path)
# print(data_fil)
id = data_fil.columns[0]
redshift = data_fil.columns[3]

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
#
# # Read BD template rows in a loop
BD_vec_flux = []
for i in range(len(data_bd_temp)):
    vec_flux = [BD_g_decam[i], BD_r_decam[i], BD_i_decam[i], BD_z_decam[i], BD_Y_vhs[i], BD_J_vhs[i], BD_H_vhs[i],
                BD_K_vhs[i], BD_W1[i], BD_W2[i]]
    BD_vec_flux.append(vec_flux)

# # Make list with BDs type and print first template as an example
BD_all_vec_flux = np.array(BD_vec_flux)
vec_flux_model_BD = BD_all_vec_flux.astype(float)
BD_type_vec = ["BD1", "BD2", "BD3", "BD4", "MD1", "MD2", "MD3", "L0", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "M4", "M5", "M6", "M7", "M8", "M9", "T0", "T1", "T2", "T3", "T4", "T6", "T8"]
BD_type_RA = ["BD1", "BD2", "BD3", "BD4", "MD1", "MD2", "MD3"]

 # -------- Read RAssef BD and MD templates for plotting and comparing the results
mdRA_data = Table(ascii.read(mdRA_path))
bdRA_data = Table(ascii.read(bdRA_path))
bdRA_l1 = bdRA_data.columns[0] * 10 ** 4  # --AA
bdRA_l2 = bdRA_data.columns[1] * 10 ** 4  # --AA
bdRA_l = bdRA_l2 + (bdRA_l1 - bdRA_l2) / 2.  # ---bin center in the  BD template
mdRA_l1 = mdRA_data.columns[0] * 10 ** 4  # --AA
mdRA_l2 = mdRA_data.columns[1] * 10 ** 4  # --AA
mdRA_l = mdRA_l2 + (mdRA_l1 - mdRA_l2) / 2.  # ---bin center in the MD template

# # ----- Fluxes vectors in nuFnu (Hz Jy)
bdRA_nf1 = bdRA_data.columns[2]
bdRA_nf2 = bdRA_data.columns[3]
bdRA_nf3 = bdRA_data.columns[4]
bdRA_nf4 = bdRA_data.columns[5]
mdRA_nf1 = mdRA_data.columns[33]
mdRA_nf2 = mdRA_data.columns[34]
mdRA_nf3 = mdRA_data.columns[35]
#
# # ----- Convert fluxes from nuFnu (Hz Jy) in (erg/s/cm2/A)
bdRA_f1 = bdRA_nf1
bdRA_f2 = bdRA_nf2
bdRA_f3 = bdRA_nf3
bdRA_f4 = bdRA_nf4
mdRA_f1 = mdRA_nf1
mdRA_f2 = mdRA_nf2
mdRA_f3 = mdRA_nf3
bd_ra_array = QTable([bdRA_f1, bdRA_f2, bdRA_f3, bdRA_f4, mdRA_f1, mdRA_f2, mdRA_f3], names=(BD_type_RA))
# List of template names
template_names = [
    "L0", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9",
    "M4", "M5", "M6", "M7", "M8", "M9",
    "T0", "T1", "T2", "T3", "T4", "T6", "T8"
]

# Get the fluxes for each template
bd_cm_fluxes = [get_template_flux(name) for name in template_names]

# Determine the minimum length among the flux arrays
min_length = len(bd_ra_array)  # Assuming bd_ra_array has consistent length across columns

# Trim all flux arrays in bd_cm_fluxes to the minimum length
bd_cm_fluxes_trimmed = [flux[:min_length] for flux in bd_cm_fluxes]

# Create bd_cm_array with trimmed fluxes
bd_cm_array = QTable(bd_cm_fluxes_trimmed, names=template_names[:min_length])

# Combine bd_ra_array and bd_cm_array into a single QTable
bd_all_array = hstack([bd_ra_array, bd_cm_array])

# Print the combined QTable
print(bd_all_array)
#
ll_g = 4798.3527009231575  # Angstrom
ll_r = 6407.493598028656  # Angstrom
ll_i = 7802.488114833454  # Angstrom
ll_z = 9144.625340022629  # Angstrom
ll_J = 12325.125694338809  # Angstrom
ll_Y = 10201.359507821942  # Angstrom
ll_H = 16473.95843628733  # Angstrom
ll_K = 22045.772662096875  # Angstrom
ll_W1 = 33791.878497259444  # Angstrom
ll_W2 = 46292.93969033106  # Angstrom # -----Define central wavelength for plot it later for check

ll_vec = np.array([ll_g, ll_r, ll_i, ll_z, ll_Y, ll_J, ll_H, ll_K, ll_W1, ll_W2])
#
# # IIb : Make vector arrays for Quasars templates
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

# # Initialize a plot
# plt.figure(figsize=(10, 6))
#
# # Plotting all L-dwarfs with a specific alpha and legend color
# folder_path_LDwarf = os.path.abspath('input_test/BD_templates_empirical_CM/BDspectra/LDwarf/specs')
# fluxes_LDwarf_empirical = os.path.abspath('input_test/BD_templates_empirical_CM/BDspectra/L_type_wise_g_r.txt')
# subtype_fluxes_and_colors_LDwarf = plot_files_from_folder(folder_path_LDwarf, plt.cm.Purples, 'LDwarf', alpha=0.3, legend_color='purple', empirical_fluxes=fluxes_LDwarf_empirical)
#
# # Plotting all T-dwarfs with a different alpha and legend color
# folder_path_TDwarf = os.path.abspath('input_test/BD_templates_empirical_CM/BDspectra/TDwarf/specs')
# fluxes_TDwarf_empirical = os.path.abspath('input_test/BD_templates_empirical_CM/BDspectra/T_type_wise_g_r.txt')
# subtype_fluxes_and_colors_TDwarf = plot_files_from_folder(folder_path_TDwarf, plt.cm.Oranges, 'TDwarf', alpha=0.3, legend_color='orange', empirical_fluxes=fluxes_TDwarf_empirical)
#
# # Plotting all M-dwarfs with a different alpha and legend color
# folder_path_MDwarf = os.path.abspath('input_test/BD_templates_empirical_CM/BDspectra/MDwarf/specs')
# fluxes_MDwarf_empirical = os.path.abspath('input_test/BD_templates_empirical_CM/BDspectra/M_type_wise_g_r.txt')
# subtype_fluxes_and_colors_MDwarf = plot_files_from_folder(folder_path_MDwarf, plt.cm.Greens, 'MDwarf', alpha=0.3, legend_color='green', empirical_fluxes=fluxes_MDwarf_empirical)
#
# # Function to plot empirical fluxes
# def plot_empirical_fluxes(subtype_fluxes_and_colors):
#     for subtype, color, g_flux, r_flux, w1_flux, w2_flux in subtype_fluxes_and_colors:
#         plt.plot(ll_g, g_flux, 'x', color=color, alpha=0.5)
#         plt.plot(ll_r, r_flux, 'x', color=color, alpha=0.5)
#         plt.plot(ll_W1, w1_flux, 'x', color=color, alpha=0.5)
#         plt.plot(ll_W2, w2_flux, 'x', color=color, alpha=0.5)
#
# # Plot empirical fluxes for all subtypes
# plot_empirical_fluxes(subtype_fluxes_and_colors_LDwarf)
# plot_empirical_fluxes(subtype_fluxes_and_colors_TDwarf)
# plot_empirical_fluxes(subtype_fluxes_and_colors_MDwarf)
#
# # Add labels and legend
# plt.xlim(0, 50000)
# plt.ylim(0,1.5)
# plt.legend(loc='upper left')
# plt.tight_layout()
#
#
# ----- Read Column of Redshift in a loop and print the first template as an example
QSO_vec_flux = []

for i in range(len(data_qso_temp)):
    vec_flux = [QSO_g_decam[i], QSO_r_decam[i], QSO_i_decam[i], QSO_z_decam[i], QSO_Y_vhs[i], QSO_J_vhs[i],
                QSO_H_vhs[i], QSO_K_vhs[i], QSO_W1[i], QSO_W2[i]]
    QSO_vec_flux.append(vec_flux)
QSO_all_vec_flux = np.array(QSO_vec_flux)
vec_flux_model_QSO = QSO_all_vec_flux.astype(float)

#
# IVa : Make array for central fluxes of the bands(each of them)
vec_flux = data_fil[
    ["g_prime_delve", "r_prime_delve", "i_prime_delve", "z_prime_delve", "vista.vircam.Y_vhs", "vista.vircam.J_vhs",
     "vista.vircam.H_vhs", "vista.vircam.Ks_vhs", "WISE1", "WISE2"]].copy()
# print(vec_flux)
#
# IVb : Make array for errors of the bands(each of them)
vec_fluxe = data_fil[
    ["g_prime_err_delve", "r_prime_err_delve", "i_prime_err_delve", "z_prime_err_delve", "vista.vircam.Y_err_vhs",
     "vista.vircam.J_err_vhs", "vista.vircam.H_err_vhs", "vista.vircam.Ks_err_vhs", "WISE1_err", "WISE2_err"]].copy()
# print(vec_fluxe)

object_name = []
data = []
object_redshift = []
vec_flux_obs = []
vec_fluxe_obs = []

print(f"CHECKING {ls_id} OBJECT IN {input_file_name} CATALOG ")

for i in range(len(id)):
    if ls_id == str(id[i]):
        print(f"{ls_id} OBJECT is FOUND in {input_file_name} CATALOG")
        object_name.append(id[i])
        data.append(data_fil[i])
        object_redshift.append(redshift[i])
        x = vec_flux[i]
        y = vec_fluxe[i]
        for j in range(len(x)):
                if y[j] < 0.1 * x[j]:
                    y[j] = 0.1 * x[j]
                if x[j] <= 10 ** -30:
                    x[j] = np.nan
                vec_flux_obs.append(x[j])
                vec_fluxe_obs.append(y[j])

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

        # ---- mask the nan value, which means it will not be used for the calculation
        mask_nan = ~np.isnan(vec_flux_obs0)
        vec_flux_obs = vec_flux_obs0[mask_nan]
        vec_fluxe_obs = vec_fluxe_obs0[mask_nan]
        ll_vec_best = ll_vec[mask_nan]
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

        # Vb : Calculate scaling factor and Chi2 for QSOs
        for k in range(len(data_qso_temp)):
            vec_flux_model_QSO_k0 = vec_flux_model_QSO[k]
            vec_flux_model_QSO_k = vec_flux_model_QSO_k0[mask_nan]
            a_QSO = a_scale(vec_flux_obs, vec_fluxe_obs, vec_flux_model_QSO_k)
            a_QSO = a_QSO.astype(float)
            Chi2_QSO = chi2_calc(vec_flux_obs, vec_fluxe_obs, vec_flux_model_QSO_k, a_QSO)
            QSO_Chi2_array.append(Chi2_QSO)
            a_QSO_array.append(a_QSO)

        # VI : Calculate and print out the template with the best (lowest) Chi2 for QSOs and BDs
        # VIa : Calculate the lowest value for the Chi2 for the BDs
        BD_Chi2_min = np.min(BD_Chi2_array)  # -- minimum value
        BD_Chi2_min_ind = np.argmin(BD_Chi2_array)  # -- position in array
        a_BD_best = a_BD_array[BD_Chi2_min_ind]  # -- corresponding scaling factor of the best chi2
        BD_Chi2_min_temp = BD_type_vec[BD_Chi2_min_ind]  # -- corresponding best template
        vec_flux_model_BD_best = BD_all_vec_flux[BD_Chi2_min_ind][mask_nan] * a_BD_best
        print("Best BD Scaling:", a_BD_best)

        #Preparing all other plots of BD's to be plotted
        templates = {
            'BD1': BD_all_vec_flux[0][mask_nan] * a_BD_array[0],
            'BD2': BD_all_vec_flux[1][mask_nan] * a_BD_array[1],
            'BD3': BD_all_vec_flux[2][mask_nan] * a_BD_array[2],
            'BD4': BD_all_vec_flux[3][mask_nan] * a_BD_array[3],
            'MD1': BD_all_vec_flux[4][mask_nan] * a_BD_array[4],
            'MD2': BD_all_vec_flux[5][mask_nan] * a_BD_array[5],
            'MD3': BD_all_vec_flux[6][mask_nan] * a_BD_array[6],
            'L0': BD_all_vec_flux[7][mask_nan] * a_BD_array[7],
            'L1': BD_all_vec_flux[8][mask_nan] * a_BD_array[8],
            'L2': BD_all_vec_flux[9][mask_nan] * a_BD_array[9],
            'L3': BD_all_vec_flux[10][mask_nan] * a_BD_array[10],
            'L4': BD_all_vec_flux[11][mask_nan] * a_BD_array[11],
            'L5': BD_all_vec_flux[12][mask_nan] * a_BD_array[12],
            'L6': BD_all_vec_flux[13][mask_nan] * a_BD_array[13],
            'L7': BD_all_vec_flux[14][mask_nan] * a_BD_array[14],
            'L8': BD_all_vec_flux[15][mask_nan] * a_BD_array[15],
            'L9': BD_all_vec_flux[16][mask_nan] * a_BD_array[16],
            'M4': BD_all_vec_flux[17][mask_nan] * a_BD_array[17],
            'M5': BD_all_vec_flux[18][mask_nan] * a_BD_array[18],
            'M6': BD_all_vec_flux[19][mask_nan] * a_BD_array[19],
            'M7': BD_all_vec_flux[20][mask_nan] * a_BD_array[20],
            'M8': BD_all_vec_flux[21][mask_nan] * a_BD_array[21],
            'M9': BD_all_vec_flux[22][mask_nan] * a_BD_array[22],
            'T0': BD_all_vec_flux[23][mask_nan] * a_BD_array[23],
            'T1': BD_all_vec_flux[24][mask_nan] * a_BD_array[24],
            'T2': BD_all_vec_flux[25][mask_nan] * a_BD_array[25],
            'T3': BD_all_vec_flux[26][mask_nan] * a_BD_array[26],
            'T4': BD_all_vec_flux[27][mask_nan] * a_BD_array[27],
            'T6': BD_all_vec_flux[28][mask_nan] * a_BD_array[28],
            'T8': BD_all_vec_flux[29][mask_nan] * a_BD_array[29]
        }

        # VIb. Calculate the lowest value for the Chi2 for the QSOs
        QSO_Chi2_min = np.min(QSO_Chi2_array)  # -- minimum value
        QSO_Chi2_min_ind = np.argmin(QSO_Chi2_array)  # -- po11036800232067sition in array
        a_QSO_best = a_QSO_array[QSO_Chi2_min_ind]  # -- corresponding scaling factor of the best chi2
        QSO_Chi2_min_z = QSO_z_vec[QSO_Chi2_min_ind]  # -- corresponding best template -> best redshift
        vec_flux_model_QSO_best = QSO_all_vec_flux[QSO_Chi2_min_ind][mask_nan] * a_QSO_best
        temp_ebv_best = QSO_temp_ebv[QSO_Chi2_min_ind]
        temp_type_best = QSO_temp_type[QSO_Chi2_min_ind]
        temp_emline_best = QSO_temp_emline[QSO_Chi2_min_ind]
    # flux_stripe82_test_03.dat
        # ----Calculate the ratio of chi2 of each QSO and BD templates
        R = QSO_Chi2_min / BD_Chi2_min
        R_Chi.append(R)
        print("Best Chi2 QSO:", QSO_Chi2_min)
        print("Best QSO Redshift:", QSO_Chi2_min_z)
        print("Best Chi2 BD:", BD_Chi2_min)
        print("Best BD Template:", BD_Chi2_min_temp)
        print("Best Chi2 Ratio (Chi2_QSO/Chi2_BD):", R)

        # ---- parameters for plotting in a loop using Matt or S15 templates
        z = QSO_Chi2_min_z
        if temp_type_best == "S15":
            w_rest = np.arange(900., 75000, 10)
            spec_f = np.array(data_qso_spec)
            wave = data_qso_spec.columns[0]
            flux = data_qso_spec.columns[1]
            sp = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux, keep_neg=True)
            sp_z_chi2 = SourceSpectrum(sp.model, z=z)
            extcurve_best = etau_madau(w_rest, z)
            sp_ext_best = sp_z_chi2 * extcurve_best
            flux_mJY_best = flux_from_cgs_to_mJy(sp_ext_best(w_rest), w_rest) * a_QSO_best
            if BD_Chi2_min_ind > 3:
                bdRA_l_fin = mdRA_l
            else:
                bdRA_l_fin = bdRA_l

            for label, vec_flux_model in templates.items():
                if label != BD_Chi2_min_temp:
                    index = list(templates.keys()).index(label)
                    bdRA_f_nonbest = (bd_all_array.columns[index] / bdRA_l_fin) * a_BD_array[index]
                    bdRA_fmJy_nonbest = flux_from_cgs_to_mJy(bdRA_f_nonbest, bdRA_l_fin)
                    # plt.plot(bdRA_l_fin, bdRA_fmJy_nonbest, alpha=0.4,
                    #           label=f'{label} Synthetic Template', linewidth=2.0,
                    #           zorder=0,linestyle='--')
            bdRA_f_best = (bd_all_array.columns[BD_Chi2_min_ind] / bdRA_l_fin) * a_BD_best
            bdRA_fmJy_best = flux_from_cgs_to_mJy(bdRA_f_best, bdRA_l_fin)
            plt.scatter(ll_vec_best, vec_flux_model_BD_best, color="#836853",
                        label=f"Best BD template ({BD_Chi2_min_temp})",
                        zorder=5, facecolor='none', linewidth=2.0 )
            plt.scatter(ll_vec_best, vec_flux_model_QSO_best, color="#004987",
                        label=f"Best QSO template (z={QSO_Chi2_min_z})",
                        zorder=5, facecolor='none', linewidth=2.0)
            plt.errorbar(ll_vec_best, vec_flux_obs, yerr=vec_fluxe_obs, fmt='o', color='#FCD12A',
                         label=f'Observed spectra with z={object_redshift}')
            plt.plot(w_rest, flux_mJY_best, color="#a7bed3", label=f'Quasar spectrum with z={QSO_Chi2_min_z}', zorder=1, linewidth=2.0)
            plt.plot(bdRA_l_fin, bdRA_fmJy_best, color="#dab894", label=f'Brown Dwarf spectrum with {BD_Chi2_min_temp}', linewidth=2.0,
                     zorder=1)
            plt.legend(loc="upper right")
            plt.xlabel("Central Wavelength (Å)")
            plt.ylabel("Flux (mJy)")
            plt.title(f"SED fitting of {ls_id} object", fontsize=9)
            plt.savefig(f"output_test/{ls_id}_BD_QSO_best_models.png")
            plt.show()
            plt.close()

        else:
            w_rest = np.arange(900., 18900, 10)
            z = round(QSO_Chi2_min_z * 100)
            ebv_dec = round((Decimal(temp_ebv_best)), 2)
            spec_ebv = str(ebv_dec).replace(".", "")
            spec_file = (f"input_test/seds2/qsogen_sed_emlines{temp_emline_best}_ebv+{spec_ebv}_z+{z}.dat")
            spec_open = os.path.abspath(spec_file)
            data_spec_Matt = ascii.read(spec_file)
            wave = data_spec_Matt.columns[0]
            flux = data_spec_Matt.columns[1]
            # sp = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux, keep_neg=True)
            # sp_z_chi2 = SourceSpectrum(sp.model, z=redshift)
            # extcurve_best = etau_madau(w_rest, redshift)
            # sp_ext_best = sp_z_chi2 * extcurve_best
            flux_mJY_best = flux_from_cgs_to_mJy(flux, wave) * a_QSO_best
            # print("THIS IS a_QSO_best",a_QSO_best)
            # print(flux_mJY_best)
            if BD_Chi2_min_ind > 3:
                bdRA_l_fin = mdRA_l
            else:
                bdRA_l_fin = bdRA_l

            for label, vec_flux_model in templates.items():
                if label != BD_Chi2_min_temp:
                    index = list(templates.keys()).index(label)
                    bdRA_f_nonbest = (bd_all_array.columns[index] / bdRA_l_fin) * a_BD_array[index]
                    bdRA_fmJy_nonbest = flux_from_cgs_to_mJy(bdRA_f_nonbest, bdRA_l_fin)
                    # plt.plot(bdRA_l_fin, bdRA_fmJy_nonbest, alpha=0.4,
                    #          label=f'{label} Synthetic Template', linewidth=2.0,
                    #          zorder=0, linestyle='--')
            bdRA_f_best = (bd_all_array.columns[BD_Chi2_min_ind] / bdRA_l_fin) * a_BD_best
            bdRA_fmJy_best = flux_from_cgs_to_mJy(bdRA_f_best, bdRA_l_fin)
            plt.scatter(ll_vec_best, vec_flux_model_BD_best, color="#836853",
                        label=f"Best BD template ({BD_Chi2_min_temp})",
                        zorder=5, facecolor='none')
            plt.scatter(ll_vec_best, vec_flux_model_QSO_best, color="#004987",
                        label=f"Best QSO template (z={QSO_Chi2_min_z})",
                        zorder=5, facecolor='none')
            plt.errorbar(ll_vec_best, vec_flux_obs, yerr=vec_fluxe_obs, fmt='o', color='#FCD12A',
                         label=f'Initial Flux')
            plt.plot(wave, flux_mJY_best*1e6, color="#a7bed3", zorder=0, linewidth=2.0)
            plt.plot(bdRA_l_fin, bdRA_fmJy_best, color="#dab894",zorder=0, linewidth=2.0)
            plt.legend(loc="upper left")
            plt.xlabel("Central Wavelength (Å)")
            plt.ylabel("Flux (mJy)")
            plt.title(f"SED fitting of {ls_id} object", fontsize=9)
            plt.savefig(f"output_test/{ls_id}_BD_QSO_best_models.png")
            plt.show()
            plt.close()





