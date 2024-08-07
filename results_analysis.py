#----------- Making histogram from χ2 ratios calculated in sed_calculation.py 4MOST project script
#----------- Analysing results of the testing and improvements

from astropy.io import ascii, fits
import matplotlib.pyplot as plt
import os
import numpy as np

object_file = 'results_4most/results_files/prioritization_delve_decals_cuts.fits'
output_path = os.path.abspath('/output_test')
objects_path = os.path.abspath(object_file)
hdu = fits.open(objects_path)
data_results = hdu[1].data
ra = data_results.field("ra")
dec = data_results.field("dec")
BD_Chi2 = data_results.field("BD_chi2_min")
BD_best_temp = data_results.field("BD chi2 template")
QSO_Chi2 = data_results.field("QSO_chi2_min")
QSO_z = data_results.field("QSO_z")
R_Chi2 = data_results.field("R_chi2_best")
F_test = data_results.field("F_test_value")
BIC = data_results.field("BIC_value")
EBV = data_results.field("ebv")
w1_AB = data_results.field("w1mpro") + 2.699
w2_AB = data_results.field("w2mpro") + 3.399
r_mag = data_results.field("mag_auto_r")
i_mag = data_results.field("mag_auto_i")
z_mag = data_results.field("mag_auto_z")
ext_z = data_results.field("extended_class_z")
data_points = data_results.field("Data_Points")

table_decals = data_results[data_results['mag_type'] == 'DECALS DR10']
mag_g_decals = table_decals.field("mag_auto_g")
snr_g_decals = table_decals.field("snr_g")
BIC_decals = table_decals.field("BIC_value")
F_test_decals = table_decals.field("F_test_value")
QSO_chi2_decals = table_decals.field("QSO_chi2_min")
R_chi2_decals = table_decals.field("R_chi2_best")
table_delve = data_results[data_results['mag_type'] == 'DELVE DR2']
mag_g_delve = table_delve.field("mag_auto_g")
BIC_delve = table_delve.field("BIC_value")
F_test_delve = table_delve.field("F_test_value")
QSO_chi2_delve = table_delve.field("QSO_chi2_min")
R_chi2_delve = table_delve.field("R_chi2_best")
print(data_results)

# ------- Scatter plot of ra and dec ratios
plt.figure(figsize=(10, 6), dpi=160)
distr = plt.scatter(QSO_chi2_delve, R_chi2_delve,  edgecolor='black', s=7, label="QSO Candidates of the Final Catalog DELVE", color="orange")
plt.xlabel("QSO chi2 minimum")
plt.ylabel("R chi2 best")
plt.axhline(y=0.3, color='gray', linestyle='--', label='Statistical cut')
# square_x = [0, 1.2, 1.2, 0, 0]
# # square_y = [0, 0, 0.3, 0.3, 0]
# # plt.plot(square_x, square_y, color='gray', linestyle='--', label='Statistical cut')
plt.legend(loc="upper right")
plt.ylim(0, 0.6)
plt.xlim(0, 25)
# # plt.yticks([0, 10, 25, 50, 75, 100, 125, 150, 175, 200])
# plt.show()
# plt.close()

# ------- Histogram of BD Best template
# plt.title("Priority 1 with DECALS")
# plt.xlabel("mag z")
# plt.ylabel("Number")
# plt.hist(z_mag, bins=150, alpha=0.5, edgecolor='black', color="gray", label="QSO Candidates of the Final Catalog")
# plt.xlim(19.5, 22.5)
# plt.legend(loc="upper right")
# plt.show()

# ra and dec distribution
# plt.figure(figsize=(10,6), dpi=160)
# plt.scatter(ra, dec, c=EBV, s=7)
# cb = plt.colorbar()
# cb.set_label('E(B-V)', fontsize=12)  # Add a label to the colorbar
# plt.xlabel("RA")
# plt.ylabel("DEC")
# plt.clim(0, 0.4)
# plt.show()