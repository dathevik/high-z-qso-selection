#----------- Making histogram from χ2 ratios calculated in sed_calculation.py 4MOST project script
#----------- Analysing results of the testing and improvements

from astropy.io import ascii, fits
import matplotlib.pyplot as plt
import os

object_file = 'results.fits'
output_path = os.path.abspath('output_test')
objects_path = os.path.abspath('output_test/' + object_file)
hdu = fits.open(objects_path)
data_results = hdu[1].data
BD_Chi2 = data_results.field(3)
QSO_Chi2 = data_results.field(5)
R_Chi2 = data_results.field(7)


# ------- Histogram for χ2 ratios
bins = [0.3, 0.6, 1, 1.5, 2.5, 4]
hist_ratio = plt.hist(R_Chi2, bins=bins,  edgecolor='black')
plt.xlabel("χ2 Ratio")
plt.ylabel("Sum of corresponding objects")
plt.title("Histogram of χ2 ratio")
plt.savefig(f"{output_path}/hist_r_chi2.png")
plt.close()


# -------------------------------------------------------------------------------------------------------------
# ------- Histogram for QSO χ2
hist_qso = plt.hist(QSO_Chi2, bins=6, edgecolor='black')
plt.xlabel("χ2 of QSO")
plt.ylabel("Sum of corresponding objects")
plt.title("Histogram of QSO χ2 ")
plt.savefig(f"{output_path}/hist_qso.png")
plt.close()

# ------- Histogram for BD χ2
hist_bd = plt.hist(BD_Chi2, bins=6, edgecolor='black')
plt.xlabel("χ2 of BD")
plt.ylabel("Sum of corresponding objects")
plt.title("Histogram of BD χ2")
plt.savefig(f"{output_path}/hist_bd.png")
plt.close()


# ------- Plot QSO_chi2 with R_chi2
plt.scatter(QSO_Chi2, R_Chi2, s=10)
plt.xlabel("χ2 of QSO")
plt.ylabel("χ2 Ratio")
# plt.title("Histogram of BD χ2")
plt.savefig(f"{output_path}/r_chi2_qso_chi2.png")
plt.close()