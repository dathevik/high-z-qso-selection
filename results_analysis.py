#----------- Making histogram from χ2 ratios calculated in sed_calculation.py 4MOST project script
#----------- Analysing results of the testing and improvements

from astropy.io import ascii
import matplotlib.pyplot as plt
import os

object_file = 'results_01.csv'
output_path = os.path.abspath('output_test')
objects_path = os.path.abspath('output_test/' + object_file)
data_results = ascii.read(objects_path)
BD_Chi2 = data_results.columns[2]
QSO_Chi2 = data_results.columns[4]
R_Chi2 = list(data_results.columns[6])

# ------- Histogram for χ2 ratios
bins = [0.3, 0.6, 1, 2.5, 4, 5]
hist_ratio = plt.hist(R_Chi2, bins=bins,  edgecolor='black')
plt.xlabel("χ2 Ratio")
plt.ylabel("Sum of corresponding objects")
plt.title("Histogram of χ2 ratio")

for i in range(5): # ---- calculating the sum of similar outputs
    if hist_ratio[0][i]!=0:
        x = float(hist_ratio[0][i])
        y = float(len(R_Chi2))
        per = round(x/y * 100)
        # plt.text(hist_bd_r[1][i],hist_bd_r[0][i],f"{per}%")
plt.savefig(f"{output_path}/hist_r_chi2.png")
plt.close()


# -------------------------------------------------------------------------------------------------------------
# ------- Histogram for QSO χ2
hist_qso = plt.hist(QSO_Chi2, bins=6, edgecolor='black')
plt.xlabel("χ2 of QSO")
plt.ylabel("Sum of corresponding objects")
plt.title("Histogram of QSO χ2 ")

for i in range(5): # ---- calculating the sum of similar outputs
    if hist_qso[0][i]!=0:
        x = float(hist_qso[0][i])
        y = float(len(R_Chi2))
        per = round(x / y * 100)
        # plt.text(hist_qso[1][i], hist_qso[0][i], f"{per}%")
plt.savefig(f"{output_path}/hist_qso.png")
plt.close()

# ------- Histogram for BD χ2
hist_bd = plt.hist(BD_Chi2, bins=6, edgecolor='black')
plt.xlabel("χ2 of BD")
plt.ylabel("Sum of corresponding objects")
plt.title("Histogram of BD χ2")

for i in range(5): # ---- calculating the sum of similar outputs
    if hist_bd[0][i]!=0:
        x = float(hist_bd[0][i])
        y = float(len(R_Chi2))
        per = round(x / y * 100)
        # plt.text(hist_qso[1][i], hist_qso[0][i], f"{per}%")
plt.savefig(f"{output_path}/hist_bd.png")
plt.close()
# ------- Plot QSO_chi2 with R_chi2
plt.scatter(QSO_Chi2, R_Chi2, s=10)
plt.xlabel("χ2 of QSO")
plt.ylabel("χ2 Ratio")
# plt.title("Histogram of BD χ2")
plt.savefig(f"{output_path}/r_chi2_qso_chi2.png")
plt.close()