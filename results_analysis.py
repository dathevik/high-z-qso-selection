#----------- Making histogram from χ2 ratios calculated in quasar-selection.py 4MOST project script
#----------- Analysing results of the testing and improvements

from astropy.io import ascii
import matplotlib.pyplot as plt


data_bd_v2 = ascii.read("output_tables/results_bd_v2.csv")
data_qso_v2 = ascii.read("output_tables/results_qso_v2.csv")
data_bd_v1 = ascii.read("output_tables/results_bd_v1.csv")
data_qso_v1 = ascii.read("output_tables/results_qso_v1.csv")

BD_Chi2_v1 = data_bd_v1.columns[0]
BD_Chi2_v2 = data_bd_v2.columns[1]
QSO_Chi2_v1 = data_qso_v1.columns[2]
QSO_Chi2_v2 = data_qso_v2.columns[3]

R_Chi2_bd_v1 = list(data_bd_v1.columns[2]/BD_Chi2_v1)
R_Chi2_bd_v2 = list(data_bd_v2.columns[5])
R_Chi2_qso_v1 = list(QSO_Chi2_v1/data_qso_v1.columns[0])
R_Chi2_qso_v2 = list(data_qso_v2.columns[5])

def subtracting_two_ratios(l1,l2):
    result = [x-y for x, y in zip(l1,l2)]
    return result

def dividing_two_ratios(l1,l2):
    result = [x/y for x, y in zip(l1,l2)]
    return result

# R_dif_bd = subtracting_two_ratios(R_Chi2_bd_v2,R_Chi2_bd_v1)
# R_dif_qso = subtracting_two_ratios(R_Chi2_qso_v2,R_Chi2_qso_v1)

# R_eff_bd = dividing_two_ratios(R_Chi2_bd_v2,R_Chi2_bd_v1)
# R_eff_qso = dividing_two_ratios(R_Chi2_qso_v2,R_Chi2_qso_v1)


# ------- Histogram for BD χ2 ratios
bins = [0.3, 0.6, 1, 2.5, 4, 5]
hist_bd_r = plt.hist(R_Chi2_bd_v1, bins=bins,  edgecolor='black')
plt.xlabel("χ2 Ratio")
plt.ylabel("Sum of corresponding objects")
plt.title("Histogram of χ2 ratio for Brown Dwarfs v.1")

for i in range(5): # ---- calculating the sum of similar outputs
    if hist_bd_r[0][i]!=0:
        x = float(hist_bd_r[0][i])
        y = float(len(R_Chi2_bd_v2))
        per = round(x/y * 100)
        plt.text(hist_bd_r[1][i],hist_bd_r[0][i],f"{per}%")
plt.savefig("hist_bd_r_chi2_v1.png")
plt.close()

hist_bd_r = plt.hist(R_Chi2_bd_v1, bins=bins, label="Brown Dwarfs",  edgecolor='black')
hist_qso_r = plt.hist(R_Chi2_qso_v1, bins=bins,  label="Quasars", edgecolor='black')
plt.xlabel("χ2 Ratio")
plt.ylabel("Sum of corresponding objects")
plt.title("Comparison of χ2 ratios v.1")
# plt.savefig("hist_qso_bd_r_chi2_v1.png")
plt.close()

# ------- Histogram for QSO χ2 ratios
hist_qso_r = plt.hist(R_Chi2_qso_v2, bins=bins, edgecolor='black')
plt.xlabel("χ2 ratio ")
plt.ylabel("Sum of corresponding objects")
plt.title("Histogram of χ2 ratio for Quasars v.2")

for i in range(5): # ---- calculating the sum of similar outputs
    if hist_qso_r[0][i]!=0:
        x = float(hist_qso_r[0][i])
        y = float(len(R_Chi2_qso_v1))
        per = round(x / y * 100)
        plt.text(hist_qso_r[1][i], hist_qso_r[0][i], f"{per}%")
plt.savefig("hist_qso_r_chi2_v2.png")
plt.close()


# -------------------------------------------------------------------------------------------------------------
# ------- Histogram for QSO χ2 V1 for QSO data
hist_qso = plt.hist(QSO_Chi2_v1, bins=6, edgecolor='black')
plt.xlabel("χ2 of QSO")
plt.ylabel("Sum of corresponding objects")
plt.title("Histogram of QSO χ2 for Quasars Data v.1")

for i in range(5): # ---- calculating the sum of similar outputs
    if hist_qso[0][i]!=0:
        x = float(hist_qso[0][i])
        y = float(len(R_Chi2_qso_v1))
        per = round(x / y * 100)
        plt.text(hist_qso[1][i], hist_qso[0][i], f"{per}%")
plt.savefig("hist_qso_temp_chi2_v1.png")
plt.close()

# ------- Histogram for QSO χ2 V2 for QSO temp
hist_qso = plt.hist(QSO_Chi2_v2, bins=6, edgecolor='black')
plt.xlabel("χ2 of QSO")
plt.ylabel("Sum of corresponding objects")
plt.title("Histogram of QSO χ2 for Quasars Data v.2")

for i in range(5): # ---- calculating the sum of similar outputs
    if hist_qso[0][i]!=0:
        x = float(hist_qso[0][i])
        y = float(len(R_Chi2_qso_v2))
        per = round(x / y * 100)
        plt.text(hist_qso[1][i], hist_qso[0][i], f"{per}%")
plt.savefig("hist_qso_temp_chi2_v2.png")
plt.close()

# ------- Histogram for BD χ2 V1 for QSO data
hist_qso_bd = plt.hist(data_qso_v1.columns[0], bins=6, edgecolor='black')
plt.xlabel("χ2 of BD")
plt.ylabel("Sum of corresponding objects")
plt.title("Histogram of BD χ2 for Quasars Data v.1")

for i in range(5): # ---- calculating the sum of similar outputs
    if hist_qso_bd[0][i]!=0:
        x = float(hist_qso_bd[0][i])
        y = float(len(R_Chi2_qso_v1))
        per = round(x / y * 100)
        plt.text(hist_qso_bd[1][i], hist_qso_bd[0][i], f"{per}%")
plt.savefig("hist_qso_bd_chi2_v1.png")
plt.close()

# ------- Histogram for BD χ2 V2 for QSO data
hist_qso_bd = plt.hist(data_qso_v2.columns[1], bins=6, edgecolor='black')
plt.xlabel("χ2 of BD")
plt.ylabel("Sum of corresponding objects")
plt.title("Histogram of BD χ2 for Quasars Data v.2")

for i in range(5): # ---- calculating the sum of similar outputs
    if hist_qso_bd[0][i]!=0:
        x = float(hist_qso_bd[0][i])
        y = float(len(R_Chi2_qso_v2))
        per = round(x / y * 100)
        plt.text(hist_qso_bd[1][i], hist_qso_bd[0][i], f"{per}%")
plt.savefig("hist_qso_bd_chi2_v2.png")
plt.close()


# -------------------------------------------------------------------------------------------------------------
# ------- Histogram for BD χ2 V1 for BD data
print(BD_Chi2_v1)
hist_bd = plt.hist(BD_Chi2_v1, bins=6, edgecolor='black')
plt.xlabel("χ2 of BD")
plt.ylabel("Sum of corresponding objects")
plt.title("Histogram of BD χ2 for Brown Dwarf's Data v.1")

for i in range(5): # ---- calculating the sum of similar outputs
    if hist_bd[0][i]!=0:
        x = float(hist_bd[0][i])
        y = float(len(R_Chi2_bd_v1))
        per = round(x / y * 100)
        plt.text(hist_bd[1][i], hist_bd[0][i], f"{per}%")
plt.savefig("hist_bd_temp_chi2_v1.png")
plt.close()

# ------- Histogram for BD χ2 V2 for BD data
hist_bd = plt.hist(BD_Chi2_v2, bins=6, edgecolor='black')
plt.xlabel("χ2 of BD")
plt.ylabel("Sum of corresponding objects")
plt.title("Histogram of BD χ2 for Brown Dwarf's Data v.2")

for i in range(5): # ---- calculating the sum of similar outputs
    if hist_bd[0][i]!=0:
        x = float(hist_bd[0][i])
        y = float(len(R_Chi2_bd_v2))
        per = round(x / y * 100)
        plt.text(hist_bd[1][i], hist_bd[0][i], f"{per}%")
plt.savefig("hist_bd_temp_chi2_v2.png")
plt.close()

# ------- Histogram for QSO χ2 V1 for BD data
hist_bd_qso = plt.hist(data_bd_v1.columns[2], bins=6, edgecolor='black')
plt.xlabel("χ2 of QSO")
plt.ylabel("Sum of corresponding objects")
plt.title("Histogram of QSO χ2 for Brown Dwarf's Data v.1")

for i in range(5): # ---- calculating the sum of similar outputs
    if hist_bd_qso[0][i]!=0:
        x = float(hist_bd_qso[0][i])
        y = float(len(R_Chi2_bd_v1))
        per = round(x / y * 100)
        plt.text(hist_bd_qso[1][i], hist_bd_qso[0][i], f"{per}%")
plt.savefig("hist_bd_qso_chi2_v1.png")
plt.close()

# ------- Histogram for QSO χ2 V2 for BD temp
hist_bd_qso = plt.hist(data_qso_v2.columns[3], bins=6, edgecolor='black')
plt.xlabel("χ2 of QSO")
plt.ylabel("Sum of corresponding objects")
plt.title("Histogram of QSO χ2 for Brown Dwarf's data v.2")

for i in range(5): # ---- calculating the sum of similar outputs
    if hist_bd_qso[0][i]!=0:
        x = float(hist_bd_qso[0][i])
        y = float(len(R_Chi2_bd_v2))
        per = round(x / y * 100)
        plt.text(hist_bd_qso[1][i], hist_bd_qso[0][i], f"{per}%")
plt.savefig("hist_bd_qso_chi2_v2.png")
plt.close()