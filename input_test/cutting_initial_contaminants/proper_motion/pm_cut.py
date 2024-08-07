import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

# Read data
data_qso_milliquas = ascii.read('proper_motion_milliquas.dat')
data_qso_yang = ascii.read('proper_motion_yang.dat')
data_bd = ascii.read('proper_motion_bd.dat')

# Miliquas proper motion
pm_ra_milliquas = data_qso_milliquas["pm_ra"]
pm_ra_err_milliquas = data_qso_milliquas["pm_ra_err"]
pm_dec_milliquas = data_qso_milliquas["pm_dec"]
pm_dec_err_milliquas = data_qso_milliquas["pm_dec_err"]

# Yang proper motion
pm_ra_yang = data_qso_yang["pm_ra"]
pm_ra_err_yang = data_qso_yang["pm_ra_err"]
pm_dec_yang = data_qso_yang["pm_dec"]
pm_dec_err_yang = data_qso_yang["pm_dec_err"]

# BD proper motion
pm_ra_bd = data_bd["pm_ra"]
pm_ra_err_bd = data_bd["pm_ra_err"]
pm_dec_bd = data_bd["pm_dec"]
pm_dec_err_bd = data_bd["pm_dec_err"]


# Scatter plots
plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(pm_ra_bd/pm_ra_err_bd, pm_dec_bd/pm_dec_err_bd, label="Brown Dwarfs", s=20, alpha=0.6, marker="o")
plt.scatter(pm_ra_milliquas/pm_ra_err_milliquas, pm_dec_milliquas/pm_dec_err_milliquas, label="Quasars Flesch+2023", s=10, alpha=0.6, marker="o")
plt.scatter(pm_ra_yang/pm_ra_err_yang, pm_dec_yang/pm_dec_err_yang, label="Quasars Yang+2023", s=10, alpha=0.6, marker="o")

# Labels
plt.xlabel("PM_ra/PM_ra_err")
plt.ylabel("PM_dec/PM_dec_err")

# Drawing the square for proper motion cut
square_x = [-2, 2, 2, -2, -2]
square_y = [-2, -2, 2, 2, -2]
plt.plot(square_x, square_y, color='gray', linestyle='--', label='Proper Motion cut')

# Legend
plt.legend(loc="upper left")

# Limits
plt.ylim(-10, 10)
plt.xlim(-10, 10)

# Show plot
# plt.show()
plt.savefig("QSO_BD_pm.png")