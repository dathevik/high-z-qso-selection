import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

# Read data
data_qso_araa = ascii.read('araa_decals.dat')
data_qso_milliquas = ascii.read('milliquas_decals.dat')
data_qso_yang = ascii.read('yang_decals.dat')
data_bd = ascii.read('bd_decals.dat', encoding='latin-1')

g_mag_yang = data_qso_yang["mag_g"]
g_mag_milliquas = data_qso_milliquas["mag_g"]
g_mag_araa = data_qso_araa["mag_g"]
g_mag_bd = data_bd["mag_g"]

g_snr_yang = data_qso_yang["snr_g"]
g_snr_milliquas = data_qso_milliquas["snr_g"]
g_snr_araa = data_qso_araa["snr_g"]
g_snr_bd = data_bd["snr_g"]

r_mag_yang = data_qso_yang["mag_r"]
r_mag_milliquas = data_qso_milliquas["mag_r"]
r_mag_araa = data_qso_araa["mag_r"]
r_mag_bd = data_bd["mag_r"]

i_mag_yang = data_qso_yang["mag_i"]
i_mag_milliquas = data_qso_milliquas["mag_i"]
i_mag_araa = data_qso_araa["mag_i"]
i_mag_bd = data_bd["mag_i"]

z_mag_yang = data_qso_yang["mag_z"]
z_mag_milliquas = data_qso_milliquas["mag_z"]
z_mag_araa = data_qso_araa["mag_z"]
z_mag_bd = data_bd["mag_z"]

redshift_yang = data_qso_yang["redshift"]
redshift_milliquas = data_qso_milliquas["Z"]
redshift_araa = data_qso_araa["redshift"]


# G band vs SNR G
plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(g_mag_bd,g_snr_bd, label="Brown Dwarfs", s=8, alpha=0.4, marker="o", color="gray")
plt.scatter(g_mag_milliquas, g_snr_milliquas, label="Quasars Flesch+2023", s=16, alpha=0.4, marker="x", color="blue")
plt.scatter(g_mag_yang, g_snr_yang, label="Quasars Yang+2023", s=16, alpha=0.4, marker="x", color="red")
plt.scatter(g_mag_araa, g_snr_araa, label="Quasars Fan+2022", s=16, alpha=0.4, marker="x", color="green")
plt.xlabel("DECALS g_mag")
plt.ylabel("Signal to Noise")
plt.legend(loc="upper left")
plt.xlim(18, 30)
plt.ylim(0, 100)
# plt.show()
plt.savefig("BD_QSO_mag_g_snr_decals.png")
plt.close()

# G-I vs redshift
plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(redshift_milliquas, g_mag_milliquas-i_mag_milliquas, label="Quasars Flesch+2023", s=8, alpha=0.4, marker="o")
plt.scatter(redshift_yang, g_mag_yang-i_mag_yang, label="Quasars Yang+2023", s=8, alpha=0.4, marker="o")
plt.scatter(redshift_araa, g_mag_araa-i_mag_araa, label="Quasars Fan+2022", s=8, alpha=0.4, marker="o")
plt.xlabel("Redshift")
plt.ylabel("g-i")
plt.legend(loc="upper left")
plt.xlim(4, 7)
plt.ylim(-2, 12)
# plt.show()
plt.savefig("BD_QSO_g_i_decals.png")
plt.close()

# G-R vs redshift
plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(redshift_milliquas, g_mag_milliquas-r_mag_milliquas, label="Quasars Flesch+2023", s=8, alpha=0.4, marker="o")
plt.scatter(redshift_yang, g_mag_yang-r_mag_yang, label="Quasars Yang+2023", s=8, alpha=0.4, marker="o")
plt.scatter(redshift_araa, g_mag_araa-r_mag_araa, label="Quasars Fan+2022", s=8, alpha=0.4, marker="o")
plt.xlabel("Redshift")
plt.ylabel("g-r")
plt.legend(loc="upper left")
plt.xlim(4, 7)
plt.ylim(-2, 12)
# plt.show()
plt.savefig("BD_QSO_g_r_decals.png")
plt.close()

# G-Z vs redshift
plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(redshift_milliquas, g_mag_milliquas-z_mag_milliquas, label="Quasars Flesch+2023", s=8, alpha=0.4, marker="o")
plt.scatter(redshift_yang, g_mag_yang-z_mag_yang, label="Quasars Yang+2023", s=8, alpha=0.4, marker="o")
plt.scatter(redshift_araa, g_mag_araa-z_mag_araa, label="Quasars Fan+2022", s=8, alpha=0.4, marker="o")
plt.xlabel("Redshift")
plt.ylabel("g-z")
plt.legend(loc="upper left")
plt.xlim(4, 7)
plt.ylim(-2, 12)
# plt.show()
plt.savefig("BD_QSO_g_z_decals.png")
plt.close()

# G-I vs G-R
plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(g_mag_bd -i_mag_bd, g_mag_bd - z_mag_bd, label="Brown Dwarfs", s=8, alpha=0.4, marker="o")
plt.scatter(g_mag_milliquas-i_mag_milliquas, g_mag_milliquas-z_mag_milliquas, label="Quasars Flesch+2023", s=8, alpha=0.4, marker="o")
plt.scatter(g_mag_yang-i_mag_yang, g_mag_yang-z_mag_yang, label="Quasars Yang+2023", s=8, alpha=0.4, marker="o")
plt.scatter(g_mag_araa-i_mag_araa, g_mag_araa-z_mag_araa, label="Quasars Fan+2022", s=8, alpha=0.4, marker="o")
plt.xlabel("g-i")
plt.ylabel("g-z")
plt.legend(loc="upper left")
plt.xlim(0, 12)
plt.ylim(0, 12)
# plt.show()
plt.savefig("BD_QSO_g_z_g_i_decals.png")
plt.close()

# G-I vs snr_g
plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(g_mag_bd -i_mag_bd, g_snr_bd, label="Brown Dwarfs", s=8, alpha=0.4, marker="o")
plt.scatter(g_mag_milliquas-i_mag_milliquas, g_snr_milliquas, label="Quasars Flesch+2023", s=8, alpha=0.4, marker="o")
plt.scatter(g_mag_yang-i_mag_yang, g_snr_yang, label="Quasars Yang+2023", s=8, alpha=0.4, marker="o")
plt.scatter(g_mag_araa-i_mag_araa, g_snr_araa, label="Quasars Fan+2022", s=8, alpha=0.4, marker="o")
plt.xlabel("g-i")
plt.ylabel("SNR_g")
plt.axhline(y=3, color='gray', linestyle='--', label='SNR < 3 cut')
plt.axvline(x=2, color='red', linestyle='--', label='SNR >3 and g-i > 2 cut')
plt.legend(loc="upper left")
plt.xlim(0, 12)
plt.ylim(0, 25)
# plt.show()
plt.savefig("BD_QSO_snr_g_i_cut.png")
plt.close()

# G-Z vs snr_g
plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(g_mag_bd - z_mag_bd, g_snr_bd, label="Brown Dwarfs", s=8, alpha=0.4, marker="o")
plt.scatter(g_mag_milliquas-z_mag_milliquas, g_snr_milliquas, label="Quasars Flesch+2023", s=8, alpha=0.4, marker="o")
plt.scatter(g_mag_yang-z_mag_yang, g_snr_yang, label="Quasars Yang+2023", s=8, alpha=0.4, marker="o")
plt.scatter(g_mag_araa-z_mag_araa, g_snr_araa, label="Quasars Fan+2022", s=8, alpha=0.4, marker="o")
plt.axhline(y=3, color='gray', linestyle='--', label='SNR < 3 cut')
plt.axvline(x=2, color='red', linestyle='--', label='SNR >3 and g-z > 2 cut')
plt.xlabel("g-z")
plt.ylabel("SNR_g")
plt.legend(loc="upper left")
plt.xlim(0, 12)
plt.ylim(0, 25)
# plt.show()
plt.savefig("BD_QSO_snr_g_z_cut.png")
plt.close()